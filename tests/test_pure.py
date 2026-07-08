"""Pure-logic tests that run on any host (no cyvcf2/scipy/container needed).

Covers the parts of the pipeline that don't need a real VCF: config resolution,
PED parsing, the annotation getters (via a tiny fake variant), genotype QC, and the
Step-6 burden helpers. Integration tests that need cyvcf2/VEP run inside the image.

Run: `python3 tests/test_pure.py`  or  `pytest tests/test_pure.py`
"""
import importlib.util
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from hprv import annotations as A  # noqa: E402
from hprv import audit  # noqa: E402
from hprv import genotype as G  # noqa: E402
from hprv.config import get  # noqa: E402
from hprv.ped import parse_ped, read_trios_file, write_ped  # noqa: E402
from hprv.selection import build_classifier  # noqa: E402


# --- tiny fakes standing in for a cyvcf2 Variant ---------------------------
class _INFO:
    def __init__(self, d):
        self.d = d

    def get(self, k):
        return self.d.get(k)


class FakeVar:
    def __init__(self, info=None, **fmt):
        self.INFO = _INFO(info or {})
        self.CHROM = fmt.get("CHROM", "chr1")
        self.POS = fmt.get("POS", 12345)
        for k in ("gt_quals", "gt_depths", "gt_ref_depths", "gt_alt_depths", "gt_types"):
            setattr(self, k, fmt.get(k))


def test_config_get():
    cfg = {"filters": {"rarity": {"dominant_max": 1e-4}}}
    assert get(cfg, "filters.rarity.dominant_max") == 1e-4
    assert get(cfg, "nope.here", "d") == "d"


def test_config_sh_skips_unresolved():
    """emit_sh must NOT export a literal ${ENV} placeholder — that would shadow the
    shell-level `:=`/`:-` defaults. It emits empty and warns instead."""
    import contextlib
    import io
    from hprv import config as C
    cfg = {"project": {"output_dir": "${NOPE_UNSET}"}, "runtime": {"tmpdir": "/real/tmp"}}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        C.emit_sh(cfg)
    out = buf.getvalue()
    assert "export HPRV_OUTPUT_DIR=''" in out       # unresolved -> empty, defaults win
    assert "${NOPE_UNSET}" not in out               # never leak the literal placeholder to stdout
    assert "HPRV_TMPDIR=/real/tmp" in out            # resolved values pass through (shlex-unquoted)


def test_ped(tmp="/tmp/_hprv_test.ped"):
    with open(tmp, "w") as fh:
        fh.write("# comment\nFAM CHILD DAD MOM 1 2\nFAM DAD 0 0 1 1\nFAM MOM 0 0 2 1\n")
    ped = parse_ped(tmp)
    os.remove(tmp)
    assert ped == {"child": "CHILD", "father": "DAD", "mother": "MOM", "sex": "1"}
    assert parse_ped(None) is None


def test_annotations_frequency_and_predictors():
    v = FakeVar({"hprv_gnomad_faf95": "0.0001", "hprv_gnomad_grpmax_af": "0.002",
                 "vep_REVEL_score": "0.8&0.9", "vep_AlphaMissense_score": "0.7",
                 "vep_SpliceAI_pred_DS_AG": "0.1", "vep_SpliceAI_pred_DS_DL": "0.45"})
    assert abs(A.faf95(v) - 0.0001) < 1e-12
    assert A.frequency(v) == 0.0001                       # faf95 preferred
    assert abs(A.revel(v) - 0.9) < 1e-9                   # max over &-joined
    assert abs(A.spliceai_max(v) - 0.45) < 1e-9           # max over the 4 DS
    assert A.frequency(FakeVar({})) is None               # absent = rarest


def test_annotations_clinvar():
    assert A.clnsig_is_plp(FakeVar({"hprv_clnsig": "Pathogenic"}))
    assert A.clnsig_is_plp(FakeVar({"hprv_clnsig": "Pathogenic/Likely_pathogenic"}))
    assert not A.clnsig_is_plp(FakeVar({"hprv_clnsig": "Conflicting_classifications_of_pathogenicity"}))
    assert not A.clnsig_is_plp(FakeVar({"hprv_clnsig": "Benign"}))
    assert A.clinvar_stars(FakeVar({"hprv_clnrevstat": "reviewed_by_expert_panel"})) == 3
    assert A.clinvar_stars(FakeVar({"hprv_clnrevstat": "no_assertion_criteria_provided"})) == 0
    assert A.clinvar_stars(FakeVar({})) == 0


def test_genotype_qc():
    thr = G.GtThresholds()
    # clean het child: DP 40, GQ 99, AB 0.5
    het = FakeVar(gt_quals=[99], gt_depths=[40], gt_ref_depths=[20], gt_alt_depths=[20])
    assert G.allele_balance(het, 0) == 0.5
    assert G.sample_qc(het, 0, thr, "het")
    # skewed AB fails het band
    skew = FakeVar(gt_quals=[99], gt_depths=[40], gt_ref_depths=[36], gt_alt_depths=[4])
    assert not G.sample_qc(skew, 0, thr, "het")
    # low GQ fails
    lowgq = FakeVar(gt_quals=[10], gt_depths=[40], gt_ref_depths=[20], gt_alt_depths=[20])
    assert not G.sample_qc(lowgq, 0, thr, "het")
    # clean parent: hom-ref, no alt reads
    par = FakeVar(gt_quals=[99], gt_depths=[30], gt_ref_depths=[30], gt_alt_depths=[0])
    assert G.sample_qc(par, 0, thr, "clean_parent")
    # parent with alt reads fails cleanliness
    dirty = FakeVar(gt_quals=[99], gt_depths=[30], gt_ref_depths=[26], gt_alt_depths=[4])
    assert not G.sample_qc(dirty, 0, thr, "clean_parent")


def test_par_x():
    assert G.is_x_nonpar(FakeVar(CHROM="chrX", POS=50_000_000))
    assert not G.is_x_nonpar(FakeVar(CHROM="chrX", POS=1_000_000))   # PAR1
    assert not G.is_x_nonpar(FakeVar(CHROM="chr1", POS=1_000_000))


def test_read_trios_file(tmp="/tmp/_hprv_trios.tsv"):
    # header order must NOT matter: dad/mom located by name, not position
    with open(tmp, "w") as fh:
        fh.write("#kid\tmom\tdad\nCH1\tMO1\tFA1\n")   # note: mom before dad
    trios = read_trios_file(tmp)
    os.remove(tmp)
    assert trios == [("CH1", "FA1", "MO1")]            # returned as (kid, dad, mom)


def test_write_ped_roundtrip(tmp="/tmp/_hprv_gen.ped"):
    write_ped(tmp, "CH1", "FA1", "MO1", kid_sex="2")
    ped = parse_ped(tmp)
    os.remove(tmp)
    assert ped == {"child": "CH1", "father": "FA1", "mother": "MO1", "sex": "2"}


def test_audit_record_and_summarize(tmpdir="/tmp/_hprv_audit"):
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)
    audit.record("01_cohort_sites", "union_sites", 1000, adir=tmpdir)
    audit.record("03_select", "sites_plausible", 120, adir=tmpdir)
    audit.record("04_subset", "candidate_genotypes", 40, scope="CH1", adir=tmpdir)
    audit.record("05_inheritance", "candidate_calls", 3, scope="CH1", adir=tmpdir)
    audit.record("05_inheritance", "mode.denovo", 1, scope="CH1", adir=tmpdir)
    md = audit.summarize(tmpdir)
    shutil.rmtree(tmpdir, ignore_errors=True)
    assert "Global variant funnel" in md and "1000" in md
    assert "plausible sites: 120" in md
    assert "CH1" in md and "denovo=1" in md


def test_step3_classifier():
    cfg = {"filters": {"rarity": {"benign_ba1": 0.05, "recessive_max": 1e-2},
                       "functional": {"revel_pp3_supporting": 0.644,
                                      "alphamissense_lp": 0.564, "spliceai_pp3": 0.2,
                                      "cadd_phred_supporting": 20.0, "mpc_strong": 2.0,
                                      "keep_impacts": ["HIGH", "MODERATE"]}}}
    classify = build_classifier(cfg)
    # BA1-common -> dropped, never rescued
    assert classify(FakeVar({"hprv_gnomad_faf95": "0.2"})) == (False, "ba1")
    # rare + HIGH impact -> kept with reason
    assert classify(FakeVar({"vep_IMPACT": "HIGH"})) == (True, "impact_high")
    # rare + LOFTEE HC -> kept
    assert classify(FakeVar({"vep_LoF": "HC"})) == (True, "loftee_hc")
    # ClinVar P/LP overrides missing function only at >= 2 stars
    keep, why = classify(FakeVar({"hprv_clnsig": "Pathogenic",
                                  "hprv_clnrevstat": "criteria_provided,_multiple_submitters,_no_conflicts"}))
    assert keep and why == "clinvar_plp"
    # a 1-star P/LP does NOT auto-override (single submitter)
    k1, _ = classify(FakeVar({"hprv_clnsig": "Pathogenic",
                              "hprv_clnrevstat": "criteria_provided,_single_submitter"}))
    assert not k1
    # rare but non-functional -> dropped
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER"})) == (False, "not_functional")
    # too common for permissive recessive gate, no P/LP -> dropped
    assert classify(FakeVar({"hprv_gnomad_faf95": "0.02", "vep_IMPACT": "HIGH"})) == (False, "too_common")
    # --- functional-predictor keep branches (reached only when impact is not HIGH/MODERATE) ---
    # SpliceAI Δ >= 0.2 (checked first)
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_SpliceAI_pred_DS_DL": "0.3"})) == (True, "spliceai")
    # REVEL >= 0.644 supporting
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_REVEL_score": "0.7"})) == (True, "revel")
    # AlphaMissense >= 0.564
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_AlphaMissense_score": "0.6"})) == (True, "alphamissense")
    # CADD PHRED >= 25.3
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_CADD_PHRED": "26"})) == (True, "cadd")
    # MPC >= 2.0
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_MPC_score": "2.5"})) == (True, "mpc")
    # sub-threshold predictor -> dropped
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_REVEL_score": "0.5"})) == (False, "not_functional")


def test_contamination():
    import tempfile
    from hprv import contamination as C
    assert C.charr(0, 40) == 0.0                       # clean hom-alt: no ref reads
    assert abs(C.charr(4, 100) - 0.04) < 1e-9          # 4% ref reads = contamination proxy
    assert C.charr(0, 0) is None
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "s.selfSM"), "w") as fh:
        fh.write("#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\tX\n")
        fh.write("SAMP1\tALL\tNA\t1000\t50000\t30\t0.031\tx\n")
    fm = C.read_selfsm(d)
    import shutil
    shutil.rmtree(d, ignore_errors=True)
    assert abs(fm.get("SAMP1", 0) - 0.031) < 1e-9
    assert C.read_selfsm("/no/such/dir") == {}


def test_recurrence_null_per_model():
    """The recurrence null must charge each inheritance model its OWN HWE probability —
    a recessive/hemizygous carrier is not a >=1-of-two-alleles (dominant) event."""
    gb = _load_gb()
    floor, q = 1e-6, 1e-3
    p_dom = gb.p_carrier_hwe([q], floor, 2)      # dominant het: 1-(1-q)^2 ~ 2q
    p_hemi = gb.p_carrier_hwe([q], floor, 1)     # X hemizygous male: ~q
    p_bi = gb.p_biallelic_hwe([q], floor)        # biallelic: ~q^2
    assert abs(p_dom - (1 - (1 - q) ** 2)) < 1e-12
    assert abs(p_hemi - q) < 1e-12
    assert abs(p_bi - q * q) < 1e-12
    # the recessive/hemizygous nulls are FAR smaller than the dominant one (the bug that was fixed)
    assert p_bi < p_hemi < p_dom
    # variants absent from gnomAD use the detection-limit floor, never zero probability
    assert gb.p_carrier_hwe([None], floor, 2) > 0
    assert gb.p_biallelic_hwe([None], floor) > 0


def _load_gb():
    spec = importlib.util.spec_from_file_location(
        "gb", os.path.join(os.path.dirname(__file__), "..", "pipeline", "06_gene_burden.py"))
    gb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gb)
    return gb


def test_burden_helpers():
    gb = _load_gb()
    assert gb.classify("stop_gained") == "lof"
    assert gb.classify("missense_variant") == "missense"
    assert gb.classify("synonymous_variant") == "other"
    q = gb.bh_fdr([0.01, 0.02, 0.03, None, 0.5])
    assert q[3] is None
    assert all(0 <= x <= 1 for x in q if x is not None)
    # monotone non-decreasing in p-value order
    ordered = [q[i] for i in sorted([0, 1, 2, 4], key=lambda i: [0.01, 0.02, 0.03, None, 0.5][i])]
    assert ordered == sorted(ordered)


def _run_all():
    import inspect
    fns = [f for n, f in sorted(globals().items())
           if n.startswith("test_") and inspect.isfunction(f)]
    for f in fns:
        f()
        print(f"PASS {f.__name__}")
    print(f"\nAll {len(fns)} pure-logic tests passed.")


if __name__ == "__main__":
    _run_all()
