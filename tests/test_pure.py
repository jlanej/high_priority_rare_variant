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
from hprv import genotype as G  # noqa: E402
from hprv.config import get  # noqa: E402
from hprv.ped import parse_ped  # noqa: E402


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
