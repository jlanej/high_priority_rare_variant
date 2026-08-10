"""Pure-logic tests that run on any host (no cyvcf2/scipy/container needed).

Covers the parts of the pipeline that don't need a real VCF: config resolution,
PED parsing, the annotation getters (via a tiny fake variant), genotype QC, and the
Step-6 burden helpers. Integration tests that need cyvcf2/VEP run inside the image.

Run: `python3 tests/test_pure.py`  or  `pytest tests/test_pure.py`
"""
import importlib.util
import math
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from hprv import annotations as A  # noqa: E402
from hprv import audit  # noqa: E402
from hprv import genotype as G  # noqa: E402
from hprv import igv  # noqa: E402
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
    v = FakeVar({"vep_gnomADe_NFE_AF": "0.002", "vep_gnomADg_EAS_AF": "0.004",
                 "vep_gnomADe_AFR_AF": "0.001", "vep_CADD_PHRED": "12.0&26.5"})
    assert abs(A.grpmax_af(v) - 0.004) < 1e-12            # max across e/g and pops
    PRX = {"resources": {"gnomad": {"oracle": "grpmax_proxy"}}}
    assert A.frequency(v, PRX) == A.grpmax_af(v)          # the proxy arm IS grpmax_af
    # ...and under the DEFAULT (faf95) a cache-only AF is NOT substituted: the oracle has no
    # value for this allele, so it reads absent/rarest. On a correct run that class is empty —
    # the joint slim is a superset of the cache (which only carries dbSNP-accessioned alleles) —
    # so a non-empty one means the transfer under-matched. Step 5 counts and warns about it.
    assert A.frequency(v) is None and A.rarity_basis(v) == "absent"
    assert abs(A.cadd(v) - 26.5) < 1e-9                   # max over &-joined
    assert A.frequency(FakeVar({})) is None               # absent = rarest


def test_frequency_uses_exactly_one_oracle_per_run():
    """ONE frequency oracle for the whole run — the arms must NEVER cross.

    An earlier design preferred faf95 and fell back to the grpmax proxy per variant. Two problems,
    both real: (a) two variants in one run were then compared on different quantities, which is
    undescribable in a methods section; and (b) the fallback direction was wrong — gnomAD emits
    fafmax as MISSING (never as 0) wherever no group's CI clears zero, and of the records with no
    faf95 but a proxy >= 1e-4, 96.5% are AC <= 2, so consulting the point estimate there filtered
    singletons on an inflated AF, exactly what faf95 exists to prevent.

    Both arms are gnomAD v4.1 — the proxy is the VEP cache's own gnomAD AFs — so this is a choice
    of QUANTITY (CI lower bound vs point estimate), not of database.
    """
    from hprv import annotations as AN
    F = AN.F
    FAF = {"resources": {"gnomad": {"oracle": "faf95"}}}
    PRX = {"resources": {"gnomad": {"oracle": "grpmax_proxy"}}}

    def _v(**info):
        class I(dict):
            def get(self, k, d=None): return dict.get(self, k, d)
        o = type("V", (), {})()
        o.INFO = I(info)
        return o

    # the oracle is RUN-level and constant; the default needs no extra resource
    assert AN.rarity_oracle(FAF) == "faf95"
    assert AN.rarity_oracle(PRX) == "grpmax_proxy"
    # faf95 is the DEFAULT — the correct, citable quantity, not the convenient one. It requires
    # the gnomAD joint slim, and run_pipeline.sh HALTS at preflight without it rather than
    # quietly running on the point estimate (the same contract as spliceai_required).
    assert AN.rarity_oracle(None) == "faf95", "faf95 must be the default oracle"
    assert AN.rarity_oracle({}) == "faf95"
    assert AN.rarity_oracle({"resources": {"gnomad": {"oracle": "GRPMAX_PROXY"}}}) == "grpmax_proxy"
    assert AN.rarity_oracle({"resources": {"gnomad": {"oracle": "nonsense"}}}) == "faf95", \
        "an unrecognised oracle must fall to the CORRECT quantity, not the convenient one"

    both = _v(**{F["faf95"]: "6e-05", F["gnomad_af_joint"]: "8e-05",
                 F["gnomade_nfe_af"]: "0.00025"})
    # each arm returns ITS OWN quantity and ignores the other entirely
    assert AN.frequency(both, FAF) == 6e-05
    assert AN.frequency(both, PRX) == 0.00025
    assert AN.rarity_basis(both, FAF) == "measured"
    assert AN.rarity_basis(both, PRX) == "measured"

    # faf95 arm: gnomAD HAS the allele but published no fafmax -> 0.0 (rarest), basis=zero_ci.
    # The proxy is NOT consulted, even though it is present and would have failed the gate.
    singleton = _v(**{F["gnomad_af_joint"]: "8.25e-05", F["gnomade_nfe_af"]: "0.00022"})
    assert AN.frequency(singleton, FAF) == 0.0
    assert AN.rarity_basis(singleton, FAF) == "zero_ci"
    assert AN.frequency(singleton, FAF) < 1.0e-4 <= 0.00022, \
        "the singleton survives the dominant gate on faf95 and would fail on the point estimate"
    # ...and on the proxy arm the SAME variant reads the point estimate, consistently
    assert AN.frequency(singleton, PRX) == 0.00022
    assert AN.rarity_basis(singleton, PRX) == "measured"

    # nothing at all -> None (rarest) on both arms, never a measured zero
    empty = _v()
    assert AN.frequency(empty, FAF) is None and AN.rarity_basis(empty, FAF) == "absent"
    assert AN.frequency(empty, PRX) is None and AN.rarity_basis(empty, PRX) == "absent"

    # the faf95 arm must not read the proxy even when faf95 is absent AND gnomAD has no record
    proxy_only = _v(**{F["gnomade_nfe_af"]: "0.02"})
    assert AN.frequency(proxy_only, FAF) is None, "the faf95 arm consulted the proxy"
    assert AN.frequency(proxy_only, PRX) == 0.02

    # MAX_AF and the global AFs are never the oracle on EITHER arm (golden rule 2)
    for cfg in (FAF, PRX):
        assert AN.frequency(_v(**{F["max_af"]: "0.03"}), cfg) is None, "MAX_AF leaked"
        assert AN.frequency(_v(**{F["gnomade_af"]: "0.03"}), cfg) is None, "a global AF leaked"
    # gnomad_AF_joint is a WITNESS on the faf95 arm, never a value
    assert AN.frequency(_v(**{F["gnomad_af_joint"]: "0.03"}), FAF) == 0.0


def test_prioritize_consumes_the_screens_resolved_rarity():
    """Step 9 must RANK on the value the SCREEN resolved, and must not re-derive it.

    The raw columns for "gnomAD published no faf95" (=> 0, rarest) and "gnomAD has no record"
    (=> absent) are IDENTICAL — faf95 blank either way. So a local re-derivation in Step 9 is
    structurally unable to tell them apart, and would rank a singleton on the point estimate that
    the screen correctly declined to gate it on. Step 5 resolves once; Step 9 consumes.
    """
    from hprv import prioritize as PR
    base = {"consequence": "missense_variant", "impact": "MODERATE", "ref": "A", "alt": "T",
            "inheritance": "dominant", "child_gt": "0/1"}
    # the screen said 0.0 via the faf95 arm; the proxy column disagrees and MUST be ignored
    r = PR.score_variant({**base, "rarity_af": "0", "rarity_oracle": "faf95",
                          "rarity_basis": "zero_ci", "grpmax_af": "0.00022"}, {}, {})
    assert r["rarity_oracle"] == "faf95" and r["rarity_basis"] == "zero_ci"
    assert r["rarity_strength"] == PR.rarity_strength("0", {}), \
        "Step 9 re-derived from the proxy instead of consuming the screen's value"
    # a measured faf95 rides through unchanged
    r = PR.score_variant({**base, "rarity_af": "6e-05", "rarity_oracle": "faf95",
                          "rarity_basis": "measured", "grpmax_af": "0.00025"}, {}, {})
    assert r["rarity_oracle"] == "faf95"
    assert r["rarity_strength"] == PR.rarity_strength("6e-05", {})
    # the proxy arm is carried verbatim too
    r = PR.score_variant({**base, "rarity_af": "9e-05", "rarity_oracle": "grpmax_proxy",
                          "rarity_basis": "measured"}, {}, {})
    assert r["rarity_oracle"] == "grpmax_proxy"
    # legacy Step-8 tables (no rarity_* columns) still score, via the documented fallback path
    r = PR.score_variant({**base, "grpmax_af": "9e-05"}, {}, {})
    assert r["rarity_oracle"] == "grpmax_proxy"
    r = PR.score_variant(dict(base), {}, {})
    assert r["rarity_strength"] == "unknown"


def test_prioritize_nhomalt_conflict_is_reported_and_costs_nothing_by_default():
    """nhomalt flags a biallelic call gnomAD already carries homozygotes for — reported, not
    penalised, because no calibration exists for how many should disqualify one.

    And None (no gnomAD transfer) must never behave like 0 (transfer ran, no homozygotes).
    """
    from hprv import prioritize as PR
    base = {"consequence": "missense_variant", "impact": "MODERATE", "ref": "A", "alt": "T",
            "child_gt": "1/1", "grpmax_af": "1e-4"}
    hit = PR.score_variant({**base, "inheritance": "hom_recessive", "nhomalt": "12"}, {}, {})
    assert hit["nhomalt_recessive_conflict"] is True
    clean = PR.score_variant({**base, "inheritance": "hom_recessive", "nhomalt": "0"}, {}, {})
    assert clean["nhomalt_recessive_conflict"] is False
    absent = PR.score_variant({**base, "inheritance": "hom_recessive"}, {}, {})
    assert absent["nhomalt_recessive_conflict"] is False, "absent nhomalt must not flag"
    # a DOMINANT call is not a recessive conflict however many homozygotes exist
    dom = PR.score_variant({**base, "inheritance": "dominant", "nhomalt": "12"}, {}, {})
    assert dom["nhomalt_recessive_conflict"] is False
    # default charges nothing; the flag is the deliverable
    assert hit["pts_quality"] == clean["pts_quality"]
    # ...but the knob works when set deliberately
    cfg = {"prioritization": {"composite": {"weights": {"quality": {"nhomalt_conflict": -1.5}}}}}
    charged = PR.score_variant({**base, "inheritance": "hom_recessive", "nhomalt": "12"}, {}, cfg)
    assert charged["pts_quality"] == clean["pts_quality"] - 1.5


def test_frequency_excludes_bottlenecked_pops():
    """The whole point of the grpmax proxy: a founder-group-only allele must NOT drive rarity.

    gnomAD's grpmax excludes ami/asj/fin/mid/remaining because their small ANs make a point AF
    unrepresentative. (The FAF group set excludes ami/asj/fin but INCLUDES mid — the one class
    where the two oracles legitimately differ; see test_frequency_uses_exactly_one_oracle_per_run.) VEP's MAX_AF does not exclude them — so if frequency() ever regressed to
    reading MAX_AF, this variant would report 1.1e-3, blow the 1e-4 dominant gate, and a real
    ultra-rare candidate would be silently dropped. It must read None (no eligible group).
    """
    PRX = {"resources": {"gnomad": {"oracle": "grpmax_proxy"}}}
    ami_only = FakeVar({"vep_gnomADg_AMI_AF": "0.0011", "vep_gnomADe_FIN_AF": "0.0009",
                        "vep_gnomADe_ASJ_AF": "0.0015", "vep_gnomADe_MID_AF": "0.002",
                        "vep_MAX_AF": "0.002", "vep_MAX_AF_POPS": "gnomADe_MID"})
    assert A.frequency(ami_only, PRX) is None
    # ...but a real NFE signal on the same variant IS counted.
    plus_nfe = FakeVar({"vep_gnomADg_AMI_AF": "0.0011", "vep_gnomADe_NFE_AF": "3e-5",
                        "vep_MAX_AF": "0.0011"})
    assert abs(A.frequency(plus_nfe, PRX) - 3e-5) < 1e-12
    # On the faf95 arm the guarantee is even stronger and needs no exclusion list: gnomAD computes
    # FAF only over afr/amr/eas/mid/nfe/sas, so an ami/asj/fin-only allele has no eligible group
    # and simply carries no fafmax. MAX_AF is not consulted by EITHER arm.
    assert A.frequency(ami_only) is None, "MAX_AF must not leak into the faf95 arm either"
    assert A.rarity_basis(ami_only) == "absent"


def test_frequency_reads_multivalued_af_tuples():
    """cyvcf2 hands back a TUPLE for a multi-valued numeric INFO field, not a comma-joined string.

    split-vep types every vep_gnomAD{e,g}_<POP>_AF column Number=.,Type=Float, so any selector
    emitting more than one CSQ block (-s all, -s mane at a two-gene locus, -s pick on a multiallelic
    site) yields a tuple. str(tuple) is "(0.001, 0.004)", whose tokens both fail float() — the old
    reader returned None, i.e. "absent from gnomAD => rarest", silently RETAINING a common
    polymorphism. Golden rule 2: frequency() is the rarity oracle, so this must never regress.
    (Every other test here feeds strings, which is exactly why this hid.)
    """
    # Asserted on grpmax_af() directly — this is about _max_float's VALUE parsing, which is the
    # same on either oracle arm, so it must not be re-pinned every time the default changes.
    tup = FakeVar({"vep_gnomADe_NFE_AF": (0.001, 0.004)})
    assert abs(A.grpmax_af(tup) - 0.004) < 1e-12
    # mixed: a tuple with a missing entry, and the max must still win across fields
    mixed = FakeVar({"vep_gnomADe_NFE_AF": (".", 2e-5), "vep_gnomADg_AFR_AF": (7e-5,)})
    assert abs(A.grpmax_af(mixed) - 7e-5) < 1e-12
    # a plain scalar float (single-value field) still works
    assert abs(A.grpmax_af(FakeVar({"vep_gnomADe_SAS_AF": 1.5e-4})) - 1.5e-4) < 1e-12
    # and the string form (what most callers see) is unchanged
    assert abs(A.grpmax_af(FakeVar({"vep_gnomADe_EAS_AF": "0.002&0.003"})) - 0.003) < 1e-12
    # the faf95 field parses the same way (it is Number=A, so a multiallelic site yields a tuple)
    assert abs(A.faf95(FakeVar({"gnomad_faf95": (1e-5, 4e-5)})) - 4e-5) < 1e-12


def test_annotations_clinvar():
    # VEP CLIN_SIG is lowercase and '&'-joined; the ClinVar VCF's CLNSIG was Capitalised and
    # '/'-joined. Both must parse, so the predicate survives either annotation source.
    assert A.clnsig_is_plp(FakeVar({"vep_CLIN_SIG": "pathogenic"}))
    assert A.clnsig_is_plp(FakeVar({"vep_CLIN_SIG": "pathogenic&likely_pathogenic"}))
    assert A.clnsig_is_plp(FakeVar({"vep_CLIN_SIG": "Pathogenic/Likely_pathogenic"}))
    assert not A.clnsig_is_plp(FakeVar({"vep_CLIN_SIG": "conflicting_classifications_of_pathogenicity"}))
    assert not A.clnsig_is_plp(FakeVar({"vep_CLIN_SIG": "conflicting_interpretations_of_pathogenicity"}))
    assert not A.clnsig_is_plp(FakeVar({"vep_CLIN_SIG": "benign"}))
    assert not A.clnsig_is_plp(FakeVar({}))


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
    # Pinned to the grpmax_proxy arm because these fixtures carry VEP-CACHE AFs. On the default
    # (faf95) arm a variant with only a cache AF and no gnomAD joint record correctly reads as
    # absent/rarest, so these BA1 fixtures would not fire — a real behavioural difference, not a
    # test artifact. The faf95 arm's Step-3 behaviour is exercised end-to-end by the integration
    # run, whose mock supplies a real gnomAD slim.
    cfg = {"resources": {"gnomad": {"oracle": "grpmax_proxy"}},
           "filters": {"rarity": {"benign_ba1": 0.05, "recessive_max": 1e-2},
                       "functional": {"cadd_phred_supporting": 20.0, "spliceai_ds_min": 0.2,
                                      "keep_impacts": ["HIGH", "MODERATE"]}}}
    classify = build_classifier(cfg)
    # BA1-common -> dropped, never rescued
    assert classify(FakeVar({"vep_gnomADe_NFE_AF": "0.2"})) == (False, "ba1")
    # rare + HIGH impact -> kept with reason
    assert classify(FakeVar({"vep_IMPACT": "HIGH"})) == (True, "impact_high")
    assert classify(FakeVar({"vep_IMPACT": "MODERATE"})) == (True, "impact_moderate")
    # ClinVar P/LP overrides missing function. No star gate: the VEP cache has no CLNREVSTAT,
    # so an unstarred assertion is all there is and it is honored (over-retains by design).
    assert classify(FakeVar({"vep_CLIN_SIG": "pathogenic"})) == (True, "clinvar_plp")
    # ...and it rescues a variant that is otherwise too common for the permissive gate
    assert classify(FakeVar({"vep_gnomADe_NFE_AF": "0.02",
                             "vep_CLIN_SIG": "likely_pathogenic"})) == (True, "clinvar_plp")
    # ...but never a BA1-common one
    assert classify(FakeVar({"vep_gnomADe_NFE_AF": "0.2",
                             "vep_CLIN_SIG": "pathogenic"})) == (False, "ba1")
    # rare but non-functional -> dropped
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER"})) == (False, "not_functional")
    # too common for permissive recessive gate, no P/LP -> dropped
    assert classify(FakeVar({"vep_gnomADe_NFE_AF": "0.02", "vep_IMPACT": "HIGH"})) == (False, "too_common")
    # --- CADD: the ONLY functional branch, and the only keep-path below MODERATE impact ---
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_CADD_PHRED": "26"})) == (True, "cadd")
    assert classify(FakeVar({"vep_IMPACT": "LOW", "vep_CADD_PHRED": "26"})) == (True, "cadd")
    # sub-threshold -> dropped
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_CADD_PHRED": "12"})) == (False, "not_functional")
    # A MODERATE variant is kept by IMPACT and never consults CADD — so a low CADD cannot
    # drop it. This is why the old missense-predictor branches were unreachable: any variant
    # carrying REVEL/AlphaMissense/MPC is missense => MODERATE => already returned here.
    assert classify(FakeVar({"vep_IMPACT": "MODERATE", "vep_CADD_PHRED": "0.1"})) == (True, "impact_moderate")
    # --- SpliceAI: the splice keep-path, checked BEFORE CADD so a splice hit is labelled 'spliceai' ---
    # deep-intronic MODIFIER with strong splice, low CADD -> kept via spliceai (not cadd, not dropped)
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_CADD_PHRED": "3",
                             "vep_SpliceAI_pred_DS_AL": "0.55"})) == (True, "spliceai")
    # spliceai_ds is the MAX over the four events (here the donor-gain field)
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_SpliceAI_pred_DS_DG": "0.30"})) == (True, "spliceai")
    # below the splice cutoff AND below CADD -> dropped
    assert classify(FakeVar({"vep_IMPACT": "MODIFIER", "vep_CADD_PHRED": "3",
                             "vep_SpliceAI_pred_DS_AL": "0.15"})) == (False, "not_functional")
    # a strong splice signal never rescues a BA1-common variant
    assert classify(FakeVar({"vep_gnomADe_NFE_AF": "0.2",
                             "vep_SpliceAI_pred_DS_AL": "0.9"})) == (False, "ba1")


def test_annotations_spliceai_ds():
    # max over the four delta-score events; None when unscored (never drops, only fails to rescue)
    assert A.spliceai_ds(FakeVar({"vep_SpliceAI_pred_DS_AG": "0.03", "vep_SpliceAI_pred_DS_AL": "0.91",
                                  "vep_SpliceAI_pred_DS_DG": "0.10", "vep_SpliceAI_pred_DS_DL": "0.00"})) == 0.91
    assert A.spliceai_ds(FakeVar({})) is None
    assert A.spliceai_ds(FakeVar({"vep_SpliceAI_pred_DS_AL": "."})) is None


def test_spliceai_backfill_parsing():
    from hprv import spliceai_backfill as SB
    # per-event max over a single gene entry
    b = SB._max_ds("A|GENE1|0.03|0.91|0.10|0.00|-5|10|3|-2")
    assert b["DS_AL"] == 0.91 and b["DS_AG"] == 0.03 and b["DS_DL"] == 0.0
    # per-event MAX across multiple gene entries
    b = SB._max_ds("A|G1|0.10|0.20|0.30|0.40|1|2|3|4,A|G2|0.50|0.05|0.00|0.90|1|2|3|4")
    assert (b["DS_AG"], b["DS_AL"], b["DS_DG"], b["DS_DL"]) == (0.5, 0.2, 0.3, 0.9)
    assert SB._max_ds("") is None and SB._max_ds(None) is None and SB._max_ds("garbage") is None

    class _V:
        def __init__(s, ref, alt): s.REF = ref; s.ALT = [alt]
    assert SB._is_indel(_V("A", "ACGT")) and SB._is_indel(_V("AT", "A"))   # ins / del
    assert not SB._is_indel(_V("A", "G"))                                   # SNV (precomputed-complete)
    assert not SB._is_indel(_V("AT", "GC")) and not SB._is_indel(_V("A", "*"))  # MNV / spanning-del


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


def test_join_constraint():
    import csv as _csv
    import tempfile
    d = tempfile.mkdtemp()
    gn, sh, ph, out = (os.path.join(d, f) for f in ("gn.txt", "sh.tsv", "ph.tsv", "c.tsv"))
    with open(gn, "w") as f:
        f.write("gene\toe_lof_upper\tpLI\nTP53\t0.21\t0.99\nBRCA2\t0.55\t0.0\n")
    with open(sh, "w") as f:
        f.write("gene\tpost_mean\nTP53\t0.35\n")            # s_het via 'post_mean' alias
    with open(ph, "w") as f:
        f.write("#gene\tpHaplo\nTP53\t0.98\nNF1\t0.91\n")   # '#gene' header + gene not in gnomAD
    spec = importlib.util.spec_from_file_location(
        "jc", os.path.join(os.path.dirname(__file__), "..", "scripts", "join_constraint.py"))
    jc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(jc)
    assert jc.main(["--gnomad", gn, "--shet", sh, "--phaplo", ph, "--out", out]) == 0
    rows = {r["gene"]: r for r in _csv.DictReader(open(out), delimiter="\t")}
    import shutil
    shutil.rmtree(d, ignore_errors=True)
    assert rows["TP53"]["oe_lof_upper"] == "0.21" and rows["TP53"]["pli"] == "0.99"
    assert rows["TP53"]["s_het"] == "0.35" and rows["TP53"]["phaplo"] == "0.98"
    assert rows["NF1"]["phaplo"] == "0.91" and rows["NF1"]["oe_lof_upper"] == ""   # left-join keeps it
    assert rows["BRCA2"]["s_het"] == "" and rows["BRCA2"]["phaplo"] == ""


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


def _write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(h, "")) for h in header) + "\n")


def _nhf_tsv(path, entries):
    """entries: list of (variant_key, supporting_reads, nonhuman_fraction)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cols = ("variant_key", "supporting_reads", "nonhuman_fraction")
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for k, reads, frac in entries:
            fh.write(f"{k}\t{reads}\t{frac}\n")


def test_igv_nhf_join_is_pos_minus_one():
    """Step-8b NHF folds into variants.tsv on the 0-based key (pos-1), with the read
    denominator beside each fraction, and nhf_flag over the min_reads floor. A decoy key at
    the WRONG offset must not match — the off-by-one is the one load-bearing join bug."""
    d = tempfile.mkdtemp(prefix="_hprv_nhf_")
    data = os.path.join(d, "igv")
    os.makedirs(data, exist_ok=True)
    manifest = os.path.join(d, "trios.resolved.tsv")
    _write_tsv(manifest, ["trio_id", "vcf", "ped", "samples"],
               [{"trio_id": "T1", "vcf": "x.vcf.gz", "ped": "x.ped", "samples": "KID,DAD,MOM"}])

    calls = os.path.join(d, "candidates.calls.tsv")
    _write_tsv(
        calls, ["chrom", "pos", "ref", "alt", "trio_id", "mode",
                "child_gt", "mother_gt", "father_gt"],
        [
            {"chrom": "chr1", "pos": "100", "ref": "A", "alt": "T", "trio_id": "T1", "mode": "dominant"},
            {"chrom": "chr2", "pos": "200", "ref": "C", "alt": "CAT", "trio_id": "T1", "mode": "dominant"},
            {"chrom": "chr3", "pos": "300", "ref": "G", "alt": "*", "trio_id": "T1", "mode": "compound_het"},
            {"chrom": "chr4", "pos": "400", "ref": "A", "alt": "T", "trio_id": "T1", "mode": "dominant"},
        ],
    )

    # CHILD table (KID): key is 0-based, so pos-100 call -> chr1:99. A decoy chr1:100 (= 1-based
    # 101, a DIFFERENT variant) must NOT bleed onto the pos-100 row.
    _nhf_tsv(os.path.join(data, "nhf", "T1", "KID.variant_nhf.tsv"), [
        ("chr1:99:A:T", 10, "0.90"),      # <- correct match for the chr1:100 call
        ("chr1:100:A:T", 99, "0.10"),     # <- decoy at the wrong offset
        ("chr2:199:C:CAT", 8, "0.75"),    # indel, pos-200 -> 199
        ("chr4:399:A:T", 3, "0.80"),      # low read count (< min_reads) -> must not flag
    ])
    # MOTHER table (MOM): clean at the chr1 locus.
    _nhf_tsv(os.path.join(data, "nhf", "T1", "MOM.variant_nhf.tsv"), [
        ("chr1:99:A:T", 8, "0.00"),
    ])
    # FATHER (DAD): no table at all -> father columns blank.

    out = os.path.join(data, "variants.tsv")
    igv.build_variants_tsv(calls, manifest, data, out, nhf_dir=os.path.join(data, "nhf"),
                           nhf_min_reads=5)
    import csv as _csv
    with open(out) as fh:
        rows = {(r["chrom"], r["pos"]): r for r in _csv.DictReader(fh, delimiter="\t")}

    r1 = rows[("chr1", "100")]
    assert r1["child_nhf"] == "0.90", f"pos-1 join wrong: {r1['child_nhf']} (decoy leaked?)"
    assert r1["child_nhf_reads"] == "10", r1["child_nhf_reads"]
    assert r1["mother_nhf"] == "0.00" and r1["mother_nhf_reads"] == "8", r1
    assert r1["father_nhf"] == "" and r1["father_nhf_reads"] == "", "unscreened father must be blank"
    assert r1["nhf_flag"] == "1", f"0.90 over 10 reads should flag: {r1['nhf_flag']}"

    r2 = rows[("chr2", "200")]
    assert r2["child_nhf"] == "0.75" and r2["child_nhf_reads"] == "8", r2  # indel join

    r3 = rows[("chr3", "300")]  # symbolic '*' — nonhuman-screen skips it -> blank, but row emitted
    assert r3["child_nhf"] == "" and r3["nhf_flag"] == "", "symbolic ALT must be blank NHF"

    r4 = rows[("chr4", "400")]  # screened but only 3 reads (< min_reads 5)
    assert r4["child_nhf"] == "0.80", r4
    assert r4["nhf_flag"] == "0", f"below-floor read count must not flag (got {r4['nhf_flag']})"


def test_igv_nhf_disabled_is_blank():
    """No nhf_dir -> every NHF column blank and nhf_flag empty (legacy behavior preserved)."""
    d = tempfile.mkdtemp(prefix="_hprv_nhf0_")
    data = os.path.join(d, "igv"); os.makedirs(data, exist_ok=True)
    manifest = os.path.join(d, "m.tsv")
    _write_tsv(manifest, ["trio_id", "samples"], [{"trio_id": "T1", "samples": "KID,DAD,MOM"}])
    calls = os.path.join(d, "c.tsv")
    _write_tsv(calls, ["chrom", "pos", "ref", "alt", "trio_id"],
               [{"chrom": "chr1", "pos": "100", "ref": "A", "alt": "T", "trio_id": "T1"}])
    out = os.path.join(data, "variants.tsv")
    n = igv.build_variants_tsv(calls, manifest, data, out)   # nhf_dir defaults to None
    assert n == 1
    import csv as _csv
    with open(out) as fh:
        row = next(_csv.DictReader(fh, delimiter="\t"))
    for col in ("child_nhf", "child_nhf_reads", "mother_nhf", "father_nhf", "nhf_flag"):
        assert row[col] == "", f"{col} should be blank when NHF disabled, got {row[col]!r}"


# =============================================================================
# Step 9 — prioritization (src/hprv/prioritize.py + pipeline/09_prioritize.py)
# =============================================================================
class _Skip(Exception):
    """A test could not run because an optional dependency is absent. See `_requires`."""


def _requires(*mods):
    """Raise `_Skip` if any of `mods` cannot be imported.

    This file is contracted to run on a BARE HOST with no heavy dependencies (CLAUDE.md, Testing:
    "Host, no heavy deps ... python3 tests/test_pure.py"). The Step-9 tests that drive
    `09_prioritize.py:main()` break that contract transitively: main() reads the config through
    `hprv.config.load_config`, which does `import yaml`. Without a guard the whole suite dies with
    a ModuleNotFoundError partway through — which is what CI did.

    Skipping is the lesser evil, but ONLY if it is loud: `_run_all` refuses to print a clean bill
    of health when anything was skipped, because a partial run that reads as a full one is exactly
    the silent coverage loss this repo guards against everywhere else. CI installs pyyaml, so
    these tests DO execute there rather than being permanently skipped.
    """
    missing = []
    for m in mods:
        try:
            importlib.import_module(m)
        except ImportError:
            missing.append(m)
    if missing:
        raise _Skip(", ".join(missing))


def _load_p9():
    # The single chokepoint for every test that drives the Step-9 CLI, so the dependency those
    # tests need transitively is declared once, here, rather than in each of them.
    _requires("yaml")
    spec = importlib.util.spec_from_file_location(
        "p9", os.path.join(os.path.dirname(__file__), "..", "pipeline", "09_prioritize.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_prioritize_nb_tail_and_midp():
    """The NB tail must be the regularized incomplete beta, not a 1-cdf complement sum.

    The artifact tail reaches ~1e-30 on the extreme loci; a complement sum returns 0.0 for every
    one of them and collapses the whole tail into one indistinguishable bin. Also checks mid-p,
    which is the CALIBRATION statistic (a discrete upper tail is conservative at any threshold,
    so P(N>=n) cannot answer "is this null calibrated?").
    """
    from hprv import prioritize as PR
    # P(N >= 0) = 1 always; P(N >= 1) = 1 - pmf(0)
    assert PR.nb_sf(0, 5.0, 0.2) == 1.0
    p1 = PR.nb_sf(1, 5.0, 0.2)
    assert abs(p1 - (1.0 - math.exp(PR.nb_logpmf(0, 5.0, 0.2)))) < 1e-12
    # the tail must AGREE with an explicit pmf sum where the sum is still accurate
    for n in (2, 3, 8):
        direct = 1.0 - sum(math.exp(PR.nb_logpmf(k, 5.0, 0.2)) for k in range(n))
        assert abs(PR.nb_sf(n, 5.0, 0.2) - direct) < 1e-9, n
    # ...and stay accurate deep in the tail, where a 1-cdf complement is pure rounding noise.
    # OR4Q3's real numbers from the validation cohort: n=229 against E=0.42. The complement sum
    # returns ~1.1e-16 (the machine epsilon left over from summing to 1.0), which is wrong by
    # >200 orders of magnitude — so every extreme locus would land in one indistinguishable bin.
    deep = PR.nb_sf(229, 0.42, 0.2143)
    assert 0.0 < deep < 1e-200, deep
    complement = 1.0 - sum(math.exp(PR.nb_logpmf(k, 0.42, 0.2143)) for k in range(229))
    # The complement sum cannot REPRESENT a value this small: summing pmf terms to ~1.0 and
    # subtracting leaves either exactly 0.0 or a residue of order machine epsilon (~1e-16),
    # depending on the platform's libm summation order. Both are the same failure — every extreme
    # locus collapses into one indistinguishable bin — and both are asserted, because which one
    # you get is not a property of the statistic.
    assert complement == 0.0 or complement > 1e-30, \
        ("the complement sum must be catastrophically wrong here — that is why nb_sf uses "
         f"betainc (complement={complement:g} vs betainc={deep:g})")
    assert complement != deep
    # monotone non-increasing in n, non-decreasing in mu
    assert PR.nb_sf(3, 5.0, 0.2) > PR.nb_sf(4, 5.0, 0.2)
    assert PR.nb_sf(4, 5.0, 0.2) < PR.nb_sf(4, 9.0, 0.2)
    # NB variance exceeds Poisson at the same mean => a FATTER tail. This is the whole reason
    # the NB is used: the Poisson p-value at the same count is anti-conservative.
    assert PR.nb_sf(12, 5.0, 0.2) > PR.pois_sf(12, 5.0)
    # mid-p sits strictly inside the two one-sided discrete tails
    n, mu, a = 4, 5.0, 0.2
    assert PR.nb_sf(n + 1, mu, a) < PR.nb_midp(n, mu, a) < PR.nb_sf(n, mu, a)
    assert PR.pois_sf(n + 1, mu) < PR.pois_midp(n, mu) < PR.pois_sf(n, mu)
    # Poisson tail against an exact sum
    assert abs(PR.pois_sf(3, 2.0) - (1.0 - sum(math.exp(-2.0) * 2.0 ** k / math.factorial(k)
                                               for k in range(3)))) < 1e-12


def test_prioritize_bh_fdr():
    """BH must match the closed form, be monotone, and pass None through untouched."""
    from hprv import prioritize as PR
    ps = [0.001, 0.008, 0.039, 0.041, 0.9]
    q = PR.bh_fdr(ps)
    m = 5
    # step-down: q_i = min over j >= i of p_j * m / (j+1)
    expect = []
    running = 1.0
    for rank in range(m, 0, -1):
        running = min(running, ps[rank - 1] * m / rank)
        expect.append(running)
    expect.reverse()
    for got, want in zip(q, expect):
        assert abs(got - want) < 1e-12, (got, want)
    assert PR.bh_fdr([0.01, None, 0.5])[1] is None
    # q >= p and q is monotone in p-order
    assert all(a >= b - 1e-12 for a, b in zip(q, ps))
    assert q == sorted(q)
    # and it agrees with Step 6's implementation on the same vector — they must not drift
    gb = _load_gb()
    assert [round(x, 12) for x in PR.bh_fdr(ps)] == [round(x, 12) for x in gb.bh_fdr(ps)]


def test_prioritize_nb_fit_recovers_known_alpha_and_C():
    """The trimmed fit must recover the (C, alpha) that generated the data.

    Synthetic NB2 counts at known C and alpha, drawn via the gamma-Poisson mixture so the test
    needs no scipy. Recovery within ~15% on 4,000 genes is the honest tolerance for a
    method-of-moments/MLE dispersion at this sample size; the point is that the estimator is
    UNBIASED, not that it is exact.
    """
    import random
    from hprv import prioritize as PR
    random.seed(20260729)
    C_true, alpha_true, n_genes = 40000.0, 0.20, 4000
    mus, counts = [], []
    theta = 1.0 / alpha_true
    for _ in range(n_genes):
        mu = 10 ** random.uniform(-6.0, -4.0)
        mus.append(mu)
        lam = random.gammavariate(theta, (C_true * mu) / theta)   # NB = gamma-Poisson mixture
        # Poisson draw by inversion (lam is small here, so this terminates immediately)
        k, p, s = 0, math.exp(-lam), math.exp(-lam)
        u = random.random()
        while u > s and k < 1000:
            k += 1
            p *= lam / k
            s += p
        counts.append(k)
    fit = PR.fit_excess_null(counts, mus)
    assert abs(fit["C"] - C_true) / C_true < 0.10, fit["C"]
    assert abs(fit["alpha"] - alpha_true) / alpha_true < 0.35, fit["alpha"]
    assert abs(fit["theta"] - 1.0 / fit["alpha"]) < 1e-9
    # the bulk dispersion of correctly-modelled data must be ~1 (that is what phi_bulk measures)
    assert 0.5 < fit["phi_bulk"] < 2.0, fit["phi_bulk"]
    # method-of-moments should land in the same neighbourhood as the MLE
    mom = PR.fit_excess_null(counts, mus, alpha_method="mom")
    assert abs(mom["alpha"] - alpha_true) / alpha_true < 0.5, mom["alpha"]


def test_prioritize_C_must_be_fit_over_the_full_universe():
    """Fitting C on the matched (non-zero-count) genes only INFLATES C and DEFLATES every ratio.

    The candidate list is a zero-truncated sample. This is the documented direction of the error
    and it is the direction that HIDES artifact loci, so it must be asserted, not commented.
    """
    from hprv import prioritize as PR
    mus = [1e-5] * 100
    counts = [0] * 60 + [3] * 40                     # 60% of genes produced no candidate
    C_full = PR.fit_scaling(counts, mus)
    matched = [(c, m) for c, m in zip(counts, mus) if c]
    C_matched = PR.fit_scaling([c for c, _ in matched], [m for _, m in matched])
    assert C_matched > C_full, (C_matched, C_full)
    assert abs(C_matched / C_full - 100.0 / 40.0) < 1e-9   # exactly the zero-truncation factor
    # ...and therefore the excess ratio of a hit gene is SMALLER under the matched-only fit
    assert 3.0 / (C_matched * 1e-5) < 3.0 / (C_full * 1e-5)


def test_prioritize_offset_mu_lof_imputation_and_cds_fallback():
    """A null mu_lof is IMPUTED and labelled, never treated as 0; no-mu genes get the CDS
    regression with its provenance recorded."""
    from hprv import prioritize as PR
    mu, src = PR.mutational_target(1e-5, 4e-6, 7e-7)
    assert src == "gnomad" and abs(mu - (1e-5 + 4e-6 + 7e-7)) < 1e-18
    mu_i, src_i = PR.mutational_target(1e-5, 4e-6, None)
    assert src_i == "imputed"
    assert abs(mu_i - 1.4e-5 * 1.0516) < 1e-18
    assert mu_i > 1.4e-5, "a null mu_lof must not be charged as 0.0"
    assert PR.mutational_target(None, None, None) == (None, "none")
    # CDS fallback: monotone in length, and lands in the right order of magnitude
    m1, m2 = PR.mu_from_cds(1000), PR.mu_from_cds(10000)
    assert m1 < m2 and 1e-8 < m1 < 1e-3, (m1, m2)
    assert PR.mu_from_cds(None) is None and PR.mu_from_cds(0) is None


def test_prioritize_trim_loop_halts_on_mis_specified_null():
    """A trim loop that removes >trim_max_fraction of genes must HALT, not proceed.

    That is a mis-specified null (usually: the offset table and the counts describe different
    gene universes), not an exome that is 5% artifact — and proceeding would mass-down-weight
    real genes.
    """
    from hprv import prioritize as PR
    # every gene wildly over its target => the trim loop eats the exome
    mus = [1e-6] * 200
    counts = [50] * 200
    try:
        PR.fit_excess_null(counts, mus, trim_max_fraction=0.05)
    except ValueError as e:
        assert "trim_max_fraction" in str(e) and "mis-specified" in str(e)
    else:
        # A degenerate fit that absorbs everything into alpha is an acceptable alternative
        # outcome; what must NOT happen is a silent mass trim. Assert that directly.
        fit = PR.fit_excess_null(counts, mus, trim_max_fraction=1.0)
        assert fit["trim_fraction"] <= 0.05, \
            "a >5% trim must raise rather than be reported as a normal fit"


def test_prioritize_artifact_signals_and_corroboration():
    """Each of the six signals fires on its own value and the count is their unweighted sum."""
    from hprv import prioritize as PR
    thr = PR.signal_thresholds({}, caf_low_cutoff=1e-5)
    clean = {"gene": "TP53", "per_trio": 0.004, "segdup98_frac": 0.0, "oe_syn": 1.02,
             "constraint_flag": "", "classic_caf": 1e-3}
    s = PR.artifact_signals(clean, thr)
    assert s["corroboration_count"] == 0, s
    assert s["signal_reasons"] == []
    dirty = {"gene": "OR4Q3", "per_trio": 1.04, "segdup98_frac": 0.55, "oe_syn": 1.63,
             "constraint_flag": "mis_too_many", "classic_caf": 0.0}
    d = PR.artifact_signals(dirty, thr)
    assert d["corroboration_count"] == 6, d
    for k in ("sig_saturation", "sig_segdup", "sig_family", "sig_oe_syn",
              "sig_constraint_flag", "sig_caf_low"):
        assert d[k] is True, k
    # every flag travels with the VALUE that set it, so the reason string is self-documenting
    joined = "; ".join(d["signal_reasons"])
    for frag in ("cohort_saturation=", "segdup98_frac=", "artifact_gene_family",
                 "oe_syn=", "gnomad_constraint_flag=mis_too_many", "classic_caf="):
        assert frag in joined, (frag, joined)
    # oe_syn departure is SYMMETRIC — under-calling (gnomAD missed variants, so mu is
    # unreliable) and over-calling (paralogue collapse) are both informative
    assert PR.artifact_signals(dict(clean, oe_syn=0.55), thr)["sig_oe_syn"]
    assert PR.artifact_signals(dict(clean, oe_syn=1.45), thr)["sig_oe_syn"]
    assert not PR.artifact_signals(dict(clean, oe_syn=1.25), thr)["sig_oe_syn"]
    # `no_exp_lof` is deliberately NOT a model-failure flag
    assert not PR.artifact_signals(dict(clean, constraint_flag="no_exp_lof"), thr)["sig_constraint_flag"]
    # a MISSING classic_caf reads as 0 (per spec) and lands inside the bottom-decile cut;
    # a missing per_trio (no --n-trios) must NOT fire saturation
    assert PR.artifact_signals(dict(clean, classic_caf=None), thr)["sig_caf_low"]
    assert not PR.artifact_signals(dict(clean, per_trio=None), thr)["sig_saturation"]
    # family membership is a REGEX list from config, and it is never sufficient alone
    assert PR.artifact_signals(dict(clean, gene="KRTAP9-1"), thr)["corroboration_count"] == 1
    assert PR.family_matcher(["^ZZZ"])("KRTAP9-1") == ""


def test_prioritize_gene_tier_rules():
    """The four tiers, each rule shape, and the never-drop contract on the penalty scale."""
    from hprv import prioritize as PR
    T = PR.assign_gene_tier
    # T0: nothing fires
    t = T(1.1, 1.0, 0, False, 1)
    assert t["gene_tier"] == "T0_no_downweight" and t["gene_artifact_penalty"] == 0.0
    # T1 needs NO corroboration (it costs -0.5 and populates a watch list)
    assert T(7.7, 0.2, 0, False, 8)["gene_tier"] == "T1_watch"
    assert T(7.7, 0.9, 0, False, 8)["gene_tier"] == "T0_no_downweight"     # q too high
    # T2 uses the RATIO plus one signal — an FDR-only rule would exempt every low-count gene
    assert T(6.4, 0.2, 1, False, 6)["gene_tier"] == "T2_downweight"
    assert T(6.4, 0.2, 0, False, 6)["gene_tier"] == "T1_watch"              # no corroboration
    # T3, all three disjuncts
    assert T(4.0, 0.01, 1, False, 20)["gene_tier"] == "T3_strong_downweight"   # q & corrob>=1
    assert T(24.8, 0.2, 2, False, 20)["gene_tier"] == "T3_strong_downweight"   # ratio>=10 & 2
    assert T(6.0, 0.2, 2, True, 20)["gene_tier"] == "T3_strong_downweight"     # saturation & 2
    assert T(6.0, 0.2, 2, False, 20)["gene_tier"] == "T2_downweight"           # no saturation
    # penalties are on the Tavtigian points scale, and the MAXIMUM penalty cannot alone demote a
    # strong variant: V4 splice (+4) + strong rarity (+2) + constraint (+1) = +7, -3 leaves +4.
    pen = T(24.8, 0.001, 3, True, 200)["gene_artifact_penalty"]
    assert pen == -3.0
    assert 4.0 + 2.0 + 1.0 + pen == 4.0, "a T3 penalty must be a re-rank, not a veto"
    # every tier keeps its reason machinery — nothing is ever dropped
    for tier in PR.GENE_TIERS:
        assert tier in ("T0_no_downweight", "T1_watch", "T2_downweight", "T3_strong_downweight")


def test_prioritize_positive_control_guard():
    """SENSITIVITY GUARD: an established predisposition gene is NEVER placed in T2/T3.

    On the validation cohort 13 control genes met a T2/T3 rule and were capped by this ceiling —
    CTSA at 148x excess and NPRL3 at 77x among them. They are established recessive-disease genes
    AND technically problematic loci; the exemption protects them from a score PENALTY, and the
    review_flag says their calls still need read-level review. If this regresses, real
    predisposition genes get their variants penalised silently, which is the one failure this
    whole scheme is built to prevent.
    """
    from hprv import prioritize as PR
    # CTSA's real numbers from the validation cohort: n=160, 148x, q~0, saturating, corroborated
    for corrob in range(0, 7):
        for sat in (False, True):
            t = PR.assign_gene_tier(148.35, 1e-40, corrob, sat, 160,
                                    established_gene_control=True)
            assert t["gene_tier"] in ("T0_no_downweight", "T1_watch"), (corrob, sat, t)
            assert t["gene_artifact_penalty"] >= -0.5, t
    # the ceiling must be VISIBLE and the gene must still be flagged for review
    t = PR.assign_gene_tier(148.35, 1e-40, 3, True, 160, established_gene_control=True)
    assert t["control_ceiling_applied"] is True
    assert t["review_flag"] == "established_gene_high_excess", t
    # ...and the identical gene WITHOUT the control flag does reach T3 (so the ceiling, not a
    # weak rule, is what protected it)
    t_no = PR.assign_gene_tier(148.35, 1e-40, 3, True, 160, established_gene_control=False)
    assert t_no["gene_tier"] == "T3_strong_downweight" and t_no["control_ceiling_applied"] is False
    # a control gene that reaches T1 ON ITS OWN MERITS is not an exemption — do not conflate them
    own = PR.assign_gene_tier(14.9, 0.2, 0, False, 21, established_gene_control=True)
    assert own["gene_tier"] == "T1_watch" and own["control_ceiling_applied"] is False


def test_prioritize_cds_fallback_ceiling():
    """A gene whose offset came from the CDS regression never reaches T3: a +/-30% offset error
    cannot support that claim."""
    from hprv import prioritize as PR
    args = dict(established_gene_control=False)
    strong = PR.assign_gene_tier(50.0, 1e-20, 3, True, 40, e_source="gnomad_mu", **args)
    assert strong["gene_tier"] == "T3_strong_downweight"
    capped = PR.assign_gene_tier(50.0, 1e-20, 3, True, 40, e_source="cds_fallback", **args)
    assert capped["gene_tier"] == "T2_downweight" and capped["cds_ceiling_applied"] is True
    # and the reason string SAYS the offset was a fallback, so the ratio is never read as exact
    reason = PR.downweight_reason(
        {"excess_ratio": 50.0, "E_expected": 0.8, "q_nb": 1e-20, "n_observed": 40,
         "E_source": "cds_fallback"}, ["cohort_saturation=0.5_variants_per_trio"], capped)
    assert "cds_fallback_ceiling_applied" in reason and "cds_length_fallback" in reason


def test_prioritize_unexplained_excess_is_flagged_not_penalised():
    """A high excess with NO mechanism identified is FLAGGED for review and carries no penalty.

    29 genes / 148 variants on the validation cohort, 2 of them established controls. With no
    mechanism and the never-drop rule in force, the correct action is to look, not to down-weight.
    """
    from hprv import prioritize as PR
    t = PR.assign_gene_tier(23.7, 0.2, 0, False, 14)
    assert t["review_flag"] == "unexplained_excess"
    assert t["gene_tier"] == "T1_watch" and t["gene_artifact_penalty"] == -0.5
    assert "no_corroborating_signature" in PR.downweight_reason(
        {"excess_ratio": 23.7, "E_expected": 0.59, "q_nb": 0.2, "n_observed": 14,
         "E_source": "gnomad_mu"}, [], t)


def test_prioritize_variant_tier_and_mechanism_gating():
    """The V0-V5 ladder, and the SVI rule that a benign prediction caps the tier regardless of
    how constrained the gene is."""
    from hprv import prioritize as PR
    V = lambda **kw: PR.assign_variant_tier(kw)["variant_tier"]          # noqa: E731
    # V4: strong splice, or an NMD-indeterminate pLoF (= every pLoF today)
    assert V(spliceai_ds="0.55", impact="MODIFIER", consequence="intron_variant") == "V4"
    assert V(impact="HIGH", consequence="stop_gained") == "V4"
    assert V(impact="HIGH", consequence="splice_donor_variant") == "V4"
    assert V(impact="HIGH", consequence="frameshift_variant") == "V4"
    # V5 is UNREACHABLE — no variants.tsv column supports the NMD-escape test
    info = PR.assign_variant_tier({"impact": "HIGH", "consequence": "stop_gained"})
    assert info["nmd_status"] == "INDETERMINATE" and "V5_unreachable" in info["variant_tier_reason"]
    assert info["plof_confidence"] == "UNAVAILABLE"          # no LOFTEE at any price
    # V3: supporting splice, or a high-CADD missense (a DISCOVERY RANK, labelled as such)
    assert V(spliceai_ds="0.25", impact="MODIFIER", consequence="intron_variant") == "V3"
    mis = PR.assign_variant_tier({"impact": "MODERATE", "consequence": "missense_variant",
                                  "cadd": "30"})
    assert mis["variant_tier"] == "V3"
    assert mis["missense_evidence_source"] == "cadd_offlabel", \
        "a CADD-based missense tier must never be presentable as calibrated PP3 evidence"
    # V2: missense below the cut, and in-frame indels
    assert V(impact="MODERATE", consequence="missense_variant", cadd="10") == "V2"
    assert V(impact="MODERATE", consequence="inframe_deletion", ref="ATCG", alt="A") == "V2"
    # V1: non-coding / synonymous discovery rank
    assert V(impact="MODIFIER", consequence="intron_variant", cadd="26") == "V1"
    # V0 requires BOTH scores PRESENT and low
    assert V(impact="MODIFIER", consequence="intron_variant", cadd="3", spliceai_ds="0.01") == "V0"
    # --- MECHANISM GATING: a benign prediction zeroes constraint AND the gene-list prior ---
    assert PR.constraint_gate("V0", "dominant") == 0.0
    assert PR.constraint_gate("V1", "dominant") == 0.5
    assert PR.constraint_gate("V2", "dominant") == 0.5
    assert PR.constraint_gate("V4", "dominant") == 1.0
    # ...and constraint is zeroed for every recessive mode: pLoF constraint measures selection
    # against HETEROZYGOTES, so it is not evidence about a biallelic candidate (either direction)
    for mode in ("compound_het", "hom_recessive", "homozygous", "x_linked_recessive"):
        assert PR.constraint_gate("V4", mode) == 0.0, mode


def test_prioritize_mechanism_gating_caps_benign_in_constrained_gene():
    """A molecularly-benign variant in a highly constrained, list-member gene is CAPPED.

    The single most important structural rule: no amount of gene-level enthusiasm — constraint,
    recurrence, or gene-list membership — can rescue a molecularly-benign prediction. Without
    this, a gene-list prior becomes confirmation bias.
    """
    from hprv import prioritize as PR
    benign = {"impact": "MODIFIER", "consequence": "intron_variant", "cadd": "2",
              "spliceai_ds": "0.01", "grpmax_af": "1e-7", "inheritance": "dominant",
              "child_GQ": "99", "child_DP": "40", "child_AB": "0.5"}
    gene = {"pLI": "0.999", "oe_lof_upper": "0.05", "gene_tier": "T0_no_downweight",
            "n_carriers": "5", "recurrence_kind": "distinct_variant"}
    sc = PR.score_variant(benign, gene, {}, gene_prior=True)
    assert sc["variant_tier"] == "V0"
    assert sc["pts_gene_constraint"] == 0.0, "constraint must not rescue a V0"
    assert sc["pts_gene_list_prior"] == 0.0, "a gene-list prior must not rescue a V0"
    assert sc["cap_applied"] == "V0_benign"
    assert sc["priority_points_agnostic"] <= 0.0 and sc["priority_points_prior"] <= 0.0
    # the SAME gene with a credible molecular effect DOES earn the constraint term
    strong = dict(benign, impact="HIGH", consequence="stop_gained", cadd="38", spliceai_ds="")
    sc2 = PR.score_variant(strong, gene, {})
    assert sc2["variant_tier"] == "V4" and sc2["pts_gene_constraint"] == 1.0
    assert sc2["cap_applied"] == "none" and sc2["priority_points_agnostic"] > 5.0
    # BA1 caps too, and neither cap removes the row
    ba1 = dict(strong, grpmax_af="0.2")
    sc3 = PR.score_variant(ba1, gene, {})
    assert sc3["rarity_strength"] == "fail" and "BA1_frequency" in sc3["cap_applied"]
    assert sc3["priority_points_agnostic"] <= -4.0


def test_prioritize_blank_vs_zero_nhf():
    """BLANK NHF IS NOT 0.0 — three states, never two, with a decoy blank that must not read
    as clean.

    blank = NOT SCREENED (not an ALT carrier, no mini-CRAM, or Step 8b never ran);
    0.0   = SCREENED and every read classified human.
    Treating blank as clean silently promotes exactly the calls nobody examined. Mirrors
    Analysis/.../inherited/prepare_igv_variants.py:nhf_status.
    """
    from hprv import prioritize as PR
    # DECOY: every member blank -> not_screened, and it must score 0 (no penalty, NO CREDIT)
    blank = {"child_nhf": "", "child_nhf_reads": "", "mother_nhf": "", "father_nhf": ""}
    assert PR.nhf_state(blank) == ("not_screened", None, None)
    # screened and clean -> a REAL 0.0, distinguishable from the blank above
    clean = {"child_nhf": "0.0", "child_nhf_reads": "30"}
    state, frac, reads = PR.nhf_state(clean)
    assert (state, frac, reads) == ("clean", 0.0, 30.0)
    # flagged: over threshold AND over the read floor
    assert PR.nhf_state({"child_nhf": "0.9", "child_nhf_reads": "10"})[0] == "flagged"
    # the min_reads floor is essential: NHF 1.0 over 2 reads is noise, not evidence
    assert PR.nhf_state({"child_nhf": "1.0", "child_nhf_reads": "2"})[0] == "clean"
    # ...and it must be reported as screened-with-a-value, not silently hidden
    assert PR.nhf_state({"child_nhf": "1.0", "child_nhf_reads": "2"})[1] == 1.0
    # any screened member counts, not just the child
    assert PR.nhf_state({"child_nhf": "", "mother_nhf": "0.8", "mother_nhf_reads": "9"})[0] == "flagged"
    # --- and the SCORE must distinguish them: not_screened scores 0, flagged is penalised ---
    base = {"impact": "HIGH", "consequence": "stop_gained", "grpmax_af": "1e-6",
            "inheritance": "dominant", "child_GQ": "99", "child_DP": "40", "child_AB": "0.5"}
    s_blank = PR.score_variant(dict(base, **blank), {}, {})
    s_clean = PR.score_variant(dict(base, **clean), {}, {})
    s_flag = PR.score_variant(dict(base, child_nhf="0.9", child_nhf_reads="10"), {}, {})
    assert s_blank["nhf_status"] == "not_screened" and s_clean["nhf_status"] == "clean"
    assert s_blank["pts_quality"] == 0.0 and s_clean["pts_quality"] == 0.0
    assert s_flag["pts_quality"] == -3.0
    assert s_blank["priority_points_agnostic"] == s_clean["priority_points_agnostic"], \
        "an unscreened call must score neither better nor worse than a clean one — it is an " \
        "explicit uncertainty flag, and the DIFFERENCE is carried by nhf_status, not the points"


def test_prioritize_absent_scores_are_not_benign_evidence():
    """Absence of evidence must never read as evidence of absence, on any axis."""
    from hprv import prioritize as PR
    # a missing SpliceAI score is "not covered by the precomputed set", not "no splice effect",
    # so it can never satisfy the V0 benign rule
    no_sai = {"impact": "MODIFIER", "consequence": "intron_variant", "cadd": "2"}
    assert PR.assign_variant_tier(no_sai)["variant_tier"] != "V0"
    assert PR.spliceai_status(None) == "not_covered"
    assert PR.spliceai_status("0.0") == "scored"          # a measured zero IS evidence
    # a missing CADD likewise
    assert PR.assign_variant_tier({"impact": "MODIFIER", "consequence": "intron_variant",
                                   "spliceai_ds": "0.01"})["variant_tier"] != "V0"
    # an absent AF reads as unknown (treated as rarest), never as a measured common frequency
    assert PR.rarity_strength(None) == "unknown"
    assert PR.rarity_strength("") == "unknown"
    assert PR.rarity_strength(3e-6) == "strong"
    assert PR.rarity_strength(5e-5) == "moderate"
    assert PR.rarity_strength(5e-4) == "supporting"
    assert PR.rarity_strength(5e-3) == "permissive"
    assert PR.rarity_strength(0.2) == "fail"
    # absent genotype metrics are "not assessed", not a QC failure
    ok, why = PR.genotype_qc({})
    assert ok is True and "absent" in why
    assert PR.genotype_qc({"child_GQ": "5", "child_DP": "40", "child_AB": "0.5"})[0] is False
    assert PR.genotype_qc({"child_GQ": "99", "child_DP": "40", "child_AB": "0.95"})[0] is False
    assert PR.genotype_qc({"child_GQ": "99", "child_DP": "40", "child_AB": "0.95",
                           "child_gt": "1/1"})[0] is True          # hom-alt band, not het
    # a missing MOI curation is EXACTLY neutral — the novel-gene case
    assert PR.moi_coherence("dominant", "") == ("unknown", "")
    assert PR.score_variant({"impact": "HIGH", "consequence": "stop_gained",
                             "inheritance": "dominant"}, {}, {})["pts_moi"] == 0.0


def test_prioritize_moi_coherence_and_long_gene_caveat():
    """MOI discordance is a small demotion + a flag, never a filter; and audit A-6's long-gene
    comp-het drift suppresses the penalty entirely."""
    from hprv import prioritize as PR
    assert PR.moi_coherence("dominant", "AD")[0] == "coherent"
    assert PR.moi_coherence("dominant", "AR")[0] == "discordant"
    assert PR.moi_coherence("hom_recessive", "AR")[0] == "coherent"
    assert PR.moi_coherence("compound_het", "AD")[0] == "discordant"
    assert PR.moi_coherence("dominant", "AD;AR")[0] == "coherent"     # multi-MOI genes
    # audit A-6: mode assignment is single-gene-keyed and long genes drift to compound_het, so a
    # discordance penalty there would punish an artifact of the mode assignment
    coh, caveat = PR.moi_coherence("compound_het", "AD", long_gene=True)
    assert coh == "discordant" and caveat == "long_gene_comphet_drift"
    row = {"impact": "HIGH", "consequence": "stop_gained", "inheritance": "compound_het"}
    long_gene = {"gene_moi": "AD", "cds_length": "100000"}
    short_gene = {"gene_moi": "AD", "cds_length": "1200"}
    assert PR.score_variant(row, long_gene, {})["pts_moi"] == 0.0
    assert PR.score_variant(row, short_gene, {})["pts_moi"] == -1.0


def test_prioritize_recurrence_scores_on_count_not_pvalue():
    """Recurrence must score on the CARRIER COUNT, capped, and same-variant gets less credit.

    Step 6's recurrence null is a case-only approximation, so p saturates: for essentially any
    gene with >= min_carriers carriers of private variants it clears the exome-wide line
    (audit A-2). A p-value that saturates cannot order anything.
    """
    from hprv import prioritize as PR
    row = {"impact": "HIGH", "consequence": "stop_gained", "inheritance": "dominant"}
    def rec(**g):
        return PR.score_variant(row, g, {})["pts_recurrence"]
    assert rec(n_carriers="1") == 0.0
    assert rec(n_carriers="2", recurrence_kind="distinct_variant") == 1.0
    assert rec(n_carriers="3", recurrence_kind="distinct_variant") == 2.0
    assert rec(n_carriers="40", recurrence_kind="distinct_variant") == 2.0, "capped at +2"
    # same-variant recurrence is as easily a mapping artifact or founder allele as a burden signal
    assert rec(n_carriers="40", recurrence_kind="same_variant") == 0.5
    # a tiny p_recurrence must NOT move the term
    assert rec(n_carriers="2", recurrence_kind="distinct_variant", p_recurrence="1e-30") == 1.0


def test_prioritize_clinvar_has_no_star_gate():
    """P/LP is honored unstarred (the cache has no CLNREVSTAT) and the column SAYS so."""
    from hprv import prioritize as PR
    assert PR.clinvar_strength("pathogenic") == "p_lp"
    assert PR.clinvar_strength("Pathogenic/Likely_pathogenic") == "p_lp"
    assert PR.clinvar_strength("conflicting_classifications_of_pathogenicity") == "conflicting"
    assert PR.clinvar_strength("benign") == "benign"
    assert PR.clinvar_strength("uncertain_significance") == "vus"
    assert PR.clinvar_strength("") == "absent"
    sc = PR.score_variant({"impact": "HIGH", "consequence": "stop_gained",
                           "clin_sig": "pathogenic"}, {}, {})
    assert sc["pts_clinical"] == 4.0
    assert sc["clinvar_review_status"] == "UNAVAILABLE", \
        "the +4 must never be mistakable for a >=2-star assertion"


def test_prioritize_every_term_is_a_separate_column():
    """Transparency is a hard requirement: a reviewer must be able to read WHY a variant ranked
    where it did, and the terms must sum to the reported total."""
    from hprv import prioritize as PR
    row = {"impact": "HIGH", "consequence": "splice_donor_variant", "spliceai_ds": "0.92",
           "grpmax_af": "2e-6", "inheritance": "dominant", "clin_sig": "pathogenic",
           "child_GQ": "99", "child_DP": "40", "child_AB": "0.5",
           "child_nhf": "0.0", "child_nhf_reads": "30"}
    gene = {"pLI": "0.99", "gene_tier": "T1_watch", "gene_artifact_penalty": "-0.5",
            "n_carriers": "3", "recurrence_kind": "distinct_variant", "gene_moi": "AD"}
    sc = PR.score_variant(row, gene, {})
    terms = ["pts_molecular", "pts_rarity", "pts_gene_constraint", "pts_recurrence",
             "pts_quality", "pts_clinical", "pts_moi", "pts_gene_artifact"]
    for t in terms + ["pts_gene_list_prior", "priority_points_agnostic",
                      "priority_points_prior", "cap_applied", "variant_tier_reason"]:
        assert t in sc, t
    assert abs(sum(sc[t] for t in terms) - sc["priority_points_agnostic"]) < 1e-12
    assert abs(sc["priority_points_prior"]
               - (sc["priority_points_agnostic"] + sc["pts_gene_list_prior"])) < 1e-12
    # and the total is what the spec's worked example predicts: +4 splice, +2 rarity,
    # +1 constraint, +2 recurrence, +4 ClinVar, -0.5 artifact
    assert sc["priority_points_agnostic"] == 4 + 2 + 1 + 2 + 4 - 0.5, sc["priority_points_agnostic"]


def test_prioritize_gene_prior_off_by_default_and_gated():
    """With no overlay the two rankings are IDENTICAL — that identity is what keeps hprv
    phenotype-agnostic, so it is asserted rather than assumed."""
    from hprv import prioritize as PR
    rows = []
    for i, (imp, cq, af) in enumerate((("HIGH", "stop_gained", "1e-6"),
                                       ("MODERATE", "missense_variant", "5e-5"),
                                       ("MODIFIER", "intron_variant", "5e-4"))):
        r = {"chrom": "chr1", "pos": str(1000 + i), "ref": "A", "alt": "T", "impact": imp,
             "consequence": cq, "grpmax_af": af, "inheritance": "dominant", "cadd": "26"}
        r.update(PR.score_variant(r, {}, {}, gene_prior=False))
        rows.append(r)
    PR.rank_rows(rows, "priority_points_agnostic", "rank_agnostic")
    PR.rank_rows(rows, "priority_points_prior", "rank_prior")
    assert all(r["rank_prior"] == r["rank_agnostic"] for r in rows)
    assert all(r["pts_gene_list_prior"] == 0.0 for r in rows)
    # with the overlay ON, the prior enters ONE term and only the prior total moves
    on = PR.score_variant(rows[0], {}, {}, gene_prior=True)
    assert on["pts_gene_list_prior"] == 2.0
    assert on["priority_points_agnostic"] == rows[0]["priority_points_agnostic"], \
        "the agnostic ranking must stay free of any gene-list influence"
    assert on["priority_points_prior"] == on["priority_points_agnostic"] + 2.0
    # ...and it is mechanism-gated: half credit at V2, none at V0
    v2 = PR.score_variant({"impact": "MODERATE", "consequence": "missense_variant", "cadd": "5"},
                          {}, {}, gene_prior=True)
    assert v2["variant_tier"] == "V2" and v2["pts_gene_list_prior"] == 1.0
    # ranks are deterministic (ties broken to a fixed order), so a re-run reproduces them
    again = [dict(r) for r in rows]
    PR.rank_rows(again, "priority_points_agnostic", "rank_agnostic")
    assert [r["rank_agnostic"] for r in again] == [r["rank_agnostic"] for r in rows]


def test_prioritize_gene_prior_overlay_accepts_csv_and_rejects_headerless_table():
    """A phenotype panel handed over as a CSV must APPLY, and a headerless table must HALT.

    The overlay reader's delimiter gate was tab-only, so a spreadsheet-exported CSV fell through
    to the bare-symbol-list branch: every line became one "symbol" ("BRCA1,0.9,GREEN") that can
    never match a gene, while the reader announced "3 genes from pheno.csv" and the run proceeded
    fully phenotype-agnostic. A plausible-looking gene count over a list that matched nothing is
    the failure class nobody catches in review — the reviewer ships an un-prioritised list
    believing their panel was applied. So: CSV works, and a table this reader cannot recognise is
    a HARD STOP rather than the usual degrade-with-a-WARN, because its absence is not honestly
    reportable the way a missing --mutrate is.
    """
    p9 = _load_p9()
    d = tempfile.mkdtemp(prefix="_hprv_p9ovl_")
    try:
        cols = ["chrom", "pos", "ref", "alt", "trio_id", "gene", "consequence", "impact",
                "inheritance", "grpmax_af", "cadd", "child_gt", "child_GQ", "child_DP", "child_AB"]

        def row(pos, gene):
            return {"chrom": "chr1", "pos": str(pos), "ref": "A", "alt": "T", "trio_id": "T1",
                    "gene": gene, "consequence": "stop_gained", "impact": "HIGH",
                    "inheritance": "dominant", "grpmax_af": "2e-6", "cadd": "38",
                    "child_gt": "0/1", "child_GQ": "99", "child_DP": "40", "child_AB": "0.5"}
        vin = os.path.join(d, "v.tsv")
        _write_tsv(vin, cols, [row(1000, "BRCA1"), row(2000, "TP53"), row(3000, "SDHB")])
        cfgp = os.path.join(d, "cfg.yaml")
        with open(cfgp, "w") as fh:
            fh.write("project: {name: t}\nprioritization:\n"
                     "  gene_downweight: {min_control_genes: 1}\n"
                     "  composite: {gene_list_prior: {enabled: true}}\n")
        outv, outg = os.path.join(d, "vp.tsv"), os.path.join(d, "gp.tsv")

        def run(overlay):
            for f in (outv + ".done",):
                if os.path.exists(f):
                    os.remove(f)
            return p9.main(["--variants", vin, "--config", cfgp, "--gene-prior", overlay,
                            "--n-trios", "1", "--out-variants", outv, "--out-genes", outg])

        def members():
            import csv
            return {r["gene"]: r for r in csv.DictReader(open(outv), delimiter="\t")
                    if r["gene_list_prior_member"] == "1"}

        # (1) a spreadsheet export: COMMA-delimited, with a header
        csvp = os.path.join(d, "panel.csv")
        with open(csvp, "w") as fh:
            fh.write("gene,prior_weight,tier\nBRCA1,0.9,GREEN\nTP53,0.5,AMBER\n")
        assert run(csvp) == 0
        m = members()
        assert set(m) == {"BRCA1", "TP53"}, f"CSV overlay did not apply: {sorted(m)}"
        assert m["BRCA1"]["gene_list_prior_tier"] == "GREEN"
        assert m["BRCA1"]["gene_list_prior_weight"] == "0.9"

        # (2) the same content with NO header row is unreadable -> hard stop, not a silent no-op
        bad = os.path.join(d, "headerless.csv")
        with open(bad, "w") as fh:
            fh.write("BRCA1,0.9,GREEN\nTP53,0.5,AMBER\n")
        assert run(bad) == 1, "a headerless table must fail loudly, not match zero genes quietly"

        # (3) a genuine bare symbol list still works
        plain = os.path.join(d, "panel.txt")
        with open(plain, "w") as fh:
            fh.write("# my phenotype panel\nBRCA1\nSDHB\n")
        assert run(plain) == 0
        assert set(members()) == {"BRCA1", "SDHB"}

        # (4) the JSON sidecar's gene_sets block, combined by MAX and never SUM
        tsvp = os.path.join(d, "panel2.tsv")
        with open(tsvp, "w") as fh:
            fh.write("gene\tprior_weight\ttier\nBRCA1\t0.9\tGREEN\n")
        with open(os.path.join(d, "panel2.json"), "w") as fh:
            fh.write('{"gene_sets": {"MYPATH": {"prior_weight": 0.6, '
                     '"members": ["BRCA1", "TP53", "SDHB"]}}}')
        assert run(tsvp) == 0
        m = members()
        assert set(m) == {"BRCA1", "TP53", "SDHB"}, f"sidecar set not applied: {sorted(m)}"
        # BRCA1 is in BOTH at 0.9 and 0.6 -> MAX, so its weight must not become 1.5
        assert float(m["BRCA1"]["gene_list_prior_weight"]) == 0.9, \
            "gene-level and set-level priors must combine by MAX, never SUM"
        assert float(m["TP53"]["gene_list_prior_weight"]) == 0.6
    finally:
        import shutil
        shutil.rmtree(d, ignore_errors=True)


def test_prioritize_igv_review_table_preserves_track_paths():
    """--out-igv-variants must ADD triage columns without disturbing ANY input column.

    The review table is what a reviewer opens in igv.js, and `variants.prioritized.tsv` cannot
    serve that role: it lives outside igv/ and its column set omits the *_file/*_index/*_vcf*
    track paths, which are RELATIVE to the igv/ data dir. Dropping them yields a sortable list
    with no mini-CRAMs and no VCF tracks — the read-level view Step 8 exists to provide. So this
    asserts the three things that make it a drop-in: every input column present, in its original
    order, byte-identical; the prioritization columns appended; and row-count conserved.

    Also pins the phenotype-overlay columns into the review table, since a --gene-prior run is
    useless for review if the reviewer cannot see WHICH prior fired and what it was worth.
    """
    p9 = _load_p9()
    d = tempfile.mkdtemp(prefix="_hprv_p9igv_")
    try:
        # A Step-8-shaped table: annotations, then the relative track paths that must survive.
        cols = ["chrom", "pos", "ref", "alt", "trio_id", "gene", "consequence", "impact",
                "inheritance", "grpmax_af", "cadd", "spliceai_ds", "clin_sig", "child_gt",
                "child_GQ", "child_DP", "child_AB", "nhf_flag",
                "child_file", "child_index", "child_vcf", "child_vcf_index", "child_vcf_id"]

        def row(pos, gene, **kw):
            r = {"chrom": "chr1", "pos": str(pos), "ref": "A", "alt": "T", "trio_id": "T1",
                 "gene": gene, "consequence": "stop_gained", "impact": "HIGH",
                 "inheritance": "dominant", "grpmax_af": "2e-6", "cadd": "38",
                 "spliceai_ds": "0.9", "clin_sig": "", "child_gt": "0/1", "child_GQ": "99",
                 "child_DP": "40", "child_AB": "0.5", "nhf_flag": "0",
                 "child_file": "crams/T1/kid.cram", "child_index": "crams/T1/kid.cram.crai",
                 "child_vcf": "vcfs/T1.vcf.gz", "child_vcf_index": "vcfs/T1.vcf.gz.tbi",
                 "child_vcf_id": "kid"}
            r.update(kw)
            return r
        # Two ALTs of one multiallelic site in one trio: a chrom/pos/ref/alt/trio_id join is
        # ambiguous here, which is why the merge carries rows by position instead.
        rows = [row(1000, "PHENO1"), row(1000, "PHENO1", alt="G"),
                row(2000, "OTHER", impact="MODIFIER", consequence="intron_variant",
                    cadd="2", spliceai_ds="0.01")]
        vin = os.path.join(d, "variants.tsv")
        _write_tsv(vin, cols, rows)
        pheno = os.path.join(d, "pheno.tsv")
        with open(pheno, "w") as fh:
            fh.write("gene\tprior_weight\ttier\tevidence_class\n"
                     "PHENO1\t0.9\tGREEN\tgermline_predisposition\n")
        cfgp = os.path.join(d, "cfg.yaml")
        with open(cfgp, "w") as fh:
            fh.write("project: {name: t}\nprioritization:\n"
                     "  gene_downweight: {min_control_genes: 1}\n"
                     "  composite: {gene_list_prior: {enabled: true}}\n")
        outv, outg = os.path.join(d, "vp.tsv"), os.path.join(d, "gp.tsv")
        outi = os.path.join(d, "variants.prioritized.tsv")
        rc = p9.main(["--variants", vin, "--config", cfgp, "--gene-prior", pheno,
                      "--n-trios", "1", "--out-variants", outv, "--out-genes", outg,
                      "--out-igv-variants", outi])
        assert rc == 0, f"Step 9 exited {rc}"

        import csv
        with open(outi) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            got = list(csv.DictReader(fh, fieldnames=header, delimiter="\t"))
        # (1) input columns kept, in order, at the FRONT — the merge only ever appends
        assert header[:len(cols)] == cols, f"input columns disturbed: {header[:len(cols)]}"
        # (2) the triage columns a reviewer sorts on, including the overlay provenance
        for c in ("variant_tier", "priority_points_agnostic", "rank_agnostic", "rank_prior",
                  "rank_delta", "gene_tier", "downweight_reason", "review_flag", "nhf_status",
                  "gene_list_prior_member", "gene_list_prior_tier", "gene_list_prior_weight",
                  "pts_gene_list_prior"):
            assert c in header, f"review table is missing {c}"
        assert len(header) == len(set(header)), "duplicate column in the review table"
        # (3) never-drop, and every input value byte-identical
        assert len(got) == len(rows), f"{len(rows)} in, {len(got)} out"
        by_alt = {}
        for r in got:
            by_alt.setdefault((r["pos"], r["alt"]), r)
        for src in rows:
            m = by_alt[(src["pos"], src["alt"])]
            for c in cols:
                assert m[c] == src[c], f"input column {c} altered: {m[c]!r} != {src[c]!r}"
        # the overlay is visible per-variant, and only for its own gene
        assert by_alt[("1000", "T")]["gene_list_prior_member"] == "1"
        assert by_alt[("1000", "T")]["gene_list_prior_tier"] == "GREEN"
        assert by_alt[("2000", "T")]["gene_list_prior_member"] == "0"
        # sorted by rank_agnostic so the file opens already prioritized
        ranks = [int(r["rank_agnostic"]) for r in got]
        assert ranks == sorted(ranks), f"review table not sorted by rank_agnostic: {ranks}"

        # adding --out-igv-variants to an ALREADY-CACHED run must still produce the file
        os.remove(outi)
        assert p9.main(["--variants", vin, "--config", cfgp, "--gene-prior", pheno,
                        "--n-trios", "1", "--out-variants", outv, "--out-genes", outg,
                        "--out-igv-variants", outi]) == 0
        assert os.path.exists(outi) and os.path.getsize(outi) > 0, \
            "cache hit skipped a requested output"
    finally:
        import shutil
        shutil.rmtree(d, ignore_errors=True)


def test_clinvar_stars_parses_every_rendering_and_absent_is_not_zero():
    """CLNREVSTAT arrives in several shapes, and absent must never read as 0 stars.

    The VCF value carries commas INSIDE it ("criteria_provided,_multiple_submitters,_no_conflicts")
    while the field is Number=., so the comma is also the array separator: cyvcf2 may hand back a
    string or a tuple and bcftools query joins with commas. A plain string compare against any one
    of those forms breaks on the others, so the parser canonicalises. Verified against a real
    bcftools transfer + cyvcf2 read before this test was written.

    The load-bearing distinction: None (the ClinVar transfer did not run — nobody looked) vs 0
    ("submitter provided no assertion criteria"). Collapsing them would let a run with no ClinVar
    resource damp every P/LP assertion as if it were unreviewed.
    """
    from hprv import annotations as AN
    canon = AN._canon_revstat
    stars = lambda v: AN._REVSTAT_STARS.get(canon(v))      # noqa: E731
    assert stars("criteria_provided,_multiple_submitters,_no_conflicts") == 2
    assert stars(("criteria_provided", "_multiple_submitters", "_no_conflicts")) == 2
    assert stars("criteria provided, multiple submitters, no conflicts") == 2
    assert stars("practice_guideline") == 4
    assert stars("reviewed_by_expert_panel") == 3
    assert stars("criteria_provided,_single_submitter") == 1
    # ClinVar renamed this in 2024; a pinned older release is still a legitimate input
    assert stars("criteria_provided,_conflicting_classifications") == 1
    assert stars("criteria_provided,_conflicting_interpretations") == 1
    assert stars("no_assertion_criteria_provided") == 0
    # absent / unknown -> None, NEVER 0
    assert stars(None) is None and stars("") is None
    assert stars("a_status_clinvar_invents_later") is None, \
        "an unrecognised status is an unknown, not a zero-star assertion"


def test_prioritize_tier_is_monotone_in_splice_evidence():
    """Adding splice evidence must never LOWER a variant's tier, and must never erase the
    missense evidence source.

    The splice rungs return early. Once the missense branch gained a V4 rung (REVEL >= Pejaver
    moderate), that early return became a SCORING INVERSION: a missense with revel=0.9 AND a
    supporting spliceai_ds=0.25 scored V3 with missense_evidence_source="none", while the SAME
    variant with ds=0.15 scored V4/revel. Strictly more evidence, strictly worse tier, and a
    source column asserting no predictor spoke when REVEL had said 0.9. Reproduced against the
    real function before the fix.

    The two signals are independent mechanisms (amino-acid effect vs splice disruption), so the
    combination takes the STRONGER tier and reports BOTH reasons.
    """
    from hprv import prioritize as PR
    base = {"consequence": "missense_variant", "impact": "MODERATE", "ref": "A", "alt": "T"}
    idx = {t: i for i, t in enumerate(PR.VARIANT_TIERS)}

    # the exact inversion, pinned
    lo = PR.assign_variant_tier({**base, "revel": "0.9", "spliceai_ds": "0.15"}, {})
    hi = PR.assign_variant_tier({**base, "revel": "0.9", "spliceai_ds": "0.25"}, {})
    assert idx[hi["variant_tier"]] >= idx[lo["variant_tier"]], \
        f"more splice evidence lowered the tier: {lo['variant_tier']} -> {hi['variant_tier']}"
    assert hi["missense_evidence_source"] == "revel", \
        f"a splice rung erased the missense source: {hi['missense_evidence_source']!r}"
    assert "revel" in hi["variant_tier_reason"] and "spliceai" in hi["variant_tier_reason"], \
        f"both signals must be reported: {hi['variant_tier_reason']!r}"

    # monotone across the whole grid, for every REVEL band
    for rev in ("", "0.1", "0.5", "0.7", "0.9"):
        best = -1
        for ds in ("", "0.05", "0.15", "0.25", "0.4", "0.6", "0.95"):
            t = PR.assign_variant_tier({**base, "revel": rev, "spliceai_ds": ds}, {})["variant_tier"]
            assert idx[t] >= best, f"revel={rev!r} ds={ds!r}: tier fell to {t}"
            best = max(best, idx[t])

    # a splice signal still wins on a BENIGN-REVEL missense (different mechanism), and says so
    r = PR.assign_variant_tier({**base, "revel": "0.05", "spliceai_ds": "0.6"}, {})
    assert r["variant_tier"] == "V4" and "revel=0.05" in r["variant_tier_reason"]

    # non-missense paths are untouched by the merge
    for cq, im in (("stop_gained", "HIGH"), ("intron_variant", "MODIFIER")):
        r = PR.assign_variant_tier({"consequence": cq, "impact": im, "ref": "A", "alt": "T",
                                    "spliceai_ds": "0.25", "cadd": "30"}, {})
        assert r["missense_evidence_source"] == "none"


def test_prioritize_missense_predictor_precedence():
    """REVEL -> AlphaMissense -> CADD(off-label) -> none, in that fixed order.

    NOT a max over whatever is available: ClinGen SVI's rule is to commit to ONE predictor chosen
    before seeing results, so taking the best of N would be an uncalibrated cherry-pick. The tier
    must always report which predictor spoke, and a variant scored by REVEL must not be re-scored
    by CADD just because CADD happens to be higher.
    """
    from hprv import prioritize as PR
    base = {"consequence": "missense_variant", "impact": "MODERATE", "ref": "A", "alt": "T"}

    def tier(**kw):
        return PR.assign_variant_tier({**base, **kw}, {})

    # REVEL wins even when CADD would give a different answer
    t = tier(revel="0.9", cadd="1.0")
    assert t["missense_evidence_source"] == "revel" and t["variant_tier"] == "V4", t
    t = tier(revel="0.70", cadd="40")
    assert t["missense_evidence_source"] == "revel" and t["variant_tier"] == "V3", t
    t = tier(revel="0.1", cadd="40")
    assert t["missense_evidence_source"] == "revel" and t["variant_tier"] == "V1", \
        f"a calibrated BENIGN REVEL must not be overridden by a high CADD: {t}"
    t = tier(revel="0.5", cadd="40")
    assert t["missense_evidence_source"] == "revel" and t["variant_tier"] == "V2", t

    # AlphaMissense only when REVEL is absent (both are missense-only; neither covers everything)
    t = tier(alphamissense="0.9", cadd="1.0")
    assert t["missense_evidence_source"] == "alphamissense" and t["variant_tier"] == "V3", t
    t = tier(alphamissense="0.1", cadd="40")
    assert t["missense_evidence_source"] == "alphamissense" and t["variant_tier"] == "V1", t

    # CADD last, and explicitly labelled off-label so no reader mistakes it for calibrated
    t = tier(cadd="40")
    assert t["missense_evidence_source"] == "cadd_offlabel" and t["variant_tier"] == "V3", t
    t = tier(cadd="2")
    assert t["missense_evidence_source"] == "none" and t["variant_tier"] == "V2", t
    # and a blank string is absent, not zero
    t = tier(revel="", alphamissense="", cadd="40")
    assert t["missense_evidence_source"] == "cadd_offlabel", t


def test_prioritize_clinvar_star_gate_damps_only_the_positive_limb():
    """Low stars damp a P/LP assertion; they must NOT damp a benign one, and absent = full weight.

    Scaling a negative (benign) term toward zero would PROMOTE a poorly-reviewed benign call —
    the opposite of the intent. And an absent star count means the ClinVar transfer did not run,
    which must leave the term at full weight rather than damping every assertion in the run.
    """
    from hprv import prioritize as PR
    cfg = {"resources": {"clinvar": {"min_review_stars": 2, "low_star_scale": 0.5}}}
    row = {"consequence": "missense_variant", "impact": "MODERATE", "ref": "A", "alt": "T",
           "inheritance": "dominant", "grpmax_af": "1e-6", "child_gt": "0/1"}

    def pts(clin_sig, stars):
        r = dict(row, clin_sig=clin_sig)
        if stars is not None:
            r["clinvar_stars"] = str(stars)
        return PR.score_variant(r, {}, cfg)

    full = pts("Pathogenic", 3)["pts_clinical"]
    assert full > 0
    assert pts("Pathogenic", 1)["pts_clinical"] == full * 0.5, "1-star P/LP must be damped"
    assert pts("Pathogenic", None)["pts_clinical"] == full, \
        "no ClinVar transfer must leave the clinical term at FULL weight, not damped"
    assert pts("Pathogenic", None)["clinvar_review_status"] == "UNAVAILABLE"
    assert pts("Pathogenic", 3)["clinvar_review_status"] == "3_star"
    ben_hi = pts("Benign", 3)["pts_clinical"]
    ben_lo = pts("Benign", 1)["pts_clinical"]
    assert ben_hi < 0 and ben_lo == ben_hi, \
        f"a low-star BENIGN call must not be shrunk toward zero (would promote it): {ben_lo} vs {ben_hi}"


def test_prioritize_cache_invalidates_on_config_and_resource_change():
    """Step 9's idempotency marker must be CONTENT-keyed over inputs + config + resources.

    Every threshold in this step comes from the config, and whether the excess statistic exists
    at all depends on which optional resources were supplied. An mtime test against --variants
    alone reported "cached" after a threshold change (reproduced: the tier table was unchanged
    because the stale output was served), which is the same silent-staleness class the
    Step-1/2/4/8b caches are content-keyed to avoid. Also guards the legacy zero-byte marker
    written by the earlier touch-only code: it must NOT be read as a cache hit.
    """
    p9 = _load_p9()
    d = tempfile.mkdtemp(prefix="_hprv_p9cache_")
    try:
        mut = os.path.join(d, "mut.tsv")
        _write_tsv(mut, ["gene", "mu_mis", "mu_syn", "mu_lof", "cds_length"],
                   [{"gene": "BG000", "mu_mis": "8e-6", "mu_syn": "3e-6", "mu_lof": "5e-7",
                     "cds_length": "1500"}])
        cols = ["chrom", "pos", "ref", "alt", "trio_id", "gene", "consequence", "impact",
                "inheritance", "grpmax_af", "cadd", "child_gt"]
        _write_tsv(os.path.join(d, "v.tsv"), cols,
                   [{"chrom": "chr1", "pos": "1000", "ref": "A", "alt": "T", "trio_id": "T1",
                     "gene": "BG000", "consequence": "missense_variant", "impact": "MODERATE",
                     "inheritance": "dominant", "grpmax_af": "2e-6", "cadd": "12",
                     "child_gt": "0/1"}])
        vin = os.path.join(d, "v.tsv")
        cfgp = os.path.join(d, "cfg.yaml")
        outv, outg = os.path.join(d, "vp.tsv"), os.path.join(d, "gp.tsv")
        marker = outv + ".done"

        def run(*extra):
            return p9.main(["--variants", vin, "--config", cfgp,
                            "--out-variants", outv, "--out-genes", outg] + list(extra))

        def write_cfg(min_n):
            with open(cfgp, "w") as fh:
                fh.write("project: {name: t}\nprioritization:\n"
                         "  gene_downweight: {min_control_genes: 1}\n"
                         f"  excess: {{min_n_for_ratio_rule: {min_n}}}\n")

        write_cfg(3)
        assert run("--force") == 0
        k1 = open(marker).read().strip()
        assert len(k1) == 64, f"marker must hold a sha256 content key, got {k1!r}"

        # unchanged -> the key is STABLE across runs (a key that never matches is useless)
        assert run() == 0
        assert open(marker).read().strip() == k1

        # config threshold changed -> key must change
        write_cfg(1)
        assert run() == 0
        k2 = open(marker).read().strip()
        assert k2 != k1, "a config threshold change must invalidate the cache"

        # a newly-supplied resource must also invalidate (it turns the excess statistic ON)
        assert run("--mutrate", mut) == 0
        k3 = open(marker).read().strip()
        assert k3 != k2, "supplying --mutrate must invalidate the cache"

        # a legacy zero-byte marker (old touch-only code) must not read as a hit
        open(marker, "w").close()
        assert run("--mutrate", mut) == 0
        assert len(open(marker).read().strip()) == 64
    finally:
        import shutil
        shutil.rmtree(d, ignore_errors=True)


def test_prioritize_never_drop_row_conservation():
    """THE NEVER-DROP INVARIANT, end to end through the CLI step: output rows == input rows.

    Every path is exercised: a BA1-common variant, a benign-prediction variant, a QC-failing
    call, an NHF-flagged call, a gene in the strongest down-weight tier, and a gene with no
    mutational-target row at all. None of them may vanish.
    """
    import shutil
    p9 = _load_p9()
    d = tempfile.mkdtemp(prefix="_hprv_p9_")
    try:
        mut = os.path.join(d, "mut.tsv")
        _write_tsv(mut, ["gene", "mu_mis", "mu_syn", "mu_lof", "oe_syn", "classic_caf",
                         "constraint_flag", "cds_length", "pli", "oe_lof_upper", "segdup98_frac"],
                   # OR4Q3 is the artifact locus (tiny target, huge count, every signal);
                   # TP53 is the clean constrained gene; NOGENE is absent from this table.
                   [{"gene": "OR4Q3", "mu_mis": "3e-6", "mu_syn": "1e-6", "mu_lof": "1e-7",
                     "oe_syn": "1.63", "classic_caf": "0", "constraint_flag": "mis_too_many",
                     "cds_length": "900", "pli": "0.01", "oe_lof_upper": "1.9",
                     "segdup98_frac": "0.55"},
                    {"gene": "TP53", "mu_mis": "9e-6", "mu_syn": "4e-6", "mu_lof": "8e-7",
                     "oe_syn": "1.01", "classic_caf": "1e-3", "constraint_flag": "",
                     "cds_length": "1182", "pli": "0.99", "oe_lof_upper": "0.2",
                     "segdup98_frac": "0.0"}]
                   + [{"gene": f"BG{i:03d}", "mu_mis": "8e-6", "mu_syn": "3e-6", "mu_lof": "5e-7",
                       "oe_syn": "1.0", "classic_caf": "2e-4", "constraint_flag": "",
                       "cds_length": "1500", "pli": "0.1", "oe_lof_upper": "1.1",
                       "segdup98_frac": "0.0"} for i in range(300)])
        cols = ["chrom", "pos", "ref", "alt", "trio_id", "gene", "consequence", "impact",
                "inheritance", "grpmax_af", "max_af", "cadd", "spliceai_ds", "clin_sig",
                "child_gt", "child_GQ", "child_DP", "child_AB",
                "child_nhf", "child_nhf_reads", "mother_nhf", "father_nhf"]
        rows = []

        decoys = {}

        def add(gene, tag=None, **kw):
            r = {"chrom": "chr1", "pos": str(1000 + len(rows) * 11), "ref": "A", "alt": "T",
                 "trio_id": f"T{len(rows) % 20:02d}", "gene": gene,
                 "consequence": "missense_variant", "impact": "MODERATE",
                 "inheritance": "dominant", "grpmax_af": "2e-6", "max_af": "2e-6",
                 "cadd": "12", "spliceai_ds": "0.02", "clin_sig": "",
                 "child_gt": "0/1", "child_GQ": "99", "child_DP": "40", "child_AB": "0.5",
                 "child_nhf": "", "child_nhf_reads": "", "mother_nhf": "", "father_nhf": ""}
            r.update(kw)
            rows.append(r)
            if tag:
                decoys[tag] = (r["chrom"], r["pos"])

        for i in range(300):
            add(f"BG{i:03d}", tag=("blank_nhf" if i == 5 else None))
        for _ in range(45):                                  # the artifact pileup
            add("OR4Q3")
        add("TP53", tag="control", consequence="stop_gained", impact="HIGH", cadd="38",
            spliceai_ds="")
        add("BG001", tag="ba1", grpmax_af="0.2", max_af="0.2")     # capped, NOT dropped
        add("BG002", tag="v0", impact="MODIFIER", consequence="intron_variant", cadd="2",
            spliceai_ds="0.01")                                    # capped, NOT dropped
        add("BG003", tag="qc_fail", child_GQ="5")                  # penalised, NOT dropped
        add("BG004", tag="nhf_flagged", child_nhf="0.95", child_nhf_reads="20")
        add("NOGENE", tag="no_offset")                             # no mutational-target row
        vin = os.path.join(d, "variants.tsv")
        _write_tsv(vin, cols, rows)
        cfgp = os.path.join(d, "cfg.yaml")
        with open(cfgp, "w") as fh:
            fh.write("project: {name: t}\nprioritization:\n"
                     "  gene_downweight: {min_control_genes: 1}\n")
        ctrl = os.path.join(d, "ctrl.txt")
        with open(ctrl, "w") as fh:
            fh.write("TP53\n")
        outv = os.path.join(d, "vp.tsv")
        outg = os.path.join(d, "gp.tsv")
        rc = p9.main(["--variants", vin, "--mutrate", mut, "--constraint", mut,
                      "--established-genes", ctrl, "--config", cfgp, "--n-trios", "20",
                      "--out-variants", outv, "--out-genes", outg])
        assert rc == 0, rc
        import csv as _csv
        with open(outv) as fh:
            got = list(_csv.DictReader(fh, delimiter="\t"))
        # THE INVARIANT
        assert len(got) == len(rows), f"never-drop VIOLATED: {len(rows)} in, {len(got)} out"
        by_gene = {}
        for r in got:
            by_gene.setdefault(r["gene"], []).append(r)
        # every input gene survives, including the one with no offset row
        for g in ("OR4Q3", "TP53", "NOGENE", "BG001", "BG002", "BG003", "BG004"):
            assert g in by_gene, g
        # the artifact locus IS down-weighted, and its reason string names the mechanisms
        art = by_gene["OR4Q3"][0]
        assert art["gene_tier"] in ("T2_downweight", "T3_strong_downweight"), art["gene_tier"]
        assert float(art["pts_gene_artifact"]) < 0.0
        for frag in ("excess_ratio=", "cohort_saturation=", "artifact_gene_family"):
            assert frag in art["downweight_reason"], (frag, art["downweight_reason"])
        # the established gene is never down-weighted past the ceiling
        assert by_gene["TP53"][0]["gene_tier"] in ("T0_no_downweight", "T1_watch")
        # a gene with NO offset gets an explicit 'none' source, not a fabricated ratio
        assert by_gene["NOGENE"][0]["E_source"] == "none"
        assert by_gene["NOGENE"][0]["excess_ratio"] == ""
        # the caps and penalties fired without removing anything. Keyed by POSITION, not gene: a
        # decoy shares its gene with a background row, and asserting on the gene would silently
        # test the wrong row.
        pos_key = {(r["chrom"], r["pos"]): r for r in got}
        assert pos_key[decoys["ba1"]]["cap_applied"].endswith("BA1_frequency")
        assert pos_key[decoys["v0"]]["cap_applied"] == "V0_benign"
        assert float(pos_key[decoys["qc_fail"]]["pts_quality"]) <= -2.0
        assert pos_key[decoys["nhf_flagged"]]["nhf_status"] == "flagged"
        assert float(pos_key[decoys["nhf_flagged"]]["pts_quality"]) <= -3.0
        # ...while a decoy BLANK-NHF row reads as not_screened, never as clean
        assert pos_key[decoys["blank_nhf"]]["nhf_status"] == "not_screened"
        assert float(pos_key[decoys["blank_nhf"]]["pts_quality"]) == 0.0
        # ranks are a permutation of 1..N in both rankings, and identical with no overlay
        assert sorted(int(r["rank_agnostic"]) for r in got) == list(range(1, len(got) + 1))
        assert all(r["rank_prior"] == r["rank_agnostic"] for r in got)
        assert all(r["rank_delta"] == "0" for r in got)
        # the gene table covers every gene that produced a call
        with open(outg) as fh:
            grows = list(_csv.DictReader(fh, delimiter="\t"))
        assert {r["gene"] for r in grows} == {r["gene"] for r in rows}
        # idempotent: a second run is a no-op (the .done marker), and --force re-runs
        assert p9.main(["--variants", vin, "--mutrate", mut, "--config", cfgp,
                        "--out-variants", outv, "--out-genes", outg]) == 0
        assert os.path.exists(outv + ".done")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_prioritize_min_count_guard_on_ratio_rules():
    """A RATIO rule needs a minimum count behind it — the spec's own `min_n_for_ratio_rule`.

    `excess_statistic_spec` §4.4 warns that "a single variant against E=0.1 is 10x excess and
    means nothing" and specifies the guard; the T2 rule as originally written never applied it.
    Independent validation measured the cost: 136 of 228 triaged genes (59.6%) had n<3 and carried
    only 196 variants (7.6% of triaged volume), and 8 of the 13 ceiling-exempted control genes had
    q_nb = 1 — no statistical evidence of excess at all. Default is now 3.
    """
    from hprv import prioritize as PR
    T = PR.assign_gene_tier
    # a single variant at 6x excess must NOT reach T2 under the default guard
    assert T(6.3, 1.0, 2, False, 1)["gene_tier"] == "T0_no_downweight"
    assert T(6.3, 1.0, 2, False, 2)["gene_tier"] == "T0_no_downweight"
    assert T(6.3, 1.0, 2, False, 3)["gene_tier"] == "T2_downweight", "n>=3 clears the floor"
    # the RATIO limbs of T3 are guarded too
    assert T(24.8, 1.0, 2, False, 2)["gene_tier"] == "T0_no_downweight"
    assert T(24.8, 1.0, 2, False, 3)["gene_tier"] == "T3_strong_downweight"
    assert T(6.0, 1.0, 2, True, 2)["gene_tier"] == "T0_no_downweight"       # saturation limb
    assert T(6.0, 1.0, 2, True, 3)["gene_tier"] == "T3_strong_downweight"
    # ...but the FDR limb is NOT guarded and needs no guard: no gene with n<5 reached q<0.05 on
    # the validation cohort, so the NB tail already handles low counts on its own. Guarding it
    # would suppress a genuinely significant low-count gene for no measured benefit.
    assert T(4.0, 0.01, 1, False, 1)["gene_tier"] == "T3_strong_downweight"
    # setting the floor to 1 reproduces the ORIGINAL spec tier table exactly (both are documented)
    off = {"prioritization": {"excess": {"min_n_for_ratio_rule": 1}}}
    assert T(6.3, 1.0, 2, False, 1, cfg=off)["gene_tier"] == "T2_downweight"
    # the four real single-variant exemptions the guard removes: SMPX/HMGA2/HAMP/PET100 were each
    # ONE variant against E~0.16 at q_nb=1, ceiling-exempted under the old default
    for ratio in (6.3, 6.1, 6.0, 6.0):
        t = T(ratio, 1.0, 2, False, 1, established_gene_control=True)
        assert t["gene_tier"] == "T0_no_downweight" and not t["control_ceiling_applied"], \
            "a 1-variant control gene needs the COUNT guard, not the ceiling"
    # ...while CTSA/NPRL3/CDH23 still require the ceiling (n=160/108/23, q<0.05)
    for ratio, q, n in ((148.35, 1e-107, 160), (77.22, 2.7e-60, 108), (7.7, 3.9e-4, 23)):
        t = T(ratio, q, 1, False, n, established_gene_control=True)
        assert t["control_ceiling_applied"] is True, (ratio, n)
        assert t["gene_tier"] == "T1_watch"


def test_prioritize_poisson_arm_gets_its_own_trim():
    """The trim loop belongs to WHICHEVER null is being fit — a fair NB-vs-Poisson comparison.

    Evaluating the Poisson at the NB's trimmed C measures a hybrid nobody would deploy. The
    published 2.51x Poisson anti-conservatism is reproducible only with the Poisson arm trimmed
    (1.82x without), so both arms must be fit the way they would be used.
    """
    from hprv import prioritize as PR
    import random
    random.seed(20260729)
    mus, counts = [], []
    for i in range(3000):
        mu = 10 ** random.uniform(-6.0, -4.0)
        mus.append(mu)
        lam = random.gammavariate(5.0, (40000.0 * mu) / 5.0)
        k, p, s, u = 0, math.exp(-lam), math.exp(-lam), random.random()
        while u > s and k < 500:
            k += 1
            p *= lam / k
            s += p
        counts.append(k)
    # a hard artifact tail, so the two fits have something to disagree about
    for _ in range(20):
        mus.append(1e-7)
        counts.append(40)
    pf = PR.fit_poisson_null(counts, mus)
    untrimmed = PR.fit_scaling(counts, mus)
    # the trimmed Poisson C is PULLED DOWN by removing the tail — that is the whole point
    assert pf["C"] < untrimmed, (pf["C"], untrimmed)
    assert pf["n_trimmed"] >= 20 and pf["iterations"] >= 2
    assert 0.0 < pf["trim_fraction"] <= 0.05
    # calibrate_null must USE it, and report which C the Poisson arm was judged at
    nb = PR.fit_excess_null(counts, mus)
    cal = PR.calibrate_null(counts, mus, nb["C"], nb["alpha"])
    assert "C_poisson" in cal and abs(cal["C_poisson"] - pf["C"]) < 1e-6, \
        "the Poisson arm must be evaluated at its OWN trimmed C, not the NB's"
    assert cal["C_poisson"] != nb["C"]
    # an explicitly supplied C_poisson wins (so a caller can reproduce a published number)
    cal2 = PR.calibrate_null(counts, mus, nb["C"], nb["alpha"], C_poisson=untrimmed)
    assert abs(cal2["C_poisson"] - untrimmed) < 1e-9
    # Under the CORRECT reading (each arm at its own trimmed C) the Poisson tail is the more
    # liberal of the two — the direction that decides the design. Asserted on the trimmed
    # comparison only, because the two readings are genuinely different measurements: that is the
    # whole finding, and on a synthetic whose untrimmed C happens to be inflated the Poisson can
    # even read conservative. Which is exactly why the number must never be quoted without saying
    # which reading produced it.
    assert cal["poisson_frac_lt_0.001"] >= cal["nb_frac_lt_0.001"], cal
    assert cal["poisson_mean_midp"] > cal["nb_mean_midp"], cal
    assert cal["poisson_mean_midp"] != cal2["poisson_mean_midp"], \
        "trimming the Poisson arm must CHANGE its calibration — otherwise the fix is inert"
    # A Poisson trim that blows the trim_max_fraction guard must not take the NB diagnostic down
    # with it: fall back to the untrimmed Poisson C, SAY SO, and still report the NB numbers.
    # 10% of genes at 300x their target — a mis-specified Poisson, which is the realistic shape
    # (the Poisson trim removes far more genes than the NB's, because the NB absorbs overdispersion
    # into alpha rather than trimming it away).
    random.seed(1)
    bmus = [1e-5] * 900 + [1e-7] * 100
    bcounts = [random.choice([0, 1, 2]) for _ in range(900)] + [30] * 100
    try:
        PR.fit_poisson_null(bcounts, bmus)
        raise AssertionError("a >5% Poisson trim must raise, not be reported as a normal fit")
    except ValueError as e:
        assert "trim_max_fraction" in str(e) and "mis-specified" in str(e)
    bad = PR.calibrate_null(bcounts, bmus,
                            PR.fit_excess_null(bcounts, bmus, trim_max_fraction=0.5)["C"], 0.3)
    assert bad.get("poisson_trim_failed") is True, bad
    assert "nb_mean_midp" in bad and bad["n_evaluated"] > 0, \
        "a failed Poisson arm must not suppress the NB diagnostic"


def test_prioritize_gene_prior_overlay_is_moi_agnostic():
    """MODEL-MISMATCH HAZARD: an overlay prior must be expressible independently of canonical MOI.

    The concrete case. The FA/HR genes (FANCA, FANCD2, SLX4, FANCE, BRCA2) carry germ-cell-tumour
    evidence about HETEROZYGOUS carriers (PMID 40906985: five-gene combined OR 10.17, 95% CI
    4.87-21.27, P=2.90e-06), while their canonical Mendelian model — and their PanelApp green
    status — is biallelic Fanconi anemia. If the overlay join were routed by canonical MOI, or if
    an observed-mode-vs-canonical mismatch were penalised, those rows would be consulted under a
    recessive model and silently miss the very hypothesis they exist to encode.
    """
    from hprv import prioritize as PR
    # a canonically-RECESSIVE gene with a het-carrier prior, observed as a dominant het
    overlay = PR.parse_gene_prior_overlay([
        {"gene": "FANCA", "tier": "T2", "prior_weight": "0.6",
         "evidence_class": "rare_variant_association_single_study_significant",
         "moi": "AR_biallelic_FA;heterozygous_carrier_risk_proposed",
         "gene_sets": "FA_HR_PATHWAY_23", "pmids": "40906985;34906449"},
    ])
    assert "FANCA" in overlay, "the join must be by SYMBOL, never routed by MOI"
    row = {"impact": "HIGH", "consequence": "stop_gained", "grpmax_af": "2e-6",
           "inheritance": "dominant", "child_gt": "0/1", "child_GQ": "99", "child_DP": "40",
           "child_AB": "0.5"}
    # canonical MOI says AR; the observation is a HET. The prior must STILL apply...
    gene = {"gene": "FANCA", "gene_moi": "AR", "cds_length": "4400"}
    sc = PR.score_variant(row, gene, {}, gene_prior=overlay["FANCA"])
    assert sc["pts_gene_list_prior"] > 0.0, \
        "a het-carrier prior must apply to a het observation in a canonically-recessive gene"
    # ...and the mismatch must be a FLAG, never a penalty
    assert sc["moi_coherence"] == "discordant"
    assert sc["moi_caveat"] == "moi_mismatch_het_in_recessive_gene"
    assert sc["pts_moi"] == 0.0, "the mismatch is informational — never a down-weight by itself"
    # with the overlay's own richer MOI vocabulary, the carrier annotation makes it COHERENT
    sc2 = PR.score_variant(row, dict(gene, gene_moi="AR_biallelic_FA;heterozygous_carrier_risk_proposed"),
                           {}, gene_prior=overlay["FANCA"])
    assert sc2["moi_coherence"] == "coherent" and sc2["pts_moi"] == 0.0
    assert sc2["pts_gene_list_prior"] == sc["pts_gene_list_prior"]
    # a genuine dominant/recessive discordance with NO carrier hypothesis still reports the flag
    # and still charges nothing (same never-drop logic on the coherence layer)
    assert PR.moi_coherence("dominant", "AR")[1] == "moi_mismatch_het_in_recessive_gene"
    # the reverse direction (a biallelic call in a curated-AD gene) keeps its ordinary penalty
    rec = dict(row, inheritance="hom_recessive")
    sc3 = PR.score_variant(rec, {"gene_moi": "AD", "cds_length": "1200"}, {})
    assert sc3["moi_coherence"] == "discordant" and sc3["pts_moi"] == -1.0
    # per-gene weights scale the term: T1 (1.0) outranks T2 (0.6) outranks T3 (0.35)
    tiered = PR.parse_gene_prior_overlay([
        {"gene": "CHEK2", "tier": "T1", "prior_weight": "1.0", "evidence_class": "replicated"},
        {"gene": "GENET3", "tier": "T3", "prior_weight": "0.35",
         "evidence_class": "gwas_common_variant_locus"}])
    p1 = PR.score_variant(row, {}, {}, gene_prior=tiered["CHEK2"])["pts_gene_list_prior"]
    p3 = PR.score_variant(row, {}, {}, gene_prior=tiered["GENET3"])["pts_gene_list_prior"]
    assert p1 > p3 > 0.0, (p1, p3)
    # the V0 mechanism gate still zeroes the prior — a list must never rescue a benign prediction
    v0 = {"impact": "MODIFIER", "consequence": "intron_variant", "cadd": "2",
          "spliceai_ds": "0.01", "grpmax_af": "1e-7", "inheritance": "dominant"}
    assert PR.score_variant(v0, {"pLI": "0.999"}, {},
                            gene_prior=tiered["CHEK2"])["pts_gene_list_prior"] == 0.0


def test_prioritize_overlay_set_prior_is_max_not_sum():
    """DOUBLE-COUNTING HAZARD: a gene-level and a set-level prior are combined by MAX, never SUM.

    The GCT resource carries both per-gene FA rows AND a pathway-collapsed FA_HR_PATHWAY_23 entry
    (23-gene set, pooled OR 4.14, 95% CI 1.98-8.66, P=.0013). Both derive from PMID 40906985, so
    summing them double-counts one study. The resource's own usage contract says so; this test is
    the code enforcing it.
    """
    from hprv import prioritize as PR
    gene_sets = {"FA_HR_PATHWAY_23": {"prior_weight": 0.6, "n_members": 3,
                                      "members": ["FANCA", "FANCE", "FANCB"]}}
    overlay = PR.parse_gene_prior_overlay([
        {"gene": "FANCA", "tier": "T2", "prior_weight": "0.6", "gene_sets": "FA_HR_PATHWAY_23"},
        {"gene": "FANCE", "tier": "T4", "prior_weight": "0.15", "gene_sets": "FA_HR_PATHWAY_23"},
    ], gene_sets=gene_sets)
    # FANCE: gene 0.15, set 0.60 -> MAX 0.60, and the SUM (0.75) must never appear
    assert abs(overlay["FANCE"]["points_scale"] - 0.60) < 1e-9, overlay["FANCE"]
    assert overlay["FANCE"]["points_scale"] < 0.15 + 0.60
    assert overlay["FANCE"]["set_id_applied"] == "FA_HR_PATHWAY_23"
    # FANCA: gene 0.60 already >= set 0.60, so nothing is added
    assert abs(overlay["FANCA"]["points_scale"] - 0.60) < 1e-9
    # a set member with no gene row of its own is admitted at the SET weight
    assert overlay["FANCB"]["source"] == "gene_set"
    assert abs(overlay["FANCB"]["points_scale"] - 0.60) < 1e-9
    # ...and the scored points obey the same bound
    row = {"impact": "HIGH", "consequence": "stop_gained", "grpmax_af": "2e-6",
           "inheritance": "dominant"}
    maxpts = PR.default_weights()["gene_list_prior"]
    for g in ("FANCA", "FANCE", "FANCB"):
        pts = PR.score_variant(row, {}, {}, gene_prior=overlay[g])["pts_gene_list_prior"]
        assert pts <= maxpts + 1e-9, (g, pts, maxpts)
        assert abs(pts - maxpts * 0.60) < 1e-9, (g, pts)


def test_prioritize_overlay_excludes_non_germline_evidence():
    """A somatic driver must contribute ZERO germline prior, whatever weight it carries.

    The GCT resource lists somatic drivers (KRAS, NRAS, CBL, MTOR, AKT1, BCORL1) deliberately, so
    a reader can see they were considered and excluded. Their prior_weight of 0.15 is a
    bookkeeping placeholder, NOT weak germline support — scoring it as the latter is exactly the
    misreading the resource's `do_not` list warns against. Reported with a flag, not dropped.
    """
    from hprv import prioritize as PR
    overlay = PR.parse_gene_prior_overlay([
        {"gene": "KRAS", "tier": "T4", "prior_weight": "0.15",
         "evidence_class": "somatic_driver_not_germline", "moi": "NA_somatic"},
        {"gene": "MNX1", "tier": "T1", "prior_weight": "1.0",
         "evidence_class": "rare_variant_syndromic"},
    ])
    assert overlay["KRAS"]["germline_excluded"] is True
    assert overlay["KRAS"]["points_scale"] == 0.0
    assert "KRAS" in overlay, "the row is REPORTED, not silently dropped — it was considered"
    assert overlay["MNX1"]["germline_excluded"] is False
    row = {"impact": "HIGH", "consequence": "stop_gained", "grpmax_af": "2e-6",
           "inheritance": "dominant"}
    sc = PR.score_variant(row, {}, {}, gene_prior=overlay["KRAS"])
    assert sc["pts_gene_list_prior"] == 0.0
    assert sc["gene_list_prior_excluded_non_germline"] is True
    assert sc["gene_list_prior_member"] is True, "membership is still visible to the reviewer"
    assert sc["gene_list_prior_tier"] == "T4"
    # a non-germline row is NOT rescued by set membership either
    ov2 = PR.parse_gene_prior_overlay(
        [{"gene": "KRAS", "prior_weight": "0.15",
          "evidence_class": "somatic_driver_not_germline", "gene_sets": "SOME_SET"}],
        gene_sets={"SOME_SET": {"prior_weight": 1.0, "members": ["KRAS"]}})
    assert ov2["KRAS"]["points_scale"] == 0.0, "set membership must not resurrect a somatic driver"
    # the germline row still scores normally
    assert PR.score_variant(row, {}, {}, gene_prior=overlay["MNX1"])["pts_gene_list_prior"] > 0.0
    # THE EXCLUSION MATCHES ON evidence_class, NEVER ON A REMEMBERED GENE LIST. KIT is the trap:
    # it is frequently named among germ-cell-tumour somatic drivers, but its row in the GCT
    # resource is tier T3 / gwas_common_variant_locus / weight 0.35 — a GWAS-locus row this
    # exclusion must NOT touch, because it carries a legitimate (if weak) T3 prior. A gene can be
    # a somatic driver in the literature and a GWAS-locus row in the file; only the file decides.
    kit = PR.parse_gene_prior_overlay([
        {"gene": "KIT", "tier": "T3", "prior_weight": "0.35",
         "evidence_class": "gwas_common_variant_locus", "moi": "complex_common_variant"}])["KIT"]
    assert kit["germline_excluded"] is False, \
        "KIT is a GWAS-locus row in this resource, not a somatic_driver_not_germline row"
    assert abs(kit["points_scale"] - 0.35) < 1e-9, kit["points_scale"]
    assert PR.score_variant(row, {}, {}, gene_prior=kit)["pts_gene_list_prior"] > 0.0
    # ...and the class name must match EXACTLY the string the resource uses — a near-miss class
    # (a plausible-looking synonym) must not silently exclude anything
    near = PR.parse_gene_prior_overlay([
        {"gene": "GENEX", "prior_weight": "0.5", "evidence_class": "somatic_driver"}])["GENEX"]
    assert near["germline_excluded"] is False, \
        "'somatic_driver' is not 'somatic_driver_not_germline' — do not fuzzy-match evidence classes"


def test_prioritize_header_row_is_never_a_gene():
    """A literal column header must never be read as a gene symbol.

    Not fastidiousness: the validation cohort's own per-gene counts file contains a literal `gene`
    line CARRYING n=1, which is why the true totals are 25,389 variants / 10,799 symbols rather
    than the 25,390 / 10,800 raw line count. A header read as a gene is a phantom gene with a real
    count in every table it reaches.
    """
    from hprv import prioritize as PR
    ov = PR.parse_gene_prior_overlay(["gene", "symbol", "#gene", "TP53", "gene_symbol", "BRCA2"])
    assert set(ov) == {"TP53", "BRCA2"}, set(ov)
    ov2 = PR.parse_gene_prior_overlay([{"gene": "gene", "prior_weight": "1.0"},
                                       {"gene": "TP53", "prior_weight": "1.0"}])
    assert set(ov2) == {"TP53"}
    for tok in ("gene", "symbol", "gene_symbol", "gene_name", "#gene"):
        assert tok in PR.GENE_KEYS_LOWER, tok


def test_prioritize_config_overlay_path_is_gated_on_enabled():
    """A config overlay PATH must be inert while `enabled: false` — that is the Class-A contract.

    Regression test for a real leak: the CLI fell back to
    prioritization.composite.gene_list_prior.path unconditionally, so a stale or leftover path
    silently promoted genes on a run the operator believed was phenotype-agnostic. An explicit
    --gene-prior IS the intent and still wins; a config path needs `enabled: true`.
    """
    import shutil
    p9 = _load_p9()
    d = tempfile.mkdtemp(prefix="_hprv_p9ov_")
    try:
        mut = os.path.join(d, "mut.tsv")
        _write_tsv(mut, ["gene", "mu_mis", "mu_syn", "mu_lof"],
                   [{"gene": f"BG{i:03d}", "mu_mis": "9e-6", "mu_syn": "4e-6", "mu_lof": "5e-7"}
                    for i in range(50)])
        # The overlay must target genes that rank LAST without it, or a rank_delta of 0 would be
        # ambiguous: every row here scores identically, so ties break by position, and promoting
        # the genes already at ranks 1-2 could not move them.
        overlay = os.path.join(d, "overlay.txt")
        with open(overlay, "w") as fh:
            fh.write("BG048\nBG049\n")
        cols = ["chrom", "pos", "ref", "alt", "trio_id", "gene", "consequence", "impact",
                "inheritance", "grpmax_af", "cadd", "spliceai_ds"]
        rows = [{"chrom": "chr1", "pos": str(100 + i), "ref": "A", "alt": "T", "trio_id": "T1",
                 "gene": f"BG{i:03d}", "consequence": "stop_gained", "impact": "HIGH",
                 "inheritance": "dominant", "grpmax_af": "2e-6", "cadd": "38",
                 "spliceai_ds": ""} for i in range(50)]
        vin = os.path.join(d, "variants.tsv")
        _write_tsv(vin, cols, rows)

        def run(enabled, extra=()):
            cfgp = os.path.join(d, f"cfg_{enabled}.yaml")
            with open(cfgp, "w") as fh:
                fh.write("project: {name: t}\nprioritization:\n  composite:\n"
                         f"    gene_list_prior: {{enabled: {str(enabled).lower()}, "
                         f"path: {overlay}}}\n")
            outv = os.path.join(d, f"vp_{enabled}_{len(extra)}.tsv")
            outg = os.path.join(d, f"gp_{enabled}_{len(extra)}.tsv")
            rc = p9.main(["--variants", vin, "--mutrate", mut, "--config", cfgp,
                          "--out-variants", outv, "--out-genes", outg] + list(extra))
            assert rc == 0, rc
            import csv as _csv
            with open(outv) as fh:
                return list(_csv.DictReader(fh, delimiter="\t"))

        # enabled: false + a config path -> the overlay is INERT and the rankings are identical
        off = run(False)
        assert all(float(r["pts_gene_list_prior"]) == 0.0 for r in off), \
            "a config overlay path must be inert while enabled: false"
        assert all(r["rank_prior"] == r["rank_agnostic"] for r in off)
        assert all(r["rank_delta"] == "0" for r in off)
        assert all(r["gene_list_prior_member"] == "0" for r in off)
        # enabled: true -> the overlay applies, but ONLY to the prior ranking
        on = run(True)
        promoted = [r for r in on if r["gene"] in ("BG048", "BG049")]
        assert len(promoted) == 2 and all(float(r["pts_gene_list_prior"]) == 2.0
                                         for r in promoted), promoted
        by_pos_off = {r["pos"]: r for r in off}
        assert all(r["priority_points_agnostic"] == by_pos_off[r["pos"]]["priority_points_agnostic"]
                   for r in on), "enabling the overlay must not move the AGNOSTIC total"
        # rank_delta must EXPOSE the promotion — a large negative delta is the confirmation-bias
        # review set, and it would be invisible under a single blended ranking
        assert all(int(r["rank_delta"]) < 0 for r in promoted), \
            f"overlay members must show a negative rank_delta ({[r['rank_delta'] for r in promoted]})"
        assert all(int(r["rank_prior"]) <= 2 for r in promoted), \
            "the +2 must lift the overlay members to the top of the PRIOR ranking"
        assert all(int(r["rank_agnostic"]) >= 49 for r in promoted), \
            "...while their AGNOSTIC ranks stay at the bottom, unmoved by list membership"
        # an EXPLICIT --gene-prior wins even with enabled: false (the flag is the intent)
        forced = run(False, ("--gene-prior", overlay))
        assert sum(1 for r in forced if float(r["pts_gene_list_prior"]) == 2.0) == 2
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_prioritize_reads_bgzipped_tables():
    """Step 9 must read a GZIPPED mutational-target table — the prepared resource IS one.

    prepare_resources.sh writes `constraint/mutational_target.by_gene.txt.bgz` (the gnomAD v2.1.1
    constraint table kept unjoined; gnomAD ships it .bgz) and emit-env exports that path straight
    into --mutrate. A plain open() there hands csv a binary stream, so a real run would crash or
    silently lose the excess statistic — while the integration mock uses a plain .tsv and would
    never notice. This is that regression test. bgzip is gzip-compatible, so stdlib gzip suffices
    and the test needs no htslib.
    """
    import gzip as _gzip
    import shutil
    p9 = _load_p9()
    d = tempfile.mkdtemp(prefix="_hprv_p9gz_")
    try:
        header = ("gene\tmu_mis\tmu_syn\tmu_lof\toe_syn\tclassic_caf\tconstraint_flag\t"
                  "cds_length\tpli\toe_lof_upper\n")
        body = "".join(f"BG{i:03d}\t9e-6\t4e-6\t5e-7\t1.0\t2e-4\t\t1500\t0.1\t1.1\n"
                       for i in range(200))
        body += "OR4Q3\t2e-7\t8e-8\t1e-8\t1.7\t0\tmis_too_many\t900\t0.01\t1.9\n"
        gz = os.path.join(d, "mutational_target.by_gene.txt.bgz")
        with _gzip.open(gz, "wt") as fh:
            fh.write(header + body)
        # a GZIPPED control union too — the same reader serves both
        ctrl_gz = os.path.join(d, "ctrl.txt.gz")
        with _gzip.open(ctrl_gz, "wt") as fh:
            fh.write("# established-gene union\nBG000\nBG001\n")
        cols = ["chrom", "pos", "ref", "alt", "trio_id", "gene", "consequence", "impact",
                "inheritance", "grpmax_af", "cadd", "spliceai_ds", "child_GQ", "child_DP",
                "child_AB"]
        rows = []
        for i in range(200):
            rows.append({"chrom": "chr1", "pos": str(1000 + i * 7), "ref": "A", "alt": "T",
                         "trio_id": f"T{i % 10:02d}", "gene": f"BG{i:03d}",
                         "consequence": "missense_variant", "impact": "MODERATE",
                         "inheritance": "dominant", "grpmax_af": "2e-6", "cadd": "12",
                         "spliceai_ds": "0.02", "child_GQ": "99", "child_DP": "40",
                         "child_AB": "0.5"})
        for k in range(12):                              # the artifact pileup
            rows.append(dict(rows[0], pos=str(50000 + k * 7), gene="OR4Q3",
                             trio_id=f"T{k % 10:02d}"))
        vin = os.path.join(d, "variants.tsv")
        _write_tsv(vin, cols, rows)
        cfgp = os.path.join(d, "cfg.yaml")
        with open(cfgp, "w") as fh:
            fh.write("project: {name: t}\nprioritization:\n"
                     "  gene_downweight: {min_control_genes: 2, max_downweight_fraction: 1.0}\n")
        outv, outg = os.path.join(d, "vp.tsv"), os.path.join(d, "gp.tsv")
        rc = p9.main(["--variants", vin, "--mutrate", gz, "--established-genes", ctrl_gz,
                      "--config", cfgp, "--n-trios", "10",
                      "--out-variants", outv, "--out-genes", outg])
        assert rc == 0, f"Step 9 failed on a bgzipped mutational-target table (rc={rc})"
        import csv as _csv
        with open(outg) as fh:
            pg = {r["gene"]: r for r in _csv.DictReader(fh, delimiter="\t")}
        # the offset was actually READ — not silently degraded to "no excess statistic"
        assert pg["OR4Q3"]["E_source"] == "gnomad_mu", pg["OR4Q3"]["E_source"]
        assert float(pg["OR4Q3"]["excess_ratio"]) > 10.0, pg["OR4Q3"]["excess_ratio"]
        assert pg["OR4Q3"]["q_nb"] not in ("", None)
        assert pg["OR4Q3"]["gene_tier"] in ("T2_downweight", "T3_strong_downweight")
        # and the gzipped control union was read too
        assert pg["BG000"]["established_gene_control"] == "1"
        assert pg["BG002"]["established_gene_control"] == "0"
        # never-drop still holds on this path
        with open(outv) as fh:
            assert len(list(_csv.DictReader(fh, delimiter="\t"))) == len(rows)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_prioritize_config_matches_canonical_defaults():
    """SINGLE SOURCE OF TRUTH: a threshold in code that disagrees with the shipped config is a
    bug. Every default prioritize.py falls back to must equal config.example.yaml's value."""
    from hprv import prioritize as PR
    try:
        from hprv.config import load_config
        cfg = load_config(os.path.join(os.path.dirname(__file__), "..", "config",
                                       "config.example.yaml"))
    except ImportError:                     # pyyaml absent on a bare host — skip, don't fail
        return
    checks = [
        ("prioritization.excess.trim_p", 1.0e-3),
        ("prioritization.excess.trim_max_fraction", 0.05),
        ("prioritization.excess.q_threshold", 0.05),
        ("prioritization.excess.min_n_for_ratio_rule", 3),
        ("prioritization.excess.offset.mu_lof_impute_factor", 0.0516),
        ("prioritization.excess.offset.cds_fallback.b0", -18.4484),
        ("prioritization.excess.offset.cds_fallback.b1", 1.0570),
        ("prioritization.signals.saturation.per_trio_min", 0.10),
        ("prioritization.signals.segdup.min_frac", 0.10),
        ("prioritization.signals.segdup.min_identity", 0.98),
        ("prioritization.signals.oe_syn.max_deviation", 0.30),
        ("prioritization.signals.caf_low.percentile", 0.10),
        ("prioritization.gene_downweight.t1_watch.min_ratio", 3.0),
        ("prioritization.gene_downweight.t1_watch.max_q", 0.25),
        ("prioritization.gene_downweight.t2_downweight.min_ratio", 5.0),
        ("prioritization.gene_downweight.t2_downweight.min_corroboration", 1),
        ("prioritization.gene_downweight.t3_strong.max_q", 0.05),
        ("prioritization.gene_downweight.t3_strong.min_ratio", 10.0),
        ("prioritization.gene_downweight.t3_strong.ratio_min_corroboration", 2),
        ("prioritization.gene_downweight.min_control_genes", 1000),
        ("prioritization.gene_downweight.max_downweight_fraction", 0.20),
        ("prioritization.variant_tier.spliceai_strong", 0.5),
        ("prioritization.variant_tier.spliceai_supporting", 0.2),
        ("prioritization.variant_tier.spliceai_benign_max", 0.1),
        ("prioritization.variant_tier.cadd_benign_max", 15.0),
        ("prioritization.variant_tier.cadd_missense_supporting", 25.3),
        # missense predictor ladder — REVEL's cuts are Pejaver 2022's, AlphaMissense's are the
        # model's own (Cheng 2023). Pinned here so a code default cannot drift from the config.
        ("prioritization.variant_tier.revel_supporting", 0.644),
        ("prioritization.variant_tier.revel_moderate", 0.773),
        ("prioritization.variant_tier.revel_benign_max", 0.290),
        ("prioritization.variant_tier.alphamissense_supporting", 0.564),
        ("prioritization.variant_tier.alphamissense_benign_max", 0.34),
        ("resources.clinvar.min_review_stars", 2),
        ("resources.clinvar.low_star_scale", 0.5),
        ("prioritization.variant_tier.rarity.ba1", 0.05),
        ("prioritization.variant_tier.rarity.moderate_max", 1.0e-4),
        ("prioritization.variant_tier.rarity.single_group_ratio", 10.0),
    ]
    for key, want in checks:
        got = get(cfg, key, "__MISSING__")
        assert got != "__MISSING__", f"{key} is absent from config.example.yaml"
        assert abs(float(got) - float(want)) < 1e-12, f"{key}: config {got} != code default {want}"
    # the SpliceAI / CADD / rarity numbers must also agree with the SCREEN's own canonical values
    assert get(cfg, "prioritization.variant_tier.spliceai_supporting") == \
        get(cfg, "filters.functional.spliceai_ds_min")
    assert get(cfg, "prioritization.variant_tier.cadd_missense_supporting") == \
        get(cfg, "filters.functional.cadd_phred_supporting")
    assert get(cfg, "prioritization.variant_tier.rarity.ba1") == \
        get(cfg, "filters.rarity.benign_ba1")
    assert get(cfg, "prioritization.variant_tier.rarity.moderate_max") == \
        get(cfg, "filters.rarity.dominant_max")
    # the config weights must match default_weights() term by term
    w_code = PR.default_weights()
    w_cfg = get(cfg, "prioritization.composite.weights")
    for tier, pts in w_cfg["molecular"].items():
        assert float(pts) == w_code["molecular"][tier], tier
    for tier, pts in w_cfg["gene_artifact"].items():
        assert float(pts) == w_code["gene_artifact"][tier], tier
    assert float(w_cfg["rarity"]["ba1"]) == w_code["rarity"]["fail"]
    # ...and the Class-B overlay must be OFF in the shipped config (that is the whole contract)
    assert get(cfg, "prioritization.composite.gene_list_prior.enabled") is False



def test_prioritize_emits_raw_source_columns():
    """Every merged source column must ride through VERBATIM, prefixed by its source.

    The point is reproducibility without hidden state: a collaborator holding
    variants.prioritized.tsv should be able to re-derive excess_ratio from the raw mu_* cells and
    filter on a source field the scoring never interprets (pLI, mis_z, the per-ancestry
    classic_caf_* columns) without re-running the join. Three failure modes this guards:
      * a source column silently dropped at write time (the pre-fix behaviour),
      * a raw value MUTATED on the way through (it must be byte-identical to the input cell),
      * a gene absent from a source reading as 0.0 rather than MISSING.
    Also asserts the per-source prefix keeps same-named columns from two tables distinct.
    """
    import csv as _csv
    import shutil
    p9 = _load_p9()
    d = tempfile.mkdtemp(prefix="_hprv_p9src_")
    try:
        # mutrate carries a column the SCORING never reads (mis_z) and one it does (pli).
        mut = os.path.join(d, "mut.tsv")
        with open(mut, "w") as fh:
            fh.write("gene\tmu_mis\tmu_syn\tmu_lof\toe_syn\tclassic_caf\tconstraint_flag\t"
                     "cds_length\tpli\toe_lof_upper\tmis_z\tclassic_caf_afr\n")
            for i in range(200):
                fh.write(f"BG{i:03d}\t9e-6\t4e-6\t5e-7\t1.0\t2e-4\t\t1500\t0.1\t1.1\t"
                         f"{i / 100.0}\t3.5e-5\n")
            fh.write("OR4Q3\t2e-7\t8e-8\t1e-8\t1.7\t0\tmis_too_many\t900\t0.01\t1.9\t"
                     "-2.75\t0\n")
        # A SECOND table with a same-named column (pli) — the prefixes must keep them apart, and
        # BG199 is deliberately ABSENT from it so the missing-source path is exercised.
        con = os.path.join(d, "con.tsv")
        with open(con, "w") as fh:
            fh.write("gene\tpli\ts_het\n")
            for i in range(199):
                fh.write(f"BG{i:03d}\t0.97\t0.08\n")
        cols = ["chrom", "pos", "ref", "alt", "trio_id", "gene", "consequence", "impact",
                "inheritance", "grpmax_af", "cadd", "spliceai_ds", "child_GQ", "child_DP",
                "child_AB"]
        rows = []
        for i in range(200):
            rows.append({"chrom": "chr1", "pos": str(1000 + i * 7), "ref": "A", "alt": "T",
                         "trio_id": f"T{i % 10:02d}", "gene": f"BG{i:03d}",
                         "consequence": "missense_variant", "impact": "MODERATE",
                         "inheritance": "dominant", "grpmax_af": "2e-6", "cadd": "12",
                         "spliceai_ds": "0.02", "child_GQ": "99", "child_DP": "40",
                         "child_AB": "0.5"})
        vin = os.path.join(d, "variants.tsv")
        _write_tsv(vin, cols, rows)
        cfgp = os.path.join(d, "cfg.yaml")
        with open(cfgp, "w") as fh:
            fh.write("project: {name: t}\nprioritization:\n"
                     "  gene_downweight: {min_control_genes: 0, max_downweight_fraction: 1.0}\n")
        outv, outg = os.path.join(d, "vp.tsv"), os.path.join(d, "gp.tsv")
        rc = p9.main(["--variants", vin, "--mutrate", mut, "--constraint", con,
                      "--config", cfgp, "--n-trios", "10",
                      "--out-variants", outv, "--out-genes", outg])
        assert rc == 0, f"Step 9 failed with source pass-through (rc={rc})"

        with open(outg) as fh:
            grows = {r["gene"]: r for r in _csv.DictReader(fh, delimiter="\t")}
        g = grows["BG000"]
        # (a) a column the scoring NEVER interprets survived, byte-identical
        assert "src_mutrate_mis_z" in g, "raw mutrate columns were dropped"
        assert g["src_mutrate_mis_z"] == "0.0", g["src_mutrate_mis_z"]
        assert g["src_mutrate_classic_caf_afr"] == "3.5e-5", g["src_mutrate_classic_caf_afr"]
        # (b) the raw mu_* cells are present verbatim, so E = C*mu is re-derivable off this file
        assert g["src_mutrate_mu_mis"] == "9e-6", g["src_mutrate_mu_mis"]
        # (c) same-named columns from two sources stay DISTINCT and keep their own values
        assert g["src_mutrate_pli"] == "0.1", g["src_mutrate_pli"]
        assert g["src_constraint_pli"] == "0.97", g["src_constraint_pli"]
        # (d) the gene KEY is not duplicated into the pass-through block
        assert "src_mutrate_gene" not in g
        # (e) a gene ABSENT from a source reads MISSING, never 0.0 — the same contract as NHF
        assert grows["BG199"]["src_constraint_pli"] == "", grows["BG199"]["src_constraint_pli"]
        assert grows["BG199"]["src_mutrate_pli"] == "0.1"
        # (f) the curated column is still the one the SCORE used, and is not clobbered
        assert g["pLI"] not in ("", None)

        # (g) the per-variant file is self-contained: the gene's raw cells ride onto its variants
        with open(outv) as fh:
            vrows = list(_csv.DictReader(fh, delimiter="\t"))
        assert len(vrows) == len(rows), "never-drop violated on the source-column path"
        v0 = next(r for r in vrows if r["gene"] == "BG000")
        assert v0["src_mutrate_mis_z"] == "0.0", v0["src_mutrate_mis_z"]
        assert v0["src_constraint_pli"] == "0.97"

        # (h) --no-source-columns actually suppresses them (and nothing else breaks)
        outv2, outg2 = os.path.join(d, "vp2.tsv"), os.path.join(d, "gp2.tsv")
        rc = p9.main(["--variants", vin, "--mutrate", mut, "--constraint", con,
                      "--config", cfgp, "--n-trios", "10", "--no-source-columns",
                      "--out-variants", outv2, "--out-genes", outg2])
        assert rc == 0
        with open(outg2) as fh:
            hdr2 = fh.readline().rstrip("\n").split("\t")
        assert not any(c.startswith("src_") for c in hdr2), "--no-source-columns did not suppress"
        assert "pLI" in hdr2, "curated columns must survive --no-source-columns"
    finally:
        shutil.rmtree(d, ignore_errors=True)



def test_prioritize_emits_raw_gene_prior_columns():
    """The --gene-prior overlay's OWN columns must ride through verbatim too.

    Tier and weight alone do not tell a reviewer WHY a gene carries a prior. The GCT resource's
    site_specificity / anatomical_transfer / replication / pmids columns are exactly what decides
    whether a phenotype prior transfers to a given cohort — e.g. a GWAS-locus row whose caveat
    says rare coding variants are NOT the implicated mechanism. Dropping them leaves a bare
    weight the reviewer cannot audit, which is the failure this guards.

    Also covers the two overlay-specific wrinkles: entries are keyed upper-cased internally, and a
    gene contributed by SET membership alone has no source row (blank raw cells) while still
    carrying its set attribution.
    """
    import csv as _csv
    import shutil
    p9 = _load_p9()
    d = tempfile.mkdtemp(prefix="_hprv_p9prior_")
    try:
        mut = os.path.join(d, "mut.tsv")
        with open(mut, "w") as fh:
            fh.write("gene\tmu_mis\tmu_syn\tmu_lof\toe_syn\tclassic_caf\tconstraint_flag\t"
                     "cds_length\tpli\toe_lof_upper\n")
            for i in range(200):
                fh.write(f"BG{i:03d}\t9e-6\t4e-6\t5e-7\t1.0\t2e-4\t\t1500\t0.1\t1.1\n")
        # An overlay shaped like the real GCT resource: the reviewer-critical caveat columns are
        # the point of the test.
        prior = os.path.join(d, "prior.tsv")
        with open(prior, "w") as fh:
            fh.write("gene\ttier\tprior_weight\tevidence_class\tsite_specificity\t"
                     "replication\tanatomical_transfer\tpmids\tnotes\n")
            fh.write("BG000\tT2\t0.75\trare_variant_association\tgonadal_tgct\t"
                     "replicated_independently_x2\tintracranial_transfer_untested\t"
                     "30676620;32451744\tthe replicated one\n")
            fh.write("BG001\tT3\t0.35\tgwas_locus\tgonadal_tgct;intracranial\t"
                     "gwas_meta_analysed\tCOMMON_VARIANT_LOCUS__weak_evidence_for_rare_coding\t"
                     "28604732\tlead variant is non-coding\n")
        cols = ["chrom", "pos", "ref", "alt", "trio_id", "gene", "consequence", "impact",
                "inheritance", "grpmax_af", "cadd", "spliceai_ds", "child_GQ", "child_DP",
                "child_AB"]
        rows = [{"chrom": "chr1", "pos": str(1000 + i * 7), "ref": "A", "alt": "T",
                 "trio_id": f"T{i % 10:02d}", "gene": f"BG{i:03d}",
                 "consequence": "missense_variant", "impact": "MODERATE",
                 "inheritance": "dominant", "grpmax_af": "2e-6", "cadd": "12",
                 "spliceai_ds": "0.02", "child_GQ": "99", "child_DP": "40",
                 "child_AB": "0.5"} for i in range(200)]
        vin = os.path.join(d, "variants.tsv")
        _write_tsv(vin, cols, rows)
        cfgp = os.path.join(d, "cfg.yaml")
        with open(cfgp, "w") as fh:
            fh.write("project: {name: t}\nprioritization:\n"
                     "  gene_downweight: {min_control_genes: 0, max_downweight_fraction: 1.0}\n"
                     "  composite:\n    gene_list_prior:\n      enabled: true\n")
        outv, outg = os.path.join(d, "vp.tsv"), os.path.join(d, "gp.tsv")
        rc = p9.main(["--variants", vin, "--mutrate", mut, "--gene-prior", prior,
                      "--config", cfgp, "--n-trios", "10",
                      "--out-variants", outv, "--out-genes", outg])
        assert rc == 0, f"Step 9 failed with a --gene-prior overlay (rc={rc})"

        with open(outg) as fh:
            g = {r["gene"]: r for r in _csv.DictReader(fh, delimiter="\t")}
        r0 = g["BG000"]
        # every overlay column present and byte-identical
        assert r0["src_prior_site_specificity"] == "gonadal_tgct", r0["src_prior_site_specificity"]
        assert r0["src_prior_replication"] == "replicated_independently_x2"
        assert r0["src_prior_pmids"] == "30676620;32451744", r0["src_prior_pmids"]
        assert r0["src_prior_notes"] == "the replicated one"
        # the caveat that makes a GWAS-locus prior interpretable must survive
        assert "COMMON_VARIANT_LOCUS" in g["BG001"]["src_prior_anatomical_transfer"]
        # the gene KEY is not duplicated into the pass-through block
        assert "src_prior_gene" not in r0
        # a gene NOT in the overlay reads MISSING, never 0.0 / a spurious tier
        assert g["BG050"]["src_prior_tier"] == "", g["BG050"]["src_prior_tier"]
        assert g["BG050"]["gene_list_prior_member"] == "0"
        # the derived columns still work and are NOT replaced by the raw ones
        assert r0["gene_list_prior_member"] == "1"

        # per-variant: the overlay context rides onto the call, so one file is self-contained
        with open(outv) as fh:
            vrows = list(_csv.DictReader(fh, delimiter="\t"))
        assert len(vrows) == len(rows), "never-drop violated on the gene-prior path"
        v0 = next(r for r in vrows if r["gene"] == "BG000")
        assert v0["src_prior_pmids"] == "30676620;32451744", v0["src_prior_pmids"]
        assert v0["gene_list_prior_tier"] == "T2", v0["gene_list_prior_tier"]
        v1 = next(r for r in vrows if r["gene"] == "BG001")
        assert "COMMON_VARIANT_LOCUS" in v1["src_prior_anatomical_transfer"]
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _run_all():
    import inspect
    fns = [f for n, f in sorted(globals().items())
           if n.startswith("test_") and inspect.isfunction(f)]
    skipped = []
    for f in fns:
        try:
            f()
        except _Skip as e:
            skipped.append((f.__name__, str(e)))
            print(f"SKIP {f.__name__} — needs {e}")
            continue
        print(f"PASS {f.__name__}")
    if skipped:
        # Deliberately NOT "all tests passed". A partial run must not read like a full one; the
        # skipped set here contains the never-drop and cache-invalidation guards.
        need = sorted({m for _, m in skipped})
        print(f"\n{len(fns) - len(skipped)} passed, {len(skipped)} SKIPPED for missing "
              f"dependencies ({', '.join(need)}) — this is NOT full coverage. Run "
              f"`pip install {' '.join(need).replace('yaml', 'pyyaml')}` to exercise the whole "
              f"suite (CI does).")
    else:
        print(f"\nAll {len(fns)} pure-logic tests passed.")


if __name__ == "__main__":
    _run_all()
