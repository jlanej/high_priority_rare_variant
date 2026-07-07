"""Read the annotations produced by pipeline Step 2 and evaluate genotype QC.

Step 2 lifts VEP/plugin fields to INFO with a ``vep_`` prefix (via ``bcftools
+split-vep``) and transfers external population/clinical annotations as
``hprv_gnomad_*`` / ``hprv_clnsig`` etc. This module is the single place that knows
those field names and how to coerce their (string) values, so the selection,
inheritance, and burden steps all read them identically.

All getters are defensive: they return ``None`` for missing/'.'/unparseable values
and take the max over ``&``/``,``-joined multi-transcript values for scores.
"""

from __future__ import annotations

from typing import Optional


# --- INFO field names written by Step 2 (single source of truth) -------------
F = {
    "consequence": "vep_Consequence",
    "impact": "vep_IMPACT",
    "symbol": "vep_SYMBOL",
    "gene": "vep_Gene",
    "feature": "vep_Feature",
    "biotype": "vep_BIOTYPE",
    "hgvsc": "vep_HGVSc",
    "hgvsp": "vep_HGVSp",
    "mane": "vep_MANE_SELECT",
    "revel": "vep_REVEL_score",
    "alphamissense": "vep_AlphaMissense_score",
    "alphamissense_pred": "vep_AlphaMissense_pred",
    "mpc": "vep_MPC_score",
    "metarnn": "vep_MetaRNN_score",
    # CADD may arrive from the CADD plugin (CADD_PHRED) or dbNSFP (CADD_phred)
    "cadd_a": "vep_CADD_PHRED",
    "cadd_b": "vep_CADD_phred",
    "spliceai_ag": "vep_SpliceAI_pred_DS_AG",
    "spliceai_al": "vep_SpliceAI_pred_DS_AL",
    "spliceai_dg": "vep_SpliceAI_pred_DS_DG",
    "spliceai_dl": "vep_SpliceAI_pred_DS_DL",
    "loftee": "vep_LoF",
    "loftee_filter": "vep_LoF_filter",
    "loftee_flags": "vep_LoF_flags",
    "gnomad_af": "hprv_gnomad_af",
    "gnomad_grpmax_af": "hprv_gnomad_grpmax_af",
    "gnomad_faf95": "hprv_gnomad_faf95",
    "gnomad_nhomalt": "hprv_gnomad_nhomalt",
    "clnsig": "hprv_clnsig",
    "clnrevstat": "hprv_clnrevstat",
    "clnsigconf": "hprv_clnsigconf",
}

# HIGH-impact / loss-of-function VEP consequence terms.
LOF_CONSEQUENCES = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification",
}

_STAR = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_single_submitter": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_assertion_provided": 0,
}


def _raw(variant, key: str):
    try:
        return variant.INFO.get(F[key])
    except KeyError:
        return None


def _str(variant, key: str) -> Optional[str]:
    v = _raw(variant, key)
    if v is None:
        return None
    s = str(v)
    return None if s in ("", ".") else s


def _max_float(variant, *keys) -> Optional[float]:
    """Max float over the given field(s), splitting &/, multi-values; None if none."""
    best = None
    for key in keys:
        v = _raw(variant, key)
        if v is None:
            continue
        for tok in str(v).replace("&", ",").split(","):
            tok = tok.strip()
            if tok in ("", "."):
                continue
            try:
                f = float(tok)
            except ValueError:
                continue
            best = f if best is None else max(best, f)
    return best


# --- population frequency ----------------------------------------------------
def faf95(variant) -> Optional[float]:
    return _max_float(variant, "gnomad_faf95")


def grpmax_af(variant) -> Optional[float]:
    return _max_float(variant, "gnomad_grpmax_af")


def nhomalt(variant) -> Optional[float]:
    return _max_float(variant, "gnomad_nhomalt")


def frequency(variant) -> Optional[float]:
    """Preferred rarity field: faf95, falling back to grpmax AF, then global AF."""
    for fn in (faf95, grpmax_af):
        v = fn(variant)
        if v is not None:
            return v
    return _max_float(variant, "gnomad_af")


# --- functional predictors ---------------------------------------------------
def revel(variant) -> Optional[float]:
    return _max_float(variant, "revel")


def alphamissense(variant) -> Optional[float]:
    return _max_float(variant, "alphamissense")


def mpc(variant) -> Optional[float]:
    return _max_float(variant, "mpc")


def cadd(variant) -> Optional[float]:
    return _max_float(variant, "cadd_a", "cadd_b")


def spliceai_max(variant) -> Optional[float]:
    return _max_float(variant, "spliceai_ag", "spliceai_al", "spliceai_dg", "spliceai_dl")


def impact(variant) -> Optional[str]:
    return _str(variant, "impact")


def symbol(variant) -> Optional[str]:
    return _str(variant, "symbol")


def consequence(variant) -> Optional[str]:
    return _str(variant, "consequence")


def is_loftee_hc(variant) -> bool:
    return (_str(variant, "loftee") or "").upper() == "HC"


def loftee_flags(variant) -> Optional[str]:
    return _str(variant, "loftee_flags")


# --- clinical ----------------------------------------------------------------
def clnsig(variant) -> Optional[str]:
    return _str(variant, "clnsig")


def clnsig_is_plp(variant) -> bool:
    s = (clnsig(variant) or "").lower()
    if not s:
        return False
    # Treat conflicting as NOT P/LP even if the token appears in the conflict string.
    if "conflicting" in s:
        return False
    return "pathogenic" in s and "likely_benign" not in s and "benign/likely" not in s


def clinvar_stars(variant) -> int:
    s = _str(variant, "clnrevstat")
    if not s:
        return 0
    return _STAR.get(s.strip().lower(), 0)
