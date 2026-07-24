"""Read the annotations produced by pipeline Step 2 and evaluate genotype QC.

**VEP-centric contract.** Every annotation this pipeline reads comes from ONE tool: VEP 115
(GRCh38) — its cache plus its score plugins — whose CSQ fields Step 2 lifts to INFO with a
``vep_`` prefix (via ``bcftools +split-vep``). No external sites VCF is bcftools-transferred in
— no gnomAD, ClinVar, dbNSFP or LOFTEE file is annotated in. The plugins are CADD (genome-wide
functional, SNV+indel) and — optionally — SpliceAI (precomputed raw genome-wide splice delta
scores); population frequency and ClinVar ride in the cache itself. This module is the single
place that knows those field names and how to coerce their (string) values, so the selection,
inheritance, and burden steps all read them identically.

What that costs is documented in docs/allele_frequency.md and docs/functional_annotation.md;
in short: no faf95 (no CI correction), no nhomalt, no LOFTEE, no ClinVar review status. SpliceAI
IS available now (when its score files are provided) — but the PRECOMPUTED set has its own limits
(a missing score is NOT evidence of no effect; see spliceai_ds()). Adding another annotation means
adding its INFO field here AND its plugin/transfer in Step 2 — nothing else reaches around this.

All getters are defensive: they return ``None`` for missing/'.'/unparseable values
and take the max over ``&``/``,``-joined multi-transcript values for scores.
"""

from __future__ import annotations

from typing import Optional


# gnomAD v4.1 genetic-ancestry groups that gnomAD's OWN grpmax includes. grpmax
# deliberately EXCLUDES asj / fin / mid / ami / remaining: bottlenecked founder groups
# whose small ANs make a point-estimate AF wildly unrepresentative of any general
# population (ami AN ~900 — a single allele reads as AF ~1e-3, ten-fold over the
# dominant gate). Mirroring that inclusion set is precisely what makes the per-population
# CSQ max a defensible grpmax proxy. VEP's own MAX_AF is NOT usable here: it maximises
# over all gnomAD groups AND the tiny 1000 Genomes phase-3 populations, reintroducing
# exactly the false-negative mode grpmax exists to prevent. See docs/allele_frequency.md.
GRPMAX_POPS = ("AFR", "AMR", "EAS", "NFE", "SAS")

# --- INFO field names written by Step 2 (single source of truth) -------------
F = {
    # --- core VEP consequence block ---
    "consequence": "vep_Consequence",
    "impact": "vep_IMPACT",
    "symbol": "vep_SYMBOL",
    "gene": "vep_Gene",
    "feature": "vep_Feature",
    "biotype": "vep_BIOTYPE",
    "hgvsc": "vep_HGVSc",
    "hgvsp": "vep_HGVSp",
    "mane": "vep_MANE_SELECT",
    # --- functional prediction ---
    # CADD from the dedicated plugin (CSQ CADD_PHRED -> vep_CADD_PHRED via split-vep):
    # genome-wide, SNV+indel, and under this contract the ONLY functional predictor.
    # It is therefore the sole keep-path for anything VEP rates below MODERATE. CADD
    # v1.6+ ingests SpliceAI/MMSplice as input features, so it carries a lossy
    # re-encoding of the splice signal the SpliceAI plugin would have supplied.
    "cadd": "vep_CADD_PHRED",
    # --- splice prediction (SpliceAI plugin; precomputed raw genome-wide scores) ---
    # Four per-event delta scores in [0,1]: Acceptor/Donor Gain/Loss. spliceai_ds() takes the MAX
    # = the standard SpliceAI "delta score" used for thresholding. DP_* are the predicted cryptic-
    # site positions (bp offset from the variant), carried for curator context. These are the CSQ
    # subfields the SpliceAI VEP plugin emits in VCF output (lifted verbatim by split-vep).
    "spliceai_ds_ag": "vep_SpliceAI_pred_DS_AG",
    "spliceai_ds_al": "vep_SpliceAI_pred_DS_AL",
    "spliceai_ds_dg": "vep_SpliceAI_pred_DS_DG",
    "spliceai_ds_dl": "vep_SpliceAI_pred_DS_DL",
    "spliceai_dp_ag": "vep_SpliceAI_pred_DP_AG",
    "spliceai_dp_al": "vep_SpliceAI_pred_DP_AL",
    "spliceai_dp_dg": "vep_SpliceAI_pred_DP_DG",
    "spliceai_dp_dl": "vep_SpliceAI_pred_DP_DL",
    "spliceai_symbol": "vep_SpliceAI_pred_SYMBOL",
    # --- clinical ---
    # ClinVar significance as cached by VEP (--check_existing, via --everything).
    # The cache exposes CLIN_SIG ONLY: there is no review status (CLNREVSTAT), so star
    # ratings are unavailable and the >=2-star auto-promote gate cannot be applied.
    # Values are lowercase, '&'-joined (e.g. "pathogenic&likely_pathogenic").
    "clnsig": "vep_CLIN_SIG",
    # --- population frequency (gnomAD v4.1, cached; --af_gnomade / --af_gnomadg) ---
    # POINT ESTIMATES. The cache carries no AC/AN, so faf95's CI correction is not
    # reconstructible from them at any cost — it is simply absent, not approximated.
    "gnomade_afr_af": "vep_gnomADe_AFR_AF",
    "gnomade_amr_af": "vep_gnomADe_AMR_AF",
    "gnomade_eas_af": "vep_gnomADe_EAS_AF",
    "gnomade_nfe_af": "vep_gnomADe_NFE_AF",
    "gnomade_sas_af": "vep_gnomADe_SAS_AF",
    "gnomadg_afr_af": "vep_gnomADg_AFR_AF",
    "gnomadg_amr_af": "vep_gnomADg_AMR_AF",
    "gnomadg_eas_af": "vep_gnomADg_EAS_AF",
    "gnomadg_nfe_af": "vep_gnomADg_NFE_AF",
    "gnomadg_sas_af": "vep_gnomADg_SAS_AF",
    # Global AFs + MAX_AF: REPORTING ONLY, never a filter field. Global AF dilutes an
    # ancestry-enriched benign variant across the whole cohort (false-positive
    # retention); MAX_AF over-counts founder groups (false-negative loss). Both failure
    # modes are real and they run in opposite directions — hence grpmax_af() below.
    "gnomade_af": "vep_gnomADe_AF",
    "gnomadg_af": "vep_gnomADg_AF",
    "max_af": "vep_MAX_AF",
    "max_af_pops": "vep_MAX_AF_POPS",
    # GATK PossibleDeNovo tags carried from the source trio VCF (value = comma-delimited
    # list of child sample IDs for which this is a candidate de novo). Not a VEP field.
    "hiconf_denovo": "hiConfDeNovo",
    "loconf_denovo": "loConfDeNovo",
}

# The grpmax-eligible AF fields, in F-key form, for grpmax_af()'s max.
_GRPMAX_KEYS = tuple(
    f"gnomad{src}_{pop.lower()}_af" for src in ("e", "g") for pop in GRPMAX_POPS
)

# HIGH-impact / loss-of-function VEP consequence terms.
LOF_CONSEQUENCES = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification",
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
def grpmax_af(variant) -> Optional[float]:
    """Max gnomAD v4.1 AF over the grpmax-ELIGIBLE ancestry groups (see GRPMAX_POPS).

    A point estimate standing in for gnomAD's published grpmax AF. It is NOT faf95:
    faf95 is the lower bound of the 95% CI, and computing it needs AC/AN, which the VEP
    cache does not carry. So this runs ~one CI-width HIGH on low-AC observations, and a
    rarity gate driven by it fires slightly more often than a faf95 gate would (i.e. it
    errs toward dropping). Excluding the bottlenecked groups removes the large half of
    that error; the residual is bounded by AC and documented in docs/allele_frequency.md.
    """
    return _max_float(variant, *_GRPMAX_KEYS)


def frequency(variant) -> Optional[float]:
    """The rarity field every gate reads — the single chokepoint for population frequency.

    None => no grpmax-eligible group reports this allele => treat as rarest. Note that
    under the VEP-cache contract, absence is weaker evidence than it was with a gnomAD
    sites VCF: the cache only carries frequencies for alleles accessioned into dbSNP, so
    an un-accessioned gnomAD variant silently returns no AF and reads as 'absent'. That
    biases toward retention (extra review), not toward missed calls.
    """
    return grpmax_af(variant)


# --- functional predictors ---------------------------------------------------
# CADD is the only one available under the VEP-only contract. REVEL / AlphaMissense /
# MPC / MetaRNN (dbNSFP) and SpliceAI are gone with their resource files; their getters
# and their Step-3 branches were removed rather than left to return None forever. Note
# the missense trio was already inert BEFORE removal: they are missense-only scores, and
# every missense is IMPACT=MODERATE, which selection.py keeps at an earlier branch.
def cadd(variant) -> Optional[float]:
    return _max_float(variant, "cadd")


def spliceai_ds(variant) -> Optional[float]:
    """Max SpliceAI delta score over the four events (acceptor/donor gain/loss), or None.

    The standard SpliceAI "delta score" (0-1): the probability that the variant alters splicing at
    the most-affected of the four possible splice-site changes. Thresholding this is the splice
    keep-path (filters.functional.spliceai_ds_min; ClinGen SVI uses >= 0.2 for PP3-supporting).
    None means SpliceAI did NOT score the variant — the precomputed raw set covers genome-wide SNVs
    + a large indel set, but not every possible indel, and a **missing score is not evidence of no
    splice effect** (it never drops a variant; it only fails to rescue one). See
    docs/functional_annotation.md.
    """
    return _max_float(variant, "spliceai_ds_ag", "spliceai_ds_al",
                      "spliceai_ds_dg", "spliceai_ds_dl")


def impact(variant) -> Optional[str]:
    return _str(variant, "impact")


def symbol(variant) -> Optional[str]:
    return _str(variant, "symbol")


def consequence(variant) -> Optional[str]:
    return _str(variant, "consequence")


# --- clinical ----------------------------------------------------------------
# ClinVar here is the VEP cache's CLIN_SIG, NOT a ClinVar VCF transfer. Two consequences
# the reader must hold onto:
#   1. NO review status. The cache has no CLNREVSTAT, so star ratings do not exist and
#      clinvar_stars()/the >=2-star auto-promote gate are gone. A 1-star single-submitter
#      P/LP assertion is now indistinguishable from an expert-panel one.
#   2. It is as stale as the cache (VEP 115 GRCh38 caches ClinVar 2025-02), whereas the
#      ClinVar VCF ships monthly. Reclassification is real; treat P/LP as a triage
#      prior, never as an answer.
# Both push toward false-positive retention (more to review), not toward missed calls.
def clnsig(variant) -> Optional[str]:
    return _str(variant, "clnsig")


def clnsig_is_plp(variant) -> bool:
    """True for a ClinVar P/LP assertion. VEP's CLIN_SIG is lowercase and '&'-joined
    (e.g. 'pathogenic&likely_pathogenic'); the ClinVar VCF's CLNSIG was Capitalised and
    '/'-or-','-joined. Both forms are matched so the predicate survives either source."""
    s = (clnsig(variant) or "").lower()
    if not s:
        return False
    # Treat conflicting as NOT P/LP even if the token appears in the conflict string.
    # VEP uses 'conflicting_interpretations_of_pathogenicity' /
    # 'conflicting_classifications_of_pathogenicity'; the substring covers both.
    if "conflicting" in s:
        return False
    return "pathogenic" in s and "likely_benign" not in s and "benign/likely" not in s


# --- GATK de novo tags (child-membership aware) ------------------------------
def hiconf_denovo_children(variant):
    """Set of child sample IDs listed in hiConfDeNovo, or None if the tag is absent."""
    s = _str(variant, "hiconf_denovo")
    if s is None:
        return None
    return {x.strip() for x in s.split(",") if x.strip()}


def is_hiconf_denovo_for(variant, child_id) -> bool:
    kids = hiconf_denovo_children(variant)
    return bool(kids) and child_id in kids
