"""Read the annotations produced by pipeline Step 2 and evaluate genotype QC.

**VEP-centric contract.** Every annotation this pipeline reads comes from ONE tool: VEP 115
(GRCh38) — its cache plus its score plugins — whose CSQ fields Step 2 lifts to INFO with a
``vep_`` prefix (via ``bcftools +split-vep``), with ONE documented exception. The plugins are CADD
(genome-wide functional, SNV+indel), SpliceAI (precomputed raw genome-wide splice deltas), and the
calibrated missense pair REVEL + AlphaMissense; population frequency and ClinVar ``CLIN_SIG`` ride
in the cache itself.

The exception is ClinVar REVIEW STATUS: the cache carries ``CLIN_SIG`` but no ``CLNREVSTAT`` at
any price, so gold stars require the ClinVar sites VCF itself, ``bcftools annotate``-transferred
in Step 2 under a ``clinvar_`` prefix. That third namespace is deliberate — it makes the different
oracle visible at a glance. No gnomAD, dbNSFP or LOFTEE file is transferred or read.

This module is the single place that knows those field names and how to coerce their (string)
values, so the selection, inheritance, and burden steps all read them identically.

What the contract still costs is documented in docs/allele_frequency.md and
docs/functional_annotation.md; in short: no faf95 (no CI correction), no nhomalt, no LOFTEE. Two
availability caveats travel with the scores: the PRECOMPUTED SpliceAI set does not cover every
indel (a missing score is NOT evidence of no effect; see spliceai_ds()), and REVEL/AlphaMissense
are missense-only, so ``None`` on any non-missense is expected rather than a gap. Adding another
annotation means adding its INFO field here AND its plugin/transfer in Step 2 — nothing else
reaches around this.

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
    # genome-wide, SNV+indel. Alongside SpliceAI it is one of only TWO keep-paths for
    # anything VEP rates below MODERATE (SpliceAI is checked first); REVEL/AlphaMissense
    # below are missense-only and so can never be a keep-path at all. CADD
    # v1.6+ ingests SpliceAI/MMSplice as input features, so it carries a lossy
    # re-encoding of the splice signal the SpliceAI plugin would have supplied.
    "cadd": "vep_CADD_PHRED",
    # --- calibrated MISSENSE predictors (REVEL + AlphaMissense plugins) ---
    # Missense-only, so they are INERT AT THE SCREEN by construction: every missense is
    # IMPACT=MODERATE and selection.py keeps it at the impact rung before any predictor is
    # consulted (docs/limitations.md #7). They exist for Step 9's missense tier, where they
    # replace an off-label CADD rank with a ClinGen-calibrated one. Adding them does not change
    # which variants are kept — only how the kept missense variants are ranked.
    # `am_pathogenicity`/`am_class` are the PLUGIN's key names; dbNSFP's name for the same
    # quantity is `AlphaMissense_score`, and using that would populate nothing.
    "revel": "vep_REVEL",
    "alphamissense": "vep_am_pathogenicity",
    "alphamissense_class": "vep_am_class",
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
    # The cache exposes CLIN_SIG ONLY — no review status. Stars come from the separate
    # ClinVar transfer below (clnrevstat). The old >=2-star auto-promote gate is NOT
    # reinstated: stars RANK in Step 9, and gating the screen would violate never-drop.
    # Values are lowercase, '&'-joined (e.g. "pathogenic&likely_pathogenic").
    "clnsig": "vep_CLIN_SIG",
    # --- ClinVar review status (the ONE bcftools transfer; see 02_annotate_sites.sh) ---
    # The cache has no CLNREVSTAT at any price, so gold stars require the ClinVar VCF itself.
    # `clinvar_` (not `vep_`) marks it as transferred rather than lifted from the CSQ — a third
    # namespace deliberately, so a reader can tell at a glance which oracle a field came from.
    # Absent when the transfer did not run; clinvar_stars() then returns None (= UNAVAILABLE),
    # which is NOT the same as 0 stars ("no assertion criteria provided"). Never conflate them.
    "clnrevstat": "clinvar_CLNREVSTAT",
    # CLNSIG from the VCF too, so a run can compare the pinned release against the cache's
    # (possibly older) CLIN_SIG. Reporting only — clnsig_is_plp still reads the cache field.
    "clnsig_clinvar": "clinvar_CLNSIG",
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
        # cyvcf2 returns a TUPLE (not a comma-joined string) for a numeric INFO field carrying more
        # than one value, and str(tuple) is "(0.001, 0.004)" — whose tokens '(0.001' / '0.004)'
        # BOTH fail float(), so the old single-str() loop returned None. `bcftools +split-vep` types
        # every vep_gnomAD{e,g}_<POP>_AF column Number=.,Type=Float, so any selector that emits more
        # than one CSQ block (-s all, -s mane at a two-gene locus, -s pick on a multiallelic site)
        # hit this — and a None frequency reads as "absent from gnomAD => rarest", RETAINING common
        # polymorphisms with no way to notice. Golden rule 2: this is the rarity oracle.
        # Iterate the container first, then split each element.
        for item in v if isinstance(v, (tuple, list)) else (v,):
            for tok in str(item).replace("&", ",").split(","):
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
# CADD (genome-wide), SpliceAI (splice), and the calibrated MISSENSE pair REVEL +
# AlphaMissense, all VEP plugins. MPC / MetaRNN remain unwired (they would need dbNSFP,
# whose pinned URL is dead upstream).
#
# The missense pair is INERT AT THE SCREEN by construction and always will be: they are
# missense-only scores, every missense is IMPACT=MODERATE, and selection.py keeps MODERATE
# at the impact rung and returns BEFORE any predictor is consulted. Their Step-3 keep-
# reasons are asserted never to fire in CI. Their only consumer is Step 9's missense tier,
# where they replace an off-label CADD rank with a ClinGen-calibrated one.
def cadd(variant) -> Optional[float]:
    return _max_float(variant, "cadd")


def revel(variant) -> Optional[float]:
    """REVEL score 0-1, or None. Missense-only — None on any non-missense is EXPECTED, not a gap.

    ClinGen SVI's calibrated thresholds (Pejaver 2022) are PP3 >= 0.644 supporting / 0.773
    moderate / 0.932 strong, BP4 <= 0.290 supporting / 0.183 moderate / 0.016 strong. hprv does
    not assign ACMG weight; Step 9 uses the same cut points to ORDER candidates.
    """
    return _max_float(variant, "revel")


def alphamissense(variant) -> Optional[float]:
    """AlphaMissense pathogenicity 0-1, or None. Missense-only, same caveat as revel().

    Note this reads the PLUGIN field (`am_pathogenicity`), not dbNSFP's `AlphaMissense_score`.
    """
    return _max_float(variant, "alphamissense")


def alphamissense_class(variant) -> Optional[str]:
    """AlphaMissense's own call: likely_benign / ambiguous / likely_pathogenic, or None."""
    return _str(variant, "alphamissense_class")


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
# TWO SOURCES, deliberately distinguishable by prefix:
#   * clnsig()  reads the VEP cache's CLIN_SIG (vep_ prefix). It is as stale as the cache
#     (VEP 115 caches ClinVar 2025-02) — reclassification is real, so treat P/LP as a
#     triage prior, never an answer.
#   * clinvar_stars() reads CLNREVSTAT from the ClinVar VCF TRANSFER (clinvar_ prefix),
#     which is also a fresher, independently version-pinned release.
# Stars RANK, they never gate: the screen is star-blind on purpose, because a keep/drop
# gate on review status would violate never-drop. A 1-star assertion is kept and reviewed,
# and ranked below a 3-star one by Step 9.
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


# ClinVar review status -> gold stars. Keys are CANONICALISED (see _canon_revstat): lowercase
# with every non-alphanumeric character stripped. That is deliberate, not lazy. The VCF value
# carries commas INSIDE it ("criteria_provided,_multiple_submitters,_no_conflicts") while the
# field is declared Number=., so the comma is also the array separator — cyvcf2 hands it back as
# a tuple, bcftools query joins it with commas, and a plain string compare against any one of
# those forms breaks on the others. Canonicalising collapses all of them to one key, and it also
# survives ClinVar's periodic renames of the separator style. The 2024 rename of
# "conflicting_interpretations" -> "conflicting_classifications" is carried as BOTH keys, because
# a pinned older release is still a legitimate input.
_REVSTAT_STARS = {
    "practiceguideline": 4,
    "reviewedbyexpertpanel": 3,
    "criteriaprovidedmultiplesubmittersnoconflicts": 2,
    "criteriaprovidedsinglesubmitter": 1,
    "criteriaprovidedconflictingclassifications": 1,
    "criteriaprovidedconflictinginterpretations": 1,      # pre-2024 spelling
    "noassertioncriteriaprovided": 0,
    "noclassificationprovided": 0,
    "noassertionprovided": 0,                             # pre-2024 spelling
    "noclassificationsfromunflaggedrecords": 0,
    "noassertionfromunflaggedrecords": 0,                 # pre-2024 spelling
    "noclassificationforthesinglevariant": 0,
    "noassertionforthesinglevariant": 0,                  # pre-2024 spelling
}


def _canon_revstat(value) -> str:
    """Collapse any CLNREVSTAT rendering to a single comparable key. See _REVSTAT_STARS."""
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        value = ",".join(str(v) for v in value)
    return "".join(ch for ch in str(value).lower() if ch.isalnum())


def clnrevstat(variant) -> Optional[str]:
    """Raw ClinVar review status, or None when the ClinVar transfer did not run."""
    return _str(variant, "clnrevstat")


def clinvar_stars(variant) -> Optional[int]:
    """ClinVar gold stars 0-4, or None when unavailable.

    **None and 0 are different facts and must never be merged.** None = the ClinVar transfer did
    not run (nobody looked); 0 = ClinVar has a record whose submitter provided no assertion
    criteria. Collapsing them would let an un-transferred run read as "every assertion is
    unreviewed", which silently damps every ClinVar-supported candidate in Step 9.

    An unrecognised status also returns None rather than 0 — a status string this table does not
    know is an unknown, not a zero-star assertion, and ClinVar has renamed these before.
    """
    return _REVSTAT_STARS.get(_canon_revstat(clnrevstat(variant)))


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
