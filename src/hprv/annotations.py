"""Read the annotations produced by pipeline Step 2 and evaluate genotype QC.

**VEP-centric contract.** Every annotation this pipeline reads comes from ONE tool: VEP 115
(GRCh38) — its cache plus its score plugins — whose CSQ fields Step 2 lifts to INFO with a
``vep_`` prefix (via ``bcftools +split-vep``), with TWO documented exceptions. The plugins are CADD
(genome-wide functional, SNV+indel), SpliceAI (precomputed raw genome-wide splice deltas), and the
calibrated missense pair REVEL + AlphaMissense; population frequency and ClinVar ``CLIN_SIG`` ride
in the cache itself.

Both exceptions exist because the cache cannot supply the field AT ANY PRICE, and both are
``bcftools annotate`` transfers in Step 2 under their own INFO namespace — deliberately, so which
oracle a field came from is visible at a glance:

  * ``clinvar_*`` — ClinVar REVIEW STATUS. The cache carries ``CLIN_SIG`` but no ``CLNREVSTAT``,
    so gold stars require the ClinVar sites VCF itself.
  * ``gnomad_*``  — **faf95** and **nhomalt**, from the gnomAD v4.1 JOINT slim — REQUIRED under
    the default oracle (``resources.gnomad.oracle: faf95``), optional only under ``grpmax_proxy``.
    faf95's CI correction needs AC/AN, which the cache omits; see ``frequency()``.

No dbNSFP or LOFTEE file is transferred or read.

This module is the single place that knows those field names and how to coerce their (string)
values, so the selection, inheritance, and burden steps all read them identically.

What the contract still costs is documented in docs/allele_frequency.md and
docs/functional_annotation.md. ONE rarity oracle per run: real faf95 (the default; the joint slim
is required) or, by explicit opt-down, the grpmax point-estimate proxy — the arms never cross.
nhomalt travels with the slim. LOFTEE remains unwired either way. Two
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

from . import splice as _splice


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
    # --- transcript geometry (VEP --numbers / --total_length) and the stock NMD plugin ---
    # EXON/INTRON are "k/n" on the picked transcript; the three positions read "pos/len" because
    # Step 2 passes --total_length. NMD is Ensembl's own plugin (no data file): it writes
    # NMD_escaping_variant on a stop_gained / frameshift / canonical-splice variant that escapes
    # nonsense-mediated decay (last exon, within 50 nt of the penultimate exon's end, first 100
    # coding bases, intronless transcript) and NOTHING otherwise — see nmd_status() for why the
    # blank must be resolved where the header is.
    "exon": "vep_EXON",
    "intron": "vep_INTRON",
    "cdna_position": "vep_cDNA_position",
    "cds_position": "vep_CDS_position",
    "protein_position": "vep_Protein_position",
    "nmd": "vep_NMD",
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
    # --- ClinVar review status (one of TWO bcftools transfers; see 02_annotate_sites.sh) ---
    # The cache has no CLNREVSTAT at any price, so gold stars require the ClinVar VCF itself.
    # `clinvar_` (not `vep_`) marks it as transferred rather than lifted from the CSQ — a third
    # namespace deliberately, so a reader can tell at a glance which oracle a field came from.
    # Absent when the transfer did not run; clinvar_stars() then returns None (= UNAVAILABLE),
    # which is NOT the same as 0 stars ("no assertion criteria provided"). Never conflate them.
    "clnrevstat": "clinvar_CLNREVSTAT",
    # CLNSIG from the VCF too, so a run can compare the pinned release against the cache's
    # (possibly older) CLIN_SIG. Reporting only — clnsig_is_plp still reads the cache field.
    "clnsig_clinvar": "clinvar_CLNSIG",
    # --- gnomAD v4.1 JOINT slim (the SECOND bcftools transfer; see 02_annotate_sites.sh) ---
    # faf95 is the real filtering allele frequency: the LOWER bound of the 95% Poisson CI, which
    # is the quantity ACMG/ClinGen specify for frequency filtering (Whiffin 2017). The VEP cache
    # cannot supply it at any price — the CI correction needs AC/AN, which the cache omits — so
    # this is a transfer, prefixed `gnomad_` as a third-party namespace beside vep_/clinvar_.
    #
    # VERIFIED against the real v4.1 joint header + data: the FAF group set is
    # afr/amr/eas/mid/nfe/sas, i.e. GRPMAX_POPS plus `mid`, EXCLUDING the bottlenecked ami/asj/fin.
    # So it does not reintroduce the MAX_AF trap golden rule 2 forbids. `mid` is the single
    # deviation, which is why the producing group rides along and is reported per variant.
    "faf95": "gnomad_faf95",
    "faf95_group": "gnomad_faf95_group",
    # Homozygote count in gnomAD. The recessive false-positive tell: a "rare" homozygous call in a
    # gene where gnomAD already carries homozygotes is usually not the diagnosis.
    "nhomalt": "gnomad_nhomalt",
    # Global + grpmax joint AF: REPORTING ONLY, never a filter field (golden rule 2). Carried so a
    # reviewer can see the point estimate beside the CI-corrected value and judge the gap.
    "gnomad_af_joint": "gnomad_AF_joint",
    "gnomad_af_grpmax": "gnomad_AF_grpmax",
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

    A point estimate standing in for gnomAD's published grpmax AF, and the **opt-down** arm of
    ``frequency()`` (``resources.gnomad.oracle: grpmax_proxy``). The default arm is real
    ``faf95()`` from the gnomAD joint slim; the two are never mixed within a run, and this
    function is never consulted under the faf95 arm — not even when gnomAD published no faf95 for
    the allele (that resolves to 0 / ``zero_ci``, see ``frequency()``).

    It is a point estimate, not a CI lower bound: computing faf95 needs AC/AN, which the VEP
    cache does not carry, so THIS field cannot be CI-corrected no matter what. It therefore runs
    ~one CI-width HIGH on low-AC observations, and a rarity gate driven by it fires slightly more
    often than a faf95 gate would (i.e. it errs toward dropping). That is also precisely why it is
    the SAFE fallback: proxy >= faf95 wherever both exist, so falling back here can only filter
    more, never silently retain what faf95 would have caught. Excluding the bottlenecked groups
    removes the large half of the error; the residual is bounded by AC. See
    docs/allele_frequency.md.
    """
    return _max_float(variant, *_GRPMAX_KEYS)


def faf95(variant) -> Optional[float]:
    """gnomAD v4.1 joint **faf95** — the 95%-CI-corrected filtering allele frequency, or None.

    The quantity ACMG/ClinGen actually specify for frequency filtering (Whiffin 2017): the LOWER
    bound of the Poisson 95% CI on the population AF, maximised over the FAF-eligible genetic
    ancestry groups. Available only when the gnomAD joint slim is transferred in Step 2.

    **None is not "AF = 0" — and it is not "use the proxy" either.** gnomAD emits fafmax only
    where some group's CI lower bound is above zero; on a chr22 sample roughly 80% of records
    carried none, and where it was present it was always > 0. So None means "no group has a
    confidently non-zero frequency". ``frequency()`` resolves that to 0.0 (rarest,
    ``rarity_basis=zero_ci``) when gnomAD has the allele at all, and to None (absent, rarest)
    when it has no record — never to the point-estimate proxy, which would filter singletons on
    an inflated AF (96.5% of that class are AC <= 2).
    """
    return _max_float(variant, "faf95")


def faf95_group(variant) -> Optional[str]:
    """Which genetic ancestry group produced faf95 (afr/amr/eas/mid/nfe/sas), or None.

    Reported for the same reason `max_af_pops` rides beside `max_af`: the FAF group set is
    GRPMAX_POPS **plus `mid`**, so this is how a reviewer sees the one case where the rarity call
    rests on a group hprv's own proxy would have excluded.
    """
    return _str(variant, "faf95_group")


def nhomalt(variant) -> Optional[int]:
    """Homozygote count in gnomAD joint, or None if the transfer did not run.

    None and 0 are DIFFERENT: None = nobody looked, 0 = gnomAD has this allele and observed no
    homozygotes. A hom-recessive call in a gene where gnomAD already carries homozygotes is
    usually not the diagnosis; that check is only meaningful when the value is actually present.
    """
    v = _max_float(variant, "nhomalt")
    return None if v is None else int(v)


def rarity_oracle(cfg=None) -> str:
    """``faf95`` | ``grpmax_proxy`` — the ONE frequency oracle this run uses, for every variant.

    Deliberately a RUN-LEVEL constant, not a per-variant choice. An earlier design preferred faf95
    and fell back to the proxy per variant; that made two variants in one run comparable on
    different quantities and impossible to describe in a methods section. It also got the
    fallback direction wrong (see ``frequency()``). One oracle, chosen in config, recorded in the
    audit, and never crossed at runtime.

    ``resources.gnomad.oracle``, default **faf95** — the CORRECT quantity (a 95% CI lower bound,
    what ACMG/ClinGen specify), not the convenient one. It requires the gnomAD joint slim, and
    ``run_pipeline.sh`` HALTS at preflight when that is missing rather than quietly running on the
    other quantity — the same contract as ``resources.vep.spliceai_required``. Opt down to
    ``grpmax_proxy`` deliberately for a run without the slim.
    """
    from .config import get as _get
    v = str(_get(cfg or {}, "resources.gnomad.oracle", "faf95")).strip().lower()
    return "grpmax_proxy" if v == "grpmax_proxy" else "faf95"


def rarity_basis(variant, cfg=None) -> str:
    """How this variant's rarity value arose WITHIN the chosen oracle: ``measured`` |
    ``zero_ci`` | ``absent``.

    Provenance, not a second oracle — the quantity is the same for every row. ``zero_ci`` is the
    faf95-arm case where gnomAD HAS the allele but published no filtering AF, meaning no ancestry
    group's 95% CI lower bound clears zero (80% of a chr22 sample). That resolves to 0.0, and it
    is worth distinguishing from ``absent`` because the two are different facts: gnomAD looked and
    could not bound the frequency, versus gnomAD has no record at all.
    """
    if rarity_oracle(cfg) == "faf95":
        if faf95(variant) is not None:
            return "measured"
        return "zero_ci" if gnomad_observed(variant) else "absent"
    return "measured" if grpmax_af(variant) is not None else "absent"


def frequency(variant, cfg=None) -> Optional[float]:
    """The rarity value every gate reads — the single chokepoint for population frequency.

    **ONE oracle per run** (``rarity_oracle``), never a per-variant blend. Both arms are gnomAD
    v4.1; they differ in the QUANTITY and in how they reach it:

    * ``faf95`` — the published filtering allele frequency from the gnomAD JOINT sites slim: the
      lower bound of the 95% Poisson CI, the quantity ACMG/ClinGen specify for frequency
      filtering (Whiffin 2017). Absent-but-present-in-gnomAD resolves to **0.0**, because gnomAD
      emits fafmax as missing rather than as 0 wherever no group's CI clears zero; of that class
      96.5% are AC <= 2, and filtering a singleton on a point estimate is the error faf95 exists
      to prevent. Absent from the slim entirely -> ``None`` (rarest).
    * ``grpmax_proxy`` — a POINT ESTIMATE from the VEP cache: max AF over the grpmax-eligible
      ancestry groups. Same underlying dataset, no CI correction available (the cache carries no
      AC/AN), so it sits ~one CI-width high on low-count alleles and errs toward dropping.

    Neither arm ever consults the other. ``None`` means the chosen oracle has no value for this
    allele and every gate treats it as rarest. MAX_AF and global AFs are never consulted by
    either arm (golden rule 2).
    """
    if rarity_oracle(cfg) == "faf95":
        f = faf95(variant)
        if f is not None:
            return f
        return 0.0 if gnomad_observed(variant) else None
    return grpmax_af(variant)


def gnomad_observed(variant) -> bool:
    """True when the gnomAD joint slim carries a record for this allele.

    The witness that distinguishes "gnomAD computed a FAF of 0" from "gnomAD has never seen this
    variant". Both leave ``faf95`` absent, and they demand opposite fallbacks — see
    ``frequency()``. False whenever the slim was not transferred at all, which correctly puts the
    whole run on the proxy.
    """
    return _max_float(variant, "gnomad_af_joint") is not None


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


# --- transcript geometry, MANE, NMD ------------------------------------------
def exon(variant) -> Optional[str]:
    """VEP EXON as 'k/n' (exon number / exon count on the picked transcript), or None."""
    return _str(variant, "exon")


def intron(variant) -> Optional[str]:
    return _str(variant, "intron")


def cds_position(variant) -> Optional[str]:
    """VEP CDS_position — 'pos/len' under --total_length (Step 2 passes it) — or None."""
    return _str(variant, "cds_position")


def mane_select(variant) -> Optional[str]:
    """The MANE Select transcript ID when the picked transcript IS MANE Select, else None."""
    return _str(variant, "mane")


# The consequences Ensembl's NMD plugin grades (the four it evaluates the escape rules on).
NMD_CONSEQUENCES = frozenset({"stop_gained", "frameshift_variant",
                              "splice_donor_variant", "splice_acceptor_variant"})


def nmd_status(variant, plugin_present: bool) -> str:
    """``escaping`` | ``triggering`` | ``not_assessed`` — resolved ONCE here, where the header is.

    The NMD plugin writes ``NMD_escaping_variant`` when a covered consequence sits in the last
    exon, within 50 nt of the penultimate exon's end, in the first 100 coding bases, or in an
    intronless transcript (Abou Tayoun 2018's PVS1 escape rules) — and writes NOTHING otherwise.
    So a blank means two different things: "assessed, predicted to trigger NMD" (the plugin ran
    and the consequence is one it covers) or "nobody looked" (no plugin, or a consequence it never
    grades). A TSV cannot tell those apart later; this function can, because it sees the header.
    Step 9 consumes the verdict and never re-derives it — the same rule as rarity_basis.
    """
    if _str(variant, "nmd"):
        return "escaping"
    cq = (consequence(variant) or "").lower()
    covered = any(tok in NMD_CONSEQUENCES for tok in cq.replace("&", ",").split(","))
    return "triggering" if (plugin_present and covered) else "not_assessed"


# --- SpliceAI event decomposition (src/hprv/splice.py) -----------------------
def spliceai_symbol(variant) -> Optional[str]:
    """The gene SpliceAI scored against (its own GENCODE annotation), or None."""
    return _str(variant, "spliceai_symbol")


def _int_field(variant, key) -> Optional[int]:
    v = _raw(variant, key)
    if isinstance(v, (tuple, list)):
        v = v[0] if v else None
    if v is None or str(v).strip() in ("", "."):
        return None
    try:
        return int(float(str(v)))
    except (TypeError, ValueError):
        return None


def spliceai_components(variant):
    """(ds, dp) dicts keyed by splice.EVENTS, straight from the four DS_* / DP_* INFO fields."""
    ds = {e: _max_float(variant, "spliceai_ds_" + _splice.EVENT_SUFFIX[e].lower())
          for e in _splice.EVENTS}
    dp = {e: _int_field(variant, "spliceai_dp_" + _splice.EVENT_SUFFIX[e].lower())
          for e in _splice.EVENTS}
    return ds, dp


def spliceai_event(variant, floor: float = 0.2) -> dict:
    """splice.decompose() over the precomputed components: WHICH event, WHERE, and the effect.

    The floor is the screen's own spliceai_ds_min, so the same number that says a variant HAS a
    splice signal decides which components count as events (a max below it names no event)."""
    ds, dp = spliceai_components(variant)
    return _splice.decompose(ds, dp, pos=getattr(variant, "POS", None), floor=floor)
