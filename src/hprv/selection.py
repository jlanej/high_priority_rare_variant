"""Step-3 plausibility classifier (pure; no VCF I/O so it is unit-testable).

Given the annotation getters and the config thresholds, decide whether a site is
biologically plausible and record WHY. Inheritance-agnostic (permissive-union rarity),
ClinVar P/LP as an override, BA1-common never rescued, gene lists NOT applied here
(never-drop rule). See docs/pipeline_design.md (Step 3).

The functional ladder is three rungs — VEP IMPACT, then SpliceAI, then CADD (an OR: any one
keeps). SpliceAI is back because it now has a real data source (the precomputed raw genome-wide
splice scores, wired as a VEP plugin in Step 2): it is the ONLY signal that reaches a deep-intronic
cryptic splice site or an exonic-synonymous splice disruption, which both VEP's positional terms
and CADD under-call. The missense predictors (REVEL/AlphaMissense/MPC) stay OUT: they are
missense-only, and a scored missense is IMPACT=MODERATE => already kept at the impact rung, so
they never fire — an OR over correlated predictors that reads as discriminative power the screen
does not have. Each rung is keep-ONLY (a None/absent score never drops a variant).
"""

from __future__ import annotations

from hprv import annotations as A
from hprv.config import get


def _f(cfg, key, default):
    v = get(cfg, key, default)
    return float(v) if v is not None else default


def build_classifier(cfg):
    """Return classify(variant) -> (keep: bool, reason: str).

    reason is a drop reason ('ba1' | 'too_common' | 'not_functional') or the specific
    keep evidence ('clinvar_plp' | 'impact_high' | 'impact_moderate' | 'cadd').
    """
    ba1 = _f(cfg, "filters.rarity.benign_ba1", 0.05)
    rec_max = _f(cfg, "filters.rarity.recessive_max", 1.0e-2)  # permissive-union cutoff
    cadd_sup = _f(cfg, "filters.functional.cadd_phred_supporting", 25.3)
    sai_min = _f(cfg, "filters.functional.spliceai_ds_min", 0.2)
    keep_impacts = set(get(cfg, "filters.functional.keep_impacts", ["HIGH", "MODERATE"]))

    def functional_reason(v):
        if (A.impact(v) or "") in keep_impacts:
            return "impact_" + (A.impact(v) or "").lower()
        # SpliceAI: the specific splice-disruption signal. It is the ONLY predictor that reaches a
        # cryptic splice site deep in an intron, or an exonic-synonymous variant that breaks
        # splicing — cases VEP's positional splice terms miss entirely and CADD only weakly
        # re-encodes. Checked BEFORE CADD so a splice hit is labelled 'spliceai' (actionable),
        # not the generic 'cadd'. Keep-only; the raw delta score rides through for reviewer tiering.
        # ClinGen SVI uses >= 0.2 for PP3-supporting; a missing score never drops (see spliceai_ds).
        if (ds := A.spliceai_ds(v)) is not None and ds >= sai_min:
            return "spliceai"
        # CADD is the general functional score and the other keep-path below MODERATE impact — for
        # intronic / synonymous / UTR / regulatory variants SpliceAI does not flag. See
        # docs/functional_annotation.md for why 25.3 is a discovery rank here and not the
        # Pejaver PP3-supporting cutoff it is named after (that calibration is missense-only).
        if (val := A.cadd(v)) is not None and val >= cadd_sup:
            return "cadd"
        return None

    def classify(v):
        fr = A.frequency(v)
        if fr is not None and fr >= ba1:            # ClinGen BA1 — never rescue
            return False, "ba1"
        # ClinVar P/LP override. Previously gated on >= 2 review stars; the VEP cache carries
        # no review status, so an unstarred assertion is all we get and the gate is gone. This
        # admits 1-star single-submitter P/LP calls — i.e. it over-retains rather than
        # over-drops, which is the safe direction for a screen but adds curation load.
        plp = A.clnsig_is_plp(v)
        rarity_ok = (fr is None) or (fr < rec_max) or plp
        if not rarity_ok:
            return False, "too_common"
        if plp:
            return True, "clinvar_plp"
        fr_reason = functional_reason(v)
        if fr_reason:
            return True, fr_reason
        return False, "not_functional"

    return classify
