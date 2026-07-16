"""Step-3 plausibility classifier (pure; no VCF I/O so it is unit-testable).

Given the annotation getters and the config thresholds, decide whether a site is
biologically plausible and record WHY. Inheritance-agnostic (permissive-union rarity),
ClinVar P/LP as an override, BA1-common never rescued, gene lists NOT applied here
(never-drop rule). See docs/pipeline_design.md (Step 3).

The functional ladder is deliberately two rungs — VEP IMPACT, then CADD. It used to try
spliceai -> revel -> alphamissense -> cadd -> mpc; under the VEP-only contract three of
those have no data source, and the missense pair were provably unreachable anyway (a
scored variant is missense => MODERATE => kept at the impact rung). Keeping them would
have meant an OR over correlated predictors that never fires — worse than useless,
because it reads as discriminative power the screen does not have.
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
    keep_impacts = set(get(cfg, "filters.functional.keep_impacts", ["HIGH", "MODERATE"]))

    def functional_reason(v):
        if (A.impact(v) or "") in keep_impacts:
            return "impact_" + (A.impact(v) or "").lower()
        # Everything below MODERATE reaches here, and CADD is the only score left that can
        # speak to it — the missense predictors (REVEL/AlphaMissense/MPC) could not, since a
        # variant carrying one is missense, hence MODERATE, hence already returned above.
        # That makes this the pipeline's ONLY keep-path for intronic / synonymous / UTR /
        # regulatory variants, so its threshold is the whole non-coding screen. See
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
