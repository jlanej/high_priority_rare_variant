"""Step-3 plausibility classifier (pure; no VCF I/O so it is unit-testable).

Given the annotation getters and the config thresholds, decide whether a site is
biologically plausible and record WHY. Inheritance-agnostic (permissive-union rarity),
ClinVar P/LP as an override, BA1-common never rescued, gene lists NOT applied here
(never-drop rule). See docs/pipeline_design.md (Step 3).
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
    keep evidence ('clinvar_plp' | 'loftee_hc' | 'impact_high' | 'revel' | ...).
    """
    ba1 = _f(cfg, "filters.rarity.benign_ba1", 0.05)
    rec_max = _f(cfg, "filters.rarity.recessive_max", 1.0e-2)  # permissive-union cutoff
    revel_sup = _f(cfg, "filters.functional.revel_pp3_supporting", 0.644)
    am_lp = _f(cfg, "filters.functional.alphamissense_lp", 0.564)
    spliceai_pp3 = _f(cfg, "filters.functional.spliceai_pp3", 0.2)
    cadd_sup = _f(cfg, "filters.functional.cadd_phred_supporting", 25.3)  # Pejaver PP3-supporting
    mpc_strong = _f(cfg, "filters.functional.mpc_strong", 2.0)
    keep_impacts = set(get(cfg, "filters.functional.keep_impacts", ["HIGH", "MODERATE"]))
    loftee_require_hc = bool(get(cfg, "filters.functional.loftee_require_hc", True))
    clinvar_min_stars = int(get(cfg, "filters.clinvar.auto_promote_min_stars", 2))

    def functional_reason(v):
        if (A.impact(v) or "") in keep_impacts:
            return "impact_" + (A.impact(v) or "").lower()
        # Step 3 is a PERMISSIVE keep/drop screen: HIGH/MODERATE impact above already keeps all
        # pLoF regardless of LOFTEE confidence (LOFTEE only labels HIGH-impact consequences). The
        # LOFTEE HC/no-flags PVS1 *strength* is applied at the planned ACMG tiering step, not here;
        # this branch only catches a LoF call VEP did not mark HIGH/MODERATE.
        if loftee_require_hc and A.is_loftee_hc(v) and not A.loftee_flags(v):
            return "loftee_hc"
        for name, val, thr in (
            ("spliceai", A.spliceai_max(v), spliceai_pp3),
            ("revel", A.revel(v), revel_sup),
            ("alphamissense", A.alphamissense(v), am_lp),
            ("cadd", A.cadd(v), cadd_sup),
            ("mpc", A.mpc(v), mpc_strong),
        ):
            if val is not None and val >= thr:
                return name
        return None

    def classify(v):
        fr = A.frequency(v)
        if fr is not None and fr >= ba1:            # ClinGen BA1 — never rescue
            return False, "ba1"
        # ClinVar P/LP override auto-promotes only at >= auto_promote_min_stars (default 2)
        plp = A.clnsig_is_plp(v) and A.clinvar_stars(v) >= clinvar_min_stars
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
