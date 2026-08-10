"""Step-9 prioritization logic (pure; no file I/O so it is unit-testable).

Step 6 nominates genes by recurrence. It cannot say whether a gene is producing more
candidate rows than its mutational target predicts, and nothing downstream re-ranks the
individual calls. This module supplies both, in three layers that compose:

  1. **Gene-level excess over mutational target.** ``E_g = C * mu_g`` from a gnomAD v2.1.1
     per-gene Samocha target, with a negative-binomial (NB2) null whose ``(C, alpha)`` are
     re-fit per cohort on an iteratively TRIMMED bulk, so the artifact tail does not
     calibrate its own null. Yields ``p_nb`` / ``q_nb`` (BH) plus a Poisson comparison and a
     mid-p calibration diagnostic — the self-check ``docs/gene_burden.md`` records as absent.
  2. **A six-signal artifact panel** collapsed to an integer ``corroboration_count``. The
     excess statistic says a gene is anomalous; the panel is the independent, mechanistically
     interpretable second witness required before any penalty is applied.
  3. **A per-variant tier + an additive Tavtigian-style points composite**, mechanism-gated:
     gene constraint counts only when the variant has a credible molecular effect, and a
     molecularly-benign prediction caps the total regardless of the gene.

Three contracts this module exists to honour, and every function below is written around them:

**Never-drop.** No layer here removes a row. A gene down-weight sets a tier, a separately
reported additive penalty, and a human-readable reason string; the maximum penalty (-3.0) is
deliberately too small to demote a variant carrying strong molecular evidence in a constrained
gene. It is a re-rank, not a veto. `09_prioritize.py` asserts row-count conservation.

**Absent evidence is never benign evidence.** A blank NHF means nobody looked, not that the
reads are clean (three states, never two — copied from
``Analysis/.../inherited/prepare_igv_variants.py:nhf_status``). A missing SpliceAI score means
the precomputed set does not cover that indel, not that splicing is unaffected — so the V0
"molecularly benign" tier requires BOTH scores to be PRESENT and below their cutoffs. A gene
with no gnomAD ``mu`` gets a CDS-length fallback offset AND a tier ceiling, because a ±30%
offset error cannot support a 5x claim.

**The excess statistic is a QUALITY question, not a biology question.** It answers *is this
gene producing more candidate rows than its mutational target predicts?* — never *is this gene
disease-associated?* High excess is evidence of a technical or population-genetic anomaly
(mismapping, paralogue collapse, a founder allele, a callability defect, an ancestry-uneven
rarity gate, a hypermutable locus). A gene can be both an established predisposition gene and
a mismapping hotspot; the correct response is read-level review, not discarding the gene.
The reason strings never assert biology, and the established-gene ceiling exists so the two
directions cannot be conflated. See docs/prioritization.md.

Numerics are pure ``math`` (regularized incomplete beta / gamma), not scipy: the tails reach
1e-30 where a ``1 - cdf`` sum would lose every significant digit, and ``tests/test_pure.py``
must run on a bare host. Every threshold arrives from the config; nothing here is law.
"""

from __future__ import annotations

import math
import re
from typing import Optional

from hprv.config import get

# --- tier / state vocabularies (single source of truth for the output columns) ------
GENE_TIERS = ("T0_no_downweight", "T1_watch", "T2_downweight", "T3_strong_downweight")
VARIANT_TIERS = ("V0", "V1", "V2", "V3", "V4", "V5")

# hprv Step-5 mode names for which pLoF constraint is NOT evidence. The canonical defaults say
# it directly: "Do not down-weight recessive candidates by pLoF constraint." Symmetrically, do
# not UP-weight them by it either — pLI/LOEUF measure selection against heterozygotes.
# `homozygous` is the spec's name for what Step 5 emits as `hom_recessive`; both are matched so
# the predicate survives either vocabulary.
RECESSIVE_MODES = frozenset({"compound_het", "hom_recessive", "homozygous", "x_linked_recessive"})

# Literal column-header tokens that must never be mistaken for a gene symbol. This is not
# theoretical fastidiousness: the validation cohort's own `gene.counts.txt` (a
# `cut -f6 | sort | uniq -c`) contains a literal `gene` line CARRYING n=1, which is why the true
# totals are 25,389 variants / 10,799 symbols rather than the 25,390 / 10,800 raw line count. A
# header row read as a gene is a phantom gene with a real count in every table it reaches.
GENE_KEYS_LOWER = frozenset({"gene", "gene_symbol", "symbol", "gene_name", "hgnc", "hgnc_symbol",
                             "#gene", "gene_id"})

# HIGH-impact terms that make a variant V5-ELIGIBLE (subject to the NMD gate, which is
# indeterminate today — see nmd_status()). Mirrors annotations.LOF_CONSEQUENCES minus the
# structural terms a short-read SNV/indel screen never emits.
PLOF_CONSEQUENCES = ("stop_gained", "frameshift_variant", "splice_acceptor_variant",
                     "splice_donor_variant", "start_lost")
CANONICAL_SPLICE = ("splice_acceptor_variant", "splice_donor_variant")

# Artifact-prone gene families (config: prioritization.signals.family.patterns). A CURATED,
# deliberately INCOMPLETE list, versioned as config rather than code. Family membership is
# never sufficient on its own — it enters only as one vote in the corroboration count, because
# olfactory receptors and zinc fingers are NOT significantly enriched among high-excess genes
# on the validation cohort (1.9x p=0.21 / 1.8x p=0.16) even though individual offenders are
# extreme. The mechanism (high paralogy, high pseudogene content) is what earns them a vote.
DEFAULT_FAMILY_PATTERNS = (
    r"^OR\d+[A-Z]*\d*$", r"^KRTAP", r"^MUC\d+", r"^PRAMEF?\d*", r"^NBPF\d+",
    r"^GOLGA\d", r"^NPIP", r"^HLA-", r"^TAS2R", r"^DEFB?\d",
    r"^LCE\d", r"^SPRR", r"^(H1-|H2A|H2B|H3-|H4-|HIST)",
)

# gnomAD's own constraint-model-failure flags. `no_exp_lof` is deliberately EXCLUDED: it is a
# coverage/annotation absence rather than a model failure, and it fires on genes where the
# family signal already catches the real problem. `mis_too_many` does the work here
# (5.1x enriched, p=6.2e-6); `syn_outlier` alone is 2.3x and not significant (p=0.074) but is
# kept for completeness.
DEFAULT_CONSTRAINT_FLAG_VALUES = ("mis_too_many", "syn_outlier")


# --- small coercions ---------------------------------------------------------
def _num(x) -> Optional[float]:
    """Float or None. A blank/'.'/'NA' cell is MISSING, never 0.0 — the distinction is
    load-bearing everywhere in this module (see the module docstring)."""
    if x is None:
        return None
    s = str(x).strip()
    if s in ("", ".", "NA", "na", "nan", "NaN", "None"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _int(x) -> Optional[int]:
    v = _num(x)
    return None if v is None else int(v)


def _s(x) -> str:
    return "" if x is None else str(x).strip()


def _f(cfg, key, default):
    v = get(cfg, key, default)
    return float(v) if v is not None else float(default)


def _i(cfg, key, default):
    v = get(cfg, key, default)
    return int(v) if v is not None else int(default)


# =============================================================================
# 1. Numerics — regularized incomplete beta / gamma, so the tails stay accurate
# =============================================================================
_EPS = 3.0e-16
_FPMIN = 1.0e-300


def _betacf(a: float, b: float, x: float) -> float:
    """Lentz continued fraction for the incomplete beta (Numerical Recipes betacf)."""
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < _FPMIN:
        d = _FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, 400):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < _FPMIN:
            d = _FPMIN
        c = 1.0 + aa / c
        if abs(c) < _FPMIN:
            c = _FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < _FPMIN:
            d = _FPMIN
        c = 1.0 + aa / c
        if abs(c) < _FPMIN:
            c = _FPMIN
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < _EPS:
            break
    return h


def betainc(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b) in [0, 1]."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(lbeta) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbeta) * _betacf(b, a, 1.0 - x) / b


def _gser(a: float, x: float) -> float:
    ap, s, d = a, 1.0 / a, 1.0 / a
    for _ in range(1000):
        ap += 1.0
        d *= x / ap
        s += d
        if abs(d) < abs(s) * _EPS:
            break
    return s * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gcf(a: float, x: float) -> float:
    """Upper regularized incomplete gamma Q(a, x) by continued fraction."""
    b = x + 1.0 - a
    c = 1.0 / _FPMIN
    d = 1.0 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < _FPMIN:
            d = _FPMIN
        c = b + an / c
        if abs(c) < _FPMIN:
            c = _FPMIN
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < _EPS:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def gammainc_lower(a: float, x: float) -> float:
    """Regularized lower incomplete gamma P(a, x) in [0, 1]."""
    if x <= 0.0:
        return 0.0
    if x < a + 1.0:
        return _gser(a, x)
    return 1.0 - _gcf(a, x)


def nb_logpmf(k: int, mu: float, alpha: float) -> float:
    """log P(N = k) for NB2: mean mu, variance mu*(1 + alpha*mu), size theta = 1/alpha."""
    theta = 1.0 / alpha
    p = theta / (theta + mu)
    return (math.lgamma(k + theta) - math.lgamma(theta) - math.lgamma(k + 1.0)
            + theta * math.log(p) + k * math.log1p(-p))


def nb_sf(n: int, mu: float, alpha: float) -> float:
    """P(N >= n) under NB2 — the primary excess p-value.

    Computed as I_{mu/(mu+theta)}(n, theta), NOT as 1 - sum(pmf): the tail reaches 1e-30 on
    the extreme loci and a complement sum would return 0.0 for every one of them, collapsing
    the whole artifact tail into a single indistinguishable bin.
    """
    if n <= 0:
        return 1.0
    if mu <= 0.0:
        return 0.0
    theta = 1.0 / alpha
    return betainc(float(n), theta, mu / (mu + theta))


def nb_midp(n: int, mu: float, alpha: float) -> float:
    """Mid-p: P(N > n) + 0.5*P(N = n). The CALIBRATION statistic, not the screening p.

    A discrete upper tail is systematically conservative at any threshold, which makes a raw
    P(N >= n) useless for asking "is this null calibrated?". Mid-p is uniform-on-average under
    a correct null, so mean(mid-p) ~ 0.5 and P(mid-p < a) ~ a are honest diagnostics.
    """
    if mu <= 0.0:
        return 1.0 if n <= 0 else 0.0
    return nb_sf(n + 1, mu, alpha) + 0.5 * math.exp(nb_logpmf(n, mu, alpha))


def pois_sf(n: int, mu: float) -> float:
    """P(N >= n) under Poisson(mu), reported for transparency only (see calibrate_null)."""
    if n <= 0:
        return 1.0
    if mu <= 0.0:
        return 0.0
    return gammainc_lower(float(n), mu)


def pois_midp(n: int, mu: float) -> float:
    if mu <= 0.0:
        return 1.0 if n <= 0 else 0.0
    logpmf = -mu + n * math.log(mu) - math.lgamma(n + 1.0)
    return pois_sf(n + 1, mu) + 0.5 * math.exp(logpmf)


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values; None passes through as None.

    Same contract and step-down algorithm as ``06_gene_burden.bh_fdr`` — deliberately
    duplicated rather than imported, because that one lives in a pipeline script. Keep them
    in step: a divergence here would make Step 6 and Step 9 disagree about the same p-vector.
    """
    idx = [i for i, p in enumerate(pvals) if p is not None]
    m = len(idx)
    q = [None] * len(pvals)
    prev = 1.0
    for rank, i in enumerate(reversed(sorted(idx, key=lambda i: pvals[i])), start=1):
        prev = min(prev, pvals[i] * m / (m - rank + 1))
        q[i] = prev
    return q


def percentile(values, frac: float) -> Optional[float]:
    """Linear-interpolated percentile of a list of numbers (frac in [0, 1]); None if empty.

    Used for the ``classic_caf`` bottom-decile cut, which is a COHORT property (the 10th
    percentile of the candidate genes in THIS run), not a constant.
    """
    xs = sorted(v for v in values if v is not None)
    if not xs:
        return None
    if len(xs) == 1:
        return xs[0]
    pos = frac * (len(xs) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (pos - lo)


# =============================================================================
# 2. The offset: mu_g and E_g
# =============================================================================
def mutational_target(mu_mis, mu_syn, mu_lof, impute_factor: float = 0.0516):
    """-> (mu_tot, mu_lof_src). ``mu_lof`` is null for ~505 of 19,643 gnomAD v2.1.1 genes.

    Those genes are neither dropped nor charged ``mu_lof = 0`` (which would inflate their
    excess ratio by ~5% for no reason). They are imputed at the MEDIAN
    ``mu_lof / (mu_mis + mu_syn)`` and the imputation is made visible in ``mu_lof_src``.
    """
    mis, syn, lof = _num(mu_mis), _num(mu_syn), _num(mu_lof)
    if mis is None and syn is None and lof is None:
        return None, "none"
    base = (mis or 0.0) + (syn or 0.0)
    if lof is None:
        if base <= 0.0:
            return None, "none"
        return base * (1.0 + impute_factor), "imputed"
    return base + lof, "gnomad"


def mu_from_cds(cds_length, b0: float = -18.4484, b1: float = 1.0570) -> Optional[float]:
    """CDS-length fallback offset: log(mu_tot) = b0 + b1*log(cds_length).

    OLS over the full gnomAD gene universe (R^2 = 0.828); the 10th-90th percentile of
    predicted/actual spans 0.712-1.283, i.e. roughly a +/-30% band. A gene offset this coarse
    cannot support a 5x claim, which is exactly why ``assign_gene_tier`` refuses to let a
    CDS-fallback gene reach the strong tier.
    """
    L = _num(cds_length)
    if L is None or L <= 0.0:
        return None
    return math.exp(b0 + b1 * math.log(L))


def fit_scaling(counts, mus) -> Optional[float]:
    """C = sum(n) / sum(mu) over the genes handed in.

    **Hand this the FULL gene universe, including the zero-count genes.** The candidate list is
    a zero-truncated sample: on the validation cohort 9,266 of 19,643 genes have n = 0, and
    those zeros carry real information about C. Fitting C on matched genes only inflates C and
    deflates every excess ratio — the direction that hides artifact loci.
    """
    sm = sum(m for m in mus if m)
    if sm <= 0.0:
        return None
    return sum(n or 0 for n in counts) / sm


def nb_alpha_mom(counts, means) -> Optional[float]:
    """Method-of-moments NB2 dispersion: alpha = (sum (n-mu)^2 - sum mu) / sum mu^2.

    The MLE's starting point, and the fallback when the likelihood is flat. Can go negative on
    under-dispersed data, in which case there is no NB2 fit to make and the caller should say
    so rather than clamp silently.
    """
    num = sum((float(n) - m) ** 2 for n, m in zip(counts, means)) - sum(means)
    den = sum(m * m for m in means)
    if den <= 0.0:
        return None
    return num / den


def nb_alpha_mle(counts, means, lo: float = 1e-4, hi: float = 20.0, iters: int = 200):
    """Golden-section maximisation of the NB2 log-likelihood in alpha (dispersion fixed mean).

    Bounded 1-D search rather than Newton: the profile likelihood is unimodal in alpha here,
    and a bracketed search cannot walk off to a negative dispersion the way a Newton step can.
    """
    def ll(a):
        return sum(nb_logpmf(int(n), m, a) for n, m in zip(counts, means) if m > 0.0)

    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - inv_phi * (b - a), a + inv_phi * (b - a)
    fc, fd = ll(c), ll(d)
    for _ in range(iters):
        if b - a < 1e-9:
            break
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - inv_phi * (b - a)
            fc = ll(c)
        else:
            a, c, fc = c, d, fd
            d = a + inv_phi * (b - a)
            fd = ll(d)
    return 0.5 * (a + b)


def pearson_dispersion(counts, means) -> Optional[float]:
    """Pearson chi-square / df about the Poisson mean. phi = 1 is Poisson."""
    pairs = [(float(n), m) for n, m in zip(counts, means) if m and m > 0.0]
    if len(pairs) <= 1:
        return None
    chi2 = sum((n - m) ** 2 / m for n, m in pairs)
    return chi2 / (len(pairs) - 1)


def fit_excess_null(counts, mus, trim_p: float = 1.0e-3, trim_max_fraction: float = 0.05,
                    max_iter: int = 20, alpha_method: str = "mle") -> dict:
    """Iteratively-trimmed NB2 null fit over the FULL gene universe.

    The tail must not calibrate its own null. Untrimmed, the artifact loci are absorbed into
    the dispersion and then hide themselves: on the validation cohort the global fit gives
    alpha = 0.804 against a trimmed alpha = 0.214, nearly 4x. So:

      1. keep = all genes
      2. C     = sum(n[keep]) / sum(mu[keep])
      3. alpha = NB dispersion on keep at mean C*mu[keep]
      4. keep  = {g : P(N >= n_g) > trim_p}
      5. repeat until keep is stable

    Raises ValueError if the loop trims more than ``trim_max_fraction`` of genes. That means
    the null is mis-specified for this cohort, not that 5% of the exome is artifact — and
    proceeding would mass-down-weight real genes.
    """
    n_all = len(counts)
    if n_all == 0 or len(mus) != n_all:
        raise ValueError("fit_excess_null: counts and mus must be equal-length and non-empty")
    keep = [i for i, m in enumerate(mus) if m and m > 0.0]
    if not keep:
        raise ValueError("fit_excess_null: no gene has a positive mutational target")

    C = alpha = None
    iterations = 0
    prev_keep = None
    for iterations in range(1, max_iter + 1):
        kc = [counts[i] or 0 for i in keep]
        km = [mus[i] for i in keep]
        C = fit_scaling(kc, km)
        if not C or C <= 0.0:
            raise ValueError("fit_excess_null: scaling constant C is not positive")
        means = [C * m for m in km]
        a_mom = nb_alpha_mom(kc, means)
        if alpha_method == "mom":
            alpha = a_mom
        else:
            alpha = nb_alpha_mle(kc, means)
        if alpha is None or alpha <= 0.0:
            raise ValueError(
                "fit_excess_null: NB dispersion is not positive (the counts are under-dispersed "
                "relative to Poisson) — use null_model: poisson for this cohort rather than "
                "forcing a negative-binomial fit")
        new_keep = [i for i, m in enumerate(mus)
                    if m and m > 0.0 and nb_sf(int(counts[i] or 0), C * m, alpha) > trim_p]
        if not new_keep:
            raise ValueError("fit_excess_null: the trim loop removed every gene")
        if new_keep == prev_keep or new_keep == keep:
            keep = new_keep
            break
        prev_keep, keep = keep, new_keep

    n_eligible = sum(1 for m in mus if m and m > 0.0)
    n_trimmed = n_eligible - len(keep)
    frac = n_trimmed / float(n_eligible) if n_eligible else 0.0
    if frac > trim_max_fraction:
        raise ValueError(
            f"fit_excess_null: the trim loop removed {n_trimmed}/{n_eligible} genes "
            f"({frac:.1%}) > trim_max_fraction={trim_max_fraction:.1%}. That is a "
            "mis-specified null, not an exome that is 5% artifact — check the mutational-target "
            "table's gene keys and that the candidate counts cover the same universe.")

    kc = [counts[i] or 0 for i in keep]
    km = [mus[i] for i in keep]
    # phi_all is the RAW overdispersion diagnostic and must be measured about the UNTRIMMED
    # POISSON mean — that is the quantity "Poisson requires phi = 1" refers to. Measuring it
    # about the trimmed-NB scaling instead conflates two questions and yields a number that
    # cannot be compared to 1. phi_bulk (about the trimmed fit) is its companion: the pair
    # `phi_all >> 1` with `phi_bulk ~ 1` is the signature of overdispersion concentrated in a
    # tiny tail rather than spread through the exome, which is what justifies the NB choice AND
    # the trim. On the validation cohort: 18.7 raw, 1.29 trimmed.
    all_counts = [c or 0 for c, m in zip(counts, mus) if m and m > 0.0]
    all_mus = [m for m in mus if m and m > 0.0]
    C_pois = fit_scaling(all_counts, all_mus)
    return {
        "C": C, "alpha": alpha, "theta": 1.0 / alpha, "alpha_mom": nb_alpha_mom(kc, [C * m for m in km]),
        "iterations": iterations, "n_trimmed": n_trimmed, "trim_fraction": frac,
        "n_genes_in_fit": len(keep), "n_genes_eligible": n_eligible,
        "used_in_null_fit": set(keep), "C_poisson": C_pois,
        "phi_bulk": pearson_dispersion(kc, [C * m for m in km]),
        "phi_all": (pearson_dispersion(all_counts, [C_pois * m for m in all_mus])
                    if C_pois else None),
    }


def fit_poisson_null(counts, mus, trim_p: float = 1.0e-3, max_iter: int = 20,
                     trim_max_fraction: float = 0.05) -> dict:
    """Iteratively-trimmed POISSON scaling constant — the Poisson arm's own fit.

    Identical fixed-point loop to ``fit_excess_null`` with the NB tail replaced by the Poisson
    tail, and it exists for one reason: **the trim loop applies to whichever null is being fit.**

    A side-by-side NB-vs-Poisson calibration is only a fair comparison if BOTH arms are fit the
    way they would be used. Evaluating the Poisson at the NB's trimmed ``C`` instead measures a
    hybrid nobody would deploy, and it changes the headline number materially: the published
    2.51x Poisson anti-conservatism at alpha = 1e-3 is reproducible only WITH the Poisson arm
    trimmed. Under a plain untrimmed reading it is 1.82x (range 1.60-2.04x over 40 random
    splits). **Never quote the 2.51x figure without that condition.** The design decision is
    robust either way — every variant tested leaves the Poisson anti-conservative (all >= 1.31x)
    and the NB conservative (0.31x) — but the two numbers are not interchangeable.
    """
    keep = [i for i, m in enumerate(mus) if m and m > 0.0]
    n_eligible = len(keep)
    if not keep:
        raise ValueError("fit_poisson_null: no gene has a positive mutational target")
    C = fit_scaling([counts[i] or 0 for i in keep], [mus[i] for i in keep])
    prev, iterations = None, 0
    for iterations in range(1, max_iter + 1):
        C = fit_scaling([counts[i] or 0 for i in keep], [mus[i] for i in keep])
        nxt = [i for i in range(len(mus))
               if mus[i] and mus[i] > 0.0 and pois_sf(counts[i] or 0, C * mus[i]) > trim_p]
        if not nxt:
            raise ValueError("fit_poisson_null: the trim removed every gene — mis-specified null")
        if nxt == keep or nxt == prev:
            keep = nxt
            break
        prev, keep = keep, nxt
    C = fit_scaling([counts[i] or 0 for i in keep], [mus[i] for i in keep])
    n_trimmed = n_eligible - len(keep)
    frac = n_trimmed / n_eligible if n_eligible else 0.0
    if frac > trim_max_fraction:
        raise ValueError(
            f"fit_poisson_null: the trim removed {n_trimmed}/{n_eligible} genes ({frac:.2%}) > "
            f"trim_max_fraction {trim_max_fraction:.2%} — the null is mis-specified for this "
            "cohort (usually the offset table and the counts describe different gene universes)")
    return {"C": C, "iterations": iterations, "n_trimmed": n_trimmed, "trim_fraction": frac,
            "n_genes_in_fit": len(keep), "n_genes_eligible": n_eligible}


def calibrate_null(counts, mus, C: float, alpha: float, ratio_max: float = 5.0,
                   C_poisson=None, trim_p: float = 1.0e-3) -> dict:
    """Mid-p calibration diagnostic for the fitted null — the A-3 self-check.

    Evaluated on the BULK (``excess_ratio < ratio_max``), a *shape* criterion chosen so the
    evaluation set is not defined by the p-value being tested. Under a correct null,
    mean(mid-p) ~ 0.5, P(mid-p < 0.05) ~ 0.05, P(mid-p < 0.001) ~ 0.001. The NB and Poisson
    numbers are reported side by side because the choice between them is exactly what this
    diagnostic decides, and the direction of error matters: an anti-conservative null nominates
    artifact genes that are not there (cost: a penalised real variant), a conservative one
    under-calls (cost: reviewer time). Under the never-drop rule, conservative is correct.

    **The Poisson arm is evaluated at its OWN iteratively-trimmed scaling constant**, not at the
    NB's — the trim loop belongs to whichever null is being fit, and evaluating the Poisson at
    the NB's ``C`` measures a hybrid nobody would deploy. ``C_poisson`` is computed by
    ``fit_poisson_null`` when not supplied; the value used is reported as ``C_poisson`` so the
    comparison is auditable. Both `ratio_max`-based evaluation sets use each arm's own ``C``, so
    the two arms are judged on the bulk each of them defines.
    """
    out = {"n_evaluated": 0}
    if C_poisson is None:
        try:
            C_poisson = fit_poisson_null(counts, mus, trim_p=trim_p)["C"]
        except ValueError:
            # A Poisson arm that cannot converge is itself informative, and it must not take the
            # NB diagnostic down with it: fall back to the untrimmed Poisson C and say so.
            C_poisson = fit_scaling([c or 0 for c, m in zip(counts, mus) if m and m > 0.0],
                                    [m for m in mus if m and m > 0.0])
            out["poisson_trim_failed"] = True
    nbp, pop = [], []
    for n, m in zip(counts, mus):
        if not m or m <= 0.0:
            continue
        k = int(n or 0)
        E = C * m
        if E > 0.0 and (k / E) < ratio_max:
            nbp.append(nb_midp(k, E, alpha))
        Ep = C_poisson * m
        if Ep > 0.0 and (k / Ep) < ratio_max:
            pop.append(pois_midp(k, Ep))
    if not nbp:
        return out
    out["n_evaluated"] = len(nbp)
    out["n_evaluated_poisson"] = len(pop)
    out["C_poisson"] = C_poisson
    for label, ps in (("nb", nbp), ("poisson", pop)):
        if not ps:
            continue
        out[f"{label}_mean_midp"] = sum(ps) / len(ps)
        out[f"{label}_frac_lt_0.05"] = sum(1 for p in ps if p < 0.05) / len(ps)
        out[f"{label}_frac_lt_0.001"] = sum(1 for p in ps if p < 0.001) / len(ps)
    return out


# =============================================================================
# 3. The six-signal artifact panel
# =============================================================================
def family_matcher(patterns=None):
    """-> match(symbol) -> matched pattern or ''. Patterns come from config, not code."""
    pats = [re.compile(p) for p in (patterns or DEFAULT_FAMILY_PATTERNS)]
    def match(symbol: str) -> str:
        s = _s(symbol)
        for p in pats:
            if p.search(s):
                return p.pattern
        return ""
    return match


def artifact_signals(row, thresholds: dict) -> dict:
    """The six orthogonal artifact signatures for one gene -> booleans + corroboration count.

    Every flag travels with the VALUE that set it, so the reason string is self-documenting:
    ``corrob=3 (segdup + oe_syn + saturation)`` is actionable, ``0.71`` is not. Unweighted sum,
    deliberately — the folds span 4x-195x, but weights fit on the same tail they are used to
    judge would not transfer to another cohort.

    ``thresholds`` keys: ``per_trio_min``, ``segdup_min_frac``, ``oe_syn_max_deviation``,
    ``constraint_flag_values``, ``family_match`` (callable), ``caf_low_cutoff``.
    """
    reasons = []
    sig = {}

    # 1. Cohort saturation. The strongest signal in the panel (195x) and the only one needing
    #    no mutation model, no external annotation and no gnomAD join — hence the most portable
    #    across cohorts and the least likely to fail silently. It is NEVER sufficient on its
    #    own: TTN/SYNE1/DNAH11 clear it at 1.9-3.7x excess.
    per_trio = _num(row.get("per_trio"))
    sig["sig_saturation"] = per_trio is not None and per_trio >= thresholds["per_trio_min"]
    if sig["sig_saturation"]:
        reasons.append(f"cohort_saturation={per_trio:.4g}_variants_per_trio")

    # 2. Segmental duplication. Coordinate build is load-bearing: gnomAD v2.1.1 constraint
    #    start/end are GRCh37, so the overlap must be computed against a track on the SAME
    #    build. Mixing builds silently produces ~0 overlap, i.e. the signal goes dark with no
    #    error. The fraction is emitted on EVERY gene regardless of tier, because it has a
    #    second, independent use: marking genes where a NEGATIVE result is uninformative
    #    (CYP21A2/SMN1/SMN2/NCF1/STRC all have zero candidates against non-trivial E).
    segdup = _num(row.get("segdup98_frac"))
    sig["sig_segdup"] = segdup is not None and segdup >= thresholds["segdup_min_frac"]
    if sig["sig_segdup"]:
        reasons.append(f"segdup98_frac={segdup:.3g}")

    # 3. Artifact-prone gene family (one vote only; see DEFAULT_FAMILY_PATTERNS).
    fam = thresholds["family_match"](row.get("gene"))
    sig["sig_family"] = bool(fam)
    if fam:
        reasons.append("artifact_gene_family")

    # 4. Synonymous o/e departure. Synonymous variation is the closest thing gnomAD offers to a
    #    neutral internal control, so oe_syn ~ 1 wherever the mutation model and the callability
    #    are both sound. Departure in EITHER direction is informative and they mean different
    #    things: < 0.7 = gnomAD under-called the gene, so mu_g (our denominator) is unreliable;
    #    > 1.3 = paralogue collapse inflating gnomAD's own counts, the same mechanism that
    #    inflates ours. Emit the raw value so the direction stays visible.
    oe_syn = _num(row.get("oe_syn"))
    dev = thresholds["oe_syn_max_deviation"]
    sig["sig_oe_syn"] = oe_syn is not None and abs(oe_syn - 1.0) > dev
    if sig["sig_oe_syn"]:
        reasons.append(f"oe_syn={oe_syn:.3g}")

    # 5. gnomAD's own model-failure flag — gnomAD telling you ITS constraint model did not fit.
    flag = _s(row.get("constraint_flag")).lower()
    hits = [v for v in thresholds["constraint_flag_values"] if v and v.lower() in flag]
    sig["sig_constraint_flag"] = bool(hits)
    if hits:
        reasons.append("gnomad_constraint_flag=" + ",".join(hits))

    # 6. Low cumulative pLoF allele frequency. The direction is counter-intuitive and must be
    #    stated: high-excess genes have FAR LOWER classic_caf, not higher (median 2.1e-5 vs
    #    1.3e-4), and the high decile is actually depleted — so this is NOT "common variants
    #    leaking through the rarity gate". A caf near zero means gnomAD reports essentially no
    #    pLoF alleles, which happens either because the gene is small/constrained or because
    #    gnomAD's own callability there is poor (median exp_syn 27.8 vs 126.7 disambiguates it:
    #    these genes are small or poorly covered). It is a LOW-INFORMATION-LOCUS flag, not a
    #    frequency flag — and the weakest mechanism in the panel, so the strong tier never
    #    rests on it alone (see assign_gene_tier).
    #
    # NULL HANDLING is an explicit choice, not an accident, because the two readings differ and
    # one is unsound. The cutoff is computed over the genes that HAVE a value (see the caller);
    # a gene with a NULL classic_caf is then FLAGGED when `include_null_as_flagged` — gnomAD
    # reporting no pLoF allele frequency at all is precisely the low-information case. The
    # alternative (fill nulls with 0, THEN take the percentile) makes the threshold applied to the
    # MEASURED data depend on how many genes gnomAD has no record for: on the validation cohort
    # the ~0.9% of genes with a null classic_caf (93 of 10,377) shift the cut from 7.96e-6 to
    # 4.17e-6 and 84 genes change flag. The tier table is identical either
    # way; the corroboration distribution is not.
    cutoff = thresholds.get("caf_low_cutoff")
    null_flagged = thresholds.get("caf_null_is_flagged", True)
    caf = _num(row.get("classic_caf"))
    if caf is None:
        sig["sig_caf_low"] = bool(null_flagged) and cutoff is not None
        if sig["sig_caf_low"]:
            reasons.append("classic_caf=absent_low_information_locus")
    else:
        sig["sig_caf_low"] = cutoff is not None and caf <= cutoff
        if sig["sig_caf_low"]:
            reasons.append("classic_caf=0_bottom_decile" if caf == 0.0
                           else f"classic_caf={caf:.3g}_bottom_decile")

    sig["corroboration_count"] = sum(1 for k, v in sig.items() if k.startswith("sig_") and v)
    sig["signal_reasons"] = reasons
    return sig


def signal_thresholds(cfg, caf_low_cutoff=None) -> dict:
    """Assemble the artifact_signals() threshold bundle from config."""
    return {
        "per_trio_min": _f(cfg, "prioritization.signals.saturation.per_trio_min", 0.10),
        "segdup_min_frac": _f(cfg, "prioritization.signals.segdup.min_frac", 0.10),
        "oe_syn_max_deviation": _f(cfg, "prioritization.signals.oe_syn.max_deviation", 0.30),
        "constraint_flag_values": list(get(cfg, "prioritization.signals.constraint_flag.values",
                                          list(DEFAULT_CONSTRAINT_FLAG_VALUES))),
        "family_match": family_matcher(get(cfg, "prioritization.signals.family.patterns",
                                           list(DEFAULT_FAMILY_PATTERNS))),
        "caf_low_cutoff": caf_low_cutoff,
        "caf_null_is_flagged": bool(get(cfg or {}, "prioritization.signals.caf_low"
                                        ".include_null_as_flagged", True)),
    }


# =============================================================================
# 4. The four-tier gene down-weight
# =============================================================================
def assign_gene_tier(excess_ratio, q_nb, corroboration_count, saturation, n_observed,
                     established_gene_control=False, e_source="gnomad_mu", cfg=None) -> dict:
    """-> {gene_tier, gene_artifact_penalty, control_ceiling_applied, cds_ceiling_applied,
    tier_rule, review_flag}.

    Graded, four-level, and NEVER a hard drop. Assignment runs in ascending severity so the
    highest matching tier wins, then the two ceilings are applied last.

    Rule shapes, each for a measured reason:
      * **T1 requires no corroboration** — it carries no real cost (-0.5) and its job is to
        populate a reviewable watch list, including the "unexplained excess" class. q < 0.25 is
        deliberately loose.
      * **T2 uses the ratio, not the FDR** — at the moderate-excess boundary the FDR is
        count-driven (no gene with n < 5 can reach q < 0.05), so an FDR-only rule would exempt
        every low-count artifact gene. Requiring one corroborating signal supplies the evidence
        the count cannot.
      * **T3 has three disjuncts** because there are three distinct ways to be confidently
        anomalous: significant with a mechanism, extreme with two mechanisms, or saturating the
        cohort with two mechanisms. The third exists because saturation is offset-free and
        therefore survives an unreliable mu_g.
      * **Established-gene ceiling** — a control-set gene caps at T1_watch. Both this and a
        stricter global threshold reach 100% control sensitivity; the exemption is preferable
        because it is AUDITABLE and REVERSIBLE (the exempted genes are named, keep their full
        statistics, and a reviewer can disagree), whereas a stricter threshold would silently
        release ~190 genes back into the review queue with no record of the trade. **It
        protects them from a score penalty; it does not mean their calls are correct** — hence
        ``review_flag = established_gene_high_excess``.
      * **CDS-fallback ceiling** — a gene whose offset came from the CDS regression never
        reaches T3, because a +/-30% offset error cannot support that claim.
    """
    cfg = cfg or {}
    pfx = "prioritization.gene_downweight"
    t1_ratio = _f(cfg, f"{pfx}.t1_watch.min_ratio", 3.0)
    t1_q = _f(cfg, f"{pfx}.t1_watch.max_q", 0.25)
    t2_ratio = _f(cfg, f"{pfx}.t2_downweight.min_ratio", 5.0)
    t2_corrob = _i(cfg, f"{pfx}.t2_downweight.min_corroboration", 1)
    t3_q = _f(cfg, f"{pfx}.t3_strong.max_q", 0.05)
    t3_q_corrob = _i(cfg, f"{pfx}.t3_strong.q_min_corroboration", 1)
    t3_ratio = _f(cfg, f"{pfx}.t3_strong.min_ratio", 10.0)
    t3_ratio_corrob = _i(cfg, f"{pfx}.t3_strong.ratio_min_corroboration", 2)
    t3_sat_ratio = _f(cfg, f"{pfx}.t3_strong.saturation_min_ratio", 5.0)
    t3_sat_corrob = _i(cfg, f"{pfx}.t3_strong.saturation_min_corroboration", 2)
    # Count floor behind a RATIO rule. Default 3 — the value `excess_statistic_spec` §4.4
    # specifies and then never applies to the T2 rule, which independent validation identified as
    # a substantive defect. Both tier tables are measured on the same real data:
    #
    #   floor=1 (spec's tier table): 228 genes / 2,565 variants (10.36%), 13 control genes need
    #           the ceiling. But 136 of those 228 genes (59.6%) have n < 3 and carry only 196
    #           variants between them (7.6% of triaged volume) — most of the down-weighted gene
    #           LIST is singletons and doubletons contributing almost none of the benefit while
    #           generating most of the exemption burden. Ten of the 13 ceiling-exempted genes have
    #           n <= 4 and EIGHT have q_nb = 1 (no statistical evidence of excess at all): SMPX,
    #           HMGA2, HAMP and PET100 are each ONE variant against E ~ 0.16. They never needed
    #           protecting — they needed a minimum-count requirement.
    #   floor=3 (DEFAULT): 92 genes / 2,369 variants (9.57%) with the ceiling — 92.4% of the
    #           triaged VOLUME using 92 genes instead of 228 — and the exemption burden falls from
    #           13 control genes to 8. Retention stays 100%.
    #
    # This guards only the RATIO rules; the FDR limb of T3 needs no floor because no gene with
    # n < 5 reached q_nb < 0.05 on the validation cohort, so the NB tail already handles low
    # counts correctly on its own. Set to 1 to reproduce the original spec tier table exactly.
    min_n_ratio = _i(cfg, "prioritization.excess.min_n_for_ratio_rule", 3)
    penalties = get(cfg, f"{pfx}.penalties", None) or {
        "T0_no_downweight": 0.0, "T1_watch": -0.5,
        "T2_downweight": -1.5, "T3_strong_downweight": -3.0,
    }
    ceiling_ctrl = _s(get(cfg, f"{pfx}.established_gene_ceiling", "T1_watch")) or "T1_watch"
    ceiling_cds = _s(get(cfg, f"{pfx}.cds_fallback_ceiling", "T2_downweight")) or "T2_downweight"

    ratio = _num(excess_ratio)
    q = _num(q_nb)
    corrob = int(corroboration_count or 0)
    n_obs = _int(n_observed) or 0

    # A ratio-driven rule needs a minimum count behind it: a single variant against E = 0.1 is
    # 10x excess and means nothing. The NB FDR handles low counts correctly on its own (the
    # minimum n reaching q < 0.05 on the validation cohort was 5), so this floor guards only
    # the ratio rules — never the q rules.
    ratio_usable = ratio is not None and n_obs >= min_n_ratio

    tier, rule = "T0_no_downweight", "no_downweight"
    if ratio is not None and q is not None and ratio >= t1_ratio and q < t1_q:
        tier, rule = "T1_watch", f"ratio>={t1_ratio:g}&q<{t1_q:g}"
    if ratio_usable and ratio >= t2_ratio and corrob >= t2_corrob:
        tier, rule = "T2_downweight", f"ratio>={t2_ratio:g}&corrob>={t2_corrob}"
    t3_rules = []
    if q is not None and q < t3_q and corrob >= t3_q_corrob:
        t3_rules.append(f"q<{t3_q:g}&corrob>={t3_q_corrob}")
    if ratio_usable and ratio >= t3_ratio and corrob >= t3_ratio_corrob:
        t3_rules.append(f"ratio>={t3_ratio:g}&corrob>={t3_ratio_corrob}")
    if ratio_usable and ratio >= t3_sat_ratio and saturation and corrob >= t3_sat_corrob:
        t3_rules.append(f"ratio>={t3_sat_ratio:g}&saturation&corrob>={t3_sat_corrob}")
    if t3_rules:
        tier, rule = "T3_strong_downweight", "|".join(t3_rules)

    order = {t: i for i, t in enumerate(GENE_TIERS)}
    ctrl_applied = cds_applied = False
    if established_gene_control and order[tier] > order.get(ceiling_ctrl, 1):
        tier, ctrl_applied = ceiling_ctrl, True
    if _s(e_source) == "cds_fallback" and order[tier] > order.get(ceiling_cds, 2):
        tier, cds_applied = ceiling_cds, True

    # Review flags. `unexplained_excess` is the class this whole panel exists for: a high
    # excess with NO mechanism identified. Under the never-drop rule the correct action there
    # is to LOOK, not to down-weight — so it is flagged and carries no penalty.
    review = ""
    if ctrl_applied:
        review = "established_gene_high_excess"
    elif ratio is not None and ratio >= t3_ratio and corrob == 0:
        review = "unexplained_excess"
    return {
        "gene_tier": tier,
        "gene_artifact_penalty": float(penalties.get(tier, 0.0)),
        "control_ceiling_applied": ctrl_applied,
        "cds_ceiling_applied": cds_applied,
        "tier_rule": rule,
        "review_flag": review,
    }


def downweight_reason(gene_row, signals_reasons, tier_info) -> str:
    """The self-documenting reason string. It reports VALUES and never asserts biology.

    Example shapes from the validation cohort:
      ``excess_ratio=545.7(n=229,E=0.42,q=0); gnomad_constraint_flag=mis_too_many; oe_syn=1.63;
        classic_caf=0_bottom_decile; artifact_gene_family; cohort_saturation=1.04_...``
      ``excess_ratio=1.1(n=1,E=0.95,q=1); no_corroborating_signature``
    """
    ratio, E, q = _num(gene_row.get("excess_ratio")), _num(gene_row.get("E_expected")), _num(gene_row.get("q_nb"))
    n = _int(gene_row.get("n_observed")) or 0
    if ratio is None or E is None:
        head = f"excess_ratio=NA(n={n},E=NA,offset_source={_s(gene_row.get('E_source')) or 'none'})"
    else:
        head = (f"excess_ratio={ratio:.4g}(n={n},E={E:.3g},"
                f"q={'NA' if q is None else format(q, '.3g')})")
    parts = [head]
    parts += list(signals_reasons) or ["no_corroborating_signature"]
    if tier_info.get("control_ceiling_applied"):
        parts.append("established_gene_control_ceiling_applied")
    if tier_info.get("cds_ceiling_applied"):
        parts.append("cds_fallback_ceiling_applied")
    if _s(gene_row.get("E_source")) == "cds_fallback":
        parts.append("offset=cds_length_fallback(+/-30%)")
    return "; ".join(parts)


# =============================================================================
# 5. Per-variant tiering
# =============================================================================
def molecular_effect_class(consequence, impact, ref="", alt="") -> str:
    c = _s(consequence).lower()
    if any(t in c for t in CANONICAL_SPLICE):
        return "canonical_splice"
    if any(t in c for t in PLOF_CONSEQUENCES):
        return "plof"
    if "missense" in c:
        return "missense"
    if "inframe_insertion" in c or "inframe_deletion" in c:
        return "inframe_indel"
    if _s(impact).upper() == "HIGH":
        return "plof"
    return "synonymous_utr_intronic"


def spliceai_status(spliceai_ds, ref="", alt="") -> str:
    """``scored`` | ``not_covered``. A blank score is NOT a zero.

    hprv runs on the PRECOMPUTED SpliceAI set with backfill off by default, which covers every
    genome-wide SNV but only 1 nt insertions and deletions <= 4 nt. So a blank on a larger indel
    means "the precomputed set does not cover this variant" — and it bites exactly the indel
    classes most likely to disrupt splicing. Same blank-vs-zero trap as NHF.
    """
    return "scored" if _num(spliceai_ds) is not None else "not_covered"


def nmd_status(row=None) -> str:
    """Always ``INDETERMINATE`` today, and that is a resource gap, not a modelling choice.

    Abou Tayoun 2018 grades PVS1 by whether the predicted truncation triggers NMD: a nonsense
    or frameshift in the last exon, or the last 50 nt of the penultimate exon, escapes NMD and
    drops from Very Strong. ``variants.tsv`` carries no exon number, CDS position or transcript
    length, so the test cannot be evaluated and **no pLoF can reach V5**. The fix is small and
    specific — carry VEP's ``EXON``, ``CDS_position`` and the transcript exon count / CDS
    length through Step 2 into ``variants.tsv`` — and needs no new resource.
    """
    return "INDETERMINATE"


def assign_variant_tier(row, cfg=None) -> dict:
    """-> {variant_tier, variant_tier_reason, molecular_effect_class, nmd_status,
    plof_confidence, spliceai_status, missense_evidence_source}.

    Screening/triage tiers, NOT an ACMG classification: hprv assigns no ACMG weight, and the
    resource gaps below make several ACMG codes unimplementable.
    """
    cfg = cfg or {}
    pfx = "prioritization.variant_tier"
    sai_strong = _f(cfg, f"{pfx}.spliceai_strong", 0.5)
    sai_sup = _f(cfg, f"{pfx}.spliceai_supporting", 0.2)
    sai_benign = _f(cfg, f"{pfx}.spliceai_benign_max", 0.1)
    cadd_mis = _f(cfg, f"{pfx}.cadd_missense_supporting", 25.3)
    cadd_benign = _f(cfg, f"{pfx}.cadd_benign_max", 15.0)
    benign_impacts = {s.upper() for s in get(cfg, f"{pfx}.benign_impacts", ["LOW", "MODIFIER"])}

    cq = _s(row.get("consequence")).lower()
    imp = _s(row.get("impact")).upper()
    ds = _num(row.get("spliceai_ds"))
    cadd = _num(row.get("cadd"))
    mec = molecular_effect_class(row.get("consequence"), imp, row.get("ref"), row.get("alt"))
    sai_state = spliceai_status(row.get("spliceai_ds"))
    nmd = nmd_status(row)

    out = {
        "molecular_effect_class": mec,
        "nmd_status": nmd,
        # No LOFTEE data is bind-mounted, so there is no HC/LC pLoF confidence at any price
        # under the VEP-only contract — and no low-confidence flags either (ancestral allele,
        # single-exon transcript, 50-bp rule, non-canonical splice).
        "plof_confidence": "UNAVAILABLE",
        "spliceai_status": sai_state,
        "missense_evidence_source": "none",
    }

    # V0 — molecularly benign PREDICTION. Requires BOTH scores PRESENT and below their cutoffs:
    # a missing SpliceAI or CADD score is absence of evidence, and letting it satisfy a benign
    # rule would cap exactly the variants nobody scored. This is the single most consequential
    # tier because it hard-caps the composite total (see score_variant).
    if imp in benign_impacts and ds is not None and cadd is not None \
            and ds < sai_benign and cadd < cadd_benign:
        out.update(variant_tier="V0",
                   variant_tier_reason=f"benign_prediction:spliceai_ds={ds:.3g}<{sai_benign:g}"
                                       f"&cadd={cadd:.3g}<{cadd_benign:g}&impact={imp}")
        return out

    # V4 — strong splice, or an NMD-INDETERMINATE pLoF. V5 is defined so the ladder is complete
    # and the implementation has a target, but it is UNREACHABLE today: every pLoF lands here.
    if ds is not None and ds >= sai_strong:
        out.update(variant_tier="V4",
                   variant_tier_reason=f"spliceai_ds={ds:.3g}>={sai_strong:g}")
        return out
    if mec in ("plof", "canonical_splice") and imp == "HIGH":
        note = ("canonical_splice_site" if mec == "canonical_splice" else "plof")
        # Walker 2023 is explicit that a canonical site with a LOW SpliceAI score deserves
        # scrutiny, not automatic Very Strong — so report the score beside the tier.
        extra = "" if ds is None else f",spliceai_ds={ds:.3g}"
        out.update(variant_tier="V4",
                   variant_tier_reason=f"{note}:nmd_status=INDETERMINATE(V5_unreachable){extra}")
        return out

    # V3 — supporting splice, or a high-CADD missense. The CADD route is a DISCOVERY RANK, not
    # PP3 evidence: 25.3 is Pejaver 2022's missense-only PP3-supporting value, but ClinGen SVI
    # says commit to ONE predictor chosen before seeing results, and the calibrated choice is
    # REVEL (PP3 0.644/0.773/0.932). So the tier carries missense_evidence_source so no
    # downstream reader mistakes it for a calibrated call.
    if ds is not None and ds >= sai_sup:
        out.update(variant_tier="V3", variant_tier_reason=f"spliceai_ds={ds:.3g}>={sai_sup:g}")
        return out
    if mec == "missense":
        # PREDICTOR PRECEDENCE, and it is deliberate. ClinGen SVI's rule is to commit to ONE
        # predictor chosen BEFORE seeing results — so this consults them in a fixed order and
        # reports which one spoke, rather than taking the max over whatever is available (that
        # would be a best-of-N cherry-pick with no calibration behind it). REVEL first because
        # its PP3/BP4 cut points are the ones Pejaver 2022 calibrated; AlphaMissense next (SVI
        # endorses it on par, and it reaches Strong where REVEL sits at Supporting); CADD last
        # and explicitly labelled off-label. A missing score falls through to the next source —
        # both are missense-only and neither covers every substitution.
        rev = _num(row.get("revel"))
        am = _num(row.get("alphamissense"))
        rev_sup = _f(cfg, f"{pfx}.revel_supporting", 0.644)
        rev_mod = _f(cfg, f"{pfx}.revel_moderate", 0.773)
        rev_benign = _f(cfg, f"{pfx}.revel_benign_max", 0.290)
        am_sup = _f(cfg, f"{pfx}.alphamissense_supporting", 0.564)
        am_benign = _f(cfg, f"{pfx}.alphamissense_benign_max", 0.34)
        if rev is not None:
            if rev >= rev_sup:
                # V4 at moderate-or-better: a calibrated predictor above Pejaver's moderate cut
                # is stronger evidence than the supporting-only rungs that share V3.
                tier = "V4" if rev >= rev_mod else "V3"
                out.update(variant_tier=tier, missense_evidence_source="revel",
                           variant_tier_reason=f"missense&revel={rev:.3g}>={rev_sup:g}"
                                               f"(ClinGen-calibrated,Pejaver2022)")
                return out
            if rev <= rev_benign:
                out.update(variant_tier="V1", missense_evidence_source="revel",
                           variant_tier_reason=f"missense&revel={rev:.3g}<={rev_benign:g}"
                                               "(calibrated_BP4-supporting_range)")
                return out
            out.update(variant_tier="V2", missense_evidence_source="revel",
                       variant_tier_reason=f"missense&revel={rev:.3g}_between_cuts"
                                           "(calibrated_but_indeterminate)")
            return out
        if am is not None:
            if am >= am_sup:
                out.update(variant_tier="V3", missense_evidence_source="alphamissense",
                           variant_tier_reason=f"missense&am_pathogenicity={am:.3g}>={am_sup:g}"
                                               "(REVEL_absent;SVI-endorsed)")
                return out
            if am <= am_benign:
                out.update(variant_tier="V1", missense_evidence_source="alphamissense",
                           variant_tier_reason=f"missense&am_pathogenicity={am:.3g}<={am_benign:g}")
                return out
            out.update(variant_tier="V2", missense_evidence_source="alphamissense",
                       variant_tier_reason=f"missense&am_pathogenicity={am:.3g}_between_cuts")
            return out
        if cadd is not None and cadd >= cadd_mis:
            out.update(variant_tier="V3", missense_evidence_source="cadd_offlabel",
                       variant_tier_reason=f"missense&cadd={cadd:.3g}>={cadd_mis:g}"
                                           "(off-label:calibrated_on_missense_but_not_the_"
                                           "ClinGen-recommended_predictor)")
            return out
        out.update(variant_tier="V2", missense_evidence_source="none",
                   variant_tier_reason="missense_below_cadd_cut(no_calibrated_predictor:"
                                       "REVEL/AlphaMissense_absent)")
        return out

    # V2 — in-frame indel. No calibrated in-frame predictor exists and the mechanism (in-frame
    # loss of a domain vs neutral) is unresolvable from annotation alone; report the length,
    # since a large in-frame deletion is a different proposition from a 3 nt one.
    if mec == "inframe_indel":
        ln = abs(len(_s(row.get("alt"))) - len(_s(row.get("ref"))))
        out.update(variant_tier="V2",
                   variant_tier_reason=f"inframe_indel(length={ln}nt,no_calibrated_predictor)")
        return out

    # V1 — non-coding / synonymous discovery rank. These reach the candidate list only via the
    # CADD rung (~top 0.3% genome-wide); there is no ClinGen-endorsed non-coding CADD threshold
    # to replace 25.3 with.
    if imp in benign_impacts or mec == "synonymous_utr_intronic":
        why = f"cadd={cadd:.3g}" if cadd is not None else "no_score"
        out.update(variant_tier="V1",
                   variant_tier_reason=f"non_coding_or_synonymous_discovery_rank({why})")
        return out

    out.update(variant_tier="V2",
               variant_tier_reason=f"uncalibrated_moderate_effect(consequence={cq or 'NA'})")
    return out


def rarity_strength(af, cfg=None) -> str:
    """``strong`` | ``moderate`` | ``supporting`` | ``permissive`` | ``fail`` | ``unknown``.

    Two caveats that must travel with this column. (1) The oracle is a **point estimate, not
    faf95** — the VEP cache carries no AC/AN, so the CI correction is unrecoverable at any
    price; do not present these bands as ACMG PM2. (2) Audit A-4: the proxy is the MAX over
    grpmax-eligible groups regardless of the cohort's ancestry composition, so effective
    stringency varies with the proband's ancestry and the loss is unevenly distributed across
    ancestry groups, in the silent false-negative direction. An absent AF reads as ``unknown``
    (treated as rarest downstream), never as a measured zero.
    """
    cfg = cfg or {}
    pfx = "prioritization.variant_tier.rarity"
    v = _num(af)
    if v is None:
        return "unknown"
    if v >= _f(cfg, f"{pfx}.ba1", 0.05):
        return "fail"
    if v < _f(cfg, f"{pfx}.strong_max", 1.0e-5):
        return "strong"
    if v < _f(cfg, f"{pfx}.moderate_max", 1.0e-4):
        return "moderate"
    if v < _f(cfg, f"{pfx}.supporting_max", 1.0e-3):
        return "supporting"
    if v < _f(cfg, f"{pfx}.permissive_max", 1.0e-2):
        return "permissive"
    return "permissive_fail"


def rarity_driven_by_single_group(grpmax_af, max_af, cfg=None) -> bool:
    """True when the rarity term was ONE population's doing (audit A-4 visibility).

    ``variants.tsv`` carries no per-group breakdown, only ``grpmax_af`` (max over eligible
    groups) and ``max_af`` (max over ALL groups, including the bottlenecked founder groups
    grpmax excludes). A ``max_af`` far above ``grpmax_af`` is the observable shadow of a single
    group driving the number, so that ratio is the available proxy for the spec's
    "second-highest group" comparison — a deviation forced by the column set, flagged rather
    than dropped so a reviewer can still see when one population set the rarity.
    """
    cfg = cfg or {}
    ratio = _f(cfg, "prioritization.variant_tier.rarity.single_group_ratio", 10.0)
    g, m = _num(grpmax_af), _num(max_af)
    if m is None:
        return False
    if g is None or g <= 0.0:
        return m > 0.0        # eligible groups say nothing; some other group carries it
    return m / g > ratio


def genotype_qc(row, cfg=None):
    """-> (pass: bool, fail_reason: str) for the PROBAND only. A failing call is flagged and
    penalised, NEVER dropped.

    Audit A-5, which belongs in the methods: these are genotype-refined VCFs, so ``GQ`` is
    derived from POSTERIOR probabilities (``PP``), not raw ``PL``. A GQ >= 20 cut on a
    posterior is a different, generally more permissive filter than on a likelihood, because
    the posterior has already been sharpened by a FAMILY prior — the same trio structure the
    inheritance model is about to use. The classification is therefore not statistically
    independent of the prior that produced the genotypes, and it biases toward
    cleaner-looking Mendelian patterns than the raw data support.

    Missing metrics are ``unknown``, not failures: a blank GQ means Step 5 did not record it.
    """
    cfg = cfg or {}
    min_gq = _f(cfg, "filters.genotype_qc.min_gq", 20.0)
    min_dp = _f(cfg, "filters.genotype_qc.min_dp", 10.0)
    ab_lo = _f(cfg, "filters.genotype_qc.het_ab_min", 0.25)
    ab_hi = _f(cfg, "filters.genotype_qc.het_ab_max", 0.75)
    ab_hom = _f(cfg, "filters.genotype_qc.homalt_ab_min", 0.90)

    gq, dp, ab = _num(row.get("child_GQ")), _num(row.get("child_DP")), _num(row.get("child_AB"))
    gt = _s(row.get("child_gt"))
    reasons = []
    if gq is not None and gq < min_gq:
        reasons.append(f"GQ={gq:.3g}<{min_gq:g}")
    if dp is not None and dp < min_dp:
        reasons.append(f"DP={dp:.3g}<{min_dp:g}")
    if ab is not None:
        hom_alt = gt in ("1/1", "1|1")
        if hom_alt and ab < ab_hom:
            reasons.append(f"homalt_AB={ab:.3g}<{ab_hom:g}")
        elif not hom_alt and not (ab_lo <= ab <= ab_hi):
            reasons.append(f"het_AB={ab:.3g}_outside[{ab_lo:g},{ab_hi:g}]")
    if reasons:
        return False, ";".join(reasons)
    if gq is None and dp is None and ab is None:
        return True, "qc_metrics_absent(not_assessed)"
    return True, ""


NHF_MEMBERS = ("child", "mother", "father")


def nhf_state(row, threshold: float = 0.5, min_reads: int = 5):
    """-> (state, max_fraction, max_reads) where state is ``clean`` | ``flagged`` |
    ``not_screened``.

    **Blank is NOT zero**, and this is the load-bearing rule (copied from
    ``prepare_igv_variants.py:nhf_status``):
      * blank  = NOT SCREENED — the member is not an ALT carrier, has no mini-CRAM, or Step 8b
        never ran (no kraken2 DB).
      * ``0.0`` = SCREENED, and every read classified human.

    Treating blank as clean silently promotes exactly the calls nobody examined, so this
    returns THREE states, never two. The ``min_reads`` floor of 5 is essential: an NHF of 1.0
    over 2 reads is noise, not evidence — always carry the ``*_nhf_reads`` denominator beside
    the fraction.
    """
    screened = False
    flag = False
    best_frac = best_reads = None
    for m in NHF_MEMBERS:
        frac = _num(row.get(f"{m}_nhf"))
        if frac is None:
            continue
        screened = True
        reads = _num(row.get(f"{m}_nhf_reads")) or 0.0
        if best_frac is None or frac > best_frac:
            best_frac, best_reads = frac, reads
        if frac >= threshold and reads >= min_reads:
            flag = True
    if not screened:
        return "not_screened", None, None
    return ("flagged" if flag else "clean"), best_frac, best_reads


def clinvar_strength(clin_sig) -> str:
    """``p_lp`` | ``conflicting`` | ``vus`` | ``benign`` | ``absent``.

    Classifies the SIGNIFICANCE string only. Review status is a separate axis and arrives on its
    own column (``clinvar_stars``, from the ClinVar VCF transfer in Step 2) — the caller applies
    the star gate. This value used to be named ``p_lp_no_star_gate`` because no star count
    existed at any price under the cache-only contract; it does now, so the name no longer
    encodes that limitation. When the transfer has not run, ``clinvar_stars`` is blank and the
    clinical term stays at full weight rather than being damped as if unreviewed.

    Matches both the VEP form (lowercase, ``&``-joined) and the ClinVar-VCF form (Capitalised,
    ``/``- or ``,``-joined), like ``annotations.clnsig_is_plp``.
    """
    s = _s(clin_sig).lower()
    if not s:
        return "absent"
    if "conflicting" in s:
        return "conflicting"
    if "pathogenic" in s and "likely_benign" not in s and "benign/likely" not in s:
        return "p_lp"
    if "benign" in s:
        return "benign"
    if "uncertain" in s:
        return "vus"
    return "vus"


def moi_coherence(inheritance, gene_moi, long_gene=False):
    """-> (coherence, caveat). ``unknown`` must be EXACTLY neutral — the novel-gene case.

    Any penalty on ``moi_unknown`` converts the score into a known-gene filter and destroys
    novel-gene discovery, violating the never-drop rule in spirit if not in letter.

    Audit A-6: mode assignment is single-gene-keyed, and a phase-confirmed mat x pat pair
    consumes both legs, so a genuinely dominant-grade variant is relabelled ``compound_het``
    whenever the child carries any other sub-1e-2 functional het in the same gene — near-certain
    in long genes (TTN, NEB, RYR1, DMD). So no discordance penalty is applied to a comp-het
    call in a long gene; the caveat is reported instead.

    **The het-in-recessive-gene case is a FLAG, never a penalty, and this is load-bearing.** A
    heterozygous observation in a canonically autosomal-recessive gene is the single most likely
    shape for a *carrier-risk* hypothesis, which is a different genetic model from the recessive
    syndrome the gene is curated for — not an incoherent call. The concrete case: the FA/HR genes
    (``FANCA``, ``FANCD2``, ``SLX4``, ``FANCE``, ``BRCA2``) carry germ-cell-tumour evidence about
    HETEROZYGOUS carriers (PMID 40906985: five-gene combined OR 10.17, 95% CI 4.87-21.27,
    P = 2.90e-06) while their canonical Mendelian model — and their PanelApp green status — is
    biallelic Fanconi anemia. Charging those rows a discordance penalty would penalise exactly
    the observation the evidence is about. So the mismatch emits
    ``moi_mismatch_het_in_recessive_gene`` and ``score_variant`` suppresses the penalty; the
    reviewer still sees the mismatch in the column. Incomplete penetrance, mosaicism, an
    undetected second hit and a genuinely novel mechanism produce the same shape, which is the
    general form of the same argument.
    """
    mode = _s(inheritance).lower()
    curated = {t.strip().upper() for t in re.split(r"[,;|/]", _s(gene_moi)) if t.strip()}
    caveat = ""
    if mode == "compound_het" and long_gene:
        caveat = "long_gene_comphet_drift"
    if not curated:
        return "unknown", caveat
    dominant_like = {"AD", "AUTOSOMAL DOMINANT", "DOMINANT", "XLD", "SD"}
    recessive_like = {"AR", "AUTOSOMAL RECESSIVE", "RECESSIVE", "XLR", "XL"}
    # The overlay's own MOI vocabulary is richer than a bare token: the FA rows read
    # `AR_biallelic_FA;heterozygous_carrier_risk_proposed`. Recognise both halves — the prefix
    # match makes `AR_BIALLELIC_FA` recessive-like, and an explicit carrier-risk annotation is
    # itself the statement that a het observation is IN model.
    carrier_risk = any("CARRIER" in t or "HETEROZYGOUS" in t for t in curated)
    has_dom = bool(curated & dominant_like) or any(t.startswith(("AD", "XLD")) for t in curated)
    has_rec = bool(curated & recessive_like) or any(t.startswith(("AR", "XLR")) for t in curated)
    if mode in ("dominant", "denovo", "denovo_x_hemi"):
        observed = "dominant"
    elif mode in RECESSIVE_MODES:
        observed = "recessive"
    else:
        return "unknown", caveat
    if observed == "dominant":
        if has_dom or carrier_risk:
            # An explicit carrier-risk annotation means a het observation is in model, so this is
            # coherent rather than a suppressed mismatch.
            return "coherent", caveat
        if has_rec:
            return "discordant", (caveat or "moi_mismatch_het_in_recessive_gene")
        return "unknown", caveat
    return ("coherent" if has_rec else ("discordant" if has_dom else "unknown")), caveat


# =============================================================================
# 5b. The optional Class-B gene-prior overlay
# =============================================================================
# Evidence classes that are NOT germline predisposition evidence and must therefore contribute
# ZERO germline prior, whatever weight the overlay assigns them. A curated overlay lists somatic
# drivers deliberately — in the GCT resource that is exactly SIX rows, all tier T4 at weight 0.15:
# AKT1, BCORL1, CBL, KRAS, MTOR, NRAS — precisely so a reader can see they were considered and
# excluded. A 0.15 weight on a somatic driver is a bookkeeping placeholder, not weak germline
# support, and scoring it as the latter would be the exact misreading the resource's own `do_not`
# list warns against.
#
# MATCH ON `evidence_class`, NEVER ON A REMEMBERED GENE LIST. KIT is the trap: it is frequently
# named among germ-cell-tumour somatic drivers, but its row in this resource is
# `tier=T3, evidence_class=gwas_common_variant_locus, prior_weight=0.35` — a GWAS-locus row that
# this exclusion does NOT touch and should not, since it contributes a legitimate (if weak) T3
# prior. A gene can be a somatic driver in the literature and a GWAS-locus row in the file; only
# the file decides what the code does.
NON_GERMLINE_EVIDENCE_CLASSES = ("somatic_driver_not_germline",)


def parse_gene_prior_overlay(rows, gene_sets=None, cfg=None) -> dict:
    """Parse an overlay into ``{GENE: {points_scale, tier, evidence_class, ...}}``.

    **The join is on gene symbol ONLY — deliberately MOI-agnostic.** Routing a prior lookup by a
    gene's canonical mode of inheritance would silently miss every hypothesis stated about a
    different genetic model than the gene is curated for. The concrete case this exists to
    prevent: the FA/HR genes carry germ-cell-tumour evidence about HETEROZYGOUS carriers while
    their canonical model is biallelic Fanconi anemia, so an MOI-routed lookup would consult
    those rows under a recessive model and never apply them to the het observations the evidence
    is actually about. A gene's prior must be expressible independently of its canonical MOI.

    ``rows`` is a list of dicts (a parsed TSV) or a list of bare symbols. Recognised optional
    columns — everything else passes through untouched, so an overlay may carry its own
    provenance fields without this parser needing to know them:

      * ``prior_weight`` — a 0..1 scale on the configured maximum. **UNCALIBRATED**: an ordering
        default, NOT a likelihood ratio, and it must never be presented as one. Absent -> 1.0.
      * ``tier`` / ``evidence_class`` / ``replication`` / ``moi`` / ``pmids`` — reported, never
        used to route the join.
      * ``gene_sets`` — set memberships; see the double-counting guard below.

    ``gene_sets`` maps a set id to ``{"prior_weight": w, "members": [...]}`` (the shape of the
    GCT resource's JSON ``gene_sets`` block). A set-level prior and a gene-level prior for the
    same variant are combined by **MAX, never SUM** — when both derive from the same study,
    summing double-counts one study. The GCT resource is exactly this case: ``FA_HR_PATHWAY_23``
    (pooled OR 4.14, 95% CI 1.98-8.66, P = .0013) and its per-gene FA rows (combined OR 10.17)
    both come from PMID 40906985, and its usage contract says so in as many words. Enforced here
    rather than documented, because a documented invariant is one nobody checks.
    """
    weight_cap = _f(cfg, "prioritization.composite.gene_list_prior.max_prior_weight", 1.0)
    honor = bool(get(cfg, "prioritization.composite.gene_list_prior.honor_prior_weight", True))
    drop_classes = tuple(get(cfg, "prioritization.composite.gene_list_prior.non_germline_classes",
                             None) or NON_GERMLINE_EVIDENCE_CLASSES)
    out = {}
    for r in rows or ():
        if not isinstance(r, dict):
            sym = _s(r)
            if sym and sym.lower() not in GENE_KEYS_LOWER:
                out[sym.upper()] = {"gene": sym, "points_scale": 1.0, "tier": "",
                                    "evidence_class": "", "prior_weight_raw": None,
                                    "germline_excluded": False, "source": "symbol_list"}
            continue
        sym = _s(r.get("gene") or r.get("symbol") or r.get("gene_symbol"))
        if not sym or sym.lower() in GENE_KEYS_LOWER:
            continue
        raw = _num(r.get("prior_weight"))
        scale = 1.0 if (raw is None or not honor) else max(0.0, min(float(raw) / weight_cap, 1.0))
        klass = _s(r.get("evidence_class"))
        excluded = any(c and c in klass for c in drop_classes)
        if excluded:
            # Reported with its own flag and a ZERO scale — visible in the output, contributing
            # nothing. Silently dropping the row instead would hide that it was considered.
            scale = 0.0
        entry = dict(r)
        entry.update(gene=sym, points_scale=scale, prior_weight_raw=raw,
                     evidence_class=klass, tier=_s(r.get("tier")),
                     germline_excluded=excluded, source="gene_row",
                     set_ids=[s for s in re.split(r"[,;|]", _s(r.get("gene_sets"))) if s.strip()])
        out[sym.upper()] = entry
    # Set-level priors: MAX against any gene-level prior, never a sum.
    for set_id, spec in (gene_sets or {}).items():
        if not isinstance(spec, dict):
            continue
        sraw = _num(spec.get("prior_weight"))
        sscale = 1.0 if (sraw is None or not honor) else max(0.0, min(float(sraw) / weight_cap, 1.0))
        for m in spec.get("members") or ():
            key = _s(m).upper()
            if not key:
                continue
            cur = out.get(key)
            if cur is None:
                out[key] = {"gene": _s(m), "points_scale": sscale, "tier": "",
                            "evidence_class": f"set:{set_id}", "prior_weight_raw": sraw,
                            "germline_excluded": False, "source": "gene_set",
                            "set_ids": [set_id], "set_id_applied": set_id}
            elif cur.get("germline_excluded"):
                continue          # a non-germline row is not rescued by set membership
            elif sscale > cur["points_scale"]:
                # MAX, not sum — see the docstring. Record that the set won, so the reason string
                # can say which prior was applied.
                cur.update(points_scale=sscale, set_id_applied=set_id)
                cur.setdefault("set_ids", []).append(set_id)
            else:
                cur.setdefault("set_ids", []).append(set_id)
    return out


def gene_prior_points(entry, max_points: float, gate: float) -> float:
    """Points contributed by an overlay hit: ``max_points * points_scale * mechanism_gate``.

    The mechanism gate is the same one constraint gets, for the same reason: a gene-list prior
    must never rescue a molecularly-benign prediction, or the prior becomes confirmation bias.
    ``entry`` may be ``None``/``False`` (no overlay hit) or ``True`` (a bare membership flag, for
    the unweighted case and for backward compatibility).
    """
    if not entry:
        return 0.0
    scale = 1.0 if entry is True else float(entry.get("points_scale", 1.0) or 0.0)
    return float(max_points) * scale * float(gate)


# =============================================================================
# 6. The additive composite
# =============================================================================
def default_weights() -> dict:
    """The Tavtigian-style point weights. Sign convention: positive raises priority.

    Borrowed from Tavtigian 2020's naturally-scaled point system for THREE reasons, each a
    design constraint: auditability (a reviewer reads ``spliceai=+3, rarity=+2, artifact=-3``,
    not ``0.71``); additivity survives a missing term honestly (under the VEP-only contract
    several terms are simply unavailable, and an additive scheme degrades to "that term
    contributed 0" rather than silently redistributing the missing evidence); and no labelled
    truth set exists here, so a fitted model would be fit on the very artifacts we are removing.

    **Where the analogy stops**: these are NOT ACMG points and the total must not be read
    against Tavtigian's P >= 10 / LP 6-9 / VUS 0-5 bands. The criteria are not ACMG criteria (a
    CADD-based term is not PP3), no phenotype/segregation/functional evidence exists, the
    ClinVar term has no review-status gate, and the artifact-penalty terms have no ACMG
    analogue at all. The column is ``priority_points``, never ``acmg_points``, and no P/LP/VUS
    label is ever emitted from it.
    """
    return {
        "molecular": {"V5": 8.0, "V4": 4.0, "V3": 2.0, "V2": 1.0, "V1": 0.5, "V0": 0.0},
        "rarity": {"strong": 2.0, "moderate": 1.5, "supporting": 1.0, "permissive": 0.5,
                   "permissive_fail": 0.0, "unknown": 2.0, "fail": -8.0},
        "gene_constraint": 1.0,
        "recurrence": {"distinct_2": 1.0, "distinct_3plus": 2.0, "same_variant": 0.5},
        "quality": {"gt_fail": -2.0, "nhf_flagged": -3.0, "nhf_not_screened": 0.0,
                    "partner_unknown": -0.5},
        "clinical": {"p_lp": 4.0, "conflicting_vus": 0.0, "benign": -4.0},
        "moi": {"discordant": -1.0, "unknown": 0.0, "coherent": 0.0},
        "gene_artifact": {"T0_no_downweight": 0.0, "T1_watch": -0.5,
                          "T2_downweight": -1.5, "T3_strong_downweight": -3.0},
        "gene_list_prior": 2.0,
    }


def _merge_weights(cfg) -> dict:
    w = default_weights()
    over = get(cfg or {}, "prioritization.composite.weights", None) or {}
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(w.get(k), dict):
            w[k].update(v)
        else:
            w[k] = v
    # The rarity band above BA1 is named `ba1` in the config (that is the criterion it
    # implements) and `fail` in rarity_strength()'s return vocabulary. Alias rather than rename
    # either: a config that sets `ba1: -8` must actually change the term, not silently leave the
    # default in place under a key nothing reads.
    if "ba1" in w["rarity"]:
        w["rarity"]["fail"] = float(w["rarity"]["ba1"])
    return w


def constraint_gate(variant_tier, inheritance, cfg=None) -> float:
    """Mechanism-gating multiplier on the gene-constraint term AND on the optional gene-list
    prior — the ACMG/ClinGen SVI principle and the single most important structural rule here.

      * ``V0`` -> 0.0: a constrained gene, or a gene-list membership, CANNOT rescue a
        molecularly-benign prediction. That is how gene-list priors turn into confirmation bias.
      * ``V1``/``V2`` -> 0.5: half credit for a discovery-rank or uncalibrated molecular effect.
      * ``V3``-``V5`` -> 1.0: a credible molecular effect earns the full term.
      * Recessive modes -> 0.0 regardless of tier: pLoF constraint measures selection against
        heterozygotes and is not evidence about a biallelic candidate.
    """
    cfg = cfg or {}
    gating = get(cfg, "prioritization.composite.gene_constraint_gating", None) or {
        "V0": 0.0, "V1": 0.5, "V2": 0.5, "V3": 1.0, "V4": 1.0, "V5": 1.0}
    zero_rec = get(cfg, "prioritization.composite.zero_constraint_for_recessive_modes", True)
    if zero_rec and _s(inheritance).lower() in RECESSIVE_MODES:
        return 0.0
    return float(gating.get(_s(variant_tier), 1.0))


def gene_is_constrained(gene_row, cfg=None) -> bool:
    """pLI >= 0.9 or LOEUF < 0.35 (hprv canonical defaults). ``s_het >= 0.1`` is specified as a
    short-gene route but is TARGET: it is not in the gnomAD constraint columns this step reads.
    """
    cfg = cfg or {}
    pli_min = _f(cfg, "filters.constraint_weighting.pli_min", 0.90)
    loeuf_max = _f(cfg, "filters.constraint_weighting.loeuf_v2_tier1", 0.35)
    shet_min = _f(cfg, "filters.constraint_weighting.shet_min", 0.10)
    pli = _num(gene_row.get("pLI") if gene_row.get("pLI") is not None else gene_row.get("pli"))
    loeuf = _num(gene_row.get("oe_lof_upper"))
    shet = _num(gene_row.get("s_het"))
    return bool((pli is not None and pli >= pli_min)
                or (loeuf is not None and loeuf < loeuf_max)
                or (shet is not None and shet >= shet_min))


def score_variant(row, gene_row=None, cfg=None, gene_prior=False) -> dict:
    """The additive composite. EVERY term is returned as its own column.

    Returns ``pts_*`` per term, both totals (``priority_points_agnostic`` /
    ``priority_points_prior``), ``cap_applied``, and the evidence columns behind each term.
    Never removes or nullifies a row.
    """
    cfg = cfg or {}
    gene_row = gene_row or {}
    w = _merge_weights(cfg)

    tier_info = assign_variant_tier(row, cfg)
    vt = tier_info["variant_tier"]
    out = dict(tier_info)

    # --- 1.1 molecular evidence ---
    out["pts_molecular"] = float(w["molecular"].get(vt, 0.0))

    # --- 1.2 rarity ---
    af_col = row.get("grpmax_af") if row.get("grpmax_af") not in (None, "") else row.get("frequency")
    rs = rarity_strength(af_col, cfg)
    out["rarity_strength"] = rs
    out["pts_rarity"] = float(w["rarity"].get(rs, 0.0))
    out["rarity_driven_by_single_group"] = rarity_driven_by_single_group(
        af_col, row.get("max_af"), cfg)

    # --- 1.3 gene constraint, MECHANISM-GATED (never unconditional) ---
    gate = constraint_gate(vt, row.get("inheritance"), cfg)
    out["constraint_gate"] = gate
    out["pts_gene_constraint"] = (float(w["gene_constraint"]) * gate
                                  if gene_is_constrained(gene_row, cfg) else 0.0)

    # --- 1.4 recurrence: a RANK contribution, capped, never presented as significance ---
    # Step 6's recurrence null is a case-only approximation built only from variants observed in
    # the cohort, so p is too small and saturates: for essentially any gene with >= min_carriers
    # carriers of private variants it clears the exome-wide line (audit A-2). Score on the
    # CARRIER COUNT, never on p_recurrence — a p-value that saturates cannot order anything.
    n_car = _int(gene_row.get("n_carriers")) or 0
    kind = _s(gene_row.get("recurrence_kind"))
    if kind == "same_variant" and n_car >= 2:
        # One recurrent site shared by many trios is as easily a mapping/caller artifact or a
        # founder allele as a burden signal, so it gets strictly less credit.
        out["pts_recurrence"] = float(w["recurrence"]["same_variant"])
    elif n_car >= 3:
        out["pts_recurrence"] = float(w["recurrence"]["distinct_3plus"])
    elif n_car >= 2:
        out["pts_recurrence"] = float(w["recurrence"]["distinct_2"])
    else:
        out["pts_recurrence"] = 0.0
    out["same_variant_recurrence"] = kind == "same_variant"

    # --- 1.5 quality penalties ---
    gt_ok, gt_reason = genotype_qc(row, cfg)
    nhf_thr = _f(cfg, "prioritization.composite.nhf.flag_fraction", 0.5)
    nhf_min_reads = _i(cfg, "outputs.igv.nonhuman_screen.min_reads", 5)
    nhf, nhf_frac, nhf_reads = nhf_state(row, nhf_thr, nhf_min_reads)
    partner_unknown = _s(row.get("inheritance")).lower() == "compound_het"
    q = 0.0
    if not gt_ok:
        q += float(w["quality"]["gt_fail"])
    if nhf == "flagged":
        q += float(w["quality"]["nhf_flagged"])
    elif nhf == "not_screened":
        # An explicit uncertainty flag: not a penalty and NOT credit. Scoring an unscreened
        # call as clean would silently promote exactly the calls nobody examined.
        q += float(w["quality"]["nhf_not_screened"])
    if partner_unknown:
        # Parental GQ/DP/AB are absent from variants.tsv, so the leg that establishes trans
        # phase for a comp-het cannot be quality-assessed. Flag it; do not assume it passed.
        q += float(w["quality"]["partner_unknown"])
    out.update(pts_quality=q, gt_qc_pass=gt_ok, gt_qc_fail_reason=gt_reason,
               nhf_status=nhf, nhf_max_fraction=nhf_frac, nhf_max_reads=nhf_reads,
               partner_leg_quality_unknown=partner_unknown)

    # --- 1.6 clinical ---
    cs = clinvar_strength(row.get("clin_sig"))
    out["clinvar_strength"] = cs
    # GOLD STARS. Blank/absent = the ClinVar VCF transfer did not run, which is NOT the same as
    # 0 stars ("submitter provided no assertion criteria"). Conflating them would let a run with
    # no ClinVar resource damp every P/LP assertion as though it were unreviewed, so an absent
    # star count leaves the term at FULL weight and merely reports UNAVAILABLE — the same
    # never-drop logic the NHF three-state rule uses.
    stars = _num(row.get("clinvar_stars"))
    min_stars = _f(cfg, "resources.clinvar.min_review_stars", 2.0)
    low_scale = _f(cfg, "resources.clinvar.low_star_scale", 0.5)
    base = {"p_lp": float(w["clinical"]["p_lp"]),
            "benign": float(w["clinical"]["benign"])}.get(
                cs, float(w["clinical"]["conflicting_vus"]))
    if stars is None:
        out["clinvar_review_status"] = "UNAVAILABLE"
        out["clinvar_stars"] = ""
        out["pts_clinical"] = base
    else:
        out["clinvar_review_status"] = f"{int(stars)}_star"
        out["clinvar_stars"] = int(stars)
        # Damp, never zero by default: a 1-star P/LP assertion is weaker evidence than a 3-star
        # one but it is still evidence, and this layer re-ranks rather than filters. Only the
        # POSITIVE limb is gated — a low-star BENIGN call should not have its (negative) weight
        # shrunk toward zero, because that would PROMOTE a poorly-reviewed benign assertion.
        out["pts_clinical"] = base * low_scale if (stars < min_stars and base > 0) else base

    # --- 1.7 inheritance-model coherence ---
    long_bp = _f(cfg, "prioritization.variant_tier.moi.long_gene_cds_bp", 10000.0)
    cds = _num(gene_row.get("cds_length"))
    coh, caveat = moi_coherence(row.get("inheritance"), gene_row.get("gene_moi"),
                                long_gene=(cds is not None and cds >= long_bp))
    out["moi_coherence"], out["moi_caveat"] = coh, caveat or "none"
    # Both caveats SUPPRESS the discordance penalty entirely and report the reason instead. They
    # are the two shapes where "discordant" is an artefact of how the mode was assigned or of
    # which genetic model the gene happens to be curated under, rather than evidence against the
    # call: `long_gene_comphet_drift` (audit A-6 — mode assignment is single-gene-keyed) and
    # `moi_mismatch_het_in_recessive_gene` (a het observation in a canonically-recessive gene is
    # the carrier-risk shape, e.g. the FA/HR germ-cell-tumour evidence, PMID 40906985). A flag,
    # never a down-weight — the same never-drop logic applied to the coherence layer.
    SUPPRESSED_MOI_CAVEATS = ("long_gene_comphet_drift", "moi_mismatch_het_in_recessive_gene")
    out["pts_moi"] = 0.0 if (coh == "discordant" and caveat in SUPPRESSED_MOI_CAVEATS) \
        else float(w["moi"].get(coh, 0.0))

    # --- 1.8 gene artifact penalty ---
    gt_tier = _s(gene_row.get("gene_tier")) or "T0_no_downweight"
    pen = _num(gene_row.get("gene_artifact_penalty"))
    out["pts_gene_artifact"] = pen if pen is not None else float(w["gene_artifact"].get(gt_tier, 0.0))

    # --- 1.9 the Class-B overlay: OFF by default, ONE term, mechanism-gated like constraint ---
    # `gene_prior` is False/None (no hit), True (a bare symbol-list membership), or an overlay
    # ENTRY dict from parse_gene_prior_overlay — which carries a `points_scale` reflecting the
    # overlay's own per-gene weight and the non-germline exclusion. The join that produced it is
    # MOI-agnostic by construction (gene symbol only), so a het-carrier hypothesis in a
    # canonically-recessive gene is applied to exactly the het observations it is about.
    out["pts_gene_list_prior"] = gene_prior_points(gene_prior, w["gene_list_prior"], gate)
    out["gene_list_prior_member"] = bool(gene_prior)
    if isinstance(gene_prior, dict):
        out["gene_list_prior_tier"] = _s(gene_prior.get("tier"))
        out["gene_list_prior_weight"] = gene_prior.get("prior_weight_raw")
        out["gene_list_prior_evidence_class"] = _s(gene_prior.get("evidence_class"))
        out["gene_list_prior_excluded_non_germline"] = bool(gene_prior.get("germline_excluded"))
        out["gene_list_prior_set_applied"] = _s(gene_prior.get("set_id_applied"))
    else:
        out["gene_list_prior_tier"] = ""
        out["gene_list_prior_weight"] = None
        out["gene_list_prior_evidence_class"] = ""
        out["gene_list_prior_excluded_non_germline"] = False
        out["gene_list_prior_set_applied"] = ""

    terms = ("pts_molecular", "pts_rarity", "pts_gene_constraint", "pts_recurrence",
             "pts_quality", "pts_clinical", "pts_moi", "pts_gene_artifact")
    total = sum(out[t] for t in terms)
    prior_total = total + out["pts_gene_list_prior"]

    # --- 2. hard caps. Both implement mechanism gating: molecular benignity and BA1-level
    # frequency are statements about the VARIANT that no amount of gene-level enthusiasm can
    # overturn. Neither cap removes the variant from the output. ---
    cap = "none"
    v0_cap = _f(cfg, "prioritization.composite.caps.v0_benign_max", 0.0)
    ba1_cap = _f(cfg, "prioritization.composite.caps.ba1_max", -4.0)
    if vt == "V0":
        total, prior_total, cap = min(total, v0_cap), min(prior_total, v0_cap), "V0_benign"
    if rs == "fail":
        total, prior_total = min(total, ba1_cap), min(prior_total, ba1_cap)
        cap = "BA1_frequency" if cap == "none" else cap + "+BA1_frequency"
    out.update(priority_points_agnostic=total, priority_points_prior=prior_total,
               cap_applied=cap)
    return out


def rank_rows(rows, points_key: str, rank_key: str):
    """Assign 1-based integer ranks in place, highest points first.

    Ties break by ``pts_molecular`` desc, then ``pts_rarity`` desc, then
    ``pts_gene_artifact`` desc (cleaner genes first), then ``chrom``/``pos`` for determinism —
    so the ordering is reproducible across runs and platforms.
    """
    def key(r):
        return (-(_num(r.get(points_key)) or 0.0),
                -(_num(r.get("pts_molecular")) or 0.0),
                -(_num(r.get("pts_rarity")) or 0.0),
                -(_num(r.get("pts_gene_artifact")) or 0.0),
                _s(r.get("chrom")), _int(r.get("pos")) or 0, _s(r.get("ref")), _s(r.get("alt")),
                _s(r.get("trio_id")))
    for i, r in enumerate(sorted(rows, key=key), start=1):
        r[rank_key] = i
    return rows
