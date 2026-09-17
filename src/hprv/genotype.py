"""Per-sample genotype QC helpers for GATK-refined trio VCFs.

Trusts the REFINED PP-derived GQ (cyvcf2 ``gt_quals`` reads the GQ field, which
CalculateGenotypePosteriors overwrites from PP). Allele balance is derived from AD
(it is not a native FORMAT field). See docs/inheritance_and_genotype_qc.md.
"""

from __future__ import annotations

from dataclasses import dataclass

# cyvcf2 gt_types encoding
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3

# GRCh38 pseudoautosomal regions (standard genome coordinates).
PAR_X = ((10001, 2781479), (155701383, 156030895))
PAR_Y = ((10001, 2781479), (56887903, 57217415))


@dataclass
class GtThresholds:
    min_gq: int = 20
    min_dp: int = 10
    denovo_min_dp: int = 20
    # the de novo depth floor on a HEMIZYGOUS site (a male's non-PAR chrX/chrY): a haploid
    # chromosome carries half the autosomal depth, so the diploid floor is twice as strict there
    denovo_min_dp_hemizygous: int = 10
    het_ab_min: float = 0.25
    het_ab_max: float = 0.75
    homalt_ab_min: float = 0.90
    homref_ab_max: float = 0.10
    parent_max_alt_ad: int = 1
    parent_min_dp: int = 10

    @classmethod
    def from_config(cls, cfg, get):
        g = "filters.genotype_qc."
        d = "filters.denovo."
        return cls(
            min_gq=int(get(cfg, g + "min_gq", 20)),
            min_dp=int(get(cfg, g + "min_dp", 10)),
            denovo_min_dp=int(get(cfg, g + "denovo_min_dp", 20)),
            denovo_min_dp_hemizygous=int(get(cfg, g + "denovo_min_dp_hemizygous", 10)),
            het_ab_min=float(get(cfg, g + "het_ab_min", 0.25)),
            het_ab_max=float(get(cfg, g + "het_ab_max", 0.75)),
            homalt_ab_min=float(get(cfg, g + "homalt_ab_min", 0.90)),
            homref_ab_max=float(get(cfg, g + "homref_ab_max", 0.10)),
            parent_max_alt_ad=int(get(cfg, d + "parent_max_alt_ad", 1)),
            parent_min_dp=int(get(cfg, d + "parent_min_dp", 10)),
        ).validated()

    def validated(self) -> "GtThresholds":
        """Refuse a band no genotype could pass. An inverted het band (`het_ab_min` above
        `het_ab_max`) or a percent-scale value (`het_ab_min: 25`) reduced Step 5 to zero calls in
        every inherited mode, with exit 0 and nothing in the audit to say why."""
        for k in ("het_ab_min", "het_ab_max", "homalt_ab_min", "homref_ab_max"):
            v = getattr(self, k)
            if not 0.0 <= v <= 1.0:
                raise ValueError(f"filters.genotype_qc.{k} = {v}: allele balance is a fraction in [0, 1], not a percentage")
        if self.het_ab_min >= self.het_ab_max:
            raise ValueError(f"filters.genotype_qc.het_ab_min ({self.het_ab_min}) must be below het_ab_max ({self.het_ab_max}); no het could pass this band")
        if (self.min_dp < 1 or self.denovo_min_dp < 1 or self.denovo_min_dp_hemizygous < 1
                or self.parent_min_dp < 1 or self.min_gq < 0):
            raise ValueError("filters.genotype_qc: min_dp / denovo_min_dp / denovo_min_dp_hemizygous / parent_min_dp must be >= 1 and min_gq >= 0")
        return self


def _int(x):
    """Coerce a cyvcf2 numpy scalar to int; treat negatives/None as missing."""
    if x is None:
        return None
    try:
        v = int(x)
    except (TypeError, ValueError):
        return None
    return None if v < 0 else v


def gq(v, i):
    return _int(v.gt_quals[i])


def dp(v, i):
    """Per-sample depth. Reads FORMAT/DP directly when cyvcf2's `gt_depths` is missing.

    cyvcf2 derives `gt_depths` from the AD sums whenever the record carries an AD field, so a
    sample with NO AD — a GATK ref-block-derived `0/0` parent, exactly the shape of a parent
    under a de novo — reads -1 even though its own DP is present. Every `sample_qc` limb tests
    depth before allele balance, so that sample failed CLOSED on depth: the "AD limbs fail open"
    analysis never got to run, and a de novo with a ref-block parent was silently dropped rather
    than called with `parent_ad_unmeasured`. Fall back to the sample's FORMAT/DP.
    """
    d = _int(v.gt_depths[i])
    if d is not None:
        return d
    try:
        arr = v.format("DP")
    except (AttributeError, KeyError, TypeError, ValueError):
        return None
    if arr is None:
        return None
    try:
        return _int(arr[i][0])
    except (IndexError, TypeError):
        return None


def alt_ad(v, i):
    return _int(v.gt_alt_depths[i])


def ref_ad(v, i):
    return _int(v.gt_ref_depths[i])


def allele_balance(v, i):
    r, a = ref_ad(v, i), alt_ad(v, i)
    if r is None or a is None:
        return None
    tot = r + a
    return (a / tot) if tot > 0 else None


def sex_contig(chrom: str):
    """'X' or 'Y' when `chrom` is a sex-chromosome contig in any spelling, else None.

    Recognises the primary contig (`chrX`, `X`, `chrY`, `Y`) AND GRCh38's unlocalised / alt
    scaffolds of it (`chrY_KI270740v1_random`, `chrX_KI270880v1_alt`): those are sex-chromosome
    sequence too, hemizygous in males, and reading them as autosomes emitted a male's het calls
    there as `dominant` (robustness audit #99). The PAR coordinates apply to the PRIMARY contig
    only; a scaffold is treated as non-PAR.
    """
    c = chrom[3:] if chrom.startswith("chr") else chrom
    head = c.split("_", 1)[0]
    return head if head in ("X", "Y") else None


def _primary_sex_contig(chrom: str):
    """The bare 'X'/'Y' for the PRIMARY contig only (a scaffold has no PAR)."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return c if c in ("X", "Y") else None


def in_par_x(v) -> bool:
    if _primary_sex_contig(v.CHROM) != "X":
        return False
    return any(lo <= v.POS <= hi for lo, hi in PAR_X)


def is_x_nonpar(v) -> bool:
    return sex_contig(v.CHROM) == "X" and not in_par_x(v)


def in_par_y(v) -> bool:
    return _primary_sex_contig(v.CHROM) == "Y" and any(lo <= v.POS <= hi for lo, hi in PAR_Y)


def is_y_nonpar(v) -> bool:
    return sex_contig(v.CHROM) == "Y" and not in_par_y(v)


def is_sex_nonpar(v) -> bool:
    """Non-PAR chrX or chrY — hemizygous in males; het calls there are a QC red flag."""
    return is_x_nonpar(v) or is_y_nonpar(v)


# --- pre-refinement genotype likelihoods (FORMAT/PL) --------------------------------------------
# GATK's genotype refinement (CalculateGenotypePosteriors) overwrites GT/GQ from a posterior that
# folds in a PEDIGREE prior and a population prior, but it leaves the original PL in place. That
# prior is DIPLOID and Mendelian: it is right on the autosomes and wrong on a haploid chromosome.
# On non-PAR chrY it has two concrete effects on a GMKF trio: a father-son hemizygous alt site
# (both truly 1/1) is "explained" by imputing the mother a 0/1 she has no reads for, and where the
# mother's own mismapped reads say 0/0 confidently the son's all-alt 1/1 is pushed to 0/1 (a 1/1
# child of a 0/0 mother is a de novo under the diploid prior, ~1e-8) — after which the male-het
# rule reads the true hemizygous call as a mapping artifact and drops it. The hemizygous call
# below therefore reads the LIKELIHOODS when the record carries them and the refined GT only as a
# fallback, and every consumer records which one decided.

HEMI_ALT, HEMI_REF, HEMI_MIXED, HEMI_NOCALL = "alt", "ref", "mixed", "nocall"
_GT_TO_HEMI = {HOM_ALT: HEMI_ALT, HOM_REF: HEMI_REF, HET: HEMI_MIXED}


def pl(v, i):
    """The sample's FORMAT/PL values as a list of ints (None where a value is missing), or None
    when the record carries no PL for that sample at all."""
    try:
        arr = v.format("PL")
    except (AttributeError, KeyError, TypeError, ValueError):
        return None
    if arr is None:
        return None
    try:
        row = arr[i]
    except (IndexError, TypeError):
        return None
    try:
        out = [_int(x) for x in row]     # cyvcf2 missing / vector-end sentinels are negative
    except TypeError:
        return None
    return out if any(x is not None for x in out) else None


def gt_from_pl(pls):
    """The biallelic diploid genotype class the PLs favour — HOM_REF / HET / HOM_ALT (the
    argmin over the first three values) — or None when the PLs are absent, incomplete, or
    tied at the minimum (flat `0,0,0` PLs carry no read evidence at all)."""
    if not pls or len(pls) < 3:
        return None
    vals = pls[:3]
    if any(x is None for x in vals):
        return None
    m = min(vals)
    if vals.count(m) != 1:
        return None
    return (HOM_REF, HET, HOM_ALT)[vals.index(m)]


def gq_from_pl(pls):
    """Pre-refinement genotype quality: the second-smallest PL minus the smallest, capped at 99
    (exactly how GATK derives GQ from PL before any prior is applied); None without >= 2 values."""
    vals = sorted(x for x in (pls or [])[:3] if x is not None)
    if len(vals) < 2:
        return None
    return min(99, vals[1] - vals[0])


def hemi_call(v, i):
    """Zygosity of a HAPLOID (hemizygous) site as the pre-refinement likelihoods see it.

    Returns ``(call, source, refined)``: `call` and `refined` are each one of `alt` / `ref` /
    `mixed` (both alleles have read support — on a haploid chromosome that is the mismapping /
    paralogue signature, never a genotype) / `nocall`; `source` is `pl` when FORMAT/PL decided
    the call and `gt` when the refined GT had to. A refined no-call stays a no-call: the
    likelihood read re-labels among CALLED genotypes only (the diploid-prior distortion moves a
    hemizygous 1/1 to 0/1 and imputes a parent; it does not manufacture data where a caller
    declined to call).
    """
    g = v.gt_types[i]
    if g not in _GT_TO_HEMI:
        return HEMI_NOCALL, "gt", HEMI_NOCALL
    refined = _GT_TO_HEMI[g]
    from_pl = gt_from_pl(pl(v, i))
    if from_pl is None:
        return refined, "gt", refined
    return _GT_TO_HEMI[from_pl], "pl", refined


def hemi_gq(v, i):
    """GQ for a decision on a hemizygous-model site: derived from PL when present (the
    prior-free quality — valid for a diploid member too, it is simply GATK's GQ before any prior),
    else the refined GQ."""
    q = gq_from_pl(pl(v, i))
    return q if q is not None else gq(v, i)


def hemi_qc_reason(v, i, thr: GtThresholds, kind: str):
    """WHY a sample fails the QC band for a HAPLOID site, or None when it passes — the same
    reason vocabulary as sample_qc_reason. kind: 'alt' (a hemizygous alternate carrier: DP >=
    min_dp, GQ >= min_gq on the pre-refinement quality, AB >= homalt_ab_min — fails CLOSED
    without AD, like `sample_qc(..., "hom_alt")`) or 'clean' (a confidently reference parent under
    a hemizygous de novo: DP >= parent_min_dp, alt AD <= parent_max_alt_ad, AB <= homref_ab_max —
    the AD limb fails OPEN, exactly like `sample_qc(..., "clean_parent")`, and
    `sample_qc_ad_measured(v, i, "clean_parent")` is the witness)."""
    q = hemi_gq(v, i)
    if q is None or q < thr.min_gq:
        return "gq"
    d = dp(v, i)
    need_dp = thr.parent_min_dp if kind == "clean" else thr.min_dp
    if d is None or d < need_dp:
        return "dp"
    ab = allele_balance(v, i)
    if kind == "alt":
        if ab is None:
            return "ab_missing"
        return None if ab >= thr.homalt_ab_min else "ab_low"
    if kind == "clean":
        a = alt_ad(v, i)
        ok = (a is None or a <= thr.parent_max_alt_ad) and (ab is None or ab <= thr.homref_ab_max)
        return None if ok else "alt_reads"
    if kind == "ref":
        # a confident hemizygous REFERENCE (the non-carrier band, mirroring the diploid `hom_ref`
        # band: AB <= homref_ab_max, the AD limb failing open) — looser than `clean`
        return None if (ab is None or ab <= thr.homref_ab_max) else "alt_reads"
    return "kind"


def hemi_qc(v, i, thr: GtThresholds, kind: str) -> bool:
    """True when the sample passes the haploid QC band for `kind` (see hemi_qc_reason)."""
    return hemi_qc_reason(v, i, thr, kind) is None


def sample_qc_reason(v, i, thr: GtThresholds, kind: str, gq_value=None):
    """WHY a sample fails the QC band for `kind`, or None when it passes.

    kind: 'het' | 'hom_alt' | 'hom_ref' | 'denovo_child' | 'clean_parent'. The reason is one of
    `gq`, `dp`, `ab_missing` (the carrier bands fail CLOSED without AD), `ab_low`, `ab_high`,
    or `alt_reads` (a hom-ref / clean-parent band failed because the sample HAS alt reads — the
    parental-mosaicism / contamination shape, which Step 5 counts apart from a quality failure).
    `gq_value` substitutes the refined GQ — Step 5 passes the PL-derived quality (`hemi_gq`) on a
    male's non-PAR chrX, where the refinement's diploid pedigree prior is invalid (see below).
    """
    q = gq(v, i) if gq_value is None else gq_value
    if q is None or q < thr.min_gq:
        return "gq"
    d = dp(v, i)
    need_dp = thr.denovo_min_dp if kind == "denovo_child" else (
        thr.parent_min_dp if kind == "clean_parent" else thr.min_dp)
    if d is None or d < need_dp:
        return "dp"
    ab = allele_balance(v, i)
    if kind in ("het", "denovo_child"):
        if ab is None:
            return "ab_missing"
        if ab < thr.het_ab_min:
            return "ab_low"
        if ab > thr.het_ab_max:
            return "ab_high"
        return None
    if kind == "hom_alt":
        if ab is None:
            return "ab_missing"
        return None if ab >= thr.homalt_ab_min else "ab_low"
    if kind == "hom_ref":
        return None if (ab is None or ab <= thr.homref_ab_max) else "alt_reads"
    if kind == "clean_parent":
        # parent must be hom-ref AND essentially free of alt reads
        a = alt_ad(v, i)
        ok = (a is None or a <= thr.parent_max_alt_ad) and (ab is None or ab <= thr.homref_ab_max)
        return None if ok else "alt_reads"
    return "kind"


def sample_qc(v, i, thr: GtThresholds, kind: str, gq_value=None) -> bool:
    """True when the sample passes the QC band for `kind` (see sample_qc_reason)."""
    return sample_qc_reason(v, i, thr, kind, gq_value=gq_value) is None


def sample_qc_ad_measured(v, i, kind: str) -> bool:
    """Was the ALLELE-DEPTH evidence behind a ``hom_ref``/``clean_parent`` pass actually present?

    **This exists because the AD limbs fail OPEN while the others fail closed**, and the asymmetry
    is the whole problem: ``het``/``hom_alt``/``denovo_child`` require ``ab is not None`` and so
    DROP a carrier when AD is missing, while ``hom_ref``/``clean_parent`` return True when it is
    missing and so AFFIRM a non-carrier. The same absent measurement decides both ways.

    It is not hypothetical: a GATK ref-block-derived ``0/0`` parent carries ``GT:DP:GQ:MIN_DP:PL``
    with no AD at all, which is exactly the shape of a parent at a site where the child is het.
    And ``clean_parent`` is the ONLY evidence that a compound-het pair is in **trans**.

    Hard-failing instead would drop the ordinary ref-block case wholesale, so the pass stands
    (never-drop) and callers mark it instead — ``origin_unverified`` previously fired only when
    the test FAILED, leaving a vacuous pass silently indistinguishable from a measured one.
    """
    if kind not in ("hom_ref", "clean_parent"):
        return True
    if allele_balance(v, i) is not None:
        return True
    return kind == "clean_parent" and alt_ad(v, i) is not None
