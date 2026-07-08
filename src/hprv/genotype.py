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
            het_ab_min=float(get(cfg, g + "het_ab_min", 0.25)),
            het_ab_max=float(get(cfg, g + "het_ab_max", 0.75)),
            homalt_ab_min=float(get(cfg, g + "homalt_ab_min", 0.90)),
            homref_ab_max=float(get(cfg, g + "homref_ab_max", 0.10)),
            parent_max_alt_ad=int(get(cfg, d + "parent_max_alt_ad", 1)),
            parent_min_dp=int(get(cfg, d + "parent_min_dp", 10)),
        )


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
    return _int(v.gt_depths[i])


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


def in_par_x(v) -> bool:
    chrom = v.CHROM.replace("chr", "")
    if chrom != "X":
        return False
    return any(lo <= v.POS <= hi for lo, hi in PAR_X)


def is_x_nonpar(v) -> bool:
    return v.CHROM.replace("chr", "") == "X" and not in_par_x(v)


def in_par_y(v) -> bool:
    return v.CHROM.replace("chr", "") == "Y" and any(lo <= v.POS <= hi for lo, hi in PAR_Y)


def is_y_nonpar(v) -> bool:
    return v.CHROM.replace("chr", "") == "Y" and not in_par_y(v)


def is_sex_nonpar(v) -> bool:
    """Non-PAR chrX or chrY — hemizygous in males; het calls there are a QC red flag."""
    return is_x_nonpar(v) or is_y_nonpar(v)


def sample_qc(v, i, thr: GtThresholds, kind: str) -> bool:
    """kind: 'het' | 'hom_alt' | 'hom_ref' | 'denovo_child' | 'clean_parent'."""
    q = gq(v, i)
    if q is None or q < thr.min_gq:
        return False
    d = dp(v, i)
    need_dp = thr.denovo_min_dp if kind == "denovo_child" else (
        thr.parent_min_dp if kind == "clean_parent" else thr.min_dp)
    if d is None or d < need_dp:
        return False
    ab = allele_balance(v, i)
    if kind in ("het", "denovo_child"):
        return ab is not None and thr.het_ab_min <= ab <= thr.het_ab_max
    if kind == "hom_alt":
        return ab is not None and ab >= thr.homalt_ab_min
    if kind == "hom_ref":
        return ab is None or ab <= thr.homref_ab_max
    if kind == "clean_parent":
        # parent must be hom-ref AND essentially free of alt reads
        a = alt_ad(v, i)
        return (a is None or a <= thr.parent_max_alt_ad) and (ab is None or ab <= thr.homref_ab_max)
    return False
