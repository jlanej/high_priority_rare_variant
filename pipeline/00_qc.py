#!/usr/bin/env python3
"""Pipeline Step 0: per-trio QC — advisory pass/flag list (garbage-in guard, human-reviewed).

Computes Mendelian-error rate, chrX-inferred sex (all three members), per-member no-call rate,
and contamination per trio and folds them into an `overall_pass` flag. This is ADVISORY: flags
are surfaced to the reviewer (xlsx QC sheet + IGV sample_qc.tsv) but no step auto-excludes a
flagged trio (see the overall_pass comment below).

Mostly self-contained (no extra resources): for each trio computes five gates and
flags trios that fail any of them:
  * Mendelian-error rate — a proxy for sample swaps / contamination. NOT for a father/mother
    label swap: the rule is symmetric under exchanging the parents, so that swap is invisible
    to it by construction — which is what the next gate is for.
  * chrX-heterozygosity sex inference for all three members. The PROBAND's is checked against
    the PED sex when the trios file stated one (`sex_match`; the PED is canonical — Step 5 keeps
    it and flags the disagreement, this step never overrides it) and is the fallback when it did
    not (`sex_source`, via ped.resolve_child_sex). BOTH PARENTS' are checked against their roles
    (the PED always states those): a male in the mother slot is the one direct detector of
    transposed parents. The inference is a HEURISTIC with a callset-dependent cutoff
    (`qc.x_het_male_max`): a diploid-called male's non-PAR chrX het ratio often sits well above
    the 0.10 default, so the report carries every member's raw ratio (`x_het_ratio`,
    `dad_x_het_ratio`, `mom_x_het_ratio`) for calibration, and a cohort where most fathers read
    female is reported as a miscalibrated cutoff, not as mass transposition.
  * chrY COVERAGE sex evidence for the proband — secondary and AUDIT-ONLY (it decides nothing:
    the trios file is canonical and the chrX inference is what fills in an unknown). Anchored
    on the FATHER's own hemizygous chrY sites (a known male) with the MOTHER as the in-trio
    female control: the proband's median depth there, autosome-normalised, is ~1x the father's
    for a son and ~0 for a daughter, which is far more robust on a diploid-called callset than
    the chrX het ratio. Reported (`y_inferred_sex`, `y_cov_ratio`, `y_covered_frac`,
    `sex_match_y`, `xy_agree`, `mom_y_cov_ratio`, `y_flag`), WARNed on disagreement with the
    pedigree or with the chrX inference, audited — and never folded into `overall_pass`. A few
    mismapped reads cannot read a female as male: the call is a median over
    >= `qc.y_min_anchor_sites` sites, and a trio whose mother herself shows chrY coverage is
    refused a call (`mother_y_coverage`) rather than trusted. Beside it, each member's RAW chrY
    coverage — the median depth over every non-PAR chrY record, as a fraction of half the
    member's autosomal depth (`*_y_cov_haploid`; a male ~1, a female ~0) — checks BOTH PARENTS'
    sex against their roles (`dad_y_sex`, `mom_y_sex`, `parent_sex_flag_y`): the one direct
    detector of a transposed pair that does not depend on a chrX cutoff, and the fallback
    basis for the proband's call when the father yields no usable anchors (`y_basis=raw`).
  * Per-member genotype no-call rate on the scanned sites. A jointly genotyped trio carries an
    affirmative `0/0` for a non-carrier parent; a merge of single-sample callsets carries `./.`
    there instead, and Step 5 then loses every de novo and marks every inherited call
    `origin_unverified` — while the MIE denominator, which skips no-call sites, reports a
    CLEANER trio. The rate makes that input shape visible before it silently empties a run.
  * Contamination — verifyBamID FREEMIX if a directory of ``*.selfSM`` files is
    configured (resources.selfsm_dir), else a VCF-only CHARR estimate (reference-read
    fraction at high-quality hom-ALT SNV sites). Mirrors the group's DNM freemix QC.

Trio kid/dad/mom roles come from upstream peddy; this step is the guard against the
less-well-curated trios. See docs/inheritance_and_genotype_qc.md.

The Mendelian-error / CHARR scan and the chrX sex scan are each CAPPED at ``qc.max_sites``
QC-passing sites (default 200,000; ``--max-sites`` overrides): on WGS the MIE rate is therefore
measured on a prefix of the genome, not genome-wide — adequate as a swap/contamination detector,
but a localised MIE cluster (UPD/CNV) outside that prefix is invisible to it.

Usage:
  00_qc.py --manifest trios.tsv --config config.yaml --out qc_report.tsv [--max-sites N]
"""
from __future__ import annotations

import argparse
import csv
import itertools
import statistics
import sys

from cyvcf2 import VCF

from hprv import audit
from hprv import contamination as C
from hprv import genotype as G
from hprv.config import get, load_config
from hprv.ped import parse_ped, resolve_child_sex

HOM_REF, HET, HOM_ALT = G.HOM_REF, G.HET, G.HOM_ALT


def mendelian_violation(gc, gd, gm) -> bool:
    """True if child genotype is impossible given parents (biallelic autosomal)."""
    if gc == HOM_ALT:      # needs an alt from each parent
        return gd == HOM_REF or gm == HOM_REF
    if gc == HOM_REF:      # needs a ref from each parent
        return gd == HOM_ALT or gm == HOM_ALT
    if gc == HET:          # needs one alt and one ref available
        return (gd == HOM_REF and gm == HOM_REF) or (gd == HOM_ALT and gm == HOM_ALT)
    return False


def infer_sex(x_het, x_hom, cutoff, min_sites):
    """chrX-heterozygosity sex call: '1' male, '2' female, None when too few informative sites."""
    total = x_het + x_hom
    if total < min_sites or total == 0:
        return None
    return "1" if (x_het / total) < cutoff else "2"


def _pick_contig(seqnames, base):
    for want in (f"chr{base}", base):
        if want in seqnames:
            return want
    return None


def _pick_x_contig(seqnames):
    return _pick_contig(seqnames, "X")


def _region_iter(vcf, contig):
    """Indexed jump to `contig`, with scan_sex's fallback: cyvcf2.VCF.__call__ is a generator
    function, so the tabix load fires on the first next() — force it here so an unindexed VCF
    falls back to a full scan instead of raising through the run."""
    try:
        it = vcf(contig)
        it = itertools.chain([next(it)], it)
    except StopIteration:
        it = iter(())                     # the contig is in the header but has no records
    except Exception:
        it = vcf                          # unindexed: full scan; callers filter by contig
    return it


def scan_y(vcf_path, kid, dad, mom, thr, max_y, min_reads):
    """chrY coverage evidence for the proband and mother, ANCHORED on the father's hemizygous
    sites — a dedicated pass over the chrY records.

    An anchor is a FILTER-passing non-PAR chrY record where the father is a confident hemizygous
    ALT carrier on the pre-refinement likelihoods (`genotype.hemi_call` + `hemi_qc('alt')`: DP,
    GQ and AB >= homalt_ab_min). Those are the sites this callset demonstrably maps Y reads to
    in a male, and they exclude the X-transposed / ampliconic sites where even the father reads
    het. At each anchor the proband's and the mother's depth (missing = 0) is what a Y
    chromosome looks like: a son matches his father, a daughter — and the mother — carries
    nothing. Genotypes are never read for the proband or mother here: on chrY the refined GT is
    the diploid pedigree prior's opinion (an imputed maternal 0/1, a daughter's imputed 1/1),
    not evidence. Returns None when the VCF has no chrY contig or a member is absent.
    """
    vcf = VCF(vcf_path, strict_gt=True)
    idx = {s: i for i, s in enumerate(vcf.samples)}
    y = _pick_contig(vcf.seqnames, "Y")
    if y is None or any(s not in idx for s in (kid, dad, mom)):
        vcf.close()
        return None
    c, d, m = idx[kid], idx[dad], idx[mom]
    out = {"anchors": 0, "kid_dp": [], "dad_dp": [], "mom_dp": [], "kid_cov": 0, "mom_cov": 0,
           # RAW coverage: every member's depth at EVERY passing non-PAR chrY record, whoever
           # drove the record. A male has reads at (nearly) all of them, a female at few — the
           # basis for checking the PARENTS' sex against their roles, which the father-anchored
           # statistic above cannot do (it presumes the father slot holds a male).
           "records": 0, "kid_all": [], "dad_all": [], "mom_all": []}
    for v in _region_iter(vcf, y):
        if not G.is_y_nonpar(v) or v.FILTER:
            continue
        kd, dd, md = (G.dp(v, i) or 0 for i in (c, d, m))
        if not max_y or out["records"] < max_y:
            out["records"] += 1
            out["kid_all"].append(kd)
            out["dad_all"].append(dd)
            out["mom_all"].append(md)
        if G.hemi_call(v, d)[0] != G.HEMI_ALT or not G.hemi_qc(v, d, thr, "alt"):
            continue
        out["anchors"] += 1
        out["kid_dp"].append(kd)
        out["dad_dp"].append(dd)
        out["mom_dp"].append(md)
        out["kid_cov"] += kd >= min_reads
        out["mom_cov"] += md >= min_reads
        if max_y and out["anchors"] >= max_y and out["records"] >= max_y:
            break
    vcf.close()
    return out


def y_cov_ratio(member_y_med, member_auto_med, dad_y_med, dad_auto_med):
    """A member's chrY depth at the father's hemizygous sites as a fraction of the father's,
    each normalised by its own autosomal median depth so library depth cancels: a son reads
    ~1, a daughter and the mother ~0. None when any denominator is unavailable."""
    if member_y_med is None or not member_auto_med or not dad_auto_med or not dad_y_med:
        return None
    return (member_y_med / member_auto_med) / (dad_y_med / dad_auto_med)


def y_cov_haploid(member_y_all_med, member_auto_med):
    """A member's median depth over every non-PAR chrY record as a fraction of HALF the member's
    autosomal median — the haploid expectation — so a male reads ~1 (his single Y at half his
    autosomal depth) and a female ~0, on the same scale as the father-normalised anchor ratio.
    None when either median is unavailable."""
    if member_y_all_med is None or not member_auto_med:
        return None
    return member_y_all_med / (0.5 * member_auto_med)


def sex_from_y_cov(ratio, male_min, female_max):
    """'1' / '2' / None from a haploid-scaled chrY coverage ratio (the same bands as the
    proband's call)."""
    if ratio is None:
        return None
    if ratio >= male_min:
        return "1"
    if ratio <= female_max:
        return "2"
    return None


def parent_y_diagnosis(dad_y_sex, mom_y_sex) -> str:
    """What a parent-sex mismatch on chrY COVERAGE can mean — the same taxonomy as the chrX
    diagnosis, without the cutoff caveat: chrY coverage is near-binary, so a father with no Y
    coverage is a female (or a Y-less) sample in the father slot, not a calibration artifact."""
    if dad_y_sex == "2" and mom_y_sex == "1":
        return ("the father slot has NO chrY coverage and the mother slot has FULL chrY coverage: "
                "the trios-file father/mother roles for this trio are transposed (every "
                "parent-of-origin call would be inverted, and Step 5's Y-linked model anchors on "
                "the wrong parent)")
    if dad_y_sex == "2":
        return ("the father slot has no chrY coverage while the mother slot reads female too: a "
                "sample swap/mislabel in the father slot, or a father with no Y (46,XX male, "
                "Y loss) — the Y-linked model has no transmitter to anchor on for this trio")
    return ("the mother slot has chrY coverage while the father slot reads male too: a sample "
            "swap/mislabel in the mother slot, or a male sample in it — her reads are not a "
            "female control, so no chrY sex call is made for the proband")


def infer_sex_y(kid_ratio, mom_ratio, n_anchor, min_anchor, male_min, female_max):
    """chrY-coverage sex call for the proband -> ('1' | '2' | None, flag).

    Refuses to call — and says why — whenever the evidence is thin or the in-trio female
    control failed: too few anchors, no ratio, a mother whose own ratio exceeds `female_max`
    (contamination, transposed parents or a sample swap — whatever it is, the trio's chrY reads
    cannot be trusted to separate the sexes), or a proband in the indeterminate band. A female
    with a handful of mismapped reads sits far below `female_max` because the statistic is a
    MEDIAN over the anchors, not a count of sites with any read.
    """
    if n_anchor < min_anchor:
        return None, "father_y_anchors_low"
    if kid_ratio is None:
        return None, "y_ratio_unavailable"
    if mom_ratio is not None and mom_ratio > female_max:
        return None, "mother_y_coverage"
    if kid_ratio >= male_min:
        return "1", ("" if mom_ratio is not None else "mother_y_control_unavailable")
    if kid_ratio <= female_max:
        return "2", ("" if mom_ratio is not None else "mother_y_control_unavailable")
    return None, "y_cov_indeterminate"


def scan_sex(vcf_path, sample_id, thr, max_x):
    """chrX non-PAR het/hom counts for one member, as a DEDICATED pass.

    chrX sorts after all autosomes, so an autosomal MIE cap in the main pass would
    otherwise starve sex inference on WGS trios. Uses the index to jump straight to
    chrX; falls back to a full scan (filtered to chrX) when the VCF is unindexed.

    Only FILTER-passing records (PASS or `.`, the same convention as Step 1 and Step 5) with a
    fully called genotype count: a filtered chrX record is enriched for exactly the mapping
    artifacts that render a hemizygous male as het, and cyvcf2's default reads a half-called
    `1/.` as HET — both inflate the male het ratio toward the female side of the cutoff.
    """
    vcf = VCF(vcf_path, strict_gt=True)      # a half-called 1/. is UNKNOWN, never HET
    ci = {s: i for i, s in enumerate(vcf.samples)}.get(sample_id)
    x = _pick_x_contig(vcf.seqnames)
    if ci is None or x is None:
        vcf.close()
        return 0, 0
    try:
        it = vcf(x)                       # indexed region jump (preferred)
        # cyvcf2.VCF.__call__ is a GENERATOR function, so `vcf(x)` returns without running its
        # body — the tabix load and its `assert self.idx != NULL` fire on the first next(), i.e.
        # inside the for-loop below and OUTSIDE this try. Without forcing the first element here
        # the fallback is dead code and an unindexed VCF raises AssertionError straight through
        # qc_trio/main into `set -euo pipefail`, killing the whole run. Unindexed trio VCFs are a
        # supported input, and Step 0 runs before anything indexes them.
        it = itertools.chain([next(it)], it)
    except StopIteration:
        it = iter(())                     # chrX in the header but no records there
    except Exception:
        it = vcf                          # unindexed: full scan, filtered below
    x_het = x_hom = 0
    for v in it:
        if v.CHROM.replace("chr", "") != "X" or G.in_par_x(v) or len(v.ALT) != 1:
            continue
        if v.FILTER:                  # cyvcf2: None for PASS and '.'; anything else is filtered
            continue
        gq, dp = G.gq(v, ci), G.dp(v, ci)
        if gq is None or gq < thr.min_gq or dp is None or dp < thr.min_dp:
            continue
        gt = v.gt_types[ci]
        if gt == HET:
            x_het += 1
        elif gt == HOM_ALT:
            x_hom += 1
        if max_x and (x_het + x_hom) >= max_x:
            break
    vcf.close()
    return x_het, x_hom


def qc_trio(vcf_path, ped, thr, max_sites, sex_cutoff=0.10, sex_min_sites=20, y_reads_min_dp=3):
    vcf = VCF(vcf_path)
    samples = {s: i for i, s in enumerate(vcf.samples)}
    for role in ("child", "father", "mother"):
        if ped[role] not in samples:
            vcf.close()
            return None
    c, d, m = samples[ped["child"]], samples[ped["father"]], samples[ped["mother"]]

    # sex inference runs as its own chrX pass so the autosomal MIE cap can't disable it
    x_het, x_hom = scan_sex(vcf_path, ped["child"], thr, max_x=(max_sites or 0))
    # ...and for the PARENTS. Their sex is the one thing the PED asserts with certainty
    # (father = 1, mother = 2), so this check is reachable even when the trios file states no
    # proband sex, and it is the one direct detector of a transposed mother/father, which the
    # Mendelian-error rate cannot see (mendelian_violation is symmetric in gd/gm). Their raw
    # het ratios are also the CALIBRATION set for qc.x_het_male_max: every father is a known
    # male and every mother a known female, so the cutoff belongs in the gap between the two.
    px = {role: scan_sex(vcf_path, ped[role], thr, max_x=(max_sites or 0))
          for role in ("father", "mother")}

    # ...and chrY: the proband's and mother's depth at the father's hemizygous sites (scan_y).
    ys = scan_y(vcf_path, ped["child"], ped["father"], ped["mother"], thr,
                max_y=(max_sites or 0), min_reads=y_reads_min_dp)

    considered = errors = 0
    n_records = 0                            # biallelic autosomal records seen (pre-QC)
    nocall = {"kid": 0, "dad": 0, "mom": 0}  # per-member `./.` among them
    auto_dp = {"kid": [], "dad": [], "mom": []}   # per-member depth at the scanned autosomal sites
    cref = {"kid": 0, "dad": 0, "mom": 0}   # CHARR: ref reads at hom-alt SNV sites
    cdp = {"kid": 0, "dad": 0, "mom": 0}    # CHARR: total (ref+alt) reads there
    for v in vcf:
        chrom = v.CHROM.replace("chr", "")
        if chrom in ("X", "Y", "MT", "M"):  # X counted in scan_sex; Y/MT out of scope
            continue
        if len(v.ALT) != 1:  # biallelic only for the simple MIE rule
            continue
        # CHARR: accumulate reference reads at each member's high-quality hom-alt SNV sites
        if len(v.REF) == 1 and len(v.ALT[0]) == 1:
            for role, i in (("kid", c), ("dad", d), ("mom", m)):
                if v.gt_types[i] == HOM_ALT:
                    gq, dp = G.gq(v, i), G.dp(v, i)
                    ra, aa = G.ref_ad(v, i), G.alt_ad(v, i)
                    if (gq and gq >= thr.min_gq and dp and dp >= thr.min_dp
                            and ra is not None and aa is not None and (ra + aa) > 0):
                        cref[role] += ra
                        cdp[role] += ra + aa
        gc, gd, gm = v.gt_types[c], v.gt_types[d], v.gt_types[m]
        n_records += 1
        for role, g in (("kid", gc), ("dad", gd), ("mom", gm)):
            if g == G.UNKNOWN:
                nocall[role] += 1
        # the autosomal depth baseline for the chrY coverage ratio (pre-QC, every member, capped
        # with the scan so it stays a bounded sample)
        if not max_sites or n_records <= max_sites:
            for role, i in (("kid", c), ("dad", d), ("mom", m)):
                _d = G.dp(v, i)
                if _d is not None:
                    auto_dp[role].append(_d)
        if G.UNKNOWN in (gc, gd, gm):
            continue                          # a no-call site is out of the MIE denominator
        if any(G.gq(v, i) is None or G.gq(v, i) < thr.min_gq for i in (c, d, m)):
            continue
        if any(G.dp(v, i) is None or G.dp(v, i) < thr.min_dp for i in (c, d, m)):
            continue
        considered += 1
        if mendelian_violation(gc, gd, gm):
            errors += 1
        if max_sites and considered >= max_sites:
            break
    vcf.close()

    mie_rate = (errors / considered) if considered else None
    x_total = x_het + x_hom
    x_het_ratio = (x_het / x_total) if x_total else None
    # only call sex with enough informative chrX sites, else leave it unknown (fail-soft)
    inferred = infer_sex(x_het, x_hom, sex_cutoff, sex_min_sites)

    def _ratio(h, a):
        return (h / (h + a)) if (h + a) else None

    def _med(xs):
        return statistics.median(xs) if xs else None

    auto_med = {r: _med(auto_dp[r]) for r in ("kid", "dad", "mom")}
    if ys is not None:
        ys = dict(ys, kid_med=_med(ys["kid_dp"]), dad_med=_med(ys["dad_dp"]),
                  mom_med=_med(ys["mom_dp"]), kid_all_med=_med(ys["kid_all"]),
                  dad_all_med=_med(ys["dad_all"]), mom_all_med=_med(ys["mom_all"]))

    return {
        "n_sites": considered, "mie_errors": errors, "mie_rate": mie_rate,
        "x_sites": x_total, "x_het_ratio": x_het_ratio, "inferred_sex": inferred,
        "dad_sex": infer_sex(*px["father"], sex_cutoff, sex_min_sites),
        "mom_sex": infer_sex(*px["mother"], sex_cutoff, sex_min_sites),
        "dad_x_sites": sum(px["father"]), "dad_x_het_ratio": _ratio(*px["father"]),
        "mom_x_sites": sum(px["mother"]), "mom_x_het_ratio": _ratio(*px["mother"]),
        "n_records": n_records,
        "nocall_rate": {r: ((nocall[r] / n_records) if n_records else None) for r in ("kid", "dad", "mom")},
        "charr": {r: C.charr(cref[r], cdp[r]) for r in ("kid", "dad", "mom")},
        "auto_dp_median": auto_med,
        "y": ys,
    }


def y_evidence(res, min_anchor, male_min, female_max):
    """Fold qc_trio's chrY scan into the report columns (blank = no evidence, never a value).

    Two bases for the proband's call: `anchor` — the father-normalised depth at his hemizygous
    sites (the cleanest, used whenever >= min_anchor anchors exist) — and `raw` — the haploid-
    scaled depth over every chrY record, the fallback when the father yields no usable anchors
    (a female or Y-less sample in the father slot, or an exome-sized callset). The same bands
    apply to both, since both read ~1 for a male and ~0 for a female; `y_basis` says which. The
    parents' raw ratios check their roles regardless (`dad_y_sex` / `mom_y_sex`).
    """
    ys, am = res.get("y"), res.get("auto_dp_median") or {}
    out = {"y_anchor_sites": "", "y_records": "", "kid_y_dp_median": "", "dad_y_dp_median": "",
           "mom_y_dp_median": "", "y_cov_ratio": "", "y_covered_frac": "", "mom_y_cov_ratio": "",
           "kid_y_cov_haploid": "", "dad_y_cov_haploid": "", "mom_y_cov_haploid": "",
           "dad_y_sex": "", "mom_y_sex": "", "parent_sex_flag_y": "",
           "y_inferred_sex": "", "y_basis": "", "y_flag": "no_chry_records",
           "_kid_ratio": None, "_mom_ratio": None, "_kid_h": None, "_dad_h": None, "_mom_h": None}
    if ys is None:
        out["y_flag"] = "no_chry_contig"
        return out
    n, nrec = ys["anchors"], ys["records"]
    out["y_anchor_sites"], out["y_records"] = n, nrec
    if nrec == 0:
        return out
    # raw haploid-scaled coverage for all three members, and the parents' Y-derived sex
    h = {r: y_cov_haploid(ys.get(f"{r}_all_med"), am.get(r)) for r in ("kid", "dad", "mom")}
    dad_sex = sex_from_y_cov(h["dad"], male_min, female_max) if nrec >= min_anchor else None
    mom_sex = sex_from_y_cov(h["mom"], male_min, female_max) if nrec >= min_anchor else None
    out.update({
        "kid_y_cov_haploid": ("" if h["kid"] is None else f"{h['kid']:.3g}"),
        "dad_y_cov_haploid": ("" if h["dad"] is None else f"{h['dad']:.3g}"),
        "mom_y_cov_haploid": ("" if h["mom"] is None else f"{h['mom']:.3g}"),
        "dad_y_sex": dad_sex or "", "mom_y_sex": mom_sex or "",
        # only a POSITIVE reading against the role counts (too few records = no expectation)
        "parent_sex_flag_y": ("1" if (dad_sex == "2" or mom_sex == "1") else
                              "0" if (dad_sex or mom_sex) else ""),
        "_kid_h": h["kid"], "_dad_h": h["dad"], "_mom_h": h["mom"],
    })
    # the proband: father-anchored when the anchors suffice, raw otherwise
    if n >= min_anchor:
        kid_r = y_cov_ratio(ys["kid_med"], am.get("kid"), ys["dad_med"], am.get("dad"))
        mom_r = y_cov_ratio(ys["mom_med"], am.get("mom"), ys["dad_med"], am.get("dad"))
        basis = "anchor"
    elif nrec >= min_anchor and h["kid"] is not None:
        kid_r, mom_r, basis = h["kid"], h["mom"], "raw"
    else:
        kid_r, mom_r, basis = None, None, ""
    sex, flag = infer_sex_y(kid_r, mom_r, (nrec if basis == "raw" else n), min_anchor,
                            male_min, female_max)
    if basis == "" and n < min_anchor:
        flag = "father_y_anchors_low"
    if n:
        out.update({"kid_y_dp_median": f"{ys['kid_med']:.4g}",
                    "dad_y_dp_median": f"{ys['dad_med']:.4g}",
                    "mom_y_dp_median": f"{ys['mom_med']:.4g}",
                    "y_covered_frac": f"{ys['kid_cov'] / n:.3g}"})
        if basis == "anchor":
            out.update({"y_cov_ratio": ("" if kid_r is None else f"{kid_r:.3g}"),
                        "mom_y_cov_ratio": ("" if mom_r is None else f"{mom_r:.3g}")})
    out.update({"y_inferred_sex": sex or "", "y_basis": (basis if sex else ""), "y_flag": flag,
                "_kid_ratio": kid_r, "_mom_ratio": mom_r})
    return out


def parent_sex_diagnosis(dad_sex, mom_sex, cutoff) -> str:
    """What a parent-sex mismatch pattern can and cannot mean (the WARN text).

    A transposed father/mother produces EXACTLY father=2 AND mother=1 (each parent reads as the
    other's sex). A father reading female while the mother does NOT read male (2/2, 2/None) is
    not that pattern — it is a sample swap/mislabel, OR the far more common case on a
    diploid-called callset: `qc.x_het_male_max` set below where genuine males' chrX het ratios
    sit, so a known male reads female. A mother reading male while the father does not read
    female (1/1, None/1) is a swap or a genuinely het-poor chrX (LOH, low coverage). Naming
    "transposed" for a pattern a transposition cannot produce sent reviewers to check the wrong
    thing.
    """
    if dad_sex == "2" and mom_sex == "1":
        return ("both parents read as the OTHER sex: the trios-file father/mother roles for this "
                "trio are probably transposed, and every parent-of-origin call would be inverted")
    if dad_sex == "2":
        return ("the father's chrX reads female but the mother's does NOT read male, so this is "
                "not a transposition (that gives 2/1): either a sample swap/mislabel, or — the "
                f"common case on a diploid-called callset — qc.x_het_male_max ({cutoff}) is below "
                "where genuine males' chrX het ratios sit (compare dad_x_het_ratio with "
                "mom_x_het_ratio across the cohort; see docs/inheritance_and_genotype_qc.md)")
    return ("the mother's chrX reads male but the father's does NOT read female, so this is not "
            "a transposition (that gives 2/1): a sample swap/mislabel, or a mother whose chrX is "
            "genuinely het-poor (loss of heterozygosity, low coverage) — check mom_x_het_ratio "
            "and mom_x_sites")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-sites", type=int, default=None,
                    help="cap on QC-passing sites scanned per trio (MIE/CHARR and the chrX sex "
                         "scan); default qc.max_sites from the config (200000)")
    ap.add_argument("--mie-threshold", type=float, default=None,
                    help="override qc.mie_max from config")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    thr = G.GtThresholds.from_config(cfg, get)
    selfsm = C.read_selfsm(get(cfg, "resources.selfsm_dir", ""))   # verifyBamID FREEMIX (optional)
    freemix_thr = float(get(cfg, "qc.freemix_threshold", 0.05))
    charr_thr = float(get(cfg, "qc.charr_threshold", 0.02))
    # thresholds are config defaults (canonical-defaults table), not hardcoded law
    mie_thr = args.mie_threshold if args.mie_threshold is not None else float(get(cfg, "qc.mie_max", 0.02))
    sex_cutoff = float(get(cfg, "qc.x_het_male_max", 0.10))
    sex_min = int(get(cfg, "qc.sex_min_sites", 20))
    max_sites = args.max_sites if args.max_sites is not None else int(get(cfg, "qc.max_sites", 200000))
    max_nocall = float(get(cfg, "qc.max_nocall_rate", 0.10))
    # chrY coverage sex evidence (audit-only): the read floor is shared with Step 5's
    # y_female_reads flag; the ratio band is where a son (~1) and a daughter (~0) separate.
    y_reads_min = int(get(cfg, "qc.y_reads_min_dp", 3))
    y_min_anchor = int(get(cfg, "qc.y_min_anchor_sites", 50))
    y_male_min = float(get(cfg, "qc.y_cov_male_min", 0.30))
    y_female_max = float(get(cfg, "qc.y_cov_female_max", 0.10))

    with open(args.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    # sex_source: which sex Step 5 will judge the proband's chrX ploidy under — `ped` (the trios
    # file stated it; canonical), `inferred` (it did not; Step 0's chrX call fills in) or `none`
    # (neither; Step 5 skips the sex chromosomes). sex_match is THREE-state: 1 = the PED sex and a
    # positive inference agree, 0 = they disagree (reported, never resolved here), blank = no
    # comparison was possible (no stated sex, or too few chrX sites). The parents' raw het ratios
    # (dad_/mom_x_het_ratio) are the calibration set for qc.x_het_male_max.
    cols = ["trio_id", "n_sites", "mie_errors", "mie_rate", "x_sites", "x_het_ratio", "inferred_sex",
            "ped_sex", "sex_source", "sex_match", "dad_inferred_sex", "dad_x_sites", "dad_x_het_ratio",
            "mom_inferred_sex", "mom_x_sites", "mom_x_het_ratio", "parent_sex_flag",
            "mie_flag", "kid_contam", "dad_contam", "mom_contam", "contam_source", "contam_flag",
            "n_records", "kid_nocall_rate", "dad_nocall_rate", "mom_nocall_rate", "nocall_flag",
            # chrY coverage evidence (audit-only; see the module docstring): the father-anchored
            # depth ratios, the resulting call, and its two comparisons — against the PED sex
            # (sex_match_y, three-state like sex_match) and against the chrX inference (xy_agree)
            "y_anchor_sites", "y_records", "kid_y_dp_median", "dad_y_dp_median", "mom_y_dp_median",
            "y_cov_ratio", "y_covered_frac", "mom_y_cov_ratio", "y_inferred_sex", "y_basis",
            "sex_match_y", "xy_agree", "y_flag",
            # ...and the RAW haploid-scaled chrY coverage of all three members, from which BOTH
            # PARENTS' sex is checked against their roles (a transposed pair reads dad=2/mom=1)
            "kid_y_cov_haploid", "dad_y_cov_haploid", "mom_y_cov_haploid", "dad_y_sex", "mom_y_sex",
            "parent_sex_flag_y",
            "overall_pass"]
    n_fail = n_done = 0
    n_flag = {"parent_sex": 0, "nocall": 0, "sex_discordant": 0, "sex_y_discordant": 0,
              "xy_disagree": 0, "mother_y_coverage": 0, "parent_sex_y": 0}
    n_source = {"ped": 0, "inferred": 0, "none": 0}
    n_y_sex = {"1": 0, "2": 0, "": 0}
    n_y_basis = {"anchor": 0, "raw": 0}
    # the chrY calibration evidence: proband ratios by PED sex, and the mothers' (known females)
    y_ratio_by_ped = {"1": [], "2": []}
    mother_y_ratios, father_y_anchors = [], []
    # ...and the parents' raw haploid coverage: fathers are known males, mothers known females
    n_father_y = {"1": 0, "2": 0}
    n_mother_y = {"1": 0, "2": 0}
    father_y_hap, mother_y_hap = [], []
    # calibration evidence for qc.x_het_male_max: fathers are known males, mothers known females
    n_father = {"1": 0, "2": 0}
    n_mother = {"1": 0, "2": 0}
    father_ratios, mother_ratios = [], []
    # fathers reading female OUTSIDE the transposition signature (father=2 AND mother=1): a
    # transposed trio's "father" is a real female and reads so under any cutoff, so it is
    # evidence about the roles, not about the cutoff — the calibration guard must not count it
    n_father_female_uncalibrated = 0

    def member_contam(sample, role, charr_map):
        """Return (value, source, flagged) for one member — verifyBamID freemix if present,
        else the VCF-only CHARR proxy."""
        if sample in selfsm:
            return selfsm[sample], "freemix", selfsm[sample] > freemix_thr
        cv = charr_map.get(role)
        return cv, "charr", (cv is not None and cv > charr_thr)

    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            tid, vcf_path, ped_path = r.get("trio_id"), r.get("vcf"), r.get("ped")
            ped = parse_ped(ped_path)
            if not ped:
                sys.stderr.write(f"WARN: no PED for {tid}; skipping QC\n")
                continue
            res = qc_trio(vcf_path, ped, thr, max_sites, sex_cutoff, sex_min, y_reads_min)
            if res is None:
                sys.stderr.write(f"WARN: {tid}: PED samples not in VCF; skipping\n")
                continue
            ped_sex = str(ped["sex"])
            # Sex-swap gate, through the ONE precedence rule (ped.resolve_child_sex) so this column
            # says exactly what Step 5 will do: a stated PED sex is canonical and only a positive
            # inference that DISAGREES is a failure (0); agreement is 1; no stated sex or no
            # inference is blank — "no comparison", not "matched".
            _sex, sex_source, discordant = resolve_child_sex(ped_sex, res["inferred_sex"])
            sex_match = ("0" if discordant else
                         "1" if (sex_source == "ped" and res["inferred_sex"] is not None) else "")
            n_source[sex_source] += 1
            n_flag["sex_discordant"] += discordant
            mie_flag = "1" if (res["mie_rate"] is not None and res["mie_rate"] > mie_thr) else "0"

            contam, flags, src = {}, False, "charr"
            for role, key in (("kid", "child"), ("dad", "father"), ("mom", "mother")):
                val, s, fl = member_contam(ped[key], role, res["charr"])
                contam[role] = val
                flags = flags or fl
                if s == "freemix":
                    src = "freemix"
            contam_flag = "1" if flags else "0"
            # "charr" with every value missing means nobody was measured, not that CHARR ran
            if src == "charr" and all(val is None for val in contam.values()):
                src = "none"
            # a father whose own chrX reads female, or a mother whose reads male: the PED roles
            # are wrong for this trio (or a sample is), and every parent-of-origin call would be
            # inverted. Only a positive inference counts — too few X sites is "no expectation".
            parent_sex_flag = "1" if (res["dad_sex"] == "2" or res["mom_sex"] == "1") else "0"
            for sx, tally in ((res["dad_sex"], n_father), (res["mom_sex"], n_mother)):
                if sx in tally:
                    tally[sx] += 1
            if res["dad_sex"] == "2" and res["mom_sex"] != "1":
                n_father_female_uncalibrated += 1
            if res["dad_x_het_ratio"] is not None:
                father_ratios.append(res["dad_x_het_ratio"])
            if res["mom_x_het_ratio"] is not None:
                mother_ratios.append(res["mom_x_het_ratio"])
            nc = res["nocall_rate"]
            nocall_flag = "1" if any(nc[r] is not None and nc[r] > max_nocall for r in nc) else "0"
            n_flag["parent_sex"] += parent_sex_flag == "1"
            n_flag["nocall"] += nocall_flag == "1"
            # chrY coverage evidence: two three-state comparisons, exactly like sex_match —
            # 1 agree / 0 disagree / blank = not comparable. Neither changes overall_pass.
            ye = y_evidence(res, y_min_anchor, y_male_min, y_female_max)
            ysex = ye["y_inferred_sex"]
            sex_match_y = ("" if not ysex or ped_sex not in ("1", "2") else
                           "1" if ysex == ped_sex else "0")
            xy_agree = ("" if not ysex or res["inferred_sex"] not in ("1", "2") else
                        "1" if ysex == res["inferred_sex"] else "0")
            n_y_sex[ysex] += 1
            if ye["y_basis"]:
                n_y_basis[ye["y_basis"]] += 1
            n_flag["sex_y_discordant"] += sex_match_y == "0"
            n_flag["xy_disagree"] += xy_agree == "0"
            n_flag["mother_y_coverage"] += ye["y_flag"] == "mother_y_coverage"
            n_flag["parent_sex_y"] += ye["parent_sex_flag_y"] == "1"
            for sx, tally in ((ye["dad_y_sex"], n_father_y), (ye["mom_y_sex"], n_mother_y)):
                if sx in tally:
                    tally[sx] += 1
            if ye["_dad_h"] is not None:
                father_y_hap.append(ye["_dad_h"])
            if ye["_mom_h"] is not None:
                mother_y_hap.append(ye["_mom_h"])
            if ye["_kid_ratio"] is not None and ped_sex in y_ratio_by_ped:
                y_ratio_by_ped[ped_sex].append(ye["_kid_ratio"])
            if ye["_mom_ratio"] is not None:
                mother_y_ratios.append(ye["_mom_ratio"])
            if ye["y_anchor_sites"] != "":
                father_y_anchors.append(int(ye["y_anchor_sites"]))

            # ADVISORY only: overall_pass folds the flags into one column that is SURFACED to the
            # mandatory human reviewer (xlsx QC sheet + IGV sample_qc.tsv) but does NOT auto-
            # exclude a trio — Steps 1/2/4 run over the full resolved manifest and Step 6 never reads
            # it, so a flagged trio still contributes calls and recurrence pending human review.
            # Automated gating is deliberately deferred (never-drop ethos + the contamination proxy's
            # limited sensitivity, see contamination.py); it would be a config-gated policy change.
            overall = "1" if (mie_flag == "0" and sex_match != "0" and contam_flag == "0"
                              and parent_sex_flag == "0" and nocall_flag == "0") else "0"
            if overall == "0":
                n_fail += 1
            n_done += 1
            row = {
                "trio_id": tid, "n_sites": res["n_sites"], "mie_errors": res["mie_errors"],
                "mie_rate": ("" if res["mie_rate"] is None else f"{res['mie_rate']:.4g}"),
                "x_sites": res["x_sites"],
                "x_het_ratio": ("" if res["x_het_ratio"] is None else f"{res['x_het_ratio']:.3g}"),
                "inferred_sex": res["inferred_sex"] or "", "ped_sex": ped_sex,
                "sex_source": sex_source, "sex_match": sex_match,
                "dad_inferred_sex": res["dad_sex"] or "", "dad_x_sites": res["dad_x_sites"],
                "dad_x_het_ratio": ("" if res["dad_x_het_ratio"] is None else f"{res['dad_x_het_ratio']:.3g}"),
                "mom_inferred_sex": res["mom_sex"] or "", "mom_x_sites": res["mom_x_sites"],
                "mom_x_het_ratio": ("" if res["mom_x_het_ratio"] is None else f"{res['mom_x_het_ratio']:.3g}"),
                "parent_sex_flag": parent_sex_flag,
                "mie_flag": mie_flag,
                "kid_contam": ("" if contam["kid"] is None else f"{contam['kid']:.4g}"),
                "dad_contam": ("" if contam["dad"] is None else f"{contam['dad']:.4g}"),
                "mom_contam": ("" if contam["mom"] is None else f"{contam['mom']:.4g}"),
                "contam_source": src, "contam_flag": contam_flag,
                "n_records": res["n_records"],
                "kid_nocall_rate": ("" if nc["kid"] is None else f"{nc['kid']:.4g}"),
                "dad_nocall_rate": ("" if nc["dad"] is None else f"{nc['dad']:.4g}"),
                "mom_nocall_rate": ("" if nc["mom"] is None else f"{nc['mom']:.4g}"),
                "nocall_flag": nocall_flag, "overall_pass": overall,
                "y_anchor_sites": ye["y_anchor_sites"], "kid_y_dp_median": ye["kid_y_dp_median"],
                "dad_y_dp_median": ye["dad_y_dp_median"], "mom_y_dp_median": ye["mom_y_dp_median"],
                "y_cov_ratio": ye["y_cov_ratio"], "y_covered_frac": ye["y_covered_frac"],
                "mom_y_cov_ratio": ye["mom_y_cov_ratio"], "y_inferred_sex": ysex,
                "y_basis": ye["y_basis"], "y_records": ye["y_records"],
                "sex_match_y": sex_match_y, "xy_agree": xy_agree, "y_flag": ye["y_flag"],
                "kid_y_cov_haploid": ye["kid_y_cov_haploid"], "dad_y_cov_haploid": ye["dad_y_cov_haploid"],
                "mom_y_cov_haploid": ye["mom_y_cov_haploid"], "dad_y_sex": ye["dad_y_sex"],
                "mom_y_sex": ye["mom_y_sex"], "parent_sex_flag_y": ye["parent_sex_flag_y"],
            }
            out.write("\t".join(str(row[c]) for c in cols) + "\n")
            # Step 0 recorded nothing in the audit before; the flags are the useful part
            for metric in ("overall_pass", "mie_flag", "sex_match", "parent_sex_flag",
                           "contam_flag", "nocall_flag", "sex_match_y", "xy_agree",
                           "parent_sex_flag_y"):
                audit.record("00_qc", metric, row[metric], scope=tid)
            audit.record("00_qc", f"sex_source.{sex_source}", 1, scope=tid)
            audit.record("00_qc", f"y_inferred_sex.{ysex or 'none'}", 1, scope=tid)
            if sex_match_y == "0":
                sys.stderr.write(
                    f"WARN: {tid}: PROBAND SEX vs chrY COVERAGE MISMATCH — the trios file says "
                    f"{ped_sex}, but the proband's depth at the father's {ye['y_anchor_sites']} "
                    f"hemizygous chrY sites reads {ysex} (1=male, 2=female; y_cov_ratio "
                    f"{ye['y_cov_ratio']} of the father's, covered fraction {ye['y_covered_frac']}; "
                    f"the mother reads {ye['mom_y_cov_ratio'] or '?'}). The PED stays canonical and "
                    "overall_pass is unchanged (this check is audit-only), but chrY coverage is "
                    "near-binary on a diploid-called callset: a genuine mismatch is a sample-swap or "
                    "mislabel tell for the whole trio"
                    + (" — and the chrX inference agrees with chrY here, so the pedigree is the odd "
                       "one out" if xy_agree == "1" else "") + ".\n")
            if xy_agree == "0":
                sys.stderr.write(
                    f"WARN: {tid}: the chrX het-ratio inference ({res['inferred_sex']}) and chrY "
                    f"coverage ({ysex}, y_cov_ratio {ye['y_cov_ratio']} over {ye['y_anchor_sites']} "
                    "father-hemizygous sites) DISAGREE on the proband's sex. On a diploid-called "
                    f"callset the chrX cutoff (qc.x_het_male_max={sex_cutoff}) is the usual culprit "
                    "and chrY coverage the more robust signal; a stated trios-file sex is unaffected, "
                    "but a proband WITHOUT one has its chrX ploidy judged on the chrX inference in "
                    "Step 5 — check this trio.\n")
            if ye["parent_sex_flag_y"] == "1":
                sys.stderr.write(
                    f"WARN: {tid}: PARENT SEX MISMATCH ON chrY COVERAGE — the father slot reads "
                    f"{ye['dad_y_sex'] or '?'} (dad_y_cov_haploid {ye['dad_y_cov_haploid'] or '?'}), "
                    f"the mother slot reads {ye['mom_y_sex'] or '?'} (mom_y_cov_haploid "
                    f"{ye['mom_y_cov_haploid'] or '?'}; 1=male, 2=female, a male ~1, a female ~0 over "
                    f"{ye['y_records']} chrY records): "
                    f"{parent_y_diagnosis(ye['dad_y_sex'], ye['mom_y_sex'])}. Unlike the chrX "
                    "parent check this does not depend on a cutoff; overall_pass is unchanged "
                    "(audit-only).\n")
            if ye["y_flag"] == "mother_y_coverage":
                sys.stderr.write(
                    f"WARN: {tid}: the MOTHER shows chrY coverage (mom_y_cov_ratio "
                    f"{ye['mom_y_cov_ratio']} of the father's at his hemizygous sites, above "
                    f"qc.y_cov_female_max={y_female_max}) — the in-trio female control failed "
                    "(contamination, transposed parents, or a sample swap), so NO chrY sex call is "
                    "made for this proband rather than a wrong one.\n")
            if discordant:
                sys.stderr.write(
                    f"WARN: {tid}: PROBAND SEX MISMATCH — the trios file says {ped_sex}, the chrX "
                    f"inference reads {res['inferred_sex']} (1=male, 2=female; het ratio "
                    f"{row['x_het_ratio']} over {res['x_sites']} sites, cutoff {sex_cutoff}). The "
                    "PED sex is canonical: Step 5 keeps it and flags every call for this trio "
                    "sex_discordant_inference. If the pedigree is wrong, fix the trios file; if the "
                    "inference is wrong for this callset, calibrate qc.x_het_male_max against the "
                    "parents' dad_x_het_ratio / mom_x_het_ratio (docs/inheritance_and_genotype_qc.md). "
                    "A genuine mismatch is a sample-swap tell for the whole trio.\n")
            if parent_sex_flag == "1":
                sys.stderr.write(f"WARN: {tid}: PARENT SEX MISMATCH — father's chrX reads "
                                 f"{res['dad_sex'] or '?'}, mother's reads {res['mom_sex'] or '?'} "
                                 "(1=male, 2=female): "
                                 f"{parent_sex_diagnosis(res['dad_sex'], res['mom_sex'], sex_cutoff)}.\n")
            if nocall_flag == "1":
                sys.stderr.write(f"WARN: {tid}: no-call rate above qc.max_nocall_rate ({max_nocall}) "
                                 f"— kid {row['kid_nocall_rate']} dad {row['dad_nocall_rate']} mom "
                                 f"{row['mom_nocall_rate']}. Was this trio jointly genotyped? Step 5 "
                                 "treats a parental ./. as uninformative, never as 0/0.\n")

    audit.record("00_qc", "trios_qc", n_done)
    audit.record("00_qc", "trios_flagged", n_fail)
    audit.record("00_qc", "trios_parent_sex_flag", n_flag["parent_sex"])
    audit.record("00_qc", "trios_sex_discordant_inference", n_flag["sex_discordant"])
    audit.record("00_qc", "trios_nocall_flag", n_flag["nocall"])
    for src, n in sorted(n_source.items()):
        audit.record("00_qc", f"trios_sex_source.{src}", n)
    # chrY coverage evidence, cohort-wide: how many probands it called each way, how often it
    # contradicted the pedigree or the chrX inference, and the calibration medians (proband
    # ratio by PED sex; the mothers, all known females; the fathers' anchor counts).
    audit.record("00_qc", "y_sex_inferred_male", n_y_sex["1"])
    audit.record("00_qc", "y_sex_inferred_female", n_y_sex["2"])
    audit.record("00_qc", "y_sex_inferred_none", n_y_sex[""])
    audit.record("00_qc", "trios_sex_match_y_discordant", n_flag["sex_y_discordant"])
    audit.record("00_qc", "trios_xy_disagree", n_flag["xy_disagree"])
    audit.record("00_qc", "trios_mother_y_coverage", n_flag["mother_y_coverage"])
    audit.record("00_qc", "trios_parent_sex_flag_y", n_flag["parent_sex_y"])
    audit.record("00_qc", "y_sex_basis_anchor", n_y_basis["anchor"])
    audit.record("00_qc", "y_sex_basis_raw", n_y_basis["raw"])
    audit.record("00_qc", "fathers_y_male", n_father_y["1"])
    audit.record("00_qc", "fathers_y_female", n_father_y["2"])
    audit.record("00_qc", "mothers_y_female", n_mother_y["2"])
    audit.record("00_qc", "mothers_y_male", n_mother_y["1"])
    if father_y_hap:
        audit.record("00_qc", "father_y_cov_haploid_median", f"{statistics.median(father_y_hap):.3g}")
    if mother_y_hap:
        audit.record("00_qc", "mother_y_cov_haploid_median", f"{statistics.median(mother_y_hap):.3g}")
    for sx, lbl in (("1", "ped_male"), ("2", "ped_female")):
        if y_ratio_by_ped[sx]:
            audit.record("00_qc", f"proband_y_cov_ratio_median.{lbl}",
                         f"{statistics.median(y_ratio_by_ped[sx]):.3g}")
    if mother_y_ratios:
        audit.record("00_qc", "mother_y_cov_ratio_median", f"{statistics.median(mother_y_ratios):.3g}")
    if father_y_anchors:
        audit.record("00_qc", "father_y_anchor_sites_median", f"{statistics.median(father_y_anchors):.4g}")
    # The calibration evidence, in the audit so a methods section can quote it: how the KNOWN
    # males (fathers) and KNOWN females (mothers) read under the configured cutoff.
    audit.record("00_qc", "fathers_inferred_male", n_father["1"])
    audit.record("00_qc", "fathers_inferred_female", n_father["2"])
    audit.record("00_qc", "mothers_inferred_female", n_mother["2"])
    audit.record("00_qc", "mothers_inferred_male", n_mother["1"])
    if father_ratios:
        audit.record("00_qc", "father_x_het_ratio_median", f"{statistics.median(father_ratios):.3g}")
    if mother_ratios:
        audit.record("00_qc", "mother_x_het_ratio_median", f"{statistics.median(mother_ratios):.3g}")
    sys.stderr.write(f"Step 0 complete: QC report -> {args.out} ({n_fail} trio(s) flagged; "
                     f"parent-sex {n_flag['parent_sex']}, proband sex discordant "
                     f"{n_flag['sex_discordant']}, no-call {n_flag['nocall']}); proband sex from "
                     f"the trios file for {n_source['ped']}, from the chrX inference for "
                     f"{n_source['inferred']}, unresolved for {n_source['none']}\n"
                     f"  chrY coverage (audit-only): probands reading male {n_y_sex['1']}, female "
                     f"{n_y_sex['2']}, no call {n_y_sex['']}; contradicting the pedigree "
                     f"{n_flag['sex_y_discordant']}, contradicting the chrX inference "
                     f"{n_flag['xy_disagree']}, mothers with chrY coverage "
                     f"{n_flag['mother_y_coverage']}; parents' roles contradicted by chrY "
                     f"coverage {n_flag['parent_sex_y']} (fathers reading male {n_father_y['1']} / "
                     f"female {n_father_y['2']}, mothers reading female {n_mother_y['2']} / male "
                     f"{n_mother_y['1']})"
                     + (f"; median proband y_cov_ratio: PED-male "
                        f"{statistics.median(y_ratio_by_ped['1']):.3g}" if y_ratio_by_ped["1"] else "")
                     + (f", PED-female {statistics.median(y_ratio_by_ped['2']):.3g}"
                        if y_ratio_by_ped["2"] else "")
                     + (f"; mothers {statistics.median(mother_y_ratios):.3g}" if mother_y_ratios else "")
                     + "\n")
    if n_flag["xy_disagree"] and n_flag["xy_disagree"] >= max(1, n_y_sex["1"] + n_y_sex["2"]) // 2:
        sys.stderr.write(
            f"WARN: chrX and chrY sex evidence disagree on {n_flag['xy_disagree']} of "
            f"{n_y_sex['1'] + n_y_sex['2']} probands with a chrY call. Systematic disagreement is "
            f"the signature of a miscalibrated qc.x_het_male_max ({sex_cutoff}) on a diploid-called "
            "callset — chrY coverage is near-binary and is the signal to trust when choosing the "
            "chrX cutoff (proband_y_cov_ratio_median.* and the parents' chrX ratios are in the "
            "audit).\n")
    # CALIBRATION GUARD. Fathers are known males, so more fathers reading female than male is
    # not mass transposition — it is the cutoff sitting below where this callset's males lie
    # (a GATK diploid-called male's non-PAR chrX het ratio routinely exceeds 0.10). Parameter-
    # free on purpose: a majority of known males mis-read needs no threshold to be wrong. Trios
    # showing the transposition signature (father=2 AND mother=1) are excluded from the count.
    if n_father_female_uncalibrated > n_father["1"]:
        sys.stderr.write(
            f"WARN: qc.x_het_male_max ({sex_cutoff}) looks MISCALIBRATED for this callset: "
            f"{n_father_female_uncalibrated} of {n_father['1'] + n_father_female_uncalibrated} "
            f"fathers with a chrX inference (transposition-pattern trios excluded) read "
            f"FEMALE (median father het ratio "
            f"{statistics.median(father_ratios) if father_ratios else float('nan'):.3g}, median "
            f"mother {statistics.median(mother_ratios) if mother_ratios else float('nan'):.3g}). "
            "Set the cutoff in the gap between the father and mother distributions "
            "(dad_x_het_ratio / mom_x_het_ratio in the report; docs/inheritance_and_genotype_qc.md) "
            "and rm qc_report.tsv.done. Until then: the parent_sex_flag rows above are NOT evidence "
            f"of transposition, a stated trios-file sex is unaffected ({n_source['ped']} trio(s)), "
            f"but {n_source['inferred']} trio(s) with NO stated sex have their chrX ploidy judged "
            "on this inference in Step 5.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
