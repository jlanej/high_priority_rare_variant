#!/usr/bin/env python3
"""Pipeline Step 0: per-trio QC gate (garbage-in guard).

Mostly self-contained (no extra resources): for each trio computes three gates and
flags trios that fail any of them:
  * Mendelian-error rate — a sensitive proxy for sample swaps / mislabeled parents.
  * chrX-heterozygosity sex inference vs. the PED sex.
  * Contamination — verifyBamID FREEMIX if a directory of ``*.selfSM`` files is
    configured (resources.selfsm_dir), else a VCF-only CHARR estimate (reference-read
    fraction at high-quality hom-ALT SNV sites). Mirrors the group's DNM freemix QC.

Trio kid/dad/mom roles come from upstream peddy; this step is the guard against the
less-well-curated trios. See docs/inheritance_and_genotype_qc.md.

Usage:
  00_qc.py --manifest trios.tsv --config config.yaml --out qc_report.tsv [--max-sites N]
"""
from __future__ import annotations

import argparse
import csv
import sys

from cyvcf2 import VCF

from hprv import contamination as C
from hprv import genotype as G
from hprv.config import get, load_config
from hprv.ped import parse_ped

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


def _pick_x_contig(seqnames):
    for want in ("chrX", "X"):
        if want in seqnames:
            return want
    return None


def scan_sex(vcf_path, child_id, thr, max_x):
    """chrX non-PAR het/hom counts in the child, as a DEDICATED pass.

    chrX sorts after all autosomes, so an autosomal MIE cap in the main pass would
    otherwise starve sex inference on WGS trios. Uses the index to jump straight to
    chrX; falls back to a full scan (filtered to chrX) when the VCF is unindexed.
    """
    vcf = VCF(vcf_path)
    ci = {s: i for i, s in enumerate(vcf.samples)}.get(child_id)
    x = _pick_x_contig(vcf.seqnames)
    if ci is None or x is None:
        vcf.close()
        return 0, 0
    try:
        it = vcf(x)                       # indexed region jump (preferred)
    except Exception:
        it = vcf                          # unindexed: full scan, filtered below
    x_het = x_hom = 0
    for v in it:
        if v.CHROM.replace("chr", "") != "X" or G.in_par_x(v) or len(v.ALT) != 1:
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


def qc_trio(vcf_path, ped, thr, max_sites, sex_cutoff=0.10, sex_min_sites=20):
    vcf = VCF(vcf_path)
    samples = {s: i for i, s in enumerate(vcf.samples)}
    for role in ("child", "father", "mother"):
        if ped[role] not in samples:
            vcf.close()
            return None
    c, d, m = samples[ped["child"]], samples[ped["father"]], samples[ped["mother"]]

    # sex inference runs as its own chrX pass so the autosomal MIE cap can't disable it
    x_het, x_hom = scan_sex(vcf_path, ped["child"], thr, max_x=(max_sites or 0))

    considered = errors = 0
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
        if G.UNKNOWN in (gc, gd, gm):
            continue
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
    inferred = None
    if x_total >= sex_min_sites and x_het_ratio is not None:
        inferred = "1" if x_het_ratio < sex_cutoff else "2"  # male vs female
    return {
        "n_sites": considered, "mie_errors": errors, "mie_rate": mie_rate,
        "x_sites": x_total, "x_het_ratio": x_het_ratio, "inferred_sex": inferred,
        "charr": {r: C.charr(cref[r], cdp[r]) for r in ("kid", "dad", "mom")},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-sites", type=int, default=200000)
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

    with open(args.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    cols = ["trio_id", "n_sites", "mie_errors", "mie_rate", "x_sites", "x_het_ratio", "inferred_sex",
            "ped_sex", "sex_match", "mie_flag",
            "kid_contam", "dad_contam", "mom_contam", "contam_source", "contam_flag",
            "overall_pass"]
    n_fail = 0

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
            res = qc_trio(vcf_path, ped, thr, args.max_sites, sex_cutoff, sex_min)
            if res is None:
                sys.stderr.write(f"WARN: {tid}: PED samples not in VCF; skipping\n")
                continue
            ped_sex = str(ped["sex"])
            sex_match = "1" if (res["inferred_sex"] == ped_sex or res["inferred_sex"] is None) else "0"
            mie_flag = "1" if (res["mie_rate"] is not None and res["mie_rate"] > mie_thr) else "0"

            contam, flags, src = {}, False, "charr"
            for role, key in (("kid", "child"), ("dad", "father"), ("mom", "mother")):
                val, s, fl = member_contam(ped[key], role, res["charr"])
                contam[role] = val
                flags = flags or fl
                if s == "freemix":
                    src = "freemix"
            contam_flag = "1" if flags else "0"

            overall = "1" if (mie_flag == "0" and sex_match == "1" and contam_flag == "0") else "0"
            if overall == "0":
                n_fail += 1
            row = {
                "trio_id": tid, "n_sites": res["n_sites"], "mie_errors": res["mie_errors"],
                "mie_rate": ("" if res["mie_rate"] is None else f"{res['mie_rate']:.4g}"),
                "x_sites": res["x_sites"],
                "x_het_ratio": ("" if res["x_het_ratio"] is None else f"{res['x_het_ratio']:.3g}"),
                "inferred_sex": res["inferred_sex"] or "", "ped_sex": ped_sex,
                "sex_match": sex_match, "mie_flag": mie_flag,
                "kid_contam": ("" if contam["kid"] is None else f"{contam['kid']:.4g}"),
                "dad_contam": ("" if contam["dad"] is None else f"{contam['dad']:.4g}"),
                "mom_contam": ("" if contam["mom"] is None else f"{contam['mom']:.4g}"),
                "contam_source": src, "contam_flag": contam_flag, "overall_pass": overall,
            }
            out.write("\t".join(str(row[c]) for c in cols) + "\n")

    sys.stderr.write(f"Step 0 complete: QC report -> {args.out} ({n_fail} trio(s) flagged)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
