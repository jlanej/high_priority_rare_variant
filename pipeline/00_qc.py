#!/usr/bin/env python3
"""Pipeline Step 0: per-trio QC gate (garbage-in guard).

Self-contained (no extra resources): for each trio computes the Mendelian-error rate
(a sensitive proxy for sample swaps / mislabeled parents) and a chrX-heterozygosity
sex inference, and flags trios that fail. Somalier is the recommended richer check
(sex + relatedness + ancestry) when a somalier sites file is available; this step is
the dependency-light baseline. See docs/inheritance_and_genotype_qc.md.

Usage:
  00_qc.py --manifest trios.tsv --config config.yaml --out qc_report.tsv [--max-sites N]
"""
from __future__ import annotations

import argparse
import csv
import sys

from cyvcf2 import VCF

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


def qc_trio(vcf_path, ped, thr, max_sites):
    vcf = VCF(vcf_path)
    samples = {s: i for i, s in enumerate(vcf.samples)}
    for role in ("child", "father", "mother"):
        if ped[role] not in samples:
            vcf.close()
            return None
    c, d, m = samples[ped["child"]], samples[ped["father"]], samples[ped["mother"]]

    considered = errors = 0
    x_het = x_hom = 0
    for v in vcf:
        chrom = v.CHROM.replace("chr", "")
        # sex inference from chrX non-PAR calls in the child
        if chrom == "X" and not G.in_par_x(v) and len(v.ALT) == 1:
            gt = v.gt_types[c]
            gq = G.gq(v, c)
            if gq is not None and gq >= thr.min_gq:
                if gt == HET:
                    x_het += 1
                elif gt == HOM_ALT:
                    x_hom += 1
            continue
        if chrom in ("Y", "MT", "M"):
            continue
        if len(v.ALT) != 1:  # biallelic only for the simple MIE rule
            continue
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
    inferred = None
    if x_het_ratio is not None:
        inferred = "1" if x_het_ratio < 0.10 else "2"  # male vs female
    return {
        "n_sites": considered, "mie_errors": errors, "mie_rate": mie_rate,
        "x_het_ratio": x_het_ratio, "inferred_sex": inferred,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-sites", type=int, default=200000)
    ap.add_argument("--mie-threshold", type=float, default=0.02)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    thr = G.GtThresholds.from_config(cfg, get)

    with open(args.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    cols = ["trio_id", "n_sites", "mie_errors", "mie_rate", "x_het_ratio",
            "inferred_sex", "ped_sex", "sex_match", "mie_flag", "overall_pass"]
    n_fail = 0
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            tid, vcf_path, ped_path = r.get("trio_id"), r.get("vcf"), r.get("ped")
            ped = parse_ped(ped_path)
            if not ped:
                sys.stderr.write(f"WARN: no PED for {tid}; skipping QC\n")
                continue
            res = qc_trio(vcf_path, ped, thr, args.max_sites)
            if res is None:
                sys.stderr.write(f"WARN: {tid}: PED samples not in VCF; skipping\n")
                continue
            ped_sex = str(ped["sex"])
            sex_match = "1" if (res["inferred_sex"] == ped_sex or res["inferred_sex"] is None) else "0"
            mie_flag = "1" if (res["mie_rate"] is not None and res["mie_rate"] > args.mie_threshold) else "0"
            overall = "1" if (mie_flag == "0" and sex_match == "1") else "0"
            if overall == "0":
                n_fail += 1
            row = {
                "trio_id": tid, "n_sites": res["n_sites"], "mie_errors": res["mie_errors"],
                "mie_rate": ("" if res["mie_rate"] is None else f"{res['mie_rate']:.4g}"),
                "x_het_ratio": ("" if res["x_het_ratio"] is None else f"{res['x_het_ratio']:.3g}"),
                "inferred_sex": res["inferred_sex"] or "", "ped_sex": ped_sex,
                "sex_match": sex_match, "mie_flag": mie_flag, "overall_pass": overall,
            }
            out.write("\t".join(str(row[c]) for c in cols) + "\n")

    sys.stderr.write(f"Step 0 complete: QC report -> {args.out} ({n_fail} trio(s) flagged)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
