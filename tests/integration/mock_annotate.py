#!/usr/bin/env python3
"""Stand in for Step 2's VEP pass: inject vep_* INFO fields from a lookup table.

VEP + its cache are too heavy for CI, so for the integration test we attach the
vep_* fields (the contract src/hprv/annotations.py reads) directly from annot.tsv.
The external gnomAD/ClinVar transfer is done separately with real `bcftools annotate`
in run_integration.sh, so only the VEP invocation itself is mocked.
"""
from __future__ import annotations

import argparse
import csv

from cyvcf2 import VCF, Writer

STR_FIELDS = [("vep_Consequence", "csq"), ("vep_IMPACT", "impact"),
              ("vep_SYMBOL", "gene"), ("vep_Gene", "gene"), ("vep_LoF", "loftee")]
FLOAT_FIELDS = [("vep_REVEL_score", "revel"), ("vep_AlphaMissense_score", "am"),
                ("vep_SpliceAI_pred_DS_AG", "spliceai")]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--lookup", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    lut = {}
    with open(args.lookup) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lut[(r["chrom"], int(r["pos"]), r["ref"], r["alt"])] = r

    vcf = VCF(args.inp)
    for fid, _ in STR_FIELDS:
        vcf.add_info_to_header({"ID": fid, "Number": "1", "Type": "String", "Description": "mock"})
    for fid, _ in FLOAT_FIELDS:
        vcf.add_info_to_header({"ID": fid, "Number": "1", "Type": "Float", "Description": "mock"})

    w = Writer(args.out, vcf)
    for v in vcf:
        row = lut.get((v.CHROM, v.POS, v.REF, v.ALT[0] if v.ALT else ""))
        if row:
            for fid, col in STR_FIELDS:
                if row.get(col):
                    v.INFO[fid] = row[col]
            for fid, col in FLOAT_FIELDS:
                if row.get(col):
                    v.INFO[fid] = float(row[col])
        w.write_record(v)
    w.close()
    vcf.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
