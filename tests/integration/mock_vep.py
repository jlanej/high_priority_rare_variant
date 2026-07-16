#!/usr/bin/env python3
"""Stand in for the `vep` binary: write a VEP-shaped CSQ INFO field from a lookup table.

VEP + its 24 GB cache are too heavy for CI, so this fakes the ONE thing the pipeline needs
from it — an `##INFO=<ID=CSQ...Format: ...>` header plus a per-record CSQ value — and lets
Step 2 do everything else for real (header/build verification, the split-vep lift, the
transcript selector, and the frequency guard). That is deliberately more than the previous
mock did: it wrote vep_* INFO fields directly, which meant Step 2 itself was never executed
and its checks were never tested.

Emitted to match `vep --everything --flag_pick --plugin CADD` output shape:
  * CSQ carries PICK, so Step 2's selector resolves to `-s pick` exactly as on real data.
  * frequency lives in the per-population gnomAD CSQ fields, not an external VCF.
  * CLIN_SIG is lowercase and '&'-joined, as VEP writes it (NOT ClinVar's Capitalised CLNSIG).

Usage: mock_vep.py --in sites.vcf.gz --lookup annot.tsv --out vep.vcf
"""
from __future__ import annotations

import argparse
import csv

from cyvcf2 import VCF, Writer

# The CSQ field order this mock emits. Superset of what Step 2 asks for, so the script also
# exercises Step 2's "lift only the fields that are present" intersection logic.
CSQ_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE",
    "HGVSc", "HGVSp", "CANONICAL", "MANE_SELECT", "PICK",
    "CADD_PHRED", "CADD_RAW", "CLIN_SIG",
    "gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF",
    "gnomADe_FIN_AF", "gnomADe_MID_AF", "gnomADe_NFE_AF", "gnomADe_SAS_AF",
    "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF", "gnomADg_ASJ_AF",
    "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF", "gnomADg_NFE_AF", "gnomADg_SAS_AF",
    "MAX_AF", "MAX_AF_POPS",
]

# VEP's own header lines. Step 2 checks these to refuse a wrong-build annotation, so the mock
# must carry them or the ingest path fails closed (which is the intended behaviour).
VEP_HEADER = (
    '##VEP="v115.0" API="v115" time="mock" cache="mock/homo_sapiens/115_GRCh38" '
    'assembly="GRCh38.p14" gnomADe="v4.1" gnomADg="v4.1" ClinVar="202502" dbSNP="156"'
)


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
    vcf.add_info_to_header({
        "ID": "CSQ", "Number": ".", "Type": "String",
        "Description": "Consequence annotations from Ensembl VEP. Format: " + "|".join(CSQ_FIELDS),
    })
    vcf.add_to_header(VEP_HEADER)

    w = Writer(args.out, vcf)
    n_csq = 0
    for v in vcf:
        row = lut.get((v.CHROM, v.POS, v.REF, v.ALT[0] if v.ALT else ""))
        if row:
            f = {k: "" for k in CSQ_FIELDS}
            f["Allele"] = v.ALT[0] if v.ALT else ""
            f["Consequence"] = row["csq"]
            f["IMPACT"] = row["impact"]
            f["SYMBOL"] = row["gene"]
            f["Gene"] = row["gene"]
            f["Feature_type"] = "Transcript"
            f["Feature"] = f"ENST_MOCK_{row['gene']}"
            f["BIOTYPE"] = "protein_coding"
            f["CANONICAL"] = "YES"
            f["MANE_SELECT"] = f"NM_MOCK_{row['gene']}"
            f["PICK"] = "1"
            if row.get("cadd"):
                f["CADD_PHRED"] = row["cadd"]
                f["CADD_RAW"] = row["cadd"]
            if row.get("clnsig"):
                f["CLIN_SIG"] = row["clnsig"]
            if row.get("af"):
                # Put the AF in the requested population only. A grpmax-eligible group drives
                # annotations.frequency(); a bottlenecked one (mid/ami/asj/fin) must not — and
                # MAX_AF sees it either way, which is exactly the trap frequency() avoids.
                pop = row.get("af_pop") or "gnomADe_NFE_AF"
                if pop not in f:
                    raise SystemExit(f"annot.tsv af_pop {pop!r} is not a CSQ field this mock emits")
                f[pop] = row["af"]
                f["MAX_AF"] = row["af"]
                f["MAX_AF_POPS"] = pop.replace("_AF", "")
            v.INFO["CSQ"] = "|".join(f[k] for k in CSQ_FIELDS)
            n_csq += 1
        w.write_record(v)
    w.close()
    vcf.close()
    print(f"mock_vep: wrote CSQ for {n_csq} records -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
