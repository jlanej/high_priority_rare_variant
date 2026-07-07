#!/usr/bin/env python3
"""Assert the mock run produced the expected resolution, funnel, and calls."""
from __future__ import annotations

import argparse
import csv
import os
import sys

from cyvcf2 import VCF

FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    a = ap.parse_args(argv)
    W = a.work

    # --- resolution ---
    res = {r["kid"]: r for r in rows(os.path.join(W, "trio_resolution.tsv"))}
    check(res.get("CH_A", {}).get("status", "").startswith("resolved"), "CH_A resolved")
    check(res.get("CH_B", {}).get("status", "").startswith("resolved"), "CH_B resolved (family VCF, extra sib)")
    check(res.get("CH_C", {}).get("status") == "unresolved", "CH_C unresolved")
    check("MO_C" in res.get("CH_C", {}).get("missing_members", ""), "CH_C reports missing mom MO_C")
    manifest = rows(os.path.join(W, "trios.resolved.tsv"))
    check(len(manifest) == 2, f"2 trios in resolved manifest (got {len(manifest)})")

    # --- Step 0: CH_B inferred male ---
    qc = {r["trio_id"]: r for r in rows(os.path.join(W, "qc_report.tsv"))}
    check(qc.get("CH_B", {}).get("inferred_sex") == "1", "CH_B inferred male (chrX)")

    # --- Step 3: plausible sites keep/drop ---
    plaus = {}
    for v in VCF(os.path.join(W, "plausible.sites.vcf.gz")):
        plaus[(v.CHROM, v.POS)] = v.INFO.get("hprv_keep_reason")
    check(("chr2", 8000) in plaus and plaus[("chr2", 8000)] == "clinvar_plp",
          "ClinVar P/LP (LOW impact) kept via clinvar_plp")
    check(("chr1", 12000) not in plaus, "BA1-common variant dropped at Step 3")
    check(("chr1", 17000) not in plaus, "non-PASS variant dropped before Step 3")
    check(("chr1", 5000) in plaus, "de novo site retained as plausible")

    # --- Step 5: inheritance calls ---
    calls = rows(os.path.join(W, "candidates.calls.tsv"))

    def has(trio, mode, chrom=None, pos=None, gene=None):
        for r in calls:
            if r["trio_id"] == trio and r["mode"] == mode \
               and (chrom is None or r["chrom"] == chrom) \
               and (pos is None or r["pos"] == str(pos)) \
               and (gene is None or r["symbol"] == gene):
                return True
        return False

    check(has("CH_A", "denovo", "chr1", 5000, "GENE1"), "CH_A de novo GENE1")
    check(has("CH_A", "hom_recessive", "chr1", 8000, "GENE2"), "CH_A hom recessive GENE2")
    ch = [r for r in calls if r["trio_id"] == "CH_A" and r["mode"] == "compound_het"]
    check(len(ch) == 2 and all(r["symbol"] == "GENE3" for r in ch), "CH_A compound het pair in GENE3")
    check(has("CH_B", "denovo", "chr1", 5100, "GENE1"), "CH_B de novo GENE1 (recurrent)")
    check(has("CH_B", "x_linked_recessive", "chrX", 2781600, "GENEX"), "CH_B X-linked recessive GENEX")
    check(not has("CH_A", "denovo", "chr1", 15000), "low-GQ pseudo-de-novo NOT called (QC gate)")

    # --- Step 6: recurrent-gene burden ---
    genes = {r["gene"]: r for r in rows(os.path.join(W, "genes.ranked.tsv"))}
    check(genes.get("GENE1", {}).get("obs_denovo") == "2", "GENE1 has 2 de novo across trios")

    # --- audit exists ---
    check(os.path.exists(os.path.join(W, "audit", "summary.md")), "audit/summary.md written")
    counts = rows(os.path.join(W, "audit", "counts.tsv"))
    check(any(r["step"] == "resolve" and r["metric"] == "trios_resolved" and r["value"] == "2"
              for r in counts), "audit records 2 resolved trios")

    if FAILS:
        sys.stderr.write(f"\n{len(FAILS)} assertion(s) FAILED\n")
        return 1
    print("\nALL INTEGRATION ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
