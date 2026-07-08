#!/usr/bin/env python3
"""Pipeline Step 3: select biologically-plausible sites from the annotated cohort union.

Inheritance-AGNOSTIC filter that shrinks the annotated cohort sites to a target list
worth genotyping per trio. Uses the permissive-union rarity gate (the looser of the
dominant/recessive cutoffs) so nothing any inheritance mode needs is dropped early,
keeps ClinVar P/LP as an override, and never rescues BA1-common variants. Gene lists
and constraint are NOT applied here (never-drop rule); they are downstream priors.

See docs/pipeline_design.md (Step 3) and docs/README.md#canonical-defaults.

Usage:
  03_select_plausible.py --in cohort.sites.annotated.vcf.gz --out plausible.sites.vcf.gz \
      --config config.yaml
"""
from __future__ import annotations

import argparse
import subprocess
import sys

from cyvcf2 import VCF, Writer

from hprv import audit
from hprv.config import load_config
from hprv.selection import build_classifier


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", required=True)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    classify = build_classifier(cfg)

    vcf = VCF(args.inp)
    vcf.add_info_to_header({
        "ID": "hprv_keep_reason", "Number": "1", "Type": "String",
        "Description": "Why this site was retained by Step 3 (evidence category)",
    })
    # cyvcf2 bgzips directly when writing a .gz path (mode "wz") — no uncompressed
    # intermediate to re-read, and no risk of aliasing the output onto its own input.
    mode = "wz" if args.out.endswith(".gz") else "w"
    w = Writer(args.out, vcf, mode=mode)
    n_in = n_out = 0
    reasons = {}
    for v in vcf:
        n_in += 1
        keep, reason = classify(v)
        reasons[reason] = reasons.get(reason, 0) + 1
        if keep:
            v.INFO["hprv_keep_reason"] = reason
            w.write_record(v)
            n_out += 1
    w.close()
    vcf.close()

    if args.out.endswith(".gz"):
        subprocess.run(["bcftools", "index", "-t", args.out], check=True)

    # audit: funnel in/out + counts by reason (drops and keeps)
    audit.record("03_select", "sites_in", n_in)
    audit.record("03_select", "sites_plausible", n_out)
    for r, c in sorted(reasons.items()):
        audit.record("03_select", f"reason.{r}", c)

    frac = (100.0 * n_out / n_in) if n_in else 0.0
    sys.stderr.write(
        f"Step 3 complete: {n_out}/{n_in} sites plausible ({frac:.1f}%) -> {args.out}\n"
        f"  reasons: {dict(sorted(reasons.items()))}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
