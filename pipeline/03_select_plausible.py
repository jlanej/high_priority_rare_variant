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

from hprv import annotations as A
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
    # HEADER GUARD for the faf95 arm. Under `oracle: faf95` the rarity value of every variant is
    # read from the gnomAD joint transfer; if that transfer never happened (a failed or skipped
    # Step 2 transfer, or an ingest run without the slim) frequency() returns None for EVERY
    # variant, every gate passes, and BA1-common alleles become candidates while the run exits 0.
    # Step 2 now halts on that too; this is the belt to its braces.
    if A.rarity_oracle(cfg) == "faf95":
        witness = A.F["gnomad_af_joint"]
        if f"##INFO=<ID={witness}," not in vcf.raw_header:
            sys.stderr.write(
                f"ERROR: resources.gnomad.oracle is 'faf95' but {args.inp} declares no INFO/"
                f"{witness} — the gnomAD joint slim was never transferred in Step 2, so every "
                "variant would read as ABSENT from gnomAD (= rarest) and no rarity gate would "
                "fire. Fix the slim (resources.gnomad.sites_slim) and re-run Step 2, or set "
                "oracle: grpmax_proxy deliberately.\n")
            return 1
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
    # `not_functional` is one bucket for two different facts — "the scores said no" and "no score
    # existed" — and the pipeline treats that distinction as load-bearing everywhere else (absence
    # is never zero). A novel non-coding indel that neither CADD's table nor the precomputed SpliceAI
    # set covers has NO keep-path by construction, and the audit used to file it beside a variant
    # both predictors scored and rejected. Sub-count it so a negative result can say which.
    n_unscored = 0
    for v in vcf:
        n_in += 1
        keep, reason = classify(v)
        reasons[reason] = reasons.get(reason, 0) + 1
        if reason == "not_functional" and A.cadd(v) is None and A.spliceai_ds(v) is None:
            n_unscored += 1
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
    if reasons.get("not_functional"):
        # both halves recorded, so a zero reads as "measured 0", never as "never emitted"
        audit.record("03_select", "reason.not_functional.unscored", n_unscored)
        audit.record("03_select", "reason.not_functional.scored",
                     reasons["not_functional"] - n_unscored)

    frac = (100.0 * n_out / n_in) if n_in else 0.0
    sys.stderr.write(
        f"Step 3 complete: {n_out}/{n_in} sites plausible ({frac:.1f}%) -> {args.out}\n"
        f"  reasons: {dict(sorted(reasons.items()))}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
