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
from hprv.config import get, load_config


def _f(cfg, key, default):
    v = get(cfg, key, default)
    return float(v) if v is not None else default


def build_predicate(cfg):
    ba1 = _f(cfg, "filters.rarity.benign_ba1", 0.05)
    rec_max = _f(cfg, "filters.rarity.recessive_max", 1.0e-2)  # permissive-union cutoff
    revel_sup = _f(cfg, "filters.functional.revel_pp3_supporting", 0.644)
    am_lp = _f(cfg, "filters.functional.alphamissense_lp", 0.564)
    spliceai_pp3 = _f(cfg, "filters.functional.spliceai_pp3", 0.2)
    cadd_sup = _f(cfg, "filters.functional.cadd_phred_supporting", 20.0)
    mpc_strong = _f(cfg, "filters.functional.mpc_strong", 2.0)
    keep_impacts = set(get(cfg, "filters.functional.keep_impacts", ["HIGH", "MODERATE"]))

    def is_functional(v) -> bool:
        if (A.impact(v) or "") in keep_impacts:
            return True
        if A.is_loftee_hc(v):
            return True
        for val, thr in (
            (A.spliceai_max(v), spliceai_pp3),
            (A.revel(v), revel_sup),
            (A.alphamissense(v), am_lp),
            (A.cadd(v), cadd_sup),
            (A.mpc(v), mpc_strong),
        ):
            if val is not None and val >= thr:
                return True
        return False

    def keep(v) -> bool:
        fr = A.frequency(v)
        if fr is not None and fr >= ba1:            # ClinGen BA1 — never rescue
            return False
        plp = A.clnsig_is_plp(v) and A.clinvar_stars(v) >= 1
        rarity_ok = (fr is None) or (fr < rec_max) or plp
        if not rarity_ok:
            return False
        return plp or is_functional(v)

    return keep


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", required=True)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    keep = build_predicate(cfg)

    vcf = VCF(args.inp)
    out = args.out[:-3] if args.out.endswith(".gz") else args.out  # write plain, bgzip after
    w = Writer(out, vcf)
    n_in = n_out = 0
    for v in vcf:
        n_in += 1
        if keep(v):
            w.write_record(v)
            n_out += 1
    w.close()
    vcf.close()

    # bgzip + index via bcftools (present in the container)
    subprocess.run(["bcftools", "view", "-Oz", "-o", args.out, out], check=True)
    subprocess.run(["bcftools", "index", "-t", args.out], check=True)
    subprocess.run(["rm", "-f", out], check=False)

    frac = (100.0 * n_out / n_in) if n_in else 0.0
    sys.stderr.write(
        f"Step 3 complete: {n_out}/{n_in} sites plausible ({frac:.1f}%) -> {args.out}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
