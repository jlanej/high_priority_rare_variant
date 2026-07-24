"""Step-2b SpliceAI live-backfill helpers (pure VCF I/O; the TensorFlow scoring is external).

The precomputed raw SpliceAI files cover every genome-wide SNV, but their INDEL set is finite and
their window is narrow (default ``-D 50``). So some candidate variants — mostly **novel indels** —
carry no precomputed delta score. This module supports scoring exactly those, live, with the stock
Illumina ``spliceai`` model (run in an isolated conda env), and folding the result into the SAME
``vep_SpliceAI_pred_DS_*`` INFO fields the SpliceAI VEP plugin produces, so ``selection.py`` reads
precomputed and backfilled scores identically.

Two steps, each a small cyvcf2 pass (the heavy part, the model inference, is the external
``spliceai`` call between them):

  1. ``select`` — write the subset of the annotated union that has NO precomputed SpliceAI score
     (``annotations.spliceai_ds`` is None), optionally restricted to indels (SNVs are fully covered).
  2. ``annotate`` — read the ``spliceai``-scored subset (INFO ``SpliceAI=…``) and emit a small
     annotation VCF carrying ``vep_SpliceAI_pred_DS_{AG,AL,DG,DL}`` (per-field MAX over gene entries),
     which Step 2b then ``bcftools annotate``s onto the union.

Keep-only semantics are preserved: a variant the model still cannot score simply stays unscored
(None), which never drops it.
"""

from __future__ import annotations

import argparse
import sys

from hprv import annotations as A

# The four SpliceAI delta-score event fields, in the CSQ order the plugin/`SpliceAI=` string uses.
DS_EVENTS = ("DS_AG", "DS_AL", "DS_DG", "DS_DL")
DS_FIELDS = tuple(f"vep_SpliceAI_pred_{e}" for e in DS_EVENTS)


def _is_indel(v) -> bool:
    """True if the (biallelic, post norm -m-) record is an indel — SNVs are precomputed-complete."""
    alt = v.ALT[0] if v.ALT else ""
    return bool(alt) and not alt.startswith(("<", "*")) and len(v.REF) != len(alt)


def select_unscored(union_vcf: str, out_vcf: str, indels_only: bool = True) -> int:
    """Write records of `union_vcf` that carry no precomputed SpliceAI score to `out_vcf`.

    Returns the number written. `indels_only` (default) skips SNVs, which the precomputed raw set
    already scores everywhere — so the backfill set is the small, tractable novel-indel gap.
    """
    from cyvcf2 import VCF, Writer

    vcf = VCF(union_vcf)
    w = Writer(out_vcf, vcf)
    n = 0
    for v in vcf:
        if A.spliceai_ds(v) is not None:      # already has a (precomputed) score — skip
            continue
        if indels_only and not _is_indel(v):  # SNVs are precomputed-complete
            continue
        w.write_record(v)
        n += 1
    w.close()
    vcf.close()
    return n


def _max_ds(spliceai_value):
    """Per-event max delta score over all `SpliceAI=` gene entries → dict event->str, or None.

    `SpliceAI=` is `ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`, comma-joined
    across genes. We keep the strongest signal per event (the standard reduction)."""
    if not spliceai_value:
        return None
    best = {e: None for e in DS_EVENTS}
    seen = False
    for entry in str(spliceai_value).split(","):
        parts = entry.split("|")
        if len(parts) < 6:
            continue
        for i, e in enumerate(DS_EVENTS):
            tok = parts[2 + i].strip()
            try:
                f = float(tok)
            except ValueError:
                continue
            best[e] = f if best[e] is None else max(best[e], f)
            seen = True
    return best if seen else None


def to_annotation(scored_vcf: str, out_vcf: str) -> int:
    """Convert a `spliceai`-scored VCF (INFO `SpliceAI=`) into an annotation VCF carrying
    `vep_SpliceAI_pred_DS_*` (per-field max), ready for `bcftools annotate -c INFO/...`.

    Returns the number of records that received a score. Records the model did not score are
    dropped from the annotation (nothing to transfer)."""
    from cyvcf2 import VCF, Writer

    vcf = VCF(scored_vcf)
    for field, ev in zip(DS_FIELDS, DS_EVENTS):
        vcf.add_info_to_header({
            "ID": field, "Number": ".", "Type": "String",
            "Description": f"SpliceAI {ev} delta score (live backfill; max over genes)",
        })
    w = Writer(out_vcf, vcf)
    n = 0
    for v in vcf:
        best = _max_ds(v.INFO.get("SpliceAI"))
        if best is None:
            continue
        for field, ev in zip(DS_FIELDS, DS_EVENTS):
            if best[ev] is not None:
                v.INFO[field] = f"{best[ev]:.4g}"
        w.write_record(v)
        n += 1
    w.close()
    vcf.close()
    return n


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    ps = sub.add_parser("select", help="write the unscored subset of the annotated union")
    ps.add_argument("--in", dest="inp", required=True)
    ps.add_argument("--out", required=True)
    ps.add_argument("--indels-only", action="store_true", default=False)
    pa = sub.add_parser("annotate", help="SpliceAI=-scored VCF -> vep_SpliceAI_pred_DS_* annotation VCF")
    pa.add_argument("--in", dest="inp", required=True)
    pa.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    if args.cmd == "select":
        n = select_unscored(args.inp, args.out, indels_only=args.indels_only)
        sys.stderr.write(f"spliceai_backfill: {n} unscored variant(s) selected -> {args.out}\n")
    else:
        n = to_annotation(args.inp, args.out)
        sys.stderr.write(f"spliceai_backfill: {n} record(s) scored -> {args.out}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
