"""Lightweight run auditing: append (step, scope, metric, value) rows.

Every step records its input/output counts and funnel tallies to a single
``counts.tsv`` under ``$HPRV_AUDIT_DIR`` so the whole run is answerable: how many
samples/trios/variants entered and left each step, and per-trio breakdowns. The
orchestrator assembles a human-readable summary from this file.

`scope` is "global" for cohort-wide metrics or a trio_id for per-trio metrics.
"""

from __future__ import annotations

import os
from datetime import datetime, timezone

_HEADER = "timestamp\tstep\tscope\tmetric\tvalue\n"


def audit_dir(explicit=None):
    return explicit or os.environ.get("HPRV_AUDIT_DIR")


def record(step, metric, value, scope="global", adir=None):
    adir = audit_dir(adir)
    if not adir:
        return
    os.makedirs(adir, exist_ok=True)
    path = os.path.join(adir, "counts.tsv")
    new = not os.path.exists(path)
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(path, "a") as fh:
        if new:
            fh.write(_HEADER)
        fh.write(f"{ts}\t{step}\t{scope}\t{metric}\t{value}\n")


def _read(adir):
    """Return {(step, scope, metric): value} (last value wins for re-runs)."""
    path = os.path.join(adir, "counts.tsv")
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == 5:
                out[(f[1], f[2], f[3])] = f[4]
    return out


def summarize(adir, out_md=None):
    """Assemble a human-readable run summary (global funnel + per-trio) from counts.tsv."""
    d = _read(adir)

    def g(step, metric, scope="global"):
        return d.get((step, scope, metric), "")

    lines = ["# Run audit summary", ""]
    lines += ["## Trio resolution",
              f"- trios in pedigree: {g('resolve','trios_input')}",
              f"- resolved to a VCF: {g('resolve','trios_resolved')}",
              f"- unresolved: {g('resolve','trios_unresolved')}  "
              f"(matched >1 VCF: {g('resolve','trios_multi_vcf')})",
              f"- VCFs scanned: {g('resolve','vcfs_scanned')}; samples indexed: {g('resolve','samples_indexed')}",
              "  (per-trio detail in trio_resolution.tsv)", ""]
    lines += ["## Global variant funnel",
              f"- cohort union sites: {g('01_cohort_sites','union_sites')}",
              f"- annotated sites: {g('02_annotate','annotated_sites')}",
              f"- plausible sites: {g('03_select','sites_plausible')} "
              f"(of {g('03_select','sites_in')} in)", ""]
    # step-3 drop/keep reasons
    reasons = sorted((m, v) for (s, sc, m), v in d.items()
                     if s == "03_select" and m.startswith("reason."))
    if reasons:
        lines.append("### Step 3 reasons (keeps + drops)")
        for m, v in reasons:
            lines.append(f"- {m[len('reason.'):]}: {v}")
        lines.append("")
    # per-trio table
    trios = sorted({sc for (s, sc, m) in d if s in ("04_subset", "05_inheritance") and sc != "global"})
    if trios:
        lines += ["## Per-trio funnel", "",
                  "| trio | candidate genotypes | candidate calls | modes |",
                  "|------|--------------------:|----------------:|-------|"]
        for t in trios:
            cg = g("04_subset", "candidate_genotypes", t)
            cc = g("05_inheritance", "candidate_calls", t)
            modes = ", ".join(f"{m[len('mode.'):]}={v}" for (s, sc, m), v in sorted(d.items())
                              if s == "05_inheritance" and sc == t and m.startswith("mode."))
            lines.append(f"| {t} | {cg} | {cc} | {modes} |")
        lines.append("")
    lines += ["## Cross-pedigree gene burden",
              f"- genes nominated: {g('06_burden','genes_nominated')}",
              f"- exome-wide significant: {g('06_burden','genes_exome_wide_sig')}; "
              f"FDR significant: {g('06_burden','genes_fdr_sig')}", ""]

    text = "\n".join(lines)
    if out_md:
        with open(out_md, "w") as fh:
            fh.write(text + "\n")
    return text


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Summarize an hprv run audit.")
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    print(summarize(a.dir, a.out or None))
