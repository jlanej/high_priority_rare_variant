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
    # Step 0: the sex evidence a methods section has to quote — who decided each proband's sex,
    # how the chrX inference calibrated, and the audit-only chrY coverage check beside it.
    lines += ["## Sample QC (Step 0)",
              f"- trios QC'd: {g('00_qc','trios_qc')}; flagged (advisory): {g('00_qc','trios_flagged')}",
              f"- proband sex source: pedigree {g('00_qc','trios_sex_source.ped')}, chrX inference "
              f"{g('00_qc','trios_sex_source.inferred')}, unresolved {g('00_qc','trios_sex_source.none')}; "
              f"pedigree-vs-chrX discordant: {g('00_qc','trios_sex_discordant_inference')}",
              f"- chrX calibration: fathers inferred male {g('00_qc','fathers_inferred_male')} / female "
              f"{g('00_qc','fathers_inferred_female')}; mothers inferred female "
              f"{g('00_qc','mothers_inferred_female')} / male {g('00_qc','mothers_inferred_male')} "
              f"(median het ratio fathers {g('00_qc','father_x_het_ratio_median')}, mothers "
              f"{g('00_qc','mother_x_het_ratio_median')})",
              f"- chrY coverage (audit-only): probands reading male {g('00_qc','y_sex_inferred_male')}, "
              f"female {g('00_qc','y_sex_inferred_female')}, no call {g('00_qc','y_sex_inferred_none')}; "
              f"contradicting the pedigree {g('00_qc','trios_sex_match_y_discordant')}, contradicting "
              f"the chrX inference {g('00_qc','trios_xy_disagree')}, mothers with chrY coverage "
              f"{g('00_qc','trios_mother_y_coverage')}; median proband y_cov_ratio PED-male "
              f"{g('00_qc','proband_y_cov_ratio_median.ped_male')} / PED-female "
              f"{g('00_qc','proband_y_cov_ratio_median.ped_female')}, mothers "
              f"{g('00_qc','mother_y_cov_ratio_median')}",
              f"- chrY parent-role check (raw haploid coverage, a male ~1 / a female ~0): fathers "
              f"reading male {g('00_qc','fathers_y_male')} / female {g('00_qc','fathers_y_female')}, "
              f"mothers reading female {g('00_qc','mothers_y_female')} / male "
              f"{g('00_qc','mothers_y_male')} (median father {g('00_qc','father_y_cov_haploid_median')}, "
              f"mother {g('00_qc','mother_y_cov_haploid_median')}); trios with a parent contradicted: "
              f"{g('00_qc','trios_parent_sex_flag_y')}", ""]
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
    # per-trio table. "call rows" can EXCEED "examined" (a compound-het leg is emitted once per
    # pair, and one variant can appear under up to three modes), which is why the input side is
    # shown beside it: examined == skipped + with a call + no row, per trio, and "no row" is the
    # count a reader of a negative result actually needs.
    trios = sorted({sc for (s, sc, m) in d if s in ("04_subset", "05_inheritance") and sc != "global"})
    if trios:
        lines += ["## Per-trio funnel", "",
                  "| trio | candidate genotypes | examined | with a call | no row | skipped | call rows | modes |",
                  "|------|--------------------:|---------:|------------:|-------:|--------:|----------:|-------|"]
        for t in trios:
            cg = g("04_subset", "candidate_genotypes", t)
            ex = g("05_inheritance", "variants_examined", t)
            wc = g("05_inheritance", "variants_with_call", t)
            nr = g("05_inheritance", "variants_no_row", t)
            sk = sum(int(v) for (s, sc, m), v in d.items()
                     if s == "05_inheritance" and sc == t and m.startswith("skipped."))
            cc = g("05_inheritance", "candidate_calls", t)
            modes = ", ".join(f"{m[len('mode.'):]}={v}" for (s, sc, m), v in sorted(d.items())
                              if s == "05_inheritance" and sc == t and m.startswith("mode."))
            lines.append(f"| {t} | {cg} | {ex} | {wc} | {nr} | {sk} | {cc} | {modes} |")
        lines.append("")
    # why examined variants produced no row (all trios) — the drop side of Step 5, by reason
    no_row = sorted((m, v) for (s, sc, m), v in d.items()
                    if s == "05_inheritance" and sc == "global" and m.startswith("no_row."))
    skipped = sorted((m, v) for (s, sc, m), v in d.items()
                     if s == "05_inheritance" and sc == "global" and m.startswith("skipped."))
    if no_row or skipped:
        lines.append("### Step 5: examined variants that produced no row (all trios)")
        for m, v in no_row:
            lines.append(f"- {m[len('no_row.'):]}: {v}")
        for m, v in skipped:
            lines.append(f"- skipped ({m[len('skipped.'):]}): {v}")
        plp = g("05_inheritance", "clinvar_plp_no_row")
        if plp != "":
            lines.append(f"- ClinVar P/LP alleles the child carried that yielded no row: {plp}")
        lines.append("")
    lines += ["## Cross-pedigree gene burden",
              f"- genes nominated: {g('06_burden','genes_nominated')}; "
              f"recurrent: {g('06_burden','genes_recurrent')}",
              f"- recurrence exome-wide significant: {g('06_burden','genes_recurrence_exome_wide_sig')}; "
              f"FDR significant: {g('06_burden','genes_recurrence_fdr_sig')} "
              "(case-only null — a RANK, not a calibrated test; see docs/gene_burden.md)",
              f"- genes ranked by the size-normalised p_carrier_excess: "
              f"{g('06_burden','genes_rank_mu_normalised')}", ""]

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
