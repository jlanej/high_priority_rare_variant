"""Build a variants.tsv for the jlanej/igv.js trio variant-review server.

The server (server/README.md) requires only chrom/pos/ref/alt; any other column
becomes a filterable annotation. Per-member alignment tracks come from
`child_/mother_/father_` `_file` (+ optional `_index`) columns, and a per-trio VCF
track from `*_vcf`/`*_vcf_index`/`*_vcf_id`. All track paths are RELATIVE to the
server's --data-dir. We emit our inheritance mode + rich annotations as extra
(filterable) columns, and point the track columns at the mini-CRAMs / VCFs that the
export step places under the data-dir (only when they exist).
"""

from __future__ import annotations

import csv
import os

# Output column order: chrom/pos/ref/alt first (required), then recommended +
# our annotations, then the per-member track columns.
COLUMNS = [
    "chrom", "pos", "ref", "alt",
    "trio_id", "gene", "consequence", "impact", "frequency", "inheritance", "origin", "pair_id",
    "child_gt", "mother_gt", "father_gt", "child_GQ", "child_DP", "child_AB",
    "faf95", "revel", "alphamissense", "spliceai_max", "clin_sig", "clinvar_stars",
    "loftee", "loftee_flags",
    "child_file", "child_index", "mother_file", "mother_index", "father_file", "father_index",
    "child_vcf", "child_vcf_index", "child_vcf_id",
    "mother_vcf", "mother_vcf_index", "mother_vcf_id",
    "father_vcf", "father_vcf_index", "father_vcf_id",
]


def _origin(flags):
    for tok in (flags or "").split(";"):
        if tok.startswith("origin="):
            return tok.split("=", 1)[1]
    return ""


def _read_samples(manifest):
    """trio_id -> (kid, dad, mom) from the resolved manifest's `samples` column."""
    out = {}
    with open(manifest) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = (r.get("samples") or "").split(",")
            if r.get("trio_id") and len(s) == 3:
                out[r["trio_id"]] = (s[0], s[1], s[2])
    return out


def build_variants_tsv(calls_tsv, manifest, data_dir, out_tsv):
    """Write variants.tsv. Track columns are populated only when the referenced
    mini-CRAM / VCF exists under data_dir (so the file always loads)."""
    samples = _read_samples(manifest)

    def rel_cram(trio, sample):
        rel = os.path.join("crams", trio, f"{sample}.cram")
        return rel if os.path.exists(os.path.join(data_dir, rel)) else ""

    def rel_vcf(trio):
        rel = os.path.join("vcfs", f"{trio}.vcf.gz")
        return rel if os.path.exists(os.path.join(data_dir, rel)) else ""

    n = 0
    with open(calls_tsv) as fh, open(out_tsv, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        for r in csv.DictReader(fh, delimiter="\t"):
            trio = r.get("trio_id", "")
            kid, dad, mom = samples.get(trio, ("", "", ""))
            cf, mf, ff = rel_cram(trio, kid), rel_cram(trio, mom), rel_cram(trio, dad)
            vcf = rel_vcf(trio)
            row = {
                "chrom": r.get("chrom"), "pos": r.get("pos"), "ref": r.get("ref"), "alt": r.get("alt"),
                "trio_id": trio, "gene": r.get("symbol") or r.get("gene"),
                "consequence": r.get("consequence"), "impact": r.get("impact"),
                "frequency": r.get("faf95") or r.get("grpmax_af"),
                "inheritance": r.get("mode"), "origin": _origin(r.get("flags")),
                "pair_id": r.get("pair_id"),
                "child_gt": r.get("child_gt"), "mother_gt": r.get("mother_gt"),
                "father_gt": r.get("father_gt"), "child_GQ": r.get("child_gq"),
                "child_DP": r.get("child_dp"), "child_AB": r.get("child_ab"),
                "faf95": r.get("faf95"), "revel": r.get("revel"),
                "alphamissense": r.get("alphamissense"), "spliceai_max": r.get("spliceai_max"),
                "clin_sig": r.get("clnsig"), "clinvar_stars": r.get("clinvar_stars"),
                "loftee": r.get("loftee"), "loftee_flags": r.get("loftee_flags"),
                "child_file": cf, "child_index": (cf + ".crai") if cf else "",
                "mother_file": mf, "mother_index": (mf + ".crai") if mf else "",
                "father_file": ff, "father_index": (ff + ".crai") if ff else "",
                "child_vcf": vcf, "child_vcf_index": (vcf + ".tbi") if vcf else "", "child_vcf_id": kid,
                "mother_vcf": vcf, "mother_vcf_index": (vcf + ".tbi") if vcf else "", "mother_vcf_id": mom,
                "father_vcf": vcf, "father_vcf_index": (vcf + ".tbi") if vcf else "", "father_vcf_id": dad,
            }
            w.writerow(row)
            n += 1
    return n


def write_sample_qc(qc_report, manifest, out_tsv):
    """Emit sample_qc.tsv (trio_id, role, sample_id, + metrics) from Step 0 QC."""
    samples = _read_samples(manifest)
    rows = []
    if qc_report and os.path.exists(qc_report):
        with open(qc_report) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                tid = r.get("trio_id")
                if tid not in samples:
                    continue
                kid, dad, mom = samples[tid]
                mie = r.get("mie_rate", "")            # trio-level metric
                cflag = r.get("contam_flag", "")       # trio-level flag
                # inferred_sex is the PROBAND's chrX inference — it does not apply to the
                # parents; contamination is per-member (kid/mom/dad columns).
                per_role = (
                    ("proband", kid, r.get("inferred_sex", ""), r.get("kid_contam", "")),
                    ("mother", mom, "", r.get("mom_contam", "")),
                    ("father", dad, "", r.get("dad_contam", "")),
                )
                for role, sid, sex, contam in per_role:
                    rows.append({
                        "trio_id": tid, "role": role, "sample_id": sid, "mie_rate": mie,
                        "inferred_sex": sex, "contam": contam, "contam_flag": cflag,
                    })
    with open(out_tsv, "w", newline="") as out:
        w = csv.DictWriter(
            out, fieldnames=["trio_id", "role", "sample_id", "mie_rate", "inferred_sex",
                             "contam", "contam_flag"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return len(rows)
