"""Build a variants.tsv for the jlanej/igv.js trio variant-review server.

The server (server/README.md) requires only chrom/pos/ref/alt; any other column
becomes a filterable annotation. Per-member alignment tracks come from
`child_/mother_/father_` `_file` (+ optional `_index`) columns, and a per-trio VCF
track from `*_vcf`/`*_vcf_index`/`*_vcf_id`. All track paths are RELATIVE to the
server's --data-dir. We emit our inheritance mode + rich annotations as extra
(filterable) columns, and point the track columns at the mini-CRAMs / VCFs that the
export step places under the data-dir (only when they exist). Every calls column not
represented by a curated column is appended verbatim after the track columns
(`passthrough_columns`), so the review table is dropless with respect to candidates.calls.tsv.

Step-8b non-human fraction (NHF): when nonhuman-screen has classified a member's
ALT-supporting reads, its per-(trio, sample) table lands at
`<data_dir>/nhf/<trio>/<sample>.variant_nhf.tsv`. We fold the fraction (and its
read-count denominator) in as `child_/mother_/father_nhf` columns. See the join
note above `build_variants_tsv`.
"""

from __future__ import annotations

import csv
import os

# Flag a call when this fraction of some screened member's ALT-supporting reads classify
# non-human (over >= min_reads reads). Both numbers come from the config
# (outputs.igv.nonhuman_screen.flag_fraction / .min_reads) so `nhf_flag` here and Step 9's
# `nhf_status` read the SAME threshold; this constant is only the fallback default.
NHF_FLAG_FRACTION = 0.5

# Output column order: chrom/pos/ref/alt first (required), then recommended +
# our annotations, then the per-member track columns.
COLUMNS = [
    "chrom", "pos", "ref", "alt",
    # `frequency` IS `rarity_af` — the value the run's oracle produced and every gate applied.
    # Duplicated into this prominent position on purpose: it is the column a reviewer reads
    # first, and it must not be a different quantity from the one that did the filtering.
    # `rarity_oracle` / `rarity_basis` say which quantity that was and how it arose; `grpmax_af`
    # and `faf95` ride further right as the raw inputs.
    "trio_id", "gene", "consequence", "impact", "hgvsc", "hgvsp", "frequency", "inheritance",
    "origin", "pair_id",
    "child_gt", "mother_gt", "father_gt", "child_GQ", "child_DP", "child_AB",
    # max_af/max_af_pops are shown next to grpmax_af so a reviewer can spot a call whose
    # frequency is driven by a founder group grpmax excludes (see annotations.GRPMAX_POPS).
    # faf95 (CI-corrected, the oracle when the gnomAD slim is transferred) beside grpmax_af (the
    # point-estimate proxy), so a reviewer can see the gap. faf95_group flags the `mid` case.
    "rarity_af", "rarity_oracle", "rarity_basis", "grpmax_af", "faf95", "faf95_group",
    "nhomalt",
    "max_af", "max_af_pops", "cadd", "spliceai_ds",
    # Calibrated missense predictors (REVEL / AlphaMissense plugins). Filterable in igv.js like
    # every other extra column. Blank on non-missense is EXPECTED — these are missense-only
    # scores, not a coverage gap.
    "revel", "alphamissense", "alphamissense_class",
    # clin_sig is the VEP cache's CLIN_SIG. clinvar_stars is the ClinVar VCF's review status
    # (0-4) and is BLANK when the transfer did not run — blank is "nobody looked", NOT 0 stars.
    # Without it a 1-star single-submitter assertion and a 3-star expert-panel one look alike.
    "clin_sig", "clinvar_stars",
    # Step-8b NHF: fraction of each member's ALT-supporting reads that classify NON-human, with
    # its read-count denominator right beside it (an NHF over few reads is noise). nhf_flag is a
    # single convenience boolean (>= NHF_FLAG_FRACTION over >= min_reads in any screened member).
    "child_nhf", "child_nhf_reads", "mother_nhf", "mother_nhf_reads",
    "father_nhf", "father_nhf_reads", "nhf_flag",
    "child_file", "child_index", "mother_file", "mother_index", "father_file", "father_index",
    "child_vcf", "child_vcf_index", "child_vcf_id",
    "mother_vcf", "mother_vcf_index", "mother_vcf_id",
    "father_vcf", "father_vcf_index", "father_vcf_id",
]

# --- dropless pass-through of the calls table ----------------------------------------------------
# build_variants_tsv used to write ONLY the fixed COLUMNS above (DictWriter extrasaction="ignore"),
# which made the review table a second silent projection on top of Step 5's: `flags` beyond
# origin (origin_unverified, transmitting_parent_qc_fail, unphased_denovo_partner...),
# hiConfDeNovo, review_prior_crosscheck and the Ensembl `gene` never reached igv.js. Now every
# calls column this module does not already represent (verbatim or renamed — the set below) is
# appended AFTER the track columns under its own name, in calls-header order. Step 9's
# igv/variants.prioritized.tsv copies every input column verbatim, so it inherits all of it.
CONSUMED_CALLS_COLUMNS = frozenset({
    "chrom", "pos", "ref", "alt", "trio_id", "symbol", "consequence", "impact", "hgvsc", "hgvsp",
    "rarity_af", "mode", "pair_id", "child_gt", "mother_gt", "father_gt",
    "child_gq", "child_dp", "child_ab", "grpmax_af", "max_af", "rarity_oracle", "rarity_basis",
    "faf95", "faf95_group", "nhomalt", "max_af_pops", "cadd", "spliceai_ds",
    "revel", "alphamissense", "alphamissense_class", "clnsig", "clinvar_stars",
})
# candidates.calls.tsv `gene` is the Ensembl gene ID; variants.tsv `gene` is the SYMBOL (the
# igv.js server's gene column), so the ID rides along under a name that says what it is.
RENAMED_CALLS_COLUMNS = {"gene": "gene_id"}


def passthrough_columns(calls_fieldnames):
    """-> [(calls_column, variants_column), ...] for every calls column not already represented.

    Calls-header order. A name that would collide with a curated COLUMNS entry is prefixed
    `calls_` rather than dropped or allowed to overwrite the curated value."""
    out, taken = [], set(COLUMNS)
    for c in calls_fieldnames or ():
        if not c or c in CONSUMED_CALLS_COLUMNS:
            continue
        name = RENAMED_CALLS_COLUMNS.get(c, c)
        if name in taken:
            name = "calls_" + name
        taken.add(name)
        out.append((c, name))
    return out


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


def _load_nhf_tsv(path):
    """variant_key -> (nonhuman_fraction, supporting_reads) for one member's Step-8b table.

    Returns None if the file does not exist (member not screened), else a dict (possibly
    empty). `variant_key` is nonhuman-screen's 0-based "{chrom}:{pos0}:{ref}:{alt}"; the
    caller joins on pos-1. See the note above build_variants_tsv."""
    if not path or not os.path.exists(path):
        return None
    m = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = r.get("variant_key")
            if k:
                m[k] = (r.get("nonhuman_fraction", ""), r.get("supporting_reads", ""))
    return m


def build_variants_tsv(calls_tsv, manifest, data_dir, out_tsv, nhf_dir=None, nhf_min_reads=5,
                       nhf_flag_fraction=NHF_FLAG_FRACTION):
    """Write variants.tsv. Track columns are populated only when the referenced
    mini-CRAM / VCF exists under data_dir (so the file always loads).

    NHF join (Step 8b): nonhuman-screen writes one row per concrete ALT allele to
    `<nhf_dir>/<trio>/<sample>.variant_nhf.tsv`, keyed 0-based
    ("{chrom}:{pos0}:{ref}:{alt}"). variants.tsv `pos` is 1-based (cyvcf2 v.POS), so the
    join is on **pos-1** — the single load-bearing off-by-one. A member with no table
    (not a carrier / no mini-CRAM / NHF disabled) gets blank NHF columns, which is distinct
    from a 0.0 fraction ("screened, all human"). Symbolic ALTs (e.g. the `*` spanning
    deletion) are skipped by nonhuman-screen, so those rows also get blanks. nhf_dir=None
    disables the whole join (all NHF columns blank), keeping legacy behavior byte-identical."""
    samples = _read_samples(manifest)
    nhf_cache = {}   # (trio, sample) -> map-or-None, so each member's table is read once

    def rel_cram(trio, sample):
        rel = os.path.join("crams", trio, f"{sample}.cram")
        return rel if os.path.exists(os.path.join(data_dir, rel)) else ""

    def rel_vcf(trio):
        rel = os.path.join("vcfs", f"{trio}.vcf.gz")
        return rel if os.path.exists(os.path.join(data_dir, rel)) else ""

    def nhf_map(trio, sample):
        if not nhf_dir or not sample:
            return None
        key = (trio, sample)
        if key not in nhf_cache:
            nhf_cache[key] = _load_nhf_tsv(
                os.path.join(nhf_dir, trio, f"{sample}.variant_nhf.tsv"))
        return nhf_cache[key]

    def nhf_flag(screened):
        """1 if any screened member's non-human fraction is >= nhf_flag_fraction over >=
        nhf_min_reads reads; else 0 if anything was screened for this variant; else '' (nothing
        screened)."""
        if not screened:
            return ""
        for frac, reads in screened:
            try:
                if float(frac) >= float(nhf_flag_fraction) and int(reads) >= nhf_min_reads:
                    return "1"
            except (TypeError, ValueError):
                continue
        return "0"

    n = 0
    with open(calls_tsv) as fh, open(out_tsv, "w", newline="") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        # Dropless: every calls column not represented below rides through verbatim after the
        # track columns (passthrough_columns). The curated COLUMNS keep their positions.
        extra = passthrough_columns(reader.fieldnames)
        w = csv.DictWriter(out, fieldnames=COLUMNS + [name for _, name in extra], delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in reader:
            trio = r.get("trio_id", "")
            kid, dad, mom = samples.get(trio, ("", "", ""))
            cf, mf, ff = rel_cram(trio, kid), rel_cram(trio, mom), rel_cram(trio, dad)
            vcf = rel_vcf(trio)

            # NHF join on the 0-based key. key0 is None (all blank) if any coordinate field is
            # missing/unparseable — never let a bad row raise.
            chrom, pos, ref, alt = r.get("chrom"), r.get("pos"), r.get("ref"), r.get("alt")
            key0 = None
            if chrom and pos not in (None, "") and ref and alt:
                try:
                    key0 = f"{chrom}:{int(pos) - 1}:{ref}:{alt}"
                except (TypeError, ValueError):
                    key0 = None

            screened = []   # (frac, reads) for members whose table CONTAINS this allele
            member_nhf = {}   # role -> (frac, reads); blank when not screened / allele absent
            for role, sample in (("child", kid), ("mother", mom), ("father", dad)):
                m = nhf_map(trio, sample)
                hit = m.get(key0) if (m is not None and key0 is not None) else None
                member_nhf[role] = hit if hit is not None else ("", "")
                if hit is not None:
                    screened.append(hit)

            row = {
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "trio_id": trio, "gene": r.get("symbol") or r.get("gene"),
                "consequence": r.get("consequence"), "impact": r.get("impact"),
                "hgvsc": r.get("hgvsc"), "hgvsp": r.get("hgvsp"),
                # THE ORACLE'S value — the same number as `rarity_af`, shown early because it
                # is the frequency a reviewer reads first. It was hardwired to `grpmax_af`, so
                # under the default faf95 oracle the headline column disagreed with the value
                # every gate actually used, and the review table showed two different
                # "frequencies" with no way to tell which was authoritative.
                # Keyed on rarity_ORACLE, not on rarity_af: the latter is legitimately EMPTY when
                # the oracle has no value for the allele (basis=absent), and keying on it would
                # send exactly those rows down the legacy branch.
                "frequency": (r.get("rarity_af") if r.get("rarity_oracle")
                              else r.get("grpmax_af")),
                "inheritance": r.get("mode"), "origin": _origin(r.get("flags")),
                "pair_id": r.get("pair_id"),
                "child_gt": r.get("child_gt"), "mother_gt": r.get("mother_gt"),
                "father_gt": r.get("father_gt"), "child_GQ": r.get("child_gq"),
                "child_DP": r.get("child_dp"), "child_AB": r.get("child_ab"),
                "grpmax_af": r.get("grpmax_af"), "max_af": r.get("max_af"),
                "rarity_af": r.get("rarity_af"), "rarity_oracle": r.get("rarity_oracle"),
                "rarity_basis": r.get("rarity_basis"),
                "faf95": r.get("faf95"), "faf95_group": r.get("faf95_group"),
                "nhomalt": r.get("nhomalt"),
                "max_af_pops": r.get("max_af_pops"), "cadd": r.get("cadd"),
                "spliceai_ds": r.get("spliceai_ds"),
                "revel": r.get("revel"), "alphamissense": r.get("alphamissense"),
                "alphamissense_class": r.get("alphamissense_class"),
                # renamed on the way through (clnsig -> clin_sig); clinvar_stars keeps its name
                "clin_sig": r.get("clnsig"), "clinvar_stars": r.get("clinvar_stars"),
                "child_nhf": member_nhf["child"][0], "child_nhf_reads": member_nhf["child"][1],
                "mother_nhf": member_nhf["mother"][0], "mother_nhf_reads": member_nhf["mother"][1],
                "father_nhf": member_nhf["father"][0], "father_nhf_reads": member_nhf["father"][1],
                "nhf_flag": nhf_flag(screened),
                "child_file": cf, "child_index": (cf + ".crai") if cf else "",
                "mother_file": mf, "mother_index": (mf + ".crai") if mf else "",
                "father_file": ff, "father_index": (ff + ".crai") if ff else "",
                "child_vcf": vcf, "child_vcf_index": (vcf + ".tbi") if vcf else "", "child_vcf_id": kid,
                "mother_vcf": vcf, "mother_vcf_index": (vcf + ".tbi") if vcf else "", "mother_vcf_id": mom,
                "father_vcf": vcf, "father_vcf_index": (vcf + ".tbi") if vcf else "", "father_vcf_id": dad,
            }
            for src, name in extra:
                row[name] = r.get(src, "")
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
