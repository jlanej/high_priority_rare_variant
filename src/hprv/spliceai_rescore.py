"""Step 5b — wide-window SpliceAI rescoring of the CALLED set (file I/O only; the model runs outside).

The screen runs on Illumina's PRECOMPUTED SpliceAI scores, which were generated with a 50 bp window
(``-D 50``): a cryptic site or pseudoexon partner further than 50 bp from the variant is invisible to
them. That is the right trade for a genome-wide gate and the wrong one for a reviewer holding a
short list. Step 5b re-scores every variant in ``candidates.calls.tsv`` LIVE with the stock model at
a wide window (default ``-D 4999``, the window Walker et al. 2023 used) and appends the result as
``spliceai_wide_*`` columns beside the precomputed ones. The screen and the Step-9 tier keep reading
the precomputed score, so results stay comparable run to run; the wide columns ADD evidence, they
never move a gate.

Three design points, each load-bearing:

* **The expensive part is cached per VARIANT, not per file.** Step 5 recomputes and rewrites
  ``candidates.calls.tsv`` on every invocation, so a marker on that file would re-score everything
  every time. Scores land in ``<outdir>/scores.tsv`` keyed by chrom/pos/ref/alt (with the window
  they were scored at in its header); ``plan`` emits only the variants that file does not already
  hold, and ``gather`` re-merges the cache into whatever calls table Step 5 just wrote. A changed
  window invalidates the whole cache — a 50 bp score and a 4999 bp score are different quantities.
* **Work is chunked so it can be a SLURM array.** ``plan`` writes ``chunks/chunk_NNNN.vcf`` plus a
  manifest; each chunk is scored independently (``05b_spliceai_rescore.sh --chunk N``), and
  ``gather`` refuses to merge while any chunk in the manifest is unscored — never a partial merge
  that looks complete.
* **Plain text, no cyvcf2.** The inputs and outputs here are tiny (the called set, not the union), so
  the VCFs are written and read as text. That keeps this module runnable — and testable — on a bare
  host with nothing but python3.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

from hprv import audit
from hprv import splice

# The window the Illumina precomputed files were generated with. A wide-window event further than
# this from the variant is one the screen could not have seen — the `spliceai_wide_distal` flag.
PRECOMPUTED_WINDOW = 50

WIDE_COLUMNS = (
    "spliceai_wide_ds", "spliceai_wide_event", "spliceai_wide_event_pos",
    "spliceai_wide_event2", "spliceai_wide_event2_ds", "spliceai_wide_event2_pos",
    "spliceai_wide_shift_nt", "spliceai_wide_shift_frame", "spliceai_wide_effect",
    "spliceai_wide_symbol", "spliceai_wide_distal", "spliceai_wide_minus_precomputed",
)
SCORE_COLUMNS = ("chrom", "pos", "ref", "alt", "symbol",
                 "ds_ag", "ds_al", "ds_dg", "ds_dl", "dp_ag", "dp_al", "dp_dg", "dp_dl", "n_genes")
_DS_COL = {e: "ds_" + splice.EVENT_SUFFIX[e].lower() for e in splice.EVENTS}
_DP_COL = {e: "dp_" + splice.EVENT_SUFFIX[e].lower() for e in splice.EVENTS}


def _fmt(x):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.4g}"
    return str(x)


def read_tsv(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return [], []
    return rows[0], rows[1:]


def variant_keys(header, rows):
    """Distinct (chrom, pos, ref, alt) in first-seen order; symbolic / missing alleles skipped."""
    idx = {c: i for i, c in enumerate(header)}
    for c in ("chrom", "pos", "ref", "alt"):
        if c not in idx:
            raise ValueError(f"calls table lacks a '{c}' column (have {header})")
    seen, out = set(), []
    for r in rows:
        try:
            k = (r[idx["chrom"]], r[idx["pos"]], r[idx["ref"]], r[idx["alt"]])
        except IndexError:
            continue
        if not all(k) or k[3].startswith(("*", "<")) or k[3] == "." or "," in k[3]:
            continue
        if k not in seen:
            seen.add(k)
            out.append(k)
    return out


# --- the per-variant score cache -------------------------------------------------------------
def load_scores(path, distance):
    """-> (scores {key: row}, status) with status in {absent, ok, distance_mismatch}.

    A cache scored at a different window is NOT reused: the header records the distance, and a
    mismatch returns an empty dict so every variant is re-planned (the old file is replaced on the
    next gather)."""
    if not path or not os.path.exists(path):
        return {}, "absent"
    with open(path) as fh:
        first = fh.readline().rstrip("\n")
        cached_d = None
        for tok in first.lstrip("#").split():
            if tok.startswith("distance="):
                cached_d = tok.split("=", 1)[1]
        if cached_d is None or str(cached_d) != str(distance):
            return {}, "distance_mismatch"
        reader = csv.DictReader(fh, delimiter="\t")
        out = {}
        for r in reader:
            out[(r["chrom"], r["pos"], r["ref"], r["alt"])] = r
    return out, "ok"


def write_scores(path, scores, distance):
    tmp = path + ".tmp"
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(tmp, "w", newline="") as fh:
        fh.write(f"#spliceai_rescore\tdistance={distance}\tmodel=spliceai\tannotation=grch38\n")
        w = csv.DictWriter(fh, fieldnames=list(SCORE_COLUMNS), delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for k in sorted(scores, key=lambda k: (k[0], int(k[1]) if str(k[1]).isdigit() else 0, k[2], k[3])):
            w.writerow(scores[k])
    os.replace(tmp, path)


# --- planning: the pending variants, chunked into minimal VCFs -------------------------------
def _contig_lines(fai):
    if not fai or not os.path.exists(fai):
        return []
    lines = []
    with open(fai) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 2:
                lines.append(f"##contig=<ID={f[0]},length={f[1]}>")
    return lines


def write_chunks(keys, outdir, chunk_size, fai=None):
    """Write `chunks/chunk_NNNN.vcf` (minimal VCFs the live model accepts) and `manifest.txt`
    (one chunk path per line). Returns the chunk paths; an empty list writes an EMPTY manifest,
    which is the documented 'nothing to do' signal for the SLURM planner."""
    cdir = os.path.join(outdir, "chunks")
    os.makedirs(cdir, exist_ok=True)
    for old in os.listdir(cdir):          # stale chunks from a previous plan must not linger
        os.remove(os.path.join(cdir, old))
    header = ["##fileformat=VCFv4.2"] + _contig_lines(fai) + \
             ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    paths = []
    size = max(1, int(chunk_size))
    for i in range(0, len(keys), size):
        p = os.path.join(cdir, f"chunk_{i // size:04d}.vcf")
        with open(p, "w") as fh:
            fh.write("\n".join(header) + "\n")
            for chrom, pos, ref, alt in keys[i:i + size]:
                fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
        paths.append(p)
    with open(os.path.join(outdir, "manifest.txt"), "w") as fh:
        fh.write("".join(p + "\n" for p in paths))
    return paths


def scored_path(chunk_path):
    return chunk_path[:-4] + ".scored.vcf" if chunk_path.endswith(".vcf") else chunk_path + ".scored.vcf"


def parse_scored_vcf(path):
    """One score row per variant that received a `SpliceAI=` value in the live tool's output."""
    out = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = dict(kv.partition("=")[::2] for kv in f[7].split(";") if kv and kv != ".")
            ent = splice.parse_spliceai_info(info.get("SpliceAI"))
            if ent is None:
                continue
            row = {"chrom": f[0], "pos": f[1], "ref": f[3], "alt": f[4], "symbol": ent["symbol"],
                   "n_genes": ent["n_genes"]}
            for e in splice.EVENTS:
                row[_DS_COL[e]] = _fmt(ent["ds"][e])
                row[_DP_COL[e]] = _fmt(ent["dp"][e])
            out.append(row)
    return out


# --- the merge into candidates.calls.tsv ----------------------------------------------------
def wide_values(score_row, pos, floor, precomputed_ds=None):
    """The `spliceai_wide_*` cells for one call from its cached score row."""
    ds = {e: score_row.get(_DS_COL[e]) for e in splice.EVENTS}
    dp = {e: score_row.get(_DP_COL[e]) for e in splice.EVENTS}
    d = splice.decompose(ds, dp, pos=pos, floor=floor)
    mx = splice.max_ds(ds)
    distal = ""
    if d["event"] is not None:
        distal = "1" if (d["dp"] is not None and abs(d["dp"]) > PRECOMPUTED_WINDOW) else "0"
    pre = splice._num(precomputed_ds)
    return {
        "spliceai_wide_ds": _fmt(mx), "spliceai_wide_event": d["event"] or "",
        "spliceai_wide_event_pos": _fmt(d["pos"]), "spliceai_wide_event2": d["event2"] or "",
        "spliceai_wide_event2_ds": _fmt(d["event2_ds"]), "spliceai_wide_event2_pos": _fmt(d["event2_pos"]),
        "spliceai_wide_shift_nt": _fmt(d["shift_nt"]), "spliceai_wide_shift_frame": d["shift_frame"] or "",
        "spliceai_wide_effect": d["effect"] or "", "spliceai_wide_symbol": score_row.get("symbol", ""),
        "spliceai_wide_distal": distal,
        "spliceai_wide_minus_precomputed": _fmt(mx - pre) if (mx is not None and pre is not None) else "",
    }


def merge(calls_path, scores, floor, out_path=None):
    """Rewrite the calls table with the `spliceai_wide_*` block filled from `scores`.

    Idempotent: an existing block is replaced, never duplicated. The block is inserted after the
    last curated `spliceai_*` column (before the `info_*` pass-through block), so the splice
    evidence reads left to right: precomputed, then wide-window. Returns stats."""
    header, rows = read_tsv(calls_path)
    if not header:
        raise ValueError(f"{calls_path} is empty")
    wide = set(WIDE_COLUMNS)
    keep = [i for i, c in enumerate(header) if c not in wide]
    header = [header[i] for i in keep]
    rows = [[r[i] if i < len(r) else "" for i in keep] for r in rows]
    idx = {c: i for i, c in enumerate(header)}
    anchors = [i for i, c in enumerate(header) if c.startswith("spliceai_")]
    if anchors:
        at = anchors[-1] + 1
    else:
        infos = [i for i, c in enumerate(header) if c.startswith("info_")]
        at = infos[0] if infos else len(header)
    new_header = header[:at] + list(WIDE_COLUMNS) + header[at:]
    n_wide = n_distal = 0
    out_rows = []
    for r in rows:
        key = (r[idx["chrom"]], r[idx["pos"]], r[idx["ref"]], r[idx["alt"]])
        sr = scores.get(key)
        if sr is None:
            cells = {c: "" for c in WIDE_COLUMNS}
        else:
            try:
                pos = int(r[idx["pos"]])
            except ValueError:
                pos = None
            pre = r[idx["spliceai_ds"]] if "spliceai_ds" in idx else None
            cells = wide_values(sr, pos, floor, precomputed_ds=pre)
            n_wide += 1
            n_distal += cells["spliceai_wide_distal"] == "1"
        out_rows.append(r[:at] + [cells[c] for c in WIDE_COLUMNS] + r[at:])
    target = out_path or calls_path
    tmp = target + ".tmp"
    with open(tmp, "w", newline="") as fh:
        fh.write("\t".join(new_header) + "\n")
        for r in out_rows:
            fh.write("\t".join(r) + "\n")
    os.replace(tmp, target)
    return {"rows": len(out_rows), "rows_with_wide_score": n_wide, "rows_distal_event": n_distal}


# --- CLI --------------------------------------------------------------------------------------
def _plan(a) -> int:
    header, rows = read_tsv(a.calls)
    keys = variant_keys(header, rows)
    scores, status = load_scores(os.path.join(a.outdir, "scores.tsv"), a.distance)
    if status == "distance_mismatch":
        sys.stderr.write(f"spliceai_rescore: cached scores were made at a different window than "
                         f"-D {a.distance}; ignoring the cache (every variant will be re-scored)\n")
    pending = [k for k in keys if k not in scores]
    paths = write_chunks(pending, a.outdir, a.chunk_size, fai=a.fai)
    audit.record("05b_spliceai_rescore", "variants_distinct", len(keys))
    audit.record("05b_spliceai_rescore", "variants_cached", len(keys) - len(pending))
    audit.record("05b_spliceai_rescore", "variants_pending", len(pending))
    audit.record("05b_spliceai_rescore", "chunks", len(paths))
    sys.stderr.write(f"spliceai_rescore: {len(keys)} distinct variants, {len(keys) - len(pending)} "
                     f"cached at -D {a.distance}, {len(pending)} pending in {len(paths)} chunk(s)\n")
    print(len(paths))
    return 0


def _gather(a) -> int:
    scores_path = os.path.join(a.outdir, "scores.tsv")
    scores, status = load_scores(scores_path, a.distance)
    manifest = os.path.join(a.outdir, "manifest.txt")
    chunks = []
    if os.path.exists(manifest):
        with open(manifest) as fh:
            chunks = [ln.strip() for ln in fh if ln.strip()]
    missing = [c for c in chunks if not os.path.exists(scored_path(c))]
    if missing:
        sys.stderr.write(f"ERROR: spliceai_rescore gather: {len(missing)} of {len(chunks)} chunk(s) "
                         f"have no scored output — the merge would be PARTIAL. Score them first "
                         f"(first missing: {scored_path(missing[0])}).\n")
        return 1
    n_new = 0
    for c in chunks:
        for row in parse_scored_vcf(scored_path(c)):
            scores[(row["chrom"], row["pos"], row["ref"], row["alt"])] = row
            n_new += 1
    if chunks or status != "ok":
        write_scores(scores_path, scores, a.distance)
    stats = merge(a.calls, scores, a.floor)
    for c in chunks:                                  # the cache now holds the results
        for p in (c, scored_path(c), scored_path(c) + ".done"):
            if os.path.exists(p):
                os.remove(p)
    if os.path.exists(manifest):
        os.remove(manifest)
    audit.record("05b_spliceai_rescore", "variants_scored_this_run", n_new)
    audit.record("05b_spliceai_rescore", "variants_in_cache", len(scores))
    for k, v in stats.items():
        audit.record("05b_spliceai_rescore", k, v)
    sys.stderr.write(f"spliceai_rescore: merged {n_new} newly scored variant(s); cache holds "
                     f"{len(scores)}; {stats['rows_with_wide_score']}/{stats['rows']} call rows carry a "
                     f"wide-window score, {stats['rows_distal_event']} with an event beyond "
                     f"+/-{PRECOMPUTED_WINDOW} bp\n")
    return 0


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "gather"):
        p = sub.add_parser(name)
        p.add_argument("--calls", required=True)
        p.add_argument("--outdir", required=True)
        p.add_argument("--distance", type=int, required=True)
        if name == "plan":
            p.add_argument("--chunk-size", type=int, default=1000)
            p.add_argument("--fai", default="")
        else:
            p.add_argument("--floor", type=float, default=0.2)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    return _plan(a) if a.cmd == "plan" else _gather(a)


if __name__ == "__main__":
    raise SystemExit(_main())
