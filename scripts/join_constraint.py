#!/usr/bin/env python3
"""Join the per-gene constraint sources into ONE TSV the pipeline's Step 6 reads.

Left-joins (on gene symbol) gnomAD v2.1.1 LOEUF/pLI + Zeng-2024 s_het + Collins-2022 pHaplo
into columns: gene, oe_lof_upper, pli, s_het, phaplo. These are priors/tiers only (never hard
filters). Column names vary across releases, so each value column is located case-insensitively
from a set of aliases. Missing inputs are skipped (join what's available); join-miss is reported.

Usage:
  join_constraint.py --gnomad lof_metrics.txt.bgz [--shet shet.tsv] [--phaplo phaplo.tsv.gz] --out constraint.by_gene.tsv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import sys


def _open(path):
    if path.endswith((".gz", ".bgz")):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path)


def _sniff_reader(fh):
    head = fh.readline()
    delim = "\t" if "\t" in head else ("," if "," in head else "\t")
    fh.seek(0)
    return csv.DictReader(fh, delimiter=delim)


def _col(fieldnames, *aliases):
    low = {c.lower().lstrip("#").strip(): c for c in (fieldnames or [])}
    for a in aliases:
        if a in low:
            return low[a]
    return None


def _load(path, gene_aliases, val_aliases):
    """Return {gene_symbol: value_str} from `path`, or {} if the file is absent/unreadable."""
    if not path:
        return {}
    try:
        with _open(path) as fh:
            r = _sniff_reader(fh)
            gk = _col(r.fieldnames, *gene_aliases)
            vk = _col(r.fieldnames, *val_aliases)
            if not gk or not vk:
                sys.stderr.write(f"WARN: {path}: gene/value column not found "
                                 f"(have {r.fieldnames}); skipping\n")
                return {}
            out = {}
            for row in r:
                g = (row.get(gk) or "").strip()
                v = (row.get(vk) or "").strip()
                if g and v not in ("", "NA", "."):
                    out.setdefault(g, v)  # first wins
            return out
    except OSError as e:
        sys.stderr.write(f"WARN: cannot read {path}: {e}\n")
        return {}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gnomad", required=True, help="gnomad.v2.1.1.lof_metrics.by_gene.txt(.bgz)")
    ap.add_argument("--shet", default="")
    ap.add_argument("--phaplo", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    loeuf = _load(a.gnomad, ("gene", "gene_symbol", "symbol"), ("oe_lof_upper", "loeuf"))
    pli = _load(a.gnomad, ("gene", "gene_symbol", "symbol"), ("pli",))
    shet = _load(a.shet, ("gene", "hgnc", "symbol", "gene_symbol"), ("s_het", "shet", "post_mean", "s_het_mean"))
    phaplo = _load(a.phaplo, ("gene", "gene_symbol", "symbol"), ("phaplo", "phaplo_score"))

    genes = set(loeuf) | set(pli) | set(shet) | set(phaplo)
    if not genes:
        sys.stderr.write("ERROR: no constraint records loaded from any source\n")
        return 1

    with open(a.out, "w") as out:
        out.write("gene\toe_lof_upper\tpli\ts_het\tphaplo\n")
        for g in sorted(genes):
            out.write("\t".join([g, loeuf.get(g, ""), pli.get(g, ""),
                                  shet.get(g, ""), phaplo.get(g, "")]) + "\n")

    sys.stderr.write(
        f"constraint join -> {a.out}: {len(genes)} genes "
        f"(LOEUF {len(loeuf)}, pLI {len(pli)}, s_het {len(shet)}, pHaplo {len(phaplo)})\n")
    # Loud WARN if a supplied source contributed ZERO genes overlapping the LOEUF gene set — the
    # usual cause is a key mismatch (e.g. GeneBayes keys s_het on ENSG, not HGNC symbol), which
    # would otherwise silently produce an all-empty column that Step 6 treats as "no constraint".
    loeuf_genes = set(loeuf)
    for name, path, table in (("s_het", a.shet, shet), ("pHaplo", a.phaplo, phaplo)):
        if path and table and loeuf_genes and not (set(table) & loeuf_genes):
            sys.stderr.write(
                f"WARN: {name} ({path}) loaded {len(table)} rows but 0 overlap the LOEUF gene set — "
                f"likely a gene-key mismatch (ENSG vs symbol). The {name} column will be empty; "
                f"map keys to HGNC symbols or point at a symbol-keyed file.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
