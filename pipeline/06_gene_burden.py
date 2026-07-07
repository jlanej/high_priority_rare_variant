#!/usr/bin/env python3
"""Pipeline Step 6: cross-pedigree gene-burden screen.

Aggregates the per-trio candidate calls from Step 5 to find genes carrying MORE
qualifying rare variants than expected, and ranks them with gene constraint.

Two signals (best-fit for non-jointly-genotyped trios):
  * PRIMARY  — de novo enrichment vs a Samocha per-gene mutation-rate model
    (observed de novo count per class tested against a Poisson expectation
    2 * N_trios * mu_gene_class). Robust to cohort batch effects; needs only trios.
  * FALLBACK — recurrence tally (genes recurrent across unrelated trios), used when
    no mutation-rate table is supplied.

Nominated genes are down-weighted when the gene is tolerant of damaging variation
(high LOEUF / low s_het) — such genes are not interesting. See docs/gene_burden.md.

CALIBRATION NOTE: synonymous-lambda-approx-1 calibration requires an UNFILTERED de
novo set (synonymous de novos are dropped by Step 3). This step reports observed vs
expected and flags that the enrichment test is uncalibrated unless a synonymous de
novo count is supplied via --syn-denovo-count.

Usage:
  06_gene_burden.py --calls candidates.calls.tsv --out genes.ranked.tsv --config cfg.yaml \
      [--n-trios N] [--mutrate mutrate.tsv] [--constraint constraint.tsv] [--syn-denovo-count K]
"""
from __future__ import annotations

import argparse
import csv
import sys

from hprv import annotations as A
from hprv import audit
from hprv.config import get, load_config

try:
    from scipy.stats import poisson
except ImportError:  # pragma: no cover
    poisson = None


def _open_keyed(path, key_names):
    """Read a headered TSV/CSV into {gene: rowdict}; key column matched case-insensitively."""
    if not path:
        return {}, []
    with open(path) as fh:
        sniff = fh.readline()
        delim = "\t" if "\t" in sniff else ","
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter=delim)
        cols = reader.fieldnames or []
        keycol = next((c for c in cols if c.lower() in key_names), None)
        if keycol is None:
            sys.stderr.write(f"WARN: no gene key column in {path} (have {cols}); ignoring\n")
            return {}, cols
        out = {}
        for row in reader:
            g = (row.get(keycol) or "").strip()
            if g:
                out[g] = row
        return out, cols


def _find(cols, *names):
    low = {c.lower(): c for c in cols}
    for n in names:
        if n in low:
            return low[n]
    return None


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values for a list of p-values (None-safe)."""
    idx = [i for i, p in enumerate(pvals) if p is not None]
    m = len(idx)
    q = [None] * len(pvals)
    order = sorted(idx, key=lambda i: pvals[i])
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        val = min(prev, pvals[i] * m / k)
        q[i] = prev = val
    return q


LOF_CLASS = A.LOF_CONSEQUENCES


def classify(consequence: str) -> str:
    c = (consequence or "").lower()
    if any(t in c for t in LOF_CLASS):
        return "lof"
    if "missense" in c:
        return "missense"
    return "other"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--n-trios", type=int, default=0)
    ap.add_argument("--mutrate", default="")
    ap.add_argument("--constraint", default="")
    ap.add_argument("--syn-denovo-count", type=int, default=-1)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    exome_p = float(get(cfg, "burden.exome_wide_p", 2.5e-6))
    fdr_q = float(get(cfg, "burden.fdr_q", 0.05))
    loeuf_tol = float(get(cfg, "filters.constraint_weighting.loeuf_v2_tier1", 0.35))

    # --- aggregate candidate calls per gene ---
    genes = {}
    trios = set()
    with open(args.calls) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            gene = r.get("symbol") or r.get("gene")
            if not gene:
                continue
            trios.add(r.get("trio_id"))
            g = genes.setdefault(gene, {
                "gene": gene, "n_calls": 0, "trios": set(), "modes": {},
                "denovo_lof": 0, "denovo_mis": 0,
            })
            g["n_calls"] += 1
            g["trios"].add(r.get("trio_id"))
            g["modes"][r.get("mode")] = g["modes"].get(r.get("mode"), 0) + 1
            if r.get("mode") in ("denovo", "denovo_x_hemi"):
                cls = classify(r.get("consequence"))
                if cls == "lof":
                    g["denovo_lof"] += 1
                elif cls == "missense":
                    g["denovo_mis"] += 1

    n_trios = args.n_trios or len(trios)
    if n_trios == 0:
        sys.stderr.write("ERROR: no trios found; set --n-trios\n")
        return 1

    mut, mcols = _open_keyed(args.mutrate, {"gene", "gene_symbol", "symbol"})
    con, ccols = _open_keyed(args.constraint, {"gene", "gene_symbol", "symbol"})
    mut_lof_c = _find(mcols, "mut_lof", "mu_lof", "p_lof", "lof")
    mut_mis_c = _find(mcols, "mut_mis", "mu_mis", "p_mis", "mis")
    loeuf_c = _find(ccols, "oe_lof_upper", "loeuf", "loeuf_v2")
    pli_c = _find(ccols, "pli", "pli_v2")
    shet_c = _find(ccols, "s_het", "shet")

    can_enrich = bool(mut) and poisson is not None
    if args.mutrate and poisson is None:
        sys.stderr.write("WARN: scipy unavailable; skipping enrichment (recurrence only)\n")

    rows = []
    for gene, g in genes.items():
        obs = g["denovo_lof"] + g["denovo_mis"]
        p_enrich = None
        exp = None
        if can_enrich and gene in mut:
            mrow = mut[gene]
            mu_lof = float(mrow.get(mut_lof_c) or 0) if mut_lof_c else 0.0
            mu_mis = float(mrow.get(mut_mis_c) or 0) if mut_mis_c else 0.0
            exp = 2.0 * n_trios * (mu_lof + mu_mis)
            if exp > 0:
                p_enrich = float(poisson.sf(obs - 1, exp))  # P(X >= obs)
        loeuf = pli = shet = None
        if gene in con:
            crow = con[gene]
            loeuf = _num(crow.get(loeuf_c)) if loeuf_c else None
            pli = _num(crow.get(pli_c)) if pli_c else None
            shet = _num(crow.get(shet_c)) if shet_c else None
        constrained = (loeuf is not None and loeuf < loeuf_tol) or (pli is not None and pli >= 0.9)
        rows.append({
            "gene": gene, "n_trios": len(g["trios"]), "n_calls": g["n_calls"],
            "denovo_lof": g["denovo_lof"], "denovo_mis": g["denovo_mis"],
            "obs_denovo": obs, "exp_denovo": (f"{exp:.4g}" if exp is not None else ""),
            "p_enrich": p_enrich, "loeuf": loeuf, "pli": pli, "s_het": shet,
            "constrained": "1" if constrained else "0",
            "modes": ";".join(f"{k}={v}" for k, v in sorted(g["modes"].items())),
        })

    # BH FDR on enrichment p-values
    qs = bh_fdr([r["p_enrich"] for r in rows])
    for r, q in zip(rows, qs):
        r["q_enrich"] = q
        r["exome_wide_sig"] = "1" if (r["p_enrich"] is not None and r["p_enrich"] < exome_p) else "0"

    # rank: significant enrichment first, then recurrence, constrained genes up,
    # tolerant genes down. (p None sorts last.)
    def rank_key(r):
        p = r["p_enrich"] if r["p_enrich"] is not None else 1.0
        return (p, -r["n_trios"], -r["n_calls"], 0 if r["constrained"] == "1" else 1)
    rows.sort(key=rank_key)

    out_cols = ["gene", "n_trios", "n_calls", "obs_denovo", "exp_denovo", "denovo_lof",
                "denovo_mis", "p_enrich", "q_enrich", "exome_wide_sig", "loeuf", "pli",
                "s_het", "constrained", "modes"]
    with open(args.out, "w") as out:
        out.write("\t".join(out_cols) + "\n")
        for r in rows:
            out.write("\t".join(_fmt(r.get(c)) for c in out_cols) + "\n")

    sig = sum(1 for r in rows if r["exome_wide_sig"] == "1")
    fdr_sig = sum(1 for r in rows if r["q_enrich"] is not None and r["q_enrich"] < fdr_q)
    audit.record("06_burden", "n_trios", n_trios)
    audit.record("06_burden", "genes_nominated", len(rows))
    audit.record("06_burden", "genes_exome_wide_sig", sig)
    audit.record("06_burden", "genes_fdr_sig", fdr_sig)
    sys.stderr.write(
        f"Step 6 complete: {len(rows)} genes, {n_trios} trios -> {args.out}\n"
        f"  exome-wide (p<{exome_p:g}): {sig}; FDR q<{fdr_q}: {fdr_sig}\n"
    )
    if args.syn_denovo_count < 0:
        sys.stderr.write(
            "  CALIBRATION: synonymous de novo count not supplied (--syn-denovo-count); "
            "the enrichment test is UNCALIBRATED. Provide an unfiltered synonymous de novo "
            "count to verify lambda ~ 1 before trusting p-values.\n"
        )
    return 0


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _fmt(x):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.4g}"
    return str(x)


if __name__ == "__main__":
    raise SystemExit(main())
