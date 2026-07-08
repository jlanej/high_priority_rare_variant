#!/usr/bin/env python3
"""Pipeline Step 6: cross-pedigree gene consolidation (recurrence-based).

Finds genes where rare, functional variants RECUR across multiple independent
individuals — the signal that a gene is of interest. Emphasis is on INHERITED
variation:

  * DOMINANT model — a rare functional HETEROZYGOUS variant is interesting when it
    stacks up across individuals. We tally the number of distinct individuals
    carrying a qualifying dominant (inherited het) variant per gene.
  * RECESSIVE model — distinct individuals with a biallelic hit (homozygous or
    compound het), and X-linked recessive.
  * De novo counts are carried as a SECONDARY column only (dedicated de novo
    filtering/review lives in separate machinery).

Genes are ranked by recurrence (>= min_carriers distinct individuals) and weighted by
gene constraint — a recurrent het in a constraint-intolerant (haploinsufficient) gene
is far more compelling than one in a tolerant gene. An OPTIONAL de novo Poisson
enrichment vs a Samocha mutation model is reported when a mutation-rate table is given.

See docs/gene_burden.md.

Usage:
  06_gene_burden.py --calls candidates.calls.tsv --out genes.ranked.tsv --config cfg.yaml \
      [--n-trios N] [--mutrate mutrate.tsv] [--constraint constraint.tsv]
"""
from __future__ import annotations

import argparse
import csv
import sys

from hprv import annotations as A
from hprv import audit
from hprv.config import get, load_config

try:
    from scipy.stats import binom, poisson
except ImportError:  # pragma: no cover
    binom = poisson = None

DOMINANT_MODES = {"dominant"}
BIALLELIC_MODES = {"hom_recessive", "compound_het"}
XLINKED_MODES = {"x_linked_recessive"}
DENOVO_MODES = {"denovo", "denovo_x_hemi"}


def _open_keyed(path, key_names):
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
        return {(r.get(keycol) or "").strip(): r for r in reader if (r.get(keycol) or "").strip()}, cols


def _find(cols, *names):
    low = {c.lower(): c for c in cols}
    for n in names:
        if n in low:
            return low[n]
    return None


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _fmt(x):
    if x is None:
        return ""
    return f"{x:.4g}" if isinstance(x, float) else str(x)


def _floored(q, floor):
    return q if (q is not None and q > 0) else floor


def p_carrier_hwe(fafs, floor, ploidy):
    """P(a random individual carries >=1 qualifying allele) under HWE.

    ploidy=2 for the autosomal dominant het model (>=1 of two alleles, ~2q per variant);
    ploidy=1 for the X-linked hemizygous-male model (a single allele, ~q per variant).
    """
    p_not = 1.0
    for q in fafs:
        p_not *= (1.0 - min(max(_floored(q, floor), 0.0), 1.0)) ** ploidy
    return 1.0 - p_not


def p_biallelic_hwe(fafs, floor):
    """P(a random individual is biallelic for this gene) ~ (sum of allele freqs)^2 under HWE.

    Approximates hom + compound-het carriage by the squared cumulative alt-allele frequency.
    Far smaller than the dominant >=1-allele probability, so recessive recurrence must NOT be
    tested against the dominant null.
    """
    s = min(sum(min(max(_floored(q, floor), 0.0), 1.0) for q in fafs), 1.0)
    return s * s


def bh_fdr(pvals):
    idx = [i for i, p in enumerate(pvals) if p is not None]
    m = len(idx)
    q = [None] * len(pvals)
    prev = 1.0
    for rank, i in enumerate(reversed(sorted(idx, key=lambda i: pvals[i])), start=1):
        prev = min(prev, pvals[i] * m / (m - rank + 1))
        q[i] = prev
    return q


def classify(consequence: str) -> str:
    c = (consequence or "").lower()
    if any(t in c for t in A.LOF_CONSEQUENCES):
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
    min_carriers = int(get(cfg, "burden.min_carriers", 2))
    exome_p = float(get(cfg, "burden.exome_wide_p", 2.5e-6))
    fdr_q = float(get(cfg, "burden.fdr_q", 0.05))
    loeuf_tol = float(get(cfg, "filters.constraint_weighting.loeuf_v2_tier1", 0.35))
    pli_min = float(get(cfg, "filters.constraint_weighting.pli_min", 0.9))
    shet_min = float(get(cfg, "filters.constraint_weighting.shet_min", 0.10))
    phaplo_min = float(get(cfg, "filters.constraint_weighting.phaplo_min", 0.86))
    weight_by_constraint = bool(get(cfg, "burden.weight_by_constraint", True))
    do_enrich = bool(get(cfg, "burden.denovo_enrichment", True))
    # For the recurrence null, an allele absent from gnomAD is floored at the detection
    # limit (~1 / 2*N_gnomAD alleles) so its expected carriers are tiny but non-zero.
    absent_floor = float(get(cfg, "burden.absent_faf95_floor", 1e-6))

    # --- aggregate distinct individuals per gene, by model ---
    genes, trios = {}, set()
    with open(args.calls) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gene = r.get("symbol") or r.get("gene")
            if not gene:
                continue
            trio, mode = r.get("trio_id"), r.get("mode")
            trios.add(trio)
            g = genes.setdefault(gene, {
                "dom": set(), "bi": set(), "x": set(), "dn": set(), "all": set(),
                "dom_faf": {}, "bi_faf": {}, "x_faf": {}, "denovo_lof": 0, "denovo_mis": 0,
            })
            g["all"].add(trio)
            # distinct qualifying variant per mode -> its gnomAD faf95, for the model-appropriate null
            key = f"{r.get('chrom')}:{r.get('pos')}:{r.get('ref')}:{r.get('alt')}"
            faf = _num(r.get("faf95"))
            if mode in DOMINANT_MODES:
                g["dom"].add(trio); g["dom_faf"][key] = faf
            elif mode in BIALLELIC_MODES:
                g["bi"].add(trio); g["bi_faf"][key] = faf
            elif mode in XLINKED_MODES:
                g["x"].add(trio); g["x_faf"][key] = faf
            elif mode in DENOVO_MODES:
                g["dn"].add(trio)
                cls = classify(r.get("consequence"))
                if cls == "lof":
                    g["denovo_lof"] += 1
                elif cls == "missense":
                    g["denovo_mis"] += 1

    # N_trios is the SCREENED population, not the subset with a call. Never infer it from the
    # calls file (that only contains trios with >=1 candidate) — that inflates significance.
    n_trios = args.n_trios or 0
    if not n_trios:
        sys.stderr.write("WARN: --n-trios not provided; recurrence null + de novo enrichment "
                         "SKIPPED (counts only). Pass the resolved-trio count for calibrated p-values.\n")

    mut, mcols = _open_keyed(args.mutrate, {"gene", "gene_symbol", "symbol"})
    con, ccols = _open_keyed(args.constraint, {"gene", "gene_symbol", "symbol"})
    mut_lof_c = _find(mcols, "mut_lof", "mu_lof", "p_lof", "lof")
    mut_mis_c = _find(mcols, "mut_mis", "mu_mis", "p_mis", "mis")
    loeuf_c = _find(ccols, "oe_lof_upper", "loeuf", "loeuf_v2")
    pli_c = _find(ccols, "pli", "pli_v2")
    shet_c = _find(ccols, "s_het", "shet")
    phaplo_c = _find(ccols, "phaplo", "phaplo_score")

    can_enrich = do_enrich and bool(mut) and poisson is not None and n_trios > 0

    rows = []
    for gene, g in genes.items():
        # Recurrence counts INHERITED models only (dominant het / biallelic / X-linked);
        # de novo is tracked separately (n_denovo) and never drives the recurrence flag.
        n_carriers = len(g["dom"] | g["bi"] | g["x"])

        # --- Calibrated recurrence null: is seeing this many distinct carriers surprising
        # given the gnomAD frequencies of the gene's qualifying variants? Each inheritance
        # model is tested against its OWN HWE null (a recessive/hemizygous carrier is NOT a
        # >=1-of-two-alleles event, so it must not be charged the dominant probability):
        #   dominant het  -> Binomial(N, 1 - prod_v (1-q_v)^2)      [PRIMARY headline signal]
        #   biallelic     -> Binomial(N, (sum_v q_v)^2)
        #   X-linked male -> Binomial(N, 1 - prod_v (1-q_v))        [hemizygous, single allele]
        # Only defined for >= min_carriers (a single observed carrier is not "recurrence" and
        # would be an ascertainment artifact). BH-FDR across genes on the primary p below.
        # (Case-only approximation using in-cohort variants; a gnomAD-derived per-gene
        # cumulative allele frequency, i.e. TRAPD/CoCoRV, is the natural upgrade.)
        def _recur(n, faf_map, prob):
            if binom is None or n_trios <= 0 or n < min_carriers or not faf_map:
                return None, None
            p = prob(list(faf_map.values()))
            if not p or p <= 0:
                return None, None
            return n_trios * p, float(binom.sf(n - 1, n_trios, p))

        exp_car, p_recurrence = _recur(len(g["dom"]), g["dom_faf"],
                                       lambda f: p_carrier_hwe(f, absent_floor, 2))
        _, p_rec_bi = _recur(len(g["bi"]), g["bi_faf"],
                             lambda f: p_biallelic_hwe(f, absent_floor))
        _, p_rec_x = _recur(len(g["x"]), g["x_faf"],
                            lambda f: p_carrier_hwe(f, absent_floor, 1))

        # optional SECONDARY de novo Poisson enrichment
        p_enrich = exp = None
        if can_enrich and gene in mut:
            mrow = mut[gene]
            # tolerate non-numeric cells ("NA", ".", "") instead of crashing on float()
            mu = ((_num(mrow.get(mut_lof_c)) or 0.0) if mut_lof_c else 0.0) + \
                 ((_num(mrow.get(mut_mis_c)) or 0.0) if mut_mis_c else 0.0)
            exp = 2.0 * n_trios * mu
            obs = g["denovo_lof"] + g["denovo_mis"]
            if exp > 0:
                p_enrich = float(poisson.sf(obs - 1, exp))
        loeuf = pli = shet = phaplo = None
        if gene in con:
            crow = con[gene]
            loeuf = _num(crow.get(loeuf_c)) if loeuf_c else None
            pli = _num(crow.get(pli_c)) if pli_c else None
            shet = _num(crow.get(shet_c)) if shet_c else None
            phaplo = _num(crow.get(phaplo_c)) if phaplo_c else None
        constrained = ((loeuf is not None and loeuf < loeuf_tol)
                       or (pli is not None and pli >= pli_min)
                       or (shet is not None and shet >= shet_min)
                       or (phaplo is not None and phaplo >= phaplo_min))
        modes = []
        if g["dom"]:
            modes.append(f"dominant={len(g['dom'])}")
        if g["bi"]:
            modes.append(f"biallelic={len(g['bi'])}")
        if g["x"]:
            modes.append(f"x_linked={len(g['x'])}")
        if g["dn"]:
            modes.append(f"denovo={len(g['dn'])}")
        rows.append({
            "gene": gene, "n_carriers": n_carriers, "n_dominant": len(g["dom"]),
            "n_biallelic": len(g["bi"]), "n_xlinked": len(g["x"]), "n_denovo": len(g["dn"]),
            "recurrent": "1" if n_carriers >= min_carriers else "0",
            "exp_carriers": exp_car, "p_recurrence": p_recurrence,
            "p_recurrence_biallelic": p_rec_bi, "p_recurrence_xlinked": p_rec_x,
            "loeuf": loeuf, "pli": pli, "s_het": shet, "phaplo": phaplo,
            "constrained": "1" if constrained else "0",
            "dn_exp": (f"{exp:.4g}" if exp is not None else ""), "dn_p_enrich": p_enrich,
            "modes": ";".join(modes),
        })

    # FDR across genes for BOTH the primary recurrence p-value and the secondary de novo one
    for r, q in zip(rows, bh_fdr([r["p_recurrence"] for r in rows])):
        r["q_recurrence"] = q
        r["recurrence_exome_wide_sig"] = "1" if (r["p_recurrence"] is not None and r["p_recurrence"] < exome_p) else "0"
    for r, q in zip(rows, bh_fdr([r["dn_p_enrich"] for r in rows])):
        r["dn_q_enrich"] = q
        r["dn_exome_wide_sig"] = "1" if (r["dn_p_enrich"] is not None and r["dn_p_enrich"] < exome_p) else "0"

    # Rank: recurrent genes first, then by the calibrated recurrence p-value (most
    # surprising first), then constraint, then counts, then (secondary) de novo enrichment.
    def rank_key(r):
        con_key = (r["constrained"] != "1") if weight_by_constraint else 0
        prec = r["p_recurrence"] if r["p_recurrence"] is not None else 1.0
        return (r["recurrent"] != "1", prec, con_key, -r["n_carriers"],
                -r["n_dominant"], -r["n_biallelic"],
                r["dn_p_enrich"] if r["dn_p_enrich"] is not None else 1.0)
    rows.sort(key=rank_key)

    out_cols = ["gene", "n_carriers", "n_dominant", "n_biallelic", "n_xlinked", "n_denovo",
                "recurrent", "exp_carriers", "p_recurrence", "q_recurrence",
                "recurrence_exome_wide_sig", "p_recurrence_biallelic", "p_recurrence_xlinked",
                "loeuf", "pli", "s_het", "phaplo", "constrained",
                "dn_exp", "dn_p_enrich", "dn_q_enrich", "dn_exome_wide_sig", "modes"]
    with open(args.out, "w") as out:
        out.write("\t".join(out_cols) + "\n")
        for r in rows:
            out.write("\t".join(_fmt(r.get(c)) for c in out_cols) + "\n")

    n_recurrent = sum(1 for r in rows if r["recurrent"] == "1")
    n_rec_con = sum(1 for r in rows if r["recurrent"] == "1" and r["constrained"] == "1")
    n_rec_sig = sum(1 for r in rows if r.get("recurrence_exome_wide_sig") == "1")
    n_rec_fdr = sum(1 for r in rows if r.get("q_recurrence") is not None and r["q_recurrence"] < fdr_q)
    audit.record("06_burden", "n_trios", n_trios)
    audit.record("06_burden", "genes_nominated", len(rows))
    audit.record("06_burden", "genes_recurrent", n_recurrent)
    audit.record("06_burden", "genes_recurrent_constrained", n_rec_con)
    audit.record("06_burden", "genes_recurrence_exome_wide_sig", n_rec_sig)
    audit.record("06_burden", "genes_recurrence_fdr_sig", n_rec_fdr)
    sys.stderr.write(
        f"Step 6 complete: {len(rows)} genes, {n_trios} trios -> {args.out}\n"
        f"  recurrent (>= {min_carriers} carriers): {n_recurrent}; recurrent+constrained: {n_rec_con}\n"
        f"  recurrence exome-wide sig (p<{exome_p:g}): {n_rec_sig}; FDR q<{fdr_q}: {n_rec_fdr}\n"
    )
    if args.mutrate and args.syn_denovo_count < 0:
        sys.stderr.write(
            "  NOTE: de novo enrichment (secondary) is uncalibrated without a synonymous "
            "de novo count (--syn-denovo-count); de novo review is handled by separate machinery.\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
