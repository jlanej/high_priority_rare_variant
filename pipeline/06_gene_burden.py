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

**The recurrence p-values are a RANK, not a calibrated test.** The null is case-only: it is built
from the frequencies of the variants OBSERVED in the cohort, so for carriers of private (absent)
variants it saturates — two carriers at N=200 give p≈3e-7, three at N=1000 give 4e-8 — and it is
monotone in the carrier count, i.e. in gene size. When the mutational-target table (`--mutrate`
with mu_mis/mu_syn/mu_lof) is supplied, Step 6 therefore also computes a SIZE-NORMALISED rank:
`exp_carriers_mu = C * mu_g` (C fit over the FULL table, zero-count genes included),
`carrier_excess_ratio` and a Poisson tail `p_carrier_excess`, and with
`burden.rank_by_mutational_target: true` (default) orders recurrent genes by that instead, so long
genes no longer lead by size. See docs/gene_burden.md.

Usage:
  06_gene_burden.py --calls candidates.calls.tsv --out genes.ranked.tsv --config cfg.yaml \
      [--n-trios N] [--n-male-trios N | --qc-report qc_report.tsv] [--mutrate mutrate.tsv]
      [--constraint constraint.tsv]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
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


def _open_text(path):
    """Open a TSV that may be bgzipped — the prepared mutational-target table IS (.bgz); a plain
    open() dies on byte 2 with UnicodeDecodeError (the Step-9 gotcha, now shared here)."""
    if path.endswith((".gz", ".bgz")):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path)


def _open_keyed(path, key_names):
    if not path:
        return {}, []
    with _open_text(path) as fh:
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
    ap.add_argument("--n-male-trios", type=int, default=0,
                    help="MALE proband count for the X-linked (hemizygous) null; defaults to the "
                         "count of inferred_sex==1 rows in --qc-report, else --n-trios")
    ap.add_argument("--qc-report", default="", help="Step 0 qc_report.tsv (male proband count)")
    ap.add_argument("--mutrate", default="",
                    help="Samocha per-gene rate table for the SECONDARY de novo enrichment")
    ap.add_argument("--mutational-target", default="",
                    help="per-gene mu_mis/mu_syn/mu_lof table (the gnomAD v2.1.1 constraint file "
                         "Step 9 also reads; .bgz ok) for the size-normalised recurrence rank. "
                         "Falls back to --mutrate's columns when omitted.")
    ap.add_argument("--constraint", default="")
    ap.add_argument("--syn-denovo-count", type=int, default=-1,
                    help="[reserved] observed synonymous de novo count for calibration "
                         "(not yet wired into the Poisson expectation)")
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
    absent_floor = float(get(cfg, "burden.absent_af_floor", 1e-6))
    rank_by_mu = bool(get(cfg, "burden.rank_by_mutational_target", True))
    # Same imputation Step 9 uses for a null mu_lof (median mu_lof / (mu_mis + mu_syn)).
    impute = float(get(cfg, "prioritization.excess.offset.mu_lof_impute_factor", 0.0516))

    # --- aggregate distinct individuals per gene, by model ---
    genes, trios = {}, set()
    # Count what came in and what this loop discards. A call with no gene attribution (an
    # intergenic variant Step 3's CADD rung kept) has nowhere to aggregate, so skipping it is
    # correct — but it used to be the only drop in the pipeline with no counter, and this step's
    # audit surface was otherwise outputs only (genes_nominated, ...).
    n_calls_in = n_no_gene = 0
    with open(args.calls) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n_calls_in += 1
            gene = r.get("symbol") or r.get("gene")
            if not gene:
                n_no_gene += 1
                continue
            trio, mode = r.get("trio_id"), r.get("mode")
            trios.add(trio)
            g = genes.setdefault(gene, {
                "dom": set(), "bi": set(), "x": set(), "dn": set(), "all": set(),
                "dom_faf": {}, "bi_faf": {}, "x_faf": {}, "denovo_lof": 0, "denovo_mis": 0,
                # per-trio SETS of variant keys, for same- vs distinct-variant recurrence
                "dom_sets": {}, "bi_sets": {}, "x_sets": {},
            })
            g["all"].add(trio)
            # distinct qualifying variant per mode -> its gnomAD frequency, for the model null.
            # rarity_af is THE run oracle's value (annotations.frequency()): faf95 by default, the
            # grpmax point-estimate proxy when the run opted down. On the proxy arm the value sits
            # ~one CI-width high on low-AC alleles, which makes the recurrence p slightly
            # CONSERVATIVE (a larger q inflates the null probability of seeing carriers).
            key = f"{r.get('chrom')}:{r.get('pos')}:{r.get('ref')}:{r.get('alt')}"
            # THE RUN'S ORACLE, not the proxy column. Step 5 writes rarity_af (the value every
            # gate actually used) alongside the raw grpmax_af; reading the raw column here meant a
            # faf95-oracle run GATED the screen on faf95 and then built the recurrence null on the
            # point estimate — two different quantities inside one result. Fall back to grpmax_af
            # only for a pre-existing calls table that predates the rarity_af column.
            faf = _num(r.get("rarity_af"))
            if faf is None and not r.get("rarity_oracle"):
                faf = _num(r.get("grpmax_af"))
            if mode in DOMINANT_MODES:
                g["dom"].add(trio); g["dom_faf"][key] = faf
                g["dom_sets"].setdefault(trio, set()).add(key)
            elif mode in BIALLELIC_MODES:
                g["bi"].add(trio); g["bi_faf"][key] = faf
                g["bi_sets"].setdefault(trio, set()).add(key)
            elif mode in XLINKED_MODES:
                g["x"].add(trio); g["x_faf"][key] = faf
                g["x_sets"].setdefault(trio, set()).add(key)
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
    # The X-linked (hemizygous-male) family is tested against the MALE proband count: a female
    # proband cannot be a hemizygous carrier, so N_trios overstates that denominator and the
    # p-value was conservative by the sex ratio. Read Step 0's inferred sex when available.
    n_male = args.n_male_trios or 0
    if not n_male and args.qc_report and __import__("os").path.exists(args.qc_report):
        with open(args.qc_report) as fh:
            n_male = sum(1 for r in csv.DictReader(fh, delimiter="\t")
                         if (r.get("inferred_sex") or "").strip() == "1")
    n_x = n_male if n_male > 0 else n_trios
    if n_trios and not n_male:
        sys.stderr.write("WARN: no male proband count (--n-male-trios / --qc-report); the X-linked "
                         "recurrence null uses N_trios as its denominator (conservative).\n")

    mut, mcols = _open_keyed(args.mutrate, {"gene", "gene_symbol", "symbol"})
    con, ccols = _open_keyed(args.constraint, {"gene", "gene_symbol", "symbol"})
    # The size-normalisation table: a dedicated --mutational-target when given (the gnomAD
    # v2.1.1 constraint table, the same file Step 9's offset reads), else --mutrate's own columns.
    if args.mutational_target:
        mtg, mtcols = _open_keyed(args.mutational_target, {"gene", "gene_symbol", "symbol"})
    else:
        mtg, mtcols = mut, mcols
    mt_mis_c = _find(mtcols, "mu_mis", "mut_mis", "p_mis", "mis")
    mt_syn_c = _find(mtcols, "mu_syn", "mut_syn", "p_syn", "syn")
    mt_lof_c = _find(mtcols, "mu_lof", "mut_lof", "p_lof", "lof")
    mut_lof_c = _find(mcols, "mut_lof", "mu_lof", "p_lof", "lof")
    mut_mis_c = _find(mcols, "mut_mis", "mu_mis", "p_mis", "mis")
    mut_syn_c = _find(mcols, "mut_syn", "mu_syn", "p_syn", "syn")
    loeuf_c = _find(ccols, "oe_lof_upper", "loeuf", "loeuf_v2")
    pli_c = _find(ccols, "pli", "pli_v2")
    shet_c = _find(ccols, "s_het", "shet")
    phaplo_c = _find(ccols, "phaplo", "phaplo_score")

    # RECORD WHICH COLUMN WON. `_find` returns the first present name over an alias chain, so a
    # table whose header differs from the canonical one is consumed silently: a bare `lof`/`mis`
    # spans denovolyzeR-shaped tables, and a gnomAD **v4** `oe_lof_upper` wins the LOEUF chain and
    # is then compared against a **v2-calibrated** cutoff (filters.constraint_weighting.
    # loeuf_v2_tier1 = 0.35). The output columns are fixed names, so nothing downstream could tell.
    for _lbl, _col in (("mut_lof", mut_lof_c), ("mut_mis", mut_mis_c), ("mut_syn", mut_syn_c),
                       ("loeuf", loeuf_c), ("pli", pli_c), ("s_het", shet_c), ("phaplo", phaplo_c)):
        if _col:
            audit.record("06_gene_burden", f"resolved_column.{_lbl}.{_col}", 1)
    _resolved = ", ".join(f"{l}={c}" for l, c in
                          (("mut_lof", mut_lof_c), ("mut_mis", mut_mis_c), ("mut_syn", mut_syn_c),
                           ("loeuf", loeuf_c), ("pli", pli_c), ("s_het", shet_c),
                           ("phaplo", phaplo_c)) if c)
    if _resolved:
        sys.stderr.write(f"Step 6: resolved constraint/mutrate columns: {_resolved}\n")
    # The LOEUF cutoff is calibrated on gnomAD v2.1.1. A v4 table uses the same column name for a
    # differently-scaled quantity, so the tier boundary would move without anything saying so.
    if loeuf_c and "v4" in (args.constraint or "").lower():
        sys.stderr.write(
            f"WARN: --constraint path names v4 but the LOEUF tier cutoff "
            f"(filters.constraint_weighting.loeuf_v2_tier1) is calibrated on gnomAD v2.1.1. "
            f"The column '{loeuf_c}' will be compared against a v2 boundary.\n")

    can_enrich = do_enrich and bool(mut) and poisson is not None and n_trios > 0

    # --- SIZE-NORMALISED recurrence rank (the mutational-target offset) -------------------------
    # mu_g = mu_mis + mu_syn + mu_lof per gene (lof imputed from mis+syn when gnomAD has none —
    # the same rule Step 9 applies; a missing mis or syn means NO target, never 0). C is fit over
    # the FULL table, zero-count genes included: the called list is a zero-truncated sample, and
    # fitting on called genes only inflates C and hides every real excess. exp_carriers_mu = C*mu_g
    # is what a gene of this mutational size is expected to collect; the Poisson tail on
    # n_carriers against it is `p_carrier_excess`, the rank the case-only p cannot provide
    # (that one is monotone in carrier count, i.e. in gene size). It is a size-normalised RANK,
    # not an association test, and Step 9's artifact panel is what separates technical excess
    # from biology.
    def _mu_tot(mrow):
        mis = _num(mrow.get(mt_mis_c)) if mt_mis_c else None
        syn = _num(mrow.get(mt_syn_c)) if mt_syn_c else None
        lof = _num(mrow.get(mt_lof_c)) if mt_lof_c else None
        if mis is None or syn is None:
            return None, "none"
        if lof is None:
            return (mis + syn) * (1.0 + impute), "imputed"
        return mis + syn + lof, "gnomad"

    mu_of = {}
    for _g, _row in mtg.items():
        _m, _src = _mu_tot(_row)
        if _m is not None and _m > 0:
            mu_of[_g] = _m
    C_mu = None
    if mu_of and n_trios > 0:
        _sum_mu = sum(mu_of.values())
        _sum_n = sum(len(g["dom"] | g["bi"] | g["x"]) for gene, g in genes.items() if gene in mu_of)
        C_mu = (_sum_n / _sum_mu) if (_sum_mu > 0 and _sum_n > 0) else None
    if (mtg or mut) and not mu_of:
        sys.stderr.write("WARN: no usable mu_mis+mu_syn(+mu_lof) columns in --mutational-target/"
                         "--mutrate; the size-normalised recurrence rank is unavailable — recurrent "
                         "genes are ordered by the case-only p (a gene-size ranking). Pass the gnomAD "
                         "v2.1.1 constraint table (prioritization.resources.mutational_target).\n")
    elif not (mtg or mut):
        sys.stderr.write("WARN: no --mutational-target (nor --mutrate): the size-normalised "
                         "recurrence rank is unavailable; recurrent genes are ordered by the "
                         "case-only p, which is a gene-size ranking.\n")

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
        def _recur(n, faf_map, prob, n_total=None):
            n_total = n_trios if n_total is None else n_total
            if binom is None or n_total <= 0 or n < min_carriers or not faf_map:
                return None, None
            p = prob(list(faf_map.values()))
            if not p or p <= 0:
                return None, None
            return n_total * p, float(binom.sf(n - 1, n_total, p))

        exp_car, p_recurrence = _recur(len(g["dom"]), g["dom_faf"],
                                       lambda f: p_carrier_hwe(f, absent_floor, 2))
        _, p_rec_bi = _recur(len(g["bi"]), g["bi_faf"],
                             lambda f: p_biallelic_hwe(f, absent_floor))
        _, p_rec_x = _recur(len(g["x"]), g["x_faf"],
                            lambda f: p_carrier_hwe(f, absent_floor, 1), n_total=n_x)

        # size-normalised rank (see above); None when the gene has no mutational target
        mu_g = mu_of.get(gene)
        exp_mu = (C_mu * mu_g) if (C_mu and mu_g) else None
        excess_mu = (n_carriers / exp_mu) if (exp_mu and exp_mu > 0) else None
        p_excess = (float(poisson.sf(n_carriers - 1, exp_mu))
                    if (poisson is not None and exp_mu and exp_mu > 0 and n_carriers > 0) else None)

        # optional SECONDARY de novo Poisson enrichment
        p_enrich = exp = None
        dn_mu_src = "none"
        if can_enrich and gene in mut:
            mrow = mut[gene]
            # A MISSING rate is never 0.0 (the rule everywhere in this repo): with no missense
            # rate there is no expectation and no test; a missing pLoF rate is imputed from
            # mis+syn exactly as Step 9 does, and the provenance rides out as dn_mu_src.
            mu_lof = _num(mrow.get(mut_lof_c)) if mut_lof_c else None
            mu_mis = _num(mrow.get(mut_mis_c)) if mut_mis_c else None
            mu_syn = _num(mrow.get(mut_syn_c)) if mut_syn_c else None
            mu = None
            if mu_mis is not None:
                if mu_lof is not None:
                    mu, dn_mu_src = mu_lof + mu_mis, "gnomad"
                elif mu_syn is not None:
                    mu, dn_mu_src = mu_mis + (mu_mis + mu_syn) * impute, "imputed"
            if mu is not None:
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
        # DISTINCT-variant recurrence (independent hits -> gene signal) vs SAME-variant recurrence
        # (one shared variant across carriers -> founder OR artifact; golden rule #2 treats internal
        # recurrence as artifact/blocklist-suspect, so it is ranked BELOW distinct-variant, not headlined).
        # Compared as per-trio SETS of variant keys, so two trios sharing one compound-het PAIR
        # read same_variant (each leg used to count as a distinct variant).
        rec_kind = ""
        for cset, per_trio in ((g["dom"], g["dom_sets"]), (g["bi"], g["bi_sets"]),
                               (g["x"], g["x_sets"])):
            if len(cset) >= min_carriers:
                distinct = {frozenset(v) for v in per_trio.values()}
                rec_kind = "same_variant" if len(distinct) <= 1 else "distinct_variant"
                break
        # rank by the STRONGEST recurrence signal across the applicable models (so recessive/X-only
        # recurrent genes are ordered by their own p, not left at 1.0)
        best_p = min([p for p in (p_recurrence, p_rec_bi, p_rec_x) if p is not None], default=None)
        rank_basis = "mu_normalised" if (rank_by_mu and p_excess is not None) else "case_only"
        rows.append({
            "gene": gene, "n_carriers": n_carriers, "n_dominant": len(g["dom"]),
            "n_biallelic": len(g["bi"]), "n_xlinked": len(g["x"]), "n_denovo": len(g["dn"]),
            "recurrent": "1" if n_carriers >= min_carriers else "0", "recurrence_kind": rec_kind,
            "exp_carriers": exp_car, "p_recurrence": p_recurrence, "best_p": best_p,
            "p_recurrence_biallelic": p_rec_bi, "p_recurrence_xlinked": p_rec_x,
            "mu_tot": mu_g, "exp_carriers_mu": exp_mu, "carrier_excess_ratio": excess_mu,
            "p_carrier_excess": p_excess, "rank_basis": rank_basis,
            "loeuf": loeuf, "pli": pli, "s_het": shet, "phaplo": phaplo,
            "constrained": "1" if constrained else "0",
            "dn_exp": (f"{exp:.4g}" if exp is not None else ""), "dn_p_enrich": p_enrich,
            "dn_mu_src": dn_mu_src,
            "modes": ";".join(modes),
        })

    # BH-FDR + exome-wide flag PER model family. The three INHERITED families correct over the
    # CALLED genes — legitimately conditional on observing >= min_carriers carriers (the test only
    # exists once a gene is nominated; see gene_burden.md) — so recessive-only recurrent genes get a
    # corrected q.
    for pcol, qcol, sigcol in (
        ("p_recurrence", "q_recurrence", "recurrence_exome_wide_sig"),
        ("p_recurrence_biallelic", "q_recurrence_biallelic", "recurrence_biallelic_exome_wide_sig"),
        ("p_recurrence_xlinked", "q_recurrence_xlinked", "recurrence_xlinked_exome_wide_sig"),
    ):
        for r, q in zip(rows, bh_fdr([r[pcol] for r in rows])):
            r[qcol] = q
            r[sigcol] = "1" if (r[pcol] is not None and r[pcol] < exome_p) else "0"

    # The secondary DE NOVO family is DIFFERENT: the denovolyzeR/Samocha model tests EVERY gene in the
    # mutation table, where a zero-DNM gene is a valid null at p = poisson.sf(-1, exp) = 1.0 — NOT an
    # absent test. Correcting over called genes only would make dn_q_enrich anti-conservative by
    # ~n_model / n_called (can be ~100x) and flag a single-DNM gene as FDR-significant when the true
    # q ≈ 1.0. So BH over the full model universe: pad the p-vector with p=1.0 for every mutation-model
    # gene that produced no candidate row, then keep the q's for the real rows. (bh_fdr's ranking is
    # driven by len(present) — pass the padded VECTOR, not just a larger m.) dn_exome_wide_sig is a
    # fixed p<exome_p Bonferroni flag and is m-independent, so it is unchanged.
    dn_ps = [r["dn_p_enrich"] for r in rows]
    n_tested = sum(1 for p in dn_ps if p is not None)
    pad = max(0, len(mut) - n_tested) if mut else 0
    dn_q_full = iter(bh_fdr([p for p in dn_ps if p is not None] + [1.0] * pad))
    for r in rows:
        r["dn_q_enrich"] = next(dn_q_full) if r["dn_p_enrich"] is not None else None
        r["dn_exome_wide_sig"] = "1" if (r["dn_p_enrich"] is not None and r["dn_p_enrich"] < exome_p) else "0"

    # Rank: recurrent first; DISTINCT-variant recurrence above SAME-variant (founder/artifact); then
    # by the strongest recurrence p across models; then constraint, counts, secondary de novo.
    # The rank p: the SIZE-NORMALISED p_carrier_excess where a mutational target exists (so long
    # genes stop leading by size), else the case-only best_p. `rank_basis` says which per gene.
    def rank_key(r):
        con_key = (r["constrained"] != "1") if weight_by_constraint else 0
        if r["rank_basis"] == "mu_normalised":
            prec = r["p_carrier_excess"]
        else:
            prec = r["best_p"] if r["best_p"] is not None else 1.0
        return (r["recurrent"] != "1", r["recurrence_kind"] == "same_variant", prec, con_key,
                -r["n_carriers"], -r["n_dominant"], -r["n_biallelic"],
                r["dn_p_enrich"] if r["dn_p_enrich"] is not None else 1.0)
    rows.sort(key=rank_key)

    out_cols = ["gene", "n_carriers", "n_dominant", "n_biallelic", "n_xlinked", "n_denovo",
                "recurrent", "recurrence_kind", "exp_carriers", "p_recurrence", "q_recurrence",
                "recurrence_exome_wide_sig",
                "p_recurrence_biallelic", "q_recurrence_biallelic", "recurrence_biallelic_exome_wide_sig",
                "p_recurrence_xlinked", "q_recurrence_xlinked", "recurrence_xlinked_exome_wide_sig",
                "mu_tot", "exp_carriers_mu", "carrier_excess_ratio", "p_carrier_excess", "rank_basis",
                "loeuf", "pli", "s_het", "phaplo", "constrained",
                "dn_exp", "dn_p_enrich", "dn_mu_src", "dn_q_enrich", "dn_exome_wide_sig", "modes"]
    with open(args.out, "w") as out:
        out.write("\t".join(out_cols) + "\n")
        for r in rows:
            out.write("\t".join(_fmt(r.get(c)) for c in out_cols) + "\n")

    n_recurrent = sum(1 for r in rows if r["recurrent"] == "1")
    n_rec_con = sum(1 for r in rows if r["recurrent"] == "1" and r["constrained"] == "1")
    # Count a gene as recurrence-significant if ANY inherited family (dominant / biallelic / X-linked)
    # is significant — matching how genes_recurrent uses the all-model carrier universe. Reading only
    # the dominant family (as before) under-counted recessive-only / X-linked-only significant genes;
    # genes.ranked.tsv already carries every per-family column, so this only fixes the summary tally.
    _sig_cols = ("recurrence_exome_wide_sig", "recurrence_biallelic_exome_wide_sig",
                 "recurrence_xlinked_exome_wide_sig")
    _q_cols = ("q_recurrence", "q_recurrence_biallelic", "q_recurrence_xlinked")
    n_rec_sig = sum(1 for r in rows if any(r.get(c) == "1" for c in _sig_cols))
    n_rec_fdr = sum(1 for r in rows
                    if any(r.get(c) is not None and r[c] < fdr_q for c in _q_cols))
    audit.record("06_burden", "calls_in", n_calls_in)
    audit.record("06_burden", "calls_no_gene", n_no_gene)
    if n_no_gene:
        sys.stderr.write(f"WARN: {n_no_gene} of {n_calls_in} calls carry no gene symbol or ID and "
                         "cannot be aggregated per gene; they survive in candidates.calls.tsv and "
                         "the variant layer but contribute to no gene row.\n")
    audit.record("06_burden", "n_trios", n_trios)
    audit.record("06_burden", "n_male_trios", n_male)
    audit.record("06_burden", "genes_rank_mu_normalised",
                 sum(1 for r in rows if r["rank_basis"] == "mu_normalised"))
    if C_mu is not None:
        audit.record("06_burden", "mu_scaling_C", f"{C_mu:.6g}")
    audit.record("06_burden", "genes_nominated", len(rows))
    audit.record("06_burden", "genes_recurrent", n_recurrent)
    audit.record("06_burden", "genes_recurrent_constrained", n_rec_con)
    audit.record("06_burden", "genes_recurrence_exome_wide_sig", n_rec_sig)
    audit.record("06_burden", "genes_recurrence_fdr_sig", n_rec_fdr)
    sys.stderr.write(
        f"Step 6 complete: {len(rows)} genes, {n_trios} trios -> {args.out}\n"
        f"  recurrent (>= {min_carriers} carriers): {n_recurrent}; recurrent+constrained: {n_rec_con}\n"
        f"  recurrence exome-wide sig (p<{exome_p:g}): {n_rec_sig}; FDR q<{fdr_q}: {n_rec_fdr}\n"
        f"  NOTE: p_recurrence is a case-only RANK (it saturates on private variants), not a "
        f"calibrated test; {sum(1 for r in rows if r['rank_basis'] == 'mu_normalised')} genes are "
        f"ordered by the size-normalised p_carrier_excess"
        + ("" if C_mu is not None else " (none — no mutational target available)") + "\n"
    )
    if args.mutrate:
        sys.stderr.write(
            "  NOTE: the de novo Poisson enrichment (secondary) is UNCALIBRATED here — the "
            "expectation is not yet scaled to an observed synonymous de novo rate "
            "(--syn-denovo-count is reserved); de novo review is handled by separate machinery.\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
