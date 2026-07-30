#!/usr/bin/env python3
"""Pipeline Step 9: prioritize the candidate calls (gene excess + artifact panel + variant tier).

Step 6 nominates GENES by recurrence across individuals. It cannot say whether a gene is
emitting more candidate rows than its mutational target predicts, and nothing downstream
re-ranks the individual CALLS — so a run leaves the reviewer a flat list. This step supplies
both:

  * **Gene layer.** ``E_g = C * mu_g`` against a gnomAD v2.1.1 per-gene Samocha target, with a
    negative-binomial null whose ``(C, alpha)`` are re-fit per cohort on an iteratively trimmed
    bulk (so the artifact tail cannot calibrate its own null) -> ``p_nb`` / BH ``q_nb``, a
    Poisson comparison, a mid-p calibration diagnostic, a six-signal artifact panel collapsed
    to an integer corroboration count, and a four-tier graded down-weight.
  * **Variant layer.** A V0-V5 screening tier under ACMG/ClinGen SVI mechanism gating, plus an
    additive Tavtigian-style ``priority_points`` composite in which EVERY term is a separate
    reported column, so a reviewer can read exactly why a call ranked where it did.

**Never-drop is a hard invariant here.** No row is ever removed: a down-weight sets a tier, a
separately-reported additive penalty, and a human-readable reason. Row-count conservation is
asserted in code before the output is written. The maximum gene penalty (-3.0) is deliberately
too small to demote a variant carrying strong molecular evidence in a constrained gene — it is
a re-rank, not a veto.

**Two rankings are always emitted**: ``rank_agnostic`` (no gene-list prior of any kind) and
``rank_prior`` (with the optional ``--gene-prior`` overlay, which defaults OFF so hprv stays
phenotype-agnostic). ``rank_delta`` exposes exactly which calls were promoted by list
membership — the set to scrutinise for confirmation bias, and the set that would be invisible
under a single blended ranking. With no overlay the two rankings are identical.

Optional resources degrade with a loud WARN, exactly as Step 6 does for ``--constraint``:
without ``--mutrate`` there is no gene excess statistic (every gene reads ``T0``, and the
variant layer still runs); without ``--established-genes`` the control exemption is empty, and
the run WARNs loudly and continues (the ceiling simply never fires). It is only when the flag IS
supplied and the file resolves to fewer than ``--min-control-genes`` entries that the run HALTS —
that case means a mis-specified or truncated control file, which would silently down-weight
established genes.

See docs/prioritization.md.

Usage:
  09_prioritize.py --variants igv/variants.tsv --genes genes.ranked.tsv \\
      --mutrate constraint.by_gene.tsv --config cfg.yaml \\
      --out-variants variants.prioritized.tsv --out-genes genes.prioritized.tsv \\
      [--n-trios N] [--established-genes ctrl.txt] [--gene-prior overlay.tsv] [--force]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import hashlib
import os
import re
import sys

from hprv import audit
from hprv import prioritize as P
from hprv.config import get, load_config

GENE_KEYS = ("gene", "gene_symbol", "symbol")

# Gene-level columns emitted for every gene in the candidate list.
GENE_COLUMNS = [
    "gene", "gene_id", "n_observed", "n_sites", "n_trios_obs", "recurrence_shape",
    "max_site_share", "per_trio",
    "mu_mis", "mu_syn", "mu_lof", "mu_tot", "mu_lof_src", "E_expected", "E_source",
    "excess_ratio", "p_nb", "q_nb", "p_pois", "q_pois", "used_in_null_fit",
    "oe_syn", "oe_lof_upper", "pLI", "s_het", "classic_caf", "constraint_flag",
    "cds_length", "segdup98_frac",
    "sig_constraint_flag", "sig_oe_syn", "sig_caf_low", "sig_segdup", "sig_family",
    "sig_saturation", "corroboration_count",
    "gene_tier", "gene_artifact_penalty", "established_gene_control",
    "control_ceiling_applied", "cds_ceiling_applied", "review_flag", "downweight_reason",
    "gene_list_prior_member",
    "n_carriers", "n_dominant", "n_biallelic", "n_xlinked", "n_denovo", "recurrent",
    "recurrence_kind", "p_recurrence", "q_recurrence",
]

# Per-variant columns. Identity first, then EVERY scoring term as its own column (the hard
# transparency requirement), then the totals/ranks, then the evidence behind each term.
VARIANT_COLUMNS = [
    "chrom", "pos", "ref", "alt", "trio_id", "gene", "inheritance", "origin", "pair_id",
    "pts_molecular", "pts_rarity", "pts_gene_constraint", "pts_recurrence", "pts_quality",
    "pts_clinical", "pts_moi", "pts_gene_artifact", "pts_gene_list_prior",
    "priority_points_agnostic", "rank_agnostic",
    "priority_points_prior", "rank_prior", "rank_delta", "cap_applied",
    "variant_tier", "variant_tier_reason", "molecular_effect_class",
    "nmd_status", "plof_confidence", "spliceai_status", "spliceai_ds",
    "missense_evidence_source", "cadd", "consequence", "impact",
    "rarity_strength", "grpmax_af", "max_af", "max_af_pops", "rarity_driven_by_single_group",
    "constraint_gate",
    "gt_qc_pass", "gt_qc_fail_reason", "partner_leg_quality_unknown",
    "child_gt", "child_GQ", "child_DP", "child_AB",
    "nhf_status", "nhf_max_fraction", "nhf_max_reads",
    "moi_coherence", "moi_caveat", "clinvar_strength", "clinvar_review_status", "clin_sig",
    "gene_tier", "gene_artifact_penalty", "excess_ratio", "E_expected", "E_source",
    "n_observed", "q_nb", "per_trio", "corroboration_count",
    "sig_constraint_flag", "sig_oe_syn", "sig_caf_low", "sig_segdup", "sig_family",
    "sig_saturation", "downweight_reason", "established_gene_control",
    "control_ceiling_applied", "review_flag",
    "pLI", "oe_lof_upper", "oe_syn", "segdup98_frac",
    "n_carriers", "n_dominant", "n_biallelic", "n_xlinked", "same_variant_recurrence",
    "p_recurrence", "q_recurrence", "gene_list_prior_member",
    # Overlay provenance, so a reviewer can see WHICH prior applied and why. `..._weight` is
    # UNCALIBRATED (an ordering default, never a likelihood ratio); `..._excluded_non_germline`
    # marks a row whose evidence class is not germline predisposition evidence (somatic drivers)
    # and which therefore contributed ZERO despite carrying a weight; `..._set_applied` names the
    # gene set when a SET-level prior beat the gene-level one (combined by MAX, never sum).
    "gene_list_prior_tier", "gene_list_prior_weight", "gene_list_prior_evidence_class",
    "gene_list_prior_excluded_non_germline", "gene_list_prior_set_applied",
]


def _fmt(x):
    """TSV cell. None -> '' (MISSING, never 0.0); bool -> 1/0; float -> 6 significant digits."""
    if x is None:
        return ""
    if isinstance(x, bool):
        return "1" if x else "0"
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def _run_key(args) -> str:
    """Content key for the idempotency marker: the input variants, the config, and every resolved
    resource input. Anything that can change the output must be in here, or a re-run silently
    serves a stale ranking (see the caller). Paths are included alongside their bytes so that
    swapping WHICH table is supplied also invalidates."""
    h = hashlib.sha256()
    for path in (args.variants, args.config, args.genes, args.mutrate, args.constraint,
                 args.segdup, args.established_genes, args.gene_moi, args.gene_prior):
        h.update(b"\x00" + str(path or "").encode())
        if path and os.path.exists(path):
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
    # scalar args that change the arithmetic
    for v in (args.n_trios, args.min_control_genes):
        h.update(b"\x00" + str(v).encode())
    return h.hexdigest()


def _open_text(path):
    """Open a TSV that may be bgzipped. Same contract as scripts/join_constraint.py:_open.

    Load-bearing: the prepared mutational-target resource IS bgzipped
    (`constraint/mutational_target.by_gene.txt.bgz` — it is the gnomAD v2.1.1 constraint table
    kept unjoined, and gnomAD ships it `.bgz`). bgzip is gzip-compatible, so stdlib gzip reads
    it. Without this branch a real run would hand `open()` a binary file and either crash or
    silently lose the excess statistic; the integration mock uses a plain `.tsv`, so the
    compressed path is exercised only here and by `test_prioritize_reads_bgzipped_tables`.
    """
    if path.endswith((".gz", ".bgz")):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path)


def _read_tsv(path):
    with _open_text(path) as fh:
        sniff = fh.readline()
        delim = "\t" if "\t" in sniff else ","
        fh.seek(0)
        return list(csv.DictReader(fh, delimiter=delim))


def _keyed_by_gene(path, label):
    """-> ({gene: row}, columns). Missing path -> ({}, []). A file with no gene key WARNs."""
    if not path:
        return {}, []
    if not os.path.exists(path):
        sys.stderr.write(f"WARN: {label} not found at {path}; ignoring\n")
        return {}, []
    rows = _read_tsv(path)
    cols = list(rows[0].keys()) if rows else []
    keycol = next((c for c in cols if c.lower().lstrip("#").strip() in GENE_KEYS), None)
    if keycol is None:
        sys.stderr.write(f"WARN: no gene key column in {label} ({path}; have {cols}); ignoring\n")
        return {}, cols
    out = {}
    for r in rows:
        g = (r.get(keycol) or "").strip()
        if g and g.lower() not in GENE_KEYS:      # guard the literal header artifact
            out.setdefault(g, r)                  # first wins
    return out, cols


def _find(cols, *names):
    low = {c.lower().lstrip("#").strip(): c for c in cols}
    for n in names:
        if n in low:
            return low[n]
    return None


def _read_gene_prior_overlay(path, cfg):
    """Read the optional Class-B overlay -> ``{GENE: entry}`` via ``prioritize.parse_gene_prior_overlay``.

    Accepts a bare symbol list OR a full TSV (the GCT resource's 12-column schema: gene, tier,
    prior_weight, evidence_class, site_specificity, moi, pmids, curated_cpg_standing, replication,
    gene_sets, anatomical_transfer, notes). **Unknown columns pass through untouched** so an
    overlay can carry its own provenance without this reader knowing the schema.

    If a sibling ``.json`` carries a ``gene_sets`` block, it is loaded too — set-level priors are
    combined with gene-level ones by MAX, never SUM (both can derive from the same study; the GCT
    resource's ``FA_HR_PATHWAY_23`` and its per-gene FA rows both come from PMID 40906985).

    **The join is on gene symbol only — never on MOI.** See ``parse_gene_prior_overlay``.
    """
    if not path:
        return {}
    if not os.path.exists(path):
        sys.stderr.write(f"WARN: --gene-prior overlay not found at {path}; ignoring\n")
        return {}
    rows = []
    with _open_text(path) as fh:
        head = fh.readline()
    if "\t" in head and any(k in head.lower() for k in ("gene", "symbol")):
        rows = _read_tsv(path)
    else:
        rows = [ln.strip() for ln in _open_text(path)
                if ln.strip() and not ln.startswith("#")]
    gene_sets = {}
    sidecar = re.sub(r"\.(tsv|txt|csv)(\.b?gz)?$", "", path) + ".json"
    if os.path.exists(sidecar):
        try:
            with _open_text(sidecar) as fh:
                gene_sets = (json.load(fh) or {}).get("gene_sets") or {}
        except (ValueError, OSError) as e:
            sys.stderr.write(f"WARN: could not read gene_sets from {sidecar} ({e}); "
                             "gene-level priors only\n")
    overlay = P.parse_gene_prior_overlay(rows, gene_sets=gene_sets, cfg=cfg)
    n_excl = sum(1 for e in overlay.values() if e.get("germline_excluded"))
    n_set = sum(1 for e in overlay.values() if e.get("source") == "gene_set")
    sys.stderr.write(
        f"--gene-prior overlay: {len(overlay)} genes from {os.path.basename(path)}"
        + (f" (+{len(gene_sets)} gene set(s) from the JSON sidecar, {n_set} members contributed "
           f"by set membership alone; set and gene priors combined by MAX, never sum)"
           if gene_sets else "")
        + (f"; {n_excl} row(s) carry a NON-GERMLINE evidence class (e.g. somatic drivers) and "
           "contribute ZERO germline prior — reported, not silently dropped" if n_excl else "")
        + "\n")
    if any(e.get("prior_weight_raw") is not None for e in overlay.values()):
        sys.stderr.write(
            "  NOTE: prior_weight is UNCALIBRATED — an ordering default, NOT a likelihood ratio. "
            "It scales the configured maximum prior; never report it as an odds ratio.\n")
    return overlay


def _read_gene_set(path, label):
    """One gene symbol per line (or the first column of a TSV); '#' comments skipped."""
    if not path:
        return set()
    if not os.path.exists(path):
        sys.stderr.write(f"WARN: {label} not found at {path}; ignoring\n")
        return set()
    genes = set()
    with _open_text(path) as fh:      # a control union / overlay may equally ship gzipped
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tok = line.split("\t")[0].split(",")[0].strip()
            if tok and tok.lower() not in GENE_KEYS:
                genes.add(tok)
    return genes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True,
                    help="igv/variants.tsv (Step 8) or candidates.calls.tsv (Step 5)")
    ap.add_argument("--genes", default="", help="genes.ranked.tsv from Step 6 (recurrence terms)")
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-variants", required=True)
    ap.add_argument("--out-genes", required=True)
    ap.add_argument("--mutrate", default="",
                    help="per-gene mutational-target table (gnomAD v2.1.1 constraint: "
                         "mu_mis/mu_syn/mu_lof + oe_syn/pLI/classic_caf/constraint_flag). "
                         "Without it there is no excess statistic (WARN, every gene T0).")
    ap.add_argument("--constraint", default="",
                    help="optional per-gene constraint TSV (oe_lof_upper/pli/s_het/phaplo) — "
                         "used when --mutrate does not already carry those columns")
    ap.add_argument("--segdup", default="",
                    help="optional gene<TAB>segdup98_frac table (build must MATCH the "
                         "mutational-target table's coordinates — see docs/prioritization.md)")
    ap.add_argument("--established-genes", default="",
                    help="version-pinned phenotype-AGNOSTIC established-gene-validity union "
                         "(one symbol per line). Enables the auditable control ceiling.")
    ap.add_argument("--gene-moi", default="",
                    help="optional gene<TAB>moi table (AD/AR/XL...) for MOI coherence; "
                         "absent => moi_unknown, which is EXACTLY neutral")
    ap.add_argument("--gene-prior", default="",
                    help="OPTIONAL Class-B phenotype overlay (one symbol per line). OFF by "
                         "default: with no overlay rank_prior == rank_agnostic.")
    ap.add_argument("--n-trios", type=int, default=0,
                    help="SCREENED trio count (never the number with a call) for per_trio")
    ap.add_argument("--min-control-genes", type=int, default=-1,
                    help="halt if the established-gene union has fewer than this many genes "
                         "(default from prioritization.gene_downweight.min_control_genes)")
    ap.add_argument("--force", action="store_true", help="ignore the .done marker and re-run")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    if not os.path.exists(args.variants):
        sys.stderr.write(f"ERROR: --variants not found: {args.variants}\n")
        return 1

    # Idempotency: CONTENT-keyed, not mtime-keyed. Every threshold in this step comes from the
    # config, and the excess statistic depends on which optional resources were supplied, so the
    # key covers the input variants, the config, AND the resolved resource inputs. An mtime test
    # against --variants alone would report "cached" after a threshold change or after a
    # --mutrate table first became available, serving a stale ranking as though it were fresh —
    # the same silent-staleness class the Step-1/2/4/8b caches were content-keyed to avoid.
    marker = args.out_variants + ".done"
    run_key = _run_key(args)   # NB distinct name: `key` is reused as a loop variable below
    if not args.force and os.path.exists(marker) \
            and os.path.getsize(args.out_variants or os.devnull) > 0 \
            and os.path.exists(args.out_genes) and os.path.getsize(args.out_genes) > 0:
        try:
            with open(marker) as fh:
                cached = fh.read().strip()
        except OSError:
            cached = ""
        if cached == run_key:
            sys.stderr.write(f"Step 9: cached — inputs, config and resources unchanged; skipping "
                             f"(rm {marker} or pass --force to re-run)\n")
            return 0
        if cached:
            sys.stderr.write("Step 9: inputs/config/resources CHANGED since the cached run — "
                             "re-prioritizing.\n")

    variants = _read_tsv(args.variants)
    n_in = len(variants)
    if n_in == 0:
        sys.stderr.write(f"ERROR: {args.variants} has no data rows — refusing to write an empty "
                         "prioritization (that would look like a successful run that found "
                         "nothing). Check Step 5/8 output.\n")
        return 1

    # variants.tsv (Step 8) names the gene column `gene`; candidates.calls.tsv (Step 5) uses
    # `symbol`, its mode column is `mode` not `inheritance`, and its QC columns are lowercase.
    # Accept BOTH so the step runs at either point in the pipeline.
    vcols = list(variants[0].keys())
    gene_c = _find(vcols, "gene", "symbol")
    if gene_c is None:
        sys.stderr.write(f"ERROR: {args.variants} has no gene/symbol column (have {vcols})\n")
        return 1
    mode_c = _find(vcols, "inheritance", "mode")
    alias = {"child_GQ": _find(vcols, "child_gq"), "child_DP": _find(vcols, "child_dp"),
             "child_AB": _find(vcols, "child_ab"), "clin_sig": _find(vcols, "clin_sig", "clnsig"),
             "grpmax_af": _find(vcols, "grpmax_af", "frequency")}

    def vget(row, canonical):
        """Read a canonical column name through the Step-5/Step-8 alias map."""
        if canonical in row and row.get(canonical) not in (None, ""):
            return row.get(canonical)
        src = alias.get(canonical)
        return row.get(src) if src else None

    # --- thresholds (all config defaults; nothing hardcoded) ---
    pfx = "prioritization"
    null_model = str(get(cfg, f"{pfx}.excess.null_model", "negative_binomial"))
    trim_p = float(get(cfg, f"{pfx}.excess.trim_p", 1.0e-3))
    trim_max = float(get(cfg, f"{pfx}.excess.trim_max_fraction", 0.05))
    q_thresh = float(get(cfg, f"{pfx}.excess.q_threshold", 0.05))
    alpha_method = str(get(cfg, f"{pfx}.excess.alpha_method", "mle"))
    impute = float(get(cfg, f"{pfx}.excess.offset.mu_lof_impute_factor", 0.0516))
    cds_b0 = float(get(cfg, f"{pfx}.excess.offset.cds_fallback.b0", -18.4484))
    cds_b1 = float(get(cfg, f"{pfx}.excess.offset.cds_fallback.b1", 1.0570))
    caf_pct = float(get(cfg, f"{pfx}.signals.caf_low.percentile", 0.10))
    max_dw_frac = float(get(cfg, f"{pfx}.gene_downweight.max_downweight_fraction", 0.20))
    min_ctrl = args.min_control_genes if args.min_control_genes >= 0 else \
        int(get(cfg, f"{pfx}.gene_downweight.min_control_genes", 1000))
    prior_enabled = bool(get(cfg, f"{pfx}.composite.gene_list_prior.enabled", False))
    # The Class-B overlay is OFF unless something says otherwise, and `enabled` is that
    # something. An explicit --gene-prior IS the intent, so it always wins; a config PATH is
    # honored only when `enabled: true`. Falling back to the config path unconditionally would
    # mean a stale or leftover path silently starts promoting genes on a run the operator
    # believes is phenotype-agnostic — which is precisely the contract this module exists to
    # keep. (Caught by the integration assertion `overlay OFF: rank_prior == rank_agnostic`.)
    prior_path = args.gene_prior or (
        (get(cfg, f"{pfx}.composite.gene_list_prior.path", "") or "") if prior_enabled else "")

    # --- observed counts per gene from the candidate rows ---
    counts, sites, trios_obs, site_trios = {}, {}, {}, {}
    for r in variants:
        g = (r.get(gene_c) or "").strip()
        if not g:
            continue
        counts[g] = counts.get(g, 0) + 1
        key = f"{r.get('chrom')}:{r.get('pos')}:{r.get('ref')}:{r.get('alt')}"
        sites.setdefault(g, set()).add(key)
        t = (r.get("trio_id") or "").strip()
        if t:
            trios_obs.setdefault(g, set()).add(t)
            site_trios.setdefault(g, {}).setdefault(key, set()).add(t)
    n_genes_in = len(counts)

    n_trios = args.n_trios or 0
    if not n_trios:
        sys.stderr.write("WARN: --n-trios not provided; per_trio and the cohort-saturation "
                         "signal are SKIPPED (saturation is the strongest signal in the panel, "
                         "195x enriched). Pass the SCREENED trio count.\n")

    # --- the mutational-target table (the one new resource) ---
    mut, mcols = _keyed_by_gene(args.mutrate, "--mutrate mutational-target table")
    con, ccols = _keyed_by_gene(args.constraint, "--constraint table")
    seg, scols = _keyed_by_gene(args.segdup, "--segdup table")
    moi, moicols = _keyed_by_gene(args.gene_moi, "--gene-moi table")
    if not mut:
        sys.stderr.write(
            "WARN: no per-gene mutational-target table (--mutrate) — the gene EXCESS statistic "
            "is UNAVAILABLE for this run. Every gene reads E_expected='', gene_tier=T0 and "
            "gene_artifact_penalty=0; the variant tier and composite still run, minus the "
            "artifact term. Fetch it with scripts/prepare_resources.sh --only mutational_target "
            "(gnomAD v2.1.1 constraint, ~3 MB). See docs/prioritization.md.\n")

    c_mis = _find(mcols, "mu_mis", "mut_mis")
    c_syn = _find(mcols, "mu_syn", "mut_syn")
    c_lof = _find(mcols, "mu_lof", "mut_lof")
    c_gid = _find(mcols, "gene_id", "ensembl_gene_id", "gene_ensembl_id")
    c_oesyn = _find(mcols, "oe_syn", "oe_syn_upper") or _find(ccols, "oe_syn")
    c_caf = _find(mcols, "classic_caf", "caf")
    c_cflag = _find(mcols, "constraint_flag", "flag")
    c_cds = _find(mcols, "cds_length", "cds_len")
    c_segf = _find(scols, "segdup98_frac", "segdup_frac", "frac") if seg else \
        _find(mcols, "segdup98_frac")
    c_loeuf = _find(ccols, "oe_lof_upper", "loeuf", "loeuf_v2") or _find(mcols, "oe_lof_upper")
    c_pli = _find(ccols, "pli", "pli_v2") or _find(mcols, "pli")
    c_shet = _find(ccols, "s_het", "shet")
    c_moi = _find(moicols, "moi", "mode_of_inheritance", "inheritance") if moi else None
    if mut and not (c_mis or c_syn or c_lof):
        sys.stderr.write(
            f"ERROR: --mutrate {args.mutrate} carries no mu_mis/mu_syn/mu_lof column "
            f"(have {mcols}). That table is the OFFSET for the excess statistic; without those "
            "columns the statistic is undefined. Use the gnomAD v2.1.1 lof_metrics table (or a "
            "join that keeps its mu_* columns), not a LOEUF/pLI-only constraint file.\n")
        return 1

    def mrow(g):
        return mut.get(g) or {}

    def crow(g):
        return con.get(g) or {}

    # --- build the per-gene universe. C MUST be fit over the FULL table (including its
    # zero-count genes): the candidate list is a zero-truncated sample, and fitting C on
    # matched genes only inflates C and deflates every ratio — the direction that hides
    # artifact loci. ---
    universe_counts, universe_mus = [], []
    universe_genes = []
    for g in (mut.keys() if mut else counts.keys()):
        r = mrow(g)
        mu_tot, _src = P.mutational_target(r.get(c_mis) if c_mis else None,
                                           r.get(c_syn) if c_syn else None,
                                           r.get(c_lof) if c_lof else None, impute)
        if mu_tot is None:
            continue
        universe_genes.append(g)
        universe_counts.append(counts.get(g, 0))
        universe_mus.append(mu_tot)
    n_universe = len(universe_genes)
    n_zero = sum(1 for c in universe_counts if not c)

    fit = None
    if universe_genes and null_model != "none":
        try:
            fit = P.fit_excess_null(universe_counts, universe_mus, trim_p=trim_p,
                                    trim_max_fraction=trim_max, alpha_method=alpha_method)
        except ValueError as e:
            sys.stderr.write(f"ERROR: excess-null fit failed: {e}\n")
            return 1
    C = fit["C"] if fit else None
    alpha = fit["alpha"] if fit else None
    in_fit = {universe_genes[i] for i in fit["used_in_null_fit"]} if fit else set()

    # The Poisson arm gets its OWN iterative trim, not the NB's C. The trim loop belongs to
    # whichever null is being fit: evaluating the Poisson at the NB's scaling constant measures a
    # hybrid nobody would deploy, and it changes the headline number materially (2.4-2.5x
    # anti-conservative when the Poisson arm is trimmed, 1.82x when it is not). Independent
    # validation showed the published 2.51x is only reproducible with the Poisson arm trimmed, so
    # both arms are now fit the way they would actually be used.
    calib = P.calibrate_null(universe_counts, universe_mus, C, alpha, trim_p=trim_p) if fit else {}

    # --- per-gene rows over the CANDIDATE genes (the universe drives the fit; the output is
    # the candidate list, so a reviewer sees every gene that produced a call) ---
    established = _read_gene_set(args.established_genes, "--established-genes union")
    if args.established_genes and len(established) < min_ctrl:
        sys.stderr.write(
            f"ERROR: the established-gene union has {len(established)} genes, below "
            f"min_control_genes={min_ctrl}. A stale or truncated control set makes the "
            "auditable control exemption silently EMPTY, which is exactly the failure mode that "
            "would down-weight established genes without anyone noticing. Fix the file or lower "
            "prioritization.gene_downweight.min_control_genes deliberately.\n")
        return 1
    if not args.established_genes:
        sys.stderr.write(
            "WARN: no --established-genes union — the auditable established-gene CEILING is "
            "INACTIVE for this run, so a real predisposition gene with a technical excess can "
            "reach T2/T3 and be penalised. The validation cohort had 13 such genes. Supply a "
            "phenotype-agnostic gene-validity union (ClinGen Definitive/Strong + dosage HI3 + "
            "actionability) to enable it.\n")

    prior_overlay = _read_gene_prior_overlay(prior_path, cfg) if prior_path else {}
    prior_genes = set(prior_overlay)
    if args.gene_prior and not prior_enabled:
        sys.stderr.write(
            f"WARN: --gene-prior was passed explicitly ({args.gene_prior}) while "
            "prioritization.composite.gene_list_prior.enabled is false. The explicit flag wins: "
            "the overlay IS applied to rank_prior (that is what the second ranking is for) and "
            "rank_agnostic stays free of it. Set enabled: true in config to record the intent.\n")

    # classic_caf's bottom-decile cut is a COHORT property (the 10th percentile among the genes
    # in THIS candidate list), not a constant. The percentile is taken over the genes that HAVE
    # a value; genes with a null classic_caf are then read as 0 and land inside the cut. Taking
    # the percentile over the nulls-as-zeros instead would drag the cut point down by however
    # many genes gnomAD has no pLoF record for (93 of 10,377 on the validation cohort) — i.e. a
    # missing value would silently change the threshold applied to the measured ones.
    caf_vals = [P._num(mrow(g).get(c_caf)) for g in counts] if c_caf else []
    caf_cutoff = P.percentile([v for v in caf_vals if v is not None], caf_pct) if caf_vals else None
    thr = P.signal_thresholds(cfg, caf_low_cutoff=caf_cutoff)

    ranked, rcols = _keyed_by_gene(args.genes, "--genes genes.ranked.tsv")
    if args.genes and not ranked:
        sys.stderr.write("WARN: no usable rows from --genes; the recurrence term will be 0 for "
                         "every variant (Step 6 output missing or unkeyed)\n")

    gene_rows = {}
    e_source_tally = {"gnomad_mu": 0, "cds_fallback": 0, "none": 0}
    for g, n_obs in counts.items():
        r, cr = mrow(g), crow(g)
        mu_tot, mu_src = P.mutational_target(r.get(c_mis) if c_mis else None,
                                             r.get(c_syn) if c_syn else None,
                                             r.get(c_lof) if c_lof else None, impute)
        cds_len = P._num(r.get(c_cds)) if c_cds else None
        e_source = "none"
        if mu_tot is not None:
            e_source = "gnomad_mu"
        else:
            # CDS-length fallback: a +/-30% offset. Emitted with its source AND a tier ceiling,
            # never silently mixed in with the gnomAD-mu genes.
            mu_tot = P.mu_from_cds(cds_len, cds_b0, cds_b1)
            e_source = "cds_fallback" if mu_tot is not None else "none"
        e_source_tally[e_source] = e_source_tally.get(e_source, 0) + 1

        E = (C * mu_tot) if (C and mu_tot) else None
        ratio = (n_obs / E) if (E and E > 0.0) else None
        p_nb = P.nb_sf(n_obs, E, alpha) if (E and alpha) else None
        p_pois = P.pois_sf(n_obs, E) if E else None
        per_trio = (n_obs / n_trios) if n_trios else None

        n_site = len(sites.get(g, ()))
        n_trio_obs = len(trios_obs.get(g, ()))
        # Same-variant vs distinct-variant excess. These have OPPOSITE interpretations: a high
        # max_site_share with many trios is a single-site excess (founder allele or a recurrent
        # mapping artifact) and should be escalated to that SITE for read-level review, not
        # charged to the whole gene; recurrence_shape ~ 1 at high excess is the distributed
        # (mismapping / paralogue-collapse) signature the gene-level down-weight is built for.
        shape = (n_obs / n_site) if n_site else None
        max_share = None
        if n_trio_obs:
            biggest = max((len(v) for v in site_trios.get(g, {}).values()), default=0)
            max_share = biggest / n_trio_obs

        row = {
            "gene": g, "gene_id": r.get(c_gid) if c_gid else "",
            "n_observed": n_obs, "n_sites": n_site or None, "n_trios_obs": n_trio_obs or None,
            "recurrence_shape": shape, "max_site_share": max_share, "per_trio": per_trio,
            "mu_mis": P._num(r.get(c_mis)) if c_mis else None,
            "mu_syn": P._num(r.get(c_syn)) if c_syn else None,
            "mu_lof": P._num(r.get(c_lof)) if c_lof else None,
            "mu_tot": mu_tot, "mu_lof_src": mu_src,
            "E_expected": E, "E_source": e_source, "excess_ratio": ratio,
            "p_nb": p_nb, "p_pois": p_pois,
            "used_in_null_fit": (g in in_fit) if fit else None,
            "oe_syn": P._num(r.get(c_oesyn) if (c_oesyn and c_oesyn in r) else cr.get(c_oesyn)) if c_oesyn else None,
            "oe_lof_upper": P._num(cr.get(c_loeuf) if (c_loeuf and c_loeuf in cr) else r.get(c_loeuf)) if c_loeuf else None,
            "pLI": P._num(cr.get(c_pli) if (c_pli and c_pli in cr) else r.get(c_pli)) if c_pli else None,
            "s_het": P._num(cr.get(c_shet)) if c_shet else None,
            "classic_caf": P._num(r.get(c_caf)) if c_caf else None,
            "constraint_flag": (r.get(c_cflag) or "") if c_cflag else "",
            "cds_length": cds_len,
            "segdup98_frac": P._num((seg.get(g) or {}).get(c_segf) if seg else r.get(c_segf)) if c_segf else None,
            "established_gene_control": g in established,
            "gene_list_prior_member": g in prior_genes,
            "gene_moi": (moi.get(g) or {}).get(c_moi) if c_moi else "",
        }
        for c in ("n_carriers", "n_dominant", "n_biallelic", "n_xlinked", "n_denovo",
                  "recurrent", "recurrence_kind", "p_recurrence", "q_recurrence"):
            row[c] = (ranked.get(g) or {}).get(c, "")
        gene_rows[g] = row

    # BH-FDR across the FULL universe (including the zero-count genes, whose p is 1) — not
    # across the candidate genes only, which would make q anti-conservative by ~n_universe/n_called.
    if fit:
        uni_p_nb = [P.nb_sf(c, C * m, alpha) for c, m in zip(universe_counts, universe_mus)]
        uni_p_po = [P.pois_sf(c, C * m) for c, m in zip(universe_counts, universe_mus)]
        q_nb_map = dict(zip(universe_genes, P.bh_fdr(uni_p_nb)))
        q_po_map = dict(zip(universe_genes, P.bh_fdr(uni_p_po)))
        # Genes whose offset came from the CDS fallback are absent from the universe (they have
        # no gnomAD mu), so their p is BH-corrected within the candidate set they belong to.
        extra = [g for g, row in gene_rows.items() if g not in q_nb_map and row["p_nb"] is not None]
        for g, q in zip(extra, P.bh_fdr([gene_rows[g]["p_nb"] for g in extra])):
            q_nb_map[g] = q
        extra_po = [g for g, row in gene_rows.items() if g not in q_po_map and row["p_pois"] is not None]
        for g, q in zip(extra_po, P.bh_fdr([gene_rows[g]["p_pois"] for g in extra_po])):
            q_po_map[g] = q
    else:
        q_nb_map = q_po_map = {}
    for g, row in gene_rows.items():
        row["q_nb"] = q_nb_map.get(g)
        row["q_pois"] = q_po_map.get(g)

    # --- artifact panel + tier per gene ---
    tier_tally = {t: 0 for t in P.GENE_TIERS}
    tier_variants = {t: 0 for t in P.GENE_TIERS}
    signal_tally = {k: 0 for k in ("sig_constraint_flag", "sig_oe_syn", "sig_caf_low",
                                   "sig_segdup", "sig_family", "sig_saturation")}
    review_tally = {}
    n_ctrl_ceiling = n_cds_ceiling = 0
    for g, row in gene_rows.items():
        sig = P.artifact_signals(row, thr)
        row.update({k: v for k, v in sig.items() if k.startswith("sig_")})
        row["corroboration_count"] = sig["corroboration_count"]
        info = P.assign_gene_tier(row["excess_ratio"], row["q_nb"], sig["corroboration_count"],
                                  row["sig_saturation"], row["n_observed"],
                                  established_gene_control=row["established_gene_control"],
                                  e_source=row["E_source"], cfg=cfg)
        row.update(info)
        row["downweight_reason"] = P.downweight_reason(row, sig["signal_reasons"], info)
        tier_tally[info["gene_tier"]] += 1
        tier_variants[info["gene_tier"]] += row["n_observed"]
        for k in signal_tally:
            if row.get(k):
                signal_tally[k] += 1
        if info["control_ceiling_applied"]:
            n_ctrl_ceiling += 1
        if info["cds_ceiling_applied"]:
            n_cds_ceiling += 1
        if info["review_flag"]:
            review_tally[info["review_flag"]] = review_tally.get(info["review_flag"], 0) + 1

    n_dw_variants = tier_variants["T2_downweight"] + tier_variants["T3_strong_downweight"]
    dw_frac = n_dw_variants / float(n_in) if n_in else 0.0
    if dw_frac > max_dw_frac:
        sys.stderr.write(
            f"ERROR: the down-weight tiers (T2+T3) cover {n_dw_variants}/{n_in} variants "
            f"({dw_frac:.1%}), above max_downweight_fraction={max_dw_frac:.0%}. On the "
            "validation cohort this was 10.4%. A figure this high means the null is "
            "mis-specified for this cohort — most likely the mutational-target table and the "
            "candidate counts do not describe the same gene universe — not that a fifth of the "
            "exome is artifact. Fix the offset table or raise the bound deliberately.\n")
        return 1

    # --- score every variant. NEVER-DROP: one output row per input row, asserted below. ---
    out_rows = []
    tier_v_tally = {t: 0 for t in P.VARIANT_TIERS}
    cap_tally, nhf_tally = {}, {}
    for r in variants:
        g = (r.get(gene_c) or "").strip()
        grow = gene_rows.get(g, {})
        norm = dict(r)
        norm["inheritance"] = (r.get(mode_c) or "") if mode_c else ""
        for canonical in ("child_GQ", "child_DP", "child_AB", "clin_sig", "grpmax_af"):
            norm[canonical] = vget(r, canonical)
        # Pass the overlay ENTRY, not a bool: it carries the per-gene weight and the non-germline
        # exclusion, so a T1 gene and a somatic-driver row do not contribute the same points. The
        # lookup is by SYMBOL ALONE — never routed by MOI (see parse_gene_prior_overlay).
        sc = P.score_variant(norm, grow, cfg,
                             gene_prior=prior_overlay.get(P._s(grow.get("gene")).upper()))
        row = {c: norm.get(c) for c in ("chrom", "pos", "ref", "alt", "trio_id", "origin",
                                        "pair_id", "consequence", "impact", "spliceai_ds",
                                        "cadd", "max_af", "max_af_pops", "child_gt",
                                        "child_GQ", "child_DP", "child_AB", "clin_sig",
                                        "grpmax_af", "inheritance")}
        row["gene"] = g
        row.update(sc)
        for c in ("gene_tier", "gene_artifact_penalty", "excess_ratio", "E_expected", "E_source",
                  "n_observed", "q_nb", "per_trio", "corroboration_count", "downweight_reason",
                  "established_gene_control", "control_ceiling_applied", "review_flag",
                  "sig_constraint_flag", "sig_oe_syn", "sig_caf_low", "sig_segdup", "sig_family",
                  "sig_saturation", "pLI", "oe_lof_upper", "oe_syn", "segdup98_frac",
                  "n_carriers", "n_dominant", "n_biallelic", "n_xlinked",
                  "p_recurrence", "q_recurrence", "gene_list_prior_member"):
            row[c] = grow.get(c, "")
        out_rows.append(row)
        tier_v_tally[sc["variant_tier"]] = tier_v_tally.get(sc["variant_tier"], 0) + 1
        cap_tally[sc["cap_applied"]] = cap_tally.get(sc["cap_applied"], 0) + 1
        nhf_tally[sc["nhf_status"]] = nhf_tally.get(sc["nhf_status"], 0) + 1

    # THE never-drop invariant. Not a comment, an assertion: if this step ever loses a row, the
    # run must fail loudly rather than hand a reviewer a silently shortened list.
    if len(out_rows) != n_in:
        sys.stderr.write(f"ERROR: never-drop invariant VIOLATED — {n_in} input rows produced "
                         f"{len(out_rows)} output rows. Step 9 must never remove a variant.\n")
        return 1

    P.rank_rows(out_rows, "priority_points_agnostic", "rank_agnostic")
    P.rank_rows(out_rows, "priority_points_prior", "rank_prior")
    for row in out_rows:
        row["rank_delta"] = row["rank_prior"] - row["rank_agnostic"]
    n_promoted = sum(1 for r in out_rows if r["rank_delta"] < 0)

    # --- write ---
    for path, cols, rows in ((args.out_genes, GENE_COLUMNS,
                              sorted(gene_rows.values(),
                                     key=lambda r: (-(r["excess_ratio"] or 0.0), r["gene"]))),
                             (args.out_variants, VARIANT_COLUMNS,
                              sorted(out_rows, key=lambda r: r["rank_agnostic"]))):
        with open(path, "w", newline="") as out:
            out.write("\t".join(cols) + "\n")
            for r in rows:
                out.write("\t".join(_fmt(r.get(c)) for c in cols) + "\n")

    # --- audit: every funnel tally, so the run is answerable ---
    A = lambda m, v: audit.record("09_prioritize", m, v)          # noqa: E731
    A("variants_in", n_in)
    A("variants_out", len(out_rows))
    A("genes_in", n_genes_in)
    A("genes_out", len(gene_rows))
    A("n_trios", n_trios)
    A("gene_universe_size", n_universe)
    A("gene_universe_zero_count", n_zero)
    if fit:
        A("null_model", null_model)
        A("null_C", f"{C:.6g}")
        A("null_alpha", f"{alpha:.6g}")
        A("null_theta", f"{fit['theta']:.6g}")
        A("null_trim_iterations", fit["iterations"])
        A("null_genes_trimmed", fit["n_trimmed"])
        A("null_trim_fraction", f"{fit['trim_fraction']:.6g}")
        if fit.get("phi_bulk") is not None:
            A("null_phi_bulk", f"{fit['phi_bulk']:.6g}")
        if fit.get("phi_all") is not None:
            A("null_phi_all_genes", f"{fit['phi_all']:.6g}")
        for k, v in calib.items():
            A(f"calibration.{k}", f"{v:.6g}" if isinstance(v, float) else v)
        A("genes_q_nb_sig", sum(1 for r in gene_rows.values()
                                if r["q_nb"] is not None and r["q_nb"] < q_thresh))
    for k, v in e_source_tally.items():
        A(f"offset_source.{k}", v)
    for t in P.GENE_TIERS:
        A(f"gene_tier.{t}", tier_tally[t])
        A(f"gene_tier_variants.{t}", tier_variants[t])
    A("genes_downweighted_T2_T3", tier_tally["T2_downweight"] + tier_tally["T3_strong_downweight"])
    A("variants_downweighted_T2_T3", n_dw_variants)
    A("variants_downweighted_fraction", f"{dw_frac:.6g}")
    A("established_control_genes_present", sum(1 for r in gene_rows.values()
                                               if r["established_gene_control"]))
    A("control_ceiling_applied", n_ctrl_ceiling)
    A("cds_fallback_ceiling_applied", n_cds_ceiling)
    for k, v in review_tally.items():
        A(f"review_flag.{k}", v)
    for k, v in signal_tally.items():
        A(f"signal.{k}", v)
    for t, v in tier_v_tally.items():
        A(f"variant_tier.{t}", v)
    for k, v in cap_tally.items():
        A(f"cap_applied.{k}", v)
    for k, v in nhf_tally.items():
        A(f"nhf_status.{k}", v)
    A("gene_prior_overlay_genes", len(prior_genes))
    A("variants_promoted_by_prior", n_promoted)

    # Stamp the CONTENT key (not a bare touch): the next run compares inputs+config+resources
    # against this, so a threshold change or a newly-supplied resource re-prioritizes.
    with open(marker, "w") as fh:
        fh.write(run_key + "\n")

    sys.stderr.write(
        f"Step 9 complete: {len(out_rows)} variants / {len(gene_rows)} genes -> "
        f"{args.out_variants}, {args.out_genes}\n")
    if fit:
        sys.stderr.write(
            f"  null: {null_model} C={C:.6g} alpha={alpha:.4g} (theta={fit['theta']:.4g}) over "
            f"{n_universe} genes ({n_zero} with zero counts), {fit['n_trimmed']} trimmed "
            f"({fit['trim_fraction']:.2%}) in {fit['iterations']} iterations\n")
        if calib.get("n_evaluated"):
            sys.stderr.write(
                f"  calibration (mid-p on the bulk, n={calib['n_evaluated']}): "
                f"NB mean={calib['nb_mean_midp']:.4f} P(p<0.001)={calib['nb_frac_lt_0.001']:.5f}; "
                f"Poisson mean={calib['poisson_mean_midp']:.4f} "
                f"P(p<0.001)={calib['poisson_frac_lt_0.001']:.5f} (ideal 0.5 / 0.001)\n")
    sys.stderr.write(
        f"  gene tiers: " + ", ".join(f"{t}={tier_tally[t]}g/{tier_variants[t]}v"
                                      for t in P.GENE_TIERS) + "\n"
        f"  down-weighted (T2+T3): "
        f"{tier_tally['T2_downweight'] + tier_tally['T3_strong_downweight']} genes / "
        f"{n_dw_variants} variants ({dw_frac:.2%}); control ceiling applied to {n_ctrl_ceiling} "
        f"gene(s), CDS-fallback ceiling to {n_cds_ceiling}\n"
        f"  variant tiers: " + ", ".join(f"{t}={tier_v_tally[t]}" for t in P.VARIANT_TIERS) + "\n")
    sys.stderr.write(
        "  NOTE: no variant reached V5 — the NMD-escape test needs VEP's EXON / CDS_position / "
        "transcript length in variants.tsv, so every pLoF is capped at V4. This is a resource "
        "gap, not a scoring choice.\n"
        "  NOTE: priority_points is NOT an ACMG score. Do not read the totals against "
        "Tavtigian's P>=10 / LP 6-9 / VUS 0-5 bands, and never emit a P/LP/VUS label from them.\n")
    if prior_genes:
        sys.stderr.write(
            f"  overlay: {len(prior_genes)} genes; {n_promoted} variants promoted by list "
            "membership (rank_delta < 0) — that set is what to scrutinise for confirmation bias\n")
    else:
        sys.stderr.write("  overlay: none — rank_prior is identical to rank_agnostic "
                         "(hprv stays phenotype-agnostic by default)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
