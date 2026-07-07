#!/usr/bin/env python3
"""Pipeline Step 5: pedigree-aware inheritance screen with genotype QC.

Reads the per-trio candidate manifest from Step 4 (trio_id, candidates_vcf, ped),
resolves child/mother/father from each PED, and classifies each candidate variant by
inheritance mode with the refined-GQ genotype-QC gates and per-mode rarity gates:

  * de novo (autosomal + X-hemizygous), with parental-cleanliness re-verification
    and (when present) the GATK hiConfDeNovo tag;
  * homozygous recessive;
  * X-linked recessive (male hemizygous);
  * compound heterozygous in TRANS (parent-of-origin: mat + pat, or inherited + de novo).

Emits one TSV of candidate calls across all trios. See
docs/inheritance_and_genotype_qc.md and docs/pipeline_design.md (Step 5).
"""
from __future__ import annotations

import argparse
import sys

from cyvcf2 import VCF

from hprv import annotations as A
from hprv import genotype as G
from hprv.config import get, load_config
from hprv.ped import parse_ped

COLS = [
    "trio_id", "mode", "pair_id", "chrom", "pos", "ref", "alt", "gene", "symbol",
    "consequence", "impact", "faf95", "grpmax_af", "nhomalt", "revel",
    "alphamissense", "spliceai_max", "loftee", "loftee_flags", "clnsig",
    "clinvar_stars", "child_gt", "child_gq", "child_dp", "child_ab",
    "mother_gt", "father_gt", "hiConfDeNovo", "review_prior_crosscheck", "flags",
]


def fmt(x):
    return "" if x is None else (f"{x:.4g}" if isinstance(x, float) else str(x))


class Trio:
    def __init__(self, vcf: VCF, ped, thr: G.GtThresholds):
        self.thr = thr
        samples = list(vcf.samples)
        idx = {s: i for i, s in enumerate(samples)}
        for role in ("child", "father", "mother"):
            if ped[role] not in idx:
                raise KeyError(f"PED {role} {ped[role]!r} not in VCF samples {samples}")
        self.c, self.d, self.m = idx[ped["child"]], idx[ped["father"]], idx[ped["mother"]]
        self.child_male = str(ped["sex"]) == "1"
        self.has_hiconf = "ID=hiConfDeNovo" in vcf.raw_header


def base_row(trio_id, v, gt, mode, pair_id=""):
    return {
        "trio_id": trio_id, "mode": mode, "pair_id": pair_id,
        "chrom": v.CHROM, "pos": v.POS, "ref": v.REF, "alt": ",".join(v.ALT),
        "gene": A._str(v, "gene") or "", "symbol": A.symbol(v) or "",
        "consequence": A.consequence(v) or "", "impact": A.impact(v) or "",
        "faf95": fmt(A.faf95(v)), "grpmax_af": fmt(A.grpmax_af(v)),
        "nhomalt": fmt(A.nhomalt(v)), "revel": fmt(A.revel(v)),
        "alphamissense": fmt(A.alphamissense(v)), "spliceai_max": fmt(A.spliceai_max(v)),
        "loftee": A._str(v, "loftee") or "", "loftee_flags": A.loftee_flags(v) or "",
        "clnsig": A.clnsig(v) or "", "clinvar_stars": A.clinvar_stars(v),
        "child_gt": (v.gt_bases[gt.c] if v.gt_bases is not None else ""),
        "child_gq": fmt(G.gq(v, gt.c)), "child_dp": fmt(G.dp(v, gt.c)),
        "child_ab": fmt(G.allele_balance(v, gt.c)),
        "mother_gt": (v.gt_bases[gt.m] if v.gt_bases is not None else ""),
        "father_gt": (v.gt_bases[gt.d] if v.gt_bases is not None else ""),
        "hiConfDeNovo": ("1" if _has_flag(v, "hiConfDeNovo") else ""),
        "review_prior_crosscheck": "", "flags": "",
    }


def _has_flag(v, name):
    try:
        return v.INFO.get(name) is not None
    except KeyError:
        return False


def screen_trio(trio_id, vcf, gt: Trio, cfg):
    thr = gt.thr
    dom_max = float(get(cfg, "filters.rarity.dominant_max", 1e-4))
    rec_max = float(get(cfg, "filters.rarity.recessive_max", 1e-2))
    require_hiconf = bool(get(cfg, "filters.denovo.use_hiconf_tag", True))
    require_absent = bool(get(cfg, "filters.denovo.require_gnomad_absent_or_singleton", True))
    crosscheck = bool(get(cfg, "filters.denovo.crosscheck_prerefinement_pl", True))

    def rare(v, limit):
        fr = A.frequency(v)
        return fr is None or fr < limit

    rows = []
    # gene -> list of (origin, variant, key) for compound-het pairing
    comphet = {}

    for v in vcf:
        c, d, m = gt.c, gt.d, gt.m
        gc, gd, gmm = v.gt_types[c], v.gt_types[d], v.gt_types[m]
        xnonpar = G.is_x_nonpar(v)

        # --- male non-PAR X het is a QC red flag: do not call het modes there ---
        male_x = xnonpar and gt.child_male

        # ---- de novo (autosomal het, or X-hemizygous in a male) ----
        denovo_hit = False
        if not male_x and gc == G.HET and gd == G.HOM_REF and gmm == G.HOM_REF:
            denovo_hit = True
        elif male_x and gc == G.HOM_ALT and gd == G.HOM_REF and gmm == G.HOM_REF:
            denovo_hit = True
        if denovo_hit:
            child_kind = "denovo_child" if not male_x else "hom_alt"
            ok = (G.sample_qc(v, c, thr, child_kind)
                  and G.sample_qc(v, d, thr, "clean_parent")
                  and G.sample_qc(v, m, thr, "clean_parent")
                  and rare(v, dom_max))
            nh = A.nhomalt(v)
            if require_absent and nh is not None and nh > 1:
                ok = False
            if require_hiconf and gt.has_hiconf and not _has_flag(v, "hiConfDeNovo"):
                ok = False  # tag exists in this callset but not on this variant
            if ok:
                r = base_row(trio_id, v, gt, "denovo_x_hemi" if male_x else "denovo")
                if crosscheck:
                    r["review_prior_crosscheck"] = "1"
                rows.append(r)

        # ---- homozygous recessive (autosomal / X-female) ----
        if not male_x and gc == G.HOM_ALT and gd == G.HET and gmm == G.HET:
            if (G.sample_qc(v, c, thr, "hom_alt")
                    and G.sample_qc(v, d, thr, "het")
                    and G.sample_qc(v, m, thr, "het")
                    and rare(v, rec_max)):
                rows.append(base_row(trio_id, v, gt, "hom_recessive"))

        # ---- X-linked recessive (male hemizygous, carrier mother) ----
        if male_x and gc == G.HOM_ALT and gmm == G.HET and gd == G.HOM_REF:
            if (G.sample_qc(v, c, thr, "hom_alt")
                    and G.sample_qc(v, m, thr, "het")
                    and rare(v, rec_max)):
                rows.append(base_row(trio_id, v, gt, "x_linked_recessive"))

        # ---- collect compound-het candidates (het child, rare, parent-of-origin) ----
        if not male_x and gc == G.HET and G.sample_qc(v, c, thr, "het") and rare(v, rec_max):
            gene = A._str(v, "gene") or A.symbol(v)
            if gene:
                origin = None
                if gd == G.HOM_REF and gmm == G.HOM_REF:
                    origin = "denovo"
                elif gmm in (G.HET, G.HOM_ALT) and gd == G.HOM_REF:
                    origin = "mat"
                elif gd in (G.HET, G.HOM_ALT) and gmm == G.HOM_REF:
                    origin = "pat"
                if origin:
                    comphet.setdefault(gene, []).append((origin, v, f"{v.CHROM}:{v.POS}:{v.REF}:{v.ALT[0]}"))

    # ---- compound-het pairing: emit trans pairs (mat×pat, mat×denovo, pat×denovo) ----
    pair_n = 0
    for gene, cands in comphet.items():
        by = {"mat": [], "pat": [], "denovo": []}
        for origin, v, key in cands:
            by[origin].append((v, key))
        pairs = []
        for a in by["mat"]:
            for b in by["pat"] + by["denovo"]:
                pairs.append((a, b))
        for a in by["pat"]:
            for b in by["denovo"]:
                pairs.append((a, b))
        for (va, ka), (vb, kb) in pairs:
            pair_n += 1
            pid = f"{trio_id}:CH{pair_n}"
            for v in (va, vb):
                rows.append(base_row(trio_id, v, gt, "compound_het", pid))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, help="trios.candidates.tsv from Step 4")
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True, help="candidate calls TSV")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    thr = G.GtThresholds.from_config(cfg, get)

    with open(args.manifest) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]

    n_trios = 0
    all_rows = []
    for r in rows:
        trio_id, vcf_path, ped_path = r.get("trio_id"), r.get("candidates_vcf"), r.get("ped")
        ped = parse_ped(ped_path)
        if not ped:
            sys.stderr.write(f"WARN: no usable PED for {trio_id} ({ped_path!r}); skipping\n")
            continue
        vcf = VCF(vcf_path)
        try:
            gt = Trio(vcf, ped, thr)
        except KeyError as e:
            sys.stderr.write(f"WARN: {trio_id}: {e}; skipping\n")
            vcf.close()
            continue
        all_rows.extend(screen_trio(trio_id, vcf, gt, cfg))
        vcf.close()
        n_trios += 1

    with open(args.out, "w") as out:
        out.write("\t".join(COLS) + "\n")
        for r in all_rows:
            out.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")

    by_mode = {}
    for r in all_rows:
        by_mode[r["mode"]] = by_mode.get(r["mode"], 0) + 1
    sys.stderr.write(
        f"Step 5 complete: {len(all_rows)} candidate calls across {n_trios} trios "
        f"-> {args.out}\n  by mode: {by_mode}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
