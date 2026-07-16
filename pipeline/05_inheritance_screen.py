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
from hprv import audit
from hprv import genotype as G
from hprv.config import get, load_config
from hprv.ped import parse_ped

COLS = [
    "trio_id", "mode", "pair_id", "chrom", "pos", "ref", "alt", "gene", "symbol",
    # grpmax_af is THE rarity field (annotations.frequency()); max_af/max_af_pops ride
    # along for review only, so a curator can see when a call is being driven by a
    # founder-group frequency that grpmax deliberately ignores.
    "consequence", "impact", "grpmax_af", "max_af", "max_af_pops", "cadd",
    "clnsig", "child_gt", "child_gq", "child_dp", "child_ab",
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
        self.child_name = ped["child"]
        self.child_male = str(ped["sex"]) == "1"
        # sex must be positively known (1/2) to apply ploidy-aware X/Y logic; unknown != female
        self.sex_known = str(ped.get("sex")) in ("1", "2")
        self.has_hiconf = "ID=hiConfDeNovo" in vcf.raw_header


def base_row(trio_id, v, gt, mode, pair_id=""):
    return {
        "trio_id": trio_id, "mode": mode, "pair_id": pair_id,
        "chrom": v.CHROM, "pos": v.POS, "ref": v.REF, "alt": ",".join(v.ALT),
        "gene": A._str(v, "gene") or "", "symbol": A.symbol(v) or "",
        "consequence": A.consequence(v) or "", "impact": A.impact(v) or "",
        "grpmax_af": fmt(A.grpmax_af(v)), "max_af": fmt(A._max_float(v, "max_af")),
        "max_af_pops": A._str(v, "max_af_pops") or "", "cadd": fmt(A.cadd(v)),
        "clnsig": A.clnsig(v) or "",
        "child_gt": (v.gt_bases[gt.c] if v.gt_bases is not None else ""),
        "child_gq": fmt(G.gq(v, gt.c)), "child_dp": fmt(G.dp(v, gt.c)),
        "child_ab": fmt(G.allele_balance(v, gt.c)),
        "mother_gt": (v.gt_bases[gt.m] if v.gt_bases is not None else ""),
        "father_gt": (v.gt_bases[gt.d] if v.gt_bases is not None else ""),
        "hiConfDeNovo": ("1" if A.is_hiconf_denovo_for(v, gt.child_name) else ""),
        "review_prior_crosscheck": "", "flags": "",
    }


def screen_trio(trio_id, vcf, gt: Trio, cfg):
    thr = gt.thr
    dom_max = float(get(cfg, "filters.rarity.dominant_max", 1e-4))
    rec_max = float(get(cfg, "filters.rarity.recessive_max", 1e-2))
    require_hiconf = bool(get(cfg, "filters.denovo.use_hiconf_tag", True))
    # NB: filters.denovo.require_gnomad_absent_or_singleton is retired — it was implemented as
    # nhomalt > 1 (a HOMOZYGOTE-count test, never the allele-count test its name promised), and
    # nhomalt does not exist in the VEP cache. Removed rather than silently no-op'd.
    crosscheck = bool(get(cfg, "filters.denovo.crosscheck_prerefinement_pl", True))
    # Focus is INHERITED variation. De novo detection is retained for cross-reference
    # only (dedicated de novo filtering/review lives in separate machinery); the
    # dominant model — recurrent inherited rare functional hets — is the new emphasis.
    emit_denovo = bool(get(cfg, "inheritance.emit_denovo", True))
    emit_dominant = bool(get(cfg, "inheritance.emit_dominant", True))
    require_pass = bool(get(cfg, "filters.genotype_qc.require_pass", True))
    rec_strict = float(get(cfg, "filters.rarity.recessive_strict", 1e-3))

    def rare(v, limit):
        fr = A.frequency(v)
        return fr is None or fr < limit

    def tag_strict(r, v):
        """Flag a recessive/X-linked call whose frequency is below the high-confidence tier.

        Reads frequency() — the same chokepoint rare() uses — rather than a field getter
        directly; the two previously disagreed (this read faf95 with no grpmax fallback,
        so a variant with only a grpmax AF silently never earned the flag).
        """
        fr = A.frequency(v)
        if fr is not None and fr < rec_strict:
            r["flags"] = (r["flags"] + ";" if r["flags"] else "") + "high_conf_rarity"
        return r

    rows = []
    # gene -> list of (origin, variant, key); origin in {mat, pat, both, denovo}.
    # Feeds both compound-het pairing (recessive) and dominant (inherited het) calls.
    hets = {}

    for v in vcf:
        if require_pass and v.FILTER:  # cyvcf2 FILTER is None for PASS/'.'
            continue
        # sex unresolved -> ploidy-aware X/Y logic cannot be applied; skip sex chromosomes
        # rather than silently assume female (fail-soft; autosomal modes still run)
        if not gt.sex_known and G.is_sex_nonpar(v):
            continue
        # non-PAR chrY only exists in males; a female chrY non-PAR call is an artifact
        if G.is_y_nonpar(v) and not gt.child_male:
            continue
        c, d, m = gt.c, gt.d, gt.m
        gc, gd, gmm = v.gt_types[c], v.gt_types[d], v.gt_types[m]

        # male non-PAR chrX/chrY = hemizygous; a het call there is a QC red flag
        male_x = G.is_sex_nonpar(v) and gt.child_male

        # ---- de novo (SECONDARY / cross-reference only; review handled elsewhere) ----
        denovo_hit = False
        if emit_denovo and not male_x and gc == G.HET and gd == G.HOM_REF and gmm == G.HOM_REF:
            denovo_hit = True
        # male-X de novo: the son's single X comes from the MOTHER; the father transmits Y, so his
        # chrX is irrelevant — require only mother hom-ref, not the father.
        elif emit_denovo and male_x and gc == G.HOM_ALT and gmm == G.HOM_REF:
            denovo_hit = True
        if denovo_hit:
            child_kind = "denovo_child" if not male_x else "hom_alt"
            # parental cleanliness: both parents for an autosomal de novo; only the transmitting
            # mother for a male-X de novo (father's chrX is not transmitted to a son).
            parents_clean = G.sample_qc(v, m, thr, "clean_parent")
            if not male_x:
                parents_clean = parents_clean and G.sample_qc(v, d, thr, "clean_parent")
            ok = (G.sample_qc(v, c, thr, child_kind) and parents_clean and rare(v, dom_max))
            if male_x and (G.dp(v, c) or 0) < thr.denovo_min_dp:
                ok = False  # X/Y-hemizygous de novo still needs the deeper de novo DP floor
            # The gnomAD-homozygote gate that used to sit here is gone with nhomalt: the VEP
            # cache carries no homozygote count. rare(v, dom_max) above still applies the
            # frequency gate, which is the bulk of what it did. De novo is secondary here
            # (dedicated machinery owns it), so this is the cheapest place to absorb the loss.
            if require_hiconf and gt.has_hiconf and not A.is_hiconf_denovo_for(v, gt.child_name):
                ok = False  # tag exists in this callset but not a hiConf de novo for THIS child
            if ok:
                r = base_row(trio_id, v, gt, "denovo_x_hemi" if male_x else "denovo")
                if crosscheck:
                    r["review_prior_crosscheck"] = "1"
                rows.append(r)

        # a transmitting/carrier parent may be HET or (consanguinity, common-ish recessive allele,
        # affected parent) HOM_ALT — both carry a transmissible alt; QC per its own genotype.
        def carrier_ok(idx, gtype):
            if gtype == G.HET:
                return G.sample_qc(v, idx, thr, "het")
            if gtype == G.HOM_ALT:
                return G.sample_qc(v, idx, thr, "hom_alt")
            return False

        # ---- autosomal homozygous recessive: HOM_ALT child, both parents carriers (HET or HOM_ALT) ----
        if not male_x and not G.is_x_nonpar(v) and gc == G.HOM_ALT \
                and gd in (G.HET, G.HOM_ALT) and gmm in (G.HET, G.HOM_ALT):
            if (G.sample_qc(v, c, thr, "hom_alt") and carrier_ok(d, gd) and carrier_ok(m, gmm)
                    and rare(v, rec_max)):
                rows.append(tag_strict(base_row(trio_id, v, gt, "hom_recessive"), v))

        # ---- X-linked recessive, affected male: hemizygous son + carrier mother. The father
        #      transmits his Y (not his X) to a son, so his chrX genotype is IRRELEVANT and is not
        #      required — an affected/carrier father or a father chrX no-call must not drop the call. ----
        if male_x and gc == G.HOM_ALT and gmm in (G.HET, G.HOM_ALT):
            if G.sample_qc(v, c, thr, "hom_alt") and carrier_ok(m, gmm) and rare(v, rec_max):
                r = base_row(trio_id, v, gt, "x_linked_recessive")
                if gd in (G.HET, G.HOM_ALT):
                    r["flags"] = (r["flags"] + ";" if r["flags"] else "") + "father_carries_x_allele"
                rows.append(tag_strict(r, v))

        # ---- X-linked recessive, affected female: HOM_ALT daughter, carrier mother, hemizygous-
        #      affected father (he DOES transmit his X to a daughter) (docs §3.4) ----
        if (not male_x and G.is_x_nonpar(v) and gc == G.HOM_ALT
                and gmm in (G.HET, G.HOM_ALT) and gd == G.HOM_ALT):
            if (G.sample_qc(v, c, thr, "hom_alt") and carrier_ok(m, gmm)
                    and G.sample_qc(v, d, thr, "hom_alt") and rare(v, rec_max)):
                rows.append(tag_strict(base_row(trio_id, v, gt, "x_linked_recessive"), v))

        # ---- collect het candidates (het child, rare, parent-of-origin) ----
        #      The transmitting parent must be a QC-confident carrier (documented rule),
        #      so dominant and compound-het calls require parent genotype QC, not just child.
        if not male_x and gc == G.HET and G.sample_qc(v, c, thr, "het") and rare(v, rec_max):
            gene = A._str(v, "gene") or A.symbol(v)
            if gene:
                mom_carries = gmm in (G.HET, G.HOM_ALT)
                dad_carries = gd in (G.HET, G.HOM_ALT)
                mom_ok = ((gmm == G.HET and G.sample_qc(v, m, thr, "het"))
                          or (gmm == G.HOM_ALT and G.sample_qc(v, m, thr, "hom_alt")))
                dad_ok = ((gd == G.HET and G.sample_qc(v, d, thr, "het"))
                          or (gd == G.HOM_ALT and G.sample_qc(v, d, thr, "hom_alt")))
                if gmm == G.HOM_REF and gd == G.HOM_REF:
                    # only a genuinely de novo het may pair in trans; require BOTH parents to
                    # pass cleanliness QC, else a dropped-out parental het masquerades as de novo
                    # and gets paired cis with a real variant from that same parent
                    origin = ("denovo" if (G.sample_qc(v, m, thr, "clean_parent")
                                           and G.sample_qc(v, d, thr, "clean_parent")) else None)
                elif mom_carries and dad_carries:
                    origin = "both" if (mom_ok and dad_ok) else None
                elif mom_carries:
                    origin = "mat" if mom_ok else None
                elif dad_carries:
                    origin = "pat" if dad_ok else None
                else:
                    origin = None                # a parent no-call — inheritance unestablished
                if origin:
                    key = f"{v.CHROM}:{v.POS}:{v.REF}:{v.ALT[0]}"
                    hets.setdefault(gene, []).append((origin, v, key))

    # ---- compound het (recessive): trans pairs with determinable parent-of-origin ----
    consumed = set()
    pair_n = 0
    for gene, cands in hets.items():
        by = {"mat": [], "pat": [], "denovo": [], "both": []}
        for origin, v, key in cands:
            by[origin].append((v, key))
        pairs = [(a, b) for a in by["mat"] for b in by["pat"] + by["denovo"]]
        pairs += [(a, b) for a in by["pat"] for b in by["denovo"]]
        for (va, ka), (vb, kb) in pairs:
            pair_n += 1
            pid = f"{trio_id}:CH{pair_n}"
            consumed.add(ka)
            consumed.add(kb)
            for v in (va, vb):
                rows.append(tag_strict(base_row(trio_id, v, gt, "compound_het", pid), v))

    # ---- dominant (inherited het): rare, functional, transmitted from >=1 parent,
    #      not part of a compound-het pair. This is the recurrence signal Step 6 tallies. ----
    if emit_dominant:
        for gene, cands in hets.items():
            for origin, v, key in cands:
                if origin in ("mat", "pat", "both") and key not in consumed and rare(v, dom_max):
                    r = base_row(trio_id, v, gt, "dominant")
                    r["flags"] = f"origin={origin}"
                    rows.append(r)
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, help="trios.candidates.tsv from Step 4")
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True, help="candidate calls TSV")
    ap.add_argument("--qc-report", default="", help="Step 0 qc_report.tsv (for inferred kid sex)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    thr = G.GtThresholds.from_config(cfg, get)

    # inferred sex per trio from Step 0 QC (used when the generated PED has sex unknown)
    sex_map = {}
    if args.qc_report and __import__("os").path.exists(args.qc_report):
        import csv as _csv
        with open(args.qc_report) as fh:
            for r in _csv.DictReader(fh, delimiter="\t"):
                if r.get("inferred_sex"):
                    sex_map[r.get("trio_id")] = r["inferred_sex"]

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
        # If the generated PED has kid sex unknown, use Step 0's inferred sex so
        # X-linked / hemizygous logic can fire correctly. If still unresolved, warn — Step 5
        # will skip X/Y modes for this trio rather than silently assume female.
        if str(ped.get("sex")) in ("0", "", "None"):
            if trio_id in sex_map:
                ped["sex"] = sex_map[trio_id]
            else:
                sys.stderr.write(f"WARN: {trio_id}: child sex unresolved (no Step-0 inference); "
                                 f"X/Y-linked modes skipped for this trio (autosomal modes still run)\n")
        vcf = VCF(vcf_path)
        try:
            gt = Trio(vcf, ped, thr)
        except KeyError as e:
            sys.stderr.write(f"WARN: {trio_id}: {e}; skipping\n")
            vcf.close()
            continue
        trio_rows = screen_trio(trio_id, vcf, gt, cfg)
        vcf.close()
        n_trios += 1
        # per-trio audit: total calls + counts by mode
        tmodes = {}
        for row in trio_rows:
            tmodes[row["mode"]] = tmodes.get(row["mode"], 0) + 1
        audit.record("05_inheritance", "candidate_calls", len(trio_rows), scope=trio_id)
        for mode, c in sorted(tmodes.items()):
            audit.record("05_inheritance", f"mode.{mode}", c, scope=trio_id)
        all_rows.extend(trio_rows)

    with open(args.out, "w") as out:
        out.write("\t".join(COLS) + "\n")
        for r in all_rows:
            out.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")

    by_mode = {}
    for r in all_rows:
        by_mode[r["mode"]] = by_mode.get(r["mode"], 0) + 1
    audit.record("05_inheritance", "trios_screened", n_trios)
    audit.record("05_inheritance", "candidate_calls_total", len(all_rows))
    for mode, c in sorted(by_mode.items()):
        audit.record("05_inheritance", f"mode.{mode}", c)
    sys.stderr.write(
        f"Step 5 complete: {len(all_rows)} candidate calls across {n_trios} trios "
        f"-> {args.out}\n  by mode: {by_mode}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
