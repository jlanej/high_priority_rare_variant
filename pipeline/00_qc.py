#!/usr/bin/env python3
"""Pipeline Step 0: per-trio QC — advisory pass/flag list (garbage-in guard, human-reviewed).

Computes Mendelian-error rate, chrX-inferred sex (all three members), per-member no-call rate,
and contamination per trio and folds them into an `overall_pass` flag. This is ADVISORY: flags
are surfaced to the reviewer (xlsx QC sheet + IGV sample_qc.tsv) but no step auto-excludes a
flagged trio (see the overall_pass comment below).

Mostly self-contained (no extra resources): for each trio computes five gates and
flags trios that fail any of them:
  * Mendelian-error rate — a proxy for sample swaps / contamination. NOT for a father/mother
    label swap: the rule is symmetric under exchanging the parents, so that swap is invisible
    to it by construction — which is what the next gate is for.
  * chrX-heterozygosity sex inference: the child's vs. the PED sex (when the PED states one),
    and BOTH PARENTS' vs. their roles (the PED always states those). A male in the mother slot
    is the one direct detector of transposed parents.
  * Per-member genotype no-call rate on the scanned sites. A jointly genotyped trio carries an
    affirmative `0/0` for a non-carrier parent; a merge of single-sample callsets carries `./.`
    there instead, and Step 5 then loses every de novo and marks every inherited call
    `origin_unverified` — while the MIE denominator, which skips no-call sites, reports a
    CLEANER trio. The rate makes that input shape visible before it silently empties a run.
  * Contamination — verifyBamID FREEMIX if a directory of ``*.selfSM`` files is
    configured (resources.selfsm_dir), else a VCF-only CHARR estimate (reference-read
    fraction at high-quality hom-ALT SNV sites). Mirrors the group's DNM freemix QC.

Trio kid/dad/mom roles come from upstream peddy; this step is the guard against the
less-well-curated trios. See docs/inheritance_and_genotype_qc.md.

The Mendelian-error / CHARR scan and the chrX sex scan are each CAPPED at ``qc.max_sites``
QC-passing sites (default 200,000; ``--max-sites`` overrides): on WGS the MIE rate is therefore
measured on a prefix of the genome, not genome-wide — adequate as a swap/contamination detector,
but a localised MIE cluster (UPD/CNV) outside that prefix is invisible to it.

Usage:
  00_qc.py --manifest trios.tsv --config config.yaml --out qc_report.tsv [--max-sites N]
"""
from __future__ import annotations

import argparse
import csv
import itertools
import sys

from cyvcf2 import VCF

from hprv import audit
from hprv import contamination as C
from hprv import genotype as G
from hprv.config import get, load_config
from hprv.ped import parse_ped

HOM_REF, HET, HOM_ALT = G.HOM_REF, G.HET, G.HOM_ALT


def mendelian_violation(gc, gd, gm) -> bool:
    """True if child genotype is impossible given parents (biallelic autosomal)."""
    if gc == HOM_ALT:      # needs an alt from each parent
        return gd == HOM_REF or gm == HOM_REF
    if gc == HOM_REF:      # needs a ref from each parent
        return gd == HOM_ALT or gm == HOM_ALT
    if gc == HET:          # needs one alt and one ref available
        return (gd == HOM_REF and gm == HOM_REF) or (gd == HOM_ALT and gm == HOM_ALT)
    return False


def infer_sex(x_het, x_hom, cutoff, min_sites):
    """chrX-heterozygosity sex call: '1' male, '2' female, None when too few informative sites."""
    total = x_het + x_hom
    if total < min_sites or total == 0:
        return None
    return "1" if (x_het / total) < cutoff else "2"


def _pick_x_contig(seqnames):
    for want in ("chrX", "X"):
        if want in seqnames:
            return want
    return None


def scan_sex(vcf_path, child_id, thr, max_x):
    """chrX non-PAR het/hom counts in the child, as a DEDICATED pass.

    chrX sorts after all autosomes, so an autosomal MIE cap in the main pass would
    otherwise starve sex inference on WGS trios. Uses the index to jump straight to
    chrX; falls back to a full scan (filtered to chrX) when the VCF is unindexed.
    """
    vcf = VCF(vcf_path)
    ci = {s: i for i, s in enumerate(vcf.samples)}.get(child_id)
    x = _pick_x_contig(vcf.seqnames)
    if ci is None or x is None:
        vcf.close()
        return 0, 0
    try:
        it = vcf(x)                       # indexed region jump (preferred)
        # cyvcf2.VCF.__call__ is a GENERATOR function, so `vcf(x)` returns without running its
        # body — the tabix load and its `assert self.idx != NULL` fire on the first next(), i.e.
        # inside the for-loop below and OUTSIDE this try. Without forcing the first element here
        # the fallback is dead code and an unindexed VCF raises AssertionError straight through
        # qc_trio/main into `set -euo pipefail`, killing the whole run. Unindexed trio VCFs are a
        # supported input, and Step 0 runs before anything indexes them.
        it = itertools.chain([next(it)], it)
    except StopIteration:
        it = iter(())                     # chrX in the header but no records there
    except Exception:
        it = vcf                          # unindexed: full scan, filtered below
    x_het = x_hom = 0
    for v in it:
        if v.CHROM.replace("chr", "") != "X" or G.in_par_x(v) or len(v.ALT) != 1:
            continue
        gq, dp = G.gq(v, ci), G.dp(v, ci)
        if gq is None or gq < thr.min_gq or dp is None or dp < thr.min_dp:
            continue
        gt = v.gt_types[ci]
        if gt == HET:
            x_het += 1
        elif gt == HOM_ALT:
            x_hom += 1
        if max_x and (x_het + x_hom) >= max_x:
            break
    vcf.close()
    return x_het, x_hom


def qc_trio(vcf_path, ped, thr, max_sites, sex_cutoff=0.10, sex_min_sites=20):
    vcf = VCF(vcf_path)
    samples = {s: i for i, s in enumerate(vcf.samples)}
    for role in ("child", "father", "mother"):
        if ped[role] not in samples:
            vcf.close()
            return None
    c, d, m = samples[ped["child"]], samples[ped["father"]], samples[ped["mother"]]

    # sex inference runs as its own chrX pass so the autosomal MIE cap can't disable it
    x_het, x_hom = scan_sex(vcf_path, ped["child"], thr, max_x=(max_sites or 0))
    # ...and for the PARENTS. Their sex is the one thing the PED asserts with certainty
    # (father = 1, mother = 2), so this is the only sex check that is reachable in the shipped
    # flow — the generated PED leaves the child's sex unknown, which made `sex_match` a constant
    # 1 — and it is the one direct detector of a transposed mother/father, which the Mendelian-
    # error rate cannot see (mendelian_violation is symmetric in gd/gm).
    px = {role: scan_sex(vcf_path, ped[role], thr, max_x=(max_sites or 0))
          for role in ("father", "mother")}

    considered = errors = 0
    n_records = 0                            # biallelic autosomal records seen (pre-QC)
    nocall = {"kid": 0, "dad": 0, "mom": 0}  # per-member `./.` among them
    cref = {"kid": 0, "dad": 0, "mom": 0}   # CHARR: ref reads at hom-alt SNV sites
    cdp = {"kid": 0, "dad": 0, "mom": 0}    # CHARR: total (ref+alt) reads there
    for v in vcf:
        chrom = v.CHROM.replace("chr", "")
        if chrom in ("X", "Y", "MT", "M"):  # X counted in scan_sex; Y/MT out of scope
            continue
        if len(v.ALT) != 1:  # biallelic only for the simple MIE rule
            continue
        # CHARR: accumulate reference reads at each member's high-quality hom-alt SNV sites
        if len(v.REF) == 1 and len(v.ALT[0]) == 1:
            for role, i in (("kid", c), ("dad", d), ("mom", m)):
                if v.gt_types[i] == HOM_ALT:
                    gq, dp = G.gq(v, i), G.dp(v, i)
                    ra, aa = G.ref_ad(v, i), G.alt_ad(v, i)
                    if (gq and gq >= thr.min_gq and dp and dp >= thr.min_dp
                            and ra is not None and aa is not None and (ra + aa) > 0):
                        cref[role] += ra
                        cdp[role] += ra + aa
        gc, gd, gm = v.gt_types[c], v.gt_types[d], v.gt_types[m]
        n_records += 1
        for role, g in (("kid", gc), ("dad", gd), ("mom", gm)):
            if g == G.UNKNOWN:
                nocall[role] += 1
        if G.UNKNOWN in (gc, gd, gm):
            continue                          # a no-call site is out of the MIE denominator
        if any(G.gq(v, i) is None or G.gq(v, i) < thr.min_gq for i in (c, d, m)):
            continue
        if any(G.dp(v, i) is None or G.dp(v, i) < thr.min_dp for i in (c, d, m)):
            continue
        considered += 1
        if mendelian_violation(gc, gd, gm):
            errors += 1
        if max_sites and considered >= max_sites:
            break
    vcf.close()

    mie_rate = (errors / considered) if considered else None
    x_total = x_het + x_hom
    x_het_ratio = (x_het / x_total) if x_total else None
    # only call sex with enough informative chrX sites, else leave it unknown (fail-soft)
    inferred = infer_sex(x_het, x_hom, sex_cutoff, sex_min_sites)
    return {
        "n_sites": considered, "mie_errors": errors, "mie_rate": mie_rate,
        "x_sites": x_total, "x_het_ratio": x_het_ratio, "inferred_sex": inferred,
        "dad_sex": infer_sex(*px["father"], sex_cutoff, sex_min_sites),
        "mom_sex": infer_sex(*px["mother"], sex_cutoff, sex_min_sites),
        "n_records": n_records,
        "nocall_rate": {r: ((nocall[r] / n_records) if n_records else None) for r in ("kid", "dad", "mom")},
        "charr": {r: C.charr(cref[r], cdp[r]) for r in ("kid", "dad", "mom")},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-sites", type=int, default=None,
                    help="cap on QC-passing sites scanned per trio (MIE/CHARR and the chrX sex "
                         "scan); default qc.max_sites from the config (200000)")
    ap.add_argument("--mie-threshold", type=float, default=None,
                    help="override qc.mie_max from config")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    thr = G.GtThresholds.from_config(cfg, get)
    selfsm = C.read_selfsm(get(cfg, "resources.selfsm_dir", ""))   # verifyBamID FREEMIX (optional)
    freemix_thr = float(get(cfg, "qc.freemix_threshold", 0.05))
    charr_thr = float(get(cfg, "qc.charr_threshold", 0.02))
    # thresholds are config defaults (canonical-defaults table), not hardcoded law
    mie_thr = args.mie_threshold if args.mie_threshold is not None else float(get(cfg, "qc.mie_max", 0.02))
    sex_cutoff = float(get(cfg, "qc.x_het_male_max", 0.10))
    sex_min = int(get(cfg, "qc.sex_min_sites", 20))
    max_sites = args.max_sites if args.max_sites is not None else int(get(cfg, "qc.max_sites", 200000))
    max_nocall = float(get(cfg, "qc.max_nocall_rate", 0.10))

    with open(args.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    cols = ["trio_id", "n_sites", "mie_errors", "mie_rate", "x_sites", "x_het_ratio", "inferred_sex",
            "ped_sex", "sex_match", "dad_inferred_sex", "mom_inferred_sex", "parent_sex_flag",
            "mie_flag", "kid_contam", "dad_contam", "mom_contam", "contam_source", "contam_flag",
            "n_records", "kid_nocall_rate", "dad_nocall_rate", "mom_nocall_rate", "nocall_flag",
            "overall_pass"]
    n_fail = n_done = 0
    n_flag = {"parent_sex": 0, "nocall": 0}

    def member_contam(sample, role, charr_map):
        """Return (value, source, flagged) for one member — verifyBamID freemix if present,
        else the VCF-only CHARR proxy."""
        if sample in selfsm:
            return selfsm[sample], "freemix", selfsm[sample] > freemix_thr
        cv = charr_map.get(role)
        return cv, "charr", (cv is not None and cv > charr_thr)

    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            tid, vcf_path, ped_path = r.get("trio_id"), r.get("vcf"), r.get("ped")
            ped = parse_ped(ped_path)
            if not ped:
                sys.stderr.write(f"WARN: no PED for {tid}; skipping QC\n")
                continue
            res = qc_trio(vcf_path, ped, thr, max_sites, sex_cutoff, sex_min)
            if res is None:
                sys.stderr.write(f"WARN: {tid}: PED samples not in VCF; skipping\n")
                continue
            ped_sex = str(ped["sex"])
            # Sex-swap gate: only a MISMATCH against a KNOWN expected sex is a failure. The generated
            # PED leaves kid sex unknown ('0'), so an unknown/absent ped_sex = "no expectation" = pass
            # (mirrors the un-inferable branch); otherwise every inferable trio would be flagged.
            sex_known = ped_sex not in ("0", "", "None")
            sex_match = "0" if (sex_known and res["inferred_sex"] is not None
                                and res["inferred_sex"] != ped_sex) else "1"
            mie_flag = "1" if (res["mie_rate"] is not None and res["mie_rate"] > mie_thr) else "0"

            contam, flags, src = {}, False, "charr"
            for role, key in (("kid", "child"), ("dad", "father"), ("mom", "mother")):
                val, s, fl = member_contam(ped[key], role, res["charr"])
                contam[role] = val
                flags = flags or fl
                if s == "freemix":
                    src = "freemix"
            contam_flag = "1" if flags else "0"
            # "charr" with every value missing means nobody was measured, not that CHARR ran
            if src == "charr" and all(val is None for val in contam.values()):
                src = "none"
            # a father whose own chrX reads female, or a mother whose reads male: the PED roles
            # are wrong for this trio (or a sample is), and every parent-of-origin call would be
            # inverted. Only a positive inference counts — too few X sites is "no expectation".
            parent_sex_flag = "1" if (res["dad_sex"] == "2" or res["mom_sex"] == "1") else "0"
            nc = res["nocall_rate"]
            nocall_flag = "1" if any(nc[r] is not None and nc[r] > max_nocall for r in nc) else "0"
            n_flag["parent_sex"] += parent_sex_flag == "1"
            n_flag["nocall"] += nocall_flag == "1"

            # ADVISORY only: overall_pass folds the flags into one column that is SURFACED to the
            # mandatory human reviewer (xlsx QC sheet + IGV sample_qc.tsv) but does NOT auto-
            # exclude a trio — Steps 1/2/4 run over the full resolved manifest and Step 6 never reads
            # it, so a flagged trio still contributes calls and recurrence pending human review.
            # Automated gating is deliberately deferred (never-drop ethos + the contamination proxy's
            # limited sensitivity, see contamination.py); it would be a config-gated policy change.
            overall = "1" if (mie_flag == "0" and sex_match == "1" and contam_flag == "0"
                              and parent_sex_flag == "0" and nocall_flag == "0") else "0"
            if overall == "0":
                n_fail += 1
            n_done += 1
            row = {
                "trio_id": tid, "n_sites": res["n_sites"], "mie_errors": res["mie_errors"],
                "mie_rate": ("" if res["mie_rate"] is None else f"{res['mie_rate']:.4g}"),
                "x_sites": res["x_sites"],
                "x_het_ratio": ("" if res["x_het_ratio"] is None else f"{res['x_het_ratio']:.3g}"),
                "inferred_sex": res["inferred_sex"] or "", "ped_sex": ped_sex,
                "sex_match": sex_match, "dad_inferred_sex": res["dad_sex"] or "",
                "mom_inferred_sex": res["mom_sex"] or "", "parent_sex_flag": parent_sex_flag,
                "mie_flag": mie_flag,
                "kid_contam": ("" if contam["kid"] is None else f"{contam['kid']:.4g}"),
                "dad_contam": ("" if contam["dad"] is None else f"{contam['dad']:.4g}"),
                "mom_contam": ("" if contam["mom"] is None else f"{contam['mom']:.4g}"),
                "contam_source": src, "contam_flag": contam_flag,
                "n_records": res["n_records"],
                "kid_nocall_rate": ("" if nc["kid"] is None else f"{nc['kid']:.4g}"),
                "dad_nocall_rate": ("" if nc["dad"] is None else f"{nc['dad']:.4g}"),
                "mom_nocall_rate": ("" if nc["mom"] is None else f"{nc['mom']:.4g}"),
                "nocall_flag": nocall_flag, "overall_pass": overall,
            }
            out.write("\t".join(str(row[c]) for c in cols) + "\n")
            # Step 0 recorded nothing in the audit before; the flags are the useful part
            for metric in ("overall_pass", "mie_flag", "sex_match", "parent_sex_flag",
                           "contam_flag", "nocall_flag"):
                audit.record("00_qc", metric, row[metric], scope=tid)
            if parent_sex_flag == "1":
                sys.stderr.write(f"WARN: {tid}: PARENT SEX MISMATCH — father's chrX reads "
                                 f"{res['dad_sex'] or '?'}, mother's reads {res['mom_sex'] or '?'} "
                                 "(1=male, 2=female). The trios-file roles for this trio are probably "
                                 "transposed; every parent-of-origin call would be inverted.\n")
            if nocall_flag == "1":
                sys.stderr.write(f"WARN: {tid}: no-call rate above qc.max_nocall_rate ({max_nocall}) "
                                 f"— kid {row['kid_nocall_rate']} dad {row['dad_nocall_rate']} mom "
                                 f"{row['mom_nocall_rate']}. Was this trio jointly genotyped? Step 5 "
                                 "treats a parental ./. as uninformative, never as 0/0.\n")

    audit.record("00_qc", "trios_qc", n_done)
    audit.record("00_qc", "trios_flagged", n_fail)
    audit.record("00_qc", "trios_parent_sex_flag", n_flag["parent_sex"])
    audit.record("00_qc", "trios_nocall_flag", n_flag["nocall"])
    sys.stderr.write(f"Step 0 complete: QC report -> {args.out} ({n_fail} trio(s) flagged; "
                     f"parent-sex {n_flag['parent_sex']}, no-call {n_flag['nocall']})\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
