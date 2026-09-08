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

Emits one TSV of candidate calls across all trios: the curated columns (COLS) first, then a
DROPLESS `info_<ID>` block carrying every INFO field of the per-trio candidate VCF verbatim (see
INFO_PREFIX). See docs/inheritance_and_genotype_qc.md and docs/pipeline_design.md (Step 5).
"""
from __future__ import annotations

import argparse
from collections import Counter
import re
import sys

from cyvcf2 import VCF

from hprv import annotations as A
from hprv import audit
from hprv import genotype as G
from hprv.config import get, get_bool, load_config, validate_filters
from hprv.ped import parse_ped

COLS = [
    "trio_id", "mode", "pair_id", "chrom", "pos", "ref", "alt", "gene", "symbol",
    # frequency() reads ONE oracle for the whole run (resources.gnomad.oracle: faf95 by default,
    # else the grpmax point-estimate proxy); the arms never cross. max_af/max_af_pops ride along
    # for review only, so a curator can see when a call would have been driven by a founder-group
    # frequency that grpmax deliberately ignores. grpmax_af (the proxy) and faf95 (the CI-corrected
    # value) BOTH ride along so a reviewer can see the gap the CI correction closed.
    # faf95_group names the ancestry group that produced it: the FAF set is GRPMAX_POPS + `mid`,
    # so this is how the one deviation stays visible. nhomalt is the recessive false-positive tell.
    # rarity_af is THE value every gate used (annotations.frequency()) and rarity_oracle is its
    # provenance — resolved ONCE here so Step 9 ranks on exactly what the screen gated on instead
    # of re-deriving it from the raw columns and risking divergence. Both raw inputs ride along.
    "consequence", "impact",
    # HGVS on the transcript split-vep selected (VEP --hgvs; CSQ HGVSc/HGVSp lifted in Step 2).
    # These sat in the per-trio VCF INFO from Step 2 onward and were dropped ONLY by this
    # projection: the VCF->TSV column list is the one place an annotation can vanish, which is
    # also why the `info_*` block exists (INFO_PREFIX below).
    "hgvsc", "hgvsp",
    "rarity_af", "rarity_oracle", "rarity_basis",
    "grpmax_af", "faf95", "faf95_group", "nhomalt",
    "max_af", "max_af_pops", "cadd", "spliceai_ds",
    # Calibrated missense predictors. Inert at the SCREEN by construction (missense is
    # IMPACT=MODERATE and selection.py returns at the impact rung), carried here purely so
    # Step 9's missense tier can be calibrated rather than an off-label CADD rank.
    "revel", "alphamissense", "alphamissense_class",
    # clnsig is the VEP cache's CLIN_SIG; clinvar_stars comes from the ClinVar VCF transfer and
    # is BLANK when that transfer did not run. Blank != 0 stars — see annotations.clinvar_stars.
    "clnsig", "clinvar_stars", "child_gt", "child_gq", "child_dp", "child_ab",
    "mother_gt", "father_gt", "hiConfDeNovo", "review_prior_crosscheck", "flags",
]

# --- dropless pass-through of the per-trio VCF's INFO -------------------------------------------
# candidates.calls.tsv used to be the ONE place in the pipeline where an annotation could vanish
# without a trace: Step 2 lifts ~45 INFO fields, Steps 3-4 carry every one of them into the per-trio
# candidate VCF, and this step then projected each record onto the fixed COLS list above — so
# vep_HGVSc / vep_HGVSp / vep_Feature / vep_MANE_SELECT / hprv_keep_reason and the rest never reached
# a TSV, the workbook or the igv.js review table. The curated columns stay exactly where they are
# (downstream readers are name-keyed); AFTER them, every INFO field declared in the candidate VCF
# header is emitted VERBATIM under `info_<ID>` — the union over trios, in header order. The prefix
# is load-bearing: the raw `hiConfDeNovo` (a comma list of children) would otherwise collide with
# the curated `hiConfDeNovo` column (this child: 1/blank). A VCF `.` reads blank, matching every
# other missing value in the table, and a Flag reads `1`; nothing else is transformed.
INFO_PREFIX = "info_"
_INFO_ID_RE = re.compile(r"^##INFO=<ID=([^,>]+)", re.MULTILINE)


def info_ids_from_header(raw_header: str):
    """Every INFO ID the VCF header declares, in header order."""
    return _INFO_ID_RE.findall(raw_header or "")


def _typed_info(x):
    """Fallback rendering of a cyvcf2-typed INFO value when the raw VCF line is unavailable."""
    if x is None:
        return ""
    if isinstance(x, bool):
        return "1" if x else ""
    if isinstance(x, (tuple, list)):
        return ",".join("" if e is None else str(e) for e in x)
    s = str(x)
    return "" if s == "." else s


def info_values(v, ids):
    """`{info_<ID>: verbatim value}` for every ID in `ids`, read from the record's own INFO column.

    The VCF LINE (`str(v)`) is parsed rather than the typed `v.INFO`, because htslib stores Float
    as float32 and cyvcf2 hands back the widened double — `0.000162583` would come out as
    `0.00016258300165`, which is not what the VCF says. A record that cannot render itself as a
    VCF line (the unit-test harness) falls back to the typed values. Percent-encoded characters
    (`%3B`, `%3D`) stay encoded: verbatim means verbatim.
    """
    if not ids:
        return {}
    raw = None
    try:
        cols = str(v).rstrip("\n").split("\t")
        if len(cols) >= 8:
            raw = {}
            for kv in cols[7].split(";"):
                if not kv or kv == ".":
                    continue
                k, sep, val = kv.partition("=")
                raw[k] = ("" if val == "." else val) if sep else "1"   # a Flag: present = 1
    except Exception:  # noqa: BLE001 — a record that cannot render falls back to typed values
        raw = None
    out = {}
    for i in ids:
        out[INFO_PREFIX + i] = raw.get(i, "") if raw is not None else _typed_info(v.INFO.get(i))
    return out


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
        # Every INFO field this candidate VCF declares, for the dropless `info_*` pass-through
        # (INFO_PREFIX). Read from the header ONCE per trio; main() unions the blocks over trios.
        self.info_ids = info_ids_from_header(vcf.raw_header)
        self.info_cols = [INFO_PREFIX + i for i in self.info_ids]


def base_row(trio_id, v, gt, mode, pair_id="", cfg=None):
    row = {
        "trio_id": trio_id, "mode": mode, "pair_id": pair_id,
        "chrom": v.CHROM, "pos": v.POS, "ref": v.REF, "alt": ",".join(v.ALT),
        "gene": A._str(v, "gene") or "", "symbol": A.symbol(v) or "",
        "consequence": A.consequence(v) or "", "impact": A.impact(v) or "",
        "hgvsc": A._str(v, "hgvsc") or "", "hgvsp": A._str(v, "hgvsp") or "",
        "rarity_af": fmt(A.frequency(v, cfg)), "rarity_oracle": A.rarity_oracle(cfg),
        "rarity_basis": A.rarity_basis(v, cfg),
        "grpmax_af": fmt(A.grpmax_af(v)), "faf95": fmt(A.faf95(v)),
        "faf95_group": A.faf95_group(v) or "", "nhomalt": fmt(A.nhomalt(v)),
        "max_af": fmt(A._max_float(v, "max_af")),
        "max_af_pops": A._str(v, "max_af_pops") or "", "cadd": fmt(A.cadd(v)),
        "spliceai_ds": fmt(A.spliceai_ds(v)),   # max SpliceAI delta score, for reviewer tiering
        "revel": fmt(A.revel(v)), "alphamissense": fmt(A.alphamissense(v)),
        "alphamissense_class": A.alphamissense_class(v) or "",
        "clnsig": A.clnsig(v) or "",
        # fmt() maps None -> "" which is exactly right: blank means the ClinVar transfer did not
        # run (nobody looked), and must never be read as 0 stars ("no assertion criteria").
        "clinvar_stars": fmt(A.clinvar_stars(v)),
        "child_gt": (v.gt_bases[gt.c] if v.gt_bases is not None else ""),
        "child_gq": fmt(G.gq(v, gt.c)), "child_dp": fmt(G.dp(v, gt.c)),
        "child_ab": fmt(G.allele_balance(v, gt.c)),
        "mother_gt": (v.gt_bases[gt.m] if v.gt_bases is not None else ""),
        "father_gt": (v.gt_bases[gt.d] if v.gt_bases is not None else ""),
        "hiConfDeNovo": ("1" if A.is_hiconf_denovo_for(v, gt.child_name) else ""),
        "review_prior_crosscheck": "", "flags": "",
    }
    # The dropless block: every INFO field the candidate VCF carries, verbatim, AFTER the curated
    # columns (INFO_PREFIX). Never shortened, never re-derived.
    row.update(info_values(v, gt.info_ids))
    return row


def screen_trio(trio_id, vcf, gt: Trio, cfg):
    thr = gt.thr
    dom_max = float(get(cfg, "filters.rarity.dominant_max", 1e-4))
    rec_max = float(get(cfg, "filters.rarity.recessive_max", 1e-2))
    # Boolean knobs are read strictly: `bool(value)` turned a quoted or ${ENV}-templated "false"
    # into True, so the gate a user had just switched off stayed on.
    require_hiconf = get_bool(cfg, "filters.denovo.use_hiconf_tag", True)
    # NB: filters.denovo.require_gnomad_absent_or_singleton is retired — it was implemented as
    # nhomalt > 1 (a HOMOZYGOTE-count test, never the allele-count test its name promised), and
    # nhomalt is absent from the VEP cache; it arrives only with the OPTIONAL gnomAD joint slim
    # (resources.gnomad.sites_slim). The old de-novo `nhomalt <= 1` condition stays removed: it
    # would silently no-op whenever the slim is not configured. nhomalt is instead REPORTED per
    # variant, and Step 9 flags biallelic calls gnomAD already carries homozygotes for.
    crosscheck = get_bool(cfg, "filters.denovo.crosscheck_prerefinement_pl", True)
    # Focus is INHERITED variation. De novo detection is retained for cross-reference
    # only (dedicated de novo filtering/review lives in separate machinery); the
    # dominant model — recurrent inherited rare functional hets — is the new emphasis.
    emit_denovo = get_bool(cfg, "inheritance.emit_denovo", True)
    emit_dominant = get_bool(cfg, "inheritance.emit_dominant", True)
    require_pass = get_bool(cfg, "filters.genotype_qc.require_pass", True)
    rec_strict = float(get(cfg, "filters.rarity.recessive_strict", 1e-3))

    def rare(v, limit):
        fr = A.frequency(v, cfg)
        return fr is None or fr < limit

    def tag_strict(r, v):
        """Flag a recessive/X-linked call whose frequency is below the high-confidence tier.

        Reads frequency() — the same chokepoint rare() uses — rather than a field getter
        directly; the two previously disagreed (this read faf95 with no grpmax fallback,
        so a variant with only a grpmax AF silently never earned the flag).
        """
        fr = A.frequency(v, cfg)
        if fr is not None and fr < rec_strict:
            r["flags"] = (r["flags"] + ";" if r["flags"] else "") + "high_conf_rarity"
        return r

    rows = []
    # gene -> list of (origin, variant, key); origin in {mat, pat, both, denovo}.
    # Feeds both compound-het pairing (recessive) and dominant (inherited het) calls.
    hets = {}
    # ACCOUNTING. This step used to record only what it EMITTED, so a candidate that reached it and
    # matched no mode left no trace anywhere — and because a compound-het leg is emitted once per
    # pair, the per-trio funnel could read "31 candidate genotypes in, 35 calls out" while 4 of the
    # 31 had produced nothing (the shipped integration fixture did exactly that). Every record is
    # now counted on the way in, every early `continue` is counted by reason, and every examined
    # variant that produced no row is classified by the FIRST condition that removed it
    # (why_no_row below), so a reviewer can reconcile, per trio,
    #     variants_examined == skipped.* + variants_with_call + variants_no_row
    # and read a negative result as "examined and rejected for <reason>" rather than "nothing was
    # there". Steps 1-4 already report a reason for every drop; this brings Step 5 to that standard.
    stats = {"examined": 0, "annotated": 0, "with_call": 0, "skipped": Counter(),
             "no_row": Counter(), "clinvar_plp_ge_recessive_max": 0, "clinvar_plp_no_row": 0}
    pending = {}     # key -> (no-row reason, is_plp): produced no row in the loop, not pooled
    collected = {}   # key -> (variant, is_plp, origin): pooled into `hets`; resolved after pairing
    # A ClinVar P/LP variant the child carries can be kept by Step 3's clinvar_plp rescue and then
    # fail every Step-5 mode. Two counters, because the band matters: `ge_recessive_max` is the
    # original metric (frequency() >= recessive_max, so no mode could ever have fired) and
    # `clinvar_plp_no_row` is the complete one — every carried P/LP allele that yielded no row for
    # ANY reason, which includes the [dominant_max, recessive_max) band the first counter missed
    # (an unpaired inherited het there fails the dominant gate and is emitted under no mode).
    n_plp_inert = 0

    def _carrier_qc(v, idx, gtype):
        if gtype == G.HET:
            return G.sample_qc(v, idx, thr, "het")
        if gtype == G.HOM_ALT:
            return G.sample_qc(v, idx, thr, "hom_alt")
        return False

    def _hiconf_blocks(v):
        return require_hiconf and gt.has_hiconf and not A.is_hiconf_denovo_for(v, gt.child_name)

    def why_no_row(v, gc, gd, gmm, male_x, male_x_chrx):
        """The FIRST condition, in the order the mode logic applies them, that left this examined
        variant with no row and not pooled for pairing. Mirrors the branches above it exactly —
        keep the two in step when a mode's conditions change."""
        c, d, m = gt.c, gt.d, gt.m
        if gc not in (G.HET, G.HOM_ALT):
            return "child_not_carrier"        # Step 4 keeps a locus any member carries
        if G.is_y_nonpar(v):
            return "chry"                     # no Y-linked model; documented
        if male_x and gc == G.HET:
            return "male_x_het"               # hemizygous het = QC red flag, never a call
        if gc == G.HET:
            if not G.sample_qc(v, c, thr, "het"):
                return "qc_child"
            if gd == G.HOM_REF and gmm == G.HOM_REF:          # de novo shape
                clean = G.sample_qc(v, m, thr, "clean_parent") and G.sample_qc(v, d, thr, "clean_parent")
                if not clean:
                    return "qc_parent"        # not a clean de novo; no inherited origin either
                if not rare(v, rec_max):
                    return "rarity"
                if not (A._str(v, "gene") or A.symbol(v)):
                    # not poolable (no gene) and the de novo row itself did not fire: why?
                    if not emit_denovo:
                        return "mode_disabled"
                    if not G.sample_qc(v, c, thr, "denovo_child"):
                        return "qc_child"
                    if not rare(v, dom_max):
                        return "rarity"
                    if _hiconf_blocks(v):
                        return "hiconf_tag"
                return "other"
            if not rare(v, rec_max):
                return "rarity"
            if not (A._str(v, "gene") or A.symbol(v)):
                return "no_gene"
            if gd == G.HOM_ALT and gmm == G.HOM_ALT:
                return "mendelian_inconsistent"           # 1/1 x 1/1 cannot make a het
            if gd not in (G.HET, G.HOM_ALT) and gmm not in (G.HET, G.HOM_ALT):
                return "parent_nocall"                    # no carrier, and not both hom-ref
            return "other"     # a carrier exists -> pooled with transmitting_parent_qc_fail, never here
        # HOM_ALT child
        if male_x_chrx:
            if gmm == G.UNKNOWN:
                return "parent_nocall"
            if gmm == G.HOM_REF:                          # male-X de novo shape
                if not emit_denovo:
                    return "mode_disabled"
                if not G.sample_qc(v, c, thr, "hom_alt"):
                    return "qc_child"
                if not G.sample_qc(v, m, thr, "clean_parent"):
                    return "qc_parent"
                if not rare(v, dom_max):
                    return "rarity"
                if (G.dp(v, c) or 0) < thr.denovo_min_dp:
                    return "qc_child"
                if _hiconf_blocks(v):
                    return "hiconf_tag"
                return "other"
            if not G.sample_qc(v, c, thr, "hom_alt"):
                return "qc_child"
            if not _carrier_qc(v, m, gmm):
                return "qc_parent"
            if not rare(v, rec_max):
                return "rarity"
            return "other"
        if G.is_x_nonpar(v):                              # hom-alt daughter
            if G.UNKNOWN in (gd, gmm):
                return "parent_nocall"
            if gd != G.HOM_ALT or gmm not in (G.HET, G.HOM_ALT):
                return "mendelian_inconsistent"           # her father must be hemizygous
            if not G.sample_qc(v, c, thr, "hom_alt"):
                return "qc_child"
            if not (_carrier_qc(v, m, gmm) and G.sample_qc(v, d, thr, "hom_alt")):
                return "qc_parent"
            if not rare(v, rec_max):
                return "rarity"
            return "other"
        if G.UNKNOWN in (gd, gmm):                        # autosomal hom-alt
            return "parent_nocall"
        if G.HOM_REF in (gd, gmm):
            return "mendelian_inconsistent"               # deletion-in-trans / UPD shape
        if not G.sample_qc(v, c, thr, "hom_alt"):
            return "qc_child"
        if not (_carrier_qc(v, d, gd) and _carrier_qc(v, m, gmm)):
            return "qc_parent"
        if not rare(v, rec_max):
            return "rarity"
        return "other"

    for v in vcf:
        stats["examined"] += 1
        if A.consequence(v) is not None or A.impact(v) is not None:
            stats["annotated"] += 1       # the Step-4 transfer landed on this record
        if require_pass and v.FILTER:  # cyvcf2 FILTER is None for PASS/'.'
            stats["skipped"]["filter"] += 1
            continue
        # sex unresolved -> ploidy-aware X/Y logic cannot be applied; skip sex chromosomes
        # rather than silently assume female (fail-soft; autosomal modes still run)
        if not gt.sex_known and G.is_sex_nonpar(v):
            stats["skipped"]["sex_unresolved"] += 1
            continue
        # non-PAR chrY only exists in males; a female chrY non-PAR call is an artifact
        if G.is_y_nonpar(v) and not gt.child_male:
            stats["skipped"]["chry_female"] += 1
            continue
        key = f"{v.CHROM}:{v.POS}:{v.REF}:{v.ALT[0]}"
        n_rows_before = len(rows)
        c, d, m = gt.c, gt.d, gt.m
        gc, gd, gmm = v.gt_types[c], v.gt_types[d], v.gt_types[m]
        is_plp = gc in (G.HET, G.HOM_ALT) and A.clnsig_is_plp(v)
        if is_plp:
            _fr = A.frequency(v, cfg)
            if _fr is not None and _fr >= rec_max:
                n_plp_inert += 1   # P/LP the child carries, above every mode's rarity gate

        # male non-PAR chrX/chrY = hemizygous; a het call there is a QC red flag
        male_x = G.is_sex_nonpar(v) and gt.child_male
        # X-SPECIFIC. The hemizygous models below are keyed on the MOTHER's genotype, which is
        # valid ONLY on chrX. On chrY the father is the sole transmitter and the mother has no Y
        # at all, so routing chrY through them (a) reports father-to-son Y transmission as a de
        # novo, since the father is never consulted, and (b) can attribute a Y allele to maternal
        # transmission via x_linked_recessive — and with no gnomAD chrY AF in the VEP cache the
        # rarity gate passes unconditionally, so recurrent Yq/X-transposed mismapping artifacts
        # would land in Step 6's X-linked tier at the floored q. That is a bounded version of the
        # chrM flood 01_make_cohort_sites.sh exists to prevent.
        # `male_x` still guards the het-suppression rule below (hemizygous is true of BOTH
        # chromosomes); `male_x_chrx` guards anything the mother's genotype drives.
        male_x_chrx = G.is_x_nonpar(v) and gt.child_male

        # ---- de novo (SECONDARY / cross-reference only; review handled elsewhere) ----
        denovo_hit = False
        if emit_denovo and not male_x and gc == G.HET and gd == G.HOM_REF and gmm == G.HOM_REF:
            denovo_hit = True
        # male-X de novo: the son's single X comes from the MOTHER; the father transmits Y, so his
        # chrX is irrelevant — require only mother hom-ref, not the father.
        elif emit_denovo and male_x_chrx and gc == G.HOM_ALT and gmm == G.HOM_REF:
            denovo_hit = True
        if denovo_hit:
            child_kind = "denovo_child" if not male_x else "hom_alt"
            # parental cleanliness: both parents for an autosomal de novo; only the transmitting
            # mother for a male-X de novo (father's chrX is not transmitted to a son).
            parents_clean = G.sample_qc(v, m, thr, "clean_parent")
            # A parent with NO allele-depth data passes clean_parent vacuously (the AD limbs fail
            # open while het/hom_alt fail closed — see genotype.sample_qc_ad_measured). Track it so
            # a de novo affirmed by an unmeasured parent is distinguishable from one affirmed by a
            # measured one; a GATK ref-block 0/0 parent has exactly this shape.
            parents_ad = G.sample_qc_ad_measured(v, m, "clean_parent")
            if not male_x:
                parents_clean = parents_clean and G.sample_qc(v, d, thr, "clean_parent")
                parents_ad = parents_ad and G.sample_qc_ad_measured(v, d, "clean_parent")
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
                r = base_row(trio_id, v, gt, "denovo_x_hemi" if male_x else "denovo", cfg=cfg)
                if not parents_ad:
                    r["flags"] = (r["flags"] + ";" if r["flags"] else "") + "parent_ad_unmeasured"
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
                rows.append(tag_strict(base_row(trio_id, v, gt, "hom_recessive", cfg=cfg), v))

        # ---- X-linked recessive, affected male: hemizygous son + carrier mother. The father
        #      transmits his Y (not his X) to a son, so his chrX genotype is IRRELEVANT and is not
        #      required — an affected/carrier father or a father chrX no-call must not drop the call. ----
        if male_x_chrx and gc == G.HOM_ALT and gmm in (G.HET, G.HOM_ALT):
            if G.sample_qc(v, c, thr, "hom_alt") and carrier_ok(m, gmm) and rare(v, rec_max):
                r = base_row(trio_id, v, gt, "x_linked_recessive", cfg=cfg)
                if gd in (G.HET, G.HOM_ALT):
                    r["flags"] = (r["flags"] + ";" if r["flags"] else "") + "father_carries_x_allele"
                rows.append(tag_strict(r, v))

        # ---- X-linked recessive, affected female: HOM_ALT daughter, carrier mother, hemizygous-
        #      affected father (he DOES transmit his X to a daughter) (docs §3.4) ----
        if (not male_x and G.is_x_nonpar(v) and gc == G.HOM_ALT
                and gmm in (G.HET, G.HOM_ALT) and gd == G.HOM_ALT):
            if (G.sample_qc(v, c, thr, "hom_alt") and carrier_ok(m, gmm)
                    and G.sample_qc(v, d, thr, "hom_alt") and rare(v, rec_max)):
                rows.append(tag_strict(base_row(trio_id, v, gt, "x_linked_recessive", cfg=cfg), v))

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
                # Trans-by-descent needs the NON-transmitting parent to be an affirmative,
                # QC-confident hom-ref — the documented rule is "mat 0/1 + pat 0/0". A no-call or
                # an unqualified 0/0 (allele dropout) there is NOT evidence of non-transmission:
                # if that parent silently carries the allele too, a "mat x pat" pair can be CIS.
                # We still emit (never-drop) but mark origin unverified, so an inferred pair is
                # distinguishable from a phase-confirmed one. NB gt_types cannot be trusted to
                # surface this on its own: with cyvcf2's default strict_gt=False a half-called
                # `0/.` parent is reported as HOM_REF, not UNKNOWN — hence strict_gt=True below.
                mom_clear = gmm == G.HOM_REF and G.sample_qc(v, m, thr, "hom_ref")
                dad_clear = gd == G.HOM_REF and G.sample_qc(v, d, thr, "hom_ref")
                # ...and whether that clearance rested on an actual measurement. hom_ref passes
                # vacuously with no AD, and mom_clear/dad_clear are the ONLY evidence a
                # compound-het pair is in TRANS — so a vacuous pass here silently manufactures
                # phase. Never-drop: the call still stands, it is flagged.
                mom_meas = G.sample_qc_ad_measured(v, m, "hom_ref")
                dad_meas = G.sample_qc_ad_measured(v, d, "hom_ref")
                unverified = False
                vacuous = False
                # The TRANSMITTING parent's own genotype QC used to DECIDE the call: a het band
                # or hom-alt band failure set origin to None, and the child's call — GQ 99, AB
                # 0.50 — vanished with no row, while the identical failure in the NON-transmitting
                # parent was emitted with `origin_unverified`. A true-het parent falls outside the
                # 0.25-0.75 band by binomial chance alone ~0.5% of the time at DP 30 and ~3.5% at
                # DP 15 (real WGS tails are heavier), so this deleted a few inherited candidates
                # per trio on ANY callset. Never-drop, applied symmetrically: the call stands and
                # carries `transmitting_parent_qc_fail`, and such a leg cannot veto the dominant
                # model (below), exactly as an unverified or unphased pair cannot. The flag keeps
                # the real ambiguity visible — a parent "het" at AB 0.05 may be a miscall, in
                # which case the child's variant is a de novo the cleanliness test also refused.
                tp_fail = False
                if gmm == G.HOM_REF and gd == G.HOM_REF:
                    # only a genuinely de novo het may pair in trans; require BOTH parents to
                    # pass cleanliness QC, else a dropped-out parental het masquerades as de novo
                    # and gets paired cis with a real variant from that same parent
                    origin = ("denovo" if (G.sample_qc(v, m, thr, "clean_parent")
                                           and G.sample_qc(v, d, thr, "clean_parent")) else None)
                elif mom_carries and dad_carries:
                    # A HOM_ALT parent transmits the alt OBLIGATELY, so origin is deterministic:
                    # a HET child of a 1/1 parent took the alt from that parent and the ref from
                    # the other. Only HET x HET is a genuine 50/50 ("both", never paired below).
                    # Collapsing these to "both" lost phase-CONFIRMED compound hets — most often
                    # on chrX, where a diploid caller renders a hemizygous carrier father as 1/1.
                    # Obligate transmission is decided by the 1/1 parent alone; the other, het
                    # parent's QC has no bearing on where the child's alt came from.
                    if gmm == G.HOM_ALT and gd == G.HET:
                        origin, tp_fail = "mat", not mom_ok
                    elif gd == G.HOM_ALT and gmm == G.HET:
                        origin, tp_fail = "pat", not dad_ok
                    elif gmm == G.HOM_ALT and gd == G.HOM_ALT:
                        origin = None            # 1/1 x 1/1 -> a HET child is a Mendelian error
                    else:
                        origin, tp_fail = "both", not (mom_ok and dad_ok)   # either may have transmitted
                elif mom_carries:
                    origin, tp_fail = "mat", not mom_ok
                    unverified = not dad_clear
                    vacuous = dad_clear and not dad_meas
                elif dad_carries:
                    origin, tp_fail = "pat", not dad_ok
                    unverified = not mom_clear
                    vacuous = mom_clear and not mom_meas
                else:
                    origin = None                # a parent no-call — inheritance unestablished
                if origin:
                    hets.setdefault(gene, []).append((origin, v, key, unverified, vacuous, tp_fail))
                    collected[key] = (v, is_plp, origin)
        # No row from any mode block and not pooled for pairing: this variant is finished, and
        # the only remaining question is why. (A pooled het is resolved after pairing below.)
        if len(rows) == n_rows_before and key not in collected:
            pending[key] = (why_no_row(v, gc, gd, gmm, male_x, male_x_chrx), is_plp)

    # ---- compound het (recessive): trans pairs with determinable parent-of-origin ----
    # A mat×pat pair is trans BY DESCENT (confirmed). A pair with a DE NOVO leg is NOT phase-
    # confirmed: trio genotypes cannot phase a de novo against an inherited variant (it is ~50/50
    # cis/trans), and read-backed phasing (WhatsHap) is not wired here. Such a pair is still emitted
    # — a de novo second hit is biologically valid — but flagged `unphased_denovo_partner` so it is
    # not read as a confirmed biallelic hit. See docs/inheritance_and_genotype_qc.md §3.3.
    consumed = set()
    pair_n = 0
    for gene, cands in hets.items():
        by = {"mat": [], "pat": [], "denovo": [], "both": []}
        for origin, v, key, unver, vac, tpf in cands:
            by[origin].append((v, key, unver, vac, tpf))
        denovo_keys = {k for _, k, _, _, _ in by["denovo"]}
        pairs = [(a, b) for a in by["mat"] for b in by["pat"] + by["denovo"]]
        pairs += [(a, b) for a in by["pat"] for b in by["denovo"]]
        for (va, ka, ua, wa, ta), (vb, kb, ub, wb, tb) in pairs:
            pair_n += 1
            pid = f"{trio_id}:CH{pair_n}"
            unphased = ka in denovo_keys or kb in denovo_keys
            # ONLY a phase-confirmed pair may veto the dominant model — and "confirmed" means the
            # trans evidence was TESTED and PASSED on QC-confident genotypes. An unphased
            # de-novo-partner pair is "a candidate to confirm, not a confirmed biallelic hit"; a
            # pair whose non-transmitting parent was never affirmatively observed hom-ref
            # (origin_unverified) has not established trans either; and a leg whose transmitting
            # parent failed its own QC has not established its origin. Letting any of them consume
            # its legs silently deleted genuine dominant calls from the recurrence tally Step 6
            # headlines. (A VACUOUS pass — a ref-block parent with no AD at all — still consumes:
            # the documented pass stands, flagged trans_evidence_unmeasured.)
            if not unphased and not (ua or ub) and not (ta or tb):
                consumed.add(ka)
                consumed.add(kb)
            pair_flags = []
            if unphased:
                pair_flags.append("unphased_denovo_partner")
            if ua or ub:
                pair_flags.append("origin_unverified")
            if ta or tb:
                pair_flags.append("transmitting_parent_qc_fail")
            # The trans evidence PASSED but rested on a parent with no allele-depth data, so the
            # phase is inferred from a genotype call alone. Distinct from origin_unverified (which
            # means the test FAILED): this one is a vacuous pass, and without the flag it is
            # indistinguishable in the output from a measured one.
            if wa or wb:
                pair_flags.append("trans_evidence_unmeasured")
            for v in (va, vb):
                r = tag_strict(base_row(trio_id, v, gt, "compound_het", pid, cfg=cfg), v)
                for fl in pair_flags:
                    r["flags"] = (r["flags"] + ";" if r["flags"] else "") + fl
                rows.append(r)

    # ---- dominant (inherited het): rare, functional, transmitted from >=1 parent,
    #      not part of a compound-het pair. This is the recurrence signal Step 6 tallies. ----
    if emit_dominant:
        for gene, cands in hets.items():
            for origin, v, key, unver, vac, tpf in cands:
                if origin in ("mat", "pat", "both") and key not in consumed and rare(v, dom_max):
                    r = base_row(trio_id, v, gt, "dominant", cfg=cfg)
                    r["flags"] = f"origin={origin}"
                    if unver:
                        r["flags"] += ";origin_unverified"
                    # the non-transmitting parent PASSED cleanliness but carried no allele-depth
                    # data, so parent-of-origin rests on a genotype call alone
                    if vac:
                        r["flags"] += ";parent_ad_unmeasured"
                    # the transmitting parent's own genotype failed GQ/DP/AB QC: the origin is a
                    # genotype call the QC would not vouch for (see the collection note above)
                    if tpf:
                        r["flags"] += ";transmitting_parent_qc_fail"
                    rows.append(r)

    # ---- accounting: every examined variant is now either skipped, called, or classified ----
    emitted = {f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}" for r in rows}
    stats["with_call"] = len(emitted)
    for key, (why, plp) in pending.items():
        stats["no_row"][why] += 1
        if plp:
            stats["clinvar_plp_no_row"] += 1
    for key, (v, plp, origin) in collected.items():
        if key in emitted:
            continue
        if origin == "denovo":
            # pooled as a possible trans partner, but the de novo row itself did not fire
            why = ("rarity" if not rare(v, dom_max) else "hiconf_tag" if _hiconf_blocks(v)
                   else "mode_disabled" if not emit_denovo else "other")
        else:
            # an inherited het pooled at recessive_max that paired with nothing and sits above
            # dominant_max: the [dominant_max, recessive_max) band, emitted under no mode. On real
            # WGS this is the default fate of most surviving hets — count it, do not hide it.
            why = ("inert_band_het" if not rare(v, dom_max)
                   else "mode_disabled" if not emit_dominant else "other")
        stats["no_row"][why] += 1
        if plp:
            stats["clinvar_plp_no_row"] += 1
    stats["clinvar_plp_ge_recessive_max"] = n_plp_inert
    return rows, stats


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, help="trios.candidates.tsv from Step 4")
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True, help="candidate calls TSV")
    ap.add_argument("--qc-report", default="", help="Step 0 qc_report.tsv (for inferred kid sex)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    problems = validate_filters(cfg)
    if problems:
        sys.stderr.write("ERROR: the screen's filter settings are incoherent:\n"
                         + "".join(f"  - {m}\n" for m in problems))
        return 1
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
    # The dropless info_* block is the UNION over trios (a trio whose source VCF lacked, say,
    # hiConfDeNovo contributes no column for it), kept in first-seen = header order.
    info_cols, _info_seen = [], set()
    n_plp_inert_total = 0
    tot = {"examined": 0, "with_call": 0, "annotated": 0, "clinvar_plp_no_row": 0,
           "skipped": Counter(), "no_row": Counter()}
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
        # strict_gt=True is LOAD-BEARING for every "is this parent a no-call?" test below. With
        # cyvcf2's default (strict_gt=False) a HALF-called genotype like `0/.` is reported as
        # HOM_REF, not UNKNOWN — see cyvcf2 helpers.c as_gts(): "if a single allele is missing
        # e.g 0/. it's still encoded as hom ref because it has no alts". That silently defeats the
        # de novo requirement that both parents be CONFIDENTLY hom-ref, and makes a half-called
        # parent look like an affirmative non-carrier when establishing compound-het trans.
        vcf = VCF(vcf_path, strict_gt=True)
        # WITNESS on the file this step actually READS. Step 3 asserts the faf95 witness on the
        # annotated union; Step 4 then re-derives per-trio genotypes from the RAW trio VCF and
        # re-attaches every annotation with one blanket `bcftools annotate -c INFO`. Nothing checked
        # that the transfer landed here — and a transfer that lifts nothing is not loud: every
        # record reads rarity_basis=absent (= rarest), the candidate list gets BIGGER, the run exits
        # 0, and every row asserts "gnomAD has no record for this allele" for alleles nobody looked
        # up. The header test catches a transfer that never ran; the value-level test after
        # screening catches one that ran and matched nothing (a stale plausible-set index).
        if A.rarity_oracle(cfg) == "faf95":
            witness = A.F["gnomad_af_joint"]
            if f"##INFO=<ID={witness}," not in vcf.raw_header:
                sys.stderr.write(
                    f"ERROR: {trio_id}: resources.gnomad.oracle is 'faf95' but {vcf_path} declares "
                    f"no INFO/{witness} — Step 4's annotation transfer never landed on this per-trio "
                    "VCF, so every variant would read as absent from gnomAD (= rarest) and no "
                    "rarity gate would fire. Re-run Step 4 (rm its .done markers) and check its "
                    "'annotated_genotypes' audit line.\n")
                vcf.close()
                return 1
        try:
            gt = Trio(vcf, ped, thr)
        except KeyError as e:
            sys.stderr.write(f"WARN: {trio_id}: {e}; skipping\n")
            vcf.close()
            continue
        for _c in gt.info_cols:
            if _c not in _info_seen:
                _info_seen.add(_c)
                info_cols.append(_c)
        trio_rows, st = screen_trio(trio_id, vcf, gt, cfg)
        vcf.close()
        if st["examined"] and not st["annotated"]:
            sys.stderr.write(
                f"ERROR: {trio_id}: none of the {st['examined']} candidate records in {vcf_path} "
                "carries a VEP consequence or impact — Step 4's annotation transfer matched nothing "
                "(a stale plausible-set index, or a contig-naming mismatch between the trio VCF and "
                "the plausible set). Every mode would gate on blank annotations. Re-run Step 4.\n")
            return 1
        n_trios += 1
        n_plp_inert = st["clinvar_plp_ge_recessive_max"]
        n_plp_inert_total += n_plp_inert
        for k in ("examined", "with_call", "annotated", "clinvar_plp_no_row"):
            tot[k] += st[k]
        tot["skipped"].update(st["skipped"])
        tot["no_row"].update(st["no_row"])
        # per-trio audit: the INPUT side first (what was examined and what became of it), then
        # the output side (calls by mode). The two reconcile:
        #     variants_examined == sum(skipped.*) + variants_with_call + variants_no_row
        n_no_row = sum(st["no_row"].values())
        audit.record("05_inheritance", "variants_examined", st["examined"], scope=trio_id)
        audit.record("05_inheritance", "variants_with_call", st["with_call"], scope=trio_id)
        audit.record("05_inheritance", "variants_no_row", n_no_row, scope=trio_id)
        for k, n in sorted(st["skipped"].items()):
            audit.record("05_inheritance", f"skipped.{k}", n, scope=trio_id)
        for k, n in sorted(st["no_row"].items()):
            audit.record("05_inheritance", f"no_row.{k}", n, scope=trio_id)
        tmodes = {}
        for row in trio_rows:
            tmodes[row["mode"]] = tmodes.get(row["mode"], 0) + 1
        audit.record("05_inheritance", "candidate_calls", len(trio_rows), scope=trio_id)
        if n_plp_inert:
            audit.record("05_inheritance", "clinvar_plp_dropped_ge_recessive_max", n_plp_inert, scope=trio_id)
        # recorded even when 0: a reviewer must be able to tell "none lost" from "never counted"
        audit.record("05_inheritance", "clinvar_plp_no_row", st["clinvar_plp_no_row"], scope=trio_id)
        for mode, c in sorted(tmodes.items()):
            audit.record("05_inheritance", f"mode.{mode}", c, scope=trio_id)
        why = ", ".join(f"{k}={n}" for k, n in sorted(st["no_row"].items()))
        sys.stderr.write(
            f"  [{trio_id}] {st['examined']} examined -> {st['with_call']} with a call, "
            f"{n_no_row} no row" + (f" ({why})" if why else "")
            + (f"; skipped {dict(st['skipped'])}" if st["skipped"] else "") + "\n")
        all_rows.extend(trio_rows)

    out_cols = COLS + info_cols     # the curated columns first, then the dropless info_* block
    with open(args.out, "w") as out:
        out.write("\t".join(out_cols) + "\n")
        for r in all_rows:
            out.write("\t".join(str(r.get(c, "")) for c in out_cols) + "\n")
    audit.record("05_inheritance", "info_passthrough_columns", len(info_cols))

    by_mode = {}
    for r in all_rows:
        by_mode[r["mode"]] = by_mode.get(r["mode"], 0) + 1
    # The rarity oracle is a RUN-LEVEL fact and belongs in the audit, not only in a column: a
    # methods section has to state which quantity every gate in the run was applied to, and a
    # reader of audit/summary.md should not have to open a TSV to find out.
    audit.record("05_inheritance", f"rarity_oracle.{A.rarity_oracle(cfg)}", 1)
    # JOIN-COVERAGE GUARD for the faf95 arm. The gnomAD joint slim is a SUPERSET of the VEP
    # cache's frequencies (the cache carries only dbSNP-accessioned alleles), so a variant with a
    # cache proxy AF but NO gnomAD record should not exist. When it does, the transfer
    # under-matched — a partially-downloaded slim, a contig-naming mismatch on some chromosomes,
    # or a normalisation difference — and those variants silently read as RAREST, which floods the
    # candidate list with genuinely common alleles. Step 2's 0-match guard catches total failure;
    # this catches PARTIAL failure, which is the more likely and more dangerous case.
    if A.rarity_oracle(cfg) == "faf95":
        n_orphan = sum(1 for r in all_rows
                       if r.get("rarity_basis") == "absent" and r.get("grpmax_af"))
        audit.record("05_inheritance", "rarity_faf95_absent_but_cache_has_af", n_orphan)
        if n_orphan:
            sys.stderr.write(
                f"WARN: {n_orphan} of {len(all_rows)} calls have NO gnomAD joint record yet DO "
                f"carry a VEP-cache gnomAD AF. The joint slim is a superset of the cache, so this "
                f"should be 0 — the Step-2 transfer likely under-matched (partial slim? contig "
                f"naming on some chromosomes?). Those calls read as RAREST and may be common "
                f"alleles. Check Step 2's 'gnomAD joint matched N / M sites' line.\n")
    for _b, _n in sorted(Counter(r.get("rarity_basis", "") for r in all_rows).items()):
        if _b:
            audit.record("05_inheritance", f"rarity_basis.{_b}", _n)
    audit.record("05_inheritance", "trios_screened", n_trios)
    audit.record("05_inheritance", "variants_examined", tot["examined"])
    audit.record("05_inheritance", "variants_with_call", tot["with_call"])
    audit.record("05_inheritance", "variants_no_row", sum(tot["no_row"].values()))
    for k, n in sorted(tot["skipped"].items()):
        audit.record("05_inheritance", f"skipped.{k}", n)
    for k, n in sorted(tot["no_row"].items()):
        audit.record("05_inheritance", f"no_row.{k}", n)
    audit.record("05_inheritance", "candidate_calls_total", len(all_rows))
    audit.record("05_inheritance", "clinvar_plp_dropped_ge_recessive_max", n_plp_inert_total)
    audit.record("05_inheritance", "clinvar_plp_no_row", tot["clinvar_plp_no_row"])
    for mode, c in sorted(by_mode.items()):
        audit.record("05_inheritance", f"mode.{mode}", c)
    sys.stderr.write(
        f"Step 5 complete: {len(all_rows)} candidate calls across {n_trios} trios "
        f"-> {args.out}\n  by mode: {by_mode}\n"
        f"  examined {tot['examined']} variants: {tot['with_call']} with a call, "
        f"{sum(tot['no_row'].values())} no row {dict(sorted(tot['no_row'].items()))}, "
        f"skipped {dict(sorted(tot['skipped'].items()))}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
