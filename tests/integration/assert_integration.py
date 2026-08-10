#!/usr/bin/env python3
"""Assert the mock run produced the expected resolution, funnel, and calls."""
from __future__ import annotations

import argparse
import csv
import os
import sys

from cyvcf2 import VCF

FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)
    return cond


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    a = ap.parse_args(argv)
    W = a.work

    # --- resolution ---
    res = {r["kid"]: r for r in rows(os.path.join(W, "trio_resolution.tsv"))}
    check(res.get("CH_A", {}).get("status", "").startswith("resolved"), "CH_A resolved")
    check(res.get("CH_B", {}).get("status", "").startswith("resolved"), "CH_B resolved (family VCF, extra sib)")
    check(res.get("CH_C", {}).get("status") == "unresolved", "CH_C unresolved")
    check("MO_C" in res.get("CH_C", {}).get("missing_members", ""), "CH_C reports missing mom MO_C")
    manifest = rows(os.path.join(W, "trios.resolved.tsv"))
    check(len(manifest) == 2, f"2 trios in resolved manifest (got {len(manifest)})")

    # --- Step 0: CH_B inferred male + contamination gate ---
    qc = {r["trio_id"]: r for r in rows(os.path.join(W, "qc_report.tsv"))}
    check(qc.get("CH_B", {}).get("inferred_sex") == "1", "CH_B inferred male (chrX)")
    # no selfSM configured -> VCF-only CHARR fallback; mock hom-alt AD is 0 ref -> ~0, unflagged
    check(qc.get("CH_A", {}).get("contam_source") == "charr", "contamination falls back to CHARR")
    check(qc.get("CH_A", {}).get("contam_flag") == "0", "CH_A not flagged contaminated (clean)")
    # CH_B's father carries 6 ref reads at a hom-alt site -> CHARR 0.15 > 0.02 -> trio flagged
    check(qc.get("CH_B", {}).get("contam_flag") == "1", "CH_B flagged contaminated (father CHARR)")
    dadc = float(qc.get("CH_B", {}).get("dad_contam") or 0)
    check(0.1 < dadc < 0.2, f"CH_B dad_contam ~0.15 (got {dadc})")
    check((qc.get("CH_B", {}).get("kid_contam") or "0") in ("", "0"), "CH_B proband CHARR clean")
    # Mendelian-error gate (sample-swap proxy): the de novo in CH_A is a Mendelian violation
    check(int(qc.get("CH_A", {}).get("mie_errors") or 0) >= 1, "CH_A Mendelian error detected (de novo)")
    check(qc.get("CH_A", {}).get("mie_flag") == "1", "CH_A MIE flag raised")

    # --- Step 3: plausible sites keep/drop ---
    plaus = {}
    for v in VCF(os.path.join(W, "plausible.sites.vcf.gz")):
        plaus[(v.CHROM, v.POS)] = v.INFO.get("hprv_keep_reason")
    check(("chr2", 8000) in plaus and plaus[("chr2", 8000)] == "clinvar_plp",
          "ClinVar P/LP (LOW impact) kept via clinvar_plp — the SCREEN gates on CLIN_SIG alone "
          "and is deliberately star-blind: stars arrive from the Step-2 ClinVar transfer and are "
          "a Step-9 RANKING input, never a keep/drop gate (never-drop). A 1-star assertion must "
          "still reach review; it is merely ranked below a 3-star one.")
    check(("chr1", 12000) not in plaus, "BA1-common variant dropped at Step 3")
    check(("chr1", 17000) not in plaus, "non-PASS variant dropped before Step 3")
    check(("chr1", 5000) in plaus, "de novo site retained as plausible")

    # --- VEP-only contract: CADD is the ONLY functional predictor, hence the ONLY way any
    # variant below MODERATE impact can survive. If this regresses the screen silently goes
    # coding-only and every intronic/synonymous candidate vanishes. ---
    check(plaus.get(("chr1", 18000)) == "cadd",
          "deep-intronic MODIFIER kept via CADD (a non-coding keep-path)")
    check(("chr1", 18500) not in plaus,
          "intronic MODIFIER with sub-threshold CADD dropped (the CADD gate really gates)")

    # --- SpliceAI keep-path: the SpliceAI plugin (now wired) is the ONLY signal that reaches a
    # deep-intronic cryptic splice site CADD misses. A MODIFIER with low CADD but DS >= 0.2 must be
    # kept via 'spliceai'; its low-DS control must drop. If this regresses, deep-intronic splice
    # detection goes silently dark again. ---
    check(plaus.get(("chr1", 18700)) == "spliceai",
          "deep-intronic MODIFIER (low CADD) kept via SpliceAI delta score >= threshold")
    check(("chr1", 18800) not in plaus,
          "intronic MODIFIER with sub-threshold SpliceAI + low CADD dropped (the SpliceAI gate gates)")

    # --- The gating-dead predictors must not reappear. REVEL/AlphaMissense/MPC are missense-only
    # scores, and every missense is IMPACT=MODERATE, which selection.py keeps at an earlier
    # branch — so these keep-reasons were unreachable even when the code still had them. This is
    # the audit's falsifiable prediction, enforced: they must be exactly 0, cohort-wide. NB
    # 'spliceai' is deliberately NOT in this list any more: SpliceAI is re-wired (a real data source
    # now) and DOES fire, asserted above at chr1:18700. ---
    for dead in ("revel", "alphamissense", "mpc", "loftee_hc"):
        check(dead not in set(plaus.values()), f"no site kept via '{dead}' (removed / unreachable)")

    # --- Step 5: inheritance calls ---
    calls = rows(os.path.join(W, "candidates.calls.tsv"))

    def has(trio, mode, chrom=None, pos=None, gene=None):
        for r in calls:
            if r["trio_id"] == trio and r["mode"] == mode \
               and (chrom is None or r["chrom"] == chrom) \
               and (pos is None or r["pos"] == str(pos)) \
               and (gene is None or r["symbol"] == gene):
                return True
        return False

    check(has("CH_A", "denovo", "chr1", 5000, "GENE1"), "CH_A de novo GENE1")
    check(has("CH_A", "hom_recessive", "chr1", 8000, "GENE2"), "CH_A hom recessive GENE2")

    # --- MULTIALLELIC trans comp-het (child 1/2). `norm -m-` splits it and, without
    # --keep-sum AD, strips the other ALT's reads -> ref_ad~0 -> allele_balance ~1.0 -> the het
    # band rejects BOTH legs and the pair vanishes silently. This is a MISSED DIAGNOSIS channel:
    # it hits exactly the loci where comp-hets concentrate, and no counter records the loss. ---
    ch2 = [r for r in calls if r["trio_id"] == "CH_A" and r["symbol"] == "GENECH2"]
    check(len([r for r in ch2 if r["mode"] == "compound_het"]) == 2,
          "multiallelic (1/2) child yields a compound_het PAIR in GENECH2 — the --keep-sum AD fix")
    check(len({r["pair_id"] for r in ch2 if r["mode"] == "compound_het"}) == 1,
          "both GENECH2 legs share one pair_id (a genuine trans pair, not two singletons)")
    # ...and the child's allele balance must be the corrected ~0.5, not the 1.0 the bug produced
    for r in ch2:
        ab = float(r["child_ab"]) if r["child_ab"] else 0.0
        check(0.25 <= ab <= 0.75,
              f"GENECH2 leg chr1:{r['pos']} child_ab={ab:.3f} inside the het band (bug gave 1.000)")
    # scoped to GENE3: CH_A now also has a legitimate multiallelic comp-het in GENECH2 (below),
    # so an unscoped count of every CH_A compound_het row would be 4, not 2.
    ch = [r for r in calls if r["trio_id"] == "CH_A" and r["mode"] == "compound_het"
          and r["symbol"] == "GENE3"]
    check(len(ch) == 2, f"CH_A compound het pair in GENE3 (got {len(ch)})")
    # comp-het requires TRANS: two cis (both maternal) hets in GENEC must NOT be a compound_het
    genec = [r for r in calls if r["trio_id"] == "CH_A" and r["symbol"] == "GENEC"]
    check(genec and all(r["mode"] != "compound_het" for r in genec),
          "cis (same-parent) GENEC pair NOT called compound_het")
    check(any(r["mode"] == "dominant" for r in genec), "cis GENEC variants emitted as dominant instead")
    # --- A HOM-ALT transmitting parent makes parent-of-origin DETERMINISTIC (a 1/1 parent transmits
    # the alt obligately, so a HET child took the alt from it and the ref from the other) — it is NOT
    # the genuine 50/50 "both" case. Collapsing it to "both" barred the pair from trans-pairing, and
    # since both legs sit at 3e-3 (inside the recessive band, above the 1e-4 dominant gate) neither
    # could fall through to a dominant call either — a phase-CONFIRMED biallelic hit vanished under
    # NO mode at all. This is the regression test for that. ---
    homalt = [r for r in calls if r["trio_id"] == "CH_A" and r["symbol"] == "GENEHOMALT"]
    check(len([r for r in homalt if r["mode"] == "compound_het"]) == 2,
          "HOM-ALT-parent comp-het: obligate maternal transmission yields a trans PAIR in GENEHOMALT")
    check(len({r["pair_id"] for r in homalt if r["mode"] == "compound_het"}) == 1,
          "both GENEHOMALT legs share one pair_id (one genuine trans pair)")
    check(all("origin_unverified" not in (r.get("flags") or "") for r in homalt),
          "GENEHOMALT pair is phase-CONFIRMED (both parents observed) — no origin_unverified flag")
    check(has("CH_B", "denovo", "chr1", 5100, "GENE1"), "CH_B de novo GENE1 (secondary)")
    check(has("CH_B", "x_linked_recessive", "chrX", 2781600, "GENEX"), "CH_B X-linked recessive GENEX")
    # X-linked recessive must fire even with an AFFECTED (hom-alt) father (father's chrX not
    # transmitted to a son) — the previously-required father-hom-ref would have wrongly dropped it
    check(has("CH_B", "x_linked_recessive", "chrX", 2782000, "GENEXAF"),
          "X-linked recessive called with an affected (hom-alt) father")
    # autosomal hom-recessive with a HOM-ALT parent (carrier rule accepts HET or HOM_ALT parents)
    check(has("CH_A", "hom_recessive", "chr1", 8500, "GENE2H"), "hom recessive called with a HOM-ALT parent")
    check(not has("CH_A", "denovo", "chr1", 15000), "low-GQ pseudo-de-novo NOT called (QC gate)")
    # dominant model: rare functional inherited het, recurrent across individuals
    check(has("CH_A", "dominant", "chr2", 10000, "GENED"), "CH_A dominant inherited het GENED")
    check(has("CH_B", "dominant", "chr2", 10000, "GENED"), "CH_B dominant inherited het GENED")
    # THE MAX_AF TRAP, and it must hold on EITHER oracle. GENEFND is at AF 0.002 in gnomAD 'ami'
    # — a bottlenecked founder group (AN~900) excluded from grpmax AND from the FAF group set, so
    # neither arm sees it: the proxy skips ami, and gnomAD publishes no fafmax because no
    # FAF-eligible group carries the allele. MAX_AF nonetheless reports 0.002, 20x over
    # dominant_max=1e-4, and reading it would silently DROP this call.
    check(has("CH_A", "dominant", "chr2", 17000, "GENEFND"),
          "founder-population-only allele (MAX_AF=0.002 in 'ami') still called dominant — "
          "neither oracle reads bottlenecked groups, and neither reads MAX_AF")
    # The `mid` case is DIFFERENT and is the one place the arms legitimately disagree: grpmax
    # excludes mid, the FAF set includes it. On the default (faf95) this allele IS gated.
    check(not has("CH_A", "dominant", "chr2", 17500, "GENEMID"),
          "a mid-enriched allele IS gated under the default faf95 oracle (the FAF group set "
          "includes mid, unlike grpmax) — the one documented difference between the arms")
    fnd = [r for r in calls if r["symbol"] == "GENEFND"]
    check(fnd and all(not r["grpmax_af"] for r in fnd),
          "GENEFND reports an EMPTY grpmax_af (no eligible group carries it)")
    check(fnd and all(r["max_af"] for r in fnd),
          "GENEFND still reports max_af for the reviewer to see why it looked common")
    # CADD-only keep survives into an actual call, not just Step 3
    check(has("CH_A", "dominant", "chr1", 18000, "GENEIN"),
          "deep-intronic CADD-kept variant becomes a dominant call")
    # SpliceAI-kept variant becomes a call AND carries its delta score into candidates.calls.tsv
    check(has("CH_A", "dominant", "chr1", 18700, "GENESAI"),
          "deep-intronic SpliceAI-kept variant becomes a dominant call")
    sai = [r for r in calls if r["symbol"] == "GENESAI"]
    check(sai and all(r.get("spliceai_ds") == "0.55" for r in sai),
          "the SpliceAI delta score (0.55) flows into candidates.calls.tsv for the curator")

    # --- Step 6: recurrence-based gene consolidation ---
    genes = {r["gene"]: r for r in rows(os.path.join(W, "genes.ranked.tsv"))}
    check(genes.get("GENED", {}).get("n_dominant") == "2", "GENED has 2 dominant carriers")
    check(genes.get("GENED", {}).get("recurrent") == "1", "GENED flagged recurrent")
    # same- vs distinct-variant recurrence: GENED shares one variant (founder/artifact-suspect);
    # GENEDD has two distinct variants across trios (the stronger gene signal)
    check(genes.get("GENED", {}).get("recurrence_kind") == "same_variant", "GENED = same-variant recurrence")
    check(genes.get("GENEDD", {}).get("recurrence_kind") == "distinct_variant", "GENEDD = distinct-variant recurrence")
    check(genes.get("GENEDD", {}).get("n_dominant") == "2", "GENEDD has 2 dominant carriers")
    check(genes.get("GENE1", {}).get("n_denovo") == "2", "GENE1 has 2 de novo carriers (secondary)")
    # calibrated recurrence null: a rare variant recurring in 2 individuals is significant
    check(float(genes.get("GENED", {}).get("p_recurrence") or 1) < 1e-4,
          "GENED has a small calibrated recurrence p-value")
    check(genes.get("GENED", {}).get("recurrence_exome_wide_sig") == "1",
          "GENED recurrence is exome-wide significant")
    # GENE1 (de novo only) must NOT get an inherited recurrence p-value
    check(not genes.get("GENE1", {}).get("p_recurrence"),
          "GENE1 (de novo only) has no inherited recurrence p-value")

    # --- chrM is OUT OF SCOPE and must never reach any output. The mock carries a near-fixed
    # rCRS haplogroup variant (m.8860A>G, whole trio hom-alt) in BOTH trios. Un-excluded it fires
    # hom_recessive everywhere and — with no gnomAD mito AF to fail the rarity gate — floors q in
    # Step 6 and lands in the recurrent, exome-wide-significant tier above real nuclear genes. ---
    check(not any(r["chrom"] in ("chrM", "chrMT", "M", "MT") for r in calls),
          "no chrM call in candidates.calls.tsv (chrM excluded at Step 1)")
    check("MT-ATP6" not in genes,
          "MT-ATP6 absent from genes.ranked.tsv — a haplogroup variant never reaches the "
          "recurrent tier")
    for v in VCF(os.path.join(W, "cohort.sites.vcf.gz")):
        if v.CHROM in ("chrM", "chrMT", "M", "MT"):
            check(False, f"chrM leaked into the cohort union at {v.CHROM}:{v.POS}")
            break
    else:
        check(True, "cohort.sites.vcf.gz contains no chrM records")

    # --- GATK's de novo tags must SURVIVE Step 4. A blanket `annotate -x INFO` stripped them,
    # which silently made Step 5's has_hiconf permanently False and filters.denovo.use_hiconf_tag
    # a no-op — invisible, because de novo is still called when the tag is simply absent. ---
    trio_vcfs = {r["trio_id"]: r["candidates_vcf"] for r in rows(os.path.join(W, "trios.candidates.tsv"))}
    if check("CH_A" in trio_vcfs, "CH_A candidate VCF in the manifest"):
        hdr = VCF(trio_vcfs["CH_A"]).raw_header
        check("ID=hiConfDeNovo" in hdr,
              "hiConfDeNovo header SURVIVES Step 4's INFO strip (else the de novo gate is dead code)")
        check(any(v.INFO.get("hiConfDeNovo") for v in VCF(trio_vcfs["CH_A"])),
              "at least one candidate record still carries a hiConfDeNovo value")
        check(any(str(k).startswith("vep_") for v in VCF(trio_vcfs["CH_A"]) for k, _ in v.INFO),
              "vep_* annotations still transfer into the per-trio VCF")

    # --- Step 9: prioritization (gene excess + artifact panel + variant tiering) ---
    vpath9 = os.path.join(W, "variants.prioritized.tsv")
    gpath9 = os.path.join(W, "genes.prioritized.tsv")
    if check(os.path.exists(vpath9), "variants.prioritized.tsv written") and \
            check(os.path.exists(gpath9), "genes.prioritized.tsv written"):
        pv = rows(vpath9)
        pg = {r["gene"]: r for r in rows(gpath9)}
        src = rows(os.path.join(W, "igv", "variants.tsv"))

        # THE NEVER-DROP INVARIANT. Step 9 re-ranks; it must never remove a call. If this
        # regresses a reviewer silently receives a shortened list with no counter recording it.
        check(len(pv) == len(src),
              f"never-drop: prioritized rows == input rows ({len(src)}; got {len(pv)})")
        check({(r["chrom"], r["pos"], r["ref"], r["alt"], r["trio_id"]) for r in pv} ==
              {(r["chrom"], r["pos"], r["ref"], r["alt"], r["trio_id"]) for r in src},
              "never-drop: the prioritized set is exactly the input set (no substitutions)")

        # EVERY scoring term ships as its own column — the transparency requirement, so a
        # reviewer can read WHY a call ranked where it did rather than trusting one number.
        for col in ("pts_molecular", "pts_rarity", "pts_gene_constraint", "pts_recurrence",
                    "pts_quality", "pts_clinical", "pts_moi", "pts_gene_artifact",
                    "pts_gene_list_prior", "priority_points_agnostic", "rank_agnostic",
                    "priority_points_prior", "rank_prior", "rank_delta", "cap_applied",
                    "variant_tier", "variant_tier_reason", "nhf_status", "gene_tier",
                    "downweight_reason", "excess_ratio", "corroboration_count"):
            check(col in pv[0], f"variants.prioritized.tsv has '{col}' column")
        # ...and the terms must SUM to the reported total, or the columns are decoration
        terms = ("pts_molecular", "pts_rarity", "pts_gene_constraint", "pts_recurrence",
                 "pts_quality", "pts_clinical", "pts_moi", "pts_gene_artifact")
        bad_sum = [r for r in pv
                   if r["cap_applied"] == "none"
                   and abs(sum(float(r[t] or 0) for t in terms)
                           - float(r["priority_points_agnostic"])) > 1e-6]
        check(not bad_sum, f"per-term columns sum to priority_points_agnostic ({len(bad_sum)} bad)")

        # Ranks are a permutation of 1..N in BOTH rankings
        check(sorted(int(r["rank_agnostic"]) for r in pv) == list(range(1, len(pv) + 1)),
              "rank_agnostic is a permutation of 1..N")
        # THE PHENOTYPE-AGNOSTIC CONTRACT: with the Class-B overlay disabled (the mock ships the
        # file but leaves enabled: false) the two rankings must be IDENTICAL. This is what keeps
        # hprv phenotype-agnostic by default, so it is asserted rather than assumed.
        check(all(r["rank_prior"] == r["rank_agnostic"] for r in pv),
              "overlay OFF: rank_prior == rank_agnostic for every row")
        check(all(float(r["pts_gene_list_prior"] or 0) == 0.0 for r in pv),
              "overlay OFF: the gene-list prior contributes 0 to every row")
        check(all(r["rank_delta"] == "0" for r in pv), "overlay OFF: rank_delta is 0 everywhere")

        # --- The ARTIFACT LOCUS is down-weighted, with a reason string naming the mechanisms.
        # OR4Q3 carries a tiny mutational target against 4 candidate rows plus every corroborating
        # signal, so it must reach a down-weight tier — and its penalty must be a real negative
        # number, not a flag nobody scores. ---
        art = pg.get("OR4Q3", {})
        check(art.get("gene_tier") in ("T2_downweight", "T3_strong_downweight"),
              f"artifact locus OR4Q3 down-weighted (got {art.get('gene_tier')})")
        check(float(art.get("gene_artifact_penalty") or 0) < 0.0,
              "OR4Q3 carries a negative artifact penalty")
        check(int(art.get("corroboration_count") or 0) >= 2,
              f"OR4Q3 has >=2 corroborating signals (got {art.get('corroboration_count')})")
        for frag in ("excess_ratio=", "artifact_gene_family", "oe_syn=",
                     "gnomad_constraint_flag=mis_too_many"):
            check(frag in (art.get("downweight_reason") or ""),
                  f"OR4Q3 reason string reports '{frag}'")
        # ...and NONE of its variants vanished — a down-weight is a re-rank, not a veto
        check(len([r for r in pv if r["gene"] == "OR4Q3"]) ==
              len([r for r in src if r["gene"] == "OR4Q3"]),
              "every OR4Q3 variant survives the down-weight")

        # --- THE POSITIVE-CONTROL GUARD. GENE1 has the SAME extreme excess shape and the same
        # corroborating signals as OR4Q3 — on the statistics alone it earns T3 — but it is in the
        # established-gene union, so the AUDITABLE ceiling must cap it at T1_watch and flag it for
        # review instead. This is the single most important assertion in Step 9: if it regresses,
        # real predisposition genes get their variants silently penalised. ---
        ctl = pg.get("GENE1", {})
        check(ctl.get("gene_tier") in ("T0_no_downweight", "T1_watch"),
              f"POSITIVE CONTROL: established gene GENE1 never reaches T2/T3 "
              f"(got {ctl.get('gene_tier')})")
        check(ctl.get("established_gene_control") == "1", "GENE1 is in the control union")
        check(ctl.get("control_ceiling_applied") == "1",
              "GENE1's tier was capped BY THE CEILING (so the ceiling, not a weak rule, "
              "is what protected it)")
        check(ctl.get("review_flag") == "established_gene_high_excess",
              "GENE1 flagged established_gene_high_excess — the exemption protects it from a "
              "score penalty, it does NOT mean its calls are correct")
        check(float(ctl.get("gene_artifact_penalty") or 0) >= -0.5,
              "GENE1's penalty is capped at the T1_watch magnitude")

        # --- CDS-FALLBACK offset + its ceiling: a gene with no gnomAD mu gets the CDS-length
        # regression, says so, and never reaches T3 (a +/-30% offset cannot support that claim). ---
        nomu = pg.get("GENENOMU", {})
        check(nomu.get("E_source") == "cds_fallback",
              f"GENENOMU uses the CDS-length fallback offset (got {nomu.get('E_source')})")
        check(nomu.get("gene_tier") != "T3_strong_downweight",
              "a CDS-fallback gene never reaches the strong down-weight tier")
        check("cds_length_fallback" in (nomu.get("downweight_reason") or ""),
              "GENENOMU's reason string discloses that its offset is a fallback")
        # every gene that produced a call is in the gene table, including the no-offset ones
        check({r["gene"] for r in pv} == set(pg), "genes.prioritized.tsv covers every called gene")

        # --- MECHANISM GATING: a molecularly-benign PREDICTION in a highly constrained gene
        # (GENE1: pLI 0.98, LOEUF 0.20) carrying a ClinVar P/LP assertion must still be CAPPED.
        # Even +4 of clinical evidence and a constrained gene cannot lift a V0 — that is the SVI
        # principle, and it is how a gene-list prior is prevented from becoming confirmation bias. ---
        # NB the fixture is emitted under MORE THAN ONE mode: the extra GENE1 rows give the
        # proband a second sub-1e-2 functional het in the same gene, so Step 5 emits the locus
        # both as `dominant` and as a `compound_het` leg (audit A-6's single-gene-keyed mode
        # assignment). Assert over ALL of them, and pick the DOMINANT one for the constraint
        # check — a recessive mode zeroes constraint on its own, which would let the V0 gate
        # pass for the wrong reason.
        v0 = [r for r in pv if r["chrom"] == "chr2" and r["pos"] == "18800"]
        if check(len(v0) >= 1, f"the V0 benign-prediction fixture reached Step 9 ({len(v0)} rows)"):
            check(all(r["variant_tier"] == "V0" for r in v0),
                  f"benign prediction scored V0 ({[r['variant_tier'] for r in v0]})")
            check(all(r["cap_applied"] == "V0_benign" for r in v0),
                  f"V0 cap applied ({[r['cap_applied'] for r in v0]})")
            check(all(float(r["priority_points_agnostic"]) <= 0.0 for r in v0),
                  "the V0 cap holds the total at <=0 even with ClinVar P/LP (+4) and a "
                  f"constrained gene ({[r['priority_points_agnostic'] for r in v0]})")
            check(all(float(r["pts_clinical"]) == 4.0 for r in v0),
                  "the ClinVar term still SCORES +4 — the cap works on the total, so the "
                  "evidence stays visible in its own column rather than being erased")
            dom0 = [r for r in v0 if r["inheritance"] == "dominant"]
            if check(bool(dom0), "the V0 fixture has a dominant-mode row (constraint not "
                                 "already zeroed by a recessive mode)"):
                check(float(dom0[0]["pts_gene_constraint"]) == 0.0,
                      "constraint is ZEROED for a V0 BY THE MECHANISM GATE — a constrained "
                      "gene cannot rescue a benign prediction")
                check(float(dom0[0]["constraint_gate"]) == 0.0,
                      "the V0 mechanism gate itself reads 0.0")
        # ...while a genuine pLoF in the same gene DOES earn the constraint term
        hi = [r for r in pv if r["gene"] == "GENE1" and r["impact"] == "HIGH"]
        check(hi and all(r["variant_tier"] == "V4" for r in hi),
              "a HIGH-impact pLoF scores V4 (V5 is unreachable — no NMD annotation)")
        check(hi and any(float(r["pts_gene_constraint"]) > 0.0 for r in hi),
              "the same constrained gene DOES earn the constraint term for a credible effect")
        check(all(r["nmd_status"] == "INDETERMINATE" for r in pv if r["impact"] == "HIGH"),
              "every pLoF carries nmd_status=INDETERMINATE (no exon/CDS columns exist)")
        check(all(r["variant_tier"] != "V5" for r in pv),
              "no variant reaches V5 — unreachable until variants.tsv carries EXON/CDS_position")
        check(all(r["plof_confidence"] == "UNAVAILABLE" for r in pv),
              "plof_confidence is UNAVAILABLE everywhere (no LOFTEE under this contract)")

        # --- BLANK NHF IS NOT 0.0. This mock has no kraken2 DB, so Step 8b never ran and every
        # NHF cell is blank — which must read as NOT SCREENED, scoring neither a penalty nor
        # credit. Reading blank as "clean" would silently promote exactly the calls nobody
        # examined, and it is the same off-by-a-semantic trap as the Step-8b pos-1 join. ---
        check(all(r["nhf_status"] == "not_screened" for r in pv),
              "no kraken2 DB => every call reads nhf_status=not_screened (blank != clean)")
        # An unscreened call must not be penalised FOR BEING UNSCREENED. Scoped to rows with no
        # OTHER quality term in play: a comp-het legitimately carries -0.5 for
        # partner_leg_quality_unknown (parental GQ/DP/AB are absent from variants.tsv), and a
        # QC-failing call carries -2. Conflating those with the NHF term would make this
        # assertion pass or fail for the wrong reason.
        no_other = [r for r in pv if r["gt_qc_pass"] == "1"
                    and r["partner_leg_quality_unknown"] != "1"]
        check(bool(no_other), "some rows have no competing quality term")
        check(all(float(r["pts_quality"] or 0) == 0.0 for r in no_other),
              "an unscreened call is not PENALISED for being unscreened (nhf_not_screened = 0)")
        # ...and the comp-het rows DO carry the partner term, so its absence above is scoping,
        # not a term that silently never fires
        comphet = [r for r in pv if r["partner_leg_quality_unknown"] == "1"]
        check(comphet and all(abs(float(r["pts_quality"]) + 0.5) < 1e-9 for r in comphet),
              "a comp-het carries -0.5 for the unassessable trans leg (parental GQ/DP/AB absent)")
        check(all(r["nhf_max_fraction"] == "" for r in pv),
              "an unscreened call reports a BLANK max NHF fraction, never 0.0")

        # --- The SpliceAI-kept deep-intronic variant must be tiered on its delta score, and the
        # off-label CADD missense route must be labelled as such wherever it fires. ---
        sai9 = [r for r in pv if r["gene"] == "GENESAI"]
        check(sai9 and all(r["variant_tier"] in ("V3", "V4") for r in sai9),
              "the SpliceAI-kept variant (DS 0.55) is tiered on its splice evidence")
        check(sai9 and all(r["spliceai_status"] == "scored" for r in sai9),
              "a scored SpliceAI variant reports spliceai_status=scored")
        unscored = [r for r in pv if r["spliceai_ds"] == ""]
        check(all(r["spliceai_status"] == "not_covered" for r in unscored),
              "an UNSCORED SpliceAI variant reports not_covered — absence is not 'no effect'")
        cadd_mis = [r for r in pv if r["missense_evidence_source"] == "cadd_offlabel"]
        check(all("off-label" in r["variant_tier_reason"] for r in cadd_mis),
              "a CADD-based missense tier is labelled off-label (never presentable as PP3)")

        # --- calibrated missense predictors: PRECEDENCE, not a max over what is available.
        # ClinGen SVI's rule is to commit to ONE predictor chosen before seeing results, so the
        # ladder is revel -> alphamissense -> cadd(off-label) -> none and always reports which
        # one spoke. The mock's GENE2 row carries REVEL 0.85 against a deliberately LOW cadd=3:
        # if the ladder ever regressed to "take the best score available", cadd would win here
        # and the source would read cadd_offlabel. ---
        rev = [r for r in pv if r["gene"] == "GENE2" and r["revel"]]
        check(rev, "the REVEL-scored missense variant reached Step 9")
        check(all(r["missense_evidence_source"] == "revel" for r in rev),
              "REVEL outranks a (low) CADD — precedence, not best-of-N")
        check(all(r["variant_tier"] == "V4" for r in rev),
              "REVEL >= 0.773 (Pejaver moderate) reaches V4, above the supporting-only rungs")
        am = [r for r in pv if r["gene"] == "GENE3" and r["alphamissense"] and not r["revel"]]
        check(am, "the AlphaMissense-only missense variant reached Step 9")
        check(all(r["missense_evidence_source"] == "alphamissense" for r in am),
              "with REVEL absent the ladder falls through to AlphaMissense, not to CADD")
        check(all(r["alphamissense_class"] in ("likely_pathogenic", "ambiguous", "likely_benign")
                  for r in am),
              "am_class rides along (the PLUGIN's key name, not dbNSFP's AlphaMissense_score)")

        # --- rarity: ONE oracle per run, both arms proven on the real transferred VCF. ---
        faf = [r for r in pv if r["faf95"]]
        check(faf, "the gnomAD joint transfer reached Step 9 (faf95 column populated)")
        check(all(float(r["faf95"]) < float(r["grpmax_af"]) for r in faf
                  if r["grpmax_af"] and float(r["grpmax_af"]) > 0),
              "faf95 (a CI LOWER bound) sits below the point estimate, as it must")
        check(all(r["faf95_group"] in ("afr", "amr", "eas", "mid", "nfe", "sas") for r in faf),
              "faf95_group names a FAF-eligible ancestry group (GRPMAX_POPS + mid)")
        # ONE ORACLE PER RUN. rarity_oracle is a run-level constant (the mock selects faf95);
        # rarity_basis is the per-variant provenance WITHIN that oracle. The arms never cross —
        # an earlier design blended them per variant, which made two rows in one run comparable
        # on different quantities and got the fallback direction wrong for singletons.
        oracles = {r["rarity_oracle"] for r in pv if r["rarity_oracle"]}
        check(oracles == {"faf95"},
              f"exactly ONE rarity oracle for the whole run (got {sorted(oracles)})")
        # BOTH ARMS, against the REAL transferred per-trio VCF. The pipeline ran on the default
        # (proxy) oracle above; this proves the faf95 arm on the same bcftools-transferred data
        # without a second full run — and proves the arms never cross.
        import glob as _glob
        from cyvcf2 import VCF as _VCF
        from hprv import annotations as _A
        FAF = {"resources": {"gnomad": {"oracle": "faf95"}}}
        PRX = {"resources": {"gnomad": {"oracle": "grpmax_proxy"}}}
        seen = {"measured": 0, "zero_ci": 0, "absent": 0}
        crossed = []
        for _tv in _glob.glob(os.path.join(W, "trios", "*.candidates.annotated.vcf.gz")):
            for _v in _VCF(_tv):
                b = _A.rarity_basis(_v, FAF)
                seen[b] = seen.get(b, 0) + 1
                f, px = _A.faf95(_v), _A.grpmax_af(_v)
                # the faf95 arm must never return the proxy's value, and vice versa
                if f is None and px is not None and _A.frequency(_v, FAF) == px:
                    crossed.append(("faf95 arm returned the proxy", _v.CHROM, _v.POS))
                if f is not None and px is not None and f != px and _A.frequency(_v, PRX) == f:
                    crossed.append(("proxy arm returned faf95", _v.CHROM, _v.POS))
                if b == "zero_ci":
                    if _A.frequency(_v, FAF) != 0.0:
                        crossed.append(("zero_ci did not resolve to 0", _v.CHROM, _v.POS))
                if b == "measured" and _A.frequency(_v, FAF) != f:
                    crossed.append(("measured != faf95", _v.CHROM, _v.POS))
        check(not crossed, f"the two oracle arms NEVER cross on real data ({crossed[:2]})")
        check(seen["measured"] > 0, "the faf95 arm has measured rows on the real transferred VCF")
        check(seen["zero_ci"] > 0,
              "the zero_ci basis is exercised (gnomAD has the allele, published no faf95)")
        check(all(r["rarity_af"] == "" and r["rarity_strength"] == "unknown"
                  for r in pv if r["rarity_basis"] == "absent"),
              "absent reads unknown — never a measured zero")

        # nhomalt: the recessive false-positive tell, reported and costing nothing by default
        nh = [r for r in pv if r["nhomalt_recessive_conflict"] == "1"]
        check(nh, "a biallelic call with gnomAD homozygotes raises nhomalt_recessive_conflict")
        check(all(r["inheritance"] in ("compound_het", "hom_recessive", "x_linked_recessive")
                  for r in nh),
              "the nhomalt conflict fires ONLY on biallelic modes")

        # --- ClinVar stars: the mock configures NO ClinVar VCF, so every row must read
        # UNAVAILABLE — and crucially NOT 0. Blank/absent means nobody looked; 0 means ClinVar
        # has a record whose submitter provided no assertion criteria. If a future change ever
        # collapses the two, every P/LP assertion in an un-transferred run gets silently damped. ---
        check(all(r["clinvar_review_status"] == "UNAVAILABLE" for r in pv),
              "with no ClinVar VCF configured, review status reads UNAVAILABLE for every call")
        check(all(r["clinvar_stars"] == "" for r in pv),
              "an absent ClinVar transfer leaves clinvar_stars BLANK, never 0")

        # --- MOI coherence: curated genes get a verdict, UNCURATED genes are EXACTLY neutral.
        # Any penalty on moi_unknown converts the score into a known-gene filter and destroys
        # novel-gene discovery. ---
        unk = [r for r in pv if r["moi_coherence"] == "unknown"]
        check(unk, "some genes are uncurated (moi_unknown) — the novel-gene case exists here")
        check(all(float(r["pts_moi"] or 0) == 0.0 for r in unk),
              "moi_unknown scores EXACTLY 0 — never a penalty (that would filter novel genes)")
        gened = [r for r in pv if r["gene"] == "GENED"]
        check(gened and all(r["moi_coherence"] == "coherent" for r in gened),
              "GENED (curated AD, dominant calls) reads moi_coherent")

        # --- Recurrence scores on the CARRIER COUNT, not on the saturating case-only p-value,
        # and same-variant recurrence gets strictly less credit than distinct-variant. ---
        gd = [r for r in pv if r["gene"] == "GENED"]        # same-variant across 2 trios
        gdd = [r for r in pv if r["gene"] == "GENEDD"]      # distinct variants across 2 trios
        check(gd and all(float(r["pts_recurrence"]) == 0.5 for r in gd),
              "same-variant recurrence earns the reduced credit (+0.5)")
        check(gdd and all(float(r["pts_recurrence"]) == 1.0 for r in gdd),
              "distinct-variant recurrence earns full 2-carrier credit (+1.0)")

        # --- The CALIBRATION DIAGNOSTIC (the A-3 gap) is recorded, and it is recorded for BOTH
        # nulls so the NB-vs-Poisson choice is auditable rather than asserted. ---
        acount = {(r["step"], r["metric"]): r["value"]
                  for r in rows(os.path.join(W, "audit", "counts.tsv"))}
        for m in ("variants_in", "variants_out", "null_model", "null_C", "null_alpha",
                  "null_trim_iterations", "null_phi_bulk", "null_phi_all_genes",
                  "calibration.nb_mean_midp", "calibration.poisson_mean_midp",
                  "genes_downweighted_T2_T3", "variants_downweighted_T2_T3",
                  "control_ceiling_applied", "gene_tier.T0_no_downweight"):
            check(("09_prioritize", m) in acount, f"audit records 09_prioritize/{m}")
        check(acount.get(("09_prioritize", "variants_in")) ==
              acount.get(("09_prioritize", "variants_out")),
              "the audit itself records never-drop (variants_in == variants_out)")
        check(int(acount.get(("09_prioritize", "control_ceiling_applied")) or 0) >= 1,
              "audit records that the established-gene ceiling actually fired")
        check(acount.get(("09_prioritize", "null_model")) == "negative_binomial",
              "the null is the negative binomial (Poisson is 2.5x anti-conservative)")
        # phi_all >> phi_bulk is the signature that justifies BOTH the NB and the trim
        check(float(acount.get(("09_prioritize", "null_phi_all_genes")) or 0) >
              float(acount.get(("09_prioritize", "null_phi_bulk")) or 0),
              "raw dispersion exceeds the trimmed bulk dispersion (overdispersion is in the tail)")

        # --- IDEMPOTENCY: the .done marker means a re-run is a no-op. ---
        check(os.path.exists(vpath9 + ".done"), "Step 9 wrote its .done marker (idempotent)")

    # --- audit exists ---
    check(os.path.exists(os.path.join(W, "audit", "summary.md")), "audit/summary.md written")
    counts = rows(os.path.join(W, "audit", "counts.tsv"))
    check(any(r["step"] == "resolve" and r["metric"] == "trios_resolved" and r["value"] == "2"
              for r in counts), "audit records 2 resolved trios")

    # --- Step 7: xlsx summary ---
    xlsx = os.path.join(W, "hprv_summary.xlsx")
    if check(os.path.exists(xlsx), "xlsx summary written"):
        from openpyxl import load_workbook
        wb = load_workbook(xlsx, read_only=True)
        for sh in ("About", "Gene consolidation", "Candidate calls"):
            check(sh in wb.sheetnames, f"xlsx has '{sh}' sheet")

    # --- Step 8: igv.js variant-review export ---
    vpath = os.path.join(W, "igv", "variants.tsv")
    if check(os.path.exists(vpath), "igv variants.tsv written"):
        vh, vrows = None, []
        with open(vpath) as fh:
            rr = list(csv.reader(fh, delimiter="\t"))
        vh, vrows = rr[0], rr[1:]
        for col in ("chrom", "pos", "ref", "alt", "inheritance", "child_file", "child_gt"):
            check(col in vh, f"variants.tsv has '{col}' column")
        check(len(vrows) == len(calls), f"variants.tsv rows == candidate calls ({len(calls)})")
        fi = vh.index("child_file")
        check(any(r[fi] for r in vrows), "at least one child_file (mini-CRAM) populated")
    check(os.path.exists(os.path.join(W, "igv", "crams", "CH_A", "CH_A.cram")),
          "CH_A mini-CRAM extracted")
    check(os.path.exists(os.path.join(W, "igv", "trios.tsv")), "igv trios.tsv written")
    check(os.path.exists(os.path.join(W, "igv", "curation.json")), "igv curation.json written")

    if FAILS:
        sys.stderr.write(f"\n{len(FAILS)} assertion(s) FAILED\n")
        return 1
    print("\nALL INTEGRATION ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
