#!/usr/bin/env python3
"""Generate a tiny, self-consistent mock dataset exercising the whole pipeline.

Writes under --out:
  reference.fa                mini GRCh38-like genome (chr1, chr2, chrX-with-nonPAR)
  vcfs/fileA.vcf              trio A (CH_A/FA_A/MO_A) — autosomal modes + filter cases
  vcfs/fileB.vcf              FAMILY VCF: [MO_B, SIB_B, CH_B, FA_B] (extra sibling,
                              shuffled order) — de novo (recurrent gene) + X-linked
  vcfs/fileC.vcf              duo CH_C/FA_C (mom absent) — resolver "unresolved" case
  trios.tsv                   #kid dad mom (A, B resolvable; C unresolvable)
  annot.tsv                   per-site lookup mock_vep.py turns into a VEP CSQ (frequency +
                              CLIN_SIG + CADD included — there is no external sites VCF)
  mutrate.tsv, constraint.tsv gene tables for Step 6
  config.mock.yaml            config pointing at the above (concrete paths, ephemeral)

All variants are SNVs (no indel left-align ambiguity). REF matches the reference at
each position. Not committed data — generated into the git-ignored work dir.
"""
from __future__ import annotations

import argparse
import json
import os

# chrM is present ON PURPOSE: Step 1 must EXCLUDE it (it is out of scope; see
# 01_make_cohort_sites.sh EXCLUDE_CONTIGS). A mock without chrM cannot prove the filter works.
CONTIGS = {"chr1": 20000, "chr2": 20000, "chrX": 2782200, "chrM": 16569}  # chrX > PAR1 end (2,781,479)
BASES = "ACGT"


def refbase(pos):
    return BASES[pos % 4]


def altbase(pos):
    return BASES[(pos + 1) % 4]


def ad(gt, dp, n_alt=1):
    """Biallelic-style AD for a genotype. n_alt > 1 pads the Number=R array to REF + n ALTs.

    A non-ref/non-ref (1/2) genotype must NOT use this — its AD is the whole point of the test
    and is supplied explicitly via `adov` (see the GENECH2 comp-het case).
    """
    pad = ",0" * (n_alt - 1)
    if gt == "0/0":
        return f"{dp},0{pad}"
    if gt == "1/1":
        return f"0,{dp}{pad}"
    h = dp // 2
    return f"{dp - h},{h}{pad}"


# Sample groupings per VCF file (note fileB order is shuffled and has an extra sib).
FILES = {
    "A": ["CH_A", "FA_A", "MO_A"],
    "B": ["MO_B", "SIB_B", "CH_B", "FA_B"],
    "C": ["CH_C", "FA_C"],
}

# Each variant: file, chrom, pos, gene, csq, impact, cadd, af (None = absent from gnomAD),
# af_pop (which gnomAD CSQ population carries `af`), clnsig, filter, hidenovo, gts.
#
# VEP-only contract: there is no gnomAD/ClinVar/dbNSFP/SpliceAI/LOFTEE file to mock, because the
# pipeline no longer reads one. Frequency and CLIN_SIG are emitted INTO the CSQ by mock_vep.py,
# exactly as a real `vep --af_gnomade --af_gnomadg --check_existing` run would.
#
# af_pop defaults to a grpmax-ELIGIBLE group, so `af` drives annotations.frequency(). Set it to a
# bottlenecked group (ami/asj/fin/mid) to model an allele grpmax deliberately ignores.
V = []


def add(**k):
    k.setdefault("cadd", ""); k.setdefault("af", None); k.setdefault("af_pop", "gnomADe_NFE_AF")
    k.setdefault("spliceai", "")
    # Calibrated missense predictors. Default ABSENT (not 0.0) so most rows exercise the
    # fall-through to CADD/none, and the few that set them exercise the calibrated limbs.
    k.setdefault("revel", ""); k.setdefault("alphamissense", "")
    # gnomAD joint transfer fields. Default ABSENT so most rows exercise the proxy fallback.
    k.setdefault("faf95", ""); k.setdefault("faf95_group", ""); k.setdefault("nhomalt", "")
    k.setdefault("clnsig", ""); k.setdefault("filter", "PASS"); k.setdefault("hidenovo", "")
    V.append(k)


# --- Trio A (autosomal), CH_A female ---------------------------------------
# 1) de novo, HIGH LoF, absent -> expect mode=denovo
add(file="A", chrom="chr1", pos=5000, gene="GENE1", csq="stop_gained", impact="HIGH",
    hidenovo="CH_A",
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/0", 99, 40)})
# 2) homozygous recessive, MODERATE missense, rare -> mode=hom_recessive
# REVEL >= 0.773 (Pejaver moderate) with a LOW cadd: the ladder must report revel, not cadd,
# and must reach V4 — this is the whole point of adding a calibrated predictor.
add(file="A", chrom="chr1", pos=8000, gene="GENE2", csq="missense_variant", impact="MODERATE",
    af=5e-4, cadd="3", revel="0.85", faf95="8e-05", faf95_group="nfe", nhomalt="7",
    gts={"CH_A": ("1/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/1", 99, 40)})
# 3+4) compound het in GENE3 (var3 maternal, var4 paternal) -> mode=compound_het
# AlphaMissense only (REVEL absent): the ladder must fall through to it rather than to CADD.
# ALSO the faf95_zero case: gnomAD HAS this allele (nhomalt present => a record exists) but
# published no faf95, so faf95 is 0 and the variant must survive a gate its inflated
# point-estimate proxy (1e-3) would fail. 96.5% of that class are AC<=2 singletons.
add(file="A", chrom="chr2", pos=5000, gene="GENE3", csq="missense_variant", impact="MODERATE",
    af=1e-3, alphamissense="0.9", nhomalt="0",
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/1", 99, 40)})
add(file="A", chrom="chr2", pos=6000, gene="GENE3", csq="missense_variant", impact="MODERATE",
    af=1e-3,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# 5) common (BA1) -> dropped at Step 3 (never a candidate)
add(file="A", chrom="chr1", pos=12000, gene="GENE4", csq="missense_variant", impact="MODERATE",
    af=0.2,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# 6) low-GQ de-novo-looking (no hiConfDeNovo) -> passes Step 3, FAILS Step 5 QC
add(file="A", chrom="chr1", pos=15000, gene="GENE5", csq="stop_gained", impact="HIGH", gts={"CH_A": ("0/1", 12, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/0", 99, 40)})
# 7) non-PASS -> dropped at Step 1
add(file="A", chrom="chr1", pos=17000, gene="GENE6", csq="stop_gained", impact="HIGH", filter="VQSRTrancheSNP99.00to99.90+", hidenovo="CH_A",
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/0", 99, 40)})
# 8) ClinVar P/LP but LOW impact (synonymous) -> kept at Step 3 via clinvar_plp override
add(file="A", chrom="chr2", pos=8000, gene="GENE7", csq="synonymous_variant", impact="LOW",
    af=2e-4, clnsig="pathogenic",  # > dominant_max: kept via ClinVar but not a dominant call
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- Trio B (CH_B male), family VCF with extra sibling ----------------------
# 9) de novo in GENE1 (recurrent across trios A+B) -> mode=denovo; drives Step-6 burden
add(file="B", chrom="chr1", pos=5100, gene="GENE1", csq="stop_gained", impact="HIGH",
    hidenovo="CH_B",
    gts={"CH_B": ("0/1", 99, 40), "FA_B": ("0/0", 99, 40), "MO_B": ("0/0", 99, 40),
         "SIB_B": ("0/0", 99, 40)})
# 10) X-linked recessive (male hemizygous), carrier mother -> mode=x_linked_recessive
add(file="B", chrom="chrX", pos=2781600, gene="GENEX", csq="missense_variant", impact="MODERATE",
    af=1e-4,
    gts={"CH_B": ("1/1", 99, 40), "FA_B": ("0/0", 99, 40), "MO_B": ("0/1", 99, 40),
         "SIB_B": ("0/0", 99, 40)})
# 11-13) chrX filler (common) so CH_B is inferred MALE (hemizygous alt -> low het ratio)
for i, p in enumerate((2781700, 2781800, 2781900)):
    add(file="B", chrom="chrX", pos=p, gene=f"XFILL{i}", csq="missense_variant", impact="MODERATE",
        af=0.3,
        gts={"CH_B": ("1/1", 99, 40), "FA_B": ("1/1", 99, 40), "MO_B": ("0/1", 99, 40),
             "SIB_B": ("0/1", 99, 40)})

# --- DOMINANT RECURRENCE: same rare functional inherited het in GENED, in BOTH trios
#     (CH_A inherits from dad, CH_B from mom) -> Step 6 nominates GENED (n_dominant=2) ---
add(file="A", chrom="chr2", pos=10000, gene="GENED", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
add(file="B", chrom="chr2", pos=10000, gene="GENED", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_B": ("0/1", 99, 40), "MO_B": ("0/1", 99, 40), "FA_B": ("0/0", 99, 40),
         "SIB_B": ("0/0", 99, 40)})

# --- CONTAMINATION: a common hom-alt site where FA_B carries reference reads at a hom-alt
#     genotype (verifyBamID-style contamination). Common (af=0.3) so Step 3 drops it and it
#     never becomes a candidate — it only exercises Step 0's CHARR gate. Mendelian-consistent
#     (all 1/1) so it adds no MIE. Expect: CH_B contam_flag=1 (dad), kid/mom clean. ---
add(file="B", chrom="chr2", pos=14000, gene="XCONTAM", csq="missense_variant", impact="MODERATE",
    af=0.3,
    gts={"CH_B": ("1/1", 99, 40), "FA_B": ("1/1", 99, 40), "MO_B": ("1/1", 99, 40),
         "SIB_B": ("1/1", 99, 40)},
    adov={"FA_B": "6,34"})   # 6 ref reads at a hom-alt site -> CHARR 0.15 > 0.02 threshold

# --- X-linked recessive with an AFFECTED FATHER (hom-alt hemizygous) + carrier mother: the son
#     must STILL be called (father transmits Y, not X, to a son) -> tests the father-genotype
#     relaxation. Flag 'father_carries_x_allele' expected. ---
add(file="B", chrom="chrX", pos=2782000, gene="GENEXAF", csq="missense_variant", impact="MODERATE",
    af=1e-4,
    gts={"CH_B": ("1/1", 99, 40), "FA_B": ("1/1", 99, 40), "MO_B": ("0/1", 99, 40),
         "SIB_B": ("0/0", 99, 40)})

# --- autosomal hom-recessive with a HOM-ALT parent (consanguinity-like): FA_A hom-alt, MO_A het,
#     CH_A hom-alt -> hom_recessive via the {HET,HOM_ALT} carrier rule (tests carrier_ok HOM_ALT). ---
add(file="A", chrom="chr1", pos=8500, gene="GENE2H", csq="missense_variant", impact="MODERATE",
    af=5e-4,
    gts={"CH_A": ("1/1", 99, 40), "FA_A": ("1/1", 99, 40), "MO_A": ("0/1", 99, 40)})

# --- DISTINCT-variant dominant recurrence across trios A+B in GENEDD (two DIFFERENT rare hets) ->
#     the stronger 'gene signal'; recurrence_kind=distinct_variant (ranks above same-variant GENED). ---
add(file="A", chrom="chr2", pos=11000, gene="GENEDD", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
add(file="B", chrom="chr2", pos=11100, gene="GENEDD", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_B": ("0/1", 99, 40), "MO_B": ("0/1", 99, 40), "FA_B": ("0/0", 99, 40),
         "SIB_B": ("0/0", 99, 40)})

# --- comp-het CIS rejection: two rare functional hets in GENEC, BOTH inherited from mom
#     (cis). compound_het requires TRANS (mat x pat), so these must NOT be paired; each is a
#     dominant (maternal-origin) call instead. Guards the trans-pairing logic. ---
add(file="A", chrom="chr2", pos=16000, gene="GENEC", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/1", 99, 40)})
add(file="A", chrom="chr2", pos=16100, gene="GENEC", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/1", 99, 40)})

# --- comp-het with a HOM-ALT transmitting parent. A 1/1 parent transmits the alt OBLIGATELY, so
#     for a HET child the parent of origin is DETERMINISTIC (the alt came from the 1/1 parent, the
#     ref from the other) — it is NOT the 50/50 "both" case. Collapsing it to "both" barred the
#     pair from trans-pairing, and because both alleles sit in the recessive band (3e-3 > the 1e-4
#     dominant gate) neither leg could fall through to a dominant call either: a phase-CONFIRMED
#     biallelic hit vanished from candidates.calls.tsv under NO mode at all. Expect a TRANS pair.
#     (Highest-yield real instance is chrX, where a diploid caller renders a hemizygous carrier
#     father as 1/1 — the autosomal form is modelled here.) ---
add(file="A", chrom="chr2", pos=12000, gene="GENEHOMALT", csq="missense_variant", impact="MODERATE",
    af=3e-3,   # mother 1/1 -> obligate maternal transmission -> origin must resolve to "mat"
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("1/1", 99, 40)})
add(file="A", chrom="chr2", pos=12500, gene="GENEHOMALT", csq="missense_variant", impact="MODERATE",
    af=3e-3,   # ordinary paternal het -> origin "pat"; pairs in trans with the maternal leg above
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- CADD-only keep: a deep-intronic MODIFIER with no impact-based evidence. CADD is the ONLY
#     functional predictor left, so this is the sole path by which any non-coding variant can
#     survive Step 3. If the CADD branch ever breaks, the screen silently goes coding-only and
#     this is the assertion that notices. Inherited het -> also a dominant call in GENEIN. ---
add(file="A", chrom="chr1", pos=18000, gene="GENEIN", csq="intron_variant", impact="MODIFIER",
    cadd="27.5", af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# ...and its control: same intronic MODIFIER, CADD BELOW the cutoff -> must be dropped.
add(file="A", chrom="chr1", pos=18500, gene="GENEINLO", csq="intron_variant", impact="MODIFIER",
    cadd="3.0", af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- SpliceAI keep: a deep-intronic MODIFIER whose CADD is BELOW the cutoff (so CADD cannot
#     rescue it) but whose SpliceAI delta score is >= spliceai_ds_min. This is the sole path by
#     which a cryptic-splice variant that both VEP's positional terms and CADD miss survives
#     Step 3. If the SpliceAI branch ever breaks, deep-intronic splice signal goes invisible and
#     this assertion notices. Inherited het -> also a dominant call in GENESAI (validates the
#     spliceai_ds column flowing into candidates.calls.tsv). ---
add(file="A", chrom="chr1", pos=18700, gene="GENESAI", csq="intron_variant", impact="MODIFIER",
    cadd="3.0", spliceai="0.55", af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# ...and its control: same, but SpliceAI BELOW the cutoff and CADD low -> must be dropped.
add(file="A", chrom="chr1", pos=18800, gene="GENESAILO", csq="intron_variant", impact="MODIFIER",
    cadd="3.0", spliceai="0.10", af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- FOUNDER-POPULATION allele: frequent ONLY in a bottlenecked group gnomAD's grpmax excludes
#     (mid, ~AN 700). MAX_AF reports 0.002 — 20x over dominant_max — but grpmax-eligible groups
#     report nothing, so annotations.frequency() must return None and the variant must SURVIVE as
#     a dominant candidate. This is the concrete false-negative that using VEP's MAX_AF as the
#     rarity field would cause, and the reason frequency() reads only GRPMAX_POPS. ---
add(file="A", chrom="chr2", pos=17000, gene="GENEFND", csq="missense_variant", impact="MODERATE",
    af=0.002, af_pop="gnomADe_MID_AF",
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- MULTIALLELIC trans compound het: child is 1/2 (one allele from each parent), which is the
#     textbook presentation of a recessive diagnosis. `bcftools norm -m-` splits this into two
#     records and, WITHOUT --keep-sum AD, discards the other ALT's reads from each leg — leaving
#     ref_ad~0, so allele_balance() reads ~1.0, the het band rejects BOTH legs, and the whole
#     compound_het vanishes with no warning and no audit counter. This case is the regression
#     test for that (Step 4's --keep-sum AD): expect a compound_het pair in GENECH2.
#     AD is explicit because it IS the thing under test:
#       CH_A 1/2 -> 0 ref, 19 for G, 20 for T  (legs must read AB 0.487 / 0.513, not 1.0)
#       FA_A 0/1 -> transmits G ;  MO_A 0/2 -> transmits T   => trans, not cis. ---
add(file="A", chrom="chr1", pos=19000, gene="GENECH2", csq="missense_variant", impact="MODERATE",
    alt2="T", af=5e-5,
    gts={"CH_A": ("1/2", 99, 39), "FA_A": ("0/1", 99, 38), "MO_A": ("0/2", 99, 38)},
    adov={"CH_A": "0,19,20", "FA_A": "18,20,0", "MO_A": "17,0,21"})

# --- chrM: MUST be excluded at Step 1 (out of scope; a dedicated mtDNA pipeline owns it).
#     Modelled on the real failure: m.8860A>G is a near-fixed rCRS haplogroup variant, so the
#     WHOLE TRIO is hom-alt. Left un-excluded it fires hom_recessive in every trio, and with no
#     gnomAD mito AF the rarity gate passes unconditionally -> Step 6 floors q -> p ~ 1e-12 ->
#     it lands in the recurrent exome-wide-significant tier above real nuclear candidates. ---
add(file="A", chrom="chrM", pos=8860, gene="MT-ATP6", csq="missense_variant", impact="MODERATE",
    gts={"CH_A": ("1/1", 99, 400), "FA_A": ("1/1", 99, 400), "MO_A": ("1/1", 99, 400)})
add(file="B", chrom="chrM", pos=8860, gene="MT-ATP6", csq="missense_variant", impact="MODERATE",
    gts={"CH_B": ("1/1", 99, 400), "FA_B": ("1/1", 99, 400), "MO_B": ("1/1", 99, 400),
         "SIB_B": ("1/1", 99, 400)})

# --- Trio C: duo only (mom MO_C absent everywhere) -> resolver unresolved ---
add(file="C", chrom="chr1", pos=9000, gene="GENE8", csq="missense_variant", impact="MODERATE",
    af=1e-3, gts={"CH_C": ("0/1", 99, 40), "FA_C": ("0/1", 99, 40)})

# --- STEP 9 (prioritization) fixtures ------------------------------------------------------
# An ARTIFACT LOCUS: several rare functional inherited hets in one gene whose mutational target
# is tiny (see mutational_target.tsv), carrying every corroborating signal. Step 9 must place it
# in a DOWN-WEIGHT tier — and must not remove a single one of its rows (never-drop). Named OR4Q3
# so the olfactory-receptor family regex is exercised on a real pattern, and it is the largest
# real T3 member from the validation cohort (n=229 at 546x its target).
for i, pos in enumerate((18000, 18100, 18200)):
    add(file="A", chrom="chr2", pos=pos, gene="OR4Q3", csq="missense_variant", impact="MODERATE",
        af=5e-5,
        gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
add(file="B", chrom="chr2", pos=18300, gene="OR4Q3", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_B": ("0/1", 99, 40), "MO_B": ("0/1", 99, 40), "FA_B": ("0/0", 99, 40),
         "SIB_B": ("0/0", 99, 40)})

# --- THE POSITIVE-CONTROL GUARD, and the reason it is a fixture rather than a comment.
#     GENE1 is given the SAME tiny mutational target and the SAME four corroborating signals as
#     the artifact locus above, plus the same pileup shape — so on the statistics alone it earns
#     T3_strong_downweight. It is in established_genes.txt, so the auditable control ceiling must
#     cap it at T1_watch and flag it established_gene_high_excess instead. If that ceiling ever
#     regresses, a real predisposition gene's variants get penalised silently; this is the one
#     failure the whole down-weight scheme exists to prevent. (GENE1 already carries the two de
#     novo calls, so these rows push it into the same extreme-excess shape.) ---
for pos in (18400, 18500):
    add(file="A", chrom="chr2", pos=pos, gene="GENE1", csq="missense_variant", impact="MODERATE",
        af=5e-5,
        gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- CDS-FALLBACK offset: GENENOMU has no mu_* values, only a cds_length, so Step 9 must fall
#     back to the CDS-length regression, mark E_source=cds_fallback, and never let the gene reach
#     T3 (a +/-30% offset cannot support that claim). ---
add(file="A", chrom="chr2", pos=18600, gene="GENENOMU", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- A V0 "molecularly benign PREDICTION" in a HIGHLY CONSTRAINED gene (GENE1: pLI 0.98,
#     LOEUF 0.20). Mechanism gating must zero the constraint term and cap the total at 0 — no
#     amount of gene-level enthusiasm may rescue a benign prediction. Kept by the CADD rung so it
#     reaches Step 9 at all (CADD 26 clears the screen), then scored benign by the TIER rule,
#     which reads spliceai_ds < 0.1 AND cadd < 15: the two thresholds are deliberately different,
#     so this row is a V1 discovery rank rather than a V0...
add(file="A", chrom="chr2", pos=18700, gene="GENE1", csq="intron_variant", impact="MODIFIER",
    cadd="26", spliceai="0.01", af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# ...and THIS one is the real V0: a ClinVar P/LP assertion carries it through the screen (the
# clinvar_plp override), so it arrives at Step 9 with a LOW impact, a sub-0.1 SpliceAI score and
# a sub-15 CADD. The ClinVar term (+4) and the constraint of GENE1 both push it up; the V0 cap
# must hold the total at 0 regardless. That is the mechanism-gating assertion in its strongest
# form: even a P/LP assertion in a constrained gene cannot lift a molecularly-benign prediction.
add(file="A", chrom="chr2", pos=18800, gene="GENE1", csq="synonymous_variant", impact="LOW",
    cadd="4", spliceai="0.01", clnsig="pathogenic", af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- MODEL-MISMATCH fixture. A single inherited HET pLoF in GENERC, a gene the mock curates as
#     canonically autosomal RECESSIVE (gene_moi.tsv) while the weighted overlay carries a
#     het-carrier hypothesis for it (the FA/HR shape: PMID 40906985 evidence is about heterozygous
#     carriers, whereas the canonical model is biallelic Fanconi anemia). Step 9 must (a) still
#     apply the gene prior — the overlay join is by SYMBOL, never routed by MOI — and (b) report
#     the mismatch as an informational flag with ZERO penalty. ---
add(file="A", chrom="chr2", pos=18900, gene="GENERC", csq="stop_gained", impact="HIGH",
    cadd="38", af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# ...and a somatic-driver row's gene, so the zero-prior path is exercised on a real call.
add(file="A", chrom="chr2", pos=19100, gene="GENESOM", csq="stop_gained", impact="HIGH",
    cadd="38", af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})


def write_reference(path):
    with open(path, "w") as fh:
        for c, n in CONTIGS.items():
            seq = bytearray(b"A" * n)
            for v in V:
                if v["chrom"] == c:
                    seq[v["pos"] - 1] = ord(refbase(v["pos"]))
            fh.write(f">{c}\n")
            s = seq.decode()
            for i in range(0, n, 60):
                fh.write(s[i:i + 60] + "\n")


VCF_HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=lowGQ,Description="GQ < 20.0">
##FILTER=<ID=VQSRTrancheSNP99.00to99.90+,Description="VQSR tranche">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=hiConfDeNovo,Number=1,Type=String,Description="High-confidence de novo child list">
##INFO=<ID=loConfDeNovo,Number=1,Type=String,Description="Low-confidence de novo child list">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
"""


def write_vcf(path, samples, variants):
    with open(path, "w") as fh:
        fh.write(VCF_HEADER)
        for c, n in CONTIGS.items():
            fh.write(f"##contig=<ID={c},length={n}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for v in sorted(variants, key=lambda x: (list(CONTIGS).index(x["chrom"]), x["pos"])):
            ref = refbase(v["pos"])
            alt = altbase(v["pos"])
            n_alt = 1
            if v.get("alt2"):            # multiallelic: ALT becomes "G,T"
                alt = f"{alt},{v['alt2']}"
                n_alt = 2
            # AC/AF are Number=A: ONE VALUE PER ALT, counted per allele index. A single summed
            # count is not just imprecise, it is malformed for a multiallelic — bcftools rejects
            # it outright ("wrong number of fields in INFO/AF ... expected 2, found 1").
            an = 2 * len(samples)
            acs = [sum(g[0].replace("|", "/").split("/").count(str(i))
                       for g in v["gts"].values())
                   for i in range(1, n_alt + 1)]
            info = ("AC=" + ",".join(str(a) for a in acs) + f";AN={an};AF="
                    + ",".join(f"{a / an:.4g}" for a in acs))
            if v["hidenovo"]:
                info += f";hiConfDeNovo={v['hidenovo']}"
            cells = []
            adov = v.get("adov", {})     # per-sample AD override (contamination; multiallelic AD)
            for s in samples:
                gt, gq, dp = v["gts"][s]
                a = adov.get(s, ad(gt, dp, n_alt))
                cells.append(f"{gt}:{a}:{dp}:{gq}")
            fh.write(f"{v['chrom']}\t{v['pos']}\t.\t{ref}\t{alt}\t100\t{v['filter']}\t"
                     f"{info}\tGT:AD:DP:GQ\t" + "\t".join(cells) + "\n")


def write_sites(path, header_info, want):
    """Write a sites-only VCF with the given INFO records (want: list of (v, infostr))."""
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for line in header_info:
            fh.write(line + "\n")
        for c, n in CONTIGS.items():
            fh.write(f"##contig=<ID={c},length={n}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        seen = set()
        for v, info in sorted(want, key=lambda x: (list(CONTIGS).index(x[0]["chrom"]), x[0]["pos"])):
            key = (v["chrom"], v["pos"])
            if key in seen:
                continue
            seen.add(key)
            ref, alt = refbase(v["pos"]), altbase(v["pos"])
            fh.write(f"{v['chrom']}\t{v['pos']}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")


def write_source_sams(W):
    """Write a tiny per-sample SAM (a few reads over each of the sample's variants,
    coordinate-sorted) and return {sample: cram_path}. run_integration.sh converts each
    SAM to a sorted+indexed CRAM so Step 8 can slice mini-CRAMs from it."""
    sam_dir = os.path.join(W, "crams_src")
    os.makedirs(sam_dir, exist_ok=True)
    L = 60
    per_sample = {}  # sample -> list of (contig_idx, pos, samline)
    counter = 0
    for fk, samples in FILES.items():
        for v in [x for x in V if x["file"] == fk]:
            cidx = list(CONTIGS).index(v["chrom"])
            start = max(1, v["pos"] - 30)
            for s in samples:
                for _ in range(4):
                    counter += 1
                    line = (f"r{counter}\t0\t{v['chrom']}\t{start}\t60\t{L}M\t*\t0\t0\t"
                            f"{'A' * L}\t{'I' * L}")
                    per_sample.setdefault(s, []).append((cidx, start, line))
    cram_map = {}
    for s, reads in per_sample.items():
        reads.sort(key=lambda x: (x[0], x[1]))
        path = os.path.join(sam_dir, f"{s}.sam")
        with open(path, "w") as fh:
            fh.write("@HD\tVN:1.6\tSO:coordinate\n")
            for c, n in CONTIGS.items():
                fh.write(f"@SQ\tSN:{c}\tLN:{n}\n")
            for _, _, line in reads:
                fh.write(line + "\n")
        cram_map[s] = os.path.join(sam_dir, f"{s}.cram")  # produced from the SAM by the runner
    with open(os.path.join(W, "cram_map.tsv"), "w") as fh:
        for s, p in sorted(cram_map.items()):
            fh.write(f"{s}\t{p}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    W = os.path.abspath(args.out)
    os.makedirs(os.path.join(W, "vcfs"), exist_ok=True)

    write_reference(os.path.join(W, "reference.fa"))
    write_source_sams(W)
    for fk, samples in FILES.items():
        write_vcf(os.path.join(W, "vcfs", f"file{fk}.vcf"), samples,
                  [v for v in V if v["file"] == fk])

    # A mock gnomAD joint slim IS written below (it exercises Step 2's real transfer). No
    # clinvar.vcf: that transfer degrades with a warning, which the run asserts. Otherwise
    # nothing is transferred from
    # an external sites VCF, so there is nothing to mock. Frequency + CLIN_SIG go into the CSQ
    # (see mock_vep.py), which is where a real `vep --af_gnomade --check_existing` puts them.

    # annotation lookup for mock_vep.py -> becomes the CSQ.
    # Keyed per (chrom, pos, ALT), not per site: Step 1's `norm -m-` splits a multiallelic into
    # one record per ALT, and mock_vep.py matches on the exact allele — so a site with alt2 needs
    # BOTH alleles here or the second leg reaches Step 3 with no CSQ and is dropped as
    # not_functional, quietly destroying the very comp-het the multiallelic case exists to test.
    with open(os.path.join(W, "annot.tsv"), "w") as fh:
        fh.write("chrom\tpos\tref\talt\tgene\tcsq\timpact\tcadd\taf\taf_pop\tclnsig\t"
                 "spliceai\trevel\talphamissense\tfaf95\tfaf95_group\tnhomalt\n")
        seen = set()
        for v in V:
            alts = [altbase(v["pos"])]
            if v.get("alt2"):
                alts.append(v["alt2"])
            for a in alts:
                key = (v["chrom"], v["pos"], a)
                if key in seen:
                    continue
                seen.add(key)
                af = "" if v["af"] is None else f"{v['af']:.6g}"
                fh.write(f"{v['chrom']}\t{v['pos']}\t{refbase(v['pos'])}\t{a}\t"
                         f"{v['gene']}\t{v['csq']}\t{v['impact']}\t{v['cadd']}\t{af}\t"
                         f"{v['af_pop']}\t{v['clnsig']}\t{v['spliceai']}\t"
                         f"{v['revel']}\t{v['alphamissense']}\t"
                         f"{v['faf95']}\t{v['faf95_group']}\t{v['nhomalt']}\n")

    # Step-6 tables
    with open(os.path.join(W, "mutrate.tsv"), "w") as fh:
        fh.write("gene\tmut_lof\tmut_mis\n")
        for g in ("GENE1", "GENE2", "GENE3", "GENE5", "GENEX"):
            fh.write(f"{g}\t1e-6\t1e-5\n")
    with open(os.path.join(W, "constraint.tsv"), "w") as fh:
        fh.write("gene\toe_lof_upper\tpli\ts_het\n")
        fh.write("GENE1\t0.2\t0.98\t0.15\n")
        fh.write("GENED\t0.25\t0.95\t0.12\n")

    # --- Step-9 tables -----------------------------------------------------------------
    # The MUTATIONAL-TARGET table: the offset (mu_mis/mu_syn/mu_lof) for the excess statistic
    # plus three of the six artifact signals (oe_syn, classic_caf, constraint_flag) and
    # cds_length for the fallback offset. Engineered so Step 9's every branch fires on a mock
    # this small:
    #   * GENEART   — an ARTIFACT locus: a tiny mutational target against several candidate rows,
    #                 with every corroborating signal set (oe_syn 1.7, caf 0, mis_too_many, a
    #                 segdup fraction, and a symbol matching the OR* family regex). Must reach a
    #                 down-weight tier. It is deliberately named OR4Q3 so the FAMILY regex is
    #                 exercised on a real pattern rather than a synthetic one.
    #   * GENE1     — the POSITIVE CONTROL: constrained, in the established-gene union, and given
    #                 the SAME extreme excess shape as the artifact locus. It must NEVER reach
    #                 T2/T3 (the control ceiling), and must carry review_flag =
    #                 established_gene_high_excess. This is the sensitivity guard that matters.
    #   * GENENOMU  — no mu_* values at all, only a cds_length => the CDS-fallback offset AND its
    #                 tier ceiling.
    #   * everything else is background at a comfortable target so the null has a bulk to fit.
    # NB the mu values are scaled for a 2-trio mock, not a real cohort: C is re-fit per run, so
    # what matters is the RATIO between GENEART/GENE1's targets and the background's.
    mut_genes = sorted({v["gene"] for v in V} | {f"BG{i:03d}" for i in range(120)})
    with open(os.path.join(W, "mutational_target.tsv"), "w") as fh:
        fh.write("gene\tgene_id\tmu_mis\tmu_syn\tmu_lof\toe_syn\tclassic_caf\t"
                 "constraint_flag\tcds_length\tpli\toe_lof_upper\tsegdup98_frac\n")
        for i, g in enumerate(mut_genes):
            if g == "OR4Q3":
                fh.write(f"{g}\tENSG{i:011d}\t2e-7\t8e-8\t1e-8\t1.70\t0\tmis_too_many"
                         f"\t900\t0.01\t1.90\t0.55\n")
            elif g == "GENE1":
                # same tiny target as the artifact locus: the ONLY thing protecting it is the
                # established-gene ceiling, which is exactly what the assertion checks
                fh.write(f"{g}\tENSG{i:011d}\t2e-7\t8e-8\t1e-8\t1.70\t0\tmis_too_many"
                         f"\t1182\t0.98\t0.20\t0.55\n")
            elif g == "GENENOMU":
                fh.write(f"{g}\tENSG{i:011d}\t\t\t\t1.00\t2e-4\t\t1500\t0.10\t1.10\t0.00\n")
            else:
                fh.write(f"{g}\tENSG{i:011d}\t9e-6\t4e-6\t\t1.00\t2e-4\t\t1500\t0.10"
                         f"\t1.10\t0.00\n")
    # The phenotype-AGNOSTIC established-gene-validity union. Small on purpose (the mock lowers
    # prioritization.gene_downweight.min_control_genes to match), and it must contain GENE1.
    with open(os.path.join(W, "established_genes.txt"), "w") as fh:
        fh.write("# phenotype-agnostic gene-validity union (mock)\nGENE1\nGENE2\nGENEX\n")
    # Curated MOI, so the moi_coherence term is exercised in all three states: GENED is curated
    # AD (coherent with its dominant calls), GENE2 is curated AR, and every other gene is
    # UNCURATED => moi_unknown, which must be EXACTLY neutral (the novel-gene case).
    with open(os.path.join(W, "gene_moi.tsv"), "w") as fh:
        # GENERC is curated as a BARE canonical AR — deliberately WITHOUT the overlay's
        # carrier-risk annotation — so the het observation there is a genuine
        # canonical-MOI mismatch and the suppression path is what gets tested.
        fh.write("gene\tmoi\nGENED\tAD\nGENE2\tAR\nGENE3\tAR\nGENERC\tAR\n")
    # The Class-B overlay. Present as a FILE but left DISABLED in the mock config, so the
    # integration run asserts the default contract: rank_prior == rank_agnostic.
    with open(os.path.join(W, "phenotype_overlay.txt"), "w") as fh:
        fh.write("# Class-B overlay (mock) — NOT enabled in config.mock.yaml\nGENED\n")
    # A WEIGHTED overlay in the real schema, for the `--gene-prior` smoke run in
    # run_integration.sh. Every hazard is represented, with a header row as a decoy:
    #   GENE2   T1, weight 1.0 — the full prior
    #   GENE3   T3, weight 0.35 — a GWAS-locus row, so a fraction of the prior
    #   GENERC  T2, weight 0.6, canonically RECESSIVE with an explicit het-carrier hypothesis —
    #           the FA/HR shape. Its prior must apply to a HET observation and the MOI mismatch
    #           must be a flag, not a penalty.
    #   GENESOM T4, weight 0.15, somatic_driver_not_germline — must contribute ZERO despite the
    #           weight, and must still be REPORTED.
    #   GENEDD  no gene row at all; admitted only via the gene set below at the SET weight.
    with open(os.path.join(W, "phenotype_overlay_weighted.tsv"), "w") as fh:
        fh.write("gene\ttier\tprior_weight\tevidence_class\tmoi\treplication\tgene_sets\tpmids\n")
        fh.write("GENE2\tT1\t1.0\trare_variant_syndromic\tAD\treplicated\t\t12345678\n")
        fh.write("GENE3\tT3\t0.35\tgwas_common_variant_locus\tcomplex_common_variant\t"
                 "unreplicated\t\t12345679\n")
        fh.write("GENERC\tT2\t0.6\trare_variant_association_single_study_significant\t"
                 "AR_biallelic;heterozygous_carrier_risk_proposed\tunreplicated\tMOCK_PATHWAY\t"
                 "40906985\n")
        fh.write("GENESOM\tT4\t0.15\tsomatic_driver_not_germline\tNA_somatic\tunreplicated\t\t"
                 "12345680\n")
    # The JSON sidecar carrying the gene set. GENERC has BOTH a gene row (0.6) and set membership
    # (0.6) from the same study — they must combine by MAX, never sum. GENEDD has only the set.
    with open(os.path.join(W, "phenotype_overlay_weighted.json"), "w") as fh:
        json.dump({"resource": "mock", "gene_sets": {"MOCK_PATHWAY": {
            "label": "mock pathway-collapsed prior (same study as the per-gene rows)",
            "members": ["GENERC", "GENEDD"], "n_members": 2,
            "prior_weight": 0.6, "source_pmid": "40906985"}}}, fh, indent=1)

    # trios file (#kid dad mom); C is unresolvable (MO_C absent everywhere)
    with open(os.path.join(W, "trios.tsv"), "w") as fh:
        fh.write("#kid\tdad\tmom\nCH_A\tFA_A\tMO_A\nCH_B\tFA_B\tMO_B\nCH_C\tFA_C\tMO_C\n")

    # minimal config (defaults fill in thresholds); concrete ephemeral paths
    # A REAL (tiny) gnomAD joint slim, so Step 2 exercises the actual `bcftools annotate`
    # transfer and its 0-match guard rather than a faked INFO field. Field names are v4.1's
    # exactly (fafmax_faf95_max_joint, ...) — the rename to gnomad_* happens in Step 2, so a
    # typo there is caught here.
    gn = os.path.join(W, "gnomad.slim.vcf")
    with open(gn, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c in sorted({v["chrom"] for v in V}):
            fh.write(f"##contig=<ID={c}>\n")
        fh.write('##INFO=<ID=AF_joint,Number=A,Type=Float,Description="x">\n'
                 '##INFO=<ID=AF_grpmax_joint,Number=A,Type=Float,Description="x">\n'
                 '##INFO=<ID=fafmax_faf95_max_joint,Number=A,Type=Float,Description="x">\n'
                 '##INFO=<ID=fafmax_faf95_max_gen_anc_joint,Number=A,Type=String,Description="x">\n'
                 '##INFO=<ID=nhomalt_joint,Number=A,Type=Integer,Description="x">\n'
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        rows = []
        for v in V:
            if not (v["faf95"] or v["nhomalt"]):
                continue
            info = [f"AF_joint={v['af'] or 1e-4:.6g}",
                    f"AF_grpmax_joint={v['af'] or 1e-4:.6g}"]
            if v["faf95"]:
                info.append(f"fafmax_faf95_max_joint={v['faf95']}")
                info.append(f"fafmax_faf95_max_gen_anc_joint={v['faf95_group'] or 'nfe'}")
            if v["nhomalt"] != "":
                info.append(f"nhomalt_joint={v['nhomalt']}")
            rows.append((v["chrom"], v["pos"], refbase(v["pos"]), altbase(v["pos"]),
                         ";".join(info)))
        for c, pos, r, a, info in sorted(rows, key=lambda x: (x[0], x[1])):
            fh.write(f"{c}\t{pos}\t.\t{r}\t{a}\t.\t.\t{info}\n")
    # bgzip/tabix are run by run_integration.sh, which owns every tool invocation here.

    with open(os.path.join(W, "config.mock.yaml"), "w") as fh:
        fh.write(f"""project: {{name: mock, genome_build: GRCh38, output_dir: {W}/work}}
runtime: {{image: none, engine: native, tmpdir: {W}/work/tmp, threads: 1}}
reference: {{fasta: {W}/reference.fa}}
resources:
  # Step 2 ingests this instead of invoking `vep` (mock_vep.py writes it). Everything else in
  # Step 2 — build checks, split-vep, selector, frequency guard — runs for real against it.
  # spliceai_backfill is off by default; pinned explicitly here so a future default flip cannot
  # silently start invoking TensorFlow in CI (host runs have no such env — it ships only in the
  # image). The precomputed SpliceAI keep-path is still exercised: mock_vep.py writes
  # vep_SpliceAI_pred_DS_* directly.
  vep: {{annotated_vcf: {W}/cohort.sites.vep.vcf.gz, version: 115,
         spliceai_backfill: {{enabled: false}}}}
  # The gnomAD joint slim -> Step 2's SECOND bcftools transfer -> real faf95 + nhomalt. Wired
  # here so the integration exercises the transfer, the 0-match guard, and the faf95-before-proxy
  # precedence in annotations.frequency() — not just the proxy fallback.
  gnomad: {{sites_slim: {W}/gnomad.slim.vcf.gz, oracle: faf95}}
  mutation_rate_table: {W}/mutrate.tsv
  constraint: {{gnomad_v2_constraint: {W}/constraint.tsv}}
  cram_map: {W}/cram_map.tsv
prioritization:
  # Step 9. Two mock-scale deviations from the shipped defaults, both deliberate and both about
  # the SIZE of the mock rather than the method:
  #  * min_control_genes 1000 -> 3. The production guard HALTS on a control union that small,
  #    because a truncated union makes the established-gene exemption silently empty. The mock's
  #    union has 3 genes by construction, so the guard has to be lowered for the ceiling to be
  #    exercised at all. Nothing else about the ceiling changes.
  #  * max_downweight_fraction 0.20 -> 1.0. The mock is engineered so a large FRACTION of its
  #    handful of candidates sit in the artifact locus; on a real cohort that figure was 10.4%
  #    and a fifth of the exome would indeed mean a mis-specified null.
  gene_downweight:
    min_control_genes: 3
    max_downweight_fraction: 1.0
    established_genes: {W}/established_genes.txt
  # The Class-B overlay FILE exists in the mock but stays DISABLED, so the integration run
  # asserts the default contract: rank_prior is identical to rank_agnostic.
  composite:
    gene_list_prior: {{enabled: false, path: {W}/phenotype_overlay.txt}}
  resources:
    mutational_target: {W}/mutational_target.tsv
    gene_moi: {W}/gene_moi.tsv
inputs:
  trios_file: {W}/trios.tsv
  vcf_dir: {W}/vcfs
  vcf_list: ""
qc:
  sex_min_sites: 2            # tiny mock has only a handful of chrX sites (prod default is 20)
outputs:
  xlsx: true
  igv: {{enabled: true, padding: 200, genome: hg38}}
""")
    print(f"mock data written to {W}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
