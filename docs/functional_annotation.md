# Functional Consequence & In-Silico Impact Prediction

Assigns molecular consequence and calibrated deleteriousness evidence to rare SNV/indels so downstream tiering can weigh loss-of-function, missense, and splice impact defensibly.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- Annotate with **VEP** (release-matched cache, MANE Select prioritized) + **LOFTEE (GRCh38 fork)** + **dbNSFP** + **SpliceAI (masked)** + **CADD v1.7 PHRED**; version-pin every component and the exact dbNSFP build.
- **VEP IMPACT** (HIGH/MODERATE/LOW/MODIFIER) is a coarse convenience tier only — never rely on it alone for splice-adjacent variants; always layer SpliceAI.
- **pLoF confidence = LOFTEE HC with no flags**; grade PVS1 via the Abou Tayoun 2018 decision tree (NMD-escape / last-exon / 3′-terminal-50 bp / single-exon downgrade), gated on ClinGen gene–disease validity ≥ Moderate and a known LoF mechanism.
- **Missense primary = REVEL** (Pejaver-2022 calibrated): PP3 supporting **≥ 0.644**, moderate **≥ 0.773**, strong **≥ 0.932**; BP4 supporting **≤ 0.290**, moderate **≤ 0.183**. **AlphaMissense** (likely_pathogenic **≥ 0.564**) as orthogonal support only.
- **Splicing = Walker-2023 calibrated masked SpliceAI**: PP3 Δ **≥ 0.2**, BP4 Δ **≤ 0.1**, 0.1–0.2 uninformative; canonical ±1,2 with Δ **≥ 0.5** = high tier. If PVS1(splice) applies, do not also apply PP3.
- **Report ONE predictor per variant — never stack correlated missense tools** as independent evidence; cap each tool at its Pejaver-calibrated strength.
- A gene-specific **ClinGen VCEP** score cutoff always overrides these generic defaults.
- **Scope limits:** SNV/indel only (no CNV/SV); pseudogene / segmental-duplication regions (e.g. PMS2/PMS2CL) are low-confidence from short reads and are flagged, not trusted.

## VEP consequence ontology → IMPACT tiers

Ensembl VEP assigns Sequence Ontology consequence terms and bins them into four IMPACT tiers. The tier is a convenience filter, not a pathogenicity call.

| IMPACT | Representative SO consequences | Interpretation |
|---|---|---|
| **HIGH** | `transcript_ablation`, `splice_acceptor_variant`, `splice_donor_variant`, `stop_gained`, `frameshift_variant`, `stop_lost`, `start_lost`, `transcript_amplification` | Assumed protein-truncating / LoF / NMD-triggering; PVS1-eligible pending LOFTEE + decision tree |
| **MODERATE** | `missense_variant`, `inframe_insertion`, `inframe_deletion`, `protein_altering_variant` | Needs a calibrated missense/impact predictor |
| **LOW** | `splice_region_variant`, `synonymous_variant`, `stop_retained_variant`, `start_retained_variant`, `incomplete_terminal_codon_variant` | Usually low priority — but `splice_region_variant` can be highly deleterious |
| **MODIFIER** | non-coding / intronic / UTR / intergenic / regulatory | Predictions difficult; never sole evidence |

Caveat: `splice_region_variant` is binned LOW yet can abolish splicing. Never use IMPACT alone for any splice-adjacent variant — layer SpliceAI on top (see [Splicing predictors](#splicing-predictors)). Prioritize MANE Select / canonical transcripts, but retain per-transcript annotations so a consequence on a clinically relevant non-canonical transcript is not lost.

## Predicted loss-of-function: LOFTEE and NMD-escape

LOFTEE is a VEP plugin that flags predicted LoF (stop-gained, frameshift, essential splice-site) as **HC** (high-confidence, passes all filters) or **LC** (low-confidence, fails ≥ 1 filter).

- **LC filters** include the **50-bp rule** (variant in the last exon or within 50 bp of the last exon–exon junction → likely NMD-escape), non-canonical splice sites of the affected intron, in-frame rescue by nearby splice sites, and small-intron GT/AG issues.
- **HC flags** (single-exon gene, weak PhyloCSF conservation, NAGNAG acceptor, non-canonical splice) warrant caution but are overridable by gene knowledge — roughly 14% of gnomAD HC pLoF carry a flag.
- Use the maintained **GRCh38-compatible fork** (Ensembl/gnomAD lineage); the original `konradjk/loftee` supports GRCh37 only.

For rare-disease and germline pediatric-cancer screening, treat **HC, no-flag pLoF** in a haploinsufficient / known disease gene as PVS1-eligible. Grade the PVS1 strength with the **Abou Tayoun 2018 / ClinGen SVI decision tree**: NMD-escape (last exon or 3′-terminal 50 bp) and single-exon context downgrade PVS1 from Very Strong. PVS1 application is gated on **ClinGen gene–disease validity ≥ Moderate** and a **known LoF disease mechanism**. See [gene_constraint.md](gene_constraint.md) for the haploinsufficiency weighting (pHaplo, ClinGen HI) that informs "is LoF the mechanism here."

## Missense / pathogenicity predictors — calibrated cutoffs

The reference standard is **Pejaver 2022 (ClinGen SVI)**, which replaced the historic "count concordant tools" approach with single-tool, strength-tiered thresholds derived from local positive predictive value. **Report one predictor per variant** — correlated predictors do not constitute independent lines of evidence, and none may exceed its Pejaver-calibrated strength.

| Predictor | PP3 strong | PP3 moderate | PP3 supporting | BP4 supporting | BP4 moderate | BP4 strong |
|---|---|---|---|---|---|---|
| **REVEL** (primary) | ≥ 0.932 | ≥ 0.773 | ≥ 0.644 | ≤ 0.290 | ≤ 0.183 | ≤ 0.016 |
| **BayesDel** (noAF) | ≥ 0.50 | ≥ 0.27 | ≥ 0.13 | ≤ −0.18 | ≤ −0.36 | — |
| **VEST4** | ≥ 0.965 | ≥ 0.861 | ≥ 0.764 | < 0.764 | — | — |
| **CADD** (PHRED) | — | ≥ 28.1 | ≥ 25.3 | ≤ 22.7 | ≤ 17.3 | — |
| **MutPred2** | ≥ 0.932 | ≥ 0.829 | ≥ 0.737 | ≤ 0.391 | ≤ 0.197 | ≤ 0.010 |

REVEL is the only tool in this set reaching both PP3-strong and BP4-very-strong, and is this pipeline's **primary** missense predictor.

**AlphaMissense** (Cheng 2023) is used as **orthogonal** support, not stacked with REVEL. Its three bins are calibrated to ~90% ClinVar precision:

| AlphaMissense bin | Score |
|---|---|
| likely_benign | ≤ 0.34 |
| ambiguous | 0.34 – 0.564 |
| likely_pathogenic | ≥ 0.564 |

Independent ClinGen-style calibration supports AlphaMissense up to roughly PP3-moderate / BP4-strong.

**CADD v1.7** (Jan 2024; adds protein language models + regulatory CNNs; PHRED-scaled; GRCh38 chr1–22,X,Y; splice-aware since CADD-Splice / v1.6) is a genome-wide fallback that covers non-missense sites, but is a weaker missense discriminator than the meta-predictors above — use it to fill gaps, not to override REVEL/AlphaMissense.

Other predictors present in dbNSFP — **PrimateAI-3D, EVE, MPC, MetaRNN, ClinPred** — are available as orthogonal support. **MPC ≥ 2** strongly up-weights a missense variant in a constrained region (missense Z > 3.09 gives gene-level support; see [gene_constraint.md](gene_constraint.md)). Do not sum multiple correlated missense predictors as independent evidence.

## Splicing predictors

- **SpliceAI** (delta score 0–1). The original heuristics were 0.2 (high recall), 0.5 (recommended), 0.8 (high precision). The **ClinGen SVI Splicing Subgroup (Walker 2023)** calibrated cutoffs for variants outside the essential ±1,2 dinucleotides:

  | SpliceAI Δ (masked) | Evidence | Note |
  |---|---|---|
  | ≥ 0.2 | **PP3** (supporting) | ~78% sensitivity |
  | 0.1 – 0.2 | uninformative | no evidence assigned |
  | ≤ 0.1 | **BP4** | ~87% specificity |

  Canonical ±1,2 splice variants with Δ ≥ 0.5 are treated as a high tier. Use **masked** scores for interpretation and raw scores for discovery. **If PVS1(splice) already applies at any strength, do NOT also apply PP3** — that would double-count the splice evidence.
- **Pangolin** (Genome Biol 2023) shows comparable-or-superior sensitivity to SpliceAI in massively-parallel-assay benchmarks; good orthogonal confirmation at an analogous ~0.2 cutoff.
- **MaxEntScan** — classic motif model, useful for quantifying donor/acceptor strength change at canonical sites; supplementary to SpliceAI.

## dbNSFP as the aggregation resource

**dbNSFP** precomputes transcript-specific scores for nonsynonymous and splice-site SNVs, so a single annotation join delivers all missense predictors without per-tool installs: REVEL, AlphaMissense, CADD, BayesDel, VEST4, MetaRNN, ClinPred, MPC, PrimateAI, EVE, ESM1b, MutPred2, plus conservation (GERP, phyloP).

- Current is **v5.3** (Oct 2025); **v4.9a** (Aug 2024) is the stable v4 academic branch.
- Use the **`a` (academic) branch** — it retains REVEL/CADD/VEST4/PolyPhen2/ClinPred, which the `c` (commercial) branch strips. AlphaMissense (CC-BY) was added at v4.7 and is present in both branches thereafter.
- **Pin the exact dbNSFP build** in the container manifest; a silent build change shifts every downstream score.

## Example annotation commands

Generic, parameterized invocations — substitute release, cache, and plugin data paths via config; never hard-code real paths.

```bash
# VEP: release-matched cache, MANE/canonical prioritized, per-transcript retained,
# LOFTEE (GRCh38 fork), dbNSFP, SpliceAI (masked), CADD.
vep \
  --offline --cache --dir_cache "${VEP_CACHE_DIR}" --cache_version "${VEP_RELEASE}" \
  --assembly GRCh38 --mane --canonical --hgvs --symbol \
  --plugin LoF,loftee_path:"${LOFTEE_DIR}",human_ancestor_fa:"${ANCESTOR_FA}" \
  --plugin dbNSFP,"${DBNSFP_BUNDLE}",REVEL_score,AlphaMissense_score,BayesDel_noAF_score,VEST4_score,MPC_score,CADD_phred \
  --plugin SpliceAI,snv="${SPLICEAI_SNV}",indel="${SPLICEAI_INDEL}" \
  --input_file "${IN_VCF}" --vcf --output_file "${OUT_VCF}" --stats_text
```

```bash
# Extract calibrated fields for tiering (illustrative; adapt CSQ index to your header).
bcftools +split-vep "${OUT_VCF}" \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%IMPACT\t%LoF\t%REVEL_score\t%AlphaMissense_score\t%SpliceAI_pred_DS_AG\n' \
  -d -A tab
```

## Recommended defaults (this pipeline)

For per-trio GATK genotype-refinement VCFs (GRCh38, not cohort-joint), rare-disease + germline pediatric cancer. All values are configurable defaults in `config/config.example.yaml`; a gene-specific ClinGen VCEP cutoff overrides any of them.

| Layer | Default | Notes |
|---|---|---|
| Annotation stack | VEP (release-matched cache, MANE prioritized) + LOFTEE GRCh38 fork + dbNSFP (v4.9a / v5.x) + SpliceAI masked + CADD v1.7 PHRED | Version-pin all; pin the exact dbNSFP build |
| pLoF confidence | LOFTEE **HC, no flags** | NMD-escape / last-exon / 3′-50 bp / single-exon → downgrade PVS1 (Abou Tayoun tree) |
| PVS1 gating | ClinGen gene–disease validity ≥ Moderate + known LoF mechanism | See [gene_constraint.md](gene_constraint.md) |
| Missense (primary) | **REVEL** PP3 supporting ≥ 0.644 / moderate ≥ 0.773 / strong ≥ 0.932; BP4 ≤ 0.290 / ≤ 0.183 | Literature/original REVEL heuristics differ; use these calibrated tiers |
| Missense (orthogonal) | **AlphaMissense** likely_pathogenic ≥ 0.564; ambiguous 0.34–0.564; likely_benign ≤ 0.34 | Do not stack with REVEL as independent evidence |
| Regional missense | **MPC ≥ 2** up-weights; missense Z > 3.09 gene support | |
| Splicing | **SpliceAI (masked)** PP3 Δ ≥ 0.2; BP4 Δ ≤ 0.1; 0.1–0.2 uninformative; canonical ±1,2 Δ ≥ 0.5 high tier | No PP3 if PVS1(splice) applies |
| Benign deprioritize (BP4) | REVEL ≤ 0.183 **AND** AlphaMissense ≤ 0.34 **AND** SpliceAI Δ ≤ 0.1 | Deprioritize, not discard, if in a critical gene |
| Evidence hygiene | One calibrated predictor per variant; cap at Pejaver strength | Never sum correlated predictors |

**Ordering:** apply the rarity gate first (grpmax `faf95`; see [allele_frequency.md](allele_frequency.md)) and genotype/QC gates (see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)) before functional tiering, so functional scoring runs only on variants that survive frequency and quality screening. Calibrated functional evidence then feeds ACMG/AMP classification (see [clinical_classification.md](clinical_classification.md)).

## Scope limitations (state honestly)

- **SNV/indel only.** This layer does not detect CNV/SV, which account for ~10–15% of pediatric-cancer and rare-disease diagnoses (single-exon RB1/SMARCB1/DICER1/NF1 deletions, PMS2 rearrangements). LOFTEE and the missense predictors cannot see these; a future GATK-gCNV / Manta / ExomeDepth module is required.
- **Pseudogene / segmental-duplication regions** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) yield low-confidence short-read calls; functional annotation on paralog-mapping variants is unreliable and those regions are flagged, not trusted.
- **Non-coding regulatory variants** fall in VEP MODIFIER, where prediction is weak; they are not scored to reportable tiers by this layer alone.
- **Predictor calibration** is anchored to ClinVar-derived truth sets (Pejaver, Walker, AlphaMissense); performance on genes/regions under-represented in those sets is not guaranteed. Pipeline-wide sensitivity/precision should be measured against GIAB/CMRG truth sets and a positive-control variant panel.

## Sources

- Pejaver et al. 2022, ClinGen SVI PP3/BP4 calibration, *Am J Hum Genet* — https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/ (DOI 10.1016/j.ajhg.2022.10.013); ClinGen summary — https://clinicalgenome.org/docs/calibration-of-computational-tools-for-missense-variant-pathogenicity-classification-and-clingen-recommendations-for-pp3-bp4-cri/
- Walker et al. 2023, ClinGen SVI Splicing Subgroup, *Am J Hum Genet* — https://pmc.ncbi.nlm.nih.gov/articles/PMC10357475/
- Cheng et al. 2023, AlphaMissense, *Science* — https://www.science.org/doi/10.1126/science.adg7492; thresholds — https://www.ebi.ac.uk/training/online/courses/alphafold/classifying-the-effects-of-missense-variants-using-alphamissense/understanding-pathogenicity-scores-from-alphamissense/
- Abou Tayoun et al. 2018, ClinGen SVI PVS1 decision tree, *Hum Mutat* — https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/
- Ensembl VEP calculated consequences & IMPACT — https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
- LOFTEE README — https://github.com/konradjk/loftee/blob/master/README.md; gnomAD flags context — https://gnomad.broadinstitute.org
- SpliceAI (Illumina) — https://github.com/Illumina/SpliceAI; lookup — https://spliceailookup.broadinstitute.org/
- Pangolin / splicing benchmark, *Genome Biol* 2023 — https://link.springer.com/article/10.1186/s13059-023-03144-z
- CADD v1.7 release notes (Jan 2024) — https://cadd.gs.washington.edu/static/ReleaseNotes_CADD_v1.7.pdf; CADD-Splice, *Genome Med* 2021 — https://link.springer.com/article/10.1186/s13073-021-00835-9
- dbNSFP releases/changelog — https://www.dbnsfp.org/releases ; https://sites.google.com/site/jpopgen/dbNSFP/changelog
