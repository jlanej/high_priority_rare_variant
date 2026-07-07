# Clinical Pathogenicity: ClinVar & ACMG/AMP

How this pipeline turns per-trio VCF variants into clinical pathogenicity evidence using ClinVar review-gated assertions and a points-based ACMG/AMP classification.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- **ClinVar** is consumed from a **dated, version-pinned release**; auto-promote **Pathogenic/Likely_pathogenic (P/LP) at ≥ 2★** (no conflicts); 1★ P/LP → prioritize + human review; `Conflicting_classifications` and `Uncertain_significance` (VUS) → flag, never auto-promote; **0★ excluded** from auto-logic.
- ClinVar star mapping: **4★** practice guideline, **3★** expert panel (a ClinGen VCEP call overrides all submitters), **2★** multiple submitters/no conflicts, **1★** single submitter or conflicting, **0★** no assertion criteria.
- **ACMG classifier backbone = AutoGVP** (CHOP/Kids First; ClinVar + modified InterVar with graded PVS1 and PP5/BP6 removed) — purpose-built for this GMKF pediatric-cancer / rare-disease GRCh38 use case.
- **Combining = Tavtigian/ClinGen Bayesian points**: **P ≥ 10, LP 6–9, VUS 0–5, LB −1…−5, B ≤ −6** (Supporting ±1, Moderate ±2, Strong ±4, Very Strong ±8).
- **PM2 at Supporting strength only** — it is ACMG *evidence*, distinct from the pipeline's rarity screening gate; do not conflate "passed the rarity filter" with "PM2 met."
- **PVS1** is applied at **graded strength** via the Abou Tayoun 2018 decision tree, gated on ClinGen gene–disease validity **≥ Moderate** and a known loss-of-function mechanism — never naive full-strength.
- **Missense in-silico:** report **one** Pejaver-2022–calibrated predictor per variant (primary **REVEL**), never stack correlated predictors — see [functional_annotation.md](functional_annotation.md).
- **Gene-level gate:** restrict high-priority auto-calls to ClinGen **Definitive / Strong / Moderate** gene–disease validity; a gene-specific VCEP threshold overrides any generic cutoff.

## ClinVar as an evidence source

### Record model and review status

ClinVar aggregates submitted records (SCVs) into variant-level (VCV) and variant/condition (RCV) records, and assigns each a **gold-star review status (0–4)** reflecting the strength of the submitting evidence. Since January 2024 ClinVar splits classification into three axes — **germline**, **somatic clinical impact**, and **oncogenicity** — each with its own VCF INFO tags (`CLNSIG`/`CLNREVSTAT`, `ONC*`/`ONCREVSTAT`, `SCI*`/`SCIREVSTAT`). This pipeline uses the **germline** axis for rare-disease and germline pediatric-cancer screening.

| Stars | `CLNREVSTAT` token | Meaning |
|------:|--------------------|---------|
| 4★ | `practice_guideline` | Professional practice guideline |
| 3★ | `reviewed_by_expert_panel` | ClinGen VCEP / expert-panel call — **overrides all submitters** |
| 2★ | `criteria_provided,_multiple_submitters,_no_conflicts` | Concordant multi-lab |
| 1★ | `criteria_provided,_single_submitter` | Single lab, criteria provided |
| 1★ | `criteria_provided,_conflicting_classifications` | Conflicting (inspect `CLNSIGCONF`) |
| 0★ | `no_assertion_criteria_provided` / `no_classification_provided` | No criteria — excluded from auto-logic |

`CLNSIG` VCF tokens use underscores: `Pathogenic`, `Likely_pathogenic`, `Pathogenic/Likely_pathogenic`, `Uncertain_significance`, `Likely_benign`, `Benign`, and `Conflicting_classifications_of_pathogenicity` (formerly `Conflicting_interpretations_of_pathogenicity`), plus low-penetrance / risk-allele and non-standard terms (`drug_response`, `association`, `protective`, `Affects`).

### Consumption rules (defaults)

- **Extract P/LP** by matching `Pathogenic` / `Likely_pathogenic` (and the combined `Pathogenic/Likely_pathogenic`) in `CLNSIG`.
- **Auto-promote at ≥ 2★** with no conflicts. **1★ P/LP** is prioritized but routed to human review.
- Treat `Conflicting_classifications_of_pathogenicity` and `Uncertain_significance` as **not** evidence of pathogenicity — **never drop them**; flag for review and inspect the `CLNSIGCONF` submitter breakdown.
- **Exclude 0★** from auto-logic. gnomAD's own guidance is to filter ClinVar to **≥ 1★ with a specified classification** to avoid 0★ noise.
- ClinVar is not a substitute for evidence: it can carry outdated or single-lab calls, and it **reclassifies continuously**. Pin a **dated release**, record the release date in run provenance, and plan **periodic re-classification** of previously reported VUS.

```bash
# Annotate a per-trio VCF with a version-pinned ClinVar release (GRCh38).
# clinvar_YYYYMMDD.vcf.gz is the dated release; record the date in run metadata.
bcftools annotate \
  -a "${CLINVAR_VCF}" \
  -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNSIGCONF \
  -Oz -o "${OUT_VCF}" "${IN_VCF}"

# High-confidence auto-promotable P/LP: P or LP in CLNSIG AND >=2-star, no conflicts.
bcftools view -i \
  'INFO/CLNSIG ~ "Pathogenic" &&
   INFO/CLNREVSTAT ~ "multiple_submitters" &&
   INFO/CLNREVSTAT !~ "conflicting"' \
  "${OUT_VCF}"
```

## ACMG/AMP framework and ClinGen SVI refinements

### The 2015 framework, points-based

The 2015 ACMG/AMP framework (Richards et al.) defines 28 criteria (16 pathogenic PVS1–PP5, 12 benign) combined by verbal rules into P / LP / VUS / LB / B. ClinGen's Sequence Variant Interpretation (SVI) working group has progressively replaced the verbal combining rules with the **Tavtigian Bayesian points system**, in which each criterion contributes exponentially-scaled points and the sum determines the class. This lets a criterion be applied at a **tunable strength** rather than a fixed weight.

| Strength | Pathogenic points | Benign points |
|----------|------------------:|--------------:|
| Supporting | +1 | −1 |
| Moderate | +2 | −2 |
| Strong | +4 | −4 |
| Very Strong | +8 | — |

| Total points | Classification |
|-------------:|----------------|
| ≥ 10 | Pathogenic |
| 6 – 9 | Likely Pathogenic |
| 0 – 5 | Uncertain Significance (VUS) |
| −1 … −5 | Likely Benign |
| ≤ −6 | Benign |

### Key SVI refinements applied here

- **PM2 → Supporting by default** (SVI Recommendation v1.0). Absence/rarity is weak evidence; do not apply PM2 at Moderate. In this pipeline PM2 is ACMG *evidence* and is deliberately kept distinct from the upstream rarity **screening** gate (see [allele_frequency.md](allele_frequency.md)) — passing the rarity filter is not the same as "PM2 met."
- **PP3 / BP4 calibration (Pejaver 2022)**: continuous predictors receive strength-stratified thresholds. Use **one** predictor per variant; do not sum correlated predictors as independent evidence. The **REVEL** primary defaults are restated below and detailed in [functional_annotation.md](functional_annotation.md). AlphaMissense is *not* part of the Pejaver 2022 calibration and lacks an SVI PP3 stratification — treat it only as orthogonal support.
- **PVS1 decision tree (Abou Tayoun 2018)**: loss-of-function variants receive graded strength (PVS1 / _Strong / _Moderate / _Supporting) by consequence, NMD escape, exon/region context (last exon, 3′-terminal 50 bp, single-exon), and only when LoF is the disease mechanism and gene–disease validity is **≥ Moderate**. Naive full-strength PVS1 is a major over-calling source. See [functional_annotation.md](functional_annotation.md) and [gene_constraint.md](gene_constraint.md).
- **PS3 / BS3 functional (Brnich 2020)**: assay-based strength is set via OddsPath from the number of validated controls (≥ 11 controls → Moderate; more → Strong).

No single consolidated 2023–2025 replacement guideline is published yet; refinements continue to arrive as individual SVI recommendations and as gene-specific VCEP specifications, which take precedence over the generic rules whenever they exist.

## Automated ACMG classifiers

Rule-based ACMG automation varies widely in rigor; the common failure is applying **PVS1 at full strength** and using the reputable-source criteria PP5/BP6, both of which inflate P/LP calls.

| Tool | Basis | Caveat |
|------|-------|--------|
| **AutoGVP** (CHOP / Kids First / NCI) | ClinVar + **modified InterVar** (graded PVS1, PP5/BP6 removed), dockerized R workflow | Purpose-built for this GMKF pediatric-cancer / rare-disease GRCh38 use case — the reference backbone adopted here |
| InterVar (Wang lab) | Rule-based on ANNOVAR | Applies **PVS1 at full strength**, uses PP5/BP6 → P/LP over-calling; no VCEP specs |
| TAPES, GeneBe | Rule-based, no phenotype integration | Lower causal-variant prioritization in benchmarks |
| Franklin (Genoox) | Proprietary | Strong benchmarks but black-box; ToS constraints for bulk/container use |

**This pipeline adopts AutoGVP** as the classifier backbone: it integrates a dated ClinVar release with a modified InterVar (PVS1-strength adjustment, PP5/BP6 removed) and is built for exactly this consortium use case. Universal caveat: automated calls are **screening aids, not diagnostic** — VUS and conflicts require human review, and none replace applicable VCEP rules.

## ClinGen gene–disease validity gate

ClinGen classifies each gene→disease relationship as **Definitive / Strong / Moderate / Limited / Disputed / Refuted / No Known Disease Relationship**. ACMG recommends diagnostic panels include only **Definitive / Strong / Moderate** genes. This pipeline uses gene–disease validity as a **gene-level gate before variant scoring**: it restricts high-priority auto-calls to Definitive/Strong/Moderate genes and also gates PVS1 applicability (PVS1 requires ≥ Moderate validity plus a known LoF mechanism). Note this is a gate on *auto-promotion*, applied consistently with the pipeline's never-drop principle — curated lists elsewhere act as priors, not hard filters (see [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md)).

## Integration pattern (per-trio VCF)

```text
Per-trio VCF (GRCh38, GATK genotype-refined)
  └─ VEP: canonical / MANE consequence            (functional_annotation.md)
     └─ join dated ClinVar VCF: CLNSIG, CLNREVSTAT, CLNSIGCONF
        └─ join gnomAD v4.1 joint grpmax faf95     (allele_frequency.md)
           └─ REVEL (single predictor)             (functional_annotation.md)
              └─ AutoGVP: ClinVar + modified InterVar
                 └─ Tavtigian points  (PVS1 graded per Abou Tayoun; PM2 → Supporting)
                    └─ gate on ClinGen gene–disease validity ≥ Moderate
                       └─ germline VCEP specs override generic cutoffs
```

Because Kids First trios are **GATK genotype-refinement output (posterior/PP-refined), not cohort joint-genotyped**, this pipeline does **not** derive internal-cohort allele frequency. Population rarity comes from external gnomAD v4.1 (see [allele_frequency.md](allele_frequency.md)); de-novo confidence comes from the trio PP/GQ, not cohort frequency (see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)). Germline pediatric-cancer specifics — dominant vs recessive CPS handling, second-hit boosts, PMS2/PMS2CL — are in [pediatric_cancer.md](pediatric_cancer.md).

## Recommended defaults (this pipeline)

| Parameter | Default | Notes |
|-----------|---------|-------|
| ClinVar release | Dated, version-pinned; recorded in provenance | ClinVar reclassifies continuously — re-classify prior VUS periodically |
| Auto-promote P/LP | `CLNSIG` P or LP **and ≥ 2★**, no conflicts | 1★ P/LP → prioritize + human review |
| Conflicting / VUS | Flag, never auto-promote; inspect `CLNSIGCONF` | Never-drop |
| 0★ records | Excluded from auto-logic | gnomAD guidance: keep ≥ 1★ with a specified classification |
| ACMG classifier | **AutoGVP** (ClinVar + modified InterVar, PP5/BP6 removed) | Graded PVS1 |
| Combining rule | Tavtigian points: **P ≥ 10, LP 6–9, VUS 0–5, LB −1…−5, B ≤ −6** | Supporting ±1 / Moderate ±2 / Strong ±4 / Very Strong ±8 |
| PM2 | **Supporting** strength only | ACMG evidence, distinct from the rarity screening gate |
| PVS1 | Graded (Abou Tayoun 2018); gene validity ≥ Moderate + known LoF | Never naive full-strength |
| PP3/BP4 missense | **REVEL** (single predictor): PP3 supporting ≥ 0.644, moderate ≥ 0.773, strong ≥ 0.932; BP4 supporting ≤ 0.290, moderate ≤ 0.183 | Do not stack; literature ranges from Pejaver 2022 |
| Gene gate | ClinGen **Definitive / Strong / Moderate** validity | VCEP threshold overrides generic cutoff |

All values are **configurable defaults** (`config/config.example.yaml`), not immutable law. A gene-specific ClinGen VCEP specification (its own PVS1 strength, PM2/BA1/BS1 thresholds, or PP3 calibration) **overrides** any generic default here.

## Scope limitations (stated honestly)

- **SNV/indel only initially.** ClinVar/ACMG logic here classifies SNVs and indels. CNV/SV pathogenic calls (e.g. single-exon *RB1*/*SMARCB1*/*DICER1*/*NF1* deletions, *PMS2* rearrangements — 10–15% of pediatric-cancer and rare-disease diagnoses) are a known blind spot pending a future GATK-gCNV / Manta / ExomeDepth module.
- **Pseudogene / segmental-duplication genes** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) are low-confidence from short reads; ClinVar/ACMG calls in those regions must be flagged regardless of star level.
- **Proband post-zygotic mosaicism** (e.g. *NF1*, overgrowth): low-VAF calls fall outside the standard heterozygous allele-balance band and can be missed before classification is even reached (see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)).
- **Genotype-refinement caveat:** gnomAD priors in GATK `CalculateGenotypePosteriors` can push a genuine ultra-rare pathogenic call toward hom-ref — for top ClinVar/ACMG candidates, cross-check the pre-refinement PL/GT.
- **Calibration/validation:** classification sensitivity/precision should be measured against GIAB/CMRG truth sets and a positive-control variant panel; no truth-set benchmark is wired in yet.

## Sources

- ClinVar review status & stars: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
- ClinVar classification / `CLNSIG`: https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
- gnomAD ClinVar review-status filter guidance: https://gnomad.broadinstitute.org/news/2023-09-clinvar-variants-filter-by-review-status/
- gnomAD v4.1: https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/ ; ClinGen guidance on gnomAD v4 for VCEPs (Mar 2024): https://clinicalgenome.org/site/assets/files/9445/clingen_guidance_to_vceps_regarding_the_use_of_gnomad_v4_march_2024.pdf
- Tavtigian points (Bayesian ACMG): https://onlinelibrary.wiley.com/doi/10.1002/humu.24088 (DOI 10.1002/humu.24088)
- ACGS 2023 UK variant-classification guidelines (points / PM2): https://www.acgs.uk.com/media/12443/uk-practice-guidelines-for-variant-classification-v1-2023.pdf
- Pejaver 2022 PP3/BP4 calibration: https://www.cell.com/ajhg/pdfExtended/S0002-9297(22)00461-X (DOI 10.1016/j.ajhg.2022.10.013)
- PVS1 decision tree, Abou Tayoun 2018: https://onlinelibrary.wiley.com/doi/abs/10.1002/humu.23626 (DOI 10.1002/humu.23626)
- PS3/BS3 functional-assay strength, Brnich 2020: https://link.springer.com/article/10.1186/s13073-019-0690-2 (DOI 10.1186/s13073-019-0690-2)
- ClinGen gene–disease validity framework: https://clinicalgenome.org/docs/evaluating-the-clinical-validity-of-gene-disease-associations-an-evidence-based-framework-developed-by-the-clinical-genome/
- AutoGVP (CHOP / Kids First / NCI): https://academic.oup.com/bioinformatics/article/40/3/btae114/7616989 (DOI 10.1093/bioinformatics/btae114)
- Automated ACMG classifier benchmark: https://academic.oup.com/bioinformatics/article/42/2/btaf623/8483023
