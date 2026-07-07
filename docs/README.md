# Methods Reference

State-of-the-art, source-cited reference for screening **GMKF Kids First per-trio VCFs**
(GRCh38, GATK Genotype-Refinement output, **not** jointly genotyped across the cohort) for
high-priority rare variants in **rare disease and germline pediatric cancer**.

These documents are the *why* behind the pipeline. The *what* (the concrete, ordered steps and
the file each produces) is in **[pipeline_design.md](pipeline_design.md)**.

## How to read this

Start with **[pipeline_design.md](pipeline_design.md)** for the vetted end-to-end flow, then dip
into the topical references as needed. Every threshold in every document is the configurable
default defined in the **[Canonical defaults](#canonical-defaults)** table below — that table is
the single source of truth; if a document ever disagrees with it, the table wins.

| Document | Covers |
|----------|--------|
| [pipeline_design.md](pipeline_design.md) | Vetted end-to-end flow; critique of the original 5-step proposal; data artifacts; scope limits |
| [cohort_construction.md](cohort_construction.md) | Why a naive `merge` of non-joint trios corrupts AC/AN; the site-only union recipe |
| [allele_frequency.md](allele_frequency.md) | gnomAD v4.1, grpmax `faf95`, maximum credible AF, why external AF not internal |
| [functional_annotation.md](functional_annotation.md) | VEP/IMPACT, LOFTEE, calibrated missense/splice predictors, dbNSFP |
| [clinical_classification.md](clinical_classification.md) | ClinVar review status, ACMG/AMP Bayesian points, AutoGVP |
| [gene_constraint.md](gene_constraint.md) | LOEUF / pLI / s_het / MPC / pHaplo as candidate-weighting priors |
| [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md) | Genotype-refinement outputs, per-mode trio logic, genotype & sample QC, trio tools |
| [pediatric_cancer.md](pediatric_cancer.md) | Germline predisposition prevalence, genes, two-hit model, cancer gene lists |
| [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md) | OMIM/PanelApp/ClinGen/COSMIC/ACMG-SF, HPO/Exomiser priors, tiering |
| [gene_burden.md](gene_burden.md) | Recurrence-based gene consolidation (dominant het + biallelic across individuals); de novo enrichment as secondary |
| [tooling_and_reproducibility.md](tooling_and_reproducibility.md) | Container/conda-lock, GHCR CI, Apptainer, PHI-safe repo |

---

## Canonical defaults

This is the **single source of truth** for every threshold in the pipeline. All values are
**configurable defaults** (see [`config/config.example.yaml`](../config/config.example.yaml)),
not immutable law. A gene-specific ClinGen VCEP value **overrides** any generic cutoff here.

### Frequency oracle
- **gnomAD v4.1** (GRCh38; 730,947 exomes + 76,215 genomes), **joint** (exome+genome) AF/AN.
- Filter field = **grpmax `faf95`** (filtering allele frequency, 95% CI lower bound), *not* the
  point-estimate popmax AF. Fall back to grpmax AF only when `faf95` is unavailable.
- **Never** use internal cohort AC/AN as population frequency; internal recurrence is valid only
  as an artifact/blocklist signal. Heed the v4.1 exome/genome **discordance flag**.

### Rarity gates (grpmax `faf95`) — a screening gate, distinct from the ACMG **PM2** criterion
| Mode | Keep candidate if `faf95` < | Notes |
|------|----------------------------|-------|
| Dominant / de novo | **1e-4** | de novo additionally requires absent-or-singleton in gnomAD, low `nhomalt` |
| Recessive / comp-het | **1e-2** per allele (permissive); **1e-3** high-confidence tier | applied per variant, not per gene |
| Benign, all modes | drop if `faf95` ≥ **0.05** (ClinGen BA1) | never rescue |

PM2 is applied at **Supporting** strength only and is *evidence*, not the rarity gate itself.

### Functional / in-silico (report **one** predictor per variant — never stack correlated tools)
| Signal | Tier / cutoff |
|--------|---------------|
| pLoF | LOFTEE **HC, no flags**; grade PVS1 by Abou-Tayoun tree (NMD-escape/last-exon/single-exon downgrade), gated on ClinGen validity ≥ Moderate + known LoF mechanism |
| Missense (primary **REVEL**, Pejaver-2022) | PP3 supporting ≥ 0.644, moderate ≥ 0.773, strong ≥ 0.932; BP4 supporting ≤ 0.290, moderate ≤ 0.183 |
| Missense (orthogonal **AlphaMissense**) | likely_pathogenic ≥ 0.564; ambiguous 0.34–0.564; likely_benign ≤ 0.34 |
| Splicing (**SpliceAI** masked, Walker-2023) | PP3 Δ ≥ 0.2; BP4 Δ ≤ 0.1; 0.1–0.2 uninformative; canonical ±1,2 with Δ ≥ 0.5 = high tier |
| Regional missense | MPC ≥ 2 up-weights; missense Z > 3.09 gene-level support |
| Benign deprioritize (BP4) | REVEL ≤ 0.183 **and** AlphaMissense ≤ 0.34 **and** SpliceAI Δ ≤ 0.1 |

### Clinical evidence
- **ClinVar**: pin a dated release; auto-promote **P/LP at ≥ 2★** (no conflicts); 1★ → prioritize
  + human review; Conflicting/VUS → flag, never auto-promote; exclude 0★ from auto-logic.
- Classifier backbone = **AutoGVP** (CHOP/Kids First). Combining = **Tavtigian/ClinGen points**
  (P ≥ 10, LP 6–9, VUS 0–5), **PM2 at Supporting**.

### Gene constraint — a **ranking weight, never a standalone exclusion filter**
- gnomAD **v2.1.1** LOEUF (established) primary; pLI ≥ 0.9 / LOEUF_v2 < 0.35 (v4 < 0.6, flagged
  experimental). Prefer **s_het (Zeng 2024) ≥ 0.1** for short genes. pHaplo ≥ 0.86 / ClinGen HI = 3.
- **Do not** down-weight recessive candidates by pLoF constraint.

### Inheritance models (Step 5) & genotype QC (GATK-refined trios)
- Trust **refined `PP`-derived GQ**. GQ ≥ 20; DP ≥ 10; het AB 0.25–0.75; hom-alt AB ≥ 0.90;
  hom-ref AB ≤ 0.10 (AB from AD); FILTER = PASS only.
- **Dominant** (inherited): rare (`faf95 < 1e-4`), functional **het** transmitted from ≥ 1 parent
  (origin recorded) — the recurrence signal Step 6 consolidates.
- **Recessive**: homozygous, or **compound het** = two rare hets, same gene, in **trans**
  (parent-of-origin; WhatsHap fallback).
- **X-linked recessive**: male hemizygous + carrier mother; sex-aware ploidy, drop male non-PAR
  het calls; kid sex inferred from chrX heterozygosity when the PED is unknown.
- **De novo** (SECONDARY / cross-reference only — filtering & review handled by separate
  machinery): `hiConfDeNovo` (child-membership checked) → re-verify DP/AB + parental cleanliness.
- **Sample QC**: Peddy/somalier (sex, relatedness); Mendelian-error < 2%.
- **Failure mode**: gnomAD priors in CalculateGenotypePosteriors can suppress genuine ultra-rare
  pathogenic calls — cross-check pre-refinement `PL` for top candidates.

### Pediatric cancer
- Version-pinned gene-list **union** (prior/tier, not hard filter): ACMG **SF v3.3** (84 genes)
  cancer subset ∪ **PanelApp GE green** (Childhood solid tumours panel 243, Adult susceptibility
  245) ∪ curated recessive CPS set (MLH1/MSH2/MSH6/PMS2, FANC\*, BRCA2, ATM, DIS3L2, BLM).
- Dominant → report het P/LP; recessive → require **biallelic**; de novo in dominant CPS = top
  tier; second hit (LOH / DICER1 hotspot) = tier **boost**, never a filter. **PMS2** needs
  pseudogene-aware handling.

### Cross-pedigree gene consolidation (recurrence across individuals)
- Tally **distinct individuals** per gene by model: **dominant** (qualifying rare functional het),
  **biallelic** (hom / comp-het), **X-linked**; de novo counted separately (secondary).
- A gene is **recurrent** at ≥ **min_carriers** (default **2**) distinct individuals; rank
  recurrent-first, **weighted by constraint** (LOEUF / pLI / s_het — a recurrent het in a
  haploinsufficient gene is the most compelling).
- OPTIONAL secondary: de novo Poisson enrichment vs the Samocha model (exome-wide **P < 2.5e-6**,
  BH **q < 0.05**) when a mutation-rate table is supplied.

### A-priori gene lists & phenotype — priors/tiers, never hard include/exclude (**never-drop rule**)
- Tier 1 known gene → lenient thresholds; Tier 2 strong candidate (constraint/expression); Tier 3
  novel → retained at lower prior. Rarity/impact/QC gating applied *before and independently of*
  list priors. Phenotype: **Exomiser** + LIRICAL as ranking priors (not hard gates); HPO per
  proband is a dependency, degrade gracefully when sparse.

### Reproducibility / tooling
- One image `FROM mambaorg/micromamba@sha256:...` + **conda-lock**. Pin bcftools/htslib/samtools
  **1.22**, bedtools **2.31.1**, vcfanno **0.3.3**, slivar **0.3.4**, **ensembl-vep 115**, Python
  cyvcf2/pysam/pandas/numpy. VEP cache **external, release-matched**. CI → **GHCR** per commit
  (buildx, provenance + SBOM, tag by SHA, amd64). Apptainer: real-disk `TMPDIR`, **no
  `--containall`**, `--cleanenv`. **No hard paths / no PHI / dbGaP-safe.**

### Scope boundaries & known limitations
**De novo** filtering/review and **mtDNA heteroplasmy** are handled by **separate dedicated
pipelines** (de novo is a cross-reference here; mtDNA is out of scope). Within scope: SNV/indel
only (CNV/SV is a real blind spot); pseudogene/seg-dup regions flagged low-confidence; phenotype
(Exomiser/HPO) prior planned; validate with GIAB/CMRG truth sets + positive controls. See
[pipeline_design.md](pipeline_design.md#known-scope-limitations-stated-honestly-not-hidden).

---

*Reference currency: gnomAD v4.1, VEP 115/gnomAD-v4.1-built-in from r113, ACMG SF v3.3 (2025),
gnomAD constraint v2.1.1 (v4 experimental), Exomiser 15.1.0 (Java 21). Every document carries its
own dated `## Sources`.*
