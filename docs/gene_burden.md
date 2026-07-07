# Cross-Pedigree Gene-Burden Screening

Finds genes carrying more qualifying rare variants than a null expectation across the trio cohort, so recurrent gene-level signal can nominate candidates beyond single-family analysis.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- **Primary signal is de novo enrichment**, not case-vs-control burden. Test observed DNM counts per gene against a **Samocha-2014** trinucleotide mutation-rate expectation (Poisson), using **denovolyzeR**. It needs only trios and is robust to cohort batch effects.
- **The non-joint per-trio design forbids an internal cohort allele frequency.** Absent genotypes are ambiguous (no-call vs hom-ref), so there is no valid internal AC/AN. This is the single fact that shapes every method choice here.
- **Corroborative signal is TRAPD** vs **ancestry-matched, coverage-intersected gnomAD v4.1 exomes** — feasible from summary counts, but stratification/coverage-sensitive, so it is secondary.
- **Calibrate with synonymous variants: expect λ ≈ 1.0.** Report the synonymous-class λ and per-gene expected-vs-observed as standing QC for both the de novo and TRAPD arms.
- **Rarity gate for qualifying variants uses gnomAD v4.1 grpmax `faf95`** (< 1e-4 for de novo candidates), never the point-estimate popmax AF.
- **Exome-wide significance P < 2.5e-6** (~0.05 / 20,000 protein-coding genes); **BH q < 0.05** for candidate discovery. Combine masks/tiers per gene via **ACAT/Cauchy** to one p-value before thresholding.
- **Rank nominated genes by constraint** — a gene tolerant of damaging variation is not interesting even if nominally enriched.
- **Defer SAIGE-GENE+ / regenie** until a true joint call set exists; both require a joint genotype matrix this pipeline does not have.

## Why the non-joint design constrains the method

These VCFs are GATK Genotype-Refinement output, refined **per family** (CalculateGenotypePosteriors CGP posteriors), and **not jointly genotyped across the cohort**. Two consequences follow:

1. **No cohort allele-frequency table.** Without a joint call set you cannot compute a defensible internal AF. An absent genotype is ambiguous — it may be a genuine hom-ref or an uncalled site with no coverage. Internal recurrence is therefore usable only as an **artifact/blocklist signal**, never as a population frequency. (See [allele_frequency.md](allele_frequency.md) and [cohort_construction.md](cohort_construction.md).)
2. **Two orthogonal, design-appropriate signal sources remain:**
   - **De novo enrichment** vs a per-gene mutation-rate model — needs only trios, insensitive to cohort-level batch effects. **This is the primary signal.**
   - **Case-vs-external-control burden** (TRAPD-style vs gnomAD) — powerful but stratification- and coverage-sensitive; **corroborative only.**

## Primary signal: de novo enrichment vs a mutation model

### The Samocha framework

Samocha et al. (*Nat Genet* 2014) derive per-gene, per-consequence-class expected DNM rates from a **trinucleotide-context mutation model** (calibrated on human–chimp intergenic divergence). The observed DNM count for a gene and class is tested against its Poisson expectation, scaled by `2 × N_trios` (two transmissible haplotypes per trio). This is the foundation for the downstream tools below.

### Tooling

| Tool | Method | Fit for this pipeline |
| --- | --- | --- |
| **denovolyzeR** (R) | Per-gene / gene-set Poisson enrichment for LoF, missense, synonymous; CCDS-level exome-wide test | **Default.** Simplest, well-suited to a modest trio cohort. |
| DeNovoWEST (Kaplanis 2020) | Unified severity-weighted simulation test adding missense-clustering | Higher power on large NDD cohorts; heavier inputs. |
| extTADA / TADA (Nguyen 2017) | Bayesian integration of de novo + case-control counts; per-gene BF/FDR | Useful only if combining both signal arms. |

Newer proteome-wide constraint models refine expectations further but require heavier inputs and are not defaults here.

### De novo variant qualification

DNMs feeding the burden test are drawn from the same genotype-QC gates as the rest of the pipeline (see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)). A qualifying DNM must satisfy:

| Gate | Default |
| --- | --- |
| GATK screen | `hiConfDeNovo` (primary); `loConfDeNovo` = lower-sensitivity tier only |
| Child GQ (refined PP-derived) | ≥ 20 |
| DP (de novo) | ≥ 20 |
| Het allele balance | 0.25–0.75 |
| Parental cleanliness | each parent alt AD ≤ 1, DP ≥ 10 |
| gnomAD v4.1 frequency (grpmax `faf95`) | < 1e-4, absent-or-singleton, low `nhomalt` |
| Region | callable-region intersect |

**Known failure mode to guard against:** gnomAD priors in CalculateGenotypePosteriors can push a genuine ultra-rare pathogenic call toward hom-ref, suppressing a real de novo. For top candidates, cross-check the pre-refinement PL/GT before counting or discarding.

## Corroborative signal: case-vs-external-control burden (TRAPD)

**TRAPD** (Guo et al., *AJHG* 2018) counts qualifying-variant carriers per gene in cases and compares to gnomAD genotype/allele counts by Fisher/binomial. It needs only summary counts, so it is feasible **without joint genotyping**. Its correctness depends entirely on controlling four confounders:

- **Ancestry match** — compare cases to a specific gnomAD ancestry group (e.g. NFE, AFR), assigned by PCA, **not** to global popmax.
- **Coverage match** — restrict to sites where both cohorts have adequate depth (both ≥ 10–20×), using gnomAD's published per-base coverage.
- **Synonymous calibration** — synonymous variants should show **no** enrichment; a genome-wide synonymous burden **λ ≈ 1** validates the qualifying filters. Report it always.
- **Capture/platform differences** — gnomAD v4.1 exomes (730,947 exomes) vs the cohort's capture differ at indels and GC-extremes; treat those regions cautiously.

**RV-EXCALIBER** (*Nat Commun* 2021) supplies individual- and gene-level correction factors that explicitly de-bias gnomAD-as-control stratification and can be layered on TRAPD counts.

Because the non-joint design already denies a valid internal AF, treat TRAPD as **corroboration of de novo nominations**, not a standalone discovery engine.

## Variant qualification (masks)

A **mask** defines which variants qualify per gene. Consequence classes come from VEP on the Ensembl/GENCODE canonical/MANE transcript (see [functional_annotation.md](functional_annotation.md)); frequencies from the gnomAD v4.1 joint (exome+genome) oracle (see [allele_frequency.md](allele_frequency.md)).

| Mask | Definition |
| --- | --- |
| **M1 — pLoF** | LOFTEE **HC, no flags** |
| **M2 — pLoF + damaging missense** | M1 ∪ (REVEL ≥ 0.5 **or** AlphaMissense likely_pathogenic **or** CADD ≥ 20) |

AAF tiers: **≤ 1e-4** (de novo candidate default) and a stricter **≤ 1e-2 / 1e-3** tier reused from the recessive discovery gates. Multiple masks × AAF tiers improve power but must be paid for in correction (see below). All qualifying variants additionally require FILTER = PASS and the genotype-QC gates above.

> The rarity field is **grpmax `faf95`** (95% CI lower bound), not the point-estimate popmax/grpmax AF. The research brief phrased the burden filter on point-estimate popmax; the pipeline standard is faf95, consistent with every other frequency filter in this repo.

## Statistical test menu (reference)

These are the standard collapsing and variance-component tests. Only the Poisson (de novo) and Fisher (TRAPD) arms are active defaults here; the rest are documented for context and for the future joint-call-set case.

| Test | Character | When it wins |
| --- | --- | --- |
| **Simple burden** | Collapse qualifying variants to a count/indicator, regress on phenotype | All variants act in one direction |
| **CMC** (Li & Leal 2008) | MAF-binned collapse + Hotelling T² | Reduces df penalty of many rare sites |
| **SKAT / SKAT-O** | Variance-component (kernel); SKAT-O is a data-adaptive ρ mix of burden + SKAT | Bidirectional / neutral-diluted effects |
| **ACAT-V** | Cauchy combination of per-variant p-values | Few causal variants; fast, robust |
| **ACAT-O / GENE_P** | Cauchy omnibus over burden + SKAT-O + ACAT-V | Absorbs mask combination into one p-value |

### Biobank-scale frameworks — deferred

- **regenie** (v4.x; GENE_P omnibus; v4.1 custom variant weights) — whole-genome-regression + Firth/SPA, handles relatedness.
- **SAIGE-GENE+** (Zhou et al., *Nat Genet* 2022) — mixed-model set tests with correct type-I error at MAF ≤ 1e-4 / 1e-5, collapsing ultra-rare variants (MAC ≤ 10).

Both **require a joint genotype matrix** and are therefore **not applicable to unmerged per-trio VCFs**. Defer until a true cohort joint call set exists.

## Multiple-testing correction

- **Exome-wide gene-based threshold: P < 2.5e-6** (≈ 0.05 / ~20,000 protein-coding genes). The literature also expresses this as a class-specific Bonferroni; the canonical default here is the single ~2.5e-6 line.
- If **not** using a single omnibus p-value, divide further by the number of masks × AAF tiers × tests. Combining masks/tiers per gene via **ACAT/Cauchy (ACAT-O)** collapses them into one p-value and **avoids that penalty** — the preferred route.
- Report **FDR (BH q < 0.05)** for candidate discovery alongside the exome-wide line.
- Apply the same ~2.5e-6 threshold to the de novo Poisson tests.

## Ranking and interpretation

A statistically enriched gene is only interesting if it is **intolerant of the class of variation showing the excess**. Weight and rank nominated genes by constraint (see [gene_constraint.md](gene_constraint.md)):

- pLoF-mask nominations → gnomAD **v2.1.1** LOEUF (established; v4 flagged experimental), pLI, `s_het` (Zeng 2024) for short genes.
- A gene tolerant of damaging variation is **de-prioritized even when nominally enriched** — it is more likely an artifact than a discovery.
- Cross-reference a-priori gene lists and phenotype priors as **tiers, never hard filters** (never-drop rule; see [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md)), and germline cancer predisposition genes via [pediatric_cancer.md](pediatric_cancer.md).

## Pitfalls

- **Population stratification** — the dominant confounder in external-control burden; matched-ancestry gnomAD subsets + PCA sample assignment are mandatory for TRAPD.
- **Differential coverage** — without cohort joint genotyping, low-coverage no-calls masquerade as hom-ref; restrict to jointly-callable high-depth intervals.
- **Batch/platform effects** — capture kit, caller, and build details differ from gnomAD; the synonymous-λ is the calibration diagnostic that catches this.
- **No internal AF** — lean on de novo (per-family, batch-insensitive) as the primary defensible signal; treat gnomAD-control burden as corroborative.

## Scope limitations (stated honestly)

- **SNV/indel only initially.** This burden module counts SNV/indel qualifying variants. CNV/SV are a real blind spot — 10–15% of pediatric-cancer and rare-disease diagnoses are CNV/SV (single-exon *RB1*/*SMARCB1*/*DICER1*/*NF1* deletions, *PMS2* rearrangements). A gene whose recurrent hits are copy-number events will be missed until a future GATK-gCNV / Manta / ExomeDepth module exists.
- **Pseudogene / segmental-duplication genes** (*PMS2*/*PMS2CL*, *CYP21A2*, *SMN1/2*, *NEB*, *GBA*) yield low-confidence short-read calls; their per-gene counts are unreliable and those regions should be flagged, not trusted, in the burden tally.
- **Proband post-zygotic mosaicism** (*NF1*, overgrowth) produces low-VAF calls that fall outside the het AB 0.25–0.75 band and are dropped before the burden count — a source of false-negative de novos.
- **Calibration / truth sets** — the synonymous λ ≈ 1 diagnostic is the burden module's built-in calibration, but pipeline-wide sensitivity/precision still needs GIAB/CMRG truth sets and a positive-control variant panel.
- **Phenotype dependency** — phenotype-driven ranking assumes per-proband HPO terms, which are frequently sparse or absent in consortium data (see [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md)); burden nomination itself is phenotype-agnostic, which is a strength for discovery.

## Example: parameterized qualifying-variant extraction

Generic, path-free extraction of pLoF + rare qualifying variants into a per-gene tally. Substitute your own annotation field names and thresholds from `config/config.example.yaml`.

```bash
# 1. Keep PASS sites below the frequency gate, on canonical transcripts.
#    grpmax faf95 is carried in an INFO field annotated upstream (see allele_frequency.md).
bcftools view -f PASS -i 'INFO/gnomad_faf95_grpmax < 1e-4' "${TRIO_VCF}" -Oz -o "${RARE_VCF}"

# 2. Restrict to the pLoF mask (LOFTEE HC, no flags) via the VEP CSQ field.
#    filter_vep operates on the VEP-annotated CSQ; adjust field names to your cache release.
filter_vep -i "${RARE_VCF}" \
  --filter "LoF is HC and LoF_flags is not defined" \
  -o "${PLOF_VCF}" --force_overwrite

# 3. Intersect with the callable-region BED so no-calls do not inflate counts.
bcftools view -R "${CALLABLE_BED}" "${PLOF_VCF}" -Oz -o "${QUAL_VCF}"

# 4. Emit gene, consequence, and sample for the per-gene tally consumed by
#    denovolyzeR (de novo arm) or the TRAPD carrier-count step.
bcftools +split-vep "${QUAL_VCF}" -f '%CHROM\t%POS\t%SYMBOL\t%Consequence[\t%SAMPLE=%GT]\n' -d
```

```r
# De novo Poisson enrichment (primary arm), denovolyzeR.
# 'dnm' has columns: gene, class (lof/mis/syn); n_trios is the trio count.
library(denovolyzeR)
denovolyzeByGene(genes = dnm$gene, classes = dnm$class, nsamples = n_trios)
# QC gate: synonymous class enrichment ratio (observed/expected) must be ~1.0.
```

## Recommended defaults (this pipeline)

| Parameter | Default | Notes |
| --- | --- | --- |
| Primary signal | **De novo enrichment** (denovolyzeR, Samocha-2014 rates, Poisson) | Trio-only, batch-insensitive |
| Corroborative signal | **TRAPD** vs gnomAD v4.1 exomes | Ancestry-matched, coverage-intersected, synonymous-calibrated |
| Frequency oracle / field | gnomAD **v4.1** joint AF, grpmax **`faf95`** | Never internal cohort AC/AN |
| Rarity gate (de novo candidate) | faf95 **< 1e-4**, absent-or-singleton, low `nhomalt` | Per variant |
| Masks | **M1** = LOFTEE HC no-flag; **M2** = M1 ∪ (REVEL ≥ 0.5 or AlphaMissense LP or CADD ≥ 20) | AAF tiers ≤ 1e-4 and ≤ 1e-2/1e-3 |
| De novo genotype QC | GQ ≥ 20, DP ≥ 20, het AB 0.25–0.75, parent alt AD ≤ 1 | `hiConfDeNovo` primary screen |
| Exome-wide significance | **P < 2.5e-6** | ~0.05 / 20,000 genes |
| Discovery FDR | **BH q < 0.05** | Reported alongside |
| Mask/tier combination | **ACAT / Cauchy (ACAT-O)** → one p-value | Avoids mask-count penalty |
| Ranking weight | Constraint (LOEUF v2.1.1 / pLI / s_het) | Tolerant genes de-prioritized |
| Standing QC | Synonymous **λ ≈ 1.0**; per-gene expected-vs-observed | Both arms |
| Biobank frameworks | **Deferred** (SAIGE-GENE+, regenie) | Require a joint call set |

All values are configurable defaults in `config/config.example.yaml`, not immutable law. A gene-specific ClinGen VCEP threshold overrides any generic cutoff here.

## Sources

- Li & Leal, CMC — *AJHG* 2008: https://pubmed.ncbi.nlm.nih.gov/18691683/
- SKAT-O / ACAT-V (regenie GENE_P): https://rgcgithub.github.io/regenie/overview/ ; SBAT/NNLS — *AJHG* 2024, DOI 10.1016/j.ajhg.2024.09.009: https://www.cell.com/ajhg/fulltext/S0002-9297(24)00307-0
- SAIGE-GENE+ — Zhou et al., *Nat Genet* 2022, DOI 10.1038/s41588-022-01178-w: https://www.nature.com/articles/s41588-022-01178-w
- TRAPD — Guo et al., *AJHG* 2018, DOI 10.1016/j.ajhg.2018.08.016: https://www.sciencedirect.com/science/article/pii/S0002929718302842
- gnomAD v4.0 (807,162 individuals; 730,947 exomes; GRCh38): https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/
- gnomAD v4.1: https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/
- RV-EXCALIBER calibration — *Nat Commun* 2021, DOI 10.1038/s41467-021-26114-0: https://www.nature.com/articles/s41467-021-26114-0
- Samocha framework — *Nat Genet* 2014, DOI 10.1038/ng.3050: https://www.nature.com/articles/ng.3050
- DeNovoWEST — Kaplanis et al., *Nature* 2020, DOI 10.1038/s41586-020-2832-5: https://www.nature.com/articles/s41586-020-2832-5
- extTADA — Nguyen et al., *Genome Med* 2017, DOI 10.1186/s13073-017-0497-y: https://pmc.ncbi.nlm.nih.gov/articles/PMC5738153/
- GATK Genotype Refinement workflow (CGP de novo posteriors): https://gatk.broadinstitute.org/hc/en-us/articles/360035531432-Genotype-Refinement-workflow-for-germline-short-variants
- GMKF Kids First trio data (gVCF, GRCh38): https://ega-archive.org/studies/phs001228
