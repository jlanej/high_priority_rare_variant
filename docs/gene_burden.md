# Cross-Pedigree Gene Consolidation (Recurrence-Based)

Finds genes where rare, functional **inherited** variants recur across multiple independent individuals in the cohort, so recurrent gene-level signal can nominate candidates beyond single-family analysis.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- **Primary signal is recurrence-based gene consolidation** across the cohort, focused on **inherited germline variation**. For each gene, tally the number of **distinct individuals** carrying a qualifying variant under each inheritance model — **dominant** (rare inherited het), **biallelic** (homozygous + compound het), and **X-linked**. A gene is **recurrent** at **≥ `min_carriers` (default 2)** distinct individuals.
- **The key new signal is recurrence of inherited heterozygous variants.** A rare (grpmax `faf95` < 1e-4), functional, inherited het is only weakly interesting in one family, but becomes compelling when it **recurs across multiple individuals in the same gene**.
- **Rank recurrent genes first, weighted by gene constraint** (LOEUF / pLI / `s_het`). A recurrent het in a **haploinsufficient** gene is the most compelling result; a gene tolerant of damaging variation is de-prioritized even when recurrent.
- **De novo enrichment is an OPTIONAL SECONDARY signal only.** De novo filtering **and review are handled by separate dedicated machinery**; here de novo is carried as a lightweight cross-reference column (GATK `hiConfDeNovo`, child-membership checked). When a mutation-rate table is supplied, an optional Samocha-2014 Poisson enrichment (denovolyzeR-style) is reported — **exome-wide P < 2.5e-6, BH q < 0.05** — but it is not the driver.
- **The non-joint per-trio design forbids an internal cohort allele frequency.** Absent genotypes are ambiguous (no-call vs hom-ref), so there is no valid internal AC/AN. Internal recurrence here means *distinct-individual carrier counts of qualifying variants*, not a population frequency.
- **TRAPD** vs **ancestry-matched, coverage-intersected gnomAD v4.1 exomes** remains an **optional corroboration** (not yet implemented) — stratification/coverage-sensitive, never a standalone discovery engine.
- **Rarity gate for qualifying variants uses gnomAD v4.1 grpmax `faf95`** (< 1e-4 for dominant/de novo candidates), never the point-estimate popmax AF.
- **Defer SAIGE-GENE+ / regenie** until a true joint call set exists; both require a joint genotype matrix this pipeline does not have. mtDNA heteroplasmy is out of scope here (handled by a separate dedicated pipeline).

## Why the non-joint design constrains the method

These VCFs are GATK Genotype-Refinement output, refined **per family** (CalculateGenotypePosteriors CGP posteriors), and **not jointly genotyped across the cohort**. Two consequences follow:

1. **No cohort allele-frequency table.** Without a joint call set you cannot compute a defensible internal AF. An absent genotype is ambiguous — it may be a genuine hom-ref or an uncalled site with no coverage. So "recurrence" here is **not** a population frequency: it is a tally of **distinct individuals carrying a qualifying variant** in a gene, each independently vetted by the same genotype-QC and rarity gates. (See [allele_frequency.md](allele_frequency.md) and [cohort_construction.md](cohort_construction.md).)
2. **Design-appropriate signal sources:**
   - **Recurrence-based gene consolidation** — distinct-individual carrier counts per gene, per inheritance model, ranked recurrent-first and weighted by constraint. **This is the primary signal**, and it centres on **inherited germline variation** (dominant het, biallelic, X-linked).
   - **De novo Poisson enrichment** vs a per-gene mutation-rate model — needs only trios, insensitive to cohort-level batch effects, but **optional and secondary** here (de novo filtering/review live in separate machinery).
   - **Case-vs-external-control burden** (TRAPD-style vs gnomAD) — powerful but stratification- and coverage-sensitive; **optional corroboration, not yet implemented.**

## Primary signal: recurrence-based gene consolidation

### Distinct-individual carrier counts per gene, by inheritance model

Step 6 aggregates the per-family candidate calls (from the inheritance screen, see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)) into a per-gene tally. For each gene it counts the number of **distinct individuals** carrying a qualifying variant under each model:

| Model | What is counted | Rarity gate |
| --- | --- | --- |
| **Dominant (inherited het)** | Rare, functional **heterozygous** variant transmitted from ≥ 1 parent (parent-of-origin recorded: maternal / paternal / both), **not** part of a compound-het pair | `faf95` < 1e-4 |
| **Biallelic** | Homozygous **or** compound het (two rare hets, same gene, in *trans*) | `faf95` < 1e-2 (permissive) / 1e-3 (high-confidence) per allele |
| **X-linked recessive** | Affected male (hemizygous) + carrier mother (father's chrX not required); or affected female (hom-alt, carrier mother, hemizygous father); sex-aware ploidy | `faf95` < 1e-2 (permissive) per allele |
| **De novo (secondary)** | GATK `hiConfDeNovo`, child-membership checked (`annotations.is_hiconf_denovo_for`) | counted in a **separate** column, not part of the recurrence driver |

The dominant-het count is the **key new signal**: individually a rare inherited het is weak evidence, but a gene that accumulates such hets across **multiple distinct individuals** is a strong nomination.

### Recurrence flag and ranking

- A gene is **recurrent** when its distinct-individual carrier count reaches **`min_carriers` (default 2)**.
- Genes are ranked **recurrent-first**, then **distinct-variant** recurrence above same-variant (founder/artifact), then by the strongest per-model recurrence p, then **constraint-weighted** (constrained = LOEUF < 0.35 **or** pLI ≥ 0.9 **or** s_het ≥ 0.1 **or** pHaplo ≥ 0.86), then by carrier / dominant-het counts. A recurrent het in a **haploinsufficient** gene ranks highest.
- De novo carrier counts and the optional de novo enrichment p-value are carried as **secondary columns** used only to break ties after the recurrence and constraint keys.

## Optional secondary signal: de novo enrichment vs a mutation model

> **De novo filtering and review are handled by separate dedicated machinery.** In this pipeline de novo is retained only as a **lightweight cross-reference** (detected via GATK `hiConfDeNovo`, child-membership checked via `annotations.is_hiconf_denovo_for`). The enrichment test below runs **only when a mutation-rate table is supplied** and is never the primary nomination signal.

### The Samocha framework

Samocha et al. (*Nat Genet* 2014) derive per-gene, per-consequence-class expected DNM rates from a **trinucleotide-context mutation model** (calibrated on human–chimp intergenic divergence). The observed DNM count for a gene and class is tested against its Poisson expectation, scaled by `2 × N_trios` (two transmissible haplotypes per trio). This is the foundation for the downstream tools below.

### Tooling

| Tool | Method | Fit for this pipeline |
| --- | --- | --- |
| **denovolyzeR** (R) | Per-gene / gene-set Poisson enrichment for LoF, missense, synonymous; CCDS-level exome-wide test | **Default for the optional de novo arm.** Simplest, well-suited to a modest trio cohort. |
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

## Optional corroboration: case-vs-external-control burden (TRAPD)

> **Not yet implemented.** TRAPD is documented here as an optional corroboration of recurrent-gene nominations, not an active default.

**TRAPD** (Guo et al., *AJHG* 2018) counts qualifying-variant carriers per gene in cases and compares to gnomAD genotype/allele counts by Fisher/binomial. It needs only summary counts, so it is feasible **without joint genotyping**. Its correctness depends entirely on controlling four confounders:

- **Ancestry match** — compare cases to a specific gnomAD ancestry group (e.g. NFE, AFR), assigned by PCA, **not** to global popmax.
- **Coverage match** — restrict to sites where both cohorts have adequate depth (both ≥ 10–20×), using gnomAD's published per-base coverage.
- **Synonymous calibration** — synonymous variants should show **no** enrichment; a genome-wide synonymous burden **λ ≈ 1** validates the qualifying filters. Report it always.
- **Capture/platform differences** — gnomAD v4.1 exomes (730,947 exomes) vs the cohort's capture differ at indels and GC-extremes; treat those regions cautiously.

**RV-EXCALIBER** (*Nat Commun* 2021) supplies individual- and gene-level correction factors that explicitly de-bias gnomAD-as-control stratification and can be layered on TRAPD counts.

Because the non-joint design already denies a valid internal AF, treat TRAPD as **optional corroboration of recurrent-gene nominations**, not a standalone discovery engine.

## Variant qualification (masks)

A **mask** defines which variants qualify per gene. Consequence classes come from VEP on the Ensembl/GENCODE canonical/MANE transcript (see [functional_annotation.md](functional_annotation.md)); frequencies from the gnomAD v4.1 joint (exome+genome) oracle (see [allele_frequency.md](allele_frequency.md)).

| Mask | Definition |
| --- | --- |
| **M1 — pLoF** | LOFTEE **HC, no flags** |
| **M2 — pLoF + damaging missense** | M1 ∪ (REVEL ≥ 0.5 **or** AlphaMissense likely_pathogenic **or** CADD ≥ 20) |

AAF tiers: **≤ 1e-4** (dominant / de novo candidate default) and a **≤ 1e-2 / 1e-3** tier reused from the recessive (biallelic) discovery gates. Multiple masks × AAF tiers improve power but must be paid for in correction where the optional statistical tests are used (see below). All qualifying variants additionally require FILTER = PASS and the genotype-QC gates above.

> The rarity field is **grpmax `faf95`** (95% CI lower bound), not the point-estimate popmax/grpmax AF. The research brief phrased the burden filter on point-estimate popmax; the pipeline standard is faf95, consistent with every other frequency filter in this repo.

## Statistical test menu (reference)

These are the standard collapsing and variance-component tests, documented for context and for the future joint-call-set case. The active primary method here is **recurrence-based gene consolidation** (distinct-individual carrier counts, not a formal association test); the optional Poisson (de novo) and Fisher (TRAPD) arms are secondary/corroborative. The remaining tests require a joint genotype matrix this pipeline does not have.

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

The primary recurrence signal is now **calibrated**: for each gene, observed distinct-individual carriers are tested against `Binomial(N_trios, p)` with `p = 1 − Π_v (1 − faf95_v)²` over the gene's qualifying inherited variants (an allele absent from gnomAD is floored at `absent_faf95_floor`, default 1e-6). This yields a per-gene `p_recurrence`, a **BH `q_recurrence`**, and an exome-wide flag — so 2 carriers of a *private* variant are genome-wide significant while 2 carriers of a common-ish variant are not. It is a **case-only approximation** using in-cohort variants; the gnomAD-derived per-gene cumulative-AF version (TRAPD/CoCoRV) is the planned upgrade. The same thresholds apply to the **optional secondary de novo enrichment**:

- **Exome-wide gene-based threshold: P < 2.5e-6** (≈ 0.05 / ~20,000 protein-coding genes). The literature also expresses this as a class-specific Bonferroni; the canonical default here is the single ~2.5e-6 line.
- If **not** using a single omnibus p-value, divide further by the number of masks × AAF tiers × tests. Combining masks/tiers per gene via **ACAT/Cauchy (ACAT-O)** collapses them into one p-value and **avoids that penalty** — the preferred route.
- Report **FDR (BH q < 0.05)** for candidate discovery alongside the exome-wide line.
- Apply the same ~2.5e-6 threshold to the de novo Poisson tests.

## Ranking and interpretation

A recurrent gene is only interesting if it is **intolerant of the class of variation showing the recurrence**. Weight and rank nominated genes by constraint (see [gene_constraint.md](gene_constraint.md)):

- pLoF-mask nominations → gnomAD **v2.1.1** LOEUF (established; v4 flagged experimental), pLI, `s_het` (Zeng 2024) for short genes. In code, "constrained" = **LOEUF < 0.35 or pLI ≥ 0.9**.
- A recurrent **dominant het in a haploinsufficient gene is the most compelling** result; a gene tolerant of damaging variation is **de-prioritized even when recurrent** — it is more likely an artifact than a discovery.
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
