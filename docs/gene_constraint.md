# Gene & Region Constraint for Candidate Weighting

How this pipeline uses gene- and region-level constraint metrics as graded priors to rank — never to exclude — candidate variants from per-trio GRCh38 screening.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- Constraint quantifies depletion of variation versus a mutation-rate-calibrated neutral expectation; it is a gene/region **prior on dominant (haploinsufficiency) relevance** and is largely uninformative for recessive genes.
- **Constraint is a ranking weight, never a standalone exclusion filter.** It is incomplete for short genes, novel genes, and recessive genes.
- Established pLoF-intolerance defaults: **gnomAD v2.1.1 LOEUF < 0.35** (or **v4.0 LOEUF < 0.6**, flagged experimental) and **pLI ≥ 0.9** (unchanged v2→v4).
- gnomAD **v4 constraint is officially "experimental"**; v2.1.1 remains the recommendation for established/clinical weighting.
- For **short genes** (LOEUF/pLoF counts underpowered, roughly the shortest quartile), prefer **s_het (Zeng 2024) ≥ 0.1**.
- Haploinsufficiency support: **pHaplo ≥ 0.86** (Collins 2022) and **ClinGen dosage HI = 3** (established).
- Regional missense: **MPC ≥ 2** strongly up-weights a missense candidate; gene-level **missense Z > 3.09** is supporting.
- **Do NOT** down-weight recessive (compound-het / homozygous) candidates using pLoF constraint — biallelic loss is rare in the population, so a tolerant pLI/LOEUF says nothing against recessive pathogenicity.

## Why constraint matters more in this design

The VCFs screened here are GMKF Kids First **per-trio** call sets (GRCh38, GATK Genotype-Refinement output) that are **not jointly genotyped across the cohort**. Internal cohort allele counts are therefore uninterpretable as population frequency (see [allele_frequency.md](allele_frequency.md)), and there is no internal cohort-level burden signal at the single-trio stage. External gene-level priors — of which constraint is the strongest and most portable — consequently carry **more** of the prioritization weight. That makes disciplined use essential: constraint informs ranking and tiering, and is always reported alongside the ClinGen dosage curation, but it never drops a variant.

## What constraint measures

For each gene, an expected count of variants (synonymous / missense / pLoF) is derived from a sequence-context (trinucleotide) mutation-rate model corrected for local coverage. Observed counts come from very rare variants (AF < 0.1%) in a reference population. **Constraint = the observed/expected (o/e) deficit**: fewer observed rare damaging variants than expected implies purifying selection, i.e. the gene is intolerant of that class of variation.

Key consequence: constraint reflects selection against the *heterozygous* effect of a damaging allele. It is directly informative for **dominant / haploinsufficient** genes and largely silent for **recessive** genes, where a single loss-of-function allele is tolerated in the population.

## Metric families and thresholds

### gnomAD constraint (pLI, LOEUF, missense o/e, missense Z)

| Metric | Meaning | Strong-constraint default | Notes |
|---|---|---|---|
| **pLI** | Probability a gene is intolerant of pLoF (dichotomous classifier) | **pLI ≥ 0.9** (tolerant ≤ 0.1) | Recommended cutoff unchanged across v2 and v4 |
| **LOEUF** | Upper bound of the 90% CI on the pLoF o/e ratio (continuous, conservative; lower = more constrained) | **v2.1.1 < 0.35**; **v4.0 < 0.6** | Preferred over pLI for ranking/weighting; thresholds NOT interchangeable across releases |
| **Missense Z** | Missense constraint z-score (higher = more constrained) | **Z > 3.09** (≈ top 15%) | Gene-level support for missense candidates |
| **Missense o/e** | Missense observed/expected ratio | graded | Continuous down-weight for missense in tolerant genes |

**v2 vs v4 — do not mix thresholds.** gnomAD v2.1.1 (125,748 exomes, GRCh37) and v4.0 (730,947 exomes, GRCh38, high-coverage bases, autosomes) differ ~6× in sample size, which shifts the LOEUF distribution. gnomAD explicitly labels **v4 constraint "experimental" and recommends v2.1.1 for established/clinical use**. Because this pipeline is GRCh38, the most defensible current stance is to annotate genes with **v2.1.1 LOEUF/pLI by gene symbol / Ensembl gene ID** (established) and optionally also carry **v4.0 LOEUF (flagged experimental)** using its own 0.6 cutoff. Do not apply the v2 (0.35) threshold to a v4 value or vice versa.

### Regional missense constraint: MPC

Whole-gene missense metrics miss **sub-genic** regional depletion — roughly 15% of genes harbor regions locally intolerant of missense variation even when the gene overall is not. **MPC** (Missense badness, PolyPhen-2, Constraint; Samocha et al. 2017) scores this per variant on a 0–5 range:

| MPC | Interpretation | Action |
|---|---|---|
| **≥ 2** | Strong regional missense constraint; enriched for pathogenic de novo missense | Up-weight the missense candidate (per canonical default) |
| ~1–2 | Moderate | Modest up-weight |
| **< 1** | No case enrichment observed | No constraint-based up-weight |

MPC is applied **per variant** (it is region-specific), complementing the per-gene missense Z. A v4-based regional missense constraint (RMC) release exists and can substitute once validated for this workflow.

### Selection-based metrics (s_het)

**s_het** is the fitness reduction in heterozygous carriers of a pLoF allele — a mechanistically interpretable dominant-constraint measure, high in Mendelian-dominant and essential genes. Estimation has been refined over successive models (Cassa 2017 on ExAC / gnomAD v2.1.1; Weghorn 2019; **Zeng 2024 GeneBayes**, a Bayesian model combining population genetics with gene features).

The key practical advantage: **s_het (Zeng 2024) outperforms LOEUF for short genes**, where pLoF counts are underpowered (roughly the shortest quartile of genes). This pipeline uses **s_het ≥ 0.1** as the strong-constraint default and prefers it over LOEUF/pLI when a gene is short. The exact 0.1 cutoff is a community convention (≳ 0.15 flags very strong selection), not a single authoritative number; treat it as a configurable default.

### Haploinsufficiency / dosage-sensitivity prediction

These predict dosage sensitivity from features orthogonal to variant counts, and are especially useful for genes where population variant counts are sparse:

| Score / source | Meaning | Established threshold |
|---|---|---|
| **pHaplo** (Collins 2022, Cell; 18,641 genes) | ML haploinsufficiency probability | **≥ 0.86** |
| **pTriplo** (Collins 2022) | ML triplosensitivity probability | **≥ 0.84** |
| **ClinGen dosage HI** | Curated haploinsufficiency evidence | **HI = 3** (established) is the strongest curated evidence |
| **Episcore** (Han 2018) / **EDS** | Epigenomic-feature HI prediction | orthogonal corroboration |

ClinGen dosage curation (`HI` / `TS`) is authoritative where present and is reported alongside any predicted score.

### Other gene-level scores (secondary)

- **DOMINO** (Quinodoz 2017): ML probability that a gene causes *dominant* disease; "very likely dominant" > 0.8, "very likely recessive" < 0.2. Useful to gate dominant-model weighting.
- **GeVIR / RVIS**: older percentile intolerance ranks, largely **superseded by LOEUF/s_het**; use only as secondary corroboration.

## Interpreting constraint by inheritance model

- **Dominant / haploinsufficient genes:** constraint is directly informative. Low LOEUF, high pLI, high s_het, high pHaplo, or ClinGen HI = 3 → strong prior for pLoF pathogenicity. A high-confidence pLoF (LOFTEE HC, no flags — see [functional_annotation.md](functional_annotation.md)) in a constrained gene is a top nomination.
- **Recessive genes:** most are **NOT** LoF-constrained (biallelic loss is rare in the population), so a tolerant pLI / high LOEUF does **not** argue against pathogenicity. **Do not down-weight recessive (compound-het / homozygous) candidates using pLoF constraint.** Missense metrics (MPC, missense Z) retain some value; rely otherwise on ClinGen recessive curation and variant-level evidence (see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)).
- **Pediatric cancer genes:** defer to curated lists (ClinGen HI = 3, COSMIC Cancer Gene Census germline, established predisposition genes — see [pediatric_cancer.md](pediatric_cancer.md)) over generic constraint, since several are recessive or tumor-suppressor genes with variable constraint.

## How constraint feeds candidate weighting

Constraint is combined with the a-priori gene tiers ([gene_lists_and_phenotype.md](gene_lists_and_phenotype.md)) as a graded, monotonic weight — never a hard gate.

**Dominant single-hit pLoF candidates — priority tiers (defaults):**

| Tier | Condition |
|---|---|
| Tier 1 | LOEUF_v2 < 0.35 (or v4 < 0.6) **OR** pLI ≥ 0.9 **OR** s_het ≥ 0.1 **OR** pHaplo ≥ 0.86 **OR** ClinGen HI = 3 |
| Tier 2 | 0.35 ≤ LOEUF_v2 < 0.6 |
| Tier 3 (down-weight) | LOEUF_v2 ≥ 0.6 **AND** s_het < 0.05 **AND** DOMINO < 0.2 |

**Missense candidates:** up-weight **MPC ≥ 2** (strong) / ≥ 1 (moderate); treat **missense Z > 3.09** as gene-level support.

**Short genes (~ shortest quartile):** prefer **s_het (Zeng 2024)** over LOEUF/pLI, which are underpowered.

**Recessive (biallelic) candidates:** apply no pLoF-constraint down-weighting.

Constraint also **ranks/weights genes nominated by cross-pedigree burden** — a gene tolerant of damaging variation is not an interesting burden hit (see [gene_burden.md](gene_burden.md)).

### Practical annotation

Constraint scores are gene- (or region-) level tables joined to the annotated VCF by gene symbol / Ensembl ID (per-variant for MPC). A typical join, keeping every record and only *adding* a weight column, avoids any exclusion:

```bash
# Annotate an already VEP-annotated VCF with gene-level constraint columns.
# Constraint tables are joined by gene/transcript ID; NO records are dropped.
vcfanno constraint.conf annotated.vcf.gz | bgzip > annotated.constraint.vcf.gz

# MPC is applied per variant (region-specific), not per gene:
bcftools annotate \
  -a mpc_scores.grch38.bed.gz \
  -c CHROM,FROM,TO,MPC \
  -h mpc_header.hdr \
  annotated.constraint.vcf.gz -Oz -o annotated.mpc.vcf.gz
```

Downstream weighting (tier assignment) is done in the scoring layer, not by filtering the VCF, so that a novel or short gene with weak constraint is retained at a lower prior rather than removed.

## Scope limitations (state honestly)

- **Constraint weights SNVs/indels only.** The pHaplo/pTriplo and ClinGen dosage material here is used to *weight* SNV/indel candidates, **not** to drive CNV/SV detection. CNV/SV are a real blind spot for this pipeline (10–15% of pediatric-cancer and rare-disease diagnoses are CNV/SV — e.g. single-exon *RB1*/*SMARCB1*/*DICER1*/*NF1* deletions, *PMS2* rearrangements); a future GATK-gCNV / Manta / ExomeDepth module is where dosage-sensitivity scores would properly drive dosage-change calling.
- **Constraint is incomplete for short genes and novel genes** — hence the s_het preference for short genes and the never-exclude rule so novel-gene discovery survives.
- **Recessive and tumor-suppressor genes** are systematically under-served by pLoF constraint; defer to curated recessive/predisposition lists.
- **No positive-control calibration of the weighting itself** is asserted here; sensitivity/precision of the prioritization logic should ultimately be measured against truth sets (GIAB/CMRG) and a positive-control variant panel.

## Recommended defaults (this pipeline)

| Parameter | Default | Notes |
|---|---|---|
| Primary pLoF-intolerance source | **gnomAD v2.1.1 LOEUF** | Established; v4.0 flagged experimental |
| LOEUF strong-constraint cutoff | **v2.1.1 < 0.35** (or **v4.0 < 0.6**) | Not interchangeable across releases |
| pLI cutoff | **≥ 0.9** | Unchanged v2 → v4 |
| Short-gene metric | **s_het (Zeng 2024) ≥ 0.1** | Preferred where LOEUF/pLoF counts are underpowered |
| Haploinsufficiency | **pHaplo ≥ 0.86**; **ClinGen HI = 3** | pTriplo ≥ 0.84 for triplosensitivity |
| Regional missense | **MPC ≥ 2** (strong) up-weight | Per variant |
| Gene-level missense | **missense Z > 3.09** support | — |
| Recessive candidates | **no pLoF-constraint down-weighting** | Biallelic loss is population-tolerated |
| Role of constraint | **ranking weight only** | Never a standalone exclusion filter |

All values are configurable defaults in `config/config.example.yaml`; a gene-specific ClinGen VCEP threshold overrides any generic cutoff where available.

## Sources

- gnomAD v4.0 gene constraint (thresholds, v2/v4 differences, "experimental" caveat): https://gnomad.broadinstitute.org/news/2024-03-gnomad-v4-0-gene-constraint/
- Karczewski et al. 2020, gnomAD constraint (pLI / LOEUF / missense Z): https://www.nature.com/articles/s41586-020-2308-7
- Samocha et al. 2017, MPC / regional missense constraint (preprint): https://www.biorxiv.org/content/10.1101/148353v1 ; code: https://github.com/broadinstitute/regional_missense_constraint
- Cassa et al. 2017, s_het (Nat Genet): https://www.nature.com/articles/ng.3831 ; scores: http://genetics.bwh.harvard.edu/genescores/selection.html
- Weghorn et al. 2019, s_het model (Nat Genet): https://www.nature.com/articles/s41588-018-0291-9
- Zeng et al. 2024, GeneBayes s_het (Nat Genet): https://www.nature.com/articles/s41588-024-01820-9 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC10245655/
- Quinodoz et al. 2017, DOMINO (AJHG): https://www.cell.com/ajhg/references/S0002-9297(17)30368-3
- Han et al. 2018, Episcore (Nat Commun): https://www.nature.com/articles/s41467-018-04552-7
- Collins et al. 2022, pHaplo / pTriplo dosage sensitivity (Cell): https://www.cell.com/cell/fulltext/S0092-8674(22)00788-7
- ClinGen Dosage Sensitivity curation: https://clinicalgenome.org/curation-activities/dosage-sensitivity/
