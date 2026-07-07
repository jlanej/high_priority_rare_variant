# Germline Pediatric Cancer Predisposition

How this pipeline recognizes and tiers germline cancer-predisposition-syndrome (CPS) variants in GMKF Kids First per-trio VCFs.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- **~8–10% of pediatric cancers carry a P/LP germline CPS variant** (Zhang 2015: 8.5%; Gröbner 2018: 7.6%), and **family history does NOT predict it** — pedigree cancer history is not a valid pre-filter.
- **Gene list is a version-pinned union used as a PRIOR/TIER, never a hard include/exclude**: ACMG SF v3.3 (84 genes) cancer subset ∪ PanelApp GE green childhood panels ∪ a curated recessive CPS set. The never-drop rule preserves novel-gene discovery.
- **Rarity default = grpmax `faf95`** (not point-estimate popmax AF): dominant/de novo `< 1e-4`; recessive per-allele `< 1e-2` (permissive), `< 1e-3` high-confidence tier; `≥ 0.05` (BA1) → drop.
- **Zygosity by mechanism**: dominant CPS → report **het** P/LP; recessive CPS (CMMRD, Fanconi, AT, DIS3L2, BLM) → require **biallelic** (hom or trans compound-het via trio phasing).
- **De novo** P/LP in a dominant CPS gene = **top tier** (trio design is a strength here); use GATK `hiConfDeNovo` then re-verify with DP/AB + parental-cleanliness.
- **Second hit** (LOH, biallelic somatic loss, or a gene-specific hotspot like the DICER1 RNase-IIIb codons) = **tier boost, never a filter requirement** — Kids First per-trio VCFs are germline-only.
- **PMS2** requires pseudogene(PMS2CL)-aware handling; short-read calls in PMS2 exons 11–15 are low-confidence.
- **Known blind spot: SNV/indel only.** 10–15% of CPS diagnoses are CNV/SV (single-exon RB1/SMARCB1/DICER1/NF1 deletions, PMS2 rearrangements) — out of initial scope.

## Why germline screening in pediatric cancer

Roughly one in ten children with cancer carries a pathogenic or likely-pathogenic (P/LP) germline variant in an established predisposition gene — far above the historical family-history-driven expectation:

- **Zhang et al. (NEJM 2015)** — P/LP variants in **8.5%** of 1,120 patients under 20 years across 60 autosomal-dominant CPS genes. Critically, **family history did not predict predisposition** in most cases. You cannot filter on pedigree cancer history alone.
- **Gröbner et al. (Nature 2018)** — **7.6%** germline P/LP in a childhood-cancer cohort, enriched in DNA-repair genes (mismatch-repair MSH2/MSH6/PMS2; double-strand-break repair TP53, BRCA2, CHEK2).

Predisposition rate is strongly **tumor-type dependent** — adrenocortical carcinoma (~50%), hypodiploid B-ALL (~28%), then high-grade glioma / medulloblastoma / retinoblastoma (~15–25%). If a reliable tumor-type diagnosis is available in Kids First metadata, it can drive **tumor-type-aware prior weighting**; the pipeline does not require it and never uses it as a hard filter.

**De novo dominant** events contribute materially (notably TP53, RB1, RASopathies). The trio design is a direct strength: a P/LP variant absent in both parents but present in the proband is simultaneously a strong pathogenicity signal (de novo) and a surveillance-relevant finding. See [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md) for the de novo calling and verification logic.

## Major genes and syndromes

Dominant with incomplete penetrance unless noted. Names below are representative, not exhaustive; the operative list is the version-pinned union assembled in [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md).

### Dominant / high-yield

| Gene(s) | Syndrome / phenotype | Notes for tiering |
|---|---|---|
| **TP53** | Li-Fraumeni (ACC, choroid plexus carcinoma, osteosarcoma, RMS, early breast) | Highest-yield single gene; strong de novo contribution |
| **RB1** | Retinoblastoma (classic two-hit) | Bilateral disease ≈ germline |
| **DICER1** | Pleuropulmonary blastoma, cystic nephroma, Sertoli-Leydig, thyroid | Germline LOF **+ somatic RNase-IIIb hotspot** is the signature (see below) |
| **APC** | FAP → hepatoblastoma, medulloblastoma (Turcot) | LOF mechanism |
| **NF1 / NF2** | Neurofibromatosis | NF1 also relevant to proband mosaicism (see limitations) |
| **WT1** | Wilms tumor | — |
| **SMARCB1 / SMARCA4** | Rhabdoid tumor / ATRT | Frequently CNV/deletion — see blind-spot note |
| **SUFU / PTCH1** | Gorlin / medulloblastoma | — |
| **PTEN, VHL, RET, SDHx, TSC1/2** | ACMG SF secondary-findings cancer core | Always-report-if-P/LP tier |
| **RASopathy set** (PTPN11, SOS1, KRAS, NRAS, RAF1, CBL, NF1) | Noonan → JMML | Mixed de novo / inherited |
| **PALB2, BRCA1/2, CHEK2, ATM (mono)** | Moderate-to-high adult-onset risk; some pediatric relevance | ATM monoallelic = moderate |

### Recessive / biallelic

| Gene(s) | Syndrome | Requirement |
|---|---|---|
| **MLH1, MSH2, MSH6, PMS2** (biallelic) | **CMMRD** — constitutional mismatch-repair deficiency; most aggressive CPS, early multi-organ tumors | Biallelic; PMS2 most common and **pseudogene-aware calling required** |
| **FANC\*, BRCA2/FANCD1** | Fanconi anemia | Biallelic |
| **ATM** (biallelic) | Ataxia-telangiectasia | Biallelic |
| **DIS3L2** | Perlman / Wilms | Biallelic |
| **BLM** | Bloom syndrome | Biallelic |

## Two-hit (Knudson) model and filtering implications

Classic tumor suppressors need **two hits**: a germline first hit plus a somatic second hit (loss of heterozygosity, a second somatic SNV, or focal deletion). Kids First is a **germline-only per-trio VCF**, so the somatic hit cannot be confirmed from these inputs. The pipeline therefore:

- **Prioritizes germline LOF/truncating and known-P/LP variants in tumor-suppressor genes even at the heterozygous state** — one germline hit is the actionable germline finding for a dominant CPS gene.
- Treats a documented somatic second hit (LOH / biallelic loss where matched tumor data exist elsewhere in Kids First, or a gene-specific hotspot such as the **DICER1 RNase-IIIb** codons E1705/D1709/E1788/D1810/E1813) as a **+1 tier boost, never a filter requirement**.

### Zygosity logic by inheritance model

- **Dominant CPS genes** → report a het P/LP variant.
- **Recessive CPS genes** (CMMRD, Fanconi, AT, DIS3L2, BLM) → require **biallelic**: homozygous, or compound-heterozygous in **trans** established by trio phasing (parent-of-origin: one allele mat-only, one pat-only) or read-backed phasing. A de novo second hit is a valid partner allele.
- Recessive-gene **hets are not discarded** — they are retained at a lower tier in case a second allele (including a CNV, currently out of scope) is present.

The genotype QC, de novo verification, and compound-het phasing mechanics live in [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md).

## Gene-panel sources to intersect

The pipeline assembles a **version-pinned union** and uses it as a prior/tier, never as a hard filter (the never-drop rule of [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md)):

- **ACMG SF v3.3 (2025)** — 84 genes; the cancer subset (TP53, RB1, APC, MMR genes, BRCA1/2, PALB2, PTEN, VHL, RET, SDHx, WT1, TSC1/2, NF2, SMARCB1, etc.) is the reportable-secondary-findings core and the "always-report-if-P/LP" tier. (v3.2 (2023) added CALM1-3; v3.3 added the non-cancer genes ABCD1/CYP27A1/PLN.)
- **Genomics England PanelApp** childhood panels — the "Childhood solid tumours cancer susceptibility" panel and haematological/childhood panels, filtered to **green** (diagnostic-grade) genes; **pin the panel `version`** in the container.
- **Curated recessive CPS set** — {MLH1, MSH2, MSH6, PMS2, FANC\*, BRCA2, ATM, DIS3L2, BLM}.
- **Jongmans et al. (2016) selection tool** and the **COG/AACR surveillance consensus** — used to annotate actionability and drive interpretation/reporting logic, not variant filtering.

Store the panel version and download date with the container image (see [tooling_and_reproducibility.md](tooling_and_reproducibility.md)).

## Rarity, consequence, and clinical evidence

Rarity and impact gating is applied **before and independently of** the CPS gene priors. All frequency filtering uses the shared frequency oracle in [allele_frequency.md](allele_frequency.md): **gnomAD v4.1** (GRCh38; 730,947 exomes + 76,215 genomes), joint (exome+genome), with **grpmax `faf95`** (the 95% CI lower bound of the group-max filtering allele frequency) as the filter field — not the point-estimate popmax/grpmax AF. Never use internal cohort AC/AN as a population frequency; the non-joint per-trio merge makes AN uninterpretable.

| Layer | Default | Notes |
|---|---|---|
| Dominant / de novo rarity | grpmax faf95 `< 1e-4` | For de novo, additionally require absent-or-singleton in gnomAD and low `nhomalt` |
| Recessive per-allele rarity | grpmax faf95 `< 1e-2` (permissive); `< 1e-3` high-confidence tier | Applied **per variant**, not per gene (literature ranges 1e-3 to 1e-2 across sources) |
| Hard benign (all modes) | grpmax faf95 `≥ 0.05` (BA1) | Drop, never rescue |
| Consequence | LOFTEE HC (no flags) pLoF; missense with ClinVar P/LP or calibrated damaging in-silico | See functional layer below |
| ClinVar | Auto-promote P/LP at **≥ 2★** (no conflicts) | 1★ P/LP → prioritize + human review; VUS/Conflicting → flag, never auto-promote |

A gene-specific ClinGen VCEP BA1/BS1 value (or a Whiffin/Ware maximum credible AF) **overrides** these generic cutoffs whenever available.

For functional consequence and calibrated in-silico impact (REVEL Pejaver-2022 tiers, AlphaMissense, SpliceAI, graded PVS1), see [functional_annotation.md](functional_annotation.md). For ClinVar handling and the AutoGVP-backed ACMG/AMP classification, see [clinical_classification.md](clinical_classification.md). Note that **PM2 is applied at Supporting strength as ACMG evidence** and is distinct from the rarity screening gate — passing the rarity filter is not the same as "PM2 met."

### PVS1 for LOF-mechanism CPS genes

Most high-yield dominant CPS genes (TP53, RB1, APC, the MMR genes, NF1, WT1, SMARCB1, PALB2, BRCA1/2) act by **loss of function**, so a novel LOFTEE-HC frameshift/stop/canonical-splice variant is a strong candidate via **PVS1**. Grade PVS1 by the Abou Tayoun 2018 decision tree (NMD-escape, last-exon, 3′-terminal-50bp, single-exon downgrades), gated on ClinGen gene–disease validity ≥ Moderate and a known LoF mechanism. Do not apply PVS1 in genes whose disease mechanism is not LOF (e.g., activating RASopathy alleles).

## Trio-aware candidate strategy

1. Rarity + impact + genotype QC gate the variant (independently of gene priors).
2. Restrict interpretation weighting to the version-pinned union gene list as a **prior/tier**.
3. ClinVar P/LP at ≥ 2★ → auto-tier; a novel LOFTEE-HC pLoF in an LOF-mechanism gene → high tier via graded PVS1.
4. Flag **de novo** (GATK `hiConfDeNovo`, re-verified) and resolve **compound-het phase** for recessive genes.
5. Apply zygosity logic: dominant → het reportable; recessive → biallelic required.
6. Apply the second-hit boost where a somatic LOH/biallelic loss or gene-specific hotspot is documented.

## Recommended defaults (this pipeline)

| Parameter | Default | Override source |
|---|---|---|
| **Gene list** | ACMG SF v3.3 cancer subset ∪ PanelApp GE green childhood panels (pin `version`) ∪ recessive CPS set {MLH1, MSH2, MSH6, PMS2, FANC\*, BRCA2, ATM, DIS3L2, BLM} — used as PRIOR/TIER, never hard include/exclude | `config/config.example.yaml`; store panel version + download date in container |
| **Population AF field** | gnomAD v4.1 joint grpmax `faf95` | [allele_frequency.md](allele_frequency.md) |
| **Dominant / de novo rarity** | faf95 `< 1e-4`; de novo also absent-or-singleton + low `nhomalt` | ClinGen VCEP / max credible AF |
| **Recessive rarity** | faf95 `< 1e-2` per allele (permissive); `< 1e-3` high-confidence | per-variant, not per-gene |
| **Hard benign** | faf95 `≥ 0.05` (BA1) → drop | — |
| **Consequence** | LOFTEE HC no-flag pLoF; missense with ClinVar P/LP or calibrated damaging in-silico; graded PVS1 only in LOF-mechanism genes | [functional_annotation.md](functional_annotation.md) |
| **ClinVar** | Auto-report P/LP at ≥ 2★, no conflicts (dated, pinned release) | [clinical_classification.md](clinical_classification.md) |
| **Genotype QC** | Refined GQ ≥ 20; DP ≥ 10 (≥ 20 for de novo); het AB 0.25–0.75; FILTER = PASS | [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md) |
| **De novo** | GATK `hiConfDeNovo` screen → re-verify DP/AB + parental cleanliness (each parent alt AD ≤ 1, DP ≥ 10) + gnomAD absent/singleton; top tier in dominant CPS genes | — |
| **Zygosity** | Dominant → het P/LP reportable; recessive → biallelic (hom or trans compound-het by trio phasing) | — |
| **Second-hit boost** | +1 tier if matched somatic LOH/biallelic loss or gene-specific hotspot (e.g., DICER1 RNase-IIIb E1705/D1709/E1788/D1810/E1813); never required | — |
| **PMS2** | Pseudogene(PMS2CL)-aware calling/annotation; flag PMS2 exons 11–15 low-confidence | — |

### Example: restrict a rarity-filtered, QC-passed VCF to the CPS union list for weighting

These commands are illustrative and parameterized — no real paths. The gene-list BED (`${CPS_UNION_BED}`) is the version-pinned union used here as a prior/annotation, not as a discovery hard filter.

```bash
# Annotate variants overlapping the CPS union region (prior/tier flag, not a drop)
bcftools annotate \
  --annotations "${CPS_UNION_BED}.gz" \
  --columns CHROM,FROM,TO,CPS_GENE \
  --header-lines "${CPS_HEADER}" \
  --output-type z --output "${OUT_VCF}" \
  "${QC_PASSED_VCF}"

# Flag (do NOT exclude) low-confidence PMS2 pseudogene region for human review
bcftools annotate \
  --annotations "${PMS2CL_LOWCONF_BED}.gz" \
  --mark-sites +PMS2_LOWCONF \
  --output-type z --output "${OUT_VCF2}" \
  "${OUT_VCF}"
```

## Known scope limitations

State these honestly — the pipeline does not pretend to cover them.

- **SNV/indel only initially.** CNV/SV are a real blind spot: **10–15% of pediatric-cancer predisposition diagnoses are CNVs/SVs** — single-exon RB1, SMARCB1, DICER1, and NF1 deletions, and PMS2 rearrangements — which an SNV/indel-only pipeline will systematically miss. A future module would add GATK-gCNV / Manta / ExomeDepth-type calling. ClinGen dosage/pHaplo material ([gene_constraint.md](gene_constraint.md)) is used to weight SNVs, not to drive CNV detection.
- **PMS2 / pseudogene-aware calling is flagged, not fully solved.** The PMS2/PMS2CL paralogy makes short-read PMS2 calls (and other segmental-duplication genes — CYP21A2, SMN1/2, NEB, GBA) low-confidence. The current mitigation is to **flag** those regions for human review; robust masking / paralog-aware calling is future work.
- **Proband post-zygotic mosaicism** (relevant for NF1 and some overgrowth/CPS phenotypes) produces low-VAF calls that fall outside the het AB 0.25–0.75 band and are filtered out — a dedicated mosaic tier is not yet implemented.
- **Calibration/validation.** There is no CPS-specific truth set wired in yet; the pipeline-wide plan (GIAB/CMRG truth sets, synonymous λ ≈ 1, a positive-control variant panel) applies here too.
- **Phenotype/HPO dependency.** Tumor-type-aware and phenotype-driven weighting depends on per-proband HPO/diagnosis annotation, which is frequently sparse or absent in consortium data; the pipeline must degrade gracefully when it is missing.

## Sources

- Zhang J et al. Germline Mutations in Predisposition Genes in Pediatric Cancer. NEJM 2015;373:2336. https://www.nejm.org/doi/full/10.1056/NEJMoa1508054
- Gröbner SN et al. The landscape of genomic alterations across childhood cancers. Nature 2018;555:321. https://www.nature.com/articles/nature25480
- ACMG SF v3.3 (2025) policy statement, Genetics in Medicine. https://www.gimjournal.org/article/S1098-3600(25)00101-7/fulltext | PubMed https://pubmed.ncbi.nlm.nih.gov/40568962/
- ACMG SF v3.2 (2023). https://pmc.ncbi.nlm.nih.gov/articles/PMC10524344/
- Jongmans MCJ et al. Recognition of genetic predisposition in pediatric cancer patients: an easy-to-use selection tool. Eur J Med Genet 2016. https://pubmed.ncbi.nlm.nih.gov/26825391/
- Genomics England PanelApp (Childhood solid tumours cancer susceptibility). https://panelapp.genomicsengland.co.uk/
- Tabori U et al. Clinical Management and Tumor Surveillance Recommendations of Inherited MMR Deficiency in Childhood (CMMRD, AACR consensus). Clin Cancer Res 2017. https://aacrjournals.org/clincancerres/article/23/11/e32/79879/ | ERN GENTURIS CMMRD guidelines https://pmc.ncbi.nlm.nih.gov/articles/PMC11607302/
- GATK Genotype Refinement workflow / CalculateGenotypePosteriors. https://gatk.broadinstitute.org/hc/en-us/articles/360035531432-Genotype-Refinement-workflow-for-germline-short-variants | https://gatk.broadinstitute.org/hc/en-us/articles/4409897225243-CalculateGenotypePosteriors
- Selection criteria for assembling a pediatric CPS gene panel. Fam Cancer 2021. https://pmc.ncbi.nlm.nih.gov/articles/PMC8484084/
