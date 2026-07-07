# A-Priori Gene Lists & Phenotype Priors

How curated disease-gene knowledge bases and HPO-driven phenotype matching are used as **ranking priors and reporting tiers** — never as hard include/exclude filters — when screening GMKF Kids First per-trio VCFs.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- Gene lists and phenotype scores are **priors/tiers, never hard gates**. Rarity, impact, and trio-QC filtering run **before and independently of** list membership, so a novel gene is never dropped for absence from curation (the **never-drop rule**).
- **Tier 1** (known disease gene: OMIM `(3)` / PanelApp Green / ClinGen Definitive–Strong / ACMG SF / CGC germline) gets **lenient** variant thresholds; **Tier 2** (constraint/expression/PanelApp Amber) is flagged for review; **Tier 3** (novel) is **retained at a lower prior**.
- Phenotype ranker default: **Exomiser 15.1.0** (June 2026; requires **Java 21**), pedigree-aware (MOI AD/AR/XR), paired with its **matching data release** for GRCh38.
- Exomiser's combined score is a **ranking prior, not a hard reporting gate** — do not gate reporting at a fixed score cutoff (that would reintroduce list-blinding).
- Tier-2 constraint up-weighting uses the pipeline constraint defaults: gnomAD **v2.1.1 LOEUF** primary (pLI ≥ 0.9 / LOEUF_v2 < 0.35, or v4 < 0.6 flagged experimental); ClinGen **HI = 3** / pHaplo ≥ 0.86 for haploinsufficiency. Constraint is never a standalone exclusion.
- Secondary-findings overlay: **ACMG SF v3.3** (84 genes) — report P/LP regardless of proband phenotype.
- **Version-pin every list** (PanelApp panel `version`, COSMIC vNN, ClinGen file date, OMIM download date, ACMG SF vN.N, Exomiser data release) in the per-run manifest; vendor/freeze frozen copies rather than fetching live at analysis time.
- **HPO per proband is a hard dependency** for the phenotype layer; consortium probands with sparse/absent phenotype must be handled gracefully (fall back to genotype-only ranking).

## Why lists are priors, not filters

Curated gene lists let you **boost sensitivity** — recover and rank known disease genes even under lenient variant filters — without **blinding discovery**. If list membership were a hard include filter, every genuinely novel diagnosis (a de novo, biallelic, or high-impact hit in a gene not yet curated for the phenotype) would be silently discarded. The pipeline therefore treats every list as a *weight* on an already-surviving candidate, and applies rarity/impact/QC gating first (see [allele_frequency.md](allele_frequency.md), [functional_annotation.md](functional_annotation.md), [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)).

## Disease-gene knowledge bases

| Resource | What it provides | Access / licensing | Reproducibility handle |
| --- | --- | --- | --- |
| **OMIM** (`genemap2.txt`, `morbidmap.txt`, `mim2gene.txt`) | Canonical gene ↔ Mendelian-phenotype map | Free registered account, per-account token URL; REST API rate-limited **250 req/day**; not redistributable | Pin download date |
| **Genomics England PanelApp** + **PanelApp Australia** (separate instance) | Traffic-light gene panels (Green/Amber/Red), MOI | Open, API + TSV | Pin panel `version` string |
| **ClinGen gene–disease validity** + **dosage sensitivity** | Definitive→Refuted validity; HI/TS scores 0–3 | Open, CSV/FTP | Pin file date |
| **COSMIC Cancer Gene Census (CGC)** | Tier 1/2 cancer genes, per-gene germline flag | Free COSMIC login; not redistributable | Pin version (vNN) |
| **ACMG SF** | Actionable secondary-findings gene list | Open (journal + ClinGen zygosity page) | Cite version (vN.N) |

### OMIM

`genemap2`/`morbidmap`/`mim2gene` are the reference gene↔phenotype map. Parse the phenotype **mapping key `(3)` = molecular basis known** to build a defensible "known Mendelian gene" set; the bracket/brace annotations flag non-disease and susceptibility entries — filter on these. OMIM is license-gated (no public mirror): store credentials as CI secrets, download-then-freeze, and pin the release **date**, not a URL.

### PanelApp (Genomics England + Australia)

Traffic-light confidence: **Green** = diagnostic-grade (≥3 unrelated families, or 2–3 with strong functional data), **Amber** = moderate, **Red** = insufficient. Fully API-driven and versioned; per-gene JSON carries `confidence_level` (3/2/1 = green/amber/red), MOI, and the panel `version`. Pediatric-cancer-relevant panels: **Childhood solid tumours (panel 243)** and **Adult solid tumours cancer susceptibility (245)**; PanelApp Australia adds KidGen and UMCCR cancer panels. Pin the panel `version` per run — see [pediatric_cancer.md](pediatric_cancer.md) for the cancer-list union.

```bash
# Fetch a panel's genes + confidence + MOI + version (open API); freeze the JSON.
curl -s "https://panelapp.genomicsengland.co.uk/api/v1/panels/${PANEL_ID}/" \
  > "panelapp_${PANEL_ID}_$(date +%F).json"
# Record the panel version for the run manifest:
jq -r '.version' "panelapp_${PANEL_ID}_$(date +%F).json"
```

### ClinGen validity and dosage

Use **gene–disease validity** (Definitive → Strong → Moderate → Limited/Disputed/Refuted) to weight priors — favour Definitive/Strong for Tier 1. **Dosage sensitivity** gives haploinsufficiency (HI) and triplosensitivity (TS) scores 0–3; files are split by mechanism and by GRCh37/**GRCh38** coordinates. **HI = 3** up-weights LoF/CNV candidates. Note the pipeline scope limitation: dosage scores currently weight **SNV/indel** calls only — CNV/SV detection is a known blind spot (see [Known scope limitations](#known-scope-limitations)). ClinGen files are open; pin the file date.

### COSMIC Cancer Gene Census

Current release **v104 (May 2026)**, which added 2 Tier 1 and 3 Tier 2 genes over v103. (Specific v104 added-gene names circulating in secondary sources are **unverified against the live release notes** and are not asserted here.) Two evidence tiers: **Tier 1** = documented reproducible cancer role; **Tier 2** = emerging. For germline pediatric cancer, subset to the **germline-flagged** genes, not the full somatic-driver census. CGC needs a free COSMIC login; pin the version.

## Actionable secondary-findings overlay (ACMG SF)

**ACMG SF v3.3** (Genetics in Medicine, June/July 2025) is the current list: **84 genes**, adding **ABCD1, CYP27A1, PLN** to v3.2 (Feb 2023, which had added CALM1/2/3). It is now paired with a ClinGen page specifying reportable zygosity (het/hemi/hom) per gene. **Report P/LP secondary-findings hits regardless of the proband's phenotype**, and cite the exact version. The cancer subset of ACMG SF feeds the pediatric-cancer gene union ([pediatric_cancer.md](pediatric_cancer.md)); classification of these hits follows [clinical_classification.md](clinical_classification.md).

## Phenotype-driven ranking (HPO)

Proband phenotypes are standardized to **HPO** terms, then genes are ranked by phenotypic match. **HPO per proband is a hard dependency** for this entire layer — GMKF Kids First probands with sparse or absent phenotype annotation are common, so the pipeline must fall back to genotype-only ranking rather than fail or silently drop candidates.

| Tool | Inputs | Role | Notes |
| --- | --- | --- | --- |
| **Exomiser** | Genotype (VCF) + HPO + pedigree | Primary phenotype+genotype ranker | v15.1.0, Java 21; hiPHIVE/PHIVE; MOI-aware (AD/AR/XR); containerized/config-driven |
| **GADO** | HPO only (gene-network) | Orthogonal phenotype-only prior | No variant input; tie-breaker |
| **Phen2Gene** | HPO only (HPO2Gene KB) | Orthogonal phenotype-only prior | Fast; tie-breaker |

**Exomiser** is the best fit for an automated trio pipeline: it consumes VCF + HPO + pedigree and does MOI-aware variant filtering plus hiPHIVE phenotype prioritization. Pin the **binary version and the matching data release together**.

> **Version/data pairing (critique correction):** the latest Exomiser is **15.1.0** (released June 2026; 15.0.0 was Feb 2026) and **15.x requires Java 21**. Do **not** pair a 15.x binary with data release **2406** — that bundle is documented as compatible with **Exomiser 14.0.0**. Pin a data release confirmed compatible with the 15.x binary you ship (consult the Exomiser data/version compatibility matrix), and record the exact `YYMM` in the run manifest.

**GADO** and **Phen2Gene** are HPO-only, variant-free rankers useful when variant-level evidence is weak; use them as orthogonal priors and tie-breakers, never as filters.

**Never-drop caveat:** Exomiser's combined score is a **ranking prior, not a hard reporting gate**. Do not gate reporting at a fixed combined-score cutoff — doing so contradicts the never-drop rule and would re-blind discovery for probands with weak or missing phenotype data. (LIRICAL's calibrated posteriors may be layered in the same spirit; see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md) for MOI and de novo handling.)

## Tiered integration strategy

Lists and phenotype scores map candidates to reporting tiers. Rarity/impact/QC gating (gnomAD v4.1 grpmax `faf95`, VEP consequence + calibrated in-silico, GATK genotype-refinement PP/GQ) is applied **independently and before** the priors below, so no tier assignment can exclude a variant that already passed.

| Tier | Definition | Threshold treatment | Outcome |
| --- | --- | --- | --- |
| **Tier 1 — Known disease gene** | Variant in OMIM `(3)` / PanelApp **Green** / ClinGen **Definitive–Strong** / ACMG SF / CGC germline gene | **Lenient** variant thresholds (recover more true positives); cross-check MOI vs observed trio segregation | Prioritized |
| **Tier 2 — Strong candidate** | Not established for this phenotype but constraint-/expression-supported: LOEUF/pLI, missense depletion (MPC/missense Z), disease-relevant tissue expression, PanelApp **Amber**, ClinGen Moderate/Limited, dosage HI/TS | Standard thresholds; keep Exomiser/phenotype-match score attached | Flag for research review |
| **Tier 3 — Novel** | Passes all variant QC/rarity/impact filters but appears in no list | Standard thresholds; lower prior | **Retained, never dropped** — preserves genome-wide discovery |

Tier-2 constraint up-weighting uses the pipeline's constraint defaults (see [gene_constraint.md](gene_constraint.md)): gnomAD **v2.1.1 LOEUF** as the established metric, **pLI ≥ 0.9 / LOEUF_v2 < 0.35** (or **v4 LOEUF < 0.6**, flagged **experimental**), **s_het (Zeng 2024) ≥ 0.1** for short genes where LOEUF is underpowered, and **ClinGen HI = 3 / pHaplo ≥ 0.86** for haploinsufficiency. Constraint is a **ranking weight only** and is **never** used to down-weight recessive (biallelic) candidates.

## Reproducible sourcing

- Prefer versioned APIs/FTP over screen-scraping; **record every release/version string** in run metadata: PanelApp panel `version`, COSMIC vNN, Exomiser data `YYMM`, ClinGen file date, OMIM download date, ACMG SF vN.N.
- **Vendor/freeze** the list files inside the container image or a pinned data bundle; do not rely on live fetch at analysis time.
- OMIM and COSMIC are license-gated (no public redistribution) — store credentials as CI secrets and download-then-freeze. PanelApp, ClinGen, HPO, and ACMG SF are open. See [tooling_and_reproducibility.md](tooling_and_reproducibility.md).

## Recommended defaults (this pipeline)

| Item | Default | Notes |
| --- | --- | --- |
| Phenotype ranker | **Exomiser 15.1.0** (Java 21), pedigree-aware, MOI = AD/AR/XR, hiPHIVE | Pair with **matching** GRCh38 data release; **not** 2406 |
| Phenotype cross-checks | **GADO** + **Phen2Gene** (HPO-only) | Orthogonal priors / tie-breakers |
| Exomiser score | **Ranking prior, NOT a hard reporting gate** | Never-drop rule |
| Tier 1 priors | OMIM `(3)` ∪ PanelApp **Green** ∪ ClinGen **Definitive–Strong** ∪ ACMG SF ∪ CGC **germline-flagged** | Lenient thresholds; MOI-checked |
| Tier 2 constraint gates | gnomAD **v2.1.1 LOEUF** primary; pLI ≥ 0.9 / LOEUF_v2 < 0.35 (v4 < 0.6 experimental); s_het ≥ 0.1 (short genes); ClinGen **HI = 3** / pHaplo ≥ 0.86 | Ranking weight only; never a standalone exclusion; never down-weights recessive |
| Secondary-findings overlay | **ACMG SF v3.3 (84 genes)** | Report P/LP regardless of phenotype |
| Cancer panels | PanelApp GE **Green** on panels **243** + **245** (+ PanelApp-AUS cancer/KidGen) | Pin panel `version`; see pediatric_cancer.md |
| COSMIC | **v104** (May 2026), germline-flagged subset | Added-gene names not asserted |
| Never-drop rule | Tier-3 (no-list) variants passing rarity + impact + trio-QC | Retained at lower prior |
| Rarity gating (applied first) | grpmax `faf95` < **1e-4** dominant/de novo; < **1e-2** recessive (**1e-3** high-conf tier); ≥ **0.05** hard benign | See allele_frequency.md; gene-specific ClinGen VCEP overrides generic cutoffs |
| Reproducibility | Pin OMIM date, COSMIC vNN, PanelApp versions, ClinGen dates, ACMG SF vN.N, Exomiser data `YYMM` | Per-run manifest; freeze copies |

## Known scope limitations

- **HPO dependency.** The entire phenotype layer requires per-proband HPO terms. Consortium probands frequently have sparse or absent phenotype annotation; the pipeline must degrade gracefully to genotype-only ranking rather than fail.
- **CNV/SV blind spot.** ClinGen dosage/HI scores currently weight **SNV/indel** candidates only. ~10–15% of pediatric-cancer and rare-disease diagnoses are CNV/SV (single-exon *RB1*/*SMARCB1*/*DICER1*/*NF1* deletions, *PMS2* rearrangements) and are not detected in the SNV/indel-only pipeline. A future GATK-gCNV / Manta / ExomeDepth module would address this.
- **Pseudogene/segmental-duplication genes** (*PMS2*/*PMS2CL*, *CYP21A2*, *SMN1/2*, *NEB*, *GBA*) are low-confidence from short reads regardless of list membership — flag those regions rather than treat list hits there as reliable.
- **Proband mosaicism** (e.g., *NF1*, overgrowth) produces low-VAF calls that fall outside the het AB band and need a dedicated mosaic tier — list membership does not rescue them.
- **Calibration/validation.** Tier assignment inherits the pipeline's need for truth-set benchmarking (GIAB/CMRG), synonymous λ ≈ 1 checks, and a positive-control variant panel to measure sensitivity/precision.

## Sources

- ACMG SF v3.3 (2025), Genetics in Medicine: https://www.gimjournal.org/article/S1098-3600(25)00101-7/fulltext ; announcement (84 genes; ABCD1/CYP27A1/PLN): https://www.eurekalert.org/news-releases/1090415
- ACMG SF v3.2 (2023): https://pmc.ncbi.nlm.nih.gov/articles/PMC10524344/
- PanelApp (GE): https://panelapp.genomicsengland.co.uk/ ; Childhood solid tumours (243): https://panelapp.genomicsengland.co.uk/panels/243/ ; Adult susceptibility (245): https://panelapp.genomicsengland.co.uk/panels/245/ ; API JSON: https://panelapp.genomicsengland.co.uk/api/v1/panels/
- PanelApp Australia: https://panelapp.agha.umccr.org/
- ClinGen dosage sensitivity + FTP: https://dosage.clinicalgenome.org/ ; https://clinicalgenome.org/docs/clingen-dosage-sensitivity-ftp-download-files-announcement/ ; validity/dosage in classification: https://pmc.ncbi.nlm.nih.gov/articles/PMC9035475/
- COSMIC Cancer Gene Census + release notes (v104, May 2026): https://cancer.sanger.ac.uk/cosmic/release_notes ; tier system: https://cosmic-blog.sanger.ac.uk/cancer-gene-census-hallmarks-and-new-tier-system/
- OMIM downloads/API (250 req/day; genemap2/morbidmap): https://www.omim.org/downloads/ ; https://www.omim.org/api
- Exomiser releases (15.1.0, Java 21) + data compatibility: https://github.com/exomiser/Exomiser/releases ; https://github.com/exomiser/Exomiser/discussions/562 ; install docs: https://exomiser.readthedocs.io/en/latest/installation.html
- Phen2Gene (NAR GAB 2020): https://academic.oup.com/nargab/article/2/2/lqaa032/5843800 ; prioritizer evaluation incl. GADO/Exomiser: https://academic.oup.com/bib/article/23/2/bbac019/6521702
