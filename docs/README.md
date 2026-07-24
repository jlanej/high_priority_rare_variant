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
| **[limitations.md](limitations.md)** | **What the first pass cannot see, why, and the cost to fix each. Read before interpreting a negative result.** |
| [pipeline_design.md](pipeline_design.md) | Vetted end-to-end flow; critique of the original 5-step proposal; data artifacts; scope limits |
| [ROADMAP.md](ROADMAP.md) | Prioritized high-priority gaps vs state of the art (from a SOTA review); dependency-ordered |
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

> ### ⚠ VEP-only contract — read this before the tables
>
> Every annotation the pipeline reads comes from **one** source: a VEP 115 GRCh38 cache plus the
> CADD plugin. No gnomAD / ClinVar / dbNSFP / SpliceAI / LOFTEE file is downloaded or
> bind-mounted. Several rows below therefore describe **targets and reference science, not
> what runs** — each is marked. The **IMPLEMENTED** column is what the code does.
>
> Not available under this contract, and what each costs:
> | Absent | Consequence |
> |--------|-------------|
> | `faf95` | The cache has no AC/AN, so the 95% CI correction **cannot be computed at any price**. Rarity uses a grpmax **point-estimate proxy**, which runs slightly stringent on low-count alleles. |
> | `nhomalt` | No gnomAD-homozygote sanity check on de novo calls. |
> | SpliceAI | **Deep-intronic and exonic-synonymous splice variants are invisible to the screen.** CADD is a lossy proxy (v1.6+ ingests SpliceAI as an input feature). The single largest loss. |
> | LOFTEE | No HC/LC pLoF confidence. Near-inert for *selection* (HIGH impact already keeps every pLoF); matters for the planned tiering step. |
> | ClinVar stars | No `CLNREVSTAT` ⇒ the ≥2★ gate is unimplementable; unstarred P/LP is honored. ClinVar is also as stale as the cache (VEP 115 ⇒ ClinVar 2025-02). |
> | REVEL / AlphaMissense / MPC | **No loss to selection** — see the note under the functional table. |

### Frequency oracle — IMPLEMENTED
- **gnomAD v4.1** (GRCh38; 730,947 exomes + 76,215 genomes), read from the **VEP cache** via
  `--af_gnomade` / `--af_gnomadg`.
- Filter field = **grpmax proxy** = max AF over the grpmax-**eligible** ancestry groups only:
  `AFR, AMR, EAS, NFE, SAS` (`src/hprv/annotations.py:GRPMAX_POPS`).
- Two things this is deliberately **not**:
  - **Not `faf95`.** It is a point estimate. faf95 needs AC/AN; the cache carries neither.
  - **Not VEP's `MAX_AF`.** MAX_AF maximises over the bottlenecked founder groups gnomAD's own
    grpmax *excludes* (`ami` AN≈900, `asj`, `fin`, `mid`) **and** the tiny 1000 Genomes
    populations. One allele in `ami` reads as AF≈1.1e-3 — ten-fold over the dominant gate — so
    using MAX_AF would **silently drop real ultra-rare candidates**. Excluding those groups is
    the entire reason the proxy is defensible. (Enforced by a test; see
    `tests/test_pure.py:test_frequency_excludes_bottlenecked_pops`.)
- **Never** use internal cohort AC/AN as population frequency; internal recurrence is valid only
  as an artifact/blocklist signal.
- Caveat: cache frequencies exist only for alleles **accessioned into dbSNP**, so an
  un-accessioned gnomAD variant returns no AF and reads as "absent ⇒ rarest". That biases toward
  retention (more review), not toward missed calls.

### Rarity gates (grpmax proxy AF) — IMPLEMENTED. A screening gate, distinct from ACMG **PM2**
| Mode | Keep candidate if AF < | Notes |
|------|----------------------|-------|
| Dominant / de novo | **1e-4** | the old `nhomalt ≤ 1` de novo condition is **removed** (no nhomalt) |
| Recessive / comp-het | **1e-2** per allele (permissive); **1e-3** high-confidence tier | applied per variant, not per gene |
| Benign, all modes | drop if AF ≥ **0.05** (ClinGen BA1) | never rescue |

PM2 is applied at **Supporting** strength only and is *evidence*, not the rarity gate itself.

### Functional / in-silico
**IMPLEMENTED — the entire ladder is two rungs**, tried in order (`src/hprv/selection.py`):

| # | Signal | Cutoff | Reaches |
|---|--------|--------|---------|
| 1 | VEP **IMPACT** | keep if `HIGH` or `MODERATE` | all pLoF + all missense + inframe indels |
| 2 | **CADD PHRED** | keep if ≥ **25.3** | everything below MODERATE — intronic / synonymous / UTR / regulatory |

CADD is therefore the **only** functional predictor, and the **only** keep-path for any
non-coding variant. Two honest caveats on that 25.3:
- **Provenance error in the name.** 25.3 is Pejaver-2022's PP3-*supporting* cutoff, calibrated on
  **missense only**. Missense never reaches rung 2 (it is MODERATE, kept at rung 1), so in
  practice 25.3 is applied *exclusively* to the non-coding variants it was **not** calibrated for.
  Read it as a discovery rank (≈ top 0.3% genome-wide), **not** as ACMG PP3 evidence.
- There is no ClinGen-endorsed non-coding CADD threshold to replace it with.

**Why REVEL / AlphaMissense / MPC are not listed — they never did anything.** They are
missense-only scores; every missense is `IMPACT=MODERATE`; rung 1 keeps it and returns *before*
any predictor is consulted. So those branches were **unreachable even when the code contained
them** — removing dbNSFP cost the screen exactly zero discrimination. This is asserted in CI
(`assert_integration.py`: no site may be kept via `revel`/`alphamissense`/`mpc`).

**On "never stack correlated tools":** the ladder is an OR, but with one live functional rung
there is nothing to stack. If you ever narrow `keep_impacts` to `[HIGH]`, missense would fall
through to CADD alone — coherent, but note ClinGen's one-tool rule governs **PP3/BP4 evidence
assignment**, and this screen assigns no ACMG weight.

*TARGET (not implemented — needs resources this contract does not have):* LOFTEE HC-no-flags +
Abou-Tayoun PVS1 grading; REVEL PP3 0.644/0.773/0.932 + BP4 ≤0.290/≤0.183; AlphaMissense
≥0.564; SpliceAI Δ≥0.2 (Walker-2023, calibrated on **raw** scores); MPC ≥2. These specify the
planned ACMG tiering step. If tiering is built, ClinGen SVI says commit to **one** predictor
(REVEL is the ClinGen-calibrated choice), chosen before seeing results.

### Clinical evidence
- **IMPLEMENTED**: ClinVar `CLIN_SIG` from the VEP cache. P/LP (excluding `conflicting`)
  overrides a failed rarity/function screen. **No star gate** — the cache has no `CLNREVSTAT`,
  so a 1★ single-submitter assertion is indistinguishable from an expert-panel one and is
  honored. Over-retention (more to review), never over-dropping. Release is pinned by the cache
  (VEP 115 ⇒ ClinVar 2025-02), not independently.
- *TARGET*: auto-promote P/LP at **≥2★** only; 1★ → prioritize + human review; Conflicting/VUS →
  flag; exclude 0★. Classifier backbone **AutoGVP**; combining via **Tavtigian/ClinGen points**
  (P ≥ 10, LP 6–9, VUS 0–5), **PM2 at Supporting**. All require a ClinVar VCF.

### Gene constraint — a **ranking weight, never a standalone exclusion filter**
- gnomAD **v2.1.1** LOEUF (established) primary; pLI ≥ 0.9 / LOEUF_v2 < 0.35 (v4 < 0.6, flagged
  experimental). Prefer **s_het (Zeng 2024) ≥ 0.1** for short genes. pHaplo ≥ 0.86 / ClinGen HI = 3.
- **Do not** down-weight recessive candidates by pLoF constraint.

### Inheritance models (Step 5) & genotype QC (GATK-refined trios)
- Trust **refined `PP`-derived GQ**. GQ ≥ 20; DP ≥ 10; het AB 0.25–0.75; hom-alt AB ≥ 0.90;
  hom-ref AB ≤ 0.10 (AB from AD); FILTER = PASS only.
- **Dominant** (inherited): rare (grpmax **proxy** AF `< 1e-4` — the field defined under
  [Frequency oracle](#frequency-oracle--implemented), **not** `faf95`), functional **het**
  transmitted from ≥ 1 parent (origin recorded) — the recurrence signal Step 6 consolidates.
- **Recessive**: homozygous, or **compound het** = two rare hets, same gene, in **trans**
  (parent-of-origin from trio genotypes; read-backed **WhatsHap** phasing is a **TARGET**, not
  wired — a pair whose second hit is a de novo is emitted but flagged `unphased_denovo_partner`,
  since trio genotypes cannot phase a de novo against an inherited variant).
- **X-linked recessive**: affected male = hemizygous + carrier mother (the **father's chrX is not
  required** — he transmits Y to a son; flagged `father_carries_x_allele` if he carries); affected
  female = `1/1` + carrier mother + hemizygous-affected father. Sex-aware ploidy, drop male non-PAR
  het calls; kid sex inferred from chrX heterozygosity when the PED is unknown. **X-dominant is not
  a separate mode** (a female's X het is emitted as `dominant`, a male hemizygote as
  `x_linked_recessive`); **chrY yields no inherited call** — the hemizygous models are keyed on the
  mother, which is meaningless on Y.
- **De novo** (SECONDARY / cross-reference only — filtering & review handled by separate
  machinery): `hiConfDeNovo` (child-membership checked) → re-verify DP/AB + parental cleanliness.
- **Sample QC (Step 0)**: trio kid/dad/mom roles come from upstream **peddy**; Step 0 guards the
  less-curated trios via genome-wide **Mendelian-error < 2%** (`qc.mie_max`), chrX-inferred sex vs.
  PED (het-ratio **< 0.10 → male**, `qc.x_het_male_max`; needs **≥ 20** informative chrX calls,
  `qc.sex_min_sites`, else sex is left unknown — a dedicated indexed chrX pass so the autosomal MIE
  cap can't starve it), and a **contamination** gate — verifyBamID **FREEMIX > 0.05**
  (`qc.freemix_threshold`) if a `*.selfSM` directory is supplied (`resources.selfsm_dir`), else a
  VCF-only **CHARR** proxy (reference-read fraction at high-quality hom-alt SNV sites) **> 0.02**
  (`qc.charr_threshold`). (Richer somalier ancestry/relatedness is a roadmap follow-on; CHARR: Lu
  et al., AJHG 2023.)
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
- **Calibrated recurrence null (primary signal):** observed carriers (≥ **min_carriers** — a lone
  carrier is not recurrence) are tested against `Binomial(N_trios, p)` with a **model-appropriate**
  per-individual carriage probability over the gene's qualifying variants (absent → floor
  `burden.absent_af_floor`, default 1e-6). `q_v` is the per-variant **grpmax proxy** AF (Step 5's
  `grpmax_af` column), **not** `faf95`: **dominant het** `p = 1 − Π_v (1 − q_v)²` (the primary,
  FDR-corrected headline); **biallelic** `p = (Σ_v q_v)²`; **X-linked male** `p = 1 − Π_v
  (1 − q_v)` (hemizygous). A recessive/hemizygous carrier is **not** a ≥1-of-two-alleles event,
  so it is never charged the dominant probability. → per-gene `p_recurrence` (+ `p_recurrence_biallelic`,
  `p_recurrence_xlinked`), BH `q`, and an exome-wide flag (`p < burden.exome_wide_p`, default 2.5e-6).
  `N_trios` is the **screened** population (passed explicitly; never inferred from the trios that
  happened to have a call). This makes 2 carriers of a *private* variant genome-wide significant
  while 2 carriers of a common-ish one are not. The proxy reads slightly **high** vs `faf95`, which
  inflates `p` — so these p-values are **conservative**, the safe direction for a discovery claim
  (the same substitution is *un*favourable at the rarity gate; see
  [gene_burden.md](gene_burden.md#multiple-testing-correction)). *(Case-only approximation from in-cohort
  variants; a gnomAD-derived per-gene cumulative-AF test — TRAPD/CoCoRV — is the upgrade.)*
- A gene is **recurrent** at ≥ **min_carriers** (default **2**) distinct individuals; rank
  recurrent-first, **then by `p_recurrence`**, then **weighted by constraint** (LOEUF / pLI /
  s_het / pHaplo — a recurrent het in a haploinsufficient gene is the most compelling).
- OPTIONAL secondary: de novo Poisson enrichment vs the Samocha model (exome-wide **P < 2.5e-6**,
  BH **q < 0.05**) when a mutation-rate table is supplied.

### A-priori gene lists & phenotype — priors/tiers, never hard include/exclude (**never-drop rule**)
- Tier 1 known gene → lenient thresholds; Tier 2 strong candidate (constraint/expression); Tier 3
  novel → retained at lower prior. Rarity/impact/QC gating applied *before and independently of*
  list priors. Phenotype: **Exomiser** + LIRICAL as ranking priors (not hard gates); HPO per
  proband is a dependency, degrade gracefully when sparse.

### Review export (Step 8) & non-human-fraction (Step 8b)
- **Step 8** (igv.js review export): `outputs.igv.padding` **1000** bp mini-CRAM flank;
  `extract_jobs` = `runtime.threads` (samtools `-@` per slice; slices run **serially**).
- **Step 8b** (non-human-fraction — *optional review aid, outside the VEP-only annotation contract;
  never a selection filter*): `outputs.igv.nonhuman_screen.enabled` **true** — but activates only
  when `resources.kraken2_db` is set (else warns + skips, NHF columns blank). `members` **carriers**
  (child always + a parent where it carries the ALT; alternatives `child_only`, `all`); kraken2
  `confidence` **0.05** (off kraken2's 0.0 default, so a lone k-mer can't call a read non-human);
  `min_reads` **5** (denominator floor for the derived `nhf_flag`; the raw `*_nhf_reads` count is
  always emitted); `memory_mapping` **true** (warm page cache across the serial invocations).
  `nhf_flag` fires at NHF ≥ **0.5** over ≥ `min_reads` reads in any screened member.

### Reproducibility / tooling
- One image `FROM ensemblorg/ensembl-vep:release_115.0` (VEP comes from the group's validated base
  image, **not** conda) with a micromamba layer (`env/environment.yml`) for the CLI/Python tools:
  bcftools/htslib/samtools **1.23** (not 1.22 — slivar 0.3.4 needs htslib ≥ 1.23.1), bedtools
  **2.31.1**, vcfanno **0.3.3**, slivar **0.3.4**, whatshap **2.3**, somalier **0.2.19**, Python
  cyvcf2/pysam/pandas/numpy/scipy. **conda-lock** (byte-identical, hash-pinned rebuilds) is a
  **TARGET** — today the env pins versions, not hashes. VEP cache **external, release-matched**. CI
  → **GHCR** per commit (buildx, provenance + SBOM, tag by SHA, amd64). Apptainer: real-disk
  `TMPDIR`, **no `--containall`**, `--cleanenv`. **No hard paths / no PHI / dbGaP-safe.**

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
