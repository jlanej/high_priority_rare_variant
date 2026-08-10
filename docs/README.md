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
| [prioritization.md](prioritization.md) | **Step 9** — gene excess over mutational target (NB2, trimmed fit, mid-p calibration), the six-signal artifact panel, the graded gene down-weight, per-variant tiering, and the additive `priority_points` composite |
| [tooling_and_reproducibility.md](tooling_and_reproducibility.md) | Container/conda-lock, GHCR CI, Apptainer, PHI-safe repo |
| [resources.md](resources.md) | **How to acquire and prepare the annotation data** (VEP cache, CADD, SpliceAI, kraken2); what `prepare_resources.sh` fetches |
| [cram_access_phase.md](cram_access_phase.md) | *Idea doc, not scheduled* — bundling the roadmap items that all need re-access to source CRAMs |
| [pipeline_review_2026-07.md](pipeline_review_2026-07.md) | *Historical findings log (point-in-time), not canonical methods* — a publication-readiness review with a per-finding resolution table |

---

## Canonical defaults

This is the **single source of truth** for every threshold in the pipeline. All values are
**configurable defaults** (see [`config/config.example.yaml`](../config/config.example.yaml)),
not immutable law. A gene-specific ClinGen VCEP value **overrides** any generic cutoff here.

> ### ⚠ VEP-only contract — read this before the tables
>
> Every annotation the pipeline reads comes from **one** tool: VEP 115 GRCh38 — its cache plus its
> plugins (CADD, and **SpliceAI** — required by default, `resources.vep.spliceai_required: true`).
> No gnomAD / ClinVar / dbNSFP / LOFTEE file is
> bcftools-transferred in. Several rows below therefore describe **targets and reference science, not
> what runs** — each is marked. The **IMPLEMENTED** column is what the code does.
>
> Not available under this contract, and what each costs:
> | Absent | Consequence |
> |--------|-------------|
> | `faf95` | The cache has no AC/AN, so the 95% CI correction **cannot be computed at any price**. Rarity uses a grpmax **point-estimate proxy**, which runs slightly stringent on low-count alleles. |
> | `nhomalt` | No gnomAD-homozygote sanity check on de novo calls. |
> | ~~ClinVar stars~~ | **RESTORED.** The ClinVar sites VCF is now transferred in Step 2 (`resources.clinvar.vcf`), so `CLNREVSTAT` -> `clinvar_stars` (0-4) is available. It is a Step-9 RANKING input, never a keep/drop gate. Absent transfer = blank, which is *not* 0 stars. |
> | ~~REVEL / AlphaMissense~~ | **RESTORED** as VEP plugins (dedicated files, not dbNSFP). Still **no effect on selection** — see the note under the functional table; they make Step 9's missense tier calibrated instead of an off-label CADD rank. |
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

**Where each gate fires** (they are not all Step 3): Step 3 applies only the PERMISSIVE union
cutoff `recessive_max` **1e-2** plus BA1 **0.05**, so a dominant-range variant is NOT dropped
there. `dominant_max` **1e-4** and `recessive_strict` **1e-3** are applied at **Step 5**
(`05_inheritance_screen.py`), the latter tagging a call `high_conf_rarity` rather than dropping
it. That ordering is deliberate: selection stays inheritance-agnostic, and the mode-specific
cutoff is only meaningful once a mode has been assigned.

PM2 is applied at **Supporting** strength only and is *evidence*, not the rarity gate itself.

### Functional / in-silico
**IMPLEMENTED — the ladder is three rungs (an OR: any one keeps)**, tried in order (`src/hprv/selection.py`):

| # | Signal | Cutoff | Reaches |
|---|--------|--------|---------|
| 1 | VEP **IMPACT** | keep if `HIGH` or `MODERATE` | all pLoF + all missense + inframe indels |
| 2 | **SpliceAI** max Δ | keep if ≥ **0.2** (`spliceai_ds_min`; ClinGen SVI PP3-supporting) | deep-intronic cryptic splice sites + exonic-synonymous splice disruption — the class VEP's positional terms and CADD both under-call. Required by default (`resources.vep.spliceai_required` **true** — missing score files HALT at Step-2 preflight; set it false to degrade with a warning; not enforced when `resources.vep.annotated_vcf` is set); keep-only |
| 3 | **CADD PHRED** | keep if ≥ **25.3** | everything else below MODERATE — intronic / synonymous / UTR / regulatory |

SpliceAI and CADD are the two keep-paths below MODERATE impact (SpliceAI checked first, so a splice
hit is labelled `spliceai`, not the generic `cadd`). Lower `spliceai_ds_min` toward 0.1/0.05 for
higher deep-intronic recall.

**SpliceAI availability defaults** (this table is the single source of truth, so they live here):
`resources.vep.spliceai_required` **true** — missing raw score files HALT at the Step-2 preflight;
not enforced when `resources.vep.annotated_vcf` is set. `resources.vep.spliceai_backfill.enabled`
**false** — the screen runs on the PRECOMPUTED scores alone; variants the precomputed set does not
cover (insertions > 1 nt, deletions > 4 nt) simply carry no splice evidence, which never drops them
but cannot rescue them either. Set it `true` to score that gap live (`indels_only` **true**,
`distance` **500** bp); an absent isolated `spliceai` env then HALTS at preflight, while a transient
scoring failure degrades to precomputed-only with the union left intact.

Two honest caveats on that CADD 25.3:
- **Provenance error in the name.** 25.3 is Pejaver-2022's PP3-*supporting* cutoff, calibrated on
  **missense only**. Missense never reaches rung 3 — the CADD rung (it is MODERATE, kept at rung 1),
  so in practice 25.3 is applied *exclusively* to the non-coding variants it was **not** calibrated for.
  Read it as a discovery rank (≈ top 0.3% genome-wide), **not** as ACMG PP3 evidence.
- There is no ClinGen-endorsed non-coding CADD threshold to replace it with.

**Why REVEL / AlphaMissense / MPC are not listed as SELECTION evidence — they cannot be.** They
are missense-only scores; every missense is `IMPACT=MODERATE`; rung 1 keeps it and returns
*before* any predictor is consulted. So those branches are **unreachable regardless of whether
the resource is configured** — this is a property of the ladder, not of the annotation contract.
Asserted in CI (`assert_integration.py`: no site may be kept via `revel`/`alphamissense`/`mpc`).

REVEL and AlphaMissense **are** wired now (VEP plugins, dedicated files — see
[resources.md](resources.md)), and that changes nothing above: they add and remove no candidates.
Their consumer is **Step 9's missense tier**, which without them can only report an off-label
CADD rank (`missense_evidence_source=cadd_offlabel`). Configuring them makes that tier calibrated.
If you expected a sharper screen, this is not it — see [prioritization.md](prioritization.md).

**On "never stack correlated tools":** the ladder is an OR, but with one live functional rung
there is nothing to stack. If you ever narrow `keep_impacts` to `[HIGH]`, missense would fall
through to CADD alone — coherent, but note ClinGen's one-tool rule governs **PP3/BP4 evidence
assignment**, and this screen assigns no ACMG weight.

*Status:* REVEL PP3 0.644/0.773/0.932 + BP4 ≤0.290 and AlphaMissense ≥0.564 are **implemented as
Step-9 tier cut points** (`prioritization.variant_tier.revel_*` / `alphamissense_*`), consulted in
a fixed precedence — REVEL, then AlphaMissense, then off-label CADD — never as a max over whatever
is available, because best-of-N is an uncalibrated cherry-pick. hprv still assigns **no ACMG
weight**: the same cut points are used to ORDER candidates.
*Still NOT implemented (needs resources this contract does not have):* LOFTEE HC-no-flags +
Abou-Tayoun PVS1 grading; MPC ≥2. These specify the
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
  VCF-only **raw (uncorrected) reference-read fraction** at high-quality hom-alt SNV sites **> 0.02**
  (`qc.charr_threshold`) — a CHARR-*like* proxy, **not** the calibrated Lu-2023 CHARR statistic (no
  `/mean(1−AF)`, no baseline subtraction; Step 0 runs pre-annotation so per-site AF is unavailable).
  It reads only ~1/3 of the true contamination fraction, so it flags **gross (≳5–8%)** contamination,
  not the 1–3% band; near 1–3%, use the FREEMIX path. (Richer somalier ancestry/relatedness is a
  roadmap follow-on; CHARR: Lu et al., AJHG 2023.)
- **Failure mode**: gnomAD priors in CalculateGenotypePosteriors can suppress genuine ultra-rare
  pathogenic calls — cross-check pre-refinement `PL` for top candidates.

### Pediatric cancer — TARGET (overlay not yet wired; `[reserved]` in config, no code reads it)
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

### Prioritization (Step 9) — IMPLEMENTED. A re-rank, **never** a drop
Full derivations, enrichment folds, and the sensitivity trade-off table: **[prioritization.md](prioritization.md)**.

**The excess statistic is a QUALITY question, not a biology question.** It answers *is this gene
producing more candidate rows than its mutational target predicts?* — never *is this gene
disease-associated?* A high excess is evidence of a technical or population-genetic anomaly
(mismapping, paralogue collapse, callability defect, founder allele, ancestry-uneven rarity gate,
hypermutable locus); a low excess is **not** evidence a gene is real. A gene can be both an
established predisposition gene *and* a mismapping hotspot — 13 were on the validation cohort.

**Gene layer** — `E_g = C·μ_g`, `μ_g = mu_mis + mu_syn + mu_lof` (gnomAD v2.1.1 Samocha targets):

| Parameter | Default | Notes |
|---|---|---|
| `prioritization.excess.null_model` | **negative_binomial** | Poisson is **2.41× anti-conservative** at α=1e-3 *with its own arm trimmed* (1.82× untrimmed — the trim loop applies to whichever null is being fit, so never quote the figure without that condition); NB is 0.31× (conservative) — the right direction under never-drop. φ = 18.7 raw / **1.29** after trimming **77 genes (0.392%)**, the CONVERGED count |
| `C` fit universe | **FULL table, zero-count genes included** | The candidate list is zero-truncated; matched-only inflates `C` and deflates every ratio |
| `excess.offset.mu_lof_impute_factor` | **0.0516** | Median `mu_lof/(mu_mis+mu_syn)`; a null `mu_lof` is imputed + labelled, never charged 0 |
| `excess.offset.cds_fallback` | `b0` **−18.4484**, `b1` **1.0570** | OLS, R² 0.828, ±30% band. Those genes never reach T3 |
| `excess.trim_p` / `trim_max_fraction` | **1e-3** / **0.05** (HALT) | Untrimmed α = 0.804 vs trimmed 0.214 — the tail hides itself |
| `excess.q_threshold` | **0.05** | BH-FDR over the full universe |
| `excess.min_n_for_ratio_rule` | **3** | The spec's own guard, applied to T2 and both RATIO limbs of T3 (never the FDR limb — no gene with n<5 reaches q<0.05). **92 genes / 2,369 variants (9.57%), 8 exemptions**; set to 1 for the original 228-gene / 2,565 (10.36%) table with 13 exemptions. 100% retention either way. Without it, 136 of 228 triaged genes had n<3 and carried just 7.6% of triaged volume, and 8 of 13 exempted control genes had q_nb=1 |
| `signals.caf_low.include_null_as_flagged` | **true** | Percentile over genes that HAVE a value (cut 7.96e-6 on the validation cohort); nulls then FLAGGED as low-information loci. The literal "nulls read as 0" reading would let missing data move the threshold applied to the measured data |
| `excess.offset.covariate_adjust` | **false** | ΔAIC −121 but Spearman 0.992, and `oe_syn` is itself a reported signal |

**Artifact panel** (six orthogonal signals → an unweighted integer `corroboration_count` 0–6; fold
enrichment at `excess_ratio ≥ 10`):

| Signal | Default | Fold |
|---|---|---|
| `signals.saturation.per_trio_min` | **0.10** (~2× the 99th pct) | **195×** — strongest, offset-free, but never sufficient alone |
| `signals.segdup.min_identity` / `min_frac` | **0.98** / **0.10** | **10.1×** — build is load-bearing (gnomAD v2.1.1 coords are **GRCh37**) |
| `signals.family.patterns` | 13 curated regexes | **6.5×** — one vote only; OR\*/ZNF\* are n.s. on real data |
| `signals.oe_syn.max_deviation` | **0.30** (symmetric) | **5.2×** — both directions informative, for different reasons |
| `signals.constraint_flag.values` | **mis_too_many, syn_outlier** | **4.0×** — `no_exp_lof` deliberately excluded |
| `signals.caf_low.percentile` | **0.10** | **3.9×** — a LOW-INFORMATION-locus flag, **not** a frequency flag |

**Gene down-weight** — graded, four-tier, and **never a hard drop**:

| Tier | Rule | Penalty |
|---|---|---|
| `T1_watch` | `ratio ≥ 3 & q_nb < 0.25` (no corroboration needed) | **−0.5** |
| `T2_downweight` | `ratio ≥ 5 & n_g ≥ 3 & corrob ≥ 1` (ratio not FDR: no gene with n<5 reaches q<0.05) | **−1.5** |
| `T3_strong_downweight` | `(q<0.05 & corrob≥1)` **or** `(ratio≥10 & n_g≥3 & corrob≥2)` **or** `(ratio≥5 & n_g≥3 & saturation & corrob≥2)` | **−3.0** |
| `established_gene_ceiling` | **T1_watch** — a control-union gene never enters T2/T3. Still needed after the count floor: CTSA, NPRL3, CDH23 require it | auditable, reversible |
| `cds_fallback_ceiling` | **T2_downweight** — a ±30% offset can't support a 5× claim | |
| `min_control_genes` / `max_downweight_fraction` | **1000** / **0.20** | both **HALT**, not degrade |

Measured: **92 genes / 2,369 variants (9.57%) down-weighted at 100% established-gene retention** with
the `n_g ≥ 3` floor (8 named exemptions); **228 genes / 2,565 (10.36%)** without it (13 exemptions).
Retention survives expanding the control union to **2,218 genes** — every added source (ACMG SF v3.3,
PanelApp 243 v5.12, PanelApp 259 v1.30) maxes below the 5× cut and triages zero, so the figure is
stress-tested rather than merely measured. Control genes are *depleted* ~2× in the down-weight region.
The maximum penalty (−3.0) **cannot alone demote** a variant with strong molecular evidence in a
constrained gene (V4 +4, rarity +2, constraint +1 = +7 → +4). Uncorroborated high excess is
**flagged (`unexplained_excess`), never penalised.** The exemption withholds a score *penalty* and
never endorses the calls — `CTSA` (0.73 variants per trio at 148× its target) is not credible as
biology, and those genes carry `review_flag = established_gene_high_excess` and sort to the **top** of
a read-level review list.

**Variant tiers** (a screening/triage tier, **not** an ACMG classification):

| Tier | Rule | Points |
|---|---|---|
| **V5** | NMD-competent pLoF in a LoF-mechanism gene | +8 — **UNREACHABLE**: `variants.tsv` has no EXON/CDS_position, so every pLoF caps at V4 |
| **V4** | `spliceai_ds ≥ 0.5`, or a HIGH-impact pLoF (NMD indeterminate) | +4 |
| **V3** | `spliceai_ds ≥ 0.2`, or missense `cadd ≥ 25.3` | +2 — the CADD route is `cadd_offlabel`, a discovery rank, **not** PP3 |
| **V2** | missense below the CADD cut; in-frame indel | +1 |
| **V1** | non-coding/synonymous kept via the CADD rung | +0.5 |
| **V0** | `spliceai < 0.1` **and** `cadd < 15` **and** LOW/MODIFIER | 0 — **caps the total**; BOTH scores must be PRESENT (absence ≠ benignity) |

**Mechanism gating** (ACMG/ClinGen SVI — the most important structural rule): the constraint term
**and** the gene-list prior are multiplied by **V0 → 0.0, V1/V2 → 0.5, V3–V5 → 1.0**, and zeroed
entirely for `compound_het`/`hom_recessive`/`x_linked_recessive` (pLoF constraint measures selection
against heterozygotes — do not up- *or* down-weight a biallelic candidate by it).

**Composite** `priority_points` — additive, every term a separate reported column:
molecular (above) + rarity (**1e-5 → +2, 1e-4 → +1.5, 1e-3 → +1, 1e-2 → +0.5, BA1 0.05 → −8 and
caps at −4**) + gene constraint (**+1**, gated) + recurrence (**+1 / +2 cap**, on the *carrier
count* — never on the saturating case-only `p_recurrence`; same-variant only **+0.5**) + quality
(GT fail **−2**, NHF flagged **−3**, **NHF not_screened 0**, comp-het partner unknown **−0.5**) +
clinical (ClinVar P/LP **+4**, benign **−4**; `review_status = UNAVAILABLE` — no `CLNREVSTAT`, so no
≥2★ gate) + MOI (discordant **−1**, **unknown exactly 0**) + gene artifact (above).

> **`priority_points` is NOT an ACMG score.** Do not read totals against Tavtigian's P ≥ 10 /
> LP 6–9 / VUS 0–5 bands: the criteria are not ACMG criteria, no phenotype/segregation/functional
> evidence exists, the ClinVar term has no star gate, and the artifact terms have no ACMG analogue.
> Never emit a P/LP/VUS label from it.

**NHF is three states, never two:** `clean` / `flagged` / **`not_screened`**. **Blank ≠ 0.0** —
blank means nobody looked (no ALT carrier, no mini-CRAM, or Step 8b never ran); `0.0` means screened
and clean. `not_screened` scores **0**: an explicit uncertainty flag, neither penalty nor credit.
Treating blank as clean would silently promote exactly the calls nobody examined.

**Two rankings always ship:** `rank_agnostic` (no gene-list prior of any kind) and `rank_prior`
(with the optional Class-B overlay), plus `rank_delta`. `composite.gene_list_prior.enabled` defaults
**false**, and a config overlay path is inert while it is false — with the overlay off the two
rankings are identical (asserted in tests). The prior is a **prior, never a filter**.

**Four hazards in the overlay layer, all guarded in code and each with a test:**

1. **The join is on gene symbol ONLY — MOI-agnostic.** Routing by canonical MOI would silently miss any hypothesis about a *different* genetic model than the gene is curated under. The FA/HR genes (FANCA, FANCD2, SLX4, FANCE, BRCA2) carry germ-cell-tumour evidence about **heterozygous carriers** (PMID 40906985) while their canonical model is biallelic Fanconi anemia; an MOI-routed lookup would never apply them to the het observations the evidence is about.
2. **An observed-mode-vs-canonical-MOI mismatch is FLAGGED, never penalised.** A het in a canonically-recessive gene emits `moi_caveat = moi_mismatch_het_in_recessive_gene` and charges **0** — that is the carrier-risk shape, not an incoherent call. Never-drop, applied to the coherence layer.
3. **Gene-level and set-level priors combine by MAX, never SUM** (`combine_gene_and_set: max`). A curated overlay may carry per-gene rows *and* a pathway-collapsed set entry derived from the **same study** (the GCT resource's `FA_HR_PATHWAY_23` and its per-gene FA rows are both PMID 40906985) — summing double-counts one study.
4. **Non-germline rows contribute ZERO** (`non_germline_classes: [somatic_driver_not_germline]`). Somatic drivers are listed in a curated overlay deliberately, so a reader can see they were excluded; their weight of 0.15 is bookkeeping, **not weak germline support**. Reported with `gene_list_prior_excluded_non_germline`, never silently dropped.

`prior_weight` is **UNCALIBRATED** — an ordering default, never a likelihood ratio or an odds ratio.
And a T3 (GWAS-locus) overlay hit means "this gene sits at a GWAS locus", not "rare coding variants
here matter": those lead variants are predominantly **non-coding**, so a rare-coding screen is not
interrogating that mechanism at all.

**Resource:** `prioritization.resources.mutational_target` — the **unjoined** gnomAD v2.1.1
constraint table (~3 MB; `prepare_resources.sh --only mutational_target`). Not interchangeable with
`resources.constraint.gnomad_v2_constraint`, which is projected down to the LOEUF/pLI priors and has
no `mu_*` columns. Absent → loud WARN, every gene reads T0, the variant layer still runs.

*TARGET:* V5 (needs three more VEP fields), calibrated missense strength (REVEL/dbNSFP), pLoF
confidence (LOFTEE), ClinVar star gate, single-site excess de-escalation, per-ancestry candidate
yield, synonymous-λ calibration.

### A-priori gene lists & phenotype — TARGET (priors/tiers not yet wired; `[reserved]` in config, no code reads gene lists or HPO/Exomiser)
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
