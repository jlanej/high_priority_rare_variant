# Cross-Pedigree Gene Consolidation (Recurrence-Based)

Finds genes where rare, functional **inherited** variants recur across multiple independent individuals in the cohort, so recurrent gene-level signal can nominate candidates beyond single-family analysis.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## Status: what of this doc is live

The pipeline runs a **VEP-only contract** — a VEP 115 GRCh38 cache plus the CADD plugin, nothing
else downloaded or bind-mounted. That reshapes this doc into three layers, marked throughout:

| Layer | What it is | Examples here |
| --- | --- | --- |
| **IMPLEMENTED** | Step 6 (`pipeline/06_gene_burden.py`) does this today | Distinct-individual carrier counts per gene per model; the per-model recurrence nulls; BH-FDR + exome-wide flag; constraint weighting and ranking |
| **TARGET** | Documented, configured, **not** running | De novo Poisson enrichment (needs a mutation-rate table; uncalibrated when it runs), TRAPD corroboration, synonymous-λ calibration |
| **REFERENCE** | Literature the roadmap is built on, not a claim about this code | The statistical-test menu, SAIGE-GENE+/regenie, mask definitions (M1/M2), DeNovoWEST/extTADA |

**The single change that touches every rarity statement below:** there is **no `faf95`** under this
contract. The VEP cache carries gnomAD v4.1 point AFs but no AC/AN, so the 95% CI lower bound is
**unrecoverable, not approximated**. `annotations.frequency()` returns a **grpmax proxy** — the max
AF over the grpmax-eligible groups (AFR/AMR/EAS/NFE/SAS). Where this doc historically said "faf95",
read "grpmax-proxy AF" — and see [§ Rarity field](#rarity-field-for-qualifying-variants) for why
that direction of error is *safe* for the recurrence p-values specifically.

The full ledger of what the first pass cannot see, and what each gap costs to close, is
[limitations.md](limitations.md) — the anchor doc. This doc links there rather than restating it.

## TL;DR

- **Primary signal is recurrence-based gene consolidation** across the cohort, focused on **inherited germline variation**. For each gene, tally the number of **distinct individuals** carrying a qualifying variant under each inheritance model — **dominant** (rare inherited het), **biallelic** (homozygous + compound het), and **X-linked**. A gene is **recurrent** at **≥ `min_carriers` (default 2)** distinct individuals.
- **The key new signal is recurrence of inherited heterozygous variants.** A rare (grpmax-proxy AF < 1e-4), functional, inherited het is only weakly interesting in one family, but becomes compelling when it **recurs across multiple individuals in the same gene**.
- **Rank recurrent genes first, weighted by gene constraint** (LOEUF / pLI / `s_het`). A recurrent het in a **haploinsufficient** gene is the most compelling result; a gene tolerant of damaging variation is de-prioritized even when recurrent.
- **De novo enrichment is an OPTIONAL SECONDARY signal only.** De novo filtering **and review are handled by separate dedicated machinery**; here de novo is carried as a lightweight cross-reference column (GATK `hiConfDeNovo`, child-membership checked). When a mutation-rate table is supplied, an optional Samocha-2014 Poisson enrichment (denovolyzeR-style) is reported — **exome-wide P < 2.5e-6, BH q < 0.05** — but it is not the driver, and it is **uncalibrated** (see [below](#optional-secondary-signal-de-novo-enrichment-vs-a-mutation-model)).
- **The non-joint per-trio design forbids an internal cohort allele frequency.** Absent genotypes are ambiguous (no-call vs hom-ref), so there is no valid internal AC/AN. Internal recurrence here means *distinct-individual carrier counts of qualifying variants*, not a population frequency.
- **TRAPD** vs **ancestry-matched, coverage-intersected gnomAD v4.1 exomes** remains an **optional corroboration** (not yet implemented) — stratification/coverage-sensitive, never a standalone discovery engine.
- **Rarity gate for qualifying variants uses a gnomAD v4.1 grpmax *proxy*** — the max point-estimate AF over the grpmax-eligible groups (AFR/AMR/EAS/NFE/SAS), < 1e-4 for dominant/de novo candidates. **Not `faf95`** (unavailable — no AC/AN in the cache) and **never** VEP's `MAX_AF` or a global AF, both of which fail, in opposite directions ([limitations.md §2a](limitations.md)).
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

| Model | What is counted | Rarity gate (grpmax-proxy AF) |
| --- | --- | --- |
| **Dominant (inherited het)** | Rare, functional **heterozygous** variant transmitted from ≥ 1 parent (parent-of-origin recorded: maternal / paternal / both), **not** part of a compound-het pair | < 1e-4 |
| **Biallelic** | Homozygous **or** compound het (two rare hets, same gene, in *trans*) | < 1e-2 (permissive) / 1e-3 (high-confidence tag) per allele |
| **X-linked recessive** | Affected male (hemizygous) + carrier mother (father's chrX not required); or affected female (hom-alt, carrier mother, hemizygous father); sex-aware ploidy | < 1e-2 (permissive) per allele |
| **De novo (secondary)** | GATK `hiConfDeNovo`, child-membership checked (`annotations.is_hiconf_denovo_for`) | counted in a **separate** column, not part of the recurrence driver |

The dominant-het count is the **key new signal**: individually a rare inherited het is weak evidence, but a gene that accumulates such hets across **multiple distinct individuals** is a strong nomination.

### Recurrence flag and ranking

- A gene is **recurrent** when its distinct-individual carrier count reaches **`min_carriers` (default 2)**.
- Genes are ranked **recurrent-first**, then **distinct-variant** recurrence above same-variant (founder/artifact), then by the strongest per-model recurrence p, then **constraint-weighted** (constrained = LOEUF < 0.35 **or** pLI ≥ 0.9 **or** s_het ≥ 0.1 **or** pHaplo ≥ 0.86), then by carrier / dominant-het counts. A recurrent het in a **haploinsufficient** gene ranks highest.
- De novo carrier counts and the optional de novo enrichment p-value are carried as **secondary columns** used only to break ties after the recurrence and constraint keys.

## Optional secondary signal: de novo enrichment vs a mutation model

> **De novo filtering and review are handled by separate dedicated machinery.** In this pipeline de novo is retained only as a **lightweight cross-reference** (detected via GATK `hiConfDeNovo`, child-membership checked via `annotations.is_hiconf_denovo_for`). The enrichment test below runs **only when a mutation-rate table is supplied** (`--mutrate`; none ships with the pipeline, so it is **off by default in practice**) and is never the primary nomination signal.
>
> **When it does run, it is UNCALIBRATED.** The Poisson expectation is `2 × N_trios × μ` straight
> from the supplied rate table — it is **not** rescaled to an observed synonymous de novo rate.
> `--syn-denovo-count` is reserved and not yet wired into the expectation, and Step 6 prints this
> warning at runtime. An uncalibrated expectation absorbs the cohort's DNM-calling sensitivity
> into the test statistic, so treat these p-values as a **rank**, not a significance claim.
>
> **Two implementation notes (distinct from the uncalibrated-mean caveat above):**
> 1. **BH denominator = the mutation-model universe.** Unlike the recurrence families (whose BH is
>    legitimately conditional on observing ≥ `min_carriers` carriers), `dn_q_enrich` is BH-corrected
>    over **every gene in the `--mutrate` table** — a zero-DNM model gene is a valid null test at
>    `p = poisson.sf(-1, exp) = 1.0`, seeded by padding the p-vector. Correcting over called genes
>    only would make `q` anti-conservative by ~`n_model / n_called`. (`dn_exome_wide_sig` is a fixed
>    `p < 2.5e-6` Bonferroni flag and is unaffected.)
> 2. **A single POOLED protein-altering test.** Step 6 sums the LoF and missense rates/counts into
>    one Poisson test (denovolyzeR's pooled "prot" class) and reports one `dn_p_enrich`; it does
>    **not** emit the per-class (LoF / missense / synonymous) breakdown of `denovolyzeByGene(classes=…)`.
>    A LoF-specific test would give more power/interpretability for haploinsufficient genes — a
>    possible enhancement, but out of scope for this off-by-default secondary arm.

### The Samocha framework

Samocha et al. (*Nat Genet* 2014) derive per-gene, per-consequence-class expected DNM rates from a **trinucleotide-context mutation model** (calibrated on human–chimp intergenic divergence). The observed DNM count for a gene and class is tested against its Poisson expectation, scaled by `2 × N_trios` (two transmissible haplotypes per trio). This is the foundation for the downstream tools below.

### Tooling

| Tool | Method | Fit for this pipeline |
| --- | --- | --- |
| **denovolyzeR** (R) | Per-gene / gene-set Poisson enrichment for LoF, missense, synonymous; CCDS-level exome-wide test | **Model for the optional de novo arm** (Step 6 implements a **single pooled** LoF+missense "prot" Poisson test, not the per-class breakdown — see the notes above). Simplest, well-suited to a modest trio cohort. |
| DeNovoWEST (Kaplanis 2020) | Unified severity-weighted simulation test adding missense-clustering | Higher power on large NDD cohorts; heavier inputs. |
| extTADA / TADA (Nguyen 2017) | Bayesian integration of de novo + case-control counts; per-gene BF/FDR | Useful only if combining both signal arms. |

Newer proteome-wide constraint models refine expectations further but require heavier inputs and are not defaults here.

### De novo variant qualification

DNMs feeding the burden test are drawn from the same genotype-QC gates as the rest of the pipeline (see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)). A qualifying DNM must satisfy:

| Gate | Default | Status |
| --- | --- | --- |
| GATK screen | `hiConfDeNovo` (primary); `loConfDeNovo` = lower-sensitivity tier only | IMPLEMENTED (only required when the tag is present in the callset header) |
| Child GQ (refined PP-derived) | ≥ 20 | IMPLEMENTED |
| DP (de novo) | ≥ 20 | IMPLEMENTED |
| Het allele balance | 0.25–0.75 | IMPLEMENTED |
| Parental cleanliness | each parent alt AD ≤ 1, DP ≥ 10 | IMPLEMENTED |
| gnomAD v4.1 frequency | grpmax-**proxy** AF < 1e-4 | IMPLEMENTED (as a proxy, not `faf95`) |
| gnomAD homozygote sanity check | `nhomalt` | **NOT IMPLEMENTED** — `nhomalt` does not exist in the VEP cache. The `filters.denovo.require_gnomad_absent_or_singleton` key is **retired**, not silently no-op'd: it was only ever implemented as `nhomalt > 1`, a homozygote-count test rather than the allele-count test its name promised ([limitations.md §3](limitations.md)). |
| Region | callable-region intersect | **NOT IMPLEMENTED** — no callable-region BED is defined or applied anywhere in the pipeline. Aspirational; see the [pitfalls](#pitfalls) below, where differential coverage is the confounder it would address. |

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

## Variant qualification

### What actually qualifies (IMPLEMENTED)

Step 6 has **no mask concept**. It consumes whatever Step 5 emitted, and Step 5's variants are whatever Step 3's classifier kept. That classifier (`src/hprv/selection.py`) is a deliberately **three-rung** ladder (IMPACT, SpliceAI, CADD):

1. VEP `IMPACT` ∈ `keep_impacts` (**HIGH, MODERATE**) → keep;
2. **else** `CADD_PHRED ≥ 25.3` → keep;
3. else drop.

Plus two overrides: grpmax-proxy AF ≥ 0.05 (ClinGen **BA1**) drops and is never rescued; a ClinVar **P/LP** assertion keeps, and also rescues a variant that would fail the rarity gate. Consequence classes come from VEP (see [functional_annotation.md](functional_annotation.md)); frequencies from the VEP cache's gnomAD v4.1 point AFs (see [allele_frequency.md](allele_frequency.md)). All qualifying variants additionally require FILTER = PASS and the genotype-QC gates above.

Two consequences a reader must hold onto before interpreting any per-gene tally:

- **The ClinVar override is unstarred.** The cache carries `CLIN_SIG` but no review status, so a 1★ single-submitter P/LP assertion enters the tally indistinguishably from an expert-panel one. It over-retains rather than over-drops.
- **CADD 25.3 is off-label.** It is Pejaver-2022's PP3-supporting cutoff, calibrated on **missense only** — and missense never reaches rung 2, because every missense is MODERATE and rung 1 returns first. So 25.3 is applied *exclusively* to the non-coding variants it was not calibrated for. It is a discovery rank (≈ top 0.3% genome-wide), **not** ACMG PP3 evidence ([limitations.md §4](limitations.md)).

### Masks (REFERENCE / TARGET — not implemented)

The standard burden masks, documented because they are what a mature version of this module would use and what the literature below assumes. **Neither is computable under the VEP-only contract** — M1 needs LOFTEE, M2 needs REVEL/AlphaMissense.

| Mask | Definition | Blocker |
| --- | --- | --- |
| **M1 — pLoF** | LOFTEE **HC, no flags** | No LOFTEE data files ([limitations.md §5](limitations.md)) |
| **M2 — pLoF + damaging missense** | M1 ∪ (REVEL ≥ 0.5 **or** AlphaMissense likely_pathogenic **or** CADD ≥ 20) | No dbNSFP/REVEL/AlphaMissense ([limitations.md §7](limitations.md)) |

AAF tiers **≤ 1e-4** (dominant / de novo) and **≤ 1e-2 / 1e-3** (recessive) *are* live as the per-model rarity gates above. Multiple masks × AAF tiers would improve power but must be paid for in correction where the optional statistical tests are used (see below).

> **Note on M2 and REVEL/AlphaMissense specifically:** their absence costs this *screen* nothing. They are missense-only scores, every missense is MODERATE, and rung 1 keeps MODERATE before any predictor is consulted — so those branches were **unreachable even when the code contained them and dbNSFP was configured**. CI now asserts the corresponding keep-reasons never fire. The genuine loss is tiering and reporting, not selection.

### Rarity field for qualifying variants

The rarity field is a **grpmax proxy**: the max gnomAD v4.1 point-estimate AF over the grpmax-eligible ancestry groups (AFR/AMR/EAS/NFE/SAS), mirroring gnomAD's own grpmax inclusion set.

It is **not `faf95`**. faf95 is the lower bound of the 95% CI, and computing it requires AC/AN, which the VEP cache does not carry — so it is unrecoverable here, not approximated. A point estimate is always ≥ its own CI lower bound, so **rarity gates fire slightly more often than a faf95 gate would**: the screen errs toward dropping on low-count alleles.

The recurrence null (below) inherits the *opposite* and more comfortable side of that same bias — see [Multiple-testing correction](#multiple-testing-correction).

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

The primary recurrence signal is **calibrated**: for each gene, observed distinct-individual carriers are tested against `Binomial(N_trios, p)` with `p = 1 − Π_v (1 − q_v)²` over the gene's qualifying inherited variants (an allele absent from gnomAD is floored at `absent_af_floor`, default 1e-6). This yields a per-gene `p_recurrence`, a **BH `q_recurrence`**, and an exome-wide flag — so 2 carriers of a *private* variant are genome-wide significant while 2 carriers of a common-ish variant are not. It is a **case-only approximation** using in-cohort variants; the gnomAD-derived per-gene cumulative-AF version (TRAPD/CoCoRV) is the planned upgrade. The same thresholds apply to the **optional secondary de novo enrichment**:

> **The per-variant `q_v` is the grpmax-proxy AF read from Step 5's `grpmax_af` column, not `faf95`** — and for *this* test the substitution errs in the safe direction. Because the proxy is a point estimate it reads slightly **high** versus faf95 on low-count alleles; a larger `q_v` **inflates** the null probability of seeing carriers by chance, which **inflates the p-value**. These recurrence p-values are therefore **conservative** — the right direction for a discovery claim, since it costs sensitivity rather than manufacturing significance. Do not read this as "equivalent to faf95": it is a bias, in a known and bounded direction, not an equality. (Note the same substitution is *un*favourable at the rarity gate, where a high point estimate drops real candidates — the direction of harm flips with how the number is used.)

> **The larger bias runs the *other* way — read "significant" as a ranking, not a calibrated claim.** The `q_v` → conservative point above is a *small* effect (roughly one CI-width on low-count alleles). It is dominated by the **case-only approximation**: `p` is built from *only the variants observed in the cohort*, not the gene's full qualifying-variant target in gnomAD. Omitting the rest of that target makes `p` **too small** (anti-conservative) for the gene-level "is this gene recurrently hit?" question — so for essentially any gene with ≥ `min_carriers` carriers of *private* variants, `p_recurrence` clears the exome-wide line (2 carriers of distinct absent→floored variants give `p ≈ 3e-7` even at only N = 200 trios). Concretely, `recurrence_exome_wide_sig` / `q_recurrence` largely **restate the recurrence flag** for private-variant genes: they *order* candidates, they do not *certify* significance. The honest fix is the gnomAD-derived per-gene cumulative-AF null (TRAPD/CoCoRV, the planned upgrade); until then, treat these columns as a **rank**, and lean on constraint + variant class for interpretation.

- **Exome-wide gene-based threshold: P < 2.5e-6** (≈ 0.05 / ~20,000 protein-coding genes). The literature also expresses this as a class-specific Bonferroni; the canonical default here is the single ~2.5e-6 line.
- If **not** using a single omnibus p-value, divide further by the number of masks × AAF tiers × tests. Combining masks/tiers per gene via **ACAT/Cauchy (ACAT-O)** collapses them into one p-value and **avoids that penalty** — the preferred route.
- Report **FDR (BH q < 0.05)** for candidate discovery alongside the exome-wide line.
- Apply the same ~2.5e-6 threshold to the de novo Poisson tests.

## Ranking and interpretation

A recurrent gene is only interesting if it is **intolerant of the class of variation showing the recurrence**. Weight and rank nominated genes by constraint (see [gene_constraint.md](gene_constraint.md)):

- pLoF nominations → gnomAD **v2.1.1** LOEUF (established; v4 flagged experimental), pLI, `s_het` (Zeng 2024) for short genes. In code, "constrained" = **LOEUF < 0.35 or pLI ≥ 0.9 or s_het ≥ 0.1 or pHaplo ≥ 0.86** — any one suffices. All four are read from the optional `--constraint` table; **with no table supplied, every gene is `constrained=0`** and the constraint key drops out of the ranking silently.
- A recurrent **dominant het in a haploinsufficient gene is the most compelling** result; a gene tolerant of damaging variation is **de-prioritized even when recurrent** — it is more likely an artifact than a discovery.
- Cross-reference a-priori gene lists and phenotype priors as **tiers, never hard filters** (never-drop rule; see [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md)), and germline cancer predisposition genes via [pediatric_cancer.md](pediatric_cancer.md).

## Pitfalls

- **Population stratification** — the dominant confounder in external-control burden; matched-ancestry gnomAD subsets + PCA sample assignment would be mandatory for TRAPD. It also reaches the *live* recurrence null: the grpmax proxy takes the max over ancestry groups regardless of the cohort's actual composition, so a variant common in one group and absent elsewhere is charged its highest AF.
- **Differential coverage** — without cohort joint genotyping, low-coverage no-calls masquerade as hom-ref. Restricting to jointly-callable high-depth intervals is the mitigation, and **it is not implemented** — there is no callable-region intersect.
- **Batch/platform effects** — capture kit, caller, and build details differ from gnomAD. The synonymous-λ is the calibration diagnostic that would catch this; it is **not implemented** (see [Scope limitations](#scope-limitations-stated-honestly)), so this confounder is currently undiagnosed rather than controlled.
- **No internal AF** — the recurrence tally is distinct-individual carrier counts against a gnomAD-derived null, never an internal AC/AN. De novo (per-family, batch-insensitive) is a **secondary cross-reference** here, not the primary signal; gnomAD-control burden (TRAPD) would be corroborative and is not implemented.

## Scope limitations (stated honestly)

- **SNV/indel only initially.** This burden module counts SNV/indel qualifying variants. CNV/SV are a real blind spot — 10–15% of pediatric-cancer and rare-disease diagnoses are CNV/SV (single-exon *RB1*/*SMARCB1*/*DICER1*/*NF1* deletions, *PMS2* rearrangements). A gene whose recurrent hits are copy-number events will be missed until a future GATK-gCNV / Manta / ExomeDepth module exists.
- **Pseudogene / segmental-duplication genes** (*PMS2*/*PMS2CL*, *CYP21A2*, *SMN1/2*, *NEB*, *GBA*) yield low-confidence short-read calls; their per-gene counts are unreliable and those regions should be flagged, not trusted, in the burden tally.
- **Proband post-zygotic mosaicism** (*NF1*, overgrowth) produces low-VAF calls that fall outside the het AB 0.25–0.75 band and are dropped before the burden count — a source of false-negative de novos.
- **No calibration diagnostic is implemented.** The synonymous λ ≈ 1 check described above is a **TARGET, not a built-in**: nothing in Step 6 computes a synonymous burden, and `--syn-denovo-count` is reserved and unwired. So the module currently has **no self-check that its qualifying filters are unbiased** — the recurrence null's calibration rests on the gnomAD frequencies being right, untested. Pipeline-wide sensitivity/precision likewise still needs GIAB/CMRG truth sets and a positive-control variant panel: currently **unmeasured**, not measured-and-acceptable.
- **Splice-aware qualification (SpliceAI, now wired).** A gene whose recurrent signal is deep-intronic or exonic-synonymous splice disruption is now nominated on its SpliceAI delta score (rung 2 of the Step-3 classifier), when the SpliceAI plugin is configured — so such a signal reaches this tally instead of relying on CADD to rescue the variant ([limitations.md §1](limitations.md)). Absence of the score files leaves this class to CADD's weak proxy again.
- **Phenotype dependency** — phenotype-driven ranking assumes per-proband HPO terms, which are frequently sparse or absent in consortium data (see [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md)); burden nomination itself is phenotype-agnostic, which is a strength for discovery.

## Example: parameterized qualifying-variant extraction

> **REFERENCE — this is not what the pipeline runs.** Steps 3–6 do this in python over the INFO
> fields `annotations.F` defines; the sketch below is the shell equivalent for ad-hoc work. Two of
> its four stages are **not reproducible under the VEP-only contract**: there is no `LoF` field
> (no LOFTEE) and no callable-region BED. They are retained to show the shape of the intended
> mature pipeline. Substitute your own field names and thresholds from `config/config.example.yaml`.

```bash
# 1. Keep PASS sites below the frequency gate.
#    LIVE, but note the field: the rarity oracle is the grpmax PROXY computed in
#    annotations.grpmax_af() as the max over vep_gnomAD{e,g}_{AFR,AMR,EAS,NFE,SAS}_AF.
#    There is no single faf95 INFO field to filter on (see allele_frequency.md), and
#    vep_MAX_AF is NOT a substitute — it over-counts founder groups.
bcftools view -f PASS -i 'INFO/vep_gnomADg_NFE_AF < 1e-4' "${TRIO_VCF}" -Oz -o "${RARE_VCF}"

# 2. TARGET (needs LOFTEE): restrict to the pLoF mask (HC, no flags) via the VEP CSQ field.
#    Under the current contract the equivalent live gate is simply IMPACT=HIGH.
filter_vep -i "${RARE_VCF}" \
  --filter "LoF is HC and LoF_flags is not defined" \
  -o "${PLOF_VCF}" --force_overwrite

# 3. TARGET (no callable-region BED exists): intersect so no-calls do not inflate counts.
bcftools view -R "${CALLABLE_BED}" "${PLOF_VCF}" -Oz -o "${QUAL_VCF}"

# 4. Emit gene, consequence, and sample for a per-gene tally.
#    NB: +split-vep does NOT accept --threads.
bcftools +split-vep "${QUAL_VCF}" -f '%CHROM\t%POS\t%SYMBOL\t%Consequence[\t%SAMPLE=%GT]\n' -d
```

```r
# De novo Poisson enrichment (SECONDARY arm), denovolyzeR. Shown as the reference
# implementation of the test Step 6 approximates in python when --mutrate is supplied.
# 'dnm' has columns: gene, class (lof/mis/syn); n_trios is the trio count.
library(denovolyzeR)
denovolyzeByGene(genes = dnm$gene, classes = dnm$class, nsamples = n_trios)
# QC gate: synonymous class enrichment ratio (observed/expected) must be ~1.0.
# NB: this synonymous check is exactly what Step 6 does NOT do — its expectation is
# raw 2 * N_trios * mu, unscaled by any observed synonymous rate.
```

## Recommended defaults (this pipeline)

| Parameter | Default | Status | Notes |
| --- | --- | --- | --- |
| **Primary signal** | **Recurrence-based gene consolidation** — distinct-individual carrier counts per gene per inheritance model | IMPLEMENTED | The driver. Inherited-germline focused |
| Recurrence flag | `burden.min_carriers` = **2** distinct individuals | IMPLEMENTED | |
| Recurrence null | `Binomial(N_trios, p)`, model-appropriate `p`; `absent_af_floor` = **1e-6** | IMPLEMENTED | Requires `--n-trios`; **skipped with a warning if omitted** |
| Secondary signal | **De novo enrichment** (denovolyzeR-style, Samocha-2014 rates, Poisson) | TARGET | Needs a `--mutrate` table; none ships. **Uncalibrated** when it runs |
| Corroborative signal | **TRAPD** vs gnomAD v4.1 exomes | NOT IMPLEMENTED | `burden.corroborative_trapd` is reserved |
| Frequency oracle / field | gnomAD **v4.1** point AF via the VEP cache, **grpmax proxy** (max over AFR/AMR/EAS/NFE/SAS) | IMPLEMENTED | **Not `faf95`** — no AC/AN in the cache. Never internal cohort AC/AN; never `MAX_AF` or global AF |
| Rarity gate (dominant / de novo candidate) | grpmax-proxy AF **< 1e-4** | IMPLEMENTED | `nhomalt` / absent-or-singleton gates are **retired** — no `nhomalt` field exists |
| Functional ladder | `IMPACT` ∈ {HIGH, MODERATE}, **else** `CADD_PHRED ≥ 25.3` | IMPLEMENTED | Two rungs, whole ladder. 25.3 is off-label on non-coding |
| Masks | **M1** = LOFTEE HC no-flag; **M2** = M1 ∪ (REVEL ≥ 0.5 or AlphaMissense LP or CADD ≥ 20) | REFERENCE | Neither computable: no LOFTEE, no REVEL/AlphaMissense |
| De novo genotype QC | GQ ≥ 20, DP ≥ 20, het AB 0.25–0.75, parent alt AD ≤ 1 | IMPLEMENTED | `hiConfDeNovo` screen when the tag is present |
| Exome-wide significance | **P < 2.5e-6** | IMPLEMENTED | ~0.05 / 20,000 genes |
| Discovery FDR | **BH q < 0.05** | IMPLEMENTED | Per model family (dominant / biallelic / X-linked / de novo) |
| Mask/tier combination | **ACAT / Cauchy (ACAT-O)** → one p-value | REFERENCE | Moot without multiple masks |
| Ranking weight | Constraint (LOEUF v2.1.1 / pLI / s_het / pHaplo) | IMPLEMENTED | Needs a `--constraint` table; without one every gene reads unconstrained |
| Standing QC | Synonymous **λ ≈ 1.0**; per-gene expected-vs-observed | NOT IMPLEMENTED | No synonymous calibration exists; `--syn-denovo-count` is reserved |
| Biobank frameworks | **Deferred** (SAIGE-GENE+, regenie) | REFERENCE | Require a joint call set |

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
