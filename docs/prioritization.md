# Prioritization (Step 9): gene excess statistic, artifact panel, variant tiering

Re-ranks the candidate calls a run produces, so a reviewer gets an ordered list with a readable
justification per row instead of a flat one. Two layers: a **gene-level test of excess over the
gene's mutational target** with a six-signal artifact panel and a graded down-weight, and a
**per-variant tier plus an additive points composite** in which every scoring term is its own
column.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the configurable
> defaults defined in [Canonical defaults](README.md#canonical-defaults).

## Status: what of this doc is live

| Layer | What it is | Examples here |
| --- | --- | --- |
| **IMPLEMENTED** | Step 9 (`pipeline/09_prioritize.py`, `src/hprv/prioritize.py`) does this today | The NB2 trimmed-fit excess statistic + mid-p calibration diagnostic; the six artifact signals; the four-tier down-weight with both ceilings; the V0–V4 variant ladder; the additive `priority_points` composite; both rankings |
| **TARGET** | Documented, configured, **not** running | V5 (needs NMD annotation), calibrated missense strength (needs REVEL), pLoF confidence (needs LOFTEE), ClinVar star gate, per-ancestry yield, single-site de-escalation |
| **REFERENCE** | Literature the design rests on, not a claim about this code | Tavtigian's point system as an ACMG classifier; PVS1 grading; the SVI one-predictor rule |

**The one thing to hold onto before reading any number below.** The excess statistic is a
**quality** question, not a biology question. It answers *is this gene producing more candidate
rows than its mutational target predicts?* It does **not** answer *is this gene
disease-associated?* A high excess is evidence of a **technical or population-genetic anomaly** —
mismapping in a segmental duplication, paralogue collapse, a callability defect, a founder allele,
an ancestry-uneven rarity gate, or a genuine hypermutable locus. A low excess is **not** evidence
that a gene is real. **The two directions must never be conflated.** A gene can be both an
established predisposition gene *and* a mismapping hotspot — 13 were, on the validation cohort —
and the correct response there is read-level review, not discarding the gene.

This mirrors, in the opposite direction, the honest caveat this repo already carries about its
Step-6 recurrence null ([gene_burden.md](gene_burden.md#multiple-testing-correction)):
*`p_recurrence` orders genes, it does not certify significance*. The excess statistic likewise
**orders genes by artifact suspicion**; it certifies nothing about disease.

## TL;DR

- **Step 6 nominates genes; Step 9 nominates variants.** Step 6 cannot say whether a gene emits
  more candidate rows than its mutational target predicts, and nothing downstream re-ranks the
  individual calls.
- **The gene statistic is excess over a mutational target:** `E_g = C · μ_g` where
  `μ_g = mu_mis + mu_syn + mu_lof` from gnomAD v2.1.1's per-gene Samocha targets, and `C` is a
  single scalar re-fit **per cohort**.
- **The null is negative binomial (NB2), not Poisson, and that is measured rather than assumed.**
  Pearson dispersion was **φ = 18.7** raw and **1.29** after the null's own trim (77 genes,
  0.392%); in a two-fold cross-fit the Poisson tail was **2.41× anti-conservative** at α = 1e-3
  *with its own arm trimmed* (1.82× untrimmed — the condition matters, see §2) while the NB was
  0.31×, i.e. conservative. Under the never-drop rule, conservative is the correct direction of
  error.
- **`C` must be fit over the FULL gene universe including the zero-count genes.** The candidate
  list is a zero-truncated sample; fitting `C` on matched genes only inflates `C` and deflates
  every ratio — the direction that hides artifact loci.
- **Six orthogonal artifact signatures**, from 195× enrichment (cohort saturation) down to 3.9×
  (low cumulative allele frequency), collapsed into an unweighted integer corroboration count.
  The statistic says a gene is anomalous; the panel is the **independent second witness** required
  before any penalty.
- **A four-tier graded down-weight, never a hard drop**, with an **auditable established-gene
  ceiling** and a **CDS-fallback ceiling**. With the `n_g ≥ 3` count floor on the ratio rules
  (default) it triaged **92 genes / 2,369 variants (9.57%) at 100% established-gene retention**,
  needing 8 named exemptions; without the floor, 228 genes / 2,565 (10.36%) and 13 exemptions. The
  floor keeps 92.4% of the volume using 40% of the genes — see §5.
- **The retention figure is now stress-tested, not just measured.** Expanding the positive-control
  union to 2,218 genes (adding ACMG SF v3.3 and two version-pinned PanelApp panels) changes it **not
  at all** — every added source maxes below the 5× cut. Control genes are *depleted* ~2× in the
  down-weight region.
- **A per-variant V0–V5 tier under ACMG/ClinGen SVI mechanism gating**: gene constraint counts
  only when the variant has a credible molecular effect, and a molecularly-benign prediction caps
  the total regardless of the gene.
- **An additive Tavtigian-style `priority_points` composite** where every term is a reported
  column — but the total is **not** an ACMG score and must never be read against Tavtigian's
  P ≥ 10 / LP 6–9 / VUS 0–5 bands.
- **Two rankings always ship**: fully agnostic, and prior-informed with an optional overlay that
  **defaults OFF**. `rank_delta` exposes exactly which calls a gene list promoted.
- **Never-drop is a hard invariant, asserted in code.** Row-count conservation is checked before
  the output is written, and the maximum penalty (−3.0) is deliberately too small to demote a
  variant carrying strong molecular evidence in a constrained gene.

## What this closes

Two open findings from `SCIENCE_AUDIT.md`:

- **A-3 (no calibration diagnostic).** Step 9 computes and audits a **mid-p calibration
  diagnostic** for the fitted null, for the NB *and* the Poisson side by side — so the choice
  between them is auditable rather than asserted. Under a correct null, mean(mid-p) ≈ 0.5 and
  P(mid-p < α) ≈ α; the numbers land in `audit/counts.tsv` as `calibration.*` on every run. This
  does **not** close the *other* half of A-3: there is still no synonymous-λ check and no
  positive-control **recovery** measurement on real data (see [Scope limitations](#scope-limitations-stated-honestly)).
- **A-2 (recurrence "significance" is a rank, not a calibrated test).** Step 9 does not repair
  Step 6's null, but it stops *compounding* it: the recurrence term scores on the **carrier
  count**, capped at +2, and never on `p_recurrence`. A p-value that saturates cannot order
  anything.

It also lands the "Gene-list tier priors" nice-to-have from [ROADMAP.md](ROADMAP.md) — as a
mechanism-gated **prior** in a **second, separately reported ranking**, never as a filter.

---

## The funnel, and what happens to a gene with no offset

![Accounting from raw emission to tier]({{artifact:art_87de9774-516f-42f2-869c-8154c366bcd7}})

**Figure 1. Accounting from raw emission to tier.** Left: 97.5% of variants reach a mutational
target; the 623 that do not are reported and kept, never dropped. Right: the four tiers of the
24,766 that do.

| stage | variants | genes |
|---|---|---|
| emitted by the screen | **25,389** | **10,799** |
| direct symbol join to gnomAD | 23,456 | 9,786 |
| recovered by Ensembl-gene-ID reconciliation | +1,310 | 591 symbols |
| **with a mutational-target offset** | **24,766 (97.54%)** | **10,377** |
| no offset — **reported, never dropped** | 623 (2.46%) | 422 |

**Note the count: 25,389 variants over 10,799 symbols, not 25,390 / 10,800.** The raw per-gene
counts file contains **one literal `gene` header line carrying n = 1**. A header row read as a gene
symbol is a phantom gene with a real count in every table it reaches, so `GENE_KEYS_LOWER` guards
every symbol-reading path in the module and a test asserts it. This is the kind of off-by-one that
survives review because every individual number still looks plausible.

**A gene with no offset gets `E_expected = ""`, `excess_ratio = ""`, `E_source = none`, `gene_tier =
T0` and a reason string saying so — never a fabricated ratio and never a dropped row.** 422 symbols
carrying 623 variants were in that class on the validation cohort, predominantly novel/renamed ORFs
and readthrough loci.

## The gene layer

### 1. The offset `E_g`

For each protein-coding gene *g* in the gnomAD v2.1.1 constraint table (19,643 genes):

```
mu_g = mu_mis_g + mu_syn_g + mu_lof_g
E_g  = C * mu_g              where  C = sum_g(n_g) / sum_g(mu_g)
```

`mu_*` are gnomAD's per-gene Samocha-model mutation-rate targets. `C` absorbs the average pass
rate of the whole hprv funnel — the rarity gate, the three-rung functional ladder, genotype QC —
and it is a **cohort property, re-fit every run**, never a constant.

**Both sums run over the FULL gene universe, not over the genes with `n_g ≥ 1`.** 9,266 of the
19,643 genes had `n_g = 0` on the validation cohort, and those zeros carry real information about
`C`. The arithmetic is exact: fitting on matched genes only multiplies `C` by
`n_universe / n_matched` and divides every excess ratio by the same factor. Poisson fit on the full
universe gave **C = 47,997**; the trimmed-NB fit gave **C = 42,712**.

**`mu_lof` is null for 505 of 19,643 genes.** They are neither dropped nor charged `mu_lof = 0`
(which would inflate their ratios ~5% for no reason). They are imputed at the median
`mu_lof / (mu_mis + mu_syn)` over the 19,138 genes that have it, and the imputation is recorded in
a `mu_lof_src` column so it stays visible.

**Genes with no constraint row at all** (418 symbols / 623 variants = 2.5% of the validation
candidate list — predominantly novel/renamed ORFs and readthrough loci) get a CDS-length
regression offset, `log(mu_tot) = −18.4484 + 1.0570·log(cds_length)`, fit by OLS over the full
universe (R² = 0.828; the 10th–90th percentile of predicted/actual spans 0.712–1.283, i.e. a
±30% band). Those genes are emitted with `E_source = cds_fallback` **and a tier ceiling**: a ±30%
offset error cannot support a 5× claim. Genes with neither `mu` nor a CDS length get
`E_expected = ""`, `excess_ratio = ""`, `gene_tier = T0` and a reason string saying so — never a
fabricated ratio.

**Symbol reconciliation, if you add it, must go through Ensembl gene IDs, not aliases.** Alias
matching produced demonstrably wrong joins on the validation data: `ACOD1` → `CAD`, `DRC3` →
`EPS8L1`, `EMSY` → `TNRC6A`, `TRDC` → `BCL11B`. Join on the versionless Ensembl gene ID and discard
any mapping whose target symbol already carries its own count, or you double-count. Reconciliation
recovered 595 symbols / 1,310 variants (5.2% of the list) on the validation cohort; Step 9 as
shipped does the direct join only and reports the unmatched count rather than guessing.

### 2. The null: negative binomial, and why not Poisson

```
n_g ~ NegBinomial(mean = C·mu_g, variance = mean·(1 + alpha·mean))
p_g = P(N >= n_g)
```

Overdispersion is large and must be modelled. Pearson dispersion about the Poisson mean over all
19,643 genes was **φ = 18.73** (Poisson requires 1). The honest reading is the **trimmed** figure:
after removing the **77 genes (0.392%)** the null itself identifies as outliers, **φ = 1.29** — so
residual overdispersion of the *bulk* is mild but real, and the raw 18.7 is dominated by ~0.4% of
genes.

> **On the trim count.** Phase 1's spec quoted 69 genes / 0.35%; that is the count at *iteration 2*.
> Iterations 3 and 4 both trim **77 (0.392%)**, which is why the loop terminates, and 77 is what
> the implementation reports. No downstream number moves — `C`, `α`, every `p_nb` and the whole
> tier table are identical — but the trim fraction feeds a **hard guard** (halt above 5%), so the
> figure a run reports has to be the converged one, not an intermediate.

The Poisson also mis-predicts the zero class: observed zero fraction **0.4717**, Poisson **0.4096**,
NB **0.4712**. The NB absorbs it, so **no zero-inflation component is needed — do not fit a ZINB.**

#### The calibration evidence

Two-fold cross-fit (seed 20260729): `(C, α)` fit on a random half of genes, mid-p evaluated on the
held-out half restricted to `excess_ratio < 5` — a *shape* criterion, chosen so the evaluation set
is not defined by the p-value being tested.

| Null | mean mid-p (ideal 0.5) | P(p < 0.05) (ideal 0.05) | P(p < 0.001) (ideal 0.001) | factor at α = 1e-3 |
|---|---|---|---|---|
| Poisson, **arm trimmed** | 0.5180 | 0.0449 | 0.00241 | **2.41× anti-conservative** |
| Poisson, arm **untrimmed** | 0.5477 | 0.0361 | 0.00182 | 1.82× anti-conservative |
| Negative binomial | 0.5122 | 0.0288 | 0.00031 | **0.31× (conservative)** |

> **The Poisson figure depends on how the Poisson arm is fit, so never quote it without that
> condition.** The trim loop belongs to **whichever null is being fit**: evaluating the Poisson at
> the NB's trimmed `C` measures a hybrid nobody would deploy, so Step 9 gives the Poisson arm its
> own iterative trim (`fit_poisson_null`) and reports the `C` it used as `calibration.C_poisson`.
> With the arm trimmed the anti-conservatism is **2.41×** (Phase 1 published 2.51× and did not
> state the condition; the residual gap is split-draw noise). Under a plain untrimmed reading it is
> **1.82×**, ranging 1.60–2.04× over 40 random splits. **The design decision is robust to the
> ambiguity** — every variant tested leaves the Poisson anti-conservative (all ≥ 1.31×) and the NB
> conservative (0.31×) — but the two numbers are different measurements and are not
> interchangeable.

**Read this in both directions.** The Poisson tail is 2.5× too liberal at α = 1e-3: it would
nominate ~2.5× as many "significant excess" genes as the false-positive rate allows. The NB is ~3×
**conservative** at the same α, i.e. it under-calls. Given the never-drop rule and the sensitivity
constraint below, **under-calling is the correct direction of error**: a missed artifact gene stays
in the review queue (cost: reviewer time), whereas an over-called one gets its variants penalised
(cost: a possible missed discovery). At BH q < 0.05 on the real data, Poisson flags 77 genes /
3,143 variants and the NB flags 52 / 2,618 — the 25 extra Poisson genes *are* the anti-conservatism.

A quasi-Poisson (φ = 1.29 on the trimmed bulk) is defensible but inferior: it has no proper
likelihood, so it cannot be compared by AIC against the covariate models below, and it imposes
variance linear in the mean whereas the data show the quadratic NB relationship. Step 9 reports φ
as a diagnostic alongside the NB rather than substituting it.

#### The trim: the tail must not calibrate its own null

Iterate to a fixed point (4 iterations on the validation cohort):

1. `keep` = all genes
2. `C = sum(n[keep]) / sum(mu[keep])`
3. `alpha` = NB dispersion MLE on `keep` at mean `C·mu[keep]`
4. `keep = { g : P(N ≥ n_g) > trim_p }`
5. repeat 2–4 until `keep` is stable

**Trimming matters enormously.** The untrimmed global fit gives `alpha = 0.804`, nearly 4× the
trimmed **0.214**, because the artifact tail is absorbed into the dispersion and then *hides
itself*. `used_in_null_fit` is emitted per gene so the fit is auditable. Guard: if the loop removes
more than `trim_max_fraction` (default 5%) of genes, Step 9 **halts** — that means the null is
mis-specified for the cohort (usually: the offset table and the candidate counts describe different
gene universes), not that 5% of the exome is artifact.

#### Covariates: significant, and deliberately not used

Tested by NB-GLM on the trimmed bulk with `log(C·mu)` as offset:

| model | ΔAIC vs offset-only | coefficient (p) |
|---|---|---|
| `+ log(cds_length)` | −5.4 | 0.029 (0.0075) |
| `+ log(num_coding_exons)` | −22.0 | 0.039 (1.1e-6) |
| `+ log(gene_length)` | −46.8 | 0.040 (4.0e-12) |
| `+ oe_syn` | −62.9 | 0.370 (1.3e-17) |
| `+ log(gene_length) + oe_syn` | **−121.1** | 0.044 (1.3e-14), 0.406 (8.0e-21) |

So there **is** real residual structure: candidate yield rises with genomic span (more intronic
territory the SpliceAI/CADD rungs can reach) and with `oe_syn` (a callability/model-fit proxy).
`covariate_adjust` nonetheless defaults **false**, for two empirical reasons:

1. **The adjustment is negligible where it matters.** Spearman(ratio, covariate-adjusted ratio)
   among candidate genes is **0.992**, and 502 of 532 genes at ratio ≥ 5 are shared. The
   disagreement is on borderline genes, not the tail: `OR4Q3` goes 545.7× → 506.9×.
2. **`oe_syn` is itself one of the six artifact signals.** Regressing it into the offset would
   *absorb* the signal we must report separately, so a callability-defective gene would look
   normal.

GC content and per-base coverage were **not testable** — neither is in the gnomAD constraint table
and no coverage file exists under the VEP-only contract. That is a gap, not a substituted proxy.

### 3. What `q_nb < 0.05` does and does not certify

**Does:** under an NB null whose dispersion was estimated from this cohort's own bulk, this gene's
candidate count is higher than its mutational target predicts, at a 5% false-discovery rate across
the exome. BH-FDR is computed across the **full universe** (including the zero-count genes, whose
p is 1) — correcting over called genes only would make `q` anti-conservative by ~`n_universe /
n_called`.

**Does not certify, explicitly:**

1. **Not disease association.** See the top of this page.
2. **Not a *cause*.** The test is agnostic between mismapping, paralogue collapse, a founder
   allele, a hypermutable locus, ancestry-uneven rarity gating (audit **A-4**) and contamination
   (audit **A-1**). The artifact panel is how you attribute a cause.
3. **Not independent of the qualifying filters.** `mu_g` is a *mutational* target, whereas `n_g`
   counts rows that survived the rarity gate, the functional ladder and genotype QC. `C` absorbs
   the average pass rate, but any gene whose *pass rate* differs from average — more intronic
   territory, so more SpliceAI/CADD rungs reachable — shows excess for that reason alone. This is
   precisely the `log(gene_length)` effect above, and it is a real limitation.
4. **Not stable at low counts.** Of the 532 genes at ratio ≥ 5, **135 had `n_g = 1`** — a single
   variant against E = 0.1 is 10× excess and means nothing. The NB FDR handles this correctly on
   its own (the minimum `n_g` reaching q < 0.05 was **5**, and the minimum ratio 7.15), which is
   why `q_nb` is preferred over a raw ratio threshold, and why `min_n_for_ratio_rule` exists as a
   knob for the ratio rules only.
5. **Not a joint-genotyped cohort test.** These are per-family VCFs, so `n_g` counts *emitted
   candidate rows*: a variant seen in three trios contributes three.

### 4. The six-signal artifact panel

High-excess reference set: `excess_ratio ≥ 10` (101 genes, 2,506 variants). Baseline: the other
10,276 candidate genes. Enrichment = (% high-excess) / (% baseline); p = one-sided Fisher exact.

| # | signal | definition | fold | Fisher p |
|---|---|---|---|---|
| 1 | `sig_saturation` | `n_g / n_trios ≥ 0.10` | **195×** | 2.8e-39 |
| 2 | `sig_segdup` | segdup ≥ 98% identity covering ≥ 10% of the gene span | **10.1×** | 1.2e-11 |
| 3 | `sig_family` | symbol matches an artifact-prone family regex | **6.5×** | 5.3e-9 |
| 4 | `sig_oe_syn` | `abs(oe_syn − 1) > 0.3` | **5.2×** | 7.0e-15 |
| 5 | `sig_constraint_flag` | gnomAD `constraint_flag` ∈ {`mis_too_many`, `syn_outlier`} | **4.0×** | 1.3e-5 |
| 6 | `sig_caf_low` | `classic_caf` below the candidate-gene 10th percentile | **3.9×** | 6.7e-15 |

All six are significant at p ≤ 1.3e-5. `corroboration_count` is their **unweighted** sum
(0–6) — deliberately unweighted, because weights fit on the same 101-gene tail they are used to
judge would not transfer to another cohort. The count is genuinely sparse: 82% of candidate genes
carry none. At ratio ≥ 10, `corrob ≥ 1` reaches **71.3%** vs 17.2% baseline (4.1×, p = 2.2e-32) and
`corrob ≥ 2` reaches 35.6% vs 5.5% (6.5×, p = 7.8e-20).

Three signals need their interpretation stated, because the obvious reading is wrong:

**`sig_saturation` is the strongest signal and the only offset-free one.** `per_trio ≥ 0.10` is ~2×
the 99th percentile of candidate genes (median 0.0091, p99 0.050), so only 35 of 10,377 genes
cleared it. Nine genes exceeded 0.5 variants per trio and `OR4Q3` reached **1.04 — more candidate
variants than there are trios**. A gene in which a tenth of all families carry a "rare" variant is
not delivering rare variants. It needs no mutation model, no external annotation and no gnomAD
join, so it is the most portable across cohorts and the least likely to fail silently — but it is
**never sufficient alone**: `TTN`, `SYNE1` and `DNAH11` clear it at only 1.9–3.7× excess, which is
why every tier rule pairs it with excess.

**`sig_oe_syn` departure is informative in BOTH directions, and they mean different things.**
`oe_syn < 0.7` (9.2×) means gnomAD saw *fewer* synonymous variants than expected → poor callability
or a mis-specified transcript model *in gnomAD*, so `mu_g` (our denominator) is unreliable and the
excess may be an artifact of the denominator. `oe_syn > 1.3` (4.5×) means *more* than expected →
paralogue collapse inflating gnomAD's own counts by the same mechanism that inflates ours. The
symmetric flag captures both at 5.2×; the raw value is emitted so the direction stays visible.

**`sig_caf_low`'s null handling is an explicit choice, and code and prose now agree.** The
percentile is taken over the genes that **have** a `classic_caf` value (cut 7.96e-6 on the
validation cohort); the genes with a **null** value are then **flagged**, because a gene for which
gnomAD reports no pLoF allele frequency at all is exactly the low-information locus this signal is
about. 1,122 genes flag. Phase 1's spec prose said "nulls read as 0", which read literally means
fill-then-percentile — cut 4.166e-6, 1,038 genes, and a threshold whose value depends on how many
genes gnomAD happens to have no record for. That is the defect: **missing data would silently move
the threshold applied to the measured data.** The 84-gene difference is the genes with `classic_caf`
between 4.169e-6 and 7.961e-6. **The tier table is identical either way** (228 genes / 2,565
variants under the old count floor) because no gene's T2/T3 status turns on those 84; the
corroboration distribution shifts slightly. Pinned as `signals.caf_low.include_null_as_flagged`.

**`sig_caf_low`'s direction is counter-intuitive and is not a frequency signal.** High-excess genes
have **far LOWER** cumulative pLoF allele frequency (median `classic_caf` 2.1e-5 vs 1.34e-4
baseline; 35.6% zero-or-null vs 7.9%), and the *high* decile is actually **depleted** (0.60×,
p = 0.94) — so this is not "common variants leaking through the rarity gate". A `caf` near zero
means gnomAD reports essentially no pLoF alleles, which happens for two very different reasons: the
gene is small/strongly constrained, or gnomAD's own coverage there is poor. The accompanying
statistic disambiguates it — median `exp_syn` is **27.8** in the high-excess set versus **126.7** at
baseline — so these genes are systematically *small or poorly covered in gnomAD*. When gnomAD has
little information about a gene, both `mu_g` (our denominator) and the grpmax rarity oracle (our
numerator's filter) degrade. It is a **low-information-locus flag**, and being partly a proxy for
gene size it is the weakest mechanism in the panel — so the strong tier never rests on it alone.

**Two things about `sig_segdup` that will bite an implementer.** First, **coordinate build is
load-bearing**: gnomAD v2.1.1 constraint `start_position`/`end_position` are **GRCh37/hg19**
(verified directly — `BRCA1` chr17:41,196,312–41,277,500 is hg19), even though the rest of this
pipeline is GRCh38-only. The validation overlap therefore used the hg19 UCSC `genomicSuperDups`
track. Mixing builds produces ~0 overlap **silently**. Second, the thresholds are **empirical, not
conventional**: any overlap is 1.9×, ≥98% identity at any coverage is 6.0×, ≥98% over ≥10% of the
gene is 10.1×, and ≥50% is 13.4× but halves the gene count — hence ≥10% as the operating point.

**And an important negative result that belongs in any methods text.** The excess statistic detects
segdup loci that produce *too many* calls; it is **blind to segdup loci that produce too few**,
because short-read callers fail to *emit* there rather than emitting garbage. `CYP21A2`, `SMN1`,
`SMN2`, `NCF1` and `STRC` all have **zero** candidates against non-trivial `E` at segdup fraction
1.00. So `segdup98_frac` is emitted on **every** gene regardless of tier, because it has a second,
independent use: marking genes where a *negative* result is uninformative. A reviewer reading "no
candidates in `CYP21A2`" needs to know not to believe it.

**`sig_family` is a curated regex list, and it is doing less work than it looks.** Per-family folds
where ≥ 2 genes hit: `LCE*` **102×**, `PRAMEF*` **41×**, `NBPF*` **34×**, `GOLGA*` **29×**,
`KRTAP*` 6.4×. But olfactory receptors (1.9×, p = 0.21) and zinc fingers (1.8×, p = 0.16) are **not
significantly enriched** on this data. They are retained because the individual offenders are
extreme (`OR4Q3` at 545×, `ZNF429` at 88×, `ZNF66` at 33×) and the mechanism — high paralogy, high
pseudogene content — is well established; but family membership is doing much less work than the
per-gene excess, which is why it enters only as **one vote**. The list is published as config,
versioned, and a family that recurs in the unexplained-excess review set is a candidate addition —
reviewed by a human, never auto-learned from the data being screened.

### 5. The four-tier down-weight

```
T1_watch:      excess_ratio >= 3  AND q_nb < 0.25
T2_downweight: excess_ratio >= 5  AND n_g >= 3 AND corroboration_count >= 1
T3_strong:     (q_nb < 0.05 AND corrob >= 1)                              # no count guard
            OR (excess_ratio >= 10 AND n_g >= 3 AND corrob >= 2)
            OR (excess_ratio >= 5 AND n_g >= 3 AND cohort_saturation AND corrob >= 2)
penalties:     T0 0.0 | T1 -0.5 | T2 -1.5 | T3 -3.0
```

**The `n_g >= 3` count floor on the ratio rules is `excess.min_n_for_ratio_rule`, and it matters.**
§4.4's own warning is that *a single variant against E = 0.1 is 10× excess and means nothing* — but
Phase 1 specified the guard and then never applied it to the T2 rule. Independent validation
measured what that omission cost, and the guard is now the default:

| rule set | genes | variants | % of 24,766 | control genes needing the ceiling | sensitivity |
|---|---|---|---|---|---|
| T2+T3 without the count floor (`min_n_for_ratio_rule: 1`) | 228 | 2,565 | 10.36% | 13 | 100% |
| **T2+T3 with `n_g ≥ 3` (DEFAULT)** | **92** | **2,369** | **9.57%** | **8** | **100%** |
| T2+T3 with `n_g ≥ 3`, ceiling disabled | 100 | 2,677 | 10.81% | — | 99.404% |

Without the floor, **136 of the 228 triaged genes (59.6%) have `n_g < 3` and carry only 196 variants
between them — 7.6% of triaged volume.** Most of the down-weighted gene *list* was singletons and
doubletons contributing almost none of the benefit while generating most of the exemption burden:
ten of the thirteen ceiling-exempted control genes had `n_g ≤ 4` and **eight had `q_nb` = 1**, i.e.
no statistical evidence of excess whatsoever. `SMPX`, `HMGA2`, `HAMP` and `PET100` were each **one
variant against E ≈ 0.16**. They never needed protecting; they needed a minimum-count requirement.

The floor keeps **92.4% of the triaged volume using 92 genes instead of 228** and cuts the exemption
burden from 13 control genes to 8, at unchanged 100% retention. It does **not** replace the ceiling —
`CTSA`, `NPRL3` and `CDH23` still require it — it shrinks what the ceiling has to excuse. Set
`min_n_for_ratio_rule: 1` to reproduce the original tier table exactly; both are measured on the same
data and both are reported here so the choice is visible rather than buried in a default.

**The FDR limb of T3 is deliberately unguarded.** No gene with `n_g < 5` reached `q_nb < 0.05` on the
validation cohort, so the NB tail already handles low counts correctly on its own; a floor there
would suppress a genuinely significant low-count gene for no measured benefit.

Assignment runs in ascending severity so the highest matching tier wins; the two ceilings apply
last. The rule shapes are each for a measured reason:

- **T1 requires no corroboration** because it carries no real cost (−0.5) and its purpose is to
  populate a reviewable watch list, including the unexplained-excess class. `q < 0.25` is
  deliberately loose.
- **T2 uses the ratio, not the FDR**, because at the moderate-excess boundary the FDR is
  count-driven: no gene with `n_g < 5` reached `q_nb < 0.05`, so an FDR-only rule would exempt
  every low-count artifact gene. Requiring one corroborating signal supplies the evidence the
  count cannot.
- **T3 has three disjuncts** because there are three distinct ways to be confidently anomalous:
  significant with a mechanism, extreme with two mechanisms, or saturating the cohort with two
  mechanisms. The third exists because saturation is offset-free and therefore survives an
  unreliable `mu_g`.

**Penalty magnitudes are on the Tavtigian-style point scale** used by the composite: −0.5 ≈ half a
supporting increment, −1.5 ≈ one supporting, −3.0 ≈ one moderate. They are sized so that **a T3
penalty cannot by itself demote a variant carrying strong molecular evidence in a constrained
gene**: V4 splice (+4) + strong rarity (+2) + constraint (+1) = +7, and −3 leaves +4, still well
above the median. That is intentional. **The artifact tier is a re-rank, not a veto**, and the
never-drop rule requires exactly that.

![Observed vs expected candidate variants per gene]({{artifact:art_5e318ad1-0911-4551-8a47-d3fda7fc3522}})

**Figure 2. Observed vs expected candidate variants per gene.** Log-log; the shaded region is the 5×
down-weight boundary. The artifact tail (upper left — high count against a small mutational target)
and the positive controls (clustered on the diagonal) occupy visibly different regions. `CTSA` and
`NPRL3` are established genes sitting **inside** the artifact tail — the case the established-gene
ceiling exists for, and the case that shows why the ceiling cannot be read as endorsing the calls.

![Artifact-signature enrichment by excess class]({{artifact:art_59ec9f6f-584f-4b43-9dfa-205609923276}})

**Figure 3. Artifact-signature enrichment by excess class.** All six signals are enriched in the
≥10× class; cohort saturation reaches 284× over the <3× baseline and needs no mutation model, making
it the most portable signal across cohorts. Enrichment is **not** monotone across all four classes
for every signal — `sig_segdup` dips to 0.96× in the 3–5× class before rising to 10.6× — so read the
figure as "enriched in the ≥10× class", which holds for all six, rather than as monotonicity, which
does not.

#### Observed effect on the validation cohort

| tier | genes | variants | % of 24,766 | median excess | max excess | control genes |
|---|---|---|---|---|---|---|
| `T0_no_downweight` | 10,102 | 21,498 | 86.80% | 1.47 | 15.15 | 1,319 |
| `T1_watch` | 47 | 703 | 2.84% | 7.66 | 148.35 | 21 |
| `T2_downweight` | 172 | 399 | 1.61% | 6.40 | 30.30 | **0** |
| `T3_strong_downweight` | 56 | 2,166 | 8.75% | 24.78 | 545.72 | **0** |

**228 genes carrying 2,565 variants (10.36% of the candidate list) are down-weighted, and zero
established-control genes are among them.** Note the shape: T3 is only 56 genes but carries 8.75%
of the list, while T2 is 172 genes carrying 1.61%. The volume reduction comes almost entirely from
a handful of extreme loci, which is the desired behaviour — the reviewer-facing benefit is
concentrated and auditable rather than diffuse. Largest T3 members: `OR4Q3` (n = 229, 546×),
`CTU2` (201, 155×), `SLC22A1` (186, 137×), `KRTAP9-1` (169, 448×), `ADCK5` (167, 102×).

### 6. Sensitivity: the established-gene ceiling

![Tier composition and control-gene distribution]({{artifact:art_b9737ba5-a9a8-4b27-821f-f9a1789d2ee8}})

**Figure 4.** (a) Every down-weighted gene carries ≥1 mechanism by construction, while 84% of the
candidate list carries none. (b) Family composition of triaged volume, with the **non-enriched**
`ZNF` class greyed — six `ZNF` genes are triaged carrying 122 variants, but `ZNF` family enrichment
is 1.0×, so those six are triaged on their own per-gene excess and not because the family is
artifact-prone. (c) Control genes concentrate below the 5× cut.

The positive-control set is a **phenotype-agnostic** union of "genes with established gene-disease
validity" — ClinGen validity **Definitive or Strong** ∪ ClinGen dosage **haploinsufficiency score
3** ∪ ClinGen **clinical actionability** (adult + pediatric) — **not** a phenotype panel. 2,218
genes, of which 1,340 (4,073 variants) appeared in the validation candidate list.

**This set is materially larger than a small cancer-gene spot-check, and it changed the
recommendation.** Within it, `excess_ratio` reaches a **maximum of 148.35** (`CTSA`, cathepsin A,
galactosialidosis, n = 160) with `NPRL3` at 77× — both established genes *and* extreme excess
outliers. The 99th percentile is 6.71, the 95th is 3.73, the median 1.220. So **a bare 5× or
FDR-only rule is not free: it triages real genes.**

| rule | genes | variants (%) | control genes triaged | sensitivity |
|---|---|---|---|---|
| `ratio ≥ 3 & corrob ≥ 1` | 539 | 3,615 (14.6%) | 30 | 97.76% |
| `q_nb < 0.05` alone | 52 | 2,618 (10.6%) | 5 | 99.63% |
| `ratio ≥ 5 & corrob ≥ 1` | 241 | 2,879 (11.6%) | 13 | 99.03% |
| `q_nb < 0.05 & corrob ≥ 1` | 36 | 2,413 (9.7%) | 3 | 99.78% |
| `q_nb < 0.05 & corrob ≥ 2` | 15 | 1,001 (4.0%) | **0** | **100%** |
| `ratio ≥ 10 & corrob ≥ 2` | 36 | 932 (3.8%) | **0** | **100%** |
| **RECOMMENDED (T2+T3 with the control ceiling)** | **228** | **2,565 (10.4%)** | **0** | **100%** |

**Note precisely which component buys what.** Requiring one corroborating signal takes sensitivity
from 99.63% to 99.78% (5 control genes triaged → 3); the ceiling closes the remaining gap to 100%.
The corroboration requirement and the ceiling are **complementary, not substitutes** — neither
alone reaches full retention at the recommended triage volume.

**Why an explicit ceiling rather than a stricter threshold.** Both routes reach 100%. The ceiling is
preferable because it is **auditable and reversible**: the exempted genes are *named*, they still
appear in the watch list with their full excess statistics, and a reviewer can disagree. A stricter
global threshold instead protects those genes by shrinking the down-weight set from 228 to 38 —
silently releasing **190 genes carrying 1,520 variants** back into the review queue, with no record
of the trade.

**Thirteen control genes were capped this way on the validation cohort** (`control_ceiling_applied
= true`): `CTSA`, `NPRL3`, `CYB5A`, `CDH23`, `DPM1`, `HNRNPDL`, `SMPX`, `HMGA2`, `HAMP`, `PET100`,
`ASAH1`, `SURF1`, `GJB2`. A further eight reached `T1_watch` **on their own merits, without the
exemption** (`PEX5` 14.9×, `DNAH12` 13.5×, `ACADM` 8.7×, `MICU1` 7.6×, `MCCC2` 6.7×, `PLG`, `SOS2`,
`DNAH11`) — do not conflate the two groups; only the 13 had a penalty withheld.

**The ceiling protects them from a score penalty; it does not mean their calls are correct.** `CTSA`
at 160 candidate variants across 220 trios (0.73 per trio, 148× its target) is not credible as
biology. The honest reading is that `CTSA` is *both* an established recessive-disease gene *and* a
technically problematic locus, and its calls need read-level review before any of them is believed.
Those genes are emitted with `review_flag = established_gene_high_excess` and belong at the top of
a dedicated review list. **Conflating "established gene" with "trustworthy call" is exactly the
failure mode the top of this page warns about.**

**Both gaps Phase 1 flagged are now closed, and closing them changed the sensitivity number not at
all — which is the useful result.** Phase 1 recorded PanelApp GE green genes as unreachable,
attributing repeated HTTP 403s to a rejection by the service. **That attribution was wrong**: the
failure signature was `CONNECT tunnel failed, response 403` — a *proxy* rejection at the sandbox
boundary, before any request reached Genomics England. With network access granted, both panels
fetched on the first attempt. The union is now **2,218 genes, 1,343 of them in the candidate list
carrying 4,085 variants**:

| added source | version pin | genes in list | max excess | T2/T3 without the ceiling |
|---|---|---|---|---|
| ACMG SF v3.3 | 84 genes; v3.2 + *ABCD1*, *CYP27A1*, *PLN* | 55 | 3.98× | **0** |
| PanelApp 243 green | Childhood solid tumours **v5.12** (2026-05-06) | 57 | 3.04× | **0** |
| PanelApp 259 green | Childhood solid tumours cancer susceptibility **v1.30** (2025-10-13) | 48 | 3.04× | **0** |
| GCT prior T1/T2 | `gct_gene_priors.tsv` | 13 | 3.51× | **0** |

**Every added source maxes out below the 5× T2 cut and triages zero genes.** The exempted set is
unchanged and sensitivity is unchanged, so **Phase 1's caveat that the number "will move when the
control set grows" is now tested and refuted — it does not move.** That was the largest stated risk
to the A-3 sensitivity claim, and it is closed. Control genes are in fact **depleted ~2× in the
down-weight region** (2.53% exceed 5× excess versus 5.13% of all candidate genes), which is the
result you would want and is not something the design assumed.

One residual provenance caveat: **ACMG SF v3.3 membership is recalled rather than machine-read**
from a versioned file. The count of 84 is verified against the policy statement (Lee 2025, PMID
40568962), but the 84 symbols themselves should be replaced with the published supplementary table
on the first production run. It contributes 55 genes at max 3.98× and triages zero, so no number
here depends on it.

### 7. The unexplained-excess class: flag, never penalise

Under the six-signal panel, **29 genes (148 variants, 0.6% of the list) had `excess_ratio ≥ 10`
with corroboration count 0**, two of them established controls (`PEX5`, `DNAH12`). They are emitted
with `gene_tier = T1_watch`, `review_flag = unexplained_excess`, and **no penalty**.

This class is the reason the panel exists and the reason the never-drop rule matters. An
unexplained 100× excess is a strong hint of an unmodelled artifact, but with no mechanism
identified and the never-drop rule in force, **the correct action is to look, not to down-weight.**

Note how much work the saturation signal does here: on a five-signal panel *without* it, the
uncorroborated set was 44 genes carrying **1,488** variants, including `SLC22A1` (186), `ADCK5`
(167), `PON1` (115) and `MBL2` (100). Adding cohort saturation moves those four — and the bulk of
the volume — into a corroborated tier, leaving a residue of 29 genes carrying 0.6% of the list.

---

## The variant layer

### 8. The tier ladder

A **screening/triage tier, not an ACMG classification.** hprv assigns no ACMG weight, and the
resource gaps below make several ACMG codes unimplementable.

| tier | name | definition | status |
|---|---|---|---|
| **V5** | high-confidence LoF | HIGH-impact pLoF + NMD-competent + LoF-mechanism gene | **UNREACHABLE** — see below |
| **V4** | strong splice / unresolved LoF | `spliceai_ds ≥ 0.5`; **or** a HIGH-impact pLoF with indeterminate NMD status (= every pLoF today) | IMPLEMENTED |
| **V3** | supporting splice / high-CADD missense | `spliceai_ds ≥ 0.2`; **or** missense with `cadd ≥ 25.3` (off-label) | IMPLEMENTED |
| **V2** | uncalibrated missense / in-frame indel | missense below the CADD cut; `inframe_insertion`/`inframe_deletion` | IMPLEMENTED |
| **V1** | non-coding / synonymous discovery rank | kept only via the CADD rung, with `spliceai_ds < 0.2` | IMPLEMENTED |
| **V0** | molecularly benign prediction | `spliceai_ds < 0.1` **and** `cadd < 15` **and** impact ∈ {LOW, MODIFIER} — **hard-caps the total** | IMPLEMENTED |

**No pLoF can currently reach V5, and that is a resource gap rather than a scoring choice.** Abou
Tayoun 2018 grades PVS1 by whether the predicted truncation triggers nonsense-mediated decay: a
nonsense or frameshift in the **last exon or the last 50 nt of the penultimate exon** escapes NMD
and drops from Very Strong, because the truncated protein may retain function or act by a different
mechanism. Singer-Berk 2023 shows empirically that applying such a framework materially reduces the
false-positive rate of predicted-LoF calls. `variants.tsv` carries **no exon number, CDS position
or transcript length**, so the test cannot be evaluated and every pLoF is capped at **V4** with
`nmd_status = INDETERMINATE`. **The fix is small, specific, and needs no new resource**: carry
VEP's `EXON`, `CDS_position` and the transcript's exon count / CDS length through Step 2 into
`variants.tsv`. This is the single highest-value addition to Step 9.

Additional unresolved pLoF caveats, flagged rather than silently ignored: no LOFTEE, so
`plof_confidence` is `UNAVAILABLE` everywhere and there is no HC/LC distinction or low-confidence
flag (ancestral allele, single-exon transcript, 50-bp rule, non-canonical splice); the last-exon
problem also affects `start_lost` (a downstream alternative start may rescue it); and a
`stop_gained` in a gene with no loss-of-function disease mechanism is not evidence at all.

**Canonical splice sites** (`splice_acceptor_variant` / `splice_donor_variant`, the ±1,2 positions,
`IMPACT=HIGH`) are treated as pLoF for tiering, subject to the same NMD caveat — a canonical-splice
variant in a final intron may produce an in-frame skip. Where `spliceai_ds` is available it is
reported alongside, because Walker 2023 is explicit that **a canonical site with a *low* SpliceAI
score deserves scrutiny, not automatic Very Strong.**

**A blank `spliceai_ds` is not evidence of no splicing effect.** hprv runs on the precomputed
SpliceAI set with backfill off by default, so insertions > 1 nt and deletions > 4 nt carry **no**
score. Step 9 emits an explicit `spliceai_status` ∈ {`scored`, `not_covered`} rather than letting a
null read as zero — the same blank-vs-zero trap as NHF, and it bites the indel classes most likely
to disrupt splicing. Consequently **V0 requires both scores to be PRESENT** and below their cutoffs;
letting a blank satisfy a benign rule would cap exactly the variants nobody scored.

**The missense CADD route is a discovery rank, not PP3 evidence, and the column says so.** 25.3 is
Pejaver 2022's missense-only PP3-supporting value, so in this *tiering* step it is at least applied
in its calibrated domain (unlike the *screen*, where it only ever reaches non-coding variants —
see [Canonical defaults](README.md#canonical-defaults)). But ClinGen SVI says commit to **one**
predictor chosen before seeing results, and the ClinGen-calibrated choice is **REVEL** (PP3
0.644 / 0.773 / 0.932; BP4 ≤ 0.290 / ≤ 0.183). So a V3-by-CADD tier carries
`missense_evidence_source = cadd_offlabel`, and **no graded missense strength is claimed**. A
calibrated missense tier requires REVEL (or AlphaMissense ≥ 0.564 / MPC ≥ 2), i.e. dbNSFP, i.e. a
change to the annotation contract — a resource question, not a design question.

In-frame indels get **V2**: no calibrated in-frame predictor is available and their mechanism
(in-frame loss of a domain vs neutral) is unresolvable from annotation alone. The indel length is
reported, since a large in-frame deletion is a different proposition from a 3 nt one.

### 9. Mechanism gating — the ACMG/ClinGen SVI principle

**Gene constraint and the gene-list prior are conditioned on molecular effect, never added
unconditionally.** This is the single most important structural rule in Step 9.

```
if variant_tier == V0:      gate = 0.0   # HARD: a constrained gene CANNOT rescue a benign
                                         #   prediction, and neither can list membership
elif variant_tier in {V1,V2}: gate = 0.5 # half credit for discovery-rank / uncalibrated effect
else:                        gate = 1.0  # V3-V5: a credible molecular effect earns the full term
if inheritance in {compound_het, hom_recessive, x_linked_recessive}: gate = 0.0
```

The recessive zeroing is the canonical defaults' own instruction: *"Do not down-weight recessive
candidates by pLoF constraint."* Symmetrically, do not **up-weight** them by it either — pLI and
LOEUF measure selection against heterozygotes and are not evidence about a biallelic candidate.

The V0 gate is what stops a gene-list prior from becoming confirmation bias. A gene-list membership
must never rescue a molecularly-benign prediction.

**One TARGET check is specified but not wired.** V5 eligibility should additionally require that
loss of function is an established mechanism for the gene (operationalised as ClinGen dosage
haploinsufficiency score 3 — a *phenotype-agnostic gene-mechanism annotation*, not a disease gene
list). It would only ever withhold a *promotion* to V5, never drop a variant. Since V5 is
unreachable anyway, this check is moot today.

### 10. Rarity, quality, and the three-state NHF rule

| strength | grpmax-proxy AF | points |
|---|---|---|
| `strong` | < 1e-5 | +2 |
| `moderate` | < 1e-4 | +1.5 |
| `supporting` | < 1e-3 | +1 |
| `permissive` | < 1e-2 | +0.5 |
| `unknown` | no eligible group reports it | +2 (treated as rarest) |
| `fail` | ≥ 0.05 (ClinGen BA1) | **−8, and caps the total** |

**Two caveats travel with this column.** (1) The oracle is a **point estimate, not `faf95`** — the
VEP cache carries no AC/AN, so the CI correction is unrecoverable at any price. Do **not** present
these bands as ACMG **PM2**; hprv applies PM2 at Supporting only, and PM2 is evidence, not the
gate. (2) **Audit A-4, a validity *and* equity problem**: the proxy is the **max** over
grpmax-eligible groups *regardless of the cohort's actual ancestry composition*, so in a
multi-ancestry cohort such as GMKF the effective stringency **varies with the proband's ancestry**
— a variant genuinely rare in the proband's own population but common in another group is charged
the other group's AF. The loss is **unevenly distributed across ancestry groups** and runs in the
silent false-negative direction. Step 9 flags `rarity_driven_by_single_group` when `max_af` exceeds
`grpmax_af` by more than 10×, so a reviewer can see when the rarity term was one population's
doing. **Per-ancestry candidate yield** is the run-level diagnostic that would expose the bias
systematically, and it is not computable without ancestry labels — a TARGET.

**Genotype QC** (proband `GQ ≥ 20`, `DP ≥ 10`, het AB 0.25–0.75, hom-alt AB ≥ 0.90) **flags and
penalises, never drops**. **Audit A-5** belongs in any methods text: these are genotype-refined
VCFs, so `GQ` derives from **posterior** probabilities (`PP`), not raw `PL`. A `GQ ≥ 20` cut on a
posterior is a different, generally more permissive filter than on a likelihood, because the
posterior has already been sharpened by a *family* prior — the same trio structure the inheritance
model is about to use. The classification is therefore **not statistically independent** of the
prior that produced the genotypes, and it biases toward cleaner-looking Mendelian patterns than the
raw data support. Spot-check top candidates against pre-refinement `PL`. **Parental GQ/DP/AB are
absent** from `variants.tsv`, so the leg that establishes trans phase for a comp-het cannot be
quality-assessed; that is flagged as `partner_leg_quality_unknown` rather than assumed to have
passed.

**Blank NHF is not 0.0 — three states, never two.**

| `nhf_status` | condition | effect |
|---|---|---|
| `clean` | screened; no member ≥ threshold over ≥ `min_reads` | no penalty |
| `flagged` | some member ≥ 0.5 over ≥ 5 reads | −3 + review flag |
| `not_screened` | **every member blank** | **0 — an explicit uncertainty flag, no penalty and no credit** |

- **blank = NOT SCREENED** — the member is not an ALT carrier, has no mini-CRAM, or Step 8b never
  ran (no kraken2 DB).
- **`0.0` = SCREENED and every read classified human.**

**Treating blank as clean silently promotes exactly the calls nobody examined.** The `min_reads`
floor of **5** is essential and is hprv's own default: an NHF of 1.0 over 2 reads is noise, not
evidence, so the `*_nhf_reads` denominator always rides beside the fraction. (Logic copied from
`prepare_igv_variants.py:nhf_status`.)

**Audit A-8, stated plainly:** NHF **annotates, it does not filter**. `hprv_summary.xlsx` and
`genes.ranked.tsv` are written *before* Step 8b and are **not** contamination-aware. Do not confuse
NHF (per-call, per-member ALT-read composition) with Step-0 FREEMIX (sample-level cross-human
contamination, audit **A-1**); neither filters the ranked output.

### 11. Inheritance-model coherence — `unknown` must be exactly neutral

`moi_coherent` / `moi_discordant` / `moi_unknown`, from the observed mode against a curated MOI.
Discordance is a **small demotion (−1) and a review flag, never a filter** — incomplete
penetrance, mosaicism, a second undetected hit and a genuinely novel mechanism all produce
discordance. **`moi_unknown` scores exactly 0**: any penalty there converts the score into a
known-gene filter and destroys novel-gene discovery, violating the never-drop rule in spirit if not
in letter.

Three audit caveats mean the observed mode is not face value:

1. **A-6: mode assignment is single-gene-keyed and long genes drift to `compound_het`.** A
   phase-confirmed mat×pat pair consumes both legs, so a genuinely dominant-grade variant is
   relabelled `compound_het` whenever the child carries *any* other sub-1e-2 functional het in the
   same gene — near-certain in long genes (`TTN`, `NEB`, `RYR1`, `DMD`). Step 9 therefore
   **suppresses the discordance penalty entirely for a `compound_het` call in a long gene** and
   reports `moi_caveat = long_gene_comphet_drift` instead.
2. **A-6, second residual:** comp-het pairing keys on the single VEP-PICK'd gene, so two damaging
   hits in one gene fail to pair when the picked block names a different overlapping gene for one
   of them — and an unpaired second hit between the 1e-2 and 1e-4 gates is emitted under **no mode
   at all**. A missing mode is not evidence of anything.
3. **A-7: `dominant` carries no penetrance evidence.** The trios file is sample IDs only, with no
   affected status, and Step 5 never reads phenotype. A `dominant` call means "the proband is het
   and a parent transmitted it" — nothing more. Since the transmitting parent is presumably
   unaffected, every such call implicitly assumes incomplete penetrance, which is never tested, and
   there is no co-segregation (PP1/BS4) evidence available. Frame the inherited arm as nominating
   **candidate moderate-penetrance alleles** that cannot, on genetic grounds alone, be separated
   from benign inherited variants.

### 12. The composite: `priority_points`

An interpretable **additive-points** score in the spirit of Tavtigian 2020. Every term ships as its
own column.

| term | condition | points |
|---|---|---|
| `pts_molecular` | V5 / V4 / V3 / V2 / V1 / V0 | +8 / +4 / +2 / +1 / +0.5 / 0 |
| `pts_rarity` | see the rarity table above | +2 … +0.5, BA1 −8 |
| `pts_gene_constraint` | `pLI ≥ 0.9` or `LOEUF < 0.35`, **× the mechanism gate** | +1 |
| `pts_recurrence` | ≥ 2 distinct-variant carriers / ≥ 3 / same-variant | +1 / +2 (cap) / +0.5 |
| `pts_quality` | GT QC fail / NHF flagged / NHF not screened / comp-het partner unknown | −2 / −3 / **0** / −0.5 |
| `pts_clinical` | ClinVar P/LP (not conflicting) / conflicting or VUS / B or LB | +4 / 0 / −4 |
| `pts_moi` | discordant / **unknown** / coherent | −1 / **0** / 0 |
| `pts_gene_artifact` | T0 / T1 / T2 / T3 | 0 / −0.5 / −1.5 / −3.0 |
| `pts_gene_list_prior` | in the optional Class-B overlay, **× the mechanism gate** | +2, **default OFF** |

**Why additive, in three points that are each a design constraint.** *Auditability*: a reviewer must
see `spliceai=+4, rarity=+2, constraint=+1, gene_artifact=−3 → +4`, not `score=0.71`. *Additivity
survives a missing term honestly*: under the VEP-only contract several terms are simply unavailable
(no REVEL, no LOFTEE, no ClinVar stars), and an additive scheme degrades to "that term contributed
0", whereas a multiplicative or learned composite would silently redistribute the missing evidence
onto the terms that remain. *No training data exists* (audit **A-3**: sensitivity/precision are
*unmeasured*, not measured-and-acceptable), so a fitted model would be fit on the very artifacts
we are trying to remove.

**Where the analogy stops, and this must be in the methods.** These are **not ACMG points** and the
total **must not** be read against Tavtigian's classification bands (P ≥ 10, LP 6–9, VUS 0–5). The
criteria here are not ACMG criteria (a CADD-based term is not PP3); no phenotype, segregation or
functional evidence exists at all; the ClinVar term has no review-status gate; and the
artifact-penalty terms have no ACMG analogue whatsoever. **The column is named `priority_points`,
never `acmg_points`, and no P/LP/VUS label is ever emitted from it.**

**The recurrence term is capped at +2 and scores on the carrier count, never on `p_recurrence`.**
`gene_burden.md` states that the recurrence null is a **case-only approximation** built only from
variants observed in the cohort, so `p` is **too small** and *"for essentially any gene with ≥
min_carriers carriers of private variants, `p_recurrence` clears the exome-wide line"*. Audit **A-2**
records that this was once overclaimed as a calibrated test. A p-value that saturates cannot order
anything. Same-variant recurrence gets strictly less credit than distinct-variant, because a single
recurrent site shared by many trios is as easily a mapping/caller artifact or a founder allele as a
burden signal. **Audit A-1 compounds here**: one contaminated proband contributes false carriers
across many genes, and two carriers of distinct private variants give p ≈ 3e-7 at N ≈ 200 — so one
contaminated sample plus one genuine carrier can manufacture an "exome-wide significant" gene.
FREEMIX-failing trios must be excluded **upstream** from `TRIOS_FILE`, not down-weighted here.

**The ClinVar term has no star gate.** The VEP cache carries no `CLNREVSTAT`, so a 1★
single-submitter assertion is indistinguishable from an expert-panel one and is honored
identically. `clinvar_review_status = UNAVAILABLE` is emitted so the +4 is never mistaken for a
≥2★ assertion. The release is pinned by the cache (VEP 115 ⇒ ClinVar 2025-02), not independently.

**Two hard caps**, both implementing mechanism gating — molecular benignity and BA1-level frequency
are statements about the *variant* that no amount of gene-level enthusiasm can overturn:

```
if variant_tier == V0:         priority_points = min(priority_points, 0)
if rarity_strength == 'fail':  priority_points = min(priority_points, -4)   # BA1
```

**Neither cap removes the variant from the output.**

### 13. The two rankings

```
priority_points_agnostic = sum of every term EXCEPT the gene-list prior
priority_points_prior    = priority_points_agnostic + pts_gene_list_prior
rank_delta               = rank_prior - rank_agnostic
```

Both ship on every run with integer ranks, ties broken by `pts_molecular` desc, then `pts_rarity`
desc, then `pts_gene_artifact` desc (cleaner genes first), then `chrom`/`pos` for determinism.

**A large negative `rank_delta` means the variant was promoted purely by gene-list membership** —
exactly the set a reviewer should scrutinise for confirmation bias, and exactly the set that would
be invisible under a single blended ranking.

**The overlay contract is what keeps this step Class A / phenotype-agnostic:** the overlay is a
**file path in config**, never a list in code; this module names no gene;
`composite.gene_list_prior.enabled` defaults **false**, and a config path is **inert** while it is
false (only an explicit `--gene-prior` overrides that, so a stale path cannot silently start
promoting genes); with the overlay off the two rankings are byte-identical (asserted in
`tests/test_pure.py` and in the integration run); the prior is mechanism-gated exactly like
constraint; and it is a **prior, never a filter** — genes off the list lose the points and nothing
else.

### 13a. Four hazards in the overlay layer, each guarded in code

**1. The join is on gene symbol ONLY — deliberately MOI-agnostic.** Routing a prior lookup by a
gene's canonical mode of inheritance would silently miss every hypothesis stated about a *different*
genetic model than the gene is curated under. The concrete case: the FA/HR genes (`FANCA`,
`FANCD2`, `SLX4`, `FANCE`, `BRCA2`) carry germ-cell-tumour evidence about **heterozygous carriers**
(five-gene combined OR 10.17, 95% CI 4.87–21.27, P = 2.90e-06; Zhang 2025, PMID 40906985) while
their canonical Mendelian model — and their PanelApp green status — is **biallelic Fanconi anemia**.
An MOI-routed lookup would consult those rows under a recessive model and never apply them to the
het observations the evidence is about.

Symmetrically, **an observed-mode-vs-canonical-MOI mismatch is reported, never penalised.** A het
observation in a canonically-recessive gene emits
`moi_caveat = moi_mismatch_het_in_recessive_gene` and charges **zero** — it is the carrier-risk
shape, not an incoherent call, and incomplete penetrance, mosaicism, an undetected second hit and a
genuinely novel mechanism produce the same shape. This is the never-drop rule applied to the
coherence layer. (Where an overlay states the carrier hypothesis explicitly in its own MOI
vocabulary — `AR_biallelic_FA;heterozygous_carrier_risk_proposed` — the observation reads
*coherent* rather than as a suppressed mismatch.)

**2. Gene-level and set-level priors combine by MAX, never SUM.** A curated overlay may carry both
per-gene rows *and* a pathway-collapsed entry: the GCT resource has per-gene FA rows and
`FA_HR_PATHWAY_23` (23-gene FA/HR set, pooled OR 4.14, 95% CI 1.98–8.66, P = .0013). **Both derive
from the same study**, so adding them would double-count one study. The resource's own usage
contract says this; Step 9 enforces it in `parse_gene_prior_overlay` and asserts it in a test,
because a documented invariant is one nobody checks. `gene_list_prior_set_applied` names the set
when a set-level prior beat the gene-level one.

**3. `prior_weight` is UNCALIBRATED and must never be presented as a likelihood ratio.** It is an
ordering default (T1 1.0; T2 0.75 replicated / 0.60 single-study; T3 0.35; T4 0.15) that scales the
configured maximum. It is not an odds ratio, not a posterior, and not calibrated against anything.
`gene_list_prior_weight` reports the raw value so a reviewer can see what scaled the term.

**4. Some overlay rows are NOT germline evidence, and must contribute zero.** A curated overlay
lists somatic drivers deliberately — `KRAS`, `NRAS`, `CBL`, `MTOR`, `AKT1`, `BCORL1` in the GCT
resource, all `evidence_class = somatic_driver_not_germline` — precisely so a reader can see they
were considered and excluded. **Their weight of 0.15 is a bookkeeping placeholder, not weak germline
support**, and scoring it as the latter is exactly the misreading the resource warns against. Those
rows contribute **0.0**, carry `gene_list_prior_excluded_non_germline`, and are **reported rather
than dropped** — silently dropping them would hide that they were considered. Set membership does
not resurrect them either.

**Two things to know about the tiers themselves before reading an overlay-informed ranking.** T3
(59 of the GCT resource's 96 rows) is **GWAS-locus membership**, where the lead variants are
predominantly **non-coding** — the intracranial `BAK1` signal is a 4-bp enhancer deletion — so **a
rare-coding screen is not interrogating that mechanism at all**, and a T3 hit means "this gene sits
at a GWAS locus", not "rare coding variants here matter". And T2 is not uniform evidence: check the
`replication` column, since only `CHEK2` is independently replicated.

---

### 13b. Activating the gene-list prior (a phenotype panel)

The overlay is **off by default and stays off unless you turn it on**, so this is the recipe. Three
steps, none of which need a rebuild.

**1. Write the list.** Three accepted shapes, all matched on **gene symbol only**:

A bare symbol list (`#` comments and blank lines are skipped):

```
# my pediatric-cancer panel, v2026-07
BRCA1
SDHB
TP53
```

Or a table with a header naming a `gene` (or `symbol` / `gene_symbol`) column — **tab or
comma delimited**, so a spreadsheet export works as-is. Every column below is optional, and
**unknown columns pass through untouched** so your list can carry its own provenance:

```
gene	prior_weight	tier	evidence_class	moi	pmids	replication
BRCA1	0.9	GREEN	germline_predisposition	AD	12345678	replicated
SDHB	0.6	AMBER	germline_predisposition	AD	23456789	single_study
KRAS	0.15	RED	somatic_driver_not_germline	NA	34567890	na
```

Optionally, a **sibling `.json`** (same path with the extension swapped: `panel.tsv` →
`panel.json`) adds set-level priors, which are combined with gene-level ones by **MAX, never SUM**:

```json
{"gene_sets": {"FA_HR_PATHWAY": {"prior_weight": 0.6,
                                 "members": ["FANCA", "FANCD2", "SLX4", "FANCE", "BRCA2"]}}}
```

**2. Point the config at it and enable it.** Both are required — a path alone is inert:

```yaml
prioritization:
  composite:
    gene_list_prior:
      enabled: true
      path: ${GENE_PRIOR_OVERLAY}      # or a literal path; export the var, or hardcode it
```

`prepare_resources.sh emit-env` emits a commented `GENE_PRIOR_OVERLAY` line for you to fill in;
the overlay itself is **not** a fetched resource — it is yours, and lives outside this repo.

**3. Confirm it fired.** Step 9 says so on stderr, and `audit/counts.tsv` records it:

```
--gene-prior overlay: 47 genes from panel.tsv (+1 gene set(s) from the JSON sidecar, 4 members
  contributed by set membership alone; set and gene priors combined by MAX, never sum)
  overlay: 47 genes; 12 variants promoted by list membership (rank_delta < 0)
```

Then sort `igv/variants.prioritized.tsv` (§14) on `rank_prior`, and read `rank_delta` to see what
the panel changed. `gene_list_prior_member` / `_tier` / `_weight` / `_evidence_class` /
`_excluded_non_germline` / `_set_applied` carry the per-variant provenance.

**Two failure modes worth knowing, because both are now loud rather than silent:**

- A **headerless table** (`BRCA1,0.9,GREEN` with no header row) is a **hard stop**, not a
  degrade-with-a-warning. It cannot be told apart from a symbol list whose symbols happen to
  contain commas, and reading it that way would produce a plausible gene *count* over a list that
  matched nothing — the reviewer would ship an un-prioritised list believing the panel applied.
  Every other optional resource degrades because its absence is honestly reportable; this one
  isn't. Add a header row, or strip the file to one symbol per line.
  (`tests/test_pure.py:test_prioritize_gene_prior_overlay_accepts_csv_and_rejects_headerless_table`)
- A row whose `evidence_class` is in `non_germline_classes` (default
  `[somatic_driver_not_germline]`, matched as a substring) contributes **exactly 0.0** and is
  reported with `gene_list_prior_excluded_non_germline`. Note the string is
  `somatic_driver_not_germline`, **not** `somatic_driver` — a near-miss class name silently gets
  the full prior instead of zero, so check the stderr line's excluded count against what you
  expect. Pediatric-cancer panels routinely list drivers (KRAS, NRAS, CBL, MTOR, AKT1, BCORL1);
  they stay visible and credited nothing.

And two things the prior deliberately **cannot** do: it cannot promote a variant with no molecular
evidence (the mechanism gate of §9 still applies), and it cannot remove anything — `rank_agnostic`
ships on every run alongside `rank_prior` precisely so the panel's effect is measurable rather than
baked in. `prior_weight` is an **uncalibrated ordering default**, never a likelihood ratio.

### 14. The igv.js review table — where the triage actually gets used

The point of the two layers above is to stop a reviewer curating thousands of variants by hand.
That only happens if the tiers reach the tool the reviewing is done in, and
`variants.prioritized.tsv` **cannot** serve that role for two independent reasons:

1. It is written to the work-dir root, while the igv.js server's data dir is `igv/`.
2. Its column set omits `child_/mother_/father_file`, `*_index` and the `*_vcf*` columns — the
   per-member mini-CRAM and VCF track paths, which are **relative to `igv/`**. Pointing the server
   at it would give a sortable list with no read-level view at all, which is the one thing Step 8
   exists to provide.

So Step 9 takes `--out-igv-variants` and writes a third table, **`igv/variants.prioritized.tsv`**:
every column of the input verbatim and in its original order, plus every prioritization column
appended. `run_pipeline.sh` passes it only when the input actually was `igv/variants.tsv` — built
from Step 5's calls there would be no tracks to point at, and putting a track-less table inside
`igv/` would merely be misleading.

Three properties make it a drop-in replacement for the reviewer's variants file:

- **Input columns are never overwritten.** The appended set is the strict complement of the input
  header, so Step 8's table stays authoritative for everything it already reports and this file
  only ever ADDS. Asserted byte-for-byte in
  `tests/test_pure.py:test_prioritize_igv_review_table_preserves_track_paths`.
- **Rows are carried by position, not re-joined.** A `chrom/pos/ref/alt/trio_id` join is ambiguous
  for two ALTs of one multiallelic site in one trio; positional identity cannot go wrong. The test
  fixture contains exactly that decoy pair.
- **Never-drop is asserted again here.** A silently shortened review list is the failure mode
  nobody notices, so row-count conservation is re-checked on this file rather than inferred from
  the loop that built it.

It is sorted by **`rank_agnostic`**, so the file opens honest even when a `--gene-prior` overlay is
loaded; `rank_prior` and `rank_delta` are columns the reviewer sorts on in the UI, which keeps the
"what did the phenotype list change?" question one click away instead of a separate run.

Everything past the required `chrom/pos/ref/alt` is filterable in igv.js, so the useful review
columns are: `variant_tier`, `priority_points_agnostic`, `rank_agnostic`, `rank_prior`,
`rank_delta`, `gene_tier`, `downweight_reason`, `review_flag`, `excess_ratio`, `q_nb`,
`nhf_status`, `moi_coherence`, `clinvar_strength`, and the overlay provenance
(`gene_list_prior_member` / `_tier` / `_weight` / `_evidence_class` /
`_excluded_non_germline` / `_set_applied`) alongside Step 8's own `gene`, `consequence`, `impact`,
`grpmax_af`, `max_af_pops`, `cadd`, `spliceai_ds`, `clin_sig` and the NHF columns.

Because a filter in igv.js is a **view**, not a deletion, filtering here does not violate
never-drop: the file on disk still carries every candidate, and the reviewer can always widen back
out. That is the whole reason the down-weight sets a tier and a reported penalty rather than
removing a row — the decision of what not to look at stays with the reviewer, and stays reversible.

## Recommended defaults (this pipeline)

| Parameter | Default | Status | Source / notes |
| --- | --- | --- | --- |
| Null model | **negative_binomial** | IMPLEMENTED | Poisson is 2.41× anti-conservative at α = 1e-3 (arm trimmed; 1.82× untrimmed — see §Calibration); NB is 0.31× (conservative), the right direction under never-drop |
| Offset | `mu_mis + mu_syn + mu_lof`, `E = C·mu` | IMPLEMENTED | Karczewski 2020 (PMID 32461654) |
| `C` fit universe | **the FULL table, zero-count genes included** | IMPLEMENTED | The candidate list is zero-truncated; matched-only inflates C |
| `mu_lof` imputation | ×(1 + **0.0516**) when null | IMPLEMENTED | Median `mu_lof/(mu_mis+mu_syn)` over 19,138 genes |
| CDS fallback | `log(mu) = −18.4484 + 1.0570·log(cds_length)` | IMPLEMENTED | OLS, R² = 0.828, ±30% band; never reaches T3 |
| `excess.trim_p` | **1e-3** | IMPLEMENTED | Untrimmed α = 0.804 vs trimmed 0.214 — the tail hides itself |
| `excess.trim_max_fraction` | **0.05** (HALT) | IMPLEMENTED | A >5% trim is a mis-specified null, not a 5%-artifact exome |
| `excess.q_threshold` | **0.05** | IMPLEMENTED | BH-FDR over the full universe |
| `excess.min_n_for_ratio_rule` | **3** | IMPLEMENTED | The spec's own guard, now applied to the T2 rule and both ratio limbs of T3. 92 genes / 9.57% triage, 8 exemptions; set to 1 for the original 228-gene / 10.36% table with 13 exemptions. 100% retention either way |
| `signals.caf_low.include_null_as_flagged` | **true** | IMPLEMENTED | Percentile over genes that HAVE a value (cut 7.96e-6); nulls flagged as low-information loci. Fill-then-percentile would let missing data move the threshold applied to measured data |
| `composite.gene_list_prior.honor_prior_weight` | **true** | IMPLEMENTED | Scales the prior by the overlay's `prior_weight`. **UNCALIBRATED** — an ordering default, never a likelihood ratio |
| `composite.gene_list_prior.non_germline_classes` | **[somatic_driver_not_germline]** | IMPLEMENTED | Those rows contribute **0.0** and are reported, not dropped — a 0.15 weight on a somatic driver is bookkeeping, not weak germline support |
| `composite.gene_list_prior.combine_gene_and_set` | **max** | IMPLEMENTED | Gene-level and set-level priors combine by MAX, never SUM (both can derive from one study) |
| `covariate_adjust` | **false** | IMPLEMENTED | ΔAIC −121 but Spearman 0.992; `oe_syn` is itself a reported signal |
| `signals.saturation.per_trio_min` | **0.10** | IMPLEMENTED | 195× enriched, p = 2.8e-39; ~2× the 99th percentile |
| `signals.segdup.min_identity` / `min_frac` | **0.98** / **0.10** | IMPLEMENTED | 10.1×; empirical operating point (≥50% is 13.4× but halves the gene count) |
| `signals.oe_syn.max_deviation` | **0.30** | IMPLEMENTED | 5.2×; symmetric — both directions are informative |
| `signals.constraint_flag.values` | **mis_too_many, syn_outlier** | IMPLEMENTED | 4.0×; `mis_too_many` alone is 5.1×, `no_exp_lof` deliberately excluded |
| `signals.caf_low.percentile` | **0.10** | IMPLEMENTED | 3.9×; a LOW-INFORMATION-locus flag, not a frequency flag |
| `signals.family.patterns` | 13 regexes | IMPLEMENTED | 6.5×; curated, versioned, incomplete by design; one vote only |
| `T1_watch` | `ratio ≥ 3 & q < 0.25`, **−0.5** | IMPLEMENTED | No corroboration required — a watch list, not a penalty |
| `T2_downweight` | `ratio ≥ 5 & corrob ≥ 1`, **−1.5** | IMPLEMENTED | Ratio not FDR: no gene with n < 5 reaches q < 0.05 |
| `T3_strong_downweight` | three disjuncts, **−3.0** | IMPLEMENTED | Three distinct ways to be confidently anomalous |
| `established_gene_ceiling` | **T1_watch** | IMPLEMENTED | 100% control retention at 10.4% triage; auditable and reversible |
| `min_control_genes` | **1000** (HALT) | IMPLEMENTED | A truncated union makes the ceiling silently empty |
| `max_downweight_fraction` | **0.20** (HALT) | IMPLEMENTED | Measured value was 10.4% |
| SpliceAI V4 / V3 bands | **0.5** / **0.2** | IMPLEMENTED | Walker 2023 (PMID 37352859), ClinGen SVI |
| V0 benign rule | `spliceai < 0.1` **and** `cadd < 15` **and** LOW/MODIFIER | IMPLEMENTED | BOTH scores must be PRESENT — absence is not benignity |
| Missense CADD | **25.3**, labelled `cadd_offlabel` | IMPLEMENTED | Pejaver 2022 (PMID 36413997); REVEL is ClinGen's calibrated choice |
| Rarity bands | 1e-5 / 1e-4 / 1e-3 / 1e-2, BA1 0.05 | IMPLEMENTED | grpmax **proxy**, not `faf95`; not ACMG PM2 |
| Mechanism gating | V0 **0.0**, V1/V2 **0.5**, V3–V5 **1.0** | IMPLEMENTED | ACMG/ClinGen SVI; also zeroed for recessive modes |
| NHF states | clean / flagged / **not_screened** | IMPLEMENTED | Blank ≠ 0.0; `not_screened` scores 0, neither penalty nor credit |
| `moi` unknown | **0 — exactly neutral** | IMPLEMENTED | Any penalty converts the score into a known-gene filter |
| Recurrence cap | **+2**, on the carrier count | IMPLEMENTED | The Step-6 null is case-only and saturates (audit A-2) |
| ClinVar P/LP | **+4**, `review_status = UNAVAILABLE` | IMPLEMENTED | No `CLNREVSTAT` in the cache — no ≥2★ gate possible |
| `gene_list_prior.enabled` | **false** | IMPLEMENTED | Class-B overlay; with it off the two rankings are identical |
| V5 (NMD-competent LoF) | +8 | **TARGET** | Needs VEP `EXON`/`CDS_position`/transcript length in `variants.tsv` |
| Calibrated missense (REVEL) | PP3 0.644/0.773/0.932 | **TARGET** | Needs dbNSFP — a contract change |
| pLoF confidence (LOFTEE HC/LC) | — | **TARGET** | LOFTEE data not bind-mounted |
| Single-site de-escalation | `max_site_share ≥ 0.5 & n_trios ≥ 5` | **TARGET** | Columns emitted; the de-escalation rule is not wired |
| Per-ancestry candidate yield | — | **TARGET** | No ancestry labels; the audit A-4 diagnostic |

All values are configurable defaults in `config/config.example.yaml`, not immutable law. A
gene-specific ClinGen VCEP threshold overrides any generic cutoff here.

---

## Scope limitations (stated honestly)

- **The gene layer is validated on real data; the variant layer is not.** The excess statistic, the
  six signals, the tier assignment and the 100%-established-gene-retention result are all measured.
  **No weight in the composite has been tested against data** — they are reasoned from Tavtigian's
  doubling scale and from which evidence classes are calibrated, not fit, and they need review once
  a real run exists.
- **Audit A-3 is only half closed.** Step 9 adds the mid-p calibration diagnostic. It does **not**
  add the **synonymous-λ** check or a positive-control **recovery** measurement. The synonymous
  burden through the identical rarity + QC path is the single highest-value remaining addition: if
  λ ≫ 1 the qualifying filters are biased and every rank in this scheme inherits the bias.
- **Same-variant vs distinct-variant excess is emitted but not acted on.** Step 9 computes
  `n_sites`, `n_trios_obs`, `recurrence_shape` and `max_site_share` per gene, but the
  **de-escalation rule is not wired**: a gene whose excess is really *one recurrent site carried by
  many trios* (a founder allele, or a recurrent mapping artifact at a single site) should have that
  **site** penalised and routed to read-level review, not the whole gene. Until it is wired, **every
  `excess_ratio` is an upper bound on distributed excess.** The two shapes have opposite
  interpretations and `E_g` is a count of mutational *opportunity*, so a per-gene `n` inflated by
  200 trios carrying one founder allele is compared against a target that never anticipated
  repeated sampling — the arithmetic is right and the interpretation is wrong.
- **Per-trio concentration is computed but not flagged.** A gene whose excess comes from 1–2 trios
  is a *sample* problem (contamination, audit **A-1**), not a *locus* problem. `n_trios_obs` is
  emitted; the `n_g ≥ 10 & n_trios ≤ 2` flag is not.
- **No NMD, no LOFTEE, no REVEL, no ClinVar stars, no parental genotype QC.** Four of these are
  resource questions with known fixes ([ROADMAP.md](ROADMAP.md) R1/R4/R5 and three extra VEP
  fields); the fifth needs columns Step 5 does not currently carry.
- **Segdup and family signals are optional inputs.** Without a segdup table matching the offset
  table's coordinate build, `sig_segdup` is silently off; Step 9 warns, but the corroboration count
  is then measured on five signals rather than six and the tier volumes will differ from the
  validated table.
- **The triaged gene LIST is coverage-confounded even though the triaged VOLUME is not.** 81% of
  triaged genes sit in the bottom two `exp_syn` deciles (a size × callability proxy) with a 42.7×
  rate gradient, yet 20 triaged genes in the *top five* deciles carry 44% of all triaged volume.
  Since the deliverable is variant triage the volume view is operative, but **the gene count must
  never be presented as N independently-suspicious loci.** This is the same finding as the
  low-count problem seen from the other side, and the `n_g ≥ 3` floor addresses both.
- **The established-gene union is now 2,218 genes and the retention figure survived expanding it**,
  which was Phase 1's largest stated risk. The residual weakness is provenance, not coverage: ACMG
  SF v3.3's 84 symbols are recalled rather than machine-read from a versioned file. It contributes
  55 genes at max 3.98× excess and triages zero, so nothing here depends on it — but replace it
  before publication.
- **The NB is ~3× conservative at α = 1e-3**, so it under-calls artifact genes. That is the correct
  direction under the never-drop rule, but it means the triaged count is a **floor, not an
  estimate**.
- **The offset table itself was not re-derived from source.** gnomAD's hosting bucket is
  denylisted in this sandbox, so `mu_*` was verified for internal consistency but not against
  gnomAD itself. Every `E_g` rests on it — the largest unchecked dependency in the whole scheme.
- **`n_trios_screened` is taken on trust from the resolved manifest.** `per_trio` and
  `sig_saturation` scale inversely with it, and saturation is the strongest signal in the panel, so
  FREEMIX-failing trios must be excluded **upstream** or the denominator and every count drift apart.
- **Homopolymer / low-complexity context and pseudogene-parent overlap are unimplemented.** Both
  need per-variant reference context or a pseudogene annotation; indel-heavy genes in
  low-complexity tracts are a distinct artifact class from segdups and are currently undetected.
- **Everything upstream still applies.** Step 9 re-ranks what Steps 3–5 produced; it cannot see
  CNV/SV, mosaic calls below the AB floor, or a comp-het that failed to pair. See
  [limitations.md](limitations.md).

---

## Sources

Retrieved and verified 2026-07-29; PMIDs and DOIs checked against PubMed metadata.

- Karczewski KJ, Francioli LC, Tiao G, et al. *The mutational constraint spectrum quantified from
  variation in 141,456 humans.* Nature 2020. PMID **32461654**, doi
  **10.1038/s41586-020-2308-7** — `mu_mis`/`mu_syn`/`mu_lof`, `oe_syn`, `constraint_flag`,
  `classic_caf`, `exp_syn`, gene coordinates (GRCh37).
  https://www.nature.com/articles/s41586-020-2308-7
- Tavtigian SV, Harrison SM, Boucher KM, Biesecker LG. *Fitting a naturally scaled point system to
  the ACMG/AMP variant classification guidelines.* Hum Mutat 2020. PMID **32720330**, doi
  **10.1002/humu.24088** — the additive-points structure this composite borrows (and the bands it
  deliberately does not use). https://onlinelibrary.wiley.com/doi/10.1002/humu.24088
- Richards S, Aziz N, Bale S, et al. *Standards and guidelines for the interpretation of sequence
  variants.* Genet Med 2015. PMID **25741868**, doi **10.1038/gim.2015.30**
- Abou Tayoun AN, Pesaran T, DiStefano MT, et al. *Recommendations for interpreting the loss of
  function PVS1 ACMG/AMP variant criterion.* Hum Mutat 2018. PMID **30192042**, doi
  **10.1002/humu.23626** — NMD-escape grading; why V5 is unreachable here.
- Singer-Berk M, Gudmundsson S, Baxter S, et al. *Advanced variant classification framework reduces
  the false positive rate of predicted loss-of-function variants in population sequencing data.*
  Am J Hum Genet 2023. PMID **37633279**, doi **10.1016/j.ajhg.2023.08.005**
- Walker LC, de la Hoya M, Wiggins GAR, et al. *Using the ACMG/AMP framework to capture evidence
  related to predicted and observed impact on splicing: Recommendations from the ClinGen SVI
  Splicing Subgroup.* Am J Hum Genet 2023. PMID **37352859**, doi
  **10.1016/j.ajhg.2023.06.002** — the 0.2 / 0.5 SpliceAI bands, and the low-score-canonical-site
  caveat.
- Pejaver V, Byrne AB, Feng B-J, et al. *Calibration of computational tools for missense variant
  pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria.* Am J Hum Genet
  2022. PMID **36413997**, doi **10.1016/j.ajhg.2022.10.013** — the 25.3 provenance and the REVEL
  bands.
- Stenton SL, Pejaver V, Bergquist T, et al. *Assessment of the evidence yield for the calibrated
  PP3/BP4 computational recommendations.* Genet Med 2024. PMID **39030733**, doi
  **10.1016/j.gim.2024.101213**
- Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. *Predicting Splicing from Primary
  Sequence with Deep Learning.* Cell 2019. PMID **30661751**, doi **10.1016/j.cell.2018.12.015**
- Guo MH, Plummer L, Chan Y-M, et al. *Burden Testing of Rare Variants Identified through Exome
  Sequencing via Publicly Available Control Data.* Am J Hum Genet 2018. PMID **30269813**, doi
  **10.1016/j.ajhg.2018.08.016** — TRAPD, the ancestry-matched case-vs-gnomAD upgrade path this
  statistic deliberately is *not*.
- Zhang C, Chen H, Deng Z, et al. *Heterozygous Germline Fanconi Anemia-Related Gene Mutations
  Increase Susceptibility to Germ Cell Tumors.* JCO Precis Oncol 2025;9:e2500435.
  PMID **40906985**, doi **10.1200/PO-25-00435** — the het-carrier evidence behind the FA/HR
  overlay rows (five-gene combined OR 10.17, 95% CI 4.87–21.27, P = 2.90e-06; 23-gene pooled OR
  4.14, 95% CI 1.98–8.66, P = .0013). **Both the per-gene and the pathway-collapsed claims come
  from this one study**, which is why the two priors are combined by MAX rather than summed.
- Lee K, Abul-Husn NS, Amendola LM, et al. *ACMG SF v3.3 list for reporting of secondary findings in
  clinical exome and genome sequencing: a policy statement of the American College of Medical
  Genetics and Genomics (ACMG).* Genet Med 2025;27:101454. PMID **40568962**, doi
  **10.1016/j.gim.2025.101454** — 84 genes, in the positive-control union. Membership recalled
  rather than machine-read; see §6.
- Genomics England PanelApp, panel **243 "Childhood solid tumours" v5.12** (version_created
  2026-05-06) and panel **259 "Childhood solid tumours cancer susceptibility" v1.30** (2025-10-13),
  `confidence_level == 3` (green). https://panelapp.genomicsengland.co.uk/
- UCSC Genome Browser `genomicSuperDups` segmental-duplication track (hg19), 50,137 records.
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
- ClinGen gene-disease validity, dosage sensitivity, and clinical actionability curations —
  the phenotype-agnostic positive-control union. https://search.clinicalgenome.org/
- This repo: [README.md#canonical-defaults](README.md#canonical-defaults),
  [gene_burden.md](gene_burden.md) (the case-only recurrence null),
  [limitations.md](limitations.md), [gene_constraint.md](gene_constraint.md),
  `SCIENCE_AUDIT.md` findings **A-1** … **A-8**.
