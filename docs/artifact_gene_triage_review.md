# Artifact-gene triage review — why `genes.ranked.tsv` is topped by TTN/MUC16, and what to do

**Status:** diagnostic review + sequenced remediation plan. Nothing here is implemented yet.
**Subject:** Step 6 (`pipeline/06_gene_burden.py`) — its recurrence null, and the `genes.ranked.tsv`
ordering that null produces. This is a review of the **code**, not of any one run's output.
**Reading order:** §1 diagnosis → §7 the fix, sequenced → §10 limits.
If you read one section, read **§10 (limits)**: the number of exome-wide-significant genes a
non-jointly-genotyped trio screen of this design can support is approximately **zero**, and every
fix below is about producing a *defensible nomination list*, not about rescuing a significance claim.

The primary evidence here is **algebraic**. Each defect is derived from the code and stated as a
closed form that holds for any cohort; a measured table then follows as confirmation. Where a
diagnostic needs an input the reviewer must supply (`candidates.calls.tsv`, `audit/counts.tsv`,
`qc_report.tsv`), it is stated as a **REQUIRED CHECK** with a runnable command, never as a result.

---

## 0. Notation, constants, and the reference run

### `N` is a symbol

`N` is the number of screened trios. Every formula below is written in `N` and holds for any `N`.
Where a formula is instantiated to make it legible, the instantiation is an explicitly
**hypothetical** worked example at `N = 200` — chosen to match the worked example already in
`docs/gene_burden.md:203` — and is never a measurement.

### `N` is audited, but absent from the output table

`n_trios` **is** audited — `06_gene_burden.py:391`,
`audit.record("06_burden", "n_trios", n_trios)` — and echoed to stderr at `:398`. It appears
nowhere in the output header (`out_cols`, `:367-373`), so a reader holding only
`genes.ranked.tsv` cannot see the denominator of every p-value in it. **Fix**: emit `n_trios` as
a column (and audit `n_males`, per §9 — that one genuinely is missing).

The gap is cosmetic; what follows is not. For any gene whose qualifying alleles are all absent
from the frequency oracle and therefore all floored at `f = burden.absent_af_floor`,
`p_carrier_hwe` (`:89-98`, a product over alleles) collapses and the emitted expectation
(`:265`) is exactly

```
exp_carriers = N · (1 − (1 − f)^(2K))          K = number of distinct qualifying alleles
             ≈ 2fNK
```

— a discrete lattice in `K` with step `N(1−f)^{2K}(1−(1−f)²) ≈ 2fN`, **independent of the gene**.
That is the point: a null that depended on anything gene-specific — length, constraint,
mutational target — could not lay its expectations on a gene-independent lattice. The lattice is
also why `N` is recoverable from an emitted table to six significant figures from the modal
`exp_carriers` alone; that is a diagnostic, never a workflow.

### Config constants (values must match `docs/README.md#canonical-defaults` — golden rule 4)

| Constant | Value | Source |
|---|---|---|
| `burden.absent_af_floor` (`f`) | 1.0e-6 | `config/config.example.yaml:383` |
| `burden.min_carriers` | 2 | `config/config.example.yaml:374` |
| `burden.exome_wide_p` | 2.5e-6 | `config/config.example.yaml:388` |
| `burden.weight_by_constraint` | true | `config/config.example.yaml:384` |
| `filters.rarity.dominant_max` | 1.0e-4 | `config/config.example.yaml:282` |
| `filters.rarity.recessive_max` | 1.0e-2 | `config/config.example.yaml:283` |
| `filters.constraint_weighting.shet_min` | 0.10 | `config/config.example.yaml:348` |
| `resources.gnomad.oracle` | `faf95` (default) | `config/config.example.yaml:194` |

### The reference run (stated once; referred back to, never restated)

Some claims below are confirmed against **one validation run of this pipeline**, reported
anonymously. Its scale, given here once so later rates have a denominator:

> **Reference run.** `N` in the low hundreds of trios. `genes.ranked.tsv` carries **12,187** gene
> rows. **6,374** of them have ≥ `min_carriers` dominant carriers and therefore carry a dominant
> `p_recurrence` — "**tested genes**" below means these. The gnomAD v2.1.1 constraint /
> mutational-target universe used for every refit is **19,658** genes, joining the emitted table
> by symbol at a **90%** rate. Refits were recomputed in pure Python from the emitted
> `genes.ranked.tsv` and the real gnomAD v2.1.1 constraint table.

Every reference-run figure below **illustrates a general property derived from the code**. None of
them is the finding. If the derivation and the measurement ever disagree, the derivation wins and
the measurement is a bug report about that run.

---

## 1. Diagnosis

### 1.1 The null is conditioned on the data it tests

`pipeline/06_gene_burden.py:259-265` (`_recur`) tests the observed carrier count against
`Binomial(N, p)` where `p` comes from `p_carrier_hwe` (`:89-98`) evaluated over `g["dom_faf"]` —
a dict populated at `:188-193` **inside the same loop that increments the carrier set**. The
alleles in the product are, by construction, exactly the alleles carried by the individuals being
counted. Conditional on "these alleles were observed segregating here", `X ≥ n` holds with
probability 1. The quantity reported as `p_recurrence` is therefore the *prior probability of the
observed allelic configuration*, not a tail probability.

**General consequence: `o/e ≥ 1` is structural.** There is no configuration of the input data that
produces a gene performing *below* its own expectation, because the expectation is assembled from
the observations. A calibrated null puts roughly half the genome on each side of 1.

Confirmation, on the reference run, over all tested genes:

| Observed / expected dominant carriers | |
|---|---|
| min | **15.6** |
| p25 / median / p75 | 71.5 / **121.2** / 294.6 |
| max | 8,080.8 |
| fraction > 100× / > 1000× | 58.6% / 10.4% |
| pooled across tested genes (Σ observed / Σ expected) | **104.9×** |

The load-bearing cell is **min = 15.6**: not one gene in the exome under-performs its own null.

### 1.2 The floor makes the excess a gene-independent constant

Take a gene with `K = n` distinct qualifying alleles, all absent from the oracle and therefore all
floored at `f`. Then `p_null = 1 − (1−f)^{2K} ≈ 2fn` and `E = 2fNn`. Both the observation and the
expectation are linear in `n`, so

> **o/e = 1/(2·f·N) — identical for every gene, in closed form.**

For a hypothetical cohort of `N = 200` trios at the default `f = 1e-6`, that is **2,500×**. The gene
contributes nothing: not its length, not its constraint, not its mutational target.

This is exact, and it is *not* an "inflation factor". A ratio of rejection counts has no calibrated
meaning when the inputs are not p-values, so phrasing this as "*x*× tail inflation" (as an earlier
draft did) is a category error. The honest empirical statement is the spread in §1.1 — real genes
mix floored and measured alleles, so the observed ratios scatter around the closed form rather than
sitting on it.

**Consequence: `min_carriers` and the significance threshold coincide by construction.** Since
`P(X ≥ n) ≈ C(N,n)·(2fn)^n ≈ (2feN)^n`, the reported *p* is an almost exactly geometric function of
`n` — `log₁₀p / n` is a constant that depends on `f` and `N` and on nothing else. Instantiated at
the hypothetical `N = 200`:

| n (all-floored, N = 200) | `p_null` | reported *p* | log₁₀p / n |
|---|---|---|---|
| 2 | 4.00e-6 | **3.182e-07** | −3.249 |
| 3 | 6.00e-6 | 2.834e-10 | −3.183 |
| 4 | 8.00e-6 | 2.646e-13 | −3.144 |
| 5 | 1.00e-5 | 2.532e-16 | −3.119 |
| 8 | 1.60e-5 | 2.360e-25 | −3.078 |
| asymptote log₁₀(2·f·e·N) | | | −2.9636 |

`n = 2` — the definition of "recurrent" — already clears the 2.5e-6 line by an order of magnitude,
and it does so for every `N` up to ≈ **560** at the default `f`. Exactly: the all-floored `K = 2`
case has `p_null = 1 − (1−f)⁴`, and the binomial `P(X ≥ 2)` is 2.492e-6 at `N` = 559 and 2.501e-6
at `N` = 560 — the crossover. (To leading order `P(X ≥ 2) ≈ C(N,2)·(4f)²`, which crosses at the
same integer; the cruder asymptotic `(2feN)² < 2.5e-6`, i.e. `2feN < (2.5e-6)^{1/2}` ≈ 1.58e-3,
gives the looser bound `N` < 291.) Below that bound the threshold is not doing any work;
above it the `n = 2` rung stops clearing the line unaided — which is a fact about `N`, not
evidence that the null is calibrated. On the reference run, **507 genes were declared exome-wide
significant on exactly two carriers**.

| Reference-run output totals | |
|---|---|
| `recurrent` (≥ 2 carriers) | 6,565 / 12,187 (53.9%) |
| `recurrence_exome_wide_sig` — dominant / biallelic / X-linked / **union** | 2,992 / 35 / 0 / **2,996** |
| max `p_recurrence` anywhere in the table | **7.566e-03** (i.e. *every* tested gene is p<0.05) |
| median `p_recurrence` | 3.623e-06 |

### 1.3 The ranking is therefore a gene-length ranking

The chain is: *p* is monotone in `n` (§1.2); `n` is the number of individuals carrying **any**
qualifying rare allele in the gene; that count scales with the gene's mutational target, i.e. with
coding length. Nothing in the sort key opposes it — `rank_key` (`:359-364`) sorts on `best_p` before
constraint. So the emitted ordering is a length ordering with noise, for any cohort.

Confirmation, on the reference run:

| Spearman ρ | value |
|---|---|
| rank vs −`n_carriers` | **+0.9174** |
| `n_carriers` vs CDS length | +0.4903 |
| rank vs CDS length | −0.4776 |
| **`exp_carriers` vs CDS length** | **+0.2340** ← the offset barely tracks the target it should absorb |

That last row is the diagnosis in one number: a correctly specified expectation would track the
mutational target nearly perfectly, because the target is what generates the count.

Median CDS length by rank band, reference run — the same length ordering, seen directly:

| Median CDS length (gnomAD v2.1.1) | bp | vs genome-wide |
|---|---|---|
| all gnomAD genes | 1,263 | 1.0× |
| top 500 of that run's order | 5,132 | 4.1× |
| top 100 | 9,489 | 7.5× |
| top 20 | 14,100 | **11.2×** |

FLAGS contamination, measured on the reference run's rank order against the published Shyr-2014
100-gene list (99 of the 100 are present in that run's table):

| K | FLAGS@K | PanelApp-green peds-cancer genes in top K |
|---|---|---|
| 20 | 0.600 | 0 |
| 50 | 0.420 | 1 |
| **100** | **0.390** | 2 |
| 500 | 0.156 | 6 |

*(Earlier drafts quoted 51/100 and 52/100. Those were measured on raw `n_carriers` or on a
hand-assembled family regex. The single headline number is **39/100 on the published FLAGS-100
list, measured on the reference run's rank order**.)*

---

## 2. Prior art: this exact failure was solved in 2013 (MutSigCV)

**Lawrence MS, Stojanov P, Polak P, et al. Mutational heterogeneity in cancer and the search for
new cancer-associated genes. Nature 2013;499(7457):214-218. PMID 23770567.
doi:10.1038/nature12213.**

This is the canonical solution and the dossier that preceded this document cited it nowhere.
The somatic-recurrence field hit precisely this wall: naive recurrence testing produced
implausibly long significant-gene lists topped by **TTN, MUC16 and the olfactory receptors**, and
the fix was to model the per-gene background rate with covariates — **gene length, DNA replication
timing, and expression level** — plus a **per-sample rate term**. The false-positive list
collapsed. The parallel here is exact down to the gene names, because the failure is a property of
the estimator, not of the assay.

Three design inputs follow directly:

| MutSigCV term | hprv analogue | Status |
|---|---|---|
| per-gene mutational target | `C·mu_g` from gnomAD v2.1.1 `mu_mis`+`mu_syn`+`mu_lof` | Available; §7 item 1 wires it |
| **per-sample rate term** | per-trio candidate yield from `audit/counts.tsv` | **Missing** — §5.2 is the required check |
| replication timing, expression | no covariate in the gnomAD constraint table | **Gap, not substitutable** — state it, do not proxy it |

**Resolving the Step-9 tension.** `prioritization.excess.offset.covariate_adjust` defaults
**false** (`config/config.example.yaml:486`; the same key in `docs/README.md:271`) — note the
`offset:` level in the path, without which `get()` silently returns the default and the switch
does nothing. `docs/prioritization.md` records a fitted covariate NB-GLM
(`+log(gene_length) + oe_syn`, ΔAIC = −121.1) that was correctly rejected *there* for two reasons
that **do not transfer to Step 6**:

- "Spearman 0.992, negligible where it matters" was measured on Step 9's *artifact-ranking* task,
  not on a gene-discovery ranking.
- "`oe_syn` is itself one of the six artifact signals" is an argument about not absorbing a signal
  Step 9 *reports separately*. Step 6 reports no such witness.

So the same covariate model that is right to reject in Step 9 is the right **offset** for Step 6.
The two steps consume the offset for opposite purposes.

---

## 3. The two inheritance arms scale differently — one offset is mis-specified

`06_gene_burden.py:246` headlines `n_carriers = len(dom | bi | x)`. That union is not a quantity any
single expectation can be written for, and the reason is dimensional:

```
E_dominant,g  = C_d · mu_g                    (LINEAR in the mutational target)
E_biallelic,g ∝ (C · mu_g)²                   (QUADRATIC: Q_g = Σ qualifying q, capped at 1;
                                               dominant = 1−(1−Q)², biallelic = Q²)
```

`p_biallelic_hwe` (`:101-109`) squares the cumulative allele frequency; `p_carrier_hwe` (`:89-98`)
does not. A gene with twice the target therefore contributes twice the dominant carriers and **four
times** the biallelic ones. Summing the two arms into one rank statistic makes the ranking's
length-sensitivity a function of the mixture, which varies by gene.

Compounding it: the two arms are gated at **different rarity thresholds** — `dominant_max` = 1e-4
against `recessive_max` = 1e-2 — so `n_carriers` sums observations taken at frequencies two orders
of magnitude apart.

Confirmation, on the reference run, against gnomAD `cds_length` over the joinable genes:

| CDS decile | median CDS (bp) | Σ n_dominant | Σ n_biallelic | bi/dom |
|---|---|---|---|---|
| 1 | 537 | 1,417 | 20 | 0.0141 |
| 5 | 1,494 | 1,940 | 68 | 0.0351 |
| 9 | 3,297 | 3,136 | 169 | 0.0539 |
| **10** | **5,466** | 5,139 | 565 | **0.1099** |

A 7.8× rise across the deciles, steeper than the pure-quadratic prediction — so a second mechanism
is present on top of the algebra. Concentration on that run: the biallelic arm is **4.72%** of all
carriers genome-wide but **15.0%** of the top 100 and **25.6%** of the top 20, with TTN supplying 40
of 79 top-20 biallelic carriers.

The two extra mechanisms are *different kinds* of defect and neither is fixed by an offset:

1. **Frequency.** `05_inheritance_screen.py:315-317` collects hets at `recessive_max` = **1e-2**,
   and the pairing at `:332-333` is an **uncapped cross product** (`m·p + m·d + p·d` rows over the
   mat / pat / denovo leg sets — quadratic in the gene's het count, with no cap). Two ≤1%
   missense in trans in a 35,991-residue protein is the expected state of a random genome, so the
   pair count grows quadratically in exactly the genes with the most hets.
2. **Phase.** Trans evidence for every pair is `clean_parent`, whose AD limb **fails open** — a
   GATK ref-block `0/0` parent carries `GT:DP:GQ:MIN_DP:PL` with no `FORMAT/AD` and passes
   vacuously, emitting `trans_evidence_unmeasured`. The number of opportunities for that vacuous
   pass is quadratic in the gene's het count, i.e. concentrated in the same long genes.

**Design consequence:** separate offsets per mode, and stop ranking on a union count.

**REQUIRED CHECK (needs `candidates.calls.tsv`):** report the fraction of `compound_het` rows
carrying `trans_evidence_unmeasured` / `origin_unverified`, overall and within the top 20. If it
is materially higher in the long genes, phase is not evidence there and the biallelic arm for
those genes must be reported as unconfirmed rather than counted. `whatshap` is already in the
image, so read-backed phasing is a wiring change, not a new dependency.

---

## 4. What a mutational-target offset actually buys — and its mirror-image artifact

The offset replaces the conditioned null of §1.1 with an expectation that does **not** see the
cohort's alleles: `E_g = C · mu_g`, with `C` fit over the **full gene universe including
zero-count genes**, and a negative-binomial tail for the overdispersion. The repo already contains
every piece (`src/hprv/prioritize.py`: `mutational_target()` :322, `fit_excess_null()` :426,
`nb_sf()` :242, `bh_fdr()` :285), and it must be applied **per mode**, on `n_dominant`, never on
`n_carriers` (§3).

The zero-count genes are not a technicality: the candidate list is a **zero-truncated** sample, and
on the reference run **45%** of the universe genes had no dominant carrier at all. Fitting `C` on
matched genes only would rescale it by `n_universe / n_matched` ≈ 1.8× and divide every excess
ratio by the same factor, which is the direction that *hides* artifact loci.

Refit on the reference run over that full universe:

**Fit: C = 4.5663e+04, α = 0.2070, 23 genes trimmed.**

`α` is a dimensionless dispersion and carries over to any cohort. **`C` is a run-scale constant
proportional to `N`** (`E_g = C·mu_g` against the public gnomAD mu table, so every `E` below is
too), quoted only so the refit is reproducible — never pair it with an exact per-proband or
per-gene rate, because the pair inverts to `N`.

| gene | n_dominant | E (= C·mu_g, reference run) | ratio | NB *p* |
|---|---|---|---|---|
| **TTN** | 39 | 65.30 | **0.60** | 0.81 — **DEPLETED** |
| OBSCN | 13 | 25.45 | 0.51 | 0.86 — depleted |
| MUC16 | 29 | 23.60 | 1.23 | 0.29 |
| RYR1 | 15 | 15.01 | 1.00 | 0.46 |
| NEB | 17 | 14.03 | 1.21 | 0.32 |
| BRCA1 | 7 | 3.06 | 2.29 | 0.078 |
| PALB2 | 5 | 2.02 | 2.47 | 0.085 |
| ATM | 5 | 5.16 | 0.97 | 0.52 |
| TP53 | 1 | 1.01 | 0.99 | 0.60 |

**TTN carries 60% of what its mutational target predicts.** That single row is the whole diagnosis:
the gene at the top of the uncorrected ordering is not merely unremarkable under a correct null, it
is *below* expectation. The same holds for OBSCN. These are statements about gene size, not about any
cohort — a gene whose target is 65 expected carriers cannot be surprising at 39 in any study.

### The catch: the offset swaps one artifact family for another

On the reference run, genes clearing *p* < 2.5e-6 under the corrected null: **6** (was 2,992).
Every one is a segdup / paralogue-collapse / VNTR locus. **Zero** are FLAGS-100. **Zero** are
cancer genes.

| gene | n | E | ratio | *p* | BH *q* |
|---|---|---|---|---|---|
| PKD1L2 | 12 | 0.48 | 25.2× | 2.80e-10 | 5.5e-06 |
| HEATR5A | 8 | 0.18 | 43.4× | 1.21e-09 | 1.2e-05 |
| TNXB | 16 | 1.13 | 14.1× | 5.09e-09 | 3.3e-05 |
| RAB44 | 10 | 0.57 | 17.6× | 9.45e-08 | 4.6e-04 |
| MUC19 | 19 | 2.17 | 8.75× | 3.83e-07 | 1.5e-03 |
| DNAH12 | 10 | 0.71 | 14.0× | 6.29e-07 | 2.1e-03 |

(BH q<0.05 over the full universe: 10 genes. The `excess.min_n_for_ratio_rule = 3` guard, config
`:461`, changes nothing here — all six have n ≥ 8.)

**FLAGS@100 = 0.000 under this ranking, and that figure is metric-gaming, not success.** FLAGS is
by construction the list of the longest, most variant-rich genes, so *any* ranking anti-correlated
with length zeroes it. The corrected head reads PKD1L2, HEATR5A, TNXB, RAB44, MUC19, DNAH12,
HMCN2, MYO15B, DNAH14, ZNRF2, NBPF12, FHAD1, SGCZ, NBPF11, STAM2 — a *different* artifact family,
and PanelApp-green genes in the top 100 = **0**. Report FLAGS@K as a **diagnostic**, never as a
success criterion, and always beside the peds-cancer count.

### Minimum detectable effect — the achievable excess is bounded by N/E

An excess ratio is `n/E`, and `n ≤ N` always. So the largest ratio a cohort can *ever* report for a
gene is `N/E`, and since `E ∝ mu_g ∝ length`, **the ceiling on detectable excess is inversely
proportional to gene length**. That is the mirror image of §1.3, and it is structural, not a
tuning problem.

Instantiated at the reference run's fitted `α` = 0.2070. The `E` column is an **illustrative
scale**, on round values; the gene names index target size only, and the exact per-gene `E` from
the refit is in the table above (BRCA1 3.06, NEB 14.03, TTN 65.30 — same conclusions, min *n* 22 /
73 / 307):

| E (illustrative scale) | min n for NB *p* < 2.5e-6 | implied ratio |
|---|---|---|
| 0.05 | 4 | 80× |
| 0.10 | 5 | 50× |
| 0.20 | 6 | 30× |
| 1.00 | 11 | 11× |
| 3.36 (BRCA1-scale) | 23 | 6.8× |
| 5.56 (APC-scale) | 34 | 6.1× |
| 15.4 (NEB-scale) | 79 | 5.1× |
| **70 (TTN-scale)** | **328** | 4.7× — exceeds `N` for any cohort in the low hundreds |

The bottom row is the one that generalizes: the required `n` grows with `E`, and `n ≤ N` always, so
for any gene whose minimum-detectable count exceeds `N` no cohort of that size can produce a
significant result at all. Conversely the top rows show the head of a corrected ranking filling with the smallest-`E`
genes. A flat `n ≥ 3` count floor does not repair this (n=3 at E=0.05 is 60× and passes). Two
mitigations, both cheap:

- Make the count floor **E-aware** (suppress the ratio rule below some `E_min`), and/or
- Rank on the **lower confidence bound** of n/E — the LOEUF construction — which shrinks small-E
  genes toward 1 automatically. See §7 item 4.

---

## 5. Cohort and batch QC — three checks that must run before any statistical conclusion

Step 6 consumes none of these signals, so none of them can be assumed to have been performed. All
three are cheap and two need no new code.

### 5.1 Contamination — the gate exists and is advisory-only *by design*

`pipeline/00_qc.py` computes per-member verifyBamID FREEMIX or a VCF-only CHARR proxy and writes
`kid_contam`/`dad_contam`/`mom_contam`/`contam_flag`/`overall_pass`. Its own comment
(`00_qc.py:227-230`) states the consequence: *"ADVISORY only … Steps 1/2/4 run over the full
resolved manifest and Step 6 never reads it, so a flagged trio still contributes calls and
recurrence pending human review."*

Contamination is the textbook generator of the observed pattern: a contaminated sample emits
spurious hets **in proportion to the number of variant sites in a gene**, i.e. in proportion to
gene length, in a **subset** of trios — which masquerades perfectly as "rare functional variants
recurring across individuals in long genes." It is indistinguishable from the §1 defect on the
output alone, and the two compound.

**REQUIRED CHECK.** From `qc_report.tsv`: report `sum(contam_flag)`, `sum(1-overall_pass)`,
`contam_source` (freemix vs the weaker charr proxy), and the max per-member contamination. Then
correlate contamination against per-trio candidate yield:

```bash
awk -F'\t' '$2=="05_inheritance" && $4=="candidate_calls" {print $3"\t"$5}' audit/counts.tsv \
  | sort -k2,2nr
```

A positive Spearman between per-trio contamination and per-trio candidate count means
contamination is generating candidate volume. If `contam_source` is `charr` rather than `freemix`,
the proxy has limited sensitivity and the question is **unresolved, not resolved negative**.

### 5.2 Per-trio outliers and batch structure

`05_inheritance_screen.py:443` already writes one `candidate_calls` row per trio into
`audit/counts.tsv`, so the distribution is already on disk. Report min / median / p90 / max and the
top-decile share. **If a handful of trios carry a disproportionate share of the calls, the cohort is
not `N` independent observations** — the effective sample size is smaller than `N`, every binomial
denominator in Step 6 is wrong, and a gene "recurrent" across three high-yield trios is a batch
artifact. Follow-through: add a per-trio rate term to the null (MutSigCV's per-sample term is the
direct analogue), and cross-tabulate outlier trios against `qc_report.tsv`'s
`contam_flag`/`mie_rate`.

Order-of-magnitude context from the reference run: the dominant arm emitted of order **10²
qualifying genes per proband**. That is call volume, not an expectation — the counts are
overdispersed and length-heterogeneous, so a flat-rate Poisson at that mean is an arithmetic
*reference* only and must never be phrased as "fewer recurrent genes than chance".

### 5.3 Platform / capture / workflow heterogeneity, and WGS-vs-WES

CLAUDE.md already records that the upstream genotype-refinement workflow **may skip**
`VariantAnnotator PossibleDeNovo` — so trios are not guaranteed to be uniformly produced, and the
pipeline is explicitly written to tolerate that. Uniformity must therefore be *measured*, not
assumed.

**REQUIRED CHECK.** `bcftools view -h` over each trio VCF; tabulate `##source`,
`##GATKCommandLine` versions, `##reference`, and `##contig` set differences per trio.
Heterogeneity means the cohort must be stratified before any recurrence claim. Resolve **WGS vs
WES definitively** from the same headers: it decides whether an external control arm should use
gnomAD v4.1's 730,947 exomes or 76,215 genomes, and using exomes against WGS candidates reports
infinite enrichment for every intronic call (the SpliceAI and CADD≥25.3 keep-paths reach
non-coding variants an exome control cannot see).

---

## 6. Reference and callability

### 6.1 Alt / decoy / HLA contigs are NOT excluded

The **only** region filter in the entire pipeline is chrM:
`pipeline/01_make_cohort_sites.sh:104` sets `EXCLUDE_CONTIGS="${HPRV_EXCLUDE_CONTIGS:-chrM,chrMT,M,MT}"`,
applied at `:137`/`:143`. `chr*_alt`, `chr*_random`, `chrUn_*`, `chrEBV` and `HLA-*` all pass.

Two consequences, the second previously unnamed:

1. An alt-contig call is a near-duplicate of a primary-contig call in the *same* individual.
   Carrier sets are keyed on `trio_id` so the carrier is not double-counted — but the alt copy adds
   a second distinct variant key to `dom_faf`, incrementing `K` and **shrinking the already-degenerate
   null further** by the exact lattice step of §0.
2. The compound-het pairer (`05_inheritance_screen.py:332`) crosses `by["mat"]` against
   `by["pat"]`. A real het and its alt-contig twin — the same physical allele, assigned opposite
   parent-of-origin by genotyping noise on a duplicated locus — is emitted as a **phase-confirmed
   in-trans compound het**. That is a concrete generator of the biallelic inflation in §3.

**REQUIRED CHECK — one command kills or confirms the hypothesis:**

```bash
bcftools query -f '%CHROM\n' <cohort_union.vcf.gz> | sort -u \
  | grep -c '_alt$\|_random$\|^chrUn\|^chrEBV\|^HLA-'
```

Zero ⇒ hypothesis dead, record it as excluded. Non-zero ⇒ extend `EXCLUDE_CONTIGS` to the
non-primary set. **This exclusion does not violate never-drop**: an alt-contig call is a duplicate
representation of a primary-contig call, not an independent observation.
Worked example to test first: **TNXB** sits in the MHC — the region with the most alt haplotypes in
GRCh38 — and has a well-known pseudogene twin, TNXA. On the reference run it is high in both the
uncorrected and the offset-corrected orderings (§4), which is exactly the signature a duplicated
locus produces.

### 6.2 No LCR / segdup / homopolymer mask anywhere in the variant path

There is no RepeatMasker, simple-repeat, ENCODE-blacklist, segdup, HLA or paralogue mask at any
step. The only segdup machinery is Step 9's **gene-level** `segdup98_frac`
(`09_prioritize.py:386-387`, loaded by `_keyed_by_gene` at `:535`, resolved at `:558`) — optional,
**symbol-keyed**, and shipping with no prep function (§7 item 10). So on any run without a
hand-built `--segdup` table the signal is simply OFF and the corroboration count is measured on
five signals instead of six. That absence is warned about only in the narrow case where the
`--mutrate` table itself carries a `segdup98_frac` column (`:559-563`); with an ordinary
mutational-target table it is silent, which is the gap worth closing.

The coordinate-build trap is **separate**, and it lives in how such a table must be *produced*,
not in the join: because the join is on gene symbol, the callsets' GRCh38 build is irrelevant to
it. gnomAD v2.1.1 gene coordinates are GRCh37/hg19 (`config/config.example.yaml:513-514`;
`scripts/prepare_resources.sh:542-543` repeats the warning), so a segdup track overlapped against
them to build the per-gene fractions **must be hg19 too** — an hg38 track yields ~0 overlap
silently. An hg19 track is the correct input here, not the defect.

Recommended, as a **variant-level FLAG feeding the never-drop down-weight, never a filter**: the
GIAB v3.6 GRCh38 stratifications (`LowComplexity`, `SegmentalDuplications`, `OtherDifficult`
incl. `collapsed_duplication_FP_regions` and the false-duplication copies) or Heng Li's LCR-hs38.
Both are GRCh38-native and chr-prefixed.

> **Contract note (golden rule 6).** A region BED is *neither* a VEP plugin *nor* one of the two
> sanctioned `bcftools annotate` transfers. Routing it through VEP `--custom` produces a CSQ column
> lifted by `+split-vep` — the same *mechanism* as a plugin field, but `--custom` is not named in
> the rule. **This document does not assume the ambiguity away.** Either (a) amend golden rule 6 in
> the same commit to admit `--custom` as a third sanctioned mechanism with its own namespace, or
> (b) route it through an actual plugin. Either way it needs an `annotations.F` entry, a getter,
> and a `*_status` provenance column, because blank ≠ zero.

---

## 7. The fix, sequenced

Each item names the **real** file, line and config key. Costs assume the reader has the repo open.

### 1. Give Step 6 the mutational-target offset — ONE atomic change with three parts

This is *not* "pass `--mutrate`". It crashes on the real file, and the config key that
`run_pipeline.sh` reads for Step 6 is not the one that points at the right table.

| Part | Location | Change |
|---|---|---|
| (a) gzip reader | `pipeline/06_gene_burden.py:48-51` — `_open_keyed` uses a bare `open(path)` | Add the gzip sniff that `pipeline/09_prioritize.py:151 _open_text()` already has. **Reproduced**: pointing `--mutrate` at the real `.bgz` raises `UnicodeDecodeError: 'utf-8' codec can't decode byte 0x8b in position 1` |
| (b) config key | `pipeline/run_pipeline.sh:290` reads `resources.mutation_rate_table`; `:371` binds `${MUTATIONAL_TARGET}` from `prioritization.resources.mutational_target` (config `:891`) | Either repoint `:290` at `prioritization.resources.mutational_target` or add a distinct `--mutational-target` flag. The two tables are **deliberately different files** — `scripts/join_constraint.py` projects the constraint file down to `gene/oe_lof_upper/pli/s_het/phaplo` and drops `mu_*`/`oe_syn`/`classic_caf`/`cds_length` |
| (c) why `dn_exp` is blank | `resources/manifest.env:129` sets `MUTRATE_URL=""`, yet `scripts/prepare_resources.sh:536` exports `MUTRATE_TABLE=$MUTRATE_OUT` unconditionally, so `resources.mutation_rate_table` resolves to a path that is never created and `run_pipeline.sh:293`'s `[[ -e "$mut" ]]` guard silently skips it | No new download needed: `prepare_resources.sh:302 prep_mutational_target` already produces `constraint/mutational_target.by_gene.txt.bgz` by copying the already-fetched constraint file |

Reuse Step 9's functions rather than reimplementing. **The names are** `fit_excess_null`
(`src/hprv/prioritize.py:426`), `fit_poisson_null` (`:513`), `calibrate_null` (`:557`), `nb_sf`
(`:242`), `nb_midp` (`:257`), `betainc` (`:178`), `bh_fdr` (`:285`), `mutational_target` (`:322`).
*(An earlier draft named `fit_nb_null`. **No such function exists** — `grep -rn fit_nb_null` returns
zero hits repo-wide.)* `prioritize.py` imports only `math`, `re`, `typing` and `hprv.config`
(`:50-54`), so Step 6 can import it without breaking the bare-`python3` test contract — and doing so
lets Step 6 drop its `scipy` try/except at `06_gene_burden.py:37-40` entirely.

Two invariants the port must carry over from Step 9, both documented in CLAUDE.md and both
direction-of-error traps: fit `C` over the **full gene universe including zero-count genes** (fitting
on matched genes only rescales `C` by `n_universe/n_matched` and divides every excess ratio by the
same factor — the direction that *hides* artifact loci), and compute the NB tail with the
regularized incomplete beta rather than a `1 − cdf` complement sum (the artifact tail reaches ~1e-241;
a complement sum bottoms out at machine epsilon and collapses every extreme locus into one bin).

**Test guard:** `tests/test_pure.py:2165 test_prioritize_reads_bgzipped_tables` covers Step 9's
reader but not Step 6's. Add the sibling.

### 2. Use a per-mode offset, never a single `n_carriers` E

Emit `n_dominant`/`E_dominant` and `n_biallelic`/`E_biallelic` as separate ranked quantities and
stop emitting a single `n_carriers`-driven rank. The dominant arm is linear in the mutational
target and the biallelic arm quadratic (§3), so no single `E` is correct for their union — and the
union additionally mixes a 1e-4-gated arm with a 1e-2-gated one.

### 3. Fix the multiple-testing treatment — the conditional p-value is the cheap correct fix

Selection is on `n` (`06_gene_burden.py:260` refuses to test below `min_carriers`) and the statistic
is strictly monotone in `n`, so the retained p-values are **not valid unconditional tail
probabilities** — they are the tails of a truncated distribution. The valid quantity is

```
P(X ≥ n | X ≥ min_carriers) = sf(n-1, N, p) / sf(min_carriers-1, N, p)
```

At `n = min_carriers` this is **exactly 1.0 for every gene, at every N** — which is the correct
answer, and the whole point: "the gene has at least 2 carriers" is precisely the event that was
selected on, so it carries no information. Instantiated at the hypothetical `N = 200`, all-floored:
n=2 → **exactly 1.0**; n=3 → 3.96e-4; n=5 → 1.27e-10; n=8 → 4.64e-20.

Applied to the reference run's emitted p-vector:

| | dominant-family exome-wide sig |
|---|---|
| unconditional (current code) | **2,992** |
| conditional `P(X≥n \| X≥2)` | **821** |
| of the 3,040 genes with n=2 | 507 sig → **0** (p_cond = 1.0 by construction) |

Three further points, each correcting an earlier draft:

- **Padding is not a substitute for conditioning.** Padding the dominant p-vector out to a full
  ~19,500-gene universe changes the rejection count **not at all** — the retained p-values are so
  far below any BH line that the denominator is irrelevant. On the reference run only the biallelic
  family moved materially (154 → 98). The "3.1× anti-conservative denominator" claim is retracted.
  *(An earlier draft simultaneously recommended padding with p=1.0 and computing the real p for n=1
  genes — mutually exclusive, and the latter would make essentially every gene significant.)*
- **Dependence.** BH needs independence or PRDS. The three families share the same variants, genes
  tested in more than one family are common (134 on the reference run), and `best_p = min(...)` is
  then taken across them with no Šidák correction and no q column. Benjamini–Yekutieli is the
  dependence-safe alternative (penalty `Σ1/i` = the harmonic number, **≈ 9.3** at m ≈ 6,000).
- **Threshold.** `exome_wide_p` = 2.5e-6 is applied three times, so the honest Bonferroni line is
  ~8.3e-7 — **irrelevant until the p-values are valid.** A stricter threshold on invalid p's just
  yields a smaller number of equally meaningless genes.

### 4. Rank on a shrunken effect size, not on *p*

Rank on the **lower confidence bound** of n/E (mirror of the LOEUF construction) or an
empirical-Bayes posterior `(n+α)/(E+β)`, keeping *p*/*q* as reported columns. This is what
disarms the inverse-length ranking of §4 without a hand-tuned floor: a gene with E=0.05 and n=3 has
a point ratio of 60× and a lower bound near 1, while a gene with E=15 and n=45 keeps most of its 3×.

### 5. Fix the two dead constraint limbs

**s_het is blank on every emitted row.** `CONSTRAINT_SHET_URL` **is** pinned
(`resources/manifest.env:104`), so the file is fetched — the failure is purely the join key.
`scripts/join_constraint.py:77` passes gene-key aliases `("gene","hgnc","symbol","gene_symbol")`;
the GeneBayes header is `ensg  hgnc  chrom  obs_lof  exp_lof  prior_mean  post_mean  …` and its
`hgnc` column holds **HGNC IDs** (`HGNC:24141`), so `_col` binds the key to an ID column and the
join yields zero overlap. (The value alias chain already includes `post_mean`; only the key is
wrong.) Fix: key on `ensg` and join via the gnomAD table's own `gene_id` column, per CLAUDE.md's
rule that symbol reconciliation must route through versionless Ensembl IDs, never aliases. Note the
repo already has the loud guard for exactly this at `join_constraint.py:98-103` — check whether it
fired and was missed. Extend `tests/test_pure.py:503 test_join_constraint`, which currently
fixtures a symbol-keyed `gene\tpost_mean` file and so cannot catch it.

Consequence today: `filters.constraint_weighting.shet_min = 0.10` is **inert config**, and the
`constrained` flag (`06_gene_burden.py:292-295`) is a 4-way OR running on three live limbs. s_het is
also the limb that rescues **short** constrained genes — exactly where LOEUF is underpowered and
exactly where the peds-cancer genes live.

### 6. Make constraint order the list, with a third state for "unmeasured"

`06_gene_burden.py:359-364` returns `(recurrent, same_variant, prec, con_key, …)`. `prec` is a
continuous float that is almost always unique — on the reference run 490 of the top 500 p-values
are distinct and the lowest rank involved in any tie is **113** — so `con_key` never fires where it
matters and `burden.weight_by_constraint: true` is effectively dead at the head of the list. This is
general: `p` is a continuous function of a floating-point product, so ties are measure-zero, and any
sort key placed *after* it is unreachable by construction.

**Present the reorder as a cosmetic mitigation, not a fix.** Promoting `con_key` above `prec`
produces a two-block partition inside which the broken *p* still does all the ordering. Three
specific defects to avoid carrying over:

- The boolean is knife-edge: TTN's LOEUF is **0.354** against a 0.35 cutoff, and on the reference
  run 162 genes lie in [0.34, 0.36] — moving the threshold by 0.01 moves the #1 gene across the
  partition. Use a **continuous** constraint score.
- **Genes with no constraint value at all** (loeuf, pli and phaplo all blank — 1,197 on the
  reference run, including all of its ENSG-keyed rows) would be forced into the "unconstrained"
  block, scoring "nobody looked" identically to "measured as tolerant": the exact blank-vs-zero
  conflation CLAUDE.md forbids for `clinvar_stars` and NHF. Add a `constraint_status` column,
  mirroring `E_source`/`nhf_status`.
- Golden rule 4: `docs/README.md:243-245` specifies the order as "rank recurrent-first, **then by
  `p_recurrence`**, then **weighted by constraint**". That table is binding — edit it in the same
  commit. Add a direct `rank_key` test next to `tests/test_pure.py:536 test_burden_helpers`
  (`_load_gb()` at `:528` already gives bare-`python3` access, and the scipy import is inside a
  try/except at `:37-40`, so the test stays dependency-free).

### 7. Resolve the ENSG-keyed rows

`06_gene_burden.py:163` keys on `r.get("symbol") or r.get("gene")`, so a call with a blank
`vep_SYMBOL` falls through to the Ensembl ID. Those rows join **zero** symbol-keyed resources
(constraint, mutational target, segdup, established genes, gene-prior overlay) and are
`constrained=0` by construction — the direction that **promotes** them, since an unconstrained gene
sinks only in a tiebreak that never fires (item 6). On the reference run 157 such rows exist and one
of them reaches rank 7 on 14 dominant + 2 biallelic carriers. Resolve via the constraint table's
`gene_id` column (versionless Ensembl ID, never alias matching), **retain** unmapped rows with an
explicit marker (never-drop applies to key resolution too), and audit the unresolved count.

### 8. `same_variant` recurrence is structurally dead where it matters

`rec_kind = "same_variant"` fires iff the **first** inheritance family to reach `min_carriers`
(dominant, then biallelic, then X-linked — `06_gene_burden.py:309-312` breaks on the first) has
exactly **one** distinct qualifying variant. It measures distinct-variant *count* within that one
family, not carrier sharing. Because a gene reaches the head of the list by accumulating many
distinct alleles (§1.3), the flag is by construction anti-correlated with rank and can never mark
the artifact head: on the reference run all 42 such genes rank 6,524 or worse, max `n_carriers` =
4, and all 100 top-ranked genes are `distinct_variant`. Yet 12 of those 42 carry
`recurrence_exome_wide_sig = 1`: two trios sharing one novel allele produce
*p* ≈ `C(N,2)·(2f)²` — **8.0e-8 at the hypothetical `N` = 200** — i.e. the single most
artifact-suspicious configuration receives the best possible score. Replace the binary with a per-gene, per-family concentration statistic
(`max carriers on one variant / n_carriers`, not short-circuited across families), and withhold
p/q/sig for K=1 rather than only demoting in the sort key.

### 9. Things that already exist — wire, do not build

| Proposal | Already in repo |
|---|---|
| NB excess statistic | `src/hprv/prioritize.py:426/513/557/242/285/322` — only the Step-6 call site is missing |
| Artifact gene-family regex | `src/hprv/prioritize.py:88-92 DEFAULT_FAMILY_PATTERNS` (`^MUC\d+`, `^OR\d+…`, `^HLA-`, `^NBPF\d+`, `^GOLGA\d`, `^PRAMEF?\d*`, `^KRTAP`), config-overridable via `prioritization.signals.family.patterns` — only a *continuous* rank column is missing |
| LOFTEE | `scripts/prepare_resources.sh:269-283 prep_loftee` fetches all four data files (URLs at `resources/manifest.env:95-99`); plugin code baked at `/plugins`. Missing: `--plugin LoF,…` in `pipeline/02_annotate_sites.sh`, the three CSQ names in the split-vep `want` list, the `annotations.F` entries |

### 10. Resource work that is NOT cheap — corrected sizing

| Proposal | Reality |
|---|---|
| "Fold gnomAD's `lcr`/`segdup`/`only_het`/`monoallelic`/`inbreeding_coeff` into the existing slim — small" | **Wrong.** `resources/manifest.env:22-26` pins **only** the JOINT release; those five fields exist only in the exomes/genomes releases. Needs a second multi-hundred-GB chromosome-wise stream (`docs/resources.md:277-278,:314` record ~877 GB). **The cheap half is real**: `fail_interval_qc`, `outside_broad/ukb_capture/calling_region`, `not_called_in_exomes/genomes` and FILTER **are** joint fields and can be appended to `GNOMAD_KEEP_INFO` (`manifest.env:35`, currently five fields) with a re-run of the existing stream. The Step-2 0-match `die` guard still holds. |
| "Populate `--established-genes`/`--segdup`; costs one script" | `scripts/prepare_resources.sh` has **no prep function** for any Step-9 optional input — they are commented-out placeholders at `:544` (`ESTABLISHED_GENES`), `:545` (`SEGDUP_TABLE`), `:546` (`GENE_MOI_TABLE`), `:548` (`GENE_PRIOR_OVERLAY`), with no manifest URLs. Budget a new `prep_*` with pinned URLs plus a `verify_extra` entry. **Guard interaction**: `prioritization.gene_downweight.min_control_genes: 1000` (config `:653`) means an under-populated control union **HALTS** — a partial list is worse than none. |
| "Make Step 6 HALT when `--mutrate` is absent" | Contradicts the documented contract (`09_prioritize.py:30-31`; CLAUDE.md names exactly three HALT resources: `spliceai_required` config `:100`, `missense_predictors_required` `:127`, `oracle: faf95` `:194`). The discriminating principle is: HALT when silent absence changes a *reported quantity* without changing whether the run appears to succeed. A missing offset **does** change every `p_recurrence`, so the case is makeable — but it must be argued, and shipped with the preflight check, the config key **and** the `docs/README.md` canonical-defaults row in the same change. Otherwise keep the WARN. |
| BP4 missense floor in `selection.py` | Contradicts a deliberate documented decision (`src/hprv/selection.py:11-15`, `docs/limitations.md #7`) and needs maintainer sign-off. **CI does assert that the predictor keep-reasons never fire** — `tests/integration/assert_integration.py:93-94` checks `revel`/`alphamissense`/`mpc`/`loftee_hc` are absent from the cohort-wide keep-reason set, and it runs on every push (`.github/workflows/build-publish.yml:171`) — so any floor must replace that assertion in the same commit. Separately, `tests/test_pure.py:406` and `:428` assert a MODERATE variant returns `(True, "impact_moderate")`, `:428` explicitly with CADD 0.1 to pin that MODERATE never consults a downstream score. That test survives **only** if the floor is specified as keep-on-missing-predictor (the FakeVar at `:428` carries no REVEL/AlphaMissense), so the missing-predictor behaviour must be written into the proposal, not left implicit. |

### 11. Internal-recurrence signals stay inside golden rule 2

Golden rule 2 permits internal recurrence as an **artifact/blocklist signal** and forbids it as a
population frequency. Each of these is inside that line, and none builds a genotype matrix:

- a **per-site called-trio count** is an artifact signal, never an AF;
- a **transmitted/untransmitted** ratio (§10.1) is not a population frequency;
- **haplotype-sharing** counts distinguish founder IBD from recurrent miscalling.

The per-site count requires changing `pipeline/01_make_cohort_sites.sh:181` from
`bcftools concat -a -D` to count-before-dedup — a change to the cohort-union contract, and
CLAUDE.md's `$(producer | head)` SIGPIPE gotcha applies to anything added there.
**Keep the negative recommendation too:** SKAT-O / regenie / SAIGE-GENE+ need the forbidden
genotype matrix. Do not attempt them.

---

## 8. External controls: what TRAPD/CoCoRV/COBT actually require here

This is the principled destination, and it is the item most likely to be implemented badly.

**A pLoF-only CAF offset does NOT fix this.** gnomAD v2.1.1 publishes `classic_caf` — the obvious
shortcut — and it is pLoF-only, while the screen's qualifying class is dominated by rare missense.
For a long gene the omitted mass is most of the target, so the offset is too small by roughly the
missense-to-pLoF ratio and the gene stays at the top. Worked on TTN with its published
`classic_caf` = 6.62e-3: `p_dom` = 0.0132, so this null expects `0.0132·N` dominant TTN carriers —
of order **3** in a cohort of `N` = 200. The observed dominant count is far above that (39 on the
reference run, §4), and the binomial tail lands **tens of orders of magnitude** below the 2.5e-6
line; substituting Roberts 2015's 2% population TTNtv rate does not move it across.
**TTN stays at rank 1.** The control carrier probability has to reach the order of **10⁻¹** — call
it a 6–7× rise on `classic_caf`, which is what including the full rare-*missense* cumulative
frequency would buy — before TTN stops clearing the line at all, and that quantity is the one
nobody has computed and the hardest part to make filter-symmetric. Also: gnomAD **v4.1's
constraint file removed `classic_caf` entirely**, so it is a v2.1.1/GRCh37 artifact as well.

> An external-control null that looks principled and leaves the reported problem intact is the
> worst outcome available here. The CAF offset works **only** if Q(g) is defined over the full
> Step-3 qualifying class on the gnomAD side — i.e. VEP 115 + CADD + SpliceAI re-run over the
> gnomAD sites file. **That, not the download, is the cost.**

**An incorrectly specified TRAPD is WORSE than the current null**, because its bias is
gene-specific rather than a diagnosable constant. Three hard prerequisites, not optional
refinements:

1. **A per-trio callable-region mask** giving a per-gene effective `N_g`. A Fisher 2×2 needs case
   non-carriers = `N − carriers`, which assumes every trio was callable at every qualifying site.
   These trios are not jointly genotyped, absent ≠ hom-ref, `01_make_cohort_sites.sh` strips INFO
   (`bcftools annotate -x INFO` at `:141`/`:147`), and no callable-region mask or gVCF ref-block
   artifact is retained anywhere. The denominator is wrong in a coverage-correlated way — worst in
   exactly the segdup/VNTR genes at the top of the corrected list. **Building this mask from gVCF
   reference blocks is the real engineering cost and must be scoped explicitly.**
2. **A group-matched control rate under a distinct column name** (e.g. `burden_control_af_group`),
   never the run's rarity oracle. `faf95` and the grpmax proxy are deliberately conservative
   maxima over ancestry groups — correct for a per-variant rarity *gate*, wrong as a control rate.
   Using them understates the control carrier rate in every non-grpmax group and manufactures case
   excess **multiplicatively across alleles** — worst in genes with many alleles, i.e. TTN and the
   mucins again. This is an independent, second mechanism driving the artifact head, and it
   survives the mutational-target fix.
3. **Q(g) defined on the gnomAD side first**, with the cohort counted within it — the data flow
   reverses relative to today.

Adopt CoCoRV's sampled-true-null inflation estimate as a **hard gate on the output**, not a
diagnostic.

**On the synonymous negative control.** It is *worthless before* the offset fix and *essential
after* it, so it is not a "decisive fix" today. The inflation is `1/(2fN)` — gene-independent
**and consequence-class-independent**, because it arises from conditioning the expectation on the
observed alleles, not from the functional filter. A synonymous run will show the same inflation for
reasons unrelated to what it is meant to test. When it is run, four constraints bind:

- Specify the **ratio** form (ProxECAT: case f/s vs control f/s), never the per-gene synonymous
  **offset** form. For *de novo* tests `mu_syn` is on the same mutation-rate scale as `mu_lof` and
  is a clean calibrator; for **inherited** counts the quantity is standing variation = mutation ×
  selection × drift, so `CAF_syn / CAF_functional` encodes the gene's selection coefficient and an
  offset deflates the statistic in exactly the constrained genes the screen is looking for.
- Restrict it to the **dominant arm** until comp-het pairing is capped: the biallelic count is
  quadratic in per-gene variant count while the dominant count is linear (§3), so a synonymous
  control (3–4× more variants per gene) inflates the biallelic arm by a gene-size-dependent
  ~10–16×, not a constant.
- Define the qualifying set for both classes on the **same measured frequency quantity** (variants
  where gnomAD published a real faf95), or the class-dependent gate-pass rate re-enters:
  `frequency()` returns **0.0** where gnomAD has the allele but published no fafmax
  (`rarity_basis = zero_ci`) and **None** where gnomAD has no record at all
  (`rarity_basis = absent`); `_floored` (`06_gene_burden.py:85-86`) charges `f` for both, and
  gnomAD AC is systematically lower for functional than synonymous alleles of the same gene.
- Exclude synonymous variants with SpliceAI ≥ `spliceai_ds_min` or CADD ≥ `cadd_phred_supporting`
  — they are not valid controls. And at cohort sizes of a few hundred trios the per-gene 2×2 cells
  are 0–3, so the LRT chi-square approximation fails; use an exact test.

---

## 9. The oracle interaction — the code comment is backwards

`06_gene_burden.py:173-178` and `docs/gene_burden.md:201` argue that reading the grpmax point
estimate makes the recurrence test **conservative**. The mechanism is right (`binom.sf` is
monotone increasing in p) but the premise inverted when `faf95` became the default oracle
(`config:194`):

- Under `oracle: faf95`, an allele gnomAD **has** but published no fafmax for resolves to **0.0**
  (`rarity_basis = zero_ci`) — CLAUDE.md records this as ~80% of a chr22 sample.
- `_floored` (`:85-86`) then charges it `f` = **1e-6**, identical to a truly novel allele, where the
  grpmax proxy would have charged its measured ~1e-5–1e-4.

Since faf95 ≤ the point estimate always, switching the default oracle moved the recurrence null
**strictly anti-conservative** while the comment kept asserting conservatism. This is the mirror
image of the correct reasoning at the rarity *gate* (where faf95 retains more) — the direction of
benefit flips with how the number is used, and only the gate side was reasoned through.

Two calibrations on the magnitude, correcting an earlier draft: the quoted 1,070×–1.3e12× figures
are the ratio between **two mis-specified nulls**, not the recoverable error. The proposed remedy
(an AC-based Poisson 95% upper bound, ~3/AN at gnomAD v4.1 joint AN ≈ 1.6e6) is ~**1.9e-6** against
the current 1e-6 floor — a ~2× change. Under a correct CAF-based null the floor question is moot,
because unobserved alleles enter the expectation directly rather than through a per-allele
detection limit.

Also note `_floored` only floors `None`/0/negative, so a genuinely measured value **below** `f`
passes through unfloored — on the reference run the minimum *p* at `n_dominant=2` back-solves to
q = 5.6e-7, below the nominal floor. **The "floor" is not a floor**, and any argument that leans on
`f` as a lower bound on `q_v` is unsound.

**Latent, not active — the X-linked null.** `p_carrier_hwe(f, floor, 1)` (`:271-272`) is a
hemizygous-*male* probability passed to `binom.sf(n-1, n_trios, p)` with **all `N`** trios, and one
blended p covers two genotype models (hemizygous male, hom female) that differ by orders of
magnitude. The family is typically inert — on the reference run 1 gene tested, 0 significant — and
both errors happen to run conservative, so it is correct by luck. It will bite when the family
populates. Step 5 already reads inferred sex from Step 0's `qc_report.tsv`, so `N_male` is available
and simply not plumbed through.

**Footnote-level, listed only so it is not re-raised.** `p_carrier_hwe` uses a product and
`p_biallelic_hwe` a sum. They agree to 0.0004% at Q=1e-5, 0.045% at Q=1e-3, 0.45% at Q=1e-2, and
diverge by 7% only at Q≈0.19 (the extreme reached by TTN on the reference run). Harmonising them
changes nothing. The two **real** biallelic defects are in §3: the missing inbreeding coefficient F
(at q=1e-3, F=0.01 the true hom probability is 11× the modelled one, and rare-disease trio cohorts
are enriched for consanguinity) and the uncapped quadratic pairing.

---

## 10. Limits — what cannot be fixed without joint genotyping or matched controls

Three hard limits. Each is a conclusion, not a caveat.

### 10.1 No internal control arm exists, by design

Golden rule 2 forbids a cohort genotype matrix. Therefore: **no internal allele frequency, no
internal HWE / excess-heterozygosity test, no internal case-control arm**, and no way to compute a
case-side allele number without a per-sample callability mask derived from gVCF reference blocks
that the pipeline does not currently produce. What *is* available, because it is computed
**within** a trio and never across the cohort, is a **transmitted/untransmitted (TDT-style)
ratio** — a within-family contrast, not a population frequency, and therefore inside golden rule 2
(§7 item 11). It is not implemented. Any external-control burden test
(TRAPD/CoCoRV/COBT) is **blocked on building that mask** (§8) — that is the real engineering cost,
and it must be scoped, not buried.

### 10.2 At any realistic `N` there is no recurrence signal to find

From the MDE relation in §4: a gene with E=1 needs **11** carriers; BRCA1-scale E=3.36 needs **23**;
a TTN-scale E=70 needs **328**, which exceeds `N` for any cohort in the hundreds and is why no gene
of that target size can be called at this scale. Meanwhile
the achievable ratio is bounded by `N/E`, so the only genes that *can* clear the line are the ones
with the smallest targets — where a handful of calls is already 30–80× (§4). The best peds-cancer
excess ratios observed on the reference run were ~2.3–2.5× at n ≤ 7 (BRCA1 2.29, PALB2 2.47) —
statistically indistinguishable from 1, and that is the expected result, not a disappointing one.

> **The defensible number of exome-wide-significant genes from a trio screen of this design and
> scale is approximately ZERO.** The handful that survive a corrected null are segdup/paralogue
> loci (§4).

Operationally: **withhold `recurrence_exome_wide_sig`, `q_recurrence` and
`q_recurrence_biallelic` from the output entirely** rather than recomputing them. Ship counts, an
excess ratio with its confidence bound, and a per-gene power flag. A column named `q_recurrence`
will be read as an FDR no matter what the documentation says.

### 10.3 An excess statistic is never evidence of disease association

There is no phenotype contrast in this design. "Excess over population expectation" is an
**artifact statistic**. CLAUDE.md already states this for Step 9 (13 established predisposition
genes were extreme excess outliers on a validation run; CTSA at 148×). It applies verbatim
here. Never phrase an excess result as evidence for or against a gene being real.

### 10.4 The corollary nobody wants to hear

**No reordering of the available columns surfaces the pediatric-cancer panel at the top**, because
`recurrent != 1` is sort key #1 and a true dominant predisposition gene in `N` trios has exactly
**one** carrier — the whole point of a rare-disease gene is that it does not recur at these sample
sizes. On the reference run TP53 ranked 9,014 and DICER1 7,292, both at n=1: the biologically
**correct** answer under a recurrence framing, which is itself the proof that the framing is wrong
at this scale. (Both are *present*, contrary to an earlier draft. Genuinely absent from that run:
RET, RUNX1, SMARCB1, PTCH1, VHL.)

The deliverable this design can support is a **ranked nomination list with an explicit artifact
axis**, not a significance claim.

---

## 11. Diagnostic: constraint-first ordering, and the ceiling it hits

This section is a **measurement, not a recommendation** — the experiment that sets the ceiling for
**§7 item 6**. It re-sorts an emitted `genes.ranked.tsv` on nothing but columns already in it
(`loeuf`, `pli`, `phaplo`), using no cohort data whatsoever:

```
score = (1 − min(loeuf, 1)) + pli + phaplo      # tie-break on n_carriers, descending
```

Two conditions the measurement required, both of which carry into any real implementation:

1. `n_carriers ≥ 1` must be filtered first, or de-novo-only genes leak in (the emitted table
   carries rows with `n_carriers = 0`).
2. Rows with no constraint value must go to a **separate unrankable bucket**, not to the bottom —
   blank is "nobody looked", not "not constrained" (§7 item 6).

Measured on the reference run:

| K | FLAGS@K before | FLAGS@K after | peds-cancer in top K before | after |
|---|---|---|---|---|
| 20 | 0.600 | **0.000** | 0 | 1 |
| 50 | 0.420 | **0.020** | 1 | 1 |
| 100 | 0.390 | **0.030** | 2 | 2 |
| 500 | 0.156 | **0.018** | 6 | **12** |

Selected rank moves on that run: TTN 1 → 2,815; MUC16 2 → 7,559; MUC19 4 → 10,723; NEB 8 → 3,474;
RB1 2,612 → 249; DICER1 7,292 → 291; TP53 9,014 → 2,323; APC 196 → 385. The *directions* are
general — FLAGS genes are long and LOEUF-tolerant, so they sink; the classic haploinsufficient
tumour suppressors rise — while the exact ranks are not.

**What it proves, and it cuts both ways.** Constraint ordering strips the artifact head
completely — and does **not** surface biology. Peds-cancer genes in the top 100 stay at 2. BRCA1
moves 49 → **7,520** on the reference run, because its LOEUF is 0.915 and it genuinely is not
LOEUF-constrained:
tumour suppressors acting through LoF-plus-somatic-second-hit are invisible to population
haploinsufficiency constraint, so a constraint-led rank actively buries the actionable class this
kind of screen exists to find.

A constraint-only ranking is a **fixed permutation of gene space** — identical for every cohort,
invariant to whether a single variant was called. That is the ceiling: constraint belongs in the
sort as a *continuous term beside a corrected excess statistic* (§7 items 4 and 6), never as the
ordering principle. Reproducing this table needs no tooling beyond an emitted `genes.ranked.tsv`;
reproducing §4's excess columns needs the mutational-target offset wired per **§7 item 1**, which
is the change that should be made anyway.

---

## 12. Doc and code edits this review makes mandatory (golden rule 4)

The canonical-defaults table is the single source of truth, and it is currently **wrong about the
recurrence null in both the field it names and the direction of the bias**. These are not optional
follow-ups; a code change that leaves them is a rule violation on its own.

| File / line | Current text | Required edit |
|---|---|---|
| `docs/README.md:230-231` | "`q_v` is the per-variant **grpmax proxy** AF (Step 5's `grpmax_af` column), **not** `faf95`" | `06_gene_burden.py:185-187` reads **`rarity_af`** first and falls back to `grpmax_af` only when `rarity_oracle` is absent (legacy tables). Name `rarity_af`; state the run-level `rarity_oracle` |
| `docs/README.md:238-240` | "these p-values are **conservative**" | **Reverse it** (§9): under the shipped `oracle: faf95` default they are strictly anti-conservative |
| `docs/README.md:243-245` | rank "recurrent-first, then by `p_recurrence`, then weighted by constraint" | Update if §7 item 6 lands |
| `docs/gene_burden.md:199` | "The primary recurrence signal is **calibrated**" | Delete "calibrated" — `:203` in the same section already says it is "too small (anti-conservative)" |
| `docs/gene_burden.md:15`, `:282` | recurrence nulls / BH-FDR / exome-wide flag listed as **IMPLEMENTED**, unqualified | Hoist `:203`'s caveat into both status tables: IMPLEMENTED-BUT-UNCALIBRATED, or remove the columns |
| `docs/gene_burden.md:19-24`, `:287` | "there is **no `faf95`** under this contract"; CI bound "unrecoverable" | Stale — contradicts CLAUDE.md golden rule 2 (faf95 is the DEFAULT with a preflight HALT) **and** line `:172` of the same document. Delete |
| `docs/gene_burden.md:288` | "Neither computable: no LOFTEE, no REVEL/AlphaMissense" | `config/config.example.yaml:127` sets `missense_predictors_required: true` |
| `docs/gene_burden.md:144-148` | `:144` announces "a deliberately **three-rung** ladder (IMPACT, SpliceAI, CADD)" and the numbered list at `:146-148` then enumerates two, omitting SpliceAI | Add the SpliceAI rung to the enumeration (`:287`'s table already lists all three, in the right order) |
| `CLAUDE.md:593` | names the integration mock `mock_annotate.py` | The file is `tests/integration/mock_vep.py` (correct at `:393` and `:442`) |

### The integration suite hard-asserts the pathology

`tests/integration/assert_integration.py:199-202` asserts
`float(genes["GENED"]["p_recurrence"]) < 1e-4` and
`genes["GENED"]["recurrence_exome_wide_sig"] == "1"`. GENED has exactly **2 dominant carriers of
one shared variant** (`recurrence_kind == "same_variant"`, asserted at `:194`) — precisely the
"2 carriers of a private allele clears the exome-wide line" behaviour this review condemns, frozen
as a CI pass condition (`.github/workflows/build-publish.yml:171`).

**Any null change must replace both assertions in the same commit**: GENED must retain its counts
and `recurrent=1` while `recurrence_exome_wide_sig` becomes `0`. Precedent for mock-scale config
overrides is already documented in `tests/integration/make_mock_data.py`
(`min_control_genes` 1000→3, `max_downweight_fraction` 0.20→1.0).

### Testing plan per change

| Change | Test location |
|---|---|
| `rank_key` reorder, offset null, `_open_keyed` gzip branch | `tests/test_pure.py`, beside `test_recurrence_null_per_model` (`:485`) and `test_burden_helpers` (`:536`) |
| s_het join key | extend `test_join_constraint` (`:503`) with an ENSG/HGNC-keyed fixture |
| output shape, never-drop | `tests/integration/assert_integration.py` |
| any threshold change | the checks list in `test_prioritize_config_matches_canonical_defaults` (`:2258`) — this is what enforces golden rule 4 mechanically |

`tests/test_pure.py` must run on a **bare `python3`**. The 9 Step-9 CLI tests already SKIP without
pyyaml via `_load_p9()` (`:678`), and `_run_all` refuses to print "all passed" when any skip
occurred. New Step-6 tests must avoid a hard yaml/scipy dependency at import time.

---

## 13. Corrections carried into this document

Recorded so they are not re-derived. Where a critic was wrong, the corrected framing is kept and
the reason given in one line.

| Claim in earlier material | Correction |
|---|---|
| An assumed `N_trios`, carried through every derivation | **`N` must be read or recovered, never assumed.** A 0.5% error in `N` propagates ~1.5% into a reproduced *p*, which is enough to make a "reproduces the emitted value" claim false. An earlier draft of *this* document said Step 6 records `N` nowhere — **wrong**: it is audited at `06_gene_burden.py:391` and echoed at `:398`. The real gap is narrower: `n_trios` is not in `out_cols` (`:367-373`), so a reader holding only `genes.ranked.tsv` has to recover it from the `exp_carriers` lattice (§0) |
| An earlier draft's `exp_carriers = N·K·(1−(1−f)²)` | **Algebraically wrong**, and asserted as exact. `p_carrier_hwe` (`:89-98`) is a *product* over alleles, so the all-floored form is `N·(1−(1−f)^{2K})`; `K` cannot be factored out of the exponent. The two agree only at `K` = 1 and the ≈ 2fNK approximation is unaffected, but the lattice claim is stated exactly and so must be exact (§0) |
| "CI asserts the missense-predictor keep-reasons never fire" — an earlier draft of *this* document said no such assertion exists | **Re-corrected: it does exist.** `tests/integration/assert_integration.py:93-94` asserts `revel`/`alphamissense`/`mpc`/`loftee_hc` are absent from the cohort-wide keep-reason set, and it runs in CI. `docs/gene_burden.md:168` and CLAUDE.md were right; the draft's `grep` was scoped wrong. The *rest* of that item stands: `tests/test_pure.py:428` also pins MODERATE-never-consults-a-score (§7 item 10) |
| "A properly-specified gnomAD CAF null demotes TTN" | Only with the **full qualifying class**. The cheap `classic_caf` form leaves TTN tens of orders of magnitude below the line, at rank 1 (§8) |
| "TTN's true carrier probability is 0.18–0.38, so p = 0.52" | **Circular** — that is a cohort's own observed rate tested against itself, which returns ~0.5 for any gene. Deleted. Only an external published rate is admissible, and it still leaves TTN significant under a CAF offset |
| "Length-aware nulls cut top-100 FLAGS from 40% to 0%" | True and **meaningless on its own** — any 1/length ranking does that. The corrected top is a different artifact family, and peds-cancer@100 = 0 (§4) |
| "Adopt constraint-score-alone as the shipped rank" | **Rejected.** It contains zero cohort data — a fixed permutation of gene space — and it buries BRCA1-class tumour suppressors (49 → 7,520 on the reference run). Retained only as the §11 diagnostic that sets the ceiling for §7 item 6 |
| "*x* significant genes = 61,374× tail inflation" | Not an inflation factor: a ratio of rejection counts has no calibrated meaning when the inputs are not p-values, and the arithmetic was unreproducible. Replaced by the exact closed form `o/e = 1/(2fN)` (§1.2) |
| "NB refit shows the tail is 9e9× heavier" | Dispersion was estimated from the **marginal** count distribution, which is dominated by the gene-size heterogeneity the offset exists to explain. Deleted; use residual dispersion around a fitted offset, as `fit_excess_null` already does |
| "Padding the BH denominator is 3.1× anti-conservative" | Verified false: padding to a full-universe *m* leaves the dominant rejection count unchanged. Only the biallelic family moves. The correct fix is the **conditional p-value** (§7 item 3) |
| "Reuse `fit_nb_null`" | **No such function.** `fit_excess_null` / `fit_poisson_null` / `calibrate_null` (§7 item 1) |
| "Pass `--mutrate ${MUTATIONAL_TARGET}` — small change" | Crashes on the real `.bgz`, and `run_pipeline.sh:290` reads a **different** config key (§7 item 1) |
| "Fold gnomAD `lcr`/`segdup` into the slim — small" | Those fields are not in the JOINT release; needs a second ~877 GB stream (§7 item 10) |
| TP53 and DICER1 "absent" | Present, both at `n_carriers = 1` and deep in the list — the correct answer under a recurrence framing (§10.4) |
| FLAGS top-100 = 51 | **39/100** on the published Shyr-2014 list measured on the reference run's rank order (§1.3). 51 was measured on raw `n_carriers` |
| TTN LOEUF = 0.8371 / MUC16 = 0.1208 | Those are the **phaplo** column, read by splitting the file on whitespace rather than tabs (7 empty fields collapse and shift by 3). Real values: TTN loeuf 0.354 / pli 2.558e-96 / phaplo 0.8371. The constraint join is correct to 99.92%; **there is no column-shift bug** and every emitted row has exactly 28 fields |
| "Excess ratio counterfactual on `n_carriers`" | Incoherent — `n_carriers` unions a 1e-4-gated arm with a 1e-2-gated arm, and the two scale differently in the mutational target. **All counterfactuals here are per-mode on `n_dominant`** (§3, §4) |
| An earlier draft's TTN CAF figures — *p* = 6.0e-15 under Roberts 2015's 2% rate, and "TTN falls out at π ≈ 5% (*p* ≈ 2e-4)" | **Neither reproduces** under the binomial the same paragraph specifies, at any cohort size: a 2% control rate is *more* extreme than 6.0e-15, and π = 5% still leaves TTN many orders of magnitude below 2.5e-6. TTN only stops clearing the line once the control carrier probability reaches the order of 10⁻¹. Significands deleted; the conclusion the paragraph exists for — TTN survives a `classic_caf` offset at rank 1 — is unchanged (§8) |
| "Two trios sharing one novel allele give *p* ≈ 6e-8" | An `N`-dependent quantity quoted as a constant, and inconsistent with this document's own hypothetical. The closed form is `C(N,2)·(2f)²` = **8.0e-8 at `N` = 200** (§7 item 8) |
| "The `n = 2` rung clears 2.5e-6 at every realistic cohort size" | Overstated. The all-floored `K = 2` binomial clears the line up to `N` = **559** and fails from 560 (the asymptotic form gives the looser `N` < 291). Above that the rung stops clearing unaided — a fact about `N`, not evidence of calibration (§1.2) |
| "E = 70 (TTN-scale) is UNREACHABLE" | Reachable in principle: at α = 0.2070 the minimum `n` is **328**. The bare "unreachable" both asserted a bound on `N` implicitly and was false above `N` ≈ 328. Now stated as the count, with the comparison to `N` explicit (§4, §10.2) |
| §4's MDE `E` labels (3.36 / 15.4 / 70) contradicting the refit table's per-gene `E` (3.06 / 14.03 / 65.30) | Both auditors flagged it. **Adjudicated by relabelling, not recomputation**: the column is now an explicitly *illustrative* round-number scale, with the refit's exact per-gene values and their min-`n` (22 / 73 / 307) quoted beside it. Recomputing the rows would have implied a per-gene precision the round `E` values never carried |
| "`segdup98_frac` is configured against an hg19 track while the callsets are GRCh38, so it silently matches ~0" | **Backwards.** The Step-9 join is on gene SYMBOL (`09_prioritize.py:386-387`, `:535`, `:558`), so the callsets' build cannot break it; hg19 is the *correct* build for the track, because the gnomAD v2.1.1 coordinates it is overlapped against are hg19. As written it would have sent an implementer to rebuild in GRCh38 — the exact error CLAUDE.md warns about. The real defect is that the table has no prep function and its absence is near-silent (§6.2) |
| Benjamini–Yekutieli penalty "≈ 8.8 at m ≈ 6,000" | `Σ1/i` is the harmonic number: **9.28** at m = 6,000, 9.34 at m = 6,374. 8.8 corresponds to m ≈ 3,700 (§7 item 3) |
| `prioritization.excess.covariate_adjust` | The key is nested one level deeper, under `offset:` — `prioritization.excess.offset.covariate_adjust` (`config:486`, named the same way at `docs/README.md:271`). The short path is silently inert, because `get()` returns the default (§2) |
| "`rec_kind = same_variant` fires iff the gene has one distinct qualifying variant in the whole cohort" | Scoped too widely: `06_gene_burden.py:309-312` **breaks** on the first inheritance family to reach `min_carriers` and tests `len(fmap) <= 1` for that family alone. The short-circuit is itself part of the defect (§7 item 8) |

---

## 14. Prior art searched, and not searched

| Area | Status |
|---|---|
| MutSigCV / somatic significantly-mutated-gene false positives (Lawrence 2013, PMID 23770567) | **Was missing entirely.** Now §2, first-class |
| Germline burden vs public controls (TRAPD PMID 30269813; CoCoRV PMID 35545612; COBT) | Covered §8 |
| Constraint / mutational target (Samocha PMID 25086666; Karczewski PMID 32461654; Zeng PMID 38977852) | Covered §7 |
| FLAGS (Shyr 2014, PMID 25466818; erratum PMID 29187224 — the raw-count table exists **only** under the erratum URL) | Covered §1.3, §4 |
| Somatic SMG false-positive literature beyond MutSigCV (expression / replication-timing covariates) | **Gap. Close before the offset is designed** — gnomAD's constraint table supplies neither covariate, and no proxy should be substituted |
| Index hopping / barcode swapping on patterned flow cells | Considered, **deferred** — pursue only if §5.1 contamination returns positive |
| Caller-specific artifacts (HaplotypeCaller vs DRAGEN) | Considered, **deferred** — pursue only if §5.3 workflow heterogeneity returns positive |

Recording the deferred items matches how `docs/limitations.md` already handles the ENCODE
blacklist and LOFTEE, and prevents them being re-proposed. (The ENCODE blacklist stays excluded:
it is derived from ChIP-seq input signal anomalies, its authors state it is unsuitable for WGS, and
it is almost entirely intergenic — it would not touch a single gene in the artifact head.)
