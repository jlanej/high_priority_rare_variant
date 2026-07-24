# Pipeline review — publication-readiness pass (2026-07)

A point-in-time review of `high_priority_rare_variant` aimed at publication quality: what is
**incorrect, mis-calibrated, silently lossy, non-reproducible, or over-claimed**, with enough
detail to assess each item and act on it. This is a findings log, **not** part of the canonical
methods reference; where it disagrees with the code, the code is the ground truth to fix.

## Resolution status — 2026-07 follow-up

Every finding below was **independently re-verified against the code** (adversarial code readers +
numeric reproduction for P1/P5), then addressed on branch `fix/pipeline-review-followups`
(commit `cd2e43e` — code/reporting fixes; `c1edf29` — honesty reframes). Eight findings are fixed;
two are honest-reframed with their heavier code fix **deferred** as a policy/calibration decision.

| # | Status | What changed (or why deferred) |
|---|--------|--------------------------------|
| **P1** | doc reframe **done**; recalibration **deferred** | Renamed to "raw (uncorrected) reference-read fraction" (CHARR-*like*, not the Lu-2023 estimator); corrected the "catches 1–3%" claim (flags only gross ≳5–8%); ROADMAP "✅ DONE" → "◐ PARTIAL". Full AF/baseline recalibration + threshold re-derivation deferred — alters which trios flag, needs spike-in calibration. |
| **P2** | **fixed** | Ingest build check now `die`s on an explicitly-declared non-GRCh38 assembly (scoped to the `##VEP=` line); an absent `assembly=` token still warns. |
| **P3** | doc reframe **done**; auto-gating **deferred** | `overall_pass` documented as advisory (surfaced to the reviewer, not auto-excluding); `pipeline_design.md` "which trios enter analysis" over-claim corrected. Automated exclusion deferred — never-drop ethos + P1's limited sensitivity make it a config-gated policy change. |
| **P4** | **fixed** | Step-5 audit counter `clinvar_plp_dropped_ge_recessive_max` + `clinical_classification.md` caveat. (Propagating P/LP into Step-5 rarity to gate at `ba1` remains a deliberate policy change — not done.) |
| **P5** | **fixed** | `dn_q_enrich` now BH-corrected over the full mutation-model universe (p=1.0 nulls padded in), not called-genes-only. `dn_exome_wide_sig` (fixed Bonferroni) unchanged. |
| **P6** | doc **fixed** | Documented the single pooled protein-altering ("prot") Poisson test; dropped the per-class advertisement. Per-class code split not done (off-by-default secondary arm). |
| **P7** | **fixed** | Recurrence sig/FDR audit counters + the xlsx About count now OR all three inherited families (`genes.ranked.tsv` was already complete). |
| **P8** | **fixed** | xlsx numeric cells coerced to typed numbers, so the auto-filter sorts numerically instead of lexicographically. |
| **P9** | **fixed** | gnomAD label derived from the `##VEP=` header (no more hardcoded "v4.1"); a Provenance block (gnomAD/VEP/cache/CADD versions + git short-SHA + verbatim `##VEP`) added to the About sheet. |
| **P10** | **fixed** | Step-1 per-trio `norm -c w` → `-c e` (fail fast, attributed to the trio) + explicit `-c e` at the union. |

**Trimmed as over-claims (verified, deliberately not acted on):** P1's framing of depth-weighting as
a co-equal defect (reproduced as near-inert at normal depth — the missing AF/baseline correction is
the real issue); P10's proposed `ref_mismatch` audit counter (unreachable under `-c e`); P5's "just
pass `m = len(mut)`" alternative (`bh_fdr` mis-ranks that way — padding the p-vector is required, and
is what was implemented).

**Deferred — genuine policy changes needing a decision / calibration data (NOT shipped):**
(1) **P1** full CHARR recalibration — move the gate post-annotation, add `mean(ref_AB)−baseline)/mean(1−AF)`,
and re-derive `qc.charr_threshold` from clean-vs-spike-in trios; keep FREEMIX as the production
default meanwhile. (2) **P3** automated exclusion of QC-flagged trios from the recurrence tally —
recommended to stay advisory, or become an **off-by-default** config knob; a safe non-policy middle
option is a non-dropping `qc_flagged_carriers` per-gene column in Step 6.

*(The sections below are the original review, preserved verbatim as the findings log.)*

## How this was produced

Nine expert lenses (statistics ×2, Step-0 QC science, normalization integrity, annotation
integrity, Step-3 selection, reporting/IGV, reproducibility, claims-vs-computation) read the
relevant code in full. **Every candidate finding was then adversarially verified twice** — once by
a code reader trying to refute it against the source, once by a statistics/genetics/reproducibility
referee judging whether it is a real defect or acceptable practice for a discovery *screen*. Only
findings that survived both verifications are listed under "Findings"; refuted candidates are in
"Considered and dismissed" for completeness. Every **major** was re-checked by hand against the
source before inclusion. No code was executed — triggers are code-path arguments, so a spike-in /
mock confirmation is the natural next step for the load-bearing ones (P1, P5).

## Scope — what is deliberately NOT here

Two earlier review passes already fixed or documented a large set of issues; those are **excluded**
and should not be re-raised:

- **Inheritance-model fixes (Step 5):** chrY/chrX conflation (`male_x_chrx`), HOM_ALT-parent
  deterministic origin, compound-het trans now requires an affirmative non-transmitting-parent
  hom-ref (`origin_unverified`), an unphased de-novo pair no longer vetoes the dominant call,
  `strict_gt=True`. All merged.
- **prepare_resources** `$?` bug + CADD/VEP-cache integrity check; **doc/tooling concordance**
  (base image `ensemblorg/ensembl-vep`, bcftools 1.23, conda-lock as TARGET); **WhatsHap** marked
  not-wired.
- **Documented residuals** in `limitations.md`: comp-het keyed on the single PICK'd gene; a male
  chrX hemizygote dropped when the mother is a no-call; the 1e-2/1e-4 permissive-partner asymmetry;
  no Y-linked model.
- **VEP-only contract losses** (no faf95/nhomalt/SpliceAI/LOFTEE/ClinVar-stars), chrM out of scope,
  CNV/SV, phenotype, co-segregation, ROH, somalier — all documented gaps. Steps 0/1/4 serial
  (documented perf). The recurrence null's **case-only** nature and "conservative" framing is
  already documented in `gene_burden.md` (P5 below is a *distinct* FDR-denominator bug, not that).

## Severity legend

- **major** — an incorrect/invalid result, a silent scientific loss, or a claim the code does not support.
- **moderate** — mis-calibration, a real best-practice gap, or a reproducibility/interpretability risk.
- **minor** — an edge case with small impact.

Effort: **S** ≈ a few lines; **M** ≈ a function or a design decision; **L** ≈ a new component.

## Summary table

| # | Title | Sev | File:line | Documented? | Effort |
|---|-------|-----|-----------|-------------|--------|
| **P1** | "CHARR" contamination gate is an uncorrected pooled ref-fraction → misses the 1–3% it claims to catch | **major** | `src/hprv/contamination.py:49`; `pipeline/00_qc.py:117,160` | partial | M |
| **P2** | Ingest-mode build/version checks only `warn` → a GRCh37 VEP VCF proceeds silently despite the "verifies the build" claim | **major** | `pipeline/02_annotate_sites.sh:170–176` | no | S |
| **P3** | `overall_pass` computed but consumed nowhere → contaminated / high-MIE / sex-mismatch trios feed Steps 5–6, not excluded | moderate | `pipeline/00_qc.py:213,226` | partial | M |
| **P4** | ClinVar P/LP rarity-rescue is void in [1e-2, 0.05): kept at Step 3, **zero** calls at Step 5, no audit counter | moderate | `src/hprv/selection.py:61`; `05_inheritance_screen.py:97` | partial | S |
| **P5** | De novo `dn_q_enrich` BH-corrected over called genes only, not the ~18–19k-gene model → q anti-conservative ~100× | moderate | `pipeline/06_gene_burden.py:306` | partial | S |
| **P6** | De novo Poisson pools LoF+missense and runs only the pooled test, never the per-class LoF test it advertises | moderate | `pipeline/06_gene_burden.py:251` | no | S |
| **P7** | Audit "recurrence sig / FDR q<" counters tally only the **dominant** family → under-report biallelic/X-linked | moderate | `pipeline/06_gene_burden.py:334–335` | no | S |
| **P8** | xlsx writes numeric cells as text → the advertised auto-filter sorts p/q/AF **lexicographically** | moderate | `src/hprv/report.py:40–41,62` | partial | S |
| **P9** | No per-run annotation provenance in any output table; xlsx hardcodes "gnomAD v4.1" regardless of the cache used | moderate | `src/hprv/report.py:145`; `src/hprv/audit.py` | no | S |
| **P10** | Union `norm` uses default `-c e` (abort) while per-trio uses `-c w` (warn) → late, trio-unattributable REF crash | minor–mod | `pipeline/01_make_cohort_sites.sh:172` vs `127,133` | partial | S |

---

## Findings

### P1 — the contamination gate does not measure contamination on a calibrated scale · major

**Files:** `src/hprv/contamination.py:47–49`; `pipeline/00_qc.py:116–117,160`

**Problem.** `charr()` returns the bare **depth-pooled** ratio `sum(ref_AD) / sum(ref_AD + alt_AD)`
over the member's high-quality hom-alt SNV sites, gated at `charr_threshold = 0.02`. Published CHARR
(Lu et al., *AJHG* 2023) is not this quantity: it is an **unweighted per-genotype mean** of
`ref_AD/DP`, then (a) divided by the mean reference-allele availability `≈ mean(1 − AF)` at those
sites and (b) baseline-corrected for reference-mapping bias / sequencing error. Two defects compound:
no AF/baseline correction, and depth-weighting instead of the per-genotype mean.

**Trigger / how to see it.** A member at true contamination α = 0.03. Hom-alt sites are HWE-weighted
toward high-AF variants, so `mean(1 − AF) ≈ 0.2–0.35`; the raw statistic ≈ `0.25·0.03 + baseline
(0.005–0.015) ≈ 0.012–0.017 < 0.02` → `contam_flag = 0`, the trio passes. Separately, a few
high-depth hom-alt loci in collapsed-repeat / segmental-duplication regions (paralogous ref reads)
dominate the depth-pooled sum and can inflate it (false positive).

**Consequence.** The gate's stated target — 1–3% contamination that "manufactures false hets /
comp-het second hits" — is not met; the raw statistic only reliably trips around ~5–8%. False
carriers then flow into Step 5 → Step 6 recurrence/FDR. The 0.02 cutoff was borrowed from CHARR
literature that thresholds a *corrected* value and does not transfer to this uncorrected quantity.

**Bounding (assess honestly).** This is the **fallback** path only: when `resources.selfsm_dir` is
set, verifyBamID2 FREEMIX is used instead and is the recommended production path. The docs do hedge
"CHARR proxy / reference-read fraction," but the specific sensitivity claim ("catches 1–3%") is
quantitatively wrong for the uncorrected statistic, and ROADMAP marks the CHARR screen "DONE."

**Suggested fix.** Step 0 runs *before* Step 2, so per-site gnomAD AF is not available there and the
full correction is infeasible at this stage. Either: **(a)** reframe honestly — rename to "raw
reference-read fraction at hom-alt sites," stop citing the 0.02 CHARR-scale threshold, and re-derive
`charr_threshold` empirically from clean-vs-spike-in trios for *this* pooled statistic; or **(b)**
move the contamination gate to post-annotation so per-site AF is available and implement
`mean(ref_AB)/(1 − AF)` with a baseline subtraction and a depth cap / winsorization. Switch the
accumulator to a per-genotype mean regardless, and correct the "1–3%" claim. **Effort M.**

**Referee angle.** A QC gate that is provably insensitive in its own stated target range, and an
implemented statistic that is not the CHARR statistic whose threshold it borrows.

### P2 — ingest-mode build/version check only warns · major

**File:** `pipeline/02_annotate_sites.sh:163–176` (only `:167` missing-`##VEP=` and the later
frequency-presence guard `die`; `:171` assembly and `:175` VEP-version only `warn`).

**Problem.** Ingest mode (`resources.vep.annotated_vcf` / `--vep-vcf`) is documented to *verify* the
build ("this script verifies that rather than trusting it"; CLAUDE.md "verifying build + frequency
presence"). But the assembly check and the VEP-version check emit `warn`, not `die`. A GRCh37 VEP VCF
built from a GRCh37 cache carries a valid `##VEP=` header and valid `vep_gnomAD*_AF`, so the only
enforcing guard (frequency presence) passes.

**Trigger.** Point `--vep-vcf` at a GRCh37 VEP VCF (or any header lacking the exact
`assembly="GRCh38"` / `assembly=GRCh38` token): one stderr warning, then execution continues with
`vep_vcf=$PRE_VEP`.

**Consequence.** Silent wrong result — Step 3 selects on GRCh37 coordinates, Step 4 transfers
annotations by position onto GRCh38 trios, so most sites fail to match and candidates are silently
lost or mis-annotated with **no hard failure**. Directly contradicts the fail-loud ethos and the
explicit "verifies the build" claim.

**Suggested fix.** Promote an assembly mismatch to `die` **when the header explicitly declares a
non-GRCh38 assembly** (and likewise a different VEP major version); keep `warn` only when the token
is *absent* — a legitimate GRCh38 VEP VCF may omit it. This delivers the documented guarantee
without false-failing valid inputs. **Effort S.**

### P3 — `overall_pass` is advisory-only; failing trios still feed the recurrence statistics · moderate

**File:** `pipeline/00_qc.py:213,226`. Consumers: `05_inheritance_screen.py:328–333` reads
`inferred_sex` only; nothing reads `overall_pass` / `mie_flag` / `contam_flag` (grep-confirmed).

**Problem.** Step 0 folds `mie_flag` / `sex_match` / `contam_flag` into `overall_pass` and the
docstring calls it a "per-trio QC gate (garbage-in guard)." But no downstream step consumes it —
Steps 1/2/4 run over the full `trios.resolved.tsv`, and Step 5 reads the QC report only for inferred
sex. `pipeline_design.md:147` claims the Step-0 list determines "which trios enter analysis," which
is false.

**Trigger.** A trio with real proband contamination or an elevated Mendelian-error rate (a subtle
sample swap): `overall_pass = 0`, yet its VCF passes through Steps 1/2/4, Step 5 emits false
inherited-het / comp-het calls, and those false carriers count in the Step 6 recurrence tally and
FDR — the exact artifact the QC was introduced to prevent.

**Consequence.** The QC "gate" gives no automated protection to the publication-critical recurrence
tier (produced *before* human review). A single contaminated or mislabeled trio can push a gene past
`burden.min_carriers`. Combined with **P1**, the recurrence tier has neither a sensitive
contamination detector nor an exclusion mechanism.

**Bounding.** The mermaid diagram says "flag suspect trios" (honest); the over-claim is in
`pipeline_design.md:147` and the `overall_pass` naming.

**Suggested fix.** Either genuinely gate — have **Step 6 at minimum** exclude `overall_pass == 0`
trios from recurrence counting (config-overridable), or quarantine them for review — **or** restate
`pipeline_design.md:147` and the README to make explicit that Step 0 is advisory-only and flagged
trios still contribute to calls and recurrence pending human review. Given P1, real gating is the
stronger fix. **Effort M** (gating) / **S** (doc restatement).

### P4 — ClinVar P/LP rarity-rescue is void in the [1e-2, 0.05) band · moderate

**Files:** `src/hprv/selection.py:60–65` (`rarity_ok = (fr is None) or (fr < rec_max) or plp`);
`pipeline/05_inheritance_screen.py:97–99` (`rare()` is frequency-only) and every mode's rarity gate.

**Problem.** Step 3 rescues a P/LP variant past the `too_common` gate for grpmax AF in [1e-2, 0.05),
tagging `clinvar_plp`. Step 5's `rare()` has no ClinVar bypass; every mode gates on `rare(v, rec_max)`
(<1e-2) or `rare(v, dom_max)` (<1e-4). So a P/LP allele at grpmax AF ≥ 1e-2 fails **all six** Step 5
modes. (The P/LP override still works for its *primary* purpose — bypassing the functional/CADD
ladder for fr<1e-2 P/LP variants; only the frequency-band rescue is nullified.)

**Trigger.** A recessive P/LP founder allele at grpmax AF 1–5% — GJB2 c.35delG (~1–1.5%), HBB HbS
(~4% AFR) — homozygous in the proband with carrier parents: kept by Step 3, then `hom_recessive`
requires fr<1e-2 and fails → **no call in any mode**.

**Consequence.** Silent loss of a homozygous known-pathogenic recessive call from
`candidates.calls.tsv`, the IGV export, and the xlsx, with **no Step 5 audit counter**.

**Bounding.** `recessive_max = 1e-2` is a deliberate recessive frequency ceiling, so dropping a >1%
recessive candidate is partly intended policy; the defect is the *contradiction* between Step 3's
uncaveated rescue promise and Step 5's silent drop. `allele_frequency.md:155–158` states Step 5 has
no ClinVar override, but only via a *dominant* 2e-4 example; the recessive corollary and the
"rescue accomplishes nothing" outcome are not surfaced, and `clinical_classification.md:76` promises
the override "rescues a variant from the rarity gate" uncaveated.

**Suggested fix (low effort).** Add a Step-5 audit counter for "clinvar_plp candidate dropped by
per-mode rarity" so the loss is not silent, and extend the `clinical_classification.md:76` caveat.
A larger, optional fix: propagate the P/LP flag into Step 5 to gate at `ba1` instead of `rec_max` for
known-pathogenic variants — but that is a conscious policy change, not a silent one. **Effort S.**

### P5 — de novo FDR denominator is called-genes, not the mutation model · moderate

**File:** `pipeline/06_gene_burden.py:306` (`bh_fdr([r['dn_p_enrich'] for r in rows])`), rows built at
`:213`; `dn_p_enrich` is non-None only when `gene in mut` and `exp > 0` (`:246–254`).

**Problem.** `dn_q_enrich` BH-corrects over `m = |genes-with-a-candidate-call ∩ mutation-table|`
(typically hundreds). The denovolyzeR/Samocha framework the docs invoke corrects over **every gene
in the mutation model** (~18–19k), where a zero-DNM gene is a valid null test at
`p = poisson.sf(-1, exp) = 1.0`, not an absent test. Those ~18k zero-DNM genes never get a row, so
they are silently dropped from `m`.

**Trigger.** `--mutrate` with 18k genes, `--n-trios 200`; gene X has one LoF DNM, expected 2e-4 →
`dn_p_enrich(X) ≈ 2e-4`. With m ≈ 250, `q ≈ 2e-4·250 ≈ 0.05` (a borderline "discovery"); the correct
exome-wide BH with m ≈ 18000 gives `q ≈ 1.0` (not significant).

**Consequence.** `dn_q_enrich`, presented as a BH FDR q-value, is anti-conservative by roughly
`n_model / n_called` (~100×) and can flag a single-DNM gene as FDR-significant when the true q ≈ 1.0.
This is **distinct** from the documented "uncalibrated Poisson mean" caveat, which concerns the
expectation `2·N·μ`, not the BH denominator.

**Bounding.** Off-by-default (no mutation-rate table ships); the fixed-threshold `dn_exome_wide_sig`
flag (`:308`, `p < 2.5e-6` Bonferroni) is m-independent and correct; docs already tell readers to
treat de novo p/q as "a rank, not a significance claim" (though they attribute it to the wrong cause).

**Suggested fix.** For the **de novo family only**, seed the BH input with `p = 1.0` for every gene
in `mut` with zero observed DNMs before calling `bh_fdr` (or pass an explicit `m = len(mut)`). Leave
the recurrence families' `m` unchanged — those are legitimately conditional on observing
≥ `min_carriers` carriers. **Effort S.**

### P6 — de novo Poisson pools LoF + missense into a single test · moderate

**File:** `pipeline/06_gene_burden.py:249–254` (`mu = mu_lof + mu_mis`, `obs = denovo_lof +
denovo_mis`, one `poisson.sf(obs-1, 2*n_trios*mu)`).

**Problem.** Pooling itself is not invalid (denovolyzeR's "prot" class is exactly a pooled
LoF+missense Poisson). The genuine defect is that the code runs **only** the pooled test and never
the per-class (especially LoF-specific) test denovolyzeR reports by default
(`denovolyzeByGene(classes=…)`, advertised at `gene_burden.md:89`).

**Trigger.** A haploinsufficient gene with 1 LoF DNM and `mu_mis ≫ mu_lof`: the single LoF
observation is diluted against a missense-dominated expectation, blunting the LoF signal relative to
a class-specific test.

**Consequence.** Power and interpretability loss for LoF-mechanism genes; the combined p cannot be
mapped back to a mechanism. Bounded: off-by-default, secondary, and per-gene DNM counts are tiny, so
no wrong headline number.

**Suggested fix.** Compute separate LoF and missense `poisson.sf` tests with separate output columns
(optionally an ACAT/Fisher combination), or explicitly document that Step 6 runs a single pooled
protein-altering ("prot") test and drop the per-class advertisement at `gene_burden.md:89`.
**Effort S.**

### P7 — recurrence "sig" / "FDR" audit counters tally the dominant family only · moderate

**File:** `pipeline/06_gene_burden.py:334–335` (`n_rec_sig` reads only `recurrence_exome_wide_sig`;
`n_rec_fdr` reads only `q_recurrence`), summary `:342–346`. The per-family columns exist (`:300–308`).

**Problem.** The biallelic and X-linked families (`recurrence_biallelic_*`, `recurrence_xlinked_*`)
are never included in the aggregate counters, yet the audit metrics
(`genes_recurrence_exome_wide_sig`, `genes_recurrence_fdr_sig`) and the stderr summary are labeled
generically.

**Trigger.** A gene with 3 distinct biallelic carriers and 0 dominant: `q_recurrence_biallelic ≪
0.05`, `recurrence_biallelic_exome_wide_sig = 1`, but `q_recurrence = None` — counted in neither
counter.

**Consequence.** `audit/summary.md` and the run log under-count significant recurrent genes; a paper
quoting `genes_recurrence_fdr_sig` omits every recessive-only / X-linked-only significant gene. The
authoritative `genes.ranked.tsv` is **complete and correct** — only the summary layer under-reports.

**Suggested fix.** OR the three inherited families per gene (any of the three `_exome_wide_sig`; any
of the three q < `fdr_q`), matching how `genes_recurrent` already uses the all-model universe — or
emit separate per-family counters and relabel the metric names. **Effort S.**

### P8 — xlsx numeric columns written as text break the advertised sort · moderate

**File:** `src/hprv/report.py:33` (`_read_tsv` → all strings), `:40–41` (`ws.append(r)` verbatim, no
coercion), `:59–62` (`freeze_panes`, `auto_filter.ref`).

**Problem.** Every cell is stored as text, and the same sheets enable interactive sort/filter. Step 6
`_fmt` emits `.4g` strings like `2.5e-06`, `0.0013`.

**Trigger.** Open `hprv_summary.xlsx`, use the Gene-consolidation filter to sort by
`q_recurrence` / `p_recurrence` / `loeuf` / `dn_q_enrich` (or Candidate calls by `grpmax_af` / `cadd`).

**Consequence.** Lexicographic sort: `0.0013` before `2.5e-06` (`'0' < '2'`), `1.2e-08` after `0.05`
— the most exome-wide-significant / rarest rows do **not** sort to the top on re-sort, and every
numeric cell shows Excel's "number stored as text" warning. No precision is lost and the default row
order is correct (pre-ranked recurrent-first); the defect manifests only on user re-sort.

**Suggested fix.** In `_style_data_sheet`, coerce each token before `append` — try `int()`, then
`float()`, fall back to the original string, preserve blanks. Typed cells sort numerically and clear
the warning. **Effort S.**

### P9 — no per-run annotation provenance; xlsx hardcodes "gnomAD v4.1" · moderate

**File:** `src/hprv/report.py:145` (literal `'gnomAD v4.1'`); `src/hprv/audit.py` records only
timestamp/step/scope/metric/value.

**Problem.** `genes.ranked.tsv`, `candidates.calls.tsv`, `hprv_summary.xlsx` and `audit/counts.tsv`
carry no VEP cache version, gnomAD version, CADD file version, or git commit / config hash.
`report.py:145` writes "gnomAD v4.1" unconditionally rather than reading it from the cache or the
`##VEP` header.

**Trigger.** Ingest a `--vep-vcf` built on a non-v4.1 gnomAD cache → the workbook still prints
"gnomAD v4.1." (Under the enforced self-run path, `--cache_version 115` + VEP r113+ ships gnomAD v4.1,
so the label is correct-by-construction there.)

**Consequence.** A reproducibility gap plus a narrow over-claim in the ingest edge case. The frequency
oracle version is load-bearing (golden rule 2) yet asserted, not recorded. Distinct from the
documented conda-lock / digest-pin TODOs (those are the *build* environment, not per-run *output*
provenance).

**Bounding.** `cohort.sites.annotated.vcf.gz` does carry the `##VEP=` header and the config records
`HPRV_VEP_VERSION`, so a run is *bindable* via its artifacts — the gap is that provenance is not
surfaced into the result tables.

**Suggested fix.** Emit one provenance record per run (git commit, VEP version, cache path/version,
CADD basenames, gnomAD version parsed from `##VEP`); surface it on the xlsx About sheet; and
**derive** the "gnomAD vX" label from the header rather than hardcoding it. **Effort S.**

### P10 — asymmetric REF-check across the two `norm` calls · minor–moderate

**File:** `pipeline/01_make_cohort_sites.sh:127,133` (per-trio `norm -m- -f REF -c w`) vs `:172`
(union `norm -d exact -f REF`, inheriting bcftools' default `-c e`).

**Problem.** Per-trio uses explicit `-c w` (warn, keep, never rewrite); the union omits `-c` and
inherits the default abort-on-mismatch. Empirically (bcftools 1.22/1.23) a REF-mismatched record
passes per-trio (exit 0, `REF_MISMATCH` warning), survives unchanged, then aborts the union with a
generic `Reference allele mismatch at chrN:POS` that names neither the trio nor the file.

**Trigger.** Any single REF-mismatched site — a trio called against a slightly different GRCh38 patch,
or an alt/decoy contig differing from `--ref`.

**Consequence.** The union fails *late* (after every per-trio site file is built) with a
locus-attributable but not trio-attributable error, and no audit counter records REF mismatches at
either step, so per-trio warnings are invisible on the happy path. (Not a silent-loss or
docs-contradiction — the `-c w` "never silently rewrite" claim is honored; this is an
asymmetry/ergonomics defect.)

**Suggested fix.** Make the REF-check mode explicit and identical at both steps — cleanest is `-c e`
at the **per-trio** step so a mismatch fails fast, attributed to the specific trio/file, before any
union work — and add a per-trio `ref_mismatch` audit counter. Document the enforcement point in
`cohort_construction.md`. **Effort S.**

---

## The three statistical / methodological surfaces a referee is most likely to attack

1. **Contamination QC is not on a calibrated scale (P1).** Needs a **code** fix (recalibrate the
   threshold or implement the AF/baseline correction post-annotation) **and** a doc correction to the
   "catches 1–3%" claim and the CHARR naming. Highest value: a reproducible, self-demonstrating
   insensitivity in the gate's own target range, feeding the recurrence null.
2. **De novo FDR corrected over the wrong gene universe (P5).** Needs a **code** fix — seed the BH
   family with p=1.0 for zero-DNM model genes (or pass `m = len(mut)`). Cheap; immediately visible to
   anyone reproducing against denovolyzeR/DeNovoWEST.
3. **QC flags never gate the statistics (P3).** Needs either a **code** fix (exclude/quarantine
   failing trios from Step 6) or a **doc / stated-assumption** restatement. A referee will ask "what
   does the contamination/MIE QC do to the results?" — the honest current answer is "nothing
   automated." Combined with P1, the recurrence tier has neither a sensitive detector nor an exclusion.

Secondary but likely: **P6** (pooled vs per-class de novo test) — a doc clarification if pooling is
kept, or a small code split into per-class Poisson tests.

---

## Considered and dismissed (refuted under verification)

Recorded so they are not re-litigated. Each was raised by a lens and **refuted** by both verifiers.

- **Cryptic relatedness / non-independent probands inflates the binomial null.** The plumbing is real
  (the only cross-trio guard is unique-proband-ID dedup, `resolve_trios.py:111`; no relatedness
  check), but it is a standard burden-test assumption, cross-sample somalier relatedness is already a
  documented roadmap item, and the recurrence p is already framed as "a rank, not significance."
  **Recommendation (not a bug):** add a one-line *stated assumption* — "probands are assumed unrelated"
  — to `gene_burden.md`, ideally alongside the somalier-relatedness roadmap note.
- **Cohort union not byte-reproducible across manifest/trio order.** Refuted: `bcftools sort` +
  `norm -d exact` make the annotated union coordinate-ordered and order-invariant; downstream calls
  are deterministic. (This is a positive — worth keeping a test.)
- **grpmax proxy takes `max` over exome AND genome, reintroducing small-AN inflation.** Refuted: it is
  the documented conservative combiner, and the bottlenecked founder groups (ami/asj/fin/mid) are not
  lifted into INFO at all, so they cannot leak in.
- **Same-variant (founder/artifact) recurrence still receives a significant p/q.** Refuted: it *is*
  scored, but transparently labeled `recurrence_kind = same_variant` and ranked below distinct-variant
  recurrence, which is the documented, intended treatment.
- **Step 3 keep-reason histogram uses first-match precedence.** Refuted as a minor interpretation
  nuance, not a defect.
- **CHARR GQ/DP site gate biases the estimate downward by excluding informative sites.** Refuted as a
  distinct mechanism (the real CHARR issue is P1's missing correction, not the QC site filter).

## Verified correct (no action)

Checked during verification and sound as designed:

- **Recurrence families' BH denominator** (`p_recurrence` / `_biallelic` / `_xlinked`) — legitimately
  conditional on observing ≥ `min_carriers` carriers; correcting over called genes only is
  appropriate here (unlike the de novo arm in P5).
- **`dn_exome_wide_sig`** fixed-threshold flag — `p < 2.5e-6` Bonferroni, m-independent and correct;
  unaffected by the P5 BH-denominator bug.
- **Every per-gene column in `genes.ranked.tsv`** — all three inherited families' sig/q columns are
  complete and correct (P7 affects only the aggregate audit counters; P8 only the xlsx re-sort).
- **The binomial tail** `binom.sf(n-1, N, p)` = P(≥ n carriers) — off-by-one correct.
- **The P/LP override's primary function** — bypassing the IMPACT+CADD functional ladder for fr<1e-2
  P/LP variants works as intended (only the frequency-band rescue in P4 is inert).
- **`p_biallelic_hwe = (Σq)²`** and **`p_carrier_hwe = 1 − Π(1−q)^ploidy`** — valid rare-variant HWE
  approximations, erring conservative; the independence-across-variants assumption is standard.

## Provenance of this review

Third of three review passes. Passes 1–2 (whole-repo audit; deep inheritance-model review) are
merged into `main`; their fixes are listed under "Scope — what is NOT here." This pass focused on the
areas those touched only lightly: Step-6 statistics, Step-0 QC science, and data-integrity in the
normalization/annotation path. Reviewed against `main` at the time of writing (post inheritance-model
merge). Findings are code-path arguments, not observed failures — validate P1 and P5 against
spike-in / mock data before acting on the calibration specifics.
