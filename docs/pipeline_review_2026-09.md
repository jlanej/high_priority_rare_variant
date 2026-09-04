# Pipeline review — scientific / logical / coding / structural pass (2026-09)

Point-in-time review of `high_priority_rare_variant` at `main` (a6359b2). Same conventions as
`docs/pipeline_review_2026-07.md`: a findings log, code is ground truth, each finding carries a
file:line and a concrete failure scenario. Items the 2026-07 pass fixed, dismissed, or listed as
documented residuals are not re-raised unless the code has since regressed.

How it was produced: the core Python (annotations, genotype, selection, Steps 0/3/5/6/9,
prioritize.py, igv.py, report.py, resolve_trios) and the Step 1/2/4 shell were read in full by the
primary reviewer; three delegated sweeps covered the remaining shell layer (Steps 8, 2b, run_pipeline,
common.sh, prepare_resources, Slurm, Dockerfile/CI), docs-vs-config-vs-code consistency, and the test
suite (which was also executed). Every finding below was verified against the source; the
load-bearing ones (R1, R3, R4, R5) were reproduced by executing the function in question.

Severity: **major** = incorrect or silently lossy result / an over-claim the outputs make;
**moderate** = mis-calibration, real best-practice gap, or a reproducibility risk; **minor** = bounded
edge case or hygiene.

## Resolution status — 2026-09 follow-up

Every finding below was re-verified against the code (the load-bearing ones by executing the
function or script in question) and then addressed on branch `fix/pipeline-review-2026-09`.
One further defect surfaced while fixing R1/T7 and is recorded as **R15** below.

| # | Status | What changed |
|---|--------|--------------|
| **R1** | **fixed** | `prioritize.is_hom_alt_call` decides zygosity from the base-form GT (`T/T`), the index form, or the mode; every hom-alt call now passes its own AB band. Unit test with base strings; integration asserts `gt_qc_pass=1` on every hom_recessive / x_linked_recessive row. |
| **R2** | **fixed** (rank) + reframed (labels) | Step 6 computes `mu_tot`, `exp_carriers_mu` (C over the FULL mutational-target table), `carrier_excess_ratio`, `p_carrier_excess`, `rank_basis`, and orders recurrent genes by the size-normalised p when `burden.rank_by_mutational_target` (default true; `--mutational-target` / `resources.mutation_rate_table`). `p_recurrence`/`q_recurrence`/`*_exome_wide_sig` keep their names but are documented everywhere as a case-only RANK; the workbook and audit summary no longer say "calibrated". |
| **R3** | **fixed** | Step 9 `n_observed` counts DISTINCT (trio, variant) observations; `n_rows` carries the raw row count; `max_downweight_fraction` and `per_trio` use the distinct unit. Integration asserts GENE1 `n_rows > n_observed = 6`. |
| **R4** | **fixed** | `rarity_driven_by_single_group` is fed `grpmax_af`, never the oracle value. |
| **R5** | **fixed** | `HEMIZYGOUS_MODES` are X-linked observations for MOI coherence: coherent with any XL* curation, neutral otherwise, never charged the autosomal discordance. |
| **R6** | **fixed** | `report.py` reads `annotations.rarity_oracle(cfg)` and describes faf95 or the proxy accordingly; provenance line records `rarity_oracle.<arm>`. Test reads the About sheet under both oracles. |
| **R7** | **fixed** | Ingest mode dies on any multi-ALT record and on union sites absent from the external VCF (`HPRV_INGEST_ALLOW_PARTIAL=1` downgrades to a WARN + audit row). |
| **R8** | **fixed** | Step 6 de novo arm: missing `mu_mis` → no test; missing `mu_lof` → imputed via `mu_lof_impute_factor`; `dn_mu_src` reports gnomad / imputed / none. |
| **R9 / D1–D9** | **fixed** (docs) | CLAUDE.md, README and every doc rewritten to the one-oracle contract; dead config keys removed; dangling references removed; numbers aligned to the canonical table (2.41×, floor=3 tier table, ~80%). |
| **R10** | **fixed** | `qc.max_sites` (default 200000) is a config key, passed through Step 0, and documented as a capped scan. |
| **R11** | **fixed** | `outputs.igv.nonhuman_screen.flag_fraction` is the ONE threshold read by Step 8 (`--nhf-flag-fraction`) and Step 9; the retired `prioritization.composite.nhf.flag_fraction` is honoured as a fallback only. |
| **R12** | **fixed** | X-linked null uses the male proband count (`--n-male-trios` / `--qc-report`); `recurrence_kind` compares per-trio variant SETS. |
| **R13** | **fixed** | `sig_caf_low` requires a mutational-target row (`in_mutrate_table`). |
| **R14** | **fixed** | Stray `.bgz.tbi` untracked, `.gitignore` covers `*.vcf.bgz` / `*.bgz.tbi` / `*.bgz.csi`, stale worktrees removed. |
| **R15** (new) | **fixed** | `genotype.dp()` falls back to FORMAT/DP: cyvcf2 derives `gt_depths` from AD sums whenever the header declares AD, so an AD-less (ref-block) parent read depth −1 and `clean_parent` failed CLOSED on depth — a de novo with a ref-block parent was silently dropped, not flagged. Found by the new integration fixture (T7). |
| **S1** | **fixed** | Per-contig shard `.done` files carry a VEP-input key (union cksum + cache dir/version + plugin file stats); a stale shard is re-annotated, and `gather` refuses stale shards. |
| **S2** | **fixed** | `HPRV_GNOMAD_ORACLE` is exported; under faf95 Step 2 dies on a missing slim or a failed transfer (also on the ingest path), run_pipeline's preflight applies the slim gate on both paths, and Step 3 asserts the `gnomad_AF_joint` witness is declared. |
| **S3** | **fixed** | Steps 1 and 4 key their per-trio caches on the source VCF (size+mtime), FILTER/contigs and samples; the union key covers every per-trio key; Step 0's cache is keyed on the manifest. |
| **S4** | **fixed** | Value-level 0-lift guards for SpliceAI and CADD after split-vep. |
| **S5** | **fixed** | `emit_sh` normalises YAML booleans; Step 2 also accepts `True`. |
| **S6** | **fixed** | Step 8 sub-tasks die without a kraken2 DB; `phase.sbatch` requires the NHF manifest to exist and allows `DOWN_TO=9`; the CI smoke probes fail the job; CI shellchecks `scripts/` and fails on a skipped pure run. |
| **S7** | **fixed** | `.done` key includes `csq_select`/cache/version/oracle; `HPRV_TMPDIR` bound for Steps 2 and 8; `--nhf-gather` no longer audits `minicrams=0`; task-private tmpdir survives `--cleanenv`; `prep_dbnsfp` SIGPIPE trap closed; `verify` reports required-by-default resources as required and `emit-env` emits the slim when present. |
| **T1–T9** | **fixed** | Base-form GT fixture; a trim-loop fixture that actually trips the HALT; code-default ↔ config probes in both directions; `_Skip` instead of a silent PASS; the cache re-run asserts "cached"; per-run temp paths; new unit harnesses for Step 5 (four tests), genotype QC (two), Step 6 counting/ranking, the Step-9 fixes and the workbook; the mock gains a ref-block parent, a half-called parent, `mut_syn`, and both alleles of the multiallelic site in the slim; integration asserts the mode-level and audit-level contracts listed under T6. |

Deliberately NOT changed: the Step-5 emission of one row per pair leg (the review table's design; Step 9 now counts distinct observations instead), the column names `p_recurrence` / `q_recurrence` / `*_exome_wide_sig` (renaming would break every downstream reader for a cosmetic gain; they are documented as a rank), the `parent_max_alt_ad` absolute count, and the unpinned `micromamba latest` / LOFTEE branch tip / `action-shellcheck@master` (pinning needs a network check of the exact versions).

*(The sections below are the original review, preserved as the findings log.)*

## Summary table

| # | Title | Sev | Where |
|---|-------|-----|-------|
| R1 | Step 9 genotype QC never recognises a hom-alt call: `child_gt` is bases (`T/T`), the test is for `1/1` → every hom_recessive / x_linked_recessive / denovo_x_hemi call is penalised −2 and flagged `het_AB … outside` | **major** | `src/hprv/prioritize.py:1179`, `pipeline/05_inheritance_screen.py:97` |
| R2 | Step 6 "recurrence significance" saturates on private variants and the headline ranking is effectively a gene-size ranking; outputs still say "calibrated" / "exome-wide significant" | **major** (over-claim; scientific) | `pipeline/06_gene_burden.py:259-273,332-339`, `src/hprv/report.py:198` |
| R3 | Step 9 `n_observed` counts rows; Step 5 emits each comp-het leg once per pair (k×m pairs → 2km rows for k+m variants) and a variant can appear in up to three modes → the excess statistic's numerator is inflated quadratically for genes with many hets per proband | **major** | `pipeline/05_inheritance_screen.py:332-360`, `pipeline/09_prioritize.py:520` |
| R4 | `rarity_driven_by_single_group` is fed `rarity_af` (a CI lower bound under the default oracle) instead of the point-estimate proxy its docstring describes; fires on every `zero_ci` row | moderate | `src/hprv/prioritize.py:1589-1590` vs `1130-1145` |
| R5 | MOI coherence treats a hemizygous male call as "recessive": `x_linked_recessive` in an XLD gene is charged −1.0; `denovo_x_hemi` in an XLR gene is labelled `moi_mismatch_het_in_recessive_gene` | moderate | `src/hprv/prioritize.py:1263-1330` |
| R6 | xlsx About sheet hardcodes the retired oracle description ("grpmax-proxy AF … POINT ESTIMATE … from the VEP cache; the sole population-frequency source") regardless of `resources.gnomad.oracle` | moderate | `src/hprv/report.py:204-210,223` |
| R7 | `--vep-vcf` ingest path never checks that the external VCF is biallelic or covers the union; Step 4's `isec -c none` then drops every uncovered site silently | moderate | `pipeline/02_annotate_sites.sh` (ingest branch), `pipeline/04_subset_and_annotate_trios.sh` (`isec`) |
| R8 | Step 6 de novo enrichment substitutes 0.0 for a missing `mu_lof`/`mu_mis` — the repo's own "a missing component is never zero" rule | moderate | `pipeline/06_gene_burden.py:279-280` |
| R9 | The frequency-oracle contract is described two incompatible ways across CLAUDE.md, README, six docs, the example config and code comments; plus a set of dead config keys | moderate (docs) | see §R9 |
| R10 | Step 0's MIE and CHARR are computed over the first 200,000 QC-passing sites (an undocumented, non-configurable `--max-sites`), not genome-wide as documented | moderate | `pipeline/00_qc.py:142,165`, `pipeline/run_pipeline.sh:213` |
| R11 | `nhf_flag` (Step 8) uses a hardcoded 0.5 while `nhf_status` (Step 9) reads `prioritization.composite.nhf.flag_fraction` → the two columns can disagree in one table | minor | `src/hprv/igv.py:25,139` |
| R12 | Step 6 X-linked null uses N = all trios for a hemizygous-male model, and `recurrence_kind` reads a shared comp-het PAIR as `distinct_variant` | minor | `pipeline/06_gene_burden.py:271-272,308-311` |
| R13 | `sig_caf_low` flags any gene with no mutational-target row as a "low-information locus", adding a corroboration vote to every unjoined gene | minor | `src/hprv/prioritize.py` (`artifact_signals`, caf branch) |
| R14 | Repo hygiene: a gnomAD `.vcf.bgz.tbi` index tracked at the repo root (the ignore pattern only covers `.vcf.gz.tbi`); three stale detached worktrees under `.claude/worktrees/` | minor | repo root, `.gitignore` |
| S1 | Step 2's per-contig VEP shard cache is existence-keyed: a changed union (trio added) or a newly supplied plugin file is annotated from stale shards and the output is re-stamped with the NEW key (reproduced) | **major** | `pipeline/02_annotate_sites.sh:166,182-197,614` |
| S2 | `oracle: faf95` is enforced by file existence only: a failed slim transfer WARNs and continues, and the `annotated_vcf` ingest path skips the slim preflight entirely → every variant reads "rarest", BA1-common alleles pass every gate, exit 0 | **major** | `02_annotate_sites.sh:566-569`, `run_pipeline.sh:105-126` |
| S3 | Per-trio caches in Steps 1 and 4 are keyed on `trio_id`, never on the source VCF; a replaced/re-called trio VCF is never re-ingested | moderate | `01_make_cohort_sites.sh:120`, `04_subset_and_annotate_trios.sh:149` |
| S4 | SpliceAI/CADD have header-only guards: a present-but-inert score file (wrong build, bad tabix) yields declared-but-empty fields and no warning, while `spliceai_required` reports satisfied | moderate | `02_annotate_sites.sh:398-409`, `run_pipeline.sh:145-148` |
| S5 | YAML `shard_by_contig: true` is exported as the string `True`, which Step 2's `case` reads as OFF → the shipped config runs VEP un-sharded and un-resumable (reproduced) | moderate | `src/hprv/config.py:107`, `02_annotate_sites.sh:335` |
| S6 | SLURM: `nhf-plan` goes green having scheduled nothing when `--kraken2-db` was omitted; the phase graph caps at Step 8 so Step 9 never runs; the CI "real VEP" smoke test cannot fail | moderate | `08_igv_export.sh:357-363`, `slurm/phase.sbatch:108,153-155`, `.github/workflows/build-publish.yml:251-271` |
| S7 | Step 2's `.done` key omits `csq_select`/cache path; Step 0's cache is bare (a new trio gets no inferred sex); `HPRV_TMPDIR` is not bound for Steps 2/8 in host-wrapped mode; `--nhf-gather` audits `minicrams=0` | minor–moderate | see §Shell |

Test-suite (T1–T9) and docs (D1–D9) findings follow the detailed sections.

## Findings

### R1 — Step 9's genotype QC cannot see a homozygous call · major

`05_inheritance_screen.py:97` writes `child_gt` from cyvcf2 `gt_bases`, i.e. allele strings such as
`A/T` or `T/T`. `prioritize.py:1179` decides zygosity with `hom_alt = gt in ("1/1", "1|1")`, which is
never true on real data. Every hom-alt call — `hom_recessive`, `x_linked_recessive`,
`denovo_x_hemi`, and any hemizygous male call — therefore takes the het branch, its AB (~1.0) falls
outside 0.25–0.75, `gt_qc_pass` reads 0 with reason `het_AB=0.98_outside[0.25,0.75]`, and
`pts_quality` carries −2.0. Reproduced:

```
genotype_qc({"child_gt": "T/T", "child_AB": "0.98", ...})  -> (False, 'het_AB=0.98_outside[0.25,0.75]')
genotype_qc({"child_gt": "1/1", "child_AB": "0.98", ...})  -> (True, '')
```

Consequence: every recessive candidate is systematically down-ranked by 2 points relative to
dominant calls and carries a false QC failure in the review table. The unit test
(`tests/test_pure.py:1182`) uses the numeric form, and `assert_integration.py:388` scopes its quality
check to rows with `gt_qc_pass == "1"`, so neither can observe it.

Fix (S): derive zygosity from the two alleles (`a == b and a != ref`, splitting on `/` or `|`) or from
`inheritance in {hom_recessive, x_linked_recessive, denovo_x_hemi}`; and add a test whose fixture
uses base strings.

### R2 — the recurrence "significance" is a gene-size ranking wearing a p-value · major (scientific)

`gene_burden.md` already says the case-only null is anti-conservative and "read significant as a
ranking". Two things go beyond what that caveat covers.

First, the magnitude. With `absent_af_floor = 1e-6`, carriers of *distinct private* variants give
(`06_gene_burden.py:259-273`):

| N trios | carriers | p_recurrence |
|---|---|---|
| 200 | 2 | 3.2e-7 |
| 1000 | 2 / 3 | 8.0e-6 / 3.6e-8 |
| 3000 | 3 / 4 | 9.6e-7 / 1.4e-8 |

So `recurrence_exome_wide_sig` is essentially "≥ 3 carriers of private functional hets", which in a
cohort of hundreds of trios is true of most large genes. Because p is monotone in carrier count for
private variants, `genes.ranked.tsv` — the workbook's headline sheet — orders genes by how many rare
functional hets they collect, i.e. by mutational target: TTN, MUC16, SYNE1 first. Nothing in Step 6
normalises by gene size; the mutational-target offset exists only in Step 9 and never feeds back.

Second, the outputs still over-claim: `report.py:198` prints "exome-wide p < 2.5e-6 on the
**calibrated** recurrence null", the audit summary reports "recurrence exome-wide significant: N",
and the column names (`recurrence_exome_wide_sig`, `q_recurrence`) invite a referee to read them as
tests. The 2026-07 review dismissed the BH-denominator question for these families as "legitimately
conditional"; that holds only if the p-values were themselves valid, which the table above shows they
are not.

Fix (M): (a) rank Step 6 by an offset-normalised statistic — the cheapest is to reuse Step 9's
`E_g = C·μ_g` (already computed from the same gnomAD v2.1.1 table) as the expected carrier count, or
to divide by CDS length as a stop-gap; (b) rename the columns to `recurrence_rank_p` or drop the
`*_exome_wide_sig` flags until a per-gene cumulative-frequency null (TRAPD/CoCoRV-style) exists;
(c) fix the "calibrated" string in `report.py`. Until then the honest headline table is
`genes.prioritized.tsv`, not `genes.ranked.tsv`.

### R3 — Step 9's numerator counts rows, and Step 5 multiplies rows · major

`05_inheritance_screen.py:332-360` emits both legs of *every* mat×pat (and mat/pat×denovo) pair, so a
gene where one proband carries k maternal and m paternal sub-1e-2 hets yields 2·k·m `compound_het`
rows for k+m distinct variants. The same variant can additionally appear as a `dominant` row (an
unphased de-novo-partner pair does not consume its legs, `:340-345`) and as a `denovo` row.
`09_prioritize.py:520` then sets `n_observed` = number of rows per gene, which drives
`excess_ratio`, `p_nb`, `per_trio` (the saturation signal), `recurrence_shape` and every tier rule.

The inflation is quadratic in hets-per-proband, i.e. it hits long genes specifically — the class
the artifact panel is calibrated on (the "TTN/SYNE1/DNAH11 clear it at 1.9–3.7x" figure in
`prioritize.py` is therefore partly this bookkeeping, not biology). A proband with 4 maternal and 4
paternal TTN hets contributes 32 rows for 8 observations.

Fix (S): count distinct `(trio_id, chrom, pos, ref, alt)` per gene in Step 9 (and use the same
distinct set for `sites`/`site_trios`); optionally emit each comp-het leg once with a
semicolon-joined `pair_id` list in Step 5. Re-derive the validation-cohort tier table afterwards —
the down-weight set will change.

### R4 — the single-group flag compares a CI lower bound to a point estimate · moderate

`rarity_driven_by_single_group` (`prioritize.py:1130-1145`) is documented and written for
`grpmax_af` (point estimate over eligible groups) vs `max_af` (point estimate over all groups).
`score_variant` passes `af_col` (`:1589`), which under the default `faf95` oracle is the CI lower
bound — 0.0 on every `zero_ci` row. Reproduced: `rarity_driven_by_single_group("0", "2.2e-4")` →
True, so the flag fires on essentially every gnomAD-observed low-count allele and says nothing about
population structure. Fix (S): pass `row.get("grpmax_af")`.

### R5 — MOI coherence mislabels hemizygous male calls · moderate

`RECESSIVE_MODES` includes `x_linked_recessive`, and `moi_coherence` maps it to `observed =
"recessive"`; `denovo_x_hemi` maps to `"dominant"`. A hemizygous male genotype is compatible with
both XLR and XLD, so: `moi_coherence("x_linked_recessive", "XLD")` → `('discordant', '')` and
charges −1.0 with no caveat; `moi_coherence("denovo_x_hemi", "XLR")` → caveat
`moi_mismatch_het_in_recessive_gene` for a call that is not a het. Fix (S): treat both hemizygous
modes as coherent with any `XL*` curation and neutral otherwise.

### R6 — the workbook documents the wrong oracle · moderate

`report.py:204-210` ("grpmax-proxy AF … A POINT ESTIMATE: faf95 … needs AC/AN, which the cache does
not carry") and `:223` ("from the VEP cache; the sole population-frequency source") are emitted
unconditionally. Under the default run every gate used faf95 from the gnomAD joint slim. Since the
About sheet is the methods record a collaborator reads, this is the 2026-07 P9 problem in a new
place. Fix (S): read `annotations.rarity_oracle(cfg)` and the audit's `rarity_oracle.*` row.

### R7 — ingest mode can silently discard union sites · moderate

With `resources.vep.annotated_vcf`, Step 2 verifies the `##VEP` build/version but never that the file
(a) is biallelic or (b) contains every site of `cohort.sites.vcf.gz`. Step 3 selects from it and
Step 4 intersects each trio with `bcftools isec -c none -n=2`, so a site absent from the external
VCF — or present only as a multiallelic record — is dropped from every trio with no counter. With
`--flag_pick` on a multiallelic record, `+split-vep -s pick` also assigns one allele's consequence to
the site. Fix (S): in ingest mode, die on any record with >1 ALT and on
`isec -n=2` count < union count (or WARN with the number lost and an audit row).

### R8 — de novo enrichment charges a missing rate as zero · moderate

`06_gene_burden.py:279-280`: `mu = (_num(mu_lof) or 0.0) + (_num(mu_mis) or 0.0)`. A gene with a
blank `mu_lof` (≈505 of 19,643 in gnomAD v2.1.1) is tested against a smaller expectation with no
provenance mark — the exact pattern CLAUDE.md records as fixed in Step 9. Fix (S): skip the test (or
impute as Step 9 does) when either rate is missing.

### R9 — one contract, two stories · moderate (docs, but the kind that changes what a user runs)

The code implements ONE oracle per run (`annotations.py:274-336`): `faf95` by default, halting
without the slim; `grpmax_proxy` only when configured; the arms never cross; `rarity_basis` carries
`measured|zero_ci|absent`. The retired per-variant design ("prefers faf95, falls back to the proxy
wherever gnomAD published none; `rarity_oracle` reports which fired per variant") is still the text
of: CLAUDE.md lines 61-74 and 405-413 (directly contradicting lines 40-48 and its own provenance
table), `docs/README.md:60,73,89,185,230`, `docs/allele_frequency.md:31-37,186-192,257-265`,
`docs/inheritance_and_genotype_qc.md:8-23`, `docs/gene_burden.md:19-24,119-120,285`,
`docs/pipeline_design.md:229-241`, `README.md:81-93`, `config/config.example.yaml:271-278,704-706`
(contradicting its own 169-194), and the docstrings of `grpmax_af`/`faf95` in `annotations.py`,
`05:33-38`, `06:173-178`, `prioritize.py:1101-1105,1557-1562`. CLAUDE.md also duplicates the MAX_AF
bullets verbatim and numbers golden rule "6" between "2" and "3".

Related config drift (none of these keys is read anywhere): `composite.gene_list_prior.points`
(the live key is `composite.weights.gene_list_prior`), `gene_list_prior.combine_gene_and_set`,
`signals.segdup.min_identity` (pinned as a "code default" by `test_pure.py:2345`),
`excess.offset.covariate_adjust`, `excess.offset.components`, `excess.refit_per_cohort`,
`composite.emit_both_rankings`, `signals.segdup.track_build`, `scope.*`. Also:
`prepare_resources.sh verify` reports gnomAD/REVEL/AlphaMissense as "not required" while
`run_pipeline.sh:124-136` halts without them; `docs/resources.md:427` says Step 2b is ON by default
(config and code say off); the Poisson anti-conservatism figure is 2.51× in README/code comments and
2.41× in the canonical table; `gene_is_constrained`'s docstring says s_het is not read while the body
reads it, and Step 9 never reads `phaplo_min` although Step 6 does.

### R10 — "genome-wide" MIE is the first 200k sites · moderate

`00_qc.py:165` defaults `--max-sites` to 200,000 and `run_pipeline.sh:213` never overrides it; the
loop breaks at `:142` once that many QC-passing biallelic autosomal sites are seen, and the CHARR
sums stop with it. On WGS that is a slice of chr1. Probably adequate as a swap detector, but it is
neither documented (`docs/README.md:195`, `inheritance_and_genotype_qc.md:34` say genome-wide) nor
configurable, and a localized MIE cluster (UPD/CNV) elsewhere is invisible. Fix (S): expose
`qc.max_sites` and say what it does; consider sampling across contigs rather than a prefix.

### R11–R14 · minor

- **R11** `igv.py:25` `NHF_FLAG_FRACTION = 0.5` vs `prioritization.composite.nhf.flag_fraction`; pass
  the config value into `build_variants_tsv`.
- **R12** `06:271-272` the X-linked binomial uses `n_trios` (all probands) for a hemizygous-male
  model, and the same column pools hom-alt daughters; `:308-311` two trios sharing one comp-het pair
  read `distinct_variant` because each leg is a key.
- **R13** `artifact_signals`: `classic_caf is None` is flagged when *any* cutoff exists, including
  genes with no mutrate row at all, so `corroboration_count ≥ 1` for every unjoined gene.
- **R14** `gnomad.joint.v4.1.sites.chr22.vcf.bgz.tbi` is tracked at the repo root (`.gitignore` has
  `*.vcf.gz.tbi` but not `*.bgz.tbi`); `git worktree list` shows three detached worktrees under
  `.claude/worktrees/` holding old copies of the tree plus `__pycache__` — prune them.

## Observations that are design choices, not defects (stated so they are not re-litigated)

- `frequency()`'s `zero_ci → 0.0` rule is the right reading of gnomAD's FAF (Hail returns 0 for
  AC ≤ 1, and `fafmax` is emitted missing when every FAF group is 0); the `gnomad_AF_joint` witness is
  the correct discriminator. The Step-5 orphan audit (`rarity_faf95_absent_but_cache_has_af`) is a
  good partial-join guard.
- The slim keeps gnomAD-filtered (AS_VQSR/AC0) records and transfers their AF; that errs toward
  treating recurrent artifacts as common, which is the conservative direction here. Worth stating.
- `dominant` rows with `origin=both` (HET×HET parents for a < 1e-4 allele) are near-impossible in
  unrelated parents and usually mean a founder allele, consanguinity or a systematic artifact; a flag
  would help review but nothing is wrong.
- `parent_max_alt_ad = 1` is an absolute count; 1 alt read in 10 and 1 in 200 pass identically.
  A fraction-based cleanliness test is the usual refinement.
- VEP already runs with `--numbers`; lifting `EXON`/`INTRON` (and `CDS_position`) through the
  split-vep `want` list would make the NMD test — and V5 — reachable with no new resource.

## Verified correct

- Numerics in `prioritize.py`: `nb_sf`/`pois_sf` agree with explicit pmf sums to 6+ digits on
  moderate inputs; the extreme tail reproduces 1.6e-241 where the Poisson complement returns 0.0;
  NB at α=1e-4 matches Poisson (the continued fraction converges at θ=1e4); `nb_midp` is ordered
  correctly; `bh_fdr` matches a hand computation with a `None` pass-through.
- `annotations.frequency()`/`rarity_basis()` implement the documented one-oracle contract; the
  `_max_float` tuple handling; `clinvar_stars` canonicalisation; `GRPMAX_POPS` excludes the founder
  groups; MAX_AF / global AF are never read by a gate.
- GRCh38 PAR coordinates in `genotype.py`; `strict_gt=True` in Step 5; chrY routed away from the
  mother-keyed models; HOM_ALT-parent deterministic origin; `--keep-sum AD` for `1/2` legs; chrM
  excluded at the union and counted in the audit.
- Step 1's manifest-keyed `.done`, Step 2's resource-keyed `.done`, Step 4's plausible-keyed
  `.done`, the two 0-match `die`s, and the three preflight halts all exist and fire as documented.
- Step 6 `binom.sf(n-1, N, p)`, `p_carrier_hwe`, `p_biallelic_hwe`, and the de novo BH padding.
- Step 9's never-drop assertions, symbol-first join, per-table constraint precedence, gene-prior
  MAX-not-sum, non-germline exclusion, CSV sniffing and headerless hard stop.

## Tests (delegated sweep; runs were executed on this host)

Runs: `py_compile` on all 24 python files and `bash -n` on all 18 shell files pass. Bare
`python3` (3.14, no yaml): `59 passed, 9 SKIPPED … NOT full coverage`, exit 0. With a yaml-capable
python (`/usr/bin/python3` 3.9): `All 68 pure-logic tests passed`. The five shell tests (Step-8
idempotency, NHF carriers, NHF scatter/gather, SLURM NHF phases, plan-resubmit) all pass. The
integration suite could not run here: no host python has cyvcf2/pysam/scipy (the micromamba
bootstrap this machine used to have is gone).

Findings, ranked:

- **T1** The suite encodes R1 as correct: `test_pure.py:1181-1182` asserts the `"1/1"` form, and
  `assert_integration.py:388` scopes the quality check to `gt_qc_pass == "1"` rows. Add a fixture with
  base-form GTs (the NHF shell tests already document that these columns are base-form).
- **T2** `test_prioritize_trim_loop_halts_on_mis_specified_null` (`test_pure.py:897-917`) cannot
  fail: its fixture is uniformly inflated, so C absorbs it and the else-branch asserts `0.0 <= 0.05`.
  Deleting the NB-arm HALT (`prioritize.py:483-488`) leaves all 68 tests green (mutation run). The
  Poisson test's own fixture (`:1974-1976`) does trip the NB HALT at 10% and can be reused.
- **T3** `test_prioritize_config_matches_canonical_defaults` (`:2322-2394`) compares config values to
  literals inside the test; only `default_weights()` is compared to code. Drifting
  `cadd_benign_max`, `t1_watch.max_q` and `spliceai_benign_max` in code leaves all 68 green
  (mutation run). Golden rule 4 is asserted in one direction only.
- **T4** Two tests (`:2308-2313`, `:2327-2331`) catch `ImportError` and `return`, printing PASS while
  asserting nothing; the bare run's "59 passed" is really 57. Use the file's own `_Skip`.
- **T5** `test_prioritize_never_drop_row_conservation` (`:1871-1874`) claims a second run is a
  cache no-op, but re-invokes without `--constraint/--established-genes/--n-trios` so the content key
  changes and it re-runs; the assertion passes either way.
- **T6** Integration assertions that pass on silently-wrong output (`assert_integration.py`):
  `:134` `any(mode == "dominant")` over two cis GENEC variants; `:170-172` never asserts GENEMID
  survived Step 3; `:204` passes if GENE1 is absent from `genes.ranked.tsv`; no assertion on the total
  call count / mode multiset (only `len(vrows)==len(calls)`, an identity); not asserted:
  `gene_universe_zero_count`, the Step-5 orphan counter `rarity_faf95_absent_but_cache_has_af`,
  or `gt_qc_pass` on any biallelic row.
- **T7** Mock-vs-real gaps in `make_mock_data.py`/`mock_vep.py`: FORMAT is always `GT:AD:DP:GQ`
  (no ref-block `0/0` parents, so `sample_qc_ad_measured` and the `*_unmeasured` flags never run
  end-to-end); GT alphabet is `{0/0,0/1,1/1,0/2,1/2}` — no `./.`, no `0/.` (the `strict_gt=True`
  claim is untested), no phased `|`, no haploid X; no chrY in the reference and no chrX variants in
  the female proband (chrY routing, PAR2, female XLR, sex-unknown skip never run); one CSQ block per
  record, no `&`-joined values, no `SpliceAI_pred_DP_*`; no ClinVar VCF (the transfer, its 0-match
  guard and CLNREVSTAT tuple parsing never run end-to-end). `run_integration.sh:9-10` says there is
  no gnomAD transfer (there is); CLAUDE.md names `mock_annotate.py`, which is `mock_vep.py`.
- **T8** Coverage holes with no test at all: the whole of Step 5's inheritance logic
  (`05_inheritance_screen.py:107-377` — origin mat/pat/both, obligate 1/1 transmission, 1/1×1/1
  Mendelian error, `consumed` only by phase-confirmed pairs, unphased de-novo partner, X-linked
  male/female, male-X de novo, chrY, PAR, sex-unknown, `require_pass`, hiConf gating); `genotype.py`
  hom_alt/hom_ref/clean_parent limbs and the AD fail-open asymmetry (`sample_qc_ad_measured` is never
  called by a test); `annotations.rarity_basis`/`gnomad_observed` with tuple witnesses; Step 6
  `main` counting (distinct individuals, comp-het pair → one carrier, `recurrence_kind`, de-novo BH
  padding); Step 9's `alpha<=0` error, `max_downweight_fraction` HALT, `min_control_genes` HALT,
  empty-input refusal, BH over the full universe. A Step-5 test needs only a stubbed `cyvcf2` module
  and FakeVar objects carrying `gt_types/gt_bases/gt_quals/gt_depths/gt_ref_depths/gt_alt_depths/INFO`.
- **T9** Tests write to fixed `/tmp/_hprv_*` paths (`:66,360,369,376`) — parallel runs collide;
  `:981-982` asserts `GENE_TIERS` against a literal copy of itself.

## Docs / config / code consistency (delegated sweep)

Every numeric default was cross-checked three ways (canonical table, `config.example.yaml`, code
fallback) and **all agree** — `filters.*`, `qc.*`, `inheritance.*`, `burden.*`, every
`prioritization.*` cut point, tier rule, penalty, ceiling, weight and cap, `resources.clinvar.*`,
`outputs.igv.*`, the three preflight halts, and the 68-test count. Every test/function/column
CLAUDE.md names exists. The discrepancies are all about *behaviour described*, not numbers:

- **D1 (see R9)** The retired per-variant oracle is the documented behaviour in ~40 locations; the
  live one-oracle contract is documented in ~6. CLAUDE.md carries both, back to back.
- **D2** Dead config keys presented as live (none read anywhere; grep-verified):
  `composite.gene_list_prior.points` (live key: `composite.weights.gene_list_prior`),
  `composite.gene_list_prior.combine_gene_and_set`, `signals.segdup.min_identity` (also pinned by
  `test_pure.py:2345` as a "code default"), `excess.offset.covariate_adjust`,
  `excess.offset.components`, `excess.refit_per_cohort`, `composite.emit_both_rankings`,
  `signals.segdup.track_build`, `scope.*` ("enforced as warnings" — nothing enforces). Other reserved
  keys carry a `[reserved]` marker; these do not.
- **D3** "VEP-only / no transfers / no stars / no REVEL" claims contradicted by code:
  `docs/pipeline_design.md:14-17,22,76-77,149,272-277,286`; `docs/functional_annotation.md:8-24,
  82-87,205-225,240-263`; `docs/clinical_classification.md:8-27,88-90,111-115,223-230` (one
  paragraph asserts stars exist and that a star gate "cannot be applied at all");
  `docs/prioritization.md:585-592` (tier table omits every REVEL/AlphaMissense route incl. V4);
  `docs/resources.md:13-36` vs `60-65`; `README.md:73-79` vs `23-31`.
- **D4** `prepare_resources.sh verify` (`:470-493`) reports gnomAD slim, REVEL and AlphaMissense as
  "extra; not required" while `run_pipeline.sh:124-136` halts without them by default.
- **D5** `docs/resources.md:427-429` says Step 2b is ON by default; config (`:138`), `run_pipeline.sh`
  and CLAUDE.md say off. `run_pipeline.sh:240-241` comment says both.
- **D6** Poisson anti-conservatism: 2.51× in `README.md:147`, `prioritize.py:523`, `09:637`; 2.41× in
  the canonical table (`docs/README.md:263`). `README.md:159-160` quotes the `min_n_for_ratio_rule: 1`
  tier table (228 genes / 10.4%) as if it were the shipped default (92 / 9.57%).
- **D7** `gene_is_constrained` docstring says `s_het` is not read; the body reads it. Step 9 never
  reads `phaplo_min`, Step 6 does — undocumented asymmetry between the two "constrained" predicates.
- **D8** `--max-sites` (R10) is undocumented; `docs/README.md:184` "FILTER = PASS only" while Step 1
  keeps `PASS,.`; UPDhmm is listed as performed QC (`inheritance_and_genotype_qc.md:34,207,273`) but
  nothing invokes it; `rarity_strength`'s `permissive_fail` band is undocumented;
  `docs/prioritization.md:724` names `moi_coherent|…` values that are actually `coherent|…` in a
  column named `moi_coherence`; REVEL 0.932 "strong" cut is documented as implemented but has no key
  and no effect; the absent-faf95 chr22 fraction is 74% in two places and 80% in five.
- **D9** Dangling references: `SCIENCE_AUDIT.md` (README, prioritization.md ×2),
  `{{artifact:…}}` image placeholders in `docs/prioritization.md:99,443,451,478`,
  `Analysis/.../prepare_igv_variants.py` (`prioritize.py:28,1199`), `mock_annotate.py` (CLAUDE.md).
  CLAUDE.md's duplicated MAX_AF bullets and the "6." numbered between "2." and "3." are cosmetic but
  signal that the file is being appended to rather than edited.

Coverage note: `test_prioritize_config_matches_canonical_defaults` pins Step-9 and ClinVar keys
only; no test compares Step 0/3/5/6 code defaults to the example config (and see T3).

## Shell layer (delegated sweep; S1, S2 and S5 reproduced against the real scripts)

- **S1 — stale VEP shards re-stamped as current · major.** `02_annotate_sites.sh:166` short-circuits a
  contig shard on a bare `.done`. The resource/union key (`:103-136`) guards only `$OUT`; when it
  mismatches, both the in-process loop (`:342`) and the SLURM scatter (`:317`) reuse every shard,
  `gather_shards` (`:182-197`) only checks that a `.done` exists per contig, and `:614` writes the NEW
  key into `$OUT.done`. Reproduced: add a site to the union → re-run logs `[chr1] cached (resume)`,
  the site is absent from the output, the marker carries the 4-site cksum, and a third run says
  "already complete". The same path means a newly supplied SpliceAI/CADD/REVEL/AlphaMissense file
  invalidates `$OUT` but is never lifted (the shards predate the plugin). This is the scenario Step 1's
  manifest key exists to close. Fix (S): stamp each shard `.done` with a VEP-input key (union cksum +
  plugin stats + cache path/version) and treat a mismatch as not-done.
- **S2 — faf95 can silently become "everything is rarest" · major.** `02:566-569` downgrades a
  non-zero `bcftools annotate` on the slim (unindexed, truncated, wrong build) to a WARN "continuing on
  the grpmax PROXY oracle" — but the oracle is config-selected (`annotations.py:290`) and Step 2 never
  reads it, so the run stays on `faf95` with no `gnomad_*` fields: `frequency()` returns `None` for
  every variant, `selection.py:70` passes everything including BA1-common alleles, `rarity_basis`
  reads `absent` on every row, the audit says `rarity_oracle.faf95`, exit 0. Independently,
  `run_pipeline.sh:105-126` puts the slim `_need` inside the non-ingest `else`, so
  `resources.vep.annotated_vcf` + default oracle + no slim runs the same way. Fix (S): export the
  oracle to Step 2 and `die` on a failed transfer under `faf95`; move the slim `_need` outside the
  branch; add a Step-3 header assertion (`gnomad_AF_joint` declared when oracle is faf95). The Step-5
  orphan counter would show the damage after the fact, as a WARN.
- **S3** `01:120` per-trio site files and `04:149` per-trio candidates resume on `trio_id` alone; a
  replaced VCF at the same path (or a `FILTER`/`EXCLUDE_CONTIGS` change) is never re-ingested. Fix:
  fold `stat` of the source VCF + filter settings into the per-trio `.done`, as Step 2 does for
  resources.
- **S4** `run_pipeline.sh:145-148` tests `-e` only and `02:405-409` tests the CSQ *header*; the plugin
  declares its keys whenever it loads, so a wrong-build/bad-tabix score file gives declared-but-empty
  `vep_SpliceAI_pred_DS_*` with no warning (`download_spliceai.sh:151-153` describes exactly this
  tabix failure). The frequency guard at `:575-601` checks values for this reason; SpliceAI/CADD have
  no equivalent. Fix: count non-missing values after split-vep and `die` on 0 when configured.
- **S5** `config.py:107` emits `export HPRV_VEP_SHARD_BY_CONTIG=True`; `02:335` matches
  `1|true|yes|on` → un-sharded single pass on the shipped config (the key absent → sharded). The
  `get` subcommand normalises bools; `emit_sh` does not. CI sets the env var directly
  (`assert_shard_equivalence.sh:40`) so it cannot see it. Fix: normalise in `emit_sh`.
- **S6** SLURM/CI: `08_igv_export.sh:357-363` places the "NHF unavailable" die inside
  `if is_set "$KRAKEN2_DB"`, so when `run_pipeline.sh:326-338` omits `--kraken2-db` the manifest is
  never written and `phase.sbatch:108` reads "no trios need NHF work" → green. `phase.sbatch:153-155`
  caps `DOWN_TO` at 8, so a SLURM run never produces `variants.prioritized.tsv`.
  `build-publish.yml:251-271` chains every probe as `cmd && echo OK` under `set -eu`, which cannot
  fail; only the final `grep -q` decides, so a conda-Perl-shadow regression would still publish
  `latest`.
- **S7** `02:103-131` key omits `HPRV_CSQ_SELECT`, cache path and VEP version (switching `pick`→`mane`
  reports "already complete"); `run_pipeline.sh:209-215` Step 0 cache is bare (a new trio gets no
  inferred sex and its X/Y modes are skipped with a WARN); `HPRV_TMPDIR` is not in `HPRV_BIND` for
  Steps 2 (`:141-149`) and 8 (`:147-159`) — in host-wrapped mode Step 8's slice/NHF fails per member
  with exit 0; `08:529` `--nhf-gather` audits `minicrams=0` and last-value-wins; `phase.sbatch:137`'s
  task-private tmpdir is stripped by `--cleanenv`; `prepare_resources.sh:392` (`prep_dbnsfp`, opt-in)
  still has the `zcat | head -n1` SIGPIPE trap and `scripts/` is not shellchecked; `prepare_resources.sh
  fetch/verify/emit-env` never prepare, and comment out, the slim the default oracle requires, so a
  fresh install "verified" by the script halts at preflight.
- Nits: `run_pipeline.sh:124-127` fires both `_need` and `_opt` for the slim with contradictory text;
  `--from/--to` unvalidated; `Dockerfile:191-202` smoke probes can't fail (`| head -1`, `|| true`);
  unpinned `micromamba latest`, LOFTEE branch tip, `action-shellcheck@master`;
  `08:48` puts CWD on `sys.path` when invoked directly; `prepare_resources.sh:263-265` records
  `ok clinvar` for the unrenamed fallback file.

Checked and found correct in the shell layer: `+split-vep` without `--threads`; `-s pick` gated on
bcftools ≥ 1.20; `annotate -c DST:=SRC` direction; `-x '^INFO/…'` keep-syntax; `concat -a -D -f`;
`sort -T` inside the tmpdir; `norm --keep-sum AD` header-gated; `isec -c none -n=2 -w1`; chrM
exclusion with a post-union count-and-die; VEP 115 option set and all four plugin strings;
`norm -m-` at Steps 1 and 4 against the same reference (so `--flag_pick` per record equals per
allele); explicit index removal before re-indexing; the SIGPIPE-guarded `head`/`grep -m1`
pipelines; SLURM afterok chain, gather re-verification, resubmit guard; the kraken2/nonhuman-screen
pins, conda-Perl removal and SpliceAI env isolation in the Dockerfile.

## Suggested order of work

1. R1, R4, R5, S5 — one-line fixes with immediate effect on every run's ranking or resumability.
2. S2, S1, S4 — the three silent-failure paths in Step 2; each is a few lines plus a test.
3. R3 — distinct-observation counting in Step 9, then re-derive the validation tier table.
4. R2 + R6 + R8 — Step 6's headline statistic and the workbook text; decide whether Step 6 ranks by
   an offset-normalised statistic or whether `genes.prioritized.tsv` becomes the headline.
5. R9/D1–D9 — collapse the docs to the one-oracle story and delete the dead config keys; CLAUDE.md
   first, since it is what every future editor reads.
6. T2–T8 — make the guard tests able to fail, and give Step 5 a unit-test harness.
