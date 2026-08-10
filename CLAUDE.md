# CLAUDE.md — implementer guide

Orientation for anyone (human or Claude) extending this repository. Read this first, then
[docs/pipeline_design.md](docs/pipeline_design.md) and
[docs/README.md#canonical-defaults](docs/README.md#canonical-defaults).

## What this is

A container-based pipeline that screens **GMKF Kids First per-trio VCFs** (GRCh38, GATK
genotype-refinement output, **not** jointly genotyped across the cohort) for high-priority
**inherited** rare variants, and consolidates **genes where rare functional variants recur across
individuals**, targeting **rare disease and germline pediatric cancer**. Runs under Apptainer on HPC.

**Scope / focus (important):** the emphasis is **inherited** germline variation —
**dominant** (a rare functional inherited het that recurs across individuals), **recessive**
(hom / compound-het-in-trans), and **X-linked**. **De novo** filtering/review and **mtDNA
heteroplasmy** are handled by **separate dedicated machinery** (the shared `.sh` orchestration and
a dedicated mtDNA pipeline). De novo is detected here only as a lightweight cross-reference
(`inheritance.emit_denovo`, default on but secondary); chrM is out of scope.

## Golden rules (do not violate)

1. **Public repo, controlled-access data.** Never commit VCF/BAM/CRAM/PED, real filesystem
   paths, sample/subject IDs, or PHI. All paths are `${ENV}` placeholders resolved at runtime.
   `.gitignore` enforces this; before every commit, sanity-check with
   `grep -rE '/Users/|/scratch|/home/[a-z]|BS_[A-Z0-9]{8}'`.
2. **gnomAD v4.1 is the ONLY population-frequency oracle.** These trios are not
   jointly genotyped, so internal cohort AC/AN is meaningless (absent ≠ hom-ref). Never `bcftools
   merge` the trios into a genotype matrix and never compute population frequency from internal
   counts; internal recurrence is valid only as an artifact/blocklist signal.
   The **fallback** rarity field is a **grpmax proxy** — max AF over the grpmax-*eligible* groups
   (`annotations.GRPMAX_POPS` = AFR/AMR/EAS/NFE/SAS) — read from the CSQ; the preferred field is
   real `faf95` (below). Two hard rules that bind BOTH arms:
   - **Never substitute VEP's `MAX_AF`.** It maxes over the bottlenecked founder groups grpmax
     deliberately excludes (ami AN≈900, asj, fin, mid) and over tiny 1000G populations; one allele
     there reads as AF≈1e-3 and silently kills dominant candidates at the 1e-4 gate.
   - **Never substitute the global AF.** It dilutes ancestry-enriched variants and fails the
     opposite way (retaining benign polymorphisms). The two wrong substitutions err in opposite
     directions — there is no single safe fallback.
   **TWO ORACLES, chosen per variant, and `annotations.frequency()` is the one chokepoint —
   never a field getter directly.** It prefers real **faf95** (`gnomad_faf95`, from the gnomAD
   v4.1 JOINT slim transferred in Step 2; `resources.gnomad.sites_slim`, opt-in ~10 GB) and falls
   back to the grpmax point-estimate proxy above wherever gnomAD published no faf95. Which one
   fired is reported per variant as `rarity_oracle`. Verified against the real v4.1 data: the FAF
   group set is afr/amr/eas/**mid**/nfe/sas — GRPMAX_POPS plus `mid`, EXCLUDING ami/asj/fin, so
   faf95 does not reintroduce the MAX_AF trap; `faf95_group` reports the producing group.
   **Absent faf95 splits in two, and the halves demand OPPOSITE fallbacks.** gnomAD emits fafmax
   as MISSING, never as 0, wherever no group's CI lower bound clears zero (80% of a chr22 sample).
   If gnomAD HAS the allele, faf95 IS 0 -> rarest, and the proxy must NOT be consulted: of the
   records with no faf95 but a proxy >= 1e-4, **96.5% are AC <= 2**, so using the point estimate
   there filters on one or two observed alleles — precisely what faf95 exists to prevent. Only
   when gnomAD has NO record is the proxy the right answer. `gnomad_AF_joint` is the witness that
   distinguishes them, which is why it is transferred. **Supplying the slim RETAINS MORE**: faf95 <= the point estimate, so
   the same cutoffs stop discarding low-count alleles the interval never justified discarding — a
   SMALLER list after enabling it means a broken join, not a better filter.
6. **VEP-centric contract.** The annotation surface is a VEP 115 GRCh38 cache + its score PLUGINS:
   CADD (required-ish) and **SpliceAI** (required by default). Nothing else is bcftools-transferred in: no
   gnomAD, ClinVar, dbNSFP or LOFTEE file. Adding an annotation means either a VEP plugin or a
   `bcftools annotate` transfer in Step 2 **and** its INFO field in `annotations.F` — never a lookup
   that reaches around that contract. **ClinVar is the one transfer** (`resources.clinvar.vcf` ->
   `clinvar_CLNREVSTAT` -> `clinvar_stars`), because the cache carries `CLIN_SIG` but no
   `CLNREVSTAT` at any price; it is prefixed `clinvar_` rather than `vep_` precisely so the
   different oracle is visible. **REVEL and AlphaMissense are plugins** (dedicated files, never
   dbNSFP — 32 GB for 5 columns and a dead URL). Remaining accepted losses (see
   [docs/README.md](docs/README.md#canonical-defaults)): no LOFTEE, and no exome/genome discordance
   flag. faf95 + nhomalt are RESTORED by the optional gnomAD joint slim (rule 2). Below MODERATE impact there are now TWO keep-paths: **SpliceAI** (max delta
   ≥ `spliceai_ds_min`, default 0.2 — the deep-intronic/synonymous splice signal; checked first) and
   **CADD** (≥ `cadd_phred_supporting`, default 25.3 — everything else non-coding). SpliceAI is a
   VEP plugin over the precomputed RAW scores (`resources.vep.spliceai_snv`/`spliceai_indel`); the
   precomputed set is not exhaustive, so a missing score is not "no effect" (keep-only, never drops).
3. **Gene lists and constraint are priors/tiers, never hard include/exclude** — the "never-drop
   rule" keeps novel-gene discovery alive. Rarity/impact/QC gating happens *before and
   independently of* any list.
4. **The canonical-defaults table is the single source of truth.** Every threshold in code and
   docs must match [docs/README.md#canonical-defaults](docs/README.md#canonical-defaults); all
   thresholds are config defaults, never hardcoded law. A gene-specific ClinGen VCEP value
   overrides a generic cutoff.
5. **Engineering ethos:** fail loudly, verify before claiming, be idempotent (`.done` files +
   integrity checks). Match the existing style.

## Architecture

- **One image** (`Dockerfile`): the group's validated `ensemblorg/ensembl-vep:release_115.0`
  base + a pinned micromamba env (`env/environment.yml`) for bcftools/bedtools/slivar/somalier/
  whatshap/python. The conda env must not shadow VEP's Perl — bioconda drags `perl` in, so the
  Dockerfile removes the conda Perl after solve (VEP's `env perl` must resolve to the base's Perl).
- **Config → env → scripts.** `config/config.example.yaml` is the contract. `src/hprv/config.py`
  loads it, expands `${ENV}`, and emits shell exports (`python -m hprv.config sh`). The bash
  steps take explicit args; `run_pipeline.sh` maps config → args.
- **Execution model:** the pipeline is meant to run **inside the container** (tools native on
  PATH). `pipeline/lib/common.sh`'s `hprv_run` auto-detects the runtime: inside an Apptainer
  container it resolves to `native` (direct calls); from a host it wraps each tool in
  apptainer/docker. So the same scripts work in both modes. `HPRV_BIND` carries the dirs that
  must be visible to wrapped calls.
- **Shared python** in `src/hprv/`: `config` (YAML+env), `annotations` (the INFO-field contract
  from Step 2), `genotype` (refined-GQ QC), `ped` (trio parsing). Steps import these so selection
  / inheritance / burden read annotations identically.

## Data contract between steps

- **User input** (git-ignored): a `trios_file` (TSV, header names kid/dad/mom in any order;
  IDs match VCF samples) + a `vcf_dir`/`vcf_list`. `pipeline/resolve_trios.py` maps each trio to
  the VCF containing all three members (exact match; picks the fewest-sample VCF on a tie; extras
  OK), generates PEDs, and writes the **internal manifest** `trios.resolved.tsv`
  (`trio_id  vcf  ped  samples`) that Steps 0/1/4 consume. Unresolved/ambiguous trios are reported
  in `trio_resolution.tsv`, never guessed. Steps 1 and 4 subset each VCF to its 3 members
  (`bcftools view -s`), so extra members and inconsistent sample order don't matter.
- **PED sex**: the generated PED leaves kid sex unknown (`0`); Step 5 reads Step 0's inferred sex
  (`qc_report.tsv`) so X-linked/hemizygous logic fires correctly.
- **Auditing**: every step calls `audit`/`hprv.audit.record` → `audit/counts.tsv`
  (timestamp, step, scope, metric, value; scope = `global` or trio_id). `python -m hprv.audit` assembles
  `audit/summary.md`. Step 3 tags kept variants with `hprv_keep_reason`.
- **VEP runs ONCE** on the cohort union (Step 2). Step 4 transfers annotations with
  `bcftools annotate` — it never re-runs VEP. Keep it that way. Already have a VEP VCF? Set
  `resources.vep.annotated_vcf` and Step 2 ingests it instead (verifying build + frequency
  presence); the rest of Step 2 is unchanged.
- **Step 2 INFO fields** (the contract `src/hprv/annotations.py` owns): all but two are CSQ
  fields lifted by `bcftools +split-vep` with a `vep_` prefix. The exceptions are
  `clinvar_*` (from the ClinVar sites VCF) and `gnomad_*` (`gnomad_faf95`, `gnomad_faf95_group`,
  `gnomad_nhomalt`, + two reporting-only AFs, from the gnomAD v4.1 joint slim) —
  `bcftools annotate` transfers under their own namespaces on purpose, so which oracle a field
  came from is readable at a glance. Those are the ONLY two transfers. `vep_Consequence`, `vep_IMPACT`, `vep_SYMBOL`, `vep_Gene`, `vep_Feature`,
  `vep_BIOTYPE`, `vep_HGVSc`, `vep_HGVSp`, `vep_MANE_SELECT`, `vep_CADD_PHRED`, `vep_CLIN_SIG`,
  `vep_SpliceAI_pred_DS_{AG,AL,DG,DL}` (+ `DP_*`, `SYMBOL`; `annotations.spliceai_ds()` = the max,
  the splice keep-path — present only when the SpliceAI plugin is configured),
  `vep_REVEL` + `vep_am_pathogenicity`/`vep_am_class` (the calibrated missense predictors —
  **inert at the screen**, consumed only by Step 9's missense tier),
  `vep_gnomAD{e,g}_{AFR,AMR,EAS,NFE,SAS}_AF` (the rarity oracle), plus `vep_MAX_AF` /
  `vep_MAX_AF_POPS` / `vep_gnomAD{e,g}_AF` for REPORTING ONLY — never as filter fields (see rule
  2). Add a new annotation by wiring it through Step 2's split-vep `want` list AND
  `annotations.F`.
- **`--pick` vs `--flag_pick`:** Step 2 runs VEP with `--flag_pick`, which keeps EVERY consequence
  block and marks the chosen one `PICK=1`, so split-vep's `-s` selector decides — and an
  externally-produced `--flag_pick` VCF takes the identical path. The selector auto-resolves to
  `pick` when the CSQ has a PICK field, else `worst`; override with `resources.vep.csq_select`.
  Watch out: a `--pick_order` starting with `rank` picks the WORST-consequence transcript, so
  `SYMBOL` can name a non-MANE/readthrough gene that Step 6 then aggregates carriers under.
- **`bcftools +split-vep` does not accept `--threads`** (it is a plugin; passing it aborts the
  step). Thread the VEP call with `--fork` instead.
- **Step 4 output**: per-trio `*.candidates.annotated.vcf.gz` + `trios.candidates.tsv`
  (`trio_id  candidates_vcf  ped`). Per-trio VCFs are the authoritative unit — no cohort genotype
  matrix is ever built.
- **Step 5 output**: `candidates.calls.tsv` (one row per candidate; `mode` ∈ `dominant`
  (inherited het; `flags=origin=mat|pat|both`), `hom_recessive`, `compound_het` (pairs share a
  `pair_id`; a pair whose second hit is a de novo is unphaseable from trio genotypes and carries
  `flags=unphased_denovo_partner`, and does NOT suppress the dominant call; a pair whose
  non-transmitting parent was never affirmatively observed hom-ref carries `origin_unverified`),
  `x_linked_recessive`, `denovo`/`denovo_x_hemi` (secondary)). A `1/1` parent transmits obligately,
  so parent-of-origin there is deterministic (`both` is reserved for HET×HET). chrY is routed away
  from the mother-keyed hemizygous models (`male_x_chrx`) and yields no inherited call. Step 5 opens
  VCFs with `strict_gt=True` — cyvcf2's default reports a half-called `0/.` as hom-ref, which would
  defeat every "parent is a confident no-call" test. Modes are configured in
  `inheritance.emit_dominant` / `inheritance.emit_denovo`. **Step 6 output**: `genes.ranked.tsv` —
  distinct-individual carrier counts per gene per model (`n_dominant`/`n_biallelic`/`n_xlinked`/
  `n_denovo`), `recurrent` flag (≥ `burden.min_carriers`), constraint columns; ranked
  recurrent-first, constraint-weighted.
- **Step 7 output**: `hprv_summary.xlsx` (openpyxl; `src/hprv/report.py`) — documented workbook:
  About/legend + Gene consolidation + Candidate calls + Trio resolution + QC + Audit counts.
- **Step 8 output**: `igv/` for the jlanej/igv.js variant-review server (`src/hprv/igv.py` +
  `08_igv_export.sh`): `variants.tsv` (only `chrom/pos/ref/alt` required; extra columns are
  filterable; per-member `*_file`/`*_index` + `*_vcf*` track paths are RELATIVE to the data-dir
  `igv/`), mini-CRAMs `crams/<trio>/<sample>.cram` sliced around candidate loci via a `sample→CRAM`
  map (`resources.cram_map`; `samtools view -C -T ref --regions-file bed`, ± `outputs.igv.padding`),
  per-trio VCF tracks `vcfs/<trio>.vcf.gz`, `trios.tsv`, `sample_qc.tsv`, empty `curation.json`, and
  `config.json` (`{"genome": ...}` from `outputs.igv.genome`, default `hg38` — the igv.js genome hint).
- **Step 8b (optional, default ON)**: non-human-fraction (NHF) annotation. If a kraken2 DB is
  provided (`resources.kraken2_db`, bind-mounted DATA — never baked), `nonhuman-screen classify`
  runs over each screened member's **mini-CRAM** (no new source-CRAM I/O) against the per-trio VCF
  → `igv/nhf/<trio>/<sample>.variant_nhf.tsv`, folded into `variants.tsv` as `child_/mother_/father_nhf`
  (+ `_reads` denominator) and a derived `nhf_flag`. Configured under `outputs.igv.nonhuman_screen`
  (`members: carriers|child_only|all`, `confidence`, `min_reads`, `memory_mapping`). Gated: activates
  only with a DB + mini-CRAMs + the `nonhuman-screen` binary (image-only) — else warns and leaves the
  NHF columns blank. **The join is on the 0-based key (`pos-1`)** — nonhuman-screen keys variants
  0-based (`{chrom}:{pos0}:{ref}:{alt}`) while `variants.tsv` `pos` is 1-based; `igv._load_nhf_tsv` +
  the `pos-1` join in `build_variants_tsv` is the single load-bearing off-by-one.
- **Step 9 output**: `variants.prioritized.tsv` + `genes.prioritized.tsv`, plus — when the input was
  Step 8's table — **`igv/variants.prioritized.tsv`** (`--out-igv-variants`), the table a reviewer
  actually opens: every input column verbatim in its original order **plus** every prioritization
  column appended. It exists because `variants.prioritized.tsv` is NOT a drop-in for
  `igv/variants.tsv` — it sits outside `igv/` and its column set omits the `*_file`/`*_index`/
  `*_vcf*` track paths, which are RELATIVE to the `igv/` data dir, so the review server would
  render a sortable list with no mini-CRAMs and no VCF tracks. The appended set is the strict
  COMPLEMENT of the input header (input values are never overwritten), rows are carried by
  POSITION rather than re-joined on `chrom/pos/ref/alt/trio_id` (ambiguous for two ALTs of one
  multiallelic site in one trio), never-drop is asserted a second time on this file, and it is
  sorted by `rank_agnostic` so it opens honest with an overlay loaded. (`src/hprv/prioritize.py`
  = ALL pure logic, no I/O; `09_prioritize.py` = the CLI). **Every merged source table is ALSO
  emitted verbatim** — one `src_<source>_<column>` per source column (`src_mutrate_` = the 76
  gnomAD constraint columns, plus `src_constraint_`/`src_segdup_`/`src_moi_`/`src_burden_`/
  `src_prior_`) —
  appended AFTER the curated set so leading columns are unchanged. Two load-bearing rules: the
  prefix is **per source**, because two tables legitimately carry the same column name (gnomAD's
  `pLI` and a LOEUF-only table's `pLI`) and one shared `src_` would silently drop one; and a gene
  ABSENT from a source yields `''` (MISSING), never `0.0`. The bare curated columns are what the
  SCORING read (possibly imputed or chosen between tables) and are NOT interchangeable with their
  `src_*` counterparts — a disagreement means a fallback fired and must stay visible.
  `--no-source-columns` suppresses the block and is part of the idempotency key. For
  `src_prior_` specifically: overlay entries are keyed UPPER-cased internally so the pass-through
  registers each gene under both casings (a case-only mismatch would otherwise blank the raw
  cells while the scoring join still fired), and a gene contributed by SET membership alone has
  no source row — blank `src_prior_*`, with `gene_list_prior_set_applied` naming the set. **Input precedence: `igv/variants.tsv`
  (Step 8) when it exists, else `candidates.calls.tsv` (Step 5)** — Step 8's table is the only one
  carrying the NHF columns the quality term reads, and with Step 5's table every call correctly
  reads `nhf_status=not_screened`. Optional inputs each degrade with a loud WARN exactly as Step 6
  does for `--constraint`: `--mutrate` (no excess statistic → every gene reads `T0`, the variant
  layer still runs), `--constraint`, `--segdup`, `--established-genes`, `--gene-moi`, `--gene-prior`.
  Idempotent via `variants.prioritized.tsv.done`; `--force` re-runs.
  **NEVER-DROP IS AN ASSERTED INVARIANT HERE**: `09_prioritize.py` checks row-count conservation
  before writing, and the down-weight sets a tier + a separately-reported additive penalty + a
  human-readable `downweight_reason` (plus `review_flag`), never a filter. The maximum penalty (−3.0) is sized so it
  cannot alone demote a variant with strong molecular evidence in a constrained gene.
  Two rankings always ship — `rank_agnostic` (no gene-list prior of any kind), `rank_prior`, and
  `rank_delta`. `prioritization.composite.gene_list_prior.enabled` defaults **false** and a config
  overlay path is inert while it is false, so hprv stays phenotype-agnostic by default.

## Gotchas that WILL bite you

- **Step 9's mutational-target table is NOT `resources.constraint.gnomad_v2_constraint`, and it
  arrives BGZIPPED.** Same gnomAD v2.1.1 download, two different prepared artifacts:
  `join_constraint.py` projects the file down to `gene/oe_lof_upper/pli/s_het/phaplo`, which DROPS
  every column the excess statistic needs (`mu_mis`/`mu_syn`/`mu_lof` = the offset itself, plus
  `oe_syn`/`classic_caf`/`constraint_flag` = three of the six artifact signals, plus `cds_length` =
  the fallback offset). So `prioritization.resources.mutational_target` points at the **unjoined**
  `constraint/mutational_target.by_gene.txt.bgz`. Because that path is `.bgz`, any reader in Step 9
  must go through `_open_text()` (the gzip branch) — a plain `open()` there dies with
  `UnicodeDecodeError` on byte 2, and the integration mock uses a plain `.tsv` so it would never
  notice. `tests/test_pure.py:test_prioritize_reads_bgzipped_tables` is the regression guard.
- **gnomAD v2.1.1 constraint coordinates are GRCh37/hg19** (verified: `BRCA1`
  chr17:41,196,312–41,277,500), even though the rest of this pipeline is GRCh38-only. Any
  interval overlap against that table — the segdup signal above all — must use an hg19 track or
  lift over first. **Mixing builds yields ~0 overlap SILENTLY**, so the signal simply never fires
  and the corroboration count is quietly measured on five signals instead of six.
- **Symbol reconciliation against the constraint table must go through Ensembl gene IDs, never
  symbol aliases.** Alias matching produced demonstrably wrong joins on the real data: `ACOD1` →
  `CAD`, `DRC3` → `EPS8L1`, `EMSY` → `TNRC6A`, `TRDC` → `BCL11B`. Join on the versionless Ensembl
  gene ID and discard any mapping whose target symbol already carries its own count, or you
  double-count.
- **The excess statistic's `C` must be fit over the FULL gene universe, zero-count genes included.**
  The candidate list is a *zero-truncated* sample: fitting on matched genes only multiplies `C` by
  `n_universe / n_matched` and divides every excess ratio by the same factor — the direction that
  HIDES artifact loci. 9,266 of 19,643 genes had zero counts on the validation cohort. Likewise the
  BH-FDR runs over the full universe (zeros at p = 1); correcting over called genes only makes `q`
  anti-conservative.
- **The NB tail must be the regularized incomplete beta, not a `1 - cdf` complement sum.** The
  artifact tail reaches ~1e-241 (`OR4Q3`: n = 229 against E = 0.42). A complement sum returns
  ~1.1e-16 there — the machine epsilon left over from summing to 1.0, wrong by 224 orders of
  magnitude — collapsing every extreme locus into one indistinguishable bin.
- **A high `excess_ratio` is a QUALITY signal, never a biology signal.** It means "this gene emits
  more candidate rows than its mutational target predicts" — mismapping, paralogue collapse, a
  founder allele, a callability defect, an ancestry-uneven rarity gate — and **not** "this gene is
  disease-associated". 13 established predisposition genes were extreme excess outliers on the
  validation cohort (`CTSA` at 148×). A gene can be both; the response is read-level review, not
  discarding the gene. Never phrase an excess result as evidence for or against a gene being real.
- **`sig_caf_low`'s direction is counter-intuitive and it is not a frequency signal.** High-excess
  genes have far LOWER cumulative pLoF allele frequency, and the *high* decile is actually depleted.
  A near-zero `classic_caf` means gnomAD reports essentially no pLoF alleles there, which makes both
  `mu_g` (our denominator) and the grpmax rarity oracle (our numerator's filter) unreliable. It is a
  **low-information-locus** flag. Do not "fix" it by inverting the comparison.
- **A RATIO rule without a count floor is not a rule.** `excess.min_n_for_ratio_rule` (default **3**)
  guards T2 and both RATIO limbs of T3 — never the FDR limb, which needs no guard because no gene
  with n<5 reached q_nb<0.05. Phase 1 specified this guard and then never applied it: without it,
  136 of 228 triaged genes had n<3 and carried just 7.6% of triaged volume, and 8 of the 13
  ceiling-exempted control genes had **q_nb = 1** — no statistical evidence of excess at all (SMPX,
  HMGA2, HAMP, PET100 were each ONE variant against E≈0.16). They needed a count floor, not an
  exemption. Set to 1 to reproduce the original tier table; both are in the docs.
- **The trim loop belongs to WHICHEVER null is being fit.** `calibrate_null` gives the Poisson arm
  its own `fit_poisson_null` and reports the C it used as `calibration.C_poisson`. Evaluating the
  Poisson at the NB's trimmed C measures a hybrid nobody would deploy — and it changes the headline
  figure (2.41× anti-conservative trimmed vs **1.82×** untrimmed). Never quote the Poisson
  anti-conservatism without saying which reading produced it. Also: report the **converged** trim
  count (77 genes / 0.392%), not an intermediate iteration's (69 / 0.35%) — the trim fraction feeds
  a hard HALT guard.
- **THE OVERLAY JOIN IS ON GENE SYMBOL ONLY — never route it by MOI.** A gene's prior must be
  expressible independently of its canonical mode of inheritance, because a prior can encode a
  hypothesis about a *different* genetic model than the gene is curated under. The concrete case:
  the FA/HR genes (FANCA, FANCD2, SLX4, FANCE, BRCA2) carry germ-cell-tumour evidence about
  **HETEROZYGOUS carriers** (PMID 40906985, five-gene combined OR 10.17) while their canonical model
  — and their PanelApp green status — is biallelic Fanconi anemia. An MOI-routed lookup consults
  those rows under a recessive model and silently misses the het observations the evidence is about.
  Symmetrically, an observed-mode-vs-canonical mismatch emits
  `moi_caveat = moi_mismatch_het_in_recessive_gene` and charges **NO penalty** — never-drop applied
  to the coherence layer.
- **Never SUM a gene-level and a set-level overlay prior — take the MAX.** A curated overlay can
  carry per-gene rows AND a pathway-collapsed set entry derived from the SAME study (the GCT
  resource's `FA_HR_PATHWAY_23`, pooled OR 4.14, and its per-gene FA rows are both PMID 40906985);
  summing double-counts one study. Enforced in `parse_gene_prior_overlay`, asserted in a test.
- **A `prior_weight` of 0.15 on a somatic driver is NOT weak germline support.** Curated overlays
  list somatic drivers (KRAS, NRAS, CBL, MTOR, AKT1, BCORL1) deliberately, so a reader can see they
  were considered and excluded. Any row whose `evidence_class` is in
  `composite.gene_list_prior.non_germline_classes` contributes **0.0** and carries
  `gene_list_prior_excluded_non_germline` — reported, never silently dropped. And `prior_weight` is
  **UNCALIBRATED** everywhere: an ordering default, never a likelihood ratio or an odds ratio.
- **The `--gene-prior` overlay reader must sniff COMMA as well as tab.** A phenotype panel handed
  over by a collaborator is usually a spreadsheet export, and a tab-only gate routed it to the
  bare-symbol-list branch: every line became one "symbol" (`BRCA1,0.9,GREEN`) that can never match
  a gene, while the reader still announced `3 genes from panel.csv` and the run went fully
  phenotype-agnostic. A plausible gene COUNT over a list that matched NOTHING is the failure nobody
  catches — the reviewer ships an un-prioritised list believing their panel applied. So: comma is
  accepted, and a HEADERLESS table is a **hard stop** (`return 1`), deliberately unlike every other
  optional resource, which degrades with a WARN because its absence is honestly reportable. Guard:
  `test_prioritize_gene_prior_overlay_accepts_csv_and_rejects_headerless_table`. Activation recipe:
  [docs/prioritization.md#13b](docs/prioritization.md).
- **A literal `gene` header row will be read as a gene symbol if you let it.** The validation
  cohort's own per-gene counts file has one, CARRYING n=1 — which is why the true totals are
  **25,389 variants / 10,799 symbols**, not the 25,390 / 10,800 raw line count. `GENE_KEYS_LOWER`
  guards every symbol-reading path. A phantom gene with a real count survives review because every
  individual number still looks plausible.
- **`priority_points` is NOT an ACMG score.** Never read a total against Tavtigian's P ≥ 10 /
  LP 6–9 / VUS 0–5 bands and never emit a P/LP/VUS label from it: the criteria are not ACMG criteria
  (a CADD-based term is not PP3), no phenotype/segregation/functional evidence exists at all, the
  ClinVar term applies only an uncalibrated review-status damp (positive limb only, not ACMG
  PP5/BP6), and the artifact-penalty terms have no ACMG analogue. The
  column is named `priority_points`, never `acmg_points`.
- **Three acquisition traps, each of which fails SILENTLY rather than loudly.** All are handled in
  `scripts/prepare_resources.sh`; the point is that none of them errors — each yields a resource
  that opens fine and matches nothing.
  1. **ClinVar ships BARE contig names** (`1`, not `chr1`) while these callsets are chr-prefixed.
     `bcftools annotate` matches on the contig STRING, so an unrenamed ClinVar transfers ZERO
     records and exits 0. `prep_clinvar` renames and then PROVES the rename landed; Step 2 also
     hard-fails on a 0-match transfer (a cohort union always overlaps ClinVar somewhere).
  2. **REVEL must be sorted on the GRCh38 column.** Its columns are
     `chr,hg19_pos,grch38_pos,...` — GRCh38 is column **3**, hence `sort -k1,1 -k3,3n` and
     `tabix -s 1 -b 3 -e 3`. Sorted/indexed on column 2 (hg19) it indexes without complaint and
     matches nothing. Its published upstream recipe also contains the `zcat | head -n1` SIGPIPE
     trap documented below — `prep_revel` guards it, and the ~8 GB sort gets an explicit `-T`.
  3. **AlphaMissense's field is `am_pathogenicity`, not `AlphaMissense_score`.** The latter is
     dbNSFP's name for the same quantity; putting it in Step 2's split-vep `want` list gives you a
     plugin that runs and a column that is never populated. `tests/integration/mock_vep.py`
     deliberately uses the plugin spelling so an upstream rename breaks the mock.
- **`clinvar_stars` blank is NOT zero, and stars RANK rather than gate.** Blank = the ClinVar
  transfer did not run (nobody looked); `0` = ClinVar has a record whose submitter provided no
  assertion criteria. Conflating them silently damps every P/LP assertion in a run with no ClinVar
  resource, so an absent star count leaves Step 9's clinical term at FULL weight. Only the
  POSITIVE limb is damped — shrinking a low-star BENIGN term toward zero would promote a
  poorly-reviewed benign call. And the screen stays star-blind: reinstating the old >=2-star
  keep/drop gate would violate never-drop.
- **Absent `faf95` means TWO different things, and conflating them is a real defect (this bit us).**
  gnomAD emits `fafmax` as MISSING, never as `0`, wherever no ancestry group's 95% CI lower bound
  clears zero — 80% of a chr22 sample. So:
  - **gnomAD HAS the allele** -> faf95 IS 0 -> **rarest, keep**. Do NOT consult the proxy. Of the
    records with no faf95 but a proxy >= 1e-4, **96.5% are AC <= 2** — gnomAD singletons, whose
    point estimate is inflated by a small group AN. Filtering those on the point estimate is
    exactly the error faf95 exists to prevent (Whiffin 2017). The first implementation did this
    and would have dropped ~8% of gnomAD-observed alleles at the dominant gate.
  - **gnomAD has NO record** -> the cache proxy is the only estimate available; use it.
  `gnomad_AF_joint` is the WITNESS that separates them — the reason a "reporting only" field is
  actually load-bearing. `rarity_oracle` reports `faf95` / `faf95_zero` / `grpmax_proxy` /
  `absent` PER VARIANT; a run-level label would be wrong, because all of them occur in one run.
  Step 5 resolves this ONCE (`rarity_af` + `rarity_oracle`) and Step 9 consumes it rather than
  re-deriving from the raw columns — a local re-derivation cannot tell the two cases apart.
- **Enabling the gnomAD slim makes the candidate list BIGGER.** faf95 <= the point estimate, so
  the same cutoffs stop discarding low-count alleles whose CI never justified the call. If the
  list got SMALLER after supplying it, the join is broken — check Step 2's
  "gnomAD joint matched N / M sites" line, which is guarded to die at 0 but not at 1.
- **`MAX_AF` is a trap, not a shortcut.** It is right there in the CSQ and looks like the rarity
  field. It is not — see golden rule 2. It maxes over founder groups (ami AN≈900) and 1000G
  populations that gnomAD's grpmax excludes on purpose, so a single allele reads as AF≈1e-3 and
  silently drops dominant candidates. Rarity comes from `annotations.frequency()` ONLY.
  `tests/test_pure.py:test_frequency_excludes_bottlenecked_pops` and the GENEFND integration case
  exist to catch a regression here.
- **Absence of a cached AF is weak evidence.** VEP caches frequencies only for alleles
  accessioned into **dbSNP** — an un-accessioned gnomAD variant returns no AF and reads as
  "absent ⇒ rarest". Ensembl itself recommends `--custom` with the gnomAD VCF over `--af_gnomad*`
  for this reason. Under the VEP-only contract we accept it; it biases toward retention.
- **REVEL/AlphaMissense are wired but STILL dead at the screen, and that is structural.** They are
  missense-only; every missense is `IMPACT=MODERATE`; `selection.py` keeps it at an earlier branch
  and returns before any predictor is consulted — true whether or not the resource is configured.
  CI asserts these keep-reasons never fire. Their only consumer is **Step 9's missense tier**,
  where the ladder is `revel -> alphamissense -> cadd(off-label) -> none`: a FIXED PRECEDENCE, not
  a max over what is available, because ClinGen SVI says commit to one predictor chosen before
  seeing results and best-of-N is an uncalibrated cherry-pick. `missense_evidence_source` always
  names which one spoke. Do not "improve" this by taking the maximum.
- **LOFTEE plugin code is baked into the image** at `/plugins` (the Dockerfile clones the
  `konradjk/loftee` **grch38** branch there — the base image ships all other VEP_plugins but
  `--skip_plugins LoF`, and master LOFTEE is GRCh37-only). The code ships; the **data**
  (human_ancestor/GERP/loftee.sql) is not fetched and LOFTEE is **not invoked** under the
  VEP-only contract. Kept in the image so re-enabling it is a config change, not a rebuild.
- **`hiConfDeNovo` may be absent.** The Kids First genotype-refinement workflow may skip
  `VariantAnnotator PossibleDeNovo`. Step 5 requires the tag only when it is present in the
  callset header; otherwise it detects de novos from genotypes + QC. Don't assume the tag exists.
- **gnomAD-prior suppression (real failure mode).** `CalculateGenotypePosteriors` uses gnomAD
  priors that can push a genuine ultra-rare pathogenic call toward hom-ref. Step 5 flags de novo
  candidates (`review_prior_crosscheck`) — for top hits, cross-check the pre-refinement `PL`/`GT`.
- **VEP cache release must match** the VEP binary (115) and be bind-mounted (never baked into the
  image).
- **Step 8b's NHF join is 0-based; the DB has two moving parts.** nonhuman-screen keys variants
  0-based (`{chrom}:{pos0}:{ref}:{alt}`), so `igv.build_variants_tsv` joins on **`pos-1`** — off by
  one and every NHF lands on the neighbouring row (`test_igv_nhf_join_is_pos_minus_one` guards it
  with a decoy key). The **kraken2 binary is source-built** in the Dockerfile at a pinned version
  whose Perl wrappers must run under VEP's Perl (the conda env is kept perl-free; the build's
  `kraken2 --version` smoke check enforces it). **nonhuman-screen is pip-pinned to a COMMIT**, not
  PyPI's latest — bumping either ref is a contract change: re-verify the CLI flags, the
  `variant_nhf.tsv` column order, and the 0-based key before merging.
- **SpliceAI live backfill (Step 2b) runs in an ISOLATED conda env.** `spliceai` needs TensorFlow,
  whose numpy/protobuf pins conflict with the main hprv env (cyvcf2 needs numpy≥2). So the Dockerfile
  builds a separate `/opt/conda/envs/spliceai` and `02b_spliceai_backfill.sh` invokes it via
  `micromamba run -n spliceai spliceai …` — never mix it into the main env. The model weights + GENCODE
  annotation are BUNDLED in the package (no data download; the build smoke-test hard-fails if the model
  won't load). Step 2b is OFF by default (`resources.vep.spliceai_backfill.enabled: false`) — the screen
  currently runs on the PRECOMPUTED scores alone. Set it `true` to also score the unscored gap; an
  absent isolated env then HALTS at the Step-2 preflight rather than silently skipping. It runs
  after Step 2 and before Step 3 (so a backfilled score is a keep-path), scores only variants with no
  precomputed value (default indels), and folds them into the same `vep_SpliceAI_pred_DS_*` fields
  (`bcftools annotate`). Idempotent via a `.spliceai_backfill.done` marker that must be NEWER than the
  union; a scoring failure warns + degrades to precomputed-only rather than aborting the run.
- **`$(producer | head ...)` under `set -euo pipefail` is a silent-abort trap.** `head` (and
  `grep -q`, `grep -m1`) exits as soon as it has what it needs and closes the pipe, so the
  producer dies of SIGPIPE = **exit 141**; `pipefail` makes that the pipeline's status and
  `set -e` aborts the assignment — before any guard can report it. It fires when the data is
  HEALTHY (enough output to make `head` exit early), so it reads as "works on toy input, dies on
  real input, no error message". `scripts/download_spliceai.sh` hit exactly this and could
  essentially never reach its success path. Two shapes, two fixes: put `|| true` INSIDE the
  `$( )` after the pipeline (the surrounding `[[ ]]` guard still catches genuinely bad data), and
  never pipe into `grep -q` when the pipeline's status is the answer — it fails *silently wrong*
  rather than aborting (a chr-prefixed reference was reported as `nochr`). Use `grep -cx` +
  a count test, which reads to EOF and cannot SIGPIPE the producer.
- **Apptainer:** point `APPTAINER_TMPDIR`/`CACHEDIR` at real disk and **do not use
  `--containall`** — a tmpfs `/tmp` OOMs heavy VEP/sort (documented failure in the group's
  original annotate script; `common.sh` already sets a disk-backed workdir).
- **Reproducibility hardening:** `environment.yml` pins versions; for byte-identical rebuilds run
  `conda-lock` and point the Dockerfile at the lockfile. Pull the image by `@sha256:` digest.

## Running

- **Resources (one-time, small):** the core surface is two things on the host — a **VEP 115 GRCh38
  cache** (~24 GB; it carries gnomAD v4.1 frequencies + ClinVar) and the **CADD** SNV+indel files —
  plus the **SpliceAI** raw score files (`resources.vep.spliceai_snv`/`spliceai_indel`), **required
  by default**: when Step 2 invokes VEP, `resources.vep.spliceai_required: true` HALTS the run at
  preflight if they are missing; set it `false` to degrade with a warning. Not enforced when
  `resources.vep.annotated_vcf` is set. The full raw set is a one-time BaseSpace (login) download
  (`scripts/download_spliceai.sh`). Nothing else: no gnomAD, ClinVar, dbNSFP or LOFTEE download.
  `prepare_resources.sh --dir DIR fetch|verify|emit-env` still fetches these, runs INSIDE the image
  (baked in at `/opt/hprv/scripts`, on PATH, with its pinned `/opt/hprv/resources/manifest.env`;
  `HPRV_RESOURCE_MANIFEST` re-pins without a rebuild), and `emit-env` writes the `${ENV}` exports
  the config expects. Never bake resource DATA into the image (`.dockerignore` keeps `resources/*`
  except the manifest out of the build context). See [docs/resources.md](docs/resources.md).
  Already have a VEP VCF? Skip all of it: set `resources.vep.annotated_vcf` and Step 2 ingests it.
  **Optional (Step 8b only):** a **kraken2 database** (`resources.kraken2_db`, e.g. PrackenDB
  `k2_NCBI_reference`, ~tens of GB, must ship `taxonomy/nodes.dmp`+`names.dmp`) enables NHF
  screening. Same rule — bind-mounted DATA, never baked; put it on **local NVMe** so the
  `--memory-mapping` page cache stays warm across the serial per-trio invocations. See
  [docs/resources.md](docs/resources.md#kraken2-database-optional-step-8b-nhf).
- **HPC (primary):** `apptainer exec --cleanenv --bind ... hprv.sif run_pipeline.sh --config
  config/config.yaml` (add `--from N --to M` for a subset).
- **Dev/host:** run individual step scripts; python steps need `PYTHONPATH=src` and the container
  env (cyvcf2/pysam/scipy/pyyaml) — easiest is to exec them inside the image.

## Testing

- **Host, no heavy deps:** `python3 -m py_compile` all scripts; `bash -n` all shell;
  `python3 tests/test_pure.py` (pure-logic: config, ped, trios-file parsing, annotation getters,
  genotype QC, selection funnel, Step-6 helpers, and the Step-9 prioritization layer — the NB
  fit/tail/BH-FDR, the never-drop invariant end-to-end through the CLI, the positive-control guard,
  both tier ceilings, blank-vs-zero NHF, mechanism gating, and a check that every default in the
  code equals `config.example.yaml`'s value). **59 tests, no network and no VCF.**
  **One documented exception to "no heavy deps":** the 6 tests that drive `09_prioritize.py:main()`
  need `yaml` transitively (`load_config` does `import yaml`). They declare it at the `_load_p9()`
  chokepoint and **SKIP** without it — and `_run_all` then refuses to print "All N passed", instead
  reporting `48 passed, 6 SKIPPED ... NOT full coverage`, because the skipped set holds the
  never-drop and cache-invalidation guards. **CI `pip install pyyaml`s** so they actually execute
  there rather than being permanently green-by-skipping. Before calling this suite green, run it the
  way CI does — a BARE `python3`, not an env that happens to carry the container's packages. (This
  bit us: the Step-9 CLI tests were reported passing from a local env with pyyaml, and CI died with
  `ModuleNotFoundError: No module named 'yaml'` partway through the run.)
- **End-to-end integration** (`tests/integration/`): `run_integration.sh` generates a tiny
  self-consistent mock genome + trios (`make_mock_data.py`) engineered to exercise every mode
  and filter path, runs resolve + Steps 0,1,3,4,5,6,8,9 with REAL bcftools + the python steps (only
  Step 2's VEP call is mocked via `mock_annotate.py`), and asserts the resolution, funnel, and
  calls (`assert_integration.py`). The Step-9 fixtures are engineered so the artifact locus and the
  established-gene positive control have the **same** extreme excess shape — the only thing
  separating them is the control ceiling, which is exactly what the assertion tests. Two mock-scale
  config deviations are deliberate and commented in `make_mock_data.py`: `min_control_genes` 1000→3
  and `max_downweight_fraction` 0.20→1.0, both because the mock has a handful of genes rather than
  an exome. Runs in CI on host bcftools — no image build needed. To run
  locally you need bcftools/samtools/bgzip/tabix + a python with cyvcf2/pysam/scipy/pyyaml on PATH.
- **Validation (TODO):** GIAB/CMRG truth sets + a positive-control variant panel to measure
  sensitivity/precision of the inheritance-model and recurrence logic on real data.

## Out of scope here (handled by separate machinery)

- **De novo** filtering/review — bespoke machinery (the shared `.sh` orchestration). De novo is a
  secondary cross-reference here (`inheritance.emit_denovo`), never the driver.
- **mtDNA heteroplasmy** — a dedicated pipeline; chrM is not analyzed. **Enforced**, not merely
  declared: Step 1 drops `EXCLUDE_CONTIGS` (`chrM,chrMT,M,MT`) from the cohort union and dies if
  any survive. It has to be enforced there rather than left un-modelled, because every inheritance
  mode is diploid and `genotype.py`'s sex predicates know only X/Y — a chrM record would route
  through `hom_recessive`/`dominant`/`compound_het` as if autosomal, and rCRS-referenced
  haplogroup variants (hom-alt in the whole trio, no gnomAD mito AF ⇒ rarity gate passes) would
  flood the recurrent tier of `genes.ranked.tsv`.

## Open TODOs

See **[docs/ROADMAP.md](docs/ROADMAP.md)** for the prioritized, dependency-ordered gap list from the
SOTA review. Done: calibrated recurrence null + FDR (Step 6), CHARR/freemix contamination gate
(Step 0). Quick wins still open: somalier ancestry/relatedness, PP1/BS4 co-segregation,
UTRannotator, UPD rescue, conda-lock. Big bets: germline CNV (GATK-gCNV), phenotype ranker
(Exomiser), read-backed/population phasing, ROH.

Lower-level items:

- **CNV/SV module** (GATK-gCNV / Manta / ExomeDepth) — the biggest coverage gap.
- **Pseudogene/seg-dup handling** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) — flag/annotate.
- **Dominant-recurrence corroboration** — add TRAPD (case-vs-gnomAD carrier frequency) as an
  optional external-control check on recurrent-gene nominations.
- **Phenotype layer** — Exomiser/LIRICAL HPO ranking as a *prior* (never-drop), plus HPO ingestion
  that degrades gracefully when phenotype is sparse.
- **Pediatric-cancer overlay** — implement the ACMG SF v3.3 / PanelApp-green tiering and
  second-hit boost as a reporting overlay.
- **conda-lock + digest-pinned base image**; **real-data validation** (GIAB/CMRG + positive controls).

## Commit conventions

Small, logical commits (`chore:` / `docs:` / `build:` / `feat:`). End messages with
`Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`. Work on a branch; open a PR into
`main`.
