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
   counts; internal recurrence is valid only as an artifact/blocklist signal. Every rarity decision
   reads `annotations.frequency()` — the single chokepoint — never a field getter directly. Two hard
   rules bind BOTH of its arms (rule 3):
   - **Never substitute VEP's `MAX_AF`.** It maxes over the bottlenecked founder groups grpmax
     deliberately excludes (ami AN≈900, asj, fin, mid) and over tiny 1000G populations; one allele
     there reads as AF≈1e-3 and silently kills dominant candidates at the 1e-4 gate.
   - **Never substitute the global AF.** It dilutes ancestry-enriched variants and fails the
     opposite way (retaining benign polymorphisms). The two wrong substitutions err in opposite
     directions — there is no single safe fallback. Both are carried for REPORTING ONLY, and
     neither is a gate field under EITHER arm.
3. **ONE frequency oracle per run — the arms NEVER cross.** `resources.gnomad.oracle` selects the
   quantity for the whole run, and `annotations.frequency()` never consults the other arm:
   - **`faf95` — the DEFAULT.** gnomAD's published filtering allele frequency
     (`fafmax_faf95_max_joint`: the lower bound of the 95% Poisson CI, the quantity ACMG/ClinGen
     specify — Whiffin 2017), transferred in Step 2 from the gnomAD v4.1 JOINT slim as
     `gnomad_faf95`. It REQUIRES `resources.gnomad.sites_slim`: `run_pipeline.sh` HALTS at
     preflight without it, and Step 2 dies on a missing or FAILED transfer (including on the
     `resources.vep.annotated_vcf` ingest path) — the same contract as `spliceai_required`.
     Within this arm, absent faf95 splits in two and the halves resolve DIFFERENTLY: gnomAD emits
     `fafmax` as MISSING, never as 0, wherever no group's CI lower bound clears zero (roughly 80%
     of a chr22 sample). An allele gnomAD HAS but published no faf95 for resolves to **0.0** —
     rarest, `rarity_basis=zero_ci` (96.5% of that class are AC ≤ 2, and filtering a singleton
     on its point estimate is precisely the error faf95 exists to prevent). An allele with NO
     gnomAD record resolves to **None** — rarest, `rarity_basis=absent`. The proxy is NEVER
     consulted under faf95; `gnomad_AF_joint` is the witness that separates the two halves, which
     is why a "reporting only" field is transferred at all.
   - **`grpmax_proxy` — a deliberate opt-down** for a run without the slim: the max VEP-cache
     POINT-ESTIMATE AF over the grpmax-eligible groups (`annotations.GRPMAX_POPS` =
     AFR/AMR/EAS/NFE/SAS), read from the CSQ. The cache carries no AC/AN, so no CI correction is
     possible on this arm; it sits ~one CI-width high on low-count alleles and errs toward
     DROPPING.
   `rarity_oracle` is a RUN-level constant, recorded once in `audit/counts.tsv`; `rarity_basis`
   (`measured` | `zero_ci` | `absent`) is the per-variant provenance WITHIN the oracle. A
   per-variant blend was tried and removed: it made two rows in one run comparable on different
   quantities. BOTH ARMS ARE gnomAD v4.1 — the proxy is the cache's own gnomAD AFs — so this is a
   choice of QUANTITY (CI lower bound vs point estimate), never of database. The FAF group set is
   afr/amr/eas/**mid**/nfe/sas, EXCLUDING the bottlenecked ami/asj/fin, so faf95 does not
   reintroduce the MAX_AF trap. **`mid` is the one behavioural difference between the arms**: the
   proxy excludes it, faf95 includes it, so a mid-enriched allele can be gated by one arm and not
   the other; `faf95_group` makes those rows identifiable. **Supplying the slim RETAINS MORE**:
   faf95 ≤ the point estimate, so the same cutoffs stop discarding low-count alleles the interval
   never justified discarding — a SMALLER list after enabling it means a broken join, not a
   better filter.
4. **VEP-centric contract.** The annotation surface is a VEP 115 GRCh38 cache + its score
   PLUGINS — CADD (warns when absent), **SpliceAI** (required by default,
   `resources.vep.spliceai_required`) and **REVEL + AlphaMissense** (required by default,
   `resources.vep.missense_predictors_required`; dedicated files, never dbNSFP — 32 GB for 5
   columns and a dead URL), the stock **NMD** plugin (no data file; behind `nmd_status` and V5)
   and **SpliceVault** (optional, `resources.vep.splicevault`, warns when absent — the 300K-RNA
   empirical mis-splicing outcome beside the SpliceAI event; review evidence that never gates or
   promotes) — plus EXACTLY TWO `bcftools annotate` transfers in Step 2, both
   because the cache cannot supply the field at any price: the **ClinVar sites VCF**
   (`resources.clinvar.vcf` -> `clinvar_CLNREVSTAT` -> `clinvar_stars`; the cache carries
   `CLIN_SIG` but no `CLNREVSTAT`; stars RANK in Step 9 and never gate the screen; absent ->
   WARN) and the **gnomAD v4.1 joint slim** (`resources.gnomad.sites_slim` -> `gnomad_faf95`,
   `gnomad_faf95_group`, `gnomad_nhomalt`, `gnomad_AF_joint`, `gnomad_AF_grpmax`; required under
   the default oracle, rule 3). Each transfer is prefixed with its own namespace (`clinvar_`,
   `gnomad_`) rather than `vep_` precisely so the different oracle is visible. Nothing else is
   transferred: no dbNSFP, no LOFTEE file. Adding an annotation means either a VEP plugin or a
   `bcftools annotate` transfer in Step 2 **and** its INFO field in `annotations.F` — never a
   lookup that reaches around that contract. Remaining accepted losses (see
   [docs/README.md](docs/README.md#canonical-defaults)): no LOFTEE (pLoF confidence) and no
   exome/genome discordance flag. Below MODERATE impact there are TWO keep-paths: **SpliceAI**
   (max delta ≥ `spliceai_ds_min`, default 0.2 — the deep-intronic/synonymous splice signal;
   checked first) and **CADD** (≥ `cadd_phred_supporting`, default 25.3 — everything else
   non-coding). SpliceAI is a VEP plugin over the precomputed RAW scores
   (`resources.vep.spliceai_snv`/`spliceai_indel`); the precomputed set is not exhaustive, so a
   missing score is not "no effect" (keep-only, never drops). Step 2b (live SpliceAI backfill) is
   OFF by default. REVEL/AlphaMissense are inert at the screen by construction and feed only
   Step 9's missense tier.
5. **Gene lists and constraint are priors/tiers, never hard include/exclude** — the "never-drop
   rule" keeps novel-gene discovery alive. Rarity/impact/QC gating happens *before and
   independently of* any list.
6. **The canonical-defaults table is the single source of truth.** Every threshold in code and
   docs must match [docs/README.md#canonical-defaults](docs/README.md#canonical-defaults); all
   thresholds are config defaults, never hardcoded law. A gene-specific ClinGen VCEP value
   overrides a generic cutoff.
7. **Engineering ethos:** fail loudly, verify before claiming, be idempotent (`.done` files +
   integrity checks). Match the existing style.

## Architecture

- **One image** (`Dockerfile`): the group's validated `ensemblorg/ensembl-vep:release_115.0`
  base + a pinned micromamba env (`env/environment.yml`) for bcftools/bedtools/slivar/somalier/
  whatshap/python. The conda env must not shadow VEP's Perl — bioconda drags `perl` in, so the
  Dockerfile removes the conda Perl after solve (VEP's `env perl` must resolve to the base's Perl).
- **Config → env → scripts.** `config/config.example.yaml` is the contract. `src/hprv/config.py`
  loads it, expands `${ENV}`, and emits shell exports (`python -m hprv.config sh`, which
  normalises YAML booleans to `true`/`false` — Python's `True` used to reach the shell verbatim and
  Step 2 read it as "sharding off"). The bash steps take explicit args; `run_pipeline.sh` maps
  config → args.
- **Execution model:** the pipeline is meant to run **inside the container** (tools native on
  PATH). `pipeline/lib/common.sh`'s `hprv_run` auto-detects the runtime: inside an Apptainer
  container it resolves to `native` (direct calls); from a host it wraps each tool in
  apptainer/docker. So the same scripts work in both modes. `HPRV_BIND` carries the dirs that
  must be visible to wrapped calls.
- **Shared python** in `src/hprv/`: `config` (YAML+env), `annotations` (the INFO-field contract
  from Step 2), `genotype` (refined-GQ QC), `ped` (trio parsing). Steps import these so selection
  / inheritance / burden read annotations identically.

## Data contract between steps

- **User input** (git-ignored): a `trios_file` (TSV, header names kid/dad/mom in any order —
  matched by NAME, an `_id` suffix tolerated, and STRICTLY: a header naming some roles but not
  all is an error, never a positional guess, because a transposed mother/father inverts every
  parent-of-origin call and the MIE gate is blind to it; IDs match VCF samples) + a `vcf_dir`/`vcf_list`. `pipeline/resolve_trios.py` maps each trio to
  the VCF containing all three members (exact match; picks the fewest-sample VCF on a tie; extras
  OK), generates PEDs, and writes the **internal manifest** `trios.resolved.tsv`
  (`trio_id  vcf  ped  samples`) that Steps 0/1/4 consume. Unresolved/ambiguous trios are reported
  in `trio_resolution.tsv`, never guessed. Steps 1 and 4 subset each VCF to its 3 members
  (`bcftools view -s`), so extra members and inconsistent sample order don't matter.
- **PED sex**: the generated PED leaves kid sex unknown (`0`); Step 5 reads Step 0's inferred sex
  (`qc_report.tsv`) so X-linked/hemizygous logic fires correctly. Step 0 also infers BOTH
  PARENTS' sex from their own chrX (`parent_sex_flag` — the only reachable sex check, and the
  only direct detector of transposed parents) and each member's no-call rate (`nocall_flag`,
  `qc.max_nocall_rate` 0.10 — a merge-shaped trio VCF reads `./.` for every non-carrier parent
  and Step 5 would silently lose every de novo).
- **Auditing**: every step calls `audit`/`hprv.audit.record` → `audit/counts.tsv`
  (timestamp, step, scope, metric, value; scope = `global` or trio_id). `python -m hprv.audit` assembles
  `audit/summary.md`. Step 3 tags kept variants with `hprv_keep_reason`.
- **VEP runs ONCE** on the cohort union (Step 2). Step 4 transfers annotations with
  `bcftools annotate` — it never re-runs VEP. Keep it that way. Already have a VEP VCF? Set
  `resources.vep.annotated_vcf` and Step 2 ingests it instead (verifying build + frequency
  presence); the rest of Step 2 is unchanged.
- **Step 2 INFO fields** (the contract `src/hprv/annotations.py` owns): all but two are CSQ
  fields lifted by `bcftools +split-vep` with a `vep_` prefix. The exceptions are
  `clinvar_*` (`clinvar_CLNREVSTAT`, `clinvar_CLNSIG`, from the ClinVar sites VCF) and `gnomad_*`
  (`gnomad_faf95`, `gnomad_faf95_group`, `gnomad_nhomalt`, `gnomad_AF_joint` — the
  `zero_ci`/`absent` witness — and `gnomad_AF_grpmax`, from the gnomAD v4.1 joint slim) —
  `bcftools annotate` transfers under their own namespaces on purpose, so which oracle a field
  came from is readable at a glance. Those are the ONLY two transfers. `vep_Consequence`, `vep_IMPACT`, `vep_SYMBOL`, `vep_Gene`, `vep_Feature`,
  `vep_BIOTYPE`, `vep_EXON`, `vep_INTRON`, `vep_HGVSc`, `vep_HGVSp`, `vep_cDNA_position`,
  `vep_CDS_position`, `vep_Protein_position` ("pos/len" under `--total_length`), `vep_MANE_SELECT`,
  `vep_NMD` (the stock NMD plugin: `NMD_escaping_variant` or absent), `vep_CADD_PHRED`, `vep_CLIN_SIG`,
  `vep_SpliceAI_pred_DS_{AG,AL,DG,DL}` (+ `DP_*`, `SYMBOL`; `annotations.spliceai_ds()` = the max,
  the splice keep-path — present only when the SpliceAI plugin is configured), `vep_STRAND` +
  `vep_SpliceVault_{top_events,out_of_frame_events,site_pos,site_type,site_sample_count,site_max_depth,SpliceAI_delta}`
  (the SpliceVault plugin: the ranked natural mis-splicing events at the site SpliceAI predicts
  lost — present only when `resources.vep.splicevault` is set; review columns, read by no gate),
  `vep_REVEL` + `vep_am_pathogenicity`/`vep_am_class` (the calibrated missense predictors —
  **inert at the screen**, consumed only by Step 9's missense tier),
  `vep_gnomAD{e,g}_{AFR,AMR,EAS,NFE,SAS}_AF` (the `grpmax_proxy` arm — read by a rarity gate ONLY
  under `resources.gnomad.oracle: grpmax_proxy`; under the default `faf95` they are reporting
  columns), plus `vep_MAX_AF` / `vep_MAX_AF_POPS` / `vep_gnomAD{e,g}_AF` for REPORTING ONLY —
  never as filter fields under EITHER arm (see rule 2). Add a new annotation by wiring it through
  Step 2's split-vep `want` list AND `annotations.F`.
- **`--pick` vs `--flag_pick`:** Step 2 runs VEP with `--flag_pick`, which keeps EVERY consequence
  block and marks the chosen one `PICK=1`, so split-vep's `-s` selector decides — and an
  externally-produced `--flag_pick` VCF takes the identical path. The selector auto-resolves to
  `pick` when the CSQ has a PICK field, else `worst`; override with `resources.vep.csq_select`.
  Watch out: a `--pick_order` starting with `rank` picks the WORST-consequence transcript, so
  `SYMBOL` can name a non-MANE/readthrough gene that Step 6 then aggregates carriers under.
- **`bcftools +split-vep` does not accept `--threads`** (it is a plugin; passing it aborts the
  step). Thread the VEP call with `--fork` instead.
- **Step 2 idempotency keys and guards.** Per-contig VEP shard `.done` files are keyed on the VEP
  INPUTS (union cksum + cache dir/version + plugin file size/mtime + the `VEP_RECIPE` flag list —
  a data-less plugin such as NMD, or `--total_length`, has no file to key on and used to be served
  from stale shards), so a changed union, a new plugin file or a changed flag re-runs the shard;
  the annotated-union `.done` additionally keys `csq_select`, the oracle and the split-vep `want`
  lift list (a widened lift re-runs split-vep on the cached shards, not VEP), and a union `.done`
  mismatch invalidates the union. SpliceAI and CADD have VALUE-level 0-lift guards
  (`die` when configured but no variant received a score). Under `oracle: faf95` a missing or
  FAILED gnomAD slim transfer is a hard stop, including on the `resources.vep.annotated_vcf`
  ingest path (Step 3 also asserts `gnomad_AF_joint` is declared in the header). Ingest mode dies
  on any multi-ALT record and on union sites absent from the external VCF (Step 4 would otherwise
  drop them silently).
- **Step 4 output**: per-trio `*.candidates.annotated.vcf.gz` + `trios.candidates.tsv`
  (`trio_id  candidates_vcf  ped`). Per-trio VCFs are the authoritative unit — no cohort genotype
  matrix is ever built. The Step 1 and Step 4 per-trio caches are keyed on the source VCF
  (size+mtime), the FILTER/contig settings and the sample list; Step 0's `qc_report.tsv.done` is
  keyed on the manifest. Step 1 keeps FILTER `PASS` **and** `.` (`--filter 'PASS,.'`), and Step 5's
  `filters.genotype_qc.require_pass` treats both as pass — it is not "PASS only".
- **Step 5 output**: `candidates.calls.tsv` (one row per candidate; `mode` ∈ `dominant`
  (inherited het; `flags=origin=mat|pat|both`), `hom_recessive`, `compound_het` (pairs share a
  `pair_id`; a pair whose second hit is a de novo is unphaseable from trio genotypes and carries
  `flags=unphased_denovo_partner`, and does NOT suppress the dominant call; a pair whose
  non-transmitting parent was never affirmatively observed hom-ref carries `origin_unverified`;
  a leg whose transmitting parent failed its own GQ/DP/AB QC carries
  `transmitting_parent_qc_fail` — emitted, never deleted — and none of those three pairs
  consumes its legs, only a confirmed trans pair does),
  `x_linked_recessive`, `denovo`/`denovo_x_hemi` (secondary)). A `1/1` parent transmits obligately,
  so parent-of-origin there is deterministic (`both` is reserved for HET×HET). chrY is routed away
  from the mother-keyed hemizygous models (`male_x_chrx`) and yields no inherited call. Step 5 opens
  VCFs with `strict_gt=True` — cyvcf2's default reports a half-called `0/.` as hom-ref, which would
  defeat every "parent is a confident no-call" test. Modes are configured in
  `inheritance.emit_dominant` / `inheritance.emit_denovo`. Every row carries `rarity_af` (the value
  the run's oracle produced and every gate applied) + `rarity_oracle` + `rarity_basis`, resolved
  ONCE here and consumed downstream, and `child_gt`/`mother_gt`/`father_gt` are cyvcf2 `gt_bases`
  — allele STRINGS such as `A/T` or `T/T`, never `0/1`. `hgvsc`/`hgvsp` (VEP HGVS on the split-vep-selected transcript),
  `exon`/`intron`/`cds_position`/`mane_select`, `nmd_status` (`escaping`/`triggering`/`not_assessed`,
  resolved HERE with the VCF header: a blank NMD-plugin value is `triggering` only for the four
  consequences the plugin grades AND only when `vep_NMD` is declared — a TSV cannot tell that apart
  later) and the SpliceAI event decomposition (`spliceai_event`, `_event_pos`, `_event2`,
  `_shift_nt`, `_shift_frame`, `_effect`, `spliceai_symbol_mismatch`; `src/hprv/splice.py`, floor =
  `spliceai_ds_min`) and the SpliceVault columns beside it (`strand`, `splicevault_top_events`
  verbatim, `_out_of_frame`, `_site_type`, `_site_samples`, `_top1_event` such as `ES` / `CD+12` /
  `CA-31`, `_top1_frame`, and `splicevault_agreement`: `cryptic_confirmed` / `cryptic_unseen` /
  `loss_outcome_supplied` / `site_type_mismatch` / `not_applicable`, blank without the plugin or
  without a SpliceAI event above the floor — SpliceAI's shift is GENOMIC, SpliceVault's offsets are
  transcript-oriented, so the comparison flips sign on the minus strand via `vep_STRAND`) are
  curated columns, and AFTER the curated columns comes a DROPLESS `info_<ID>` block: every INFO field the
  per-trio candidate VCF header declares, verbatim — including the raw multi-transcript `CSQ` as
  `info_CSQ`, i.e. every transcript's HGVS, not only the picked one (the union over trios, in
  header order; a VCF
  `.` reads blank, a Flag reads `1`; `05_inheritance_screen.info_values` parses the VCF LINE so a
  Float keeps the VCF's digits). The prefix is load-bearing: the raw `hiConfDeNovo` (a comma list
  of children) would otherwise collide with the curated flag. This projection was the ONE place an
  annotation could vanish — HGVSc/HGVSp were lifted in Step 2 and never reached a TSV — so keep
  both blocks intact. **Step 6 output**: `genes.ranked.tsv` —
  distinct-individual carrier counts per gene per model (`n_dominant`/`n_biallelic`/`n_xlinked`/
  `n_denovo`), `recurrent` flag (≥ `burden.min_carriers`), the case-only recurrence null
  (`p_recurrence`/`q_recurrence`/`*_exome_wide_sig` — a RANK, never a calibrated test: 2 carriers
  of private variants at N=200 give p≈3e-7 and 3 carriers at N=1000 give 4e-8, i.e. essentially
  "≥3 carriers of private hets"; `p_recurrence_xlinked` uses `--n-male-trios`, which
  `run_pipeline.sh` counts from Step 0's `inferred_sex == 1` among resolved trios, else N_trios;
  `recurrence_kind` for biallelic carriers compares the per-trio SET of variant keys, so two trios
  sharing one comp-het pair read `same_variant`), constraint columns, and — when `--mutrate`
  carries `mu_mis`/`mu_syn`/`mu_lof` — `mu_tot`, `exp_carriers_mu` (= C·μ_g, with C = Σ n_carriers
  / Σ μ over the FULL mutational-target universe, zero-count genes included),
  `carrier_excess_ratio` (n_carriers / exp_carriers_mu) and `p_carrier_excess` (Poisson upper
  tail). `burden.rank_by_mutational_target: true` orders recurrent genes by `p_carrier_excess`
  instead of `best_p`, so long genes no longer lead by size. De novo enrichment: a missing
  `mu_mis` -> no test; a missing `mu_lof` -> imputed as (mu_mis+mu_syn) ×
  `prioritization.excess.offset.mu_lof_impute_factor`, recorded in `dn_mu_src`
  (`gnomad`|`imputed`|`none`).
- **Step 7 output**: `hprv_summary.xlsx` (openpyxl; `src/hprv/report.py`) — documented workbook:
  About/legend + Gene consolidation + Candidate calls + Trio resolution + QC + Audit counts.
- **Step 8 output**: `igv/` for the jlanej/igv.js variant-review server (`src/hprv/igv.py` +
  `08_igv_export.sh`): `variants.tsv` — its headline `frequency` column IS `rarity_af`, the value
  the run's oracle produced and every gate applied; it must never be a different quantity from the
  one that did the filtering, with `grpmax_af`/`faf95` alongside as the raw inputs. `hgvsc`/`hgvsp` sit beside `impact`, and every calls column Step 8
  does not already represent is appended verbatim after the track columns
  (`igv.passthrough_columns`: `flags`, `hiConfDeNovo`, `review_prior_crosscheck`, the Ensembl
  `gene` as `gene_id`, the whole `info_*` block) — the review table is dropless with respect to
  `candidates.calls.tsv`, and it used to be a second silent projection. (Only
  `chrom/pos/ref/alt` required; extra columns are
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
  (`members: carriers|child_only|all`, `confidence`, `min_reads`, `flag_fraction` — default 0.5,
  the ONE threshold read by both Step 8's `nhf_flag` and Step 9's `nhf_status`; the old
  `prioritization.composite.nhf.flag_fraction` key is gone — `memory_mapping`). Gated: activates
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
  Gene-layer counts: `n_observed` = DISTINCT (trio_id, chrom, pos, ref, alt) observations per gene
  and `n_rows` = raw candidate rows (a comp-het leg appears once per pair in
  `candidates.calls.tsv`, and one variant can appear under up to three modes);
  `max_downweight_fraction` is measured over distinct observations. `genotype_qc` derives zygosity
  from the base-form `child_gt` (`T/T`) or the mode, so a hom-alt call is judged on the hom-alt AB
  band, not the het band. `rarity_driven_by_single_group` compares `grpmax_af` (the point
  estimate) with `max_af`, never the oracle value. `gene_is_constrained` reads pLI / LOEUF / s_het
  / `phaplo` (≥ `filters.constraint_weighting.phaplo_min`), matching Step 6. MOI coherence: the
  hemizygous modes `x_linked_recessive` and `denovo_x_hemi` are `coherent` with any XL/XLR/XLD
  curation and `unknown` otherwise — never charged the dominant/recessive discordance penalty.
  `sig_caf_low` flags only genes gnomAD looked at and reported no pLoF CAF, never genes ABSENT
  from the mutational-target table.
  **NEVER-DROP IS AN ASSERTED INVARIANT HERE**: `09_prioritize.py` checks row-count conservation
  before writing, and the down-weight sets a tier + a separately-reported additive penalty + a
  human-readable `downweight_reason` (plus `review_flag`), never a filter. The maximum penalty (−3.0) is sized so it
  cannot alone demote a variant with strong molecular evidence in a constrained gene.
  Two rankings always ship — `rank_agnostic` (no gene-list prior of any kind), `rank_prior`, and
  `rank_delta`. `prioritization.composite.gene_list_prior.enabled` defaults **false** and a config
  overlay path is inert while it is false, so hprv stays phenotype-agnostic by default.

## The annotation surface (authoritative — check here before adding or reading one)

Every annotation the pipeline reads, its provenance, and what happens when it is missing. If you
are adding one, it must appear here, in `annotations.F`, and in a Step-2 producer (a VEP plugin or
one of the two `bcftools annotate` transfers) — nothing may reach around that.

| Annotation | Source | INFO field | Getter | Consumer | Missing ⇒ |
|---|---|---|---|---|---|
| consequence / IMPACT / SYMBOL / Gene / … | VEP cache CSQ | `vep_*` | `consequence()` etc. | Step 3 keep-ladder rung 1 | Step 2 **dies** (core field) |
| CADD | CADD plugin | `vep_CADD_PHRED` | `cadd()` | Step 3 rung 3; Step 9 off-label missense | WARN, screen goes impact+splice-only |
| SpliceAI | SpliceAI plugin | `vep_SpliceAI_pred_DS_{AG,AL,DG,DL}` | `spliceai_ds()` (max of 4) | Step 3 rung 2; Step 9 V4/V3/V0 | **HALT** (`spliceai_required: true`) |
| REVEL | REVEL plugin | `vep_REVEL` | `revel()` | Step 9 missense tier **only** | **HALT** (`missense_predictors_required: true`) |
| AlphaMissense | AlphaMissense plugin | `vep_am_pathogenicity` / `vep_am_class` | `alphamissense()` | Step 9 missense tier **only** | **HALT** (same knob) |
| EXON / INTRON / cDNA·CDS·Protein position | VEP cache CSQ (`--numbers`, `--total_length`) | `vep_EXON`, `vep_INTRON`, `vep_CDS_position`, … | `exon()` / `intron()` / `cds_position()` | Step 5 columns (`exon`, `intron`, `cds_position`); reviewer geometry | blank |
| NMD escape | **NMD plugin** (stock Ensembl, no data file) | `vep_NMD` | `nmd_status(v, plugin_present)` — resolved in Step 5 WITH the header | Step 9's V5 rung (`nmd_status=triggering` only) | `not_assessed` — never a promotion; Step 2 WARNs |
| SpliceAI event | the four DS/DP components (SpliceAI plugin) | `vep_SpliceAI_pred_D{S,P}_*` | `spliceai_event()` → `splice.decompose()` | Step 5 `spliceai_event`/`_effect`/… columns; Step 9 reason strings | no event named (max still reported) |
| SpliceVault | **SpliceVault plugin** (300K-RNA table; optional) | `vep_SpliceVault_*` (+ `vep_STRAND`) | `splicevault_events()` → `splice.parse_splicevault_events()`; `splice.splicevault_agreement()` | Step 5 `splicevault_*` columns (rank-1 event/frame, agreement with the SpliceAI event); Step 9 reason suffix | WARN, columns blank — never a gate, never a promotion |
| ClinVar significance | VEP cache CSQ | `vep_CLIN_SIG` | `clnsig()` / `clnsig_is_plp()` | Step 3 P/LP override; Step 9 clinical | cache always has it |
| ClinVar review status | **transfer** (ClinVar VCF) | `clinvar_CLNREVSTAT` | `clinvar_stars()` | Step 9 clinical damp | WARN, stars `UNAVAILABLE` = FULL weight |
| faf95 | **transfer** (gnomAD joint slim) | `gnomad_faf95` (+ `_group`) | `faf95()` | **every rarity gate** | **HALT** when `oracle: faf95` (the default) |
| nhomalt | **transfer** (same slim) | `gnomad_nhomalt` | `nhomalt()` | Step 9 recessive flag (`nhomalt_recessive_conflict`, 0 points by default) | flag never fires (absent ≠ 0) |
| joint AF / grpmax AF | **transfer** (same slim) | `gnomad_AF_joint`, `gnomad_AF_grpmax` | `gnomad_observed()` | the `zero_ci`/`absent` witness under `oracle: faf95`; otherwise REPORTING ONLY | absent record ⇒ `rarity_basis=absent` (None = rarest) |
| grpmax proxy | VEP cache CSQ | `vep_gnomAD{e,g}_{POP}_AF` | `grpmax_af()` | rarity gates **only** under `oracle: grpmax_proxy` | — |
| MAX_AF / global AF | VEP cache CSQ | `vep_MAX_AF`, `vep_gnomAD{e,g}_AF` | — | **REPORTING ONLY** | never a filter field (rule 2) |

**Exactly TWO `bcftools annotate` transfers exist** — ClinVar and the gnomAD joint slim — both
because the cache cannot supply the field at any price, both under their own INFO namespace, both
guarded by a 0-match `die`. Everything else is a CSQ field lifted by `+split-vep`.

**Three resources HALT at preflight when required-and-missing** (`spliceai_required`,
`missense_predictors_required`, and `oracle: faf95` needing the slim). The pattern is deliberate:
each one silently absent would change a *reported quantity* — the splice keep-path, the
calibration of the missense tier, or which frequency quantity every gate used — without changing
whether the run appears to succeed. Contrast ClinVar, which only damps a score and so warns.

**Provenance columns you must keep populated** when touching any of this — each exists because
two things that look identical in the output are not the same fact:

| Column | Says |
|---|---|
| `rarity_oracle` | run-level: `faf95` or `grpmax_proxy`. ONE per run; the arms never cross |
| `rarity_basis` | per variant: `measured` / `zero_ci` / `absent` **within** that oracle |
| `missense_evidence_source` | `revel` / `alphamissense` / `cadd_offlabel` / `none` — fixed precedence, never a max |
| `E_source` | `gnomad_mu` / `cds_fallback` / `none` (also selects the FDR pool) |
| `mu_lof_src` | `gnomad` / `imputed` / `none` (Step 9 offset) |
| `dn_mu_src` | `gnomad` / `imputed` / `none` (Step 6 de novo expectation) |
| `constraint_source` | which table (`--constraint` first, then `--mutrate`) supplied pLI/LOEUF/oe_syn |
| `moi_coherence` | `coherent` / `discordant` / `unknown` — `unknown` is EXACTLY neutral; hemizygous modes are never `discordant` |
| `spliceai_status` | `scored` / `not_covered` — absence is not "no effect" |
| `nhf_status` | `clean` / `flagged` / `not_screened` — three states, never two |
| `nmd_status` | `escaping` / `triggering` / `not_assessed` — resolved in Step 5 (plugin declared + graded consequence); V5 ONLY via `triggering`, canonical splice never promotes |
| `spliceai_effect` | `site_loss` / `site_gain` / `cryptic_shift_in_frame` / `cryptic_shift_frameshift` / `cryptic_shift` / `paired_gains` / `paired_losses` / `complex` — WHICH SpliceAI event, not just its size |
| `splicevault_agreement` | `cryptic_confirmed` / `cryptic_unseen` / `loss_outcome_supplied` / `site_type_mismatch` / `not_applicable` / blank — how 300K-RNA's observed mis-splicing at the lost site relates to the SpliceAI event; blank = no SpliceVault data or no SpliceAI event, never "disagrees" |
| `clinvar_review_status` | `N_star` or `UNAVAILABLE` (= full weight, NOT 0 stars) |

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
  `mu_g` (our denominator) and the gnomAD rarity oracle (our numerator's filter) unreliable. It is a
  **low-information-locus** flag. Do not "fix" it by inverting the comparison. It fires only for
  genes gnomAD looked at (present in the mutational-target table) and reported no pLoF CAF; a gene
  ABSENT from the table is not flagged — absence from the table is not a measurement.
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
  clears zero — roughly 80% of a chr22 sample. So:
  - **gnomAD HAS the allele** -> faf95 IS 0 -> **rarest, keep**. Do NOT consult the proxy. Of the
    records with no faf95 but a proxy >= 1e-4, **96.5% are AC <= 2** — gnomAD singletons, whose
    point estimate is inflated by a small group AN. Filtering those on the point estimate is
    exactly the error faf95 exists to prevent (Whiffin 2017). The first implementation did this
    and would have dropped ~8% of gnomAD-observed alleles at the dominant gate.
  - **gnomAD has NO record** -> `None` = rarest, `rarity_basis=absent`. The cache proxy is NOT
    consulted — under `oracle: faf95` it is never consulted for any row; it is the rarity field
    only under `oracle: grpmax_proxy`, and then for every row.
  `gnomad_AF_joint` is the WITNESS that separates them — the reason a "reporting only" field is
  actually load-bearing. `rarity_oracle` is the RUN-level constant; `rarity_basis` is the
  per-variant provenance (`measured` / `zero_ci` / `absent`).
  Step 5 resolves this ONCE (`rarity_af` + `rarity_oracle` + `rarity_basis`) and Step 9 consumes
  it rather than re-deriving — a local re-derivation cannot tell `zero_ci` from `absent`, because
  the raw columns are identical. **Step 9 keys on `rarity_oracle`, NOT on `rarity_af`**: the value
  is legitimately empty when the oracle has none for the allele, and keying on it made those rows
  look like a legacy table and fall through to a per-row re-derivation — reintroducing the mixed
  output the design removes.
- **Enabling the gnomAD slim makes the candidate list BIGGER.** faf95 <= the point estimate, so
  the same cutoffs stop discarding low-count alleles whose CI never justified the call. If the
  list got SMALLER after supplying it, the join is broken — check Step 2's
  "gnomAD joint matched N / M sites" line, which is guarded to die at 0 but not at 1.
- **The AD limbs of `sample_qc` fail OPEN while the others fail closed — and that asymmetry
  decided both ways on one missing measurement.** `het`/`hom_alt`/`denovo_child` require
  `ab is not None` and so DROP a carrier when FORMAT/AD is absent; `hom_ref`/`clean_parent`
  returned True when it was absent and so AFFIRMED a non-carrier. A GATK ref-block-derived `0/0`
  parent carries `GT:DP:GQ:MIN_DP:PL` with **no AD** — exactly the shape of a parent at a site
  where the child is het — and `clean_parent` is the ONLY evidence a compound-het pair is in
  TRANS. Hard-failing would drop the ordinary ref-block case wholesale, so the pass stands
  (never-drop) and is now MARKED: `genotype.sample_qc_ad_measured()` is the witness, and Step 5
  emits `parent_ad_unmeasured` / `trans_evidence_unmeasured`. Distinct from `origin_unverified`,
  which means the test FAILED — this one is a vacuous pass.
- **Step 9 joins the gene layer on SYMBOL, not on `gene`.** `candidates.calls.tsv` carries BOTH
  `gene` (an Ensembl ID) and `symbol`, and `_find` returns the FIRST name present. Preferring
  `gene` keyed the whole gene layer on ENSG while every joined resource stayed symbol-keyed
  (`GENE_KEYS_LOWER`): mutational target, constraint, segdup, established genes and the
  gene-prior overlay all matched ZERO rows, `mu_lof` was imputed rather than measured, and the
  run reported success. Step 6 was already symbol-first, so the two now agree, and the resolved
  key is logged + audited. The integration suite CANNOT catch this — its input is Step 8's table,
  which has one `gene` column already holding the symbol — so the guard is
  `tests/test_pure.py:test_prioritize_joins_on_symbol_not_ensembl_id`, and `mock_vep.py` now
  emits a distinct ENSG-style `Gene` (via `zlib.crc32`, NOT `hash()`, which is per-process random
  and would break shard equivalence).
- **A missing component is never zero.** Three separate places substituted 0.0 for an absent
  measurement and kept an unchanged provenance label: a blank NHF read denominator made a 0.9
  non-human fraction read `clean` (0.0 instead of -3.0 points); a missing `mu_mis`/`mu_syn` was
  summed as 0.0 and inflated every `excess_ratio` in that gene by 3.5x while `E_source` still
  said `gnomad`; and `frequency()`'s old proxy fallback. All three are fixed, and the rule is
  general: if you write `_num(x) or 0.0`, you have almost certainly just made absence into
  evidence.
- **Two tables can carry the same constraint column, and the precedence must be ONE direction.**
  `--constraint` and `--mutrate` both carry `pLI`/`oe_lof_upper`/`oe_syn`. The old code resolved
  them with `_find(ccols,…) or _find(mcols,…)`, which gave oe_syn mutrate-first and pLI/LOEUF
  constraint-first — so a gene's rank could be decided by which table it was in. Worse, it
  resolved ONE column name and used it to index the OTHER file: the prepared constraint file is
  lowercase (`pli`, via `join_constraint.py`) while gnomAD's own table is `pLI`, so the fallback
  read BLANK despite the value sitting in `src_mutrate_pLI`. Resolve PER TABLE, constraint-first
  throughout, and record the winner in `constraint_source`.
- **Alias chains consume the wrong table silently.** `_find` returns the first present name, so a
  bare `lof`/`mis` matches denovolyzeR-shaped tables and a gnomAD **v4** `oe_lof_upper` wins the
  LOEUF chain and is then compared against the **v2-calibrated** `loeuf_v2_tier1` (0.35). Step 6
  now logs and audits the resolved column for every field and warns on a v4-named constraint path.
- **A same-named column is not the same measurement.** `segdup98_frac` used to fall back from
  `--segdup` to a same-named column in `--mutrate` — contradicting the documented "absent =>
  signal off (WARN)" contract, leaving the signal silently ON, and crossing coordinate builds
  (the segdup track is hg19; the mutrate table is the gnomAD v2.1.1 projection). Now `--segdup`
  only, with a warn if `--mutrate` carries the column.
- **The MANE-only SpliceAI mirror must never land at the full file's filename.** It used to be
  written to `spliceai_scores.raw.snv.hg38.vcf.gz` — the BaseSpace name — so `ls` could not tell
  the builds apart, every check downstream is existence-only, and `spliceai_required: true` halts
  on ABSENCE, which made the wrong file look exactly like the right one. A MANE-only file blanks
  non-MANE transcripts and `selection.py` reads a blank as "no rescue". It now lands under its own
  name and must be symlinked in deliberately.
- **Step 8b stamps the mask it APPLIED, not the one it requested.** When the `-F` filter fails it
  classifies the unfiltered mini-CRAM (deliberately — a missing NHF row is worse). It used to
  stamp the requested mask in the `.done` key, so a re-run reused the unfiltered row forever, and
  one member's `*_nhf_reads` could count PCR duplicates while its siblings' did not inside one
  `variants.tsv`.
- **`MAX_AF` is a trap, not a shortcut.** It is right there in the CSQ and looks like the rarity
  field. It is not — see golden rule 2. It maxes over founder groups (ami AN≈900) and 1000G
  populations that gnomAD's grpmax excludes on purpose, so a single allele reads as AF≈1e-3 and
  silently drops dominant candidates. Rarity comes from `annotations.frequency()` ONLY.
  `tests/test_pure.py:test_frequency_excludes_bottlenecked_pops` and the GENEFND integration case
  exist to catch a regression here.
- **Absence of a cached AF is weak evidence.** VEP caches frequencies only for alleles
  accessioned into **dbSNP** — an un-accessioned gnomAD variant returns no AF and reads as
  "absent ⇒ rarest". Ensembl itself recommends `--custom` with the gnomAD VCF over `--af_gnomad*`
  for this reason. This is a `grpmax_proxy`-arm property only: the joint slim carries every gnomAD
  allele regardless of dbSNP accession, so under the default `faf95` oracle the witness is the
  slim (`gnomad_AF_joint`), not the cache. On the proxy arm we accept it; it biases toward
  retention.
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
  VEP-centric contract. Kept in the image so re-enabling it is a config change, not a rebuild.
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
- **The screen's knobs are validated at Step 3 and Step 5 start-up, and a boolean is read
  strictly.** `hprv.config.validate_filters` halts on the settings that used to switch a rung
  OFF with exit 0: a scalar `keep_impacts: HIGH` (`set()` of a string is `{'H','I','G'}` — every
  stop_gained then filed under `not_functional`), a lower-cased or unknown IMPACT, an inverted
  rarity ladder (`dominant_max` above `recessive_max` is silently inert), a percent-scale or
  inverted allele-balance band (`het_ab_min: 25` reduced Step 5 to zero calls), and a quoted or
  `${ENV}`-templated `"false"`, which `bool()` read as True. Read booleans through
  `config.get_bool` / `as_bool`, never `bool(get(...))`; read `keep_impacts` through
  `config.keep_impacts`.
- **`child_gt`/`mother_gt`/`father_gt` are cyvcf2 `gt_bases`, never `0/1`.** Step 5 writes allele
  STRINGS (`A/T`, `T/T`, `./.`). Any consumer that tests `gt in ("1/1", "1|1")` never matches, so
  Step 9's `genotype_qc` used to judge every hom-alt call on the het AB band and fail it; zygosity
  is now derived from the base-form genotype (both alleles equal and non-ref) or from the mode.
- **Step 9's `n_observed` counts DISTINCT observations, not rows.** One candidate row is not one
  observation: a comp-het leg appears once per pair in `candidates.calls.tsv`, and a variant can be
  emitted under up to three modes. `n_observed` is the number of distinct (trio_id, chrom, pos,
  ref, alt) per gene; `n_rows` is the raw row count, kept beside it so the inflation is visible.
  `max_downweight_fraction` is measured over distinct observations. Counting rows inflated the
  excess statistic exactly where comp-het pairing is most active (long genes).
- **Step 0's Mendelian-error rate is measured on a CAPPED scan, not genome-wide.** `qc.max_sites`
  (default 200000) caps the MIE/CHARR scan at the first `max_sites` QC-passing autosomal biallelic
  sites, and the chrX sex scan is capped the same way. Quote `mie_rate` as a rate over `n_sites`,
  not as a genome-wide figure.
- **The SLURM graph has a `calls` phase between `gather` and `downstream`.** It runs Steps 3-5 and
  PLANS the Step-5b SpliceAI rescoring array (`--rescore-emit-manifest`), then submits
  `rescore-scatter` (one chunk per task, `--rescore-chunk N`) → `rescore-gather` → `downstream`,
  or `downstream` alone on an empty manifest. `DOWN_FROM` therefore defaults to **6**. A sub-task
  never re-runs Step 5 proper (concurrent tasks must not rewrite `candidates.calls.tsv`); the
  scores are cached PER VARIANT in `spliceai_rescore/scores.tsv` because Step 5 rewrites the
  calls table every run, and a changed `-D` invalidates the cache. Guards:
  `tests/test_slurm_rescore_phases.sh`, `tests/test_spliceai_rescore.sh`.
- **SLURM `DOWN_TO` may be 9** (Step 9 runs on the downstream node); `nhf-plan` now dies when no
  kraken2 DB was supplied instead of scheduling nothing, and the CI smoke test's probes fail the
  job rather than merely printing.
- **`prepare_resources.sh verify` mirrors the DEFAULT config.** `gnomad_sites`, `revel` and
  `alphamissense` are reported as REQUIRED (missing -> non-zero exit) unless the config opts down
  (`oracle: grpmax_proxy`, `missense_predictors_required: false`), and `emit-env` emits
  `GNOMAD_SITES` uncommented when the slim is present. `.gitignore` covers `*.vcf.bgz` /
  `*.bgz.tbi` so a prepared slim or a bgzipped constraint table cannot be committed by accident.
- **Apptainer:** point `APPTAINER_TMPDIR`/`CACHEDIR` at real disk and **do not use
  `--containall`** — a tmpfs `/tmp` OOMs heavy VEP/sort (documented failure in the group's
  original annotate script; `common.sh` already sets a disk-backed workdir).
- **Reproducibility hardening:** `environment.yml` pins versions; for byte-identical rebuilds run
  `conda-lock` and point the Dockerfile at the lockfile. Pull the image by `@sha256:` digest.

## Running

- **Resources (one-time):** on the host — a **VEP 115 GRCh38 cache** (~24 GB; it carries gnomAD
  v4.1 point AFs + ClinVar `CLIN_SIG`), the **CADD** SNV+indel files (warn if absent), the
  **SpliceAI** raw score files (`resources.vep.spliceai_snv`/`spliceai_indel`; a one-time BaseSpace
  login download via `scripts/download_spliceai.sh`), **REVEL + AlphaMissense** (dedicated files,
  `resources.vep.revel`/`alphamissense`), the **gnomAD v4.1 joint slim** (~10 GB,
  `resources.gnomad.sites_slim`, `prepare_resources.sh --only gnomad_sites fetch`) and the
  **ClinVar sites VCF** (~0.18 GB, `resources.clinvar.vcf`, warns if absent). Three of those HALT
  the run at preflight when Step 2 invokes VEP and they are missing: SpliceAI
  (`resources.vep.spliceai_required: true`), REVEL/AlphaMissense
  (`resources.vep.missense_predictors_required: true`) and the gnomAD slim (required by the default
  `resources.gnomad.oracle: faf95`); each has a documented opt-down (`false` / `grpmax_proxy`).
  Not enforced when `resources.vep.annotated_vcf` is set — except the gnomAD slim, which the
  ingest path still needs under `faf95`. No dbNSFP or LOFTEE download.
  `prepare_resources.sh --dir DIR fetch|verify|emit-env` fetches these, runs INSIDE the image
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
  genotype QC incl. the AD fail-open asymmetry and the FORMAT/DP fallback, the selection funnel,
  **Step 5's whole inheritance model** (a stubbed-cyvcf2 `FakeVar` harness: origins, obligate
  transmission, comp-het phasing, the trans-evidence flags, X/Y/PAR and the sex gates), Step-6
  helpers **and its counting/ranking through `main()` with a stubbed scipy**, and the Step-9
  prioritization layer — the NB fit/tail/BH-FDR, the never-drop invariant end-to-end through the
  CLI, the positive-control guard, both tier ceilings, blank-vs-zero NHF, mechanism gating).
  **89 tests, no network and no VCF.**
  **Two documented exceptions to "no heavy deps":** the tests that drive `09_prioritize.py:main()`
  or `06_gene_burden.py:main()` need `yaml` transitively (`load_config` does `import yaml`), and the
  workbook test needs `openpyxl`. They declare it at the `_requires()` chokepoint and **SKIP**
  without it — and `_run_all` then refuses to print "All N passed", instead reporting
  `76 passed, 13 SKIPPED ... NOT full coverage`, because the skipped set holds the never-drop and
  cache-invalidation guards. **CI installs pyyaml + openpyxl AND fails the job on that
  "NOT full coverage" line**, so a skipped test can never read as a green run. Before calling this
  suite green, run it the way CI does — with those two installed, not an env that happens to carry
  the container's packages. (This bit us twice: the Step-9 CLI tests were reported passing from a
  local env with pyyaml while CI died with `ModuleNotFoundError`, and two tests later
  caught `ImportError` and `return`ed, printing PASS while asserting nothing.)
  **A guard test must be able to FAIL.** Two could not, and both were found by mutation: the
  trim-loop HALT test used a uniformly-inflated fixture (C absorbs it, nothing is trimmed), and the
  config-vs-code check compared the config to literals *in the test*. Deleting the HALT, or drifting
  three code defaults, left the whole suite green. The fixtures now trip the guard, and
  `test_prioritize_config_matches_canonical_defaults` probes **both directions** — every code
  default is evaluated with `cfg={}` and with the shipped config on inputs either side of each cut
  point, and the dead keys are asserted ABSENT.
- **End-to-end integration** (`tests/integration/`): `run_integration.sh` generates a tiny
  self-consistent mock genome + trios (`make_mock_data.py`) engineered to exercise every mode
  and filter path, runs resolve + Steps 0,1,3,4,5,6,8,9 with REAL bcftools + the python steps (only
  Step 2's VEP call is mocked via `mock_vep.py`; the gnomAD joint-slim transfer is a REAL
  `bcftools annotate` onto a tiny mock slim), and asserts the resolution, funnel, and
  calls (`assert_integration.py`). The mock deliberately includes a **GATK ref-block parent**
  (`GT:DP:GQ`, no AD) and a **half-called `0/.` parent**, which is how the `parent_ad_unmeasured`
  flag, the AD fail-open asymmetry and the `strict_gt=True` contract are exercised end-to-end —
  and how the `genotype.dp()` FORMAT/DP fallback was found. The Step-9 fixtures are engineered so the artifact locus and the
  established-gene positive control have the **same** extreme excess shape — the only thing
  separating them is the control ceiling, which is exactly what the assertion tests. Two mock-scale
  config deviations are deliberate and commented in `make_mock_data.py`: `min_control_genes` 1000→3
  and `max_downweight_fraction` 0.20→1.0, both because the mock has a handful of genes rather than
  an exome. Runs in CI on host bcftools — no image build needed. To run
  locally you need bcftools/samtools/bgzip/tabix + a python with cyvcf2/pysam/scipy/pyyaml on PATH.
  `assert_shard_equivalence.sh` drives Step 2 DIRECTLY with no gnomAD slim, so it sets
  `HPRV_GNOMAD_ORACLE=grpmax_proxy` explicitly — under the default faf95 oracle Step 2 now halts
  rather than run with no oracle value at all.
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
SOTA review. Done: the case-only recurrence null + FDR (Step 6 — a rank, never a calibrated
test) plus the mutational-target carrier expectation (`p_carrier_excess`), CHARR/freemix
contamination gate (Step 0). Quick wins still open: somalier ancestry/relatedness, PP1/BS4 co-segregation,
UTRannotator, UPD rescue, conda-lock. Big bets: germline CNV (GATK-gCNV), phenotype ranker
(Exomiser), read-backed/population phasing, ROH.

See also **[docs/robustness_audit_2026-09.md](docs/robustness_audit_2026-09.md)** — the consolidated
robustness audit of the plausibility filter (Step 3) and the genotype-to-meaning assignment
(Step 5): a drop ledger showing where the pipeline loses a variant without a trace, 106 verified
findings with the finder's and the verifiers' severities side by side, and a tiered work order.
Its resolution table at the top records what landed: Tiers 1–2 (the input-side counters, the
Step-4/5 transfer witness, the ClinVar counter, strict parent columns + parental sex check, the
transmitting-parent flag and the comp-het veto rule, config validation, Step-9 freshness) are done;
Tier 3 (unassigned-call rows, PGT/PID phasing, all-gene pairing, multiallelic AB, hom-rec with one
non-carrier parent) is open and each item changes which variants are recovered — decide them from
the `no_row.*` counts a real run now produces.

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
