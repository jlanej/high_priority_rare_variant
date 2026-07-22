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
2. **gnomAD v4.1 from the VEP cache is the ONLY population-frequency oracle.** These trios are not
   jointly genotyped, so internal cohort AC/AN is meaningless (absent ≠ hom-ref). Never `bcftools
   merge` the trios into a genotype matrix and never compute population frequency from internal
   counts; internal recurrence is valid only as an artifact/blocklist signal.
   The rarity field is a **grpmax proxy** — max AF over the grpmax-*eligible* groups
   (`annotations.GRPMAX_POPS` = AFR/AMR/EAS/NFE/SAS) — read from the CSQ. Two hard rules:
   - **Never substitute VEP's `MAX_AF`.** It maxes over the bottlenecked founder groups grpmax
     deliberately excludes (ami AN≈900, asj, fin, mid) and over tiny 1000G populations; one allele
     there reads as AF≈1e-3 and silently kills dominant candidates at the 1e-4 gate.
   - **Never substitute the global AF.** It dilutes ancestry-enriched variants and fails the
     opposite way (retaining benign polymorphisms). The two wrong substitutions err in opposite
     directions — there is no single safe fallback.
   This is a **point estimate, not `faf95`**: faf95's CI correction needs AC/AN, which the cache
   does not carry, so it is unrecoverable rather than approximated. It is the deliberate cost of
   the VEP-only contract (rule 6). Everything reads it through `annotations.frequency()` — one
   chokepoint, never a field getter directly.
6. **VEP-only contract.** The annotation surface is a VEP 115 GRCh38 cache + the CADD plugin.
   Nothing else is downloaded or bind-mounted: no gnomAD, ClinVar, dbNSFP, SpliceAI or LOFTEE
   file. Adding an annotation means adding a `bcftools annotate` transfer in Step 2 **and** its
   INFO field in `annotations.F` — never a lookup that reaches around that contract. Known,
   accepted losses (see [docs/README.md](docs/README.md#canonical-defaults)): no faf95 CI, no
   nhomalt, **no SpliceAI ⇒ deep-intronic/synonymous splice variants are invisible**, no LOFTEE,
   no ClinVar stars. CADD is consequently the ONLY functional predictor and the ONLY keep-path
   below MODERATE impact.
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
  (step, scope, metric, value; scope = `global` or trio_id). `python -m hprv.audit` assembles
  `audit/summary.md`. Step 3 tags kept variants with `hprv_keep_reason`.
- **VEP runs ONCE** on the cohort union (Step 2). Step 4 transfers annotations with
  `bcftools annotate` — it never re-runs VEP. Keep it that way. Already have a VEP VCF? Set
  `resources.vep.annotated_vcf` and Step 2 ingests it instead (verifying build + frequency
  presence); the rest of Step 2 is unchanged.
- **Step 2 INFO fields** (the contract `src/hprv/annotations.py` owns): ALL of them are CSQ
  fields lifted by `bcftools +split-vep` with a `vep_` prefix — there are no `hprv_*` transfers
  any more. `vep_Consequence`, `vep_IMPACT`, `vep_SYMBOL`, `vep_Gene`, `vep_Feature`,
  `vep_BIOTYPE`, `vep_HGVSc`, `vep_HGVSp`, `vep_MANE_SELECT`, `vep_CADD_PHRED`, `vep_CLIN_SIG`,
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
  per-trio VCF tracks `vcfs/<trio>.vcf.gz`, `trios.tsv`, `sample_qc.tsv`, empty `curation.json`.

## Gotchas that WILL bite you

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
- **The predictor branches you may be tempted to re-add were dead.** REVEL/AlphaMissense/MPC are
  missense-only; every missense is `IMPACT=MODERATE`; `selection.py` keeps it at an earlier branch
  and returns. Re-adding dbNSFP without also narrowing `keep_impacts` buys exactly nothing at the
  screen (it is only worth reporting/tiering value). CI asserts these reasons never fire.
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
- **Apptainer:** point `APPTAINER_TMPDIR`/`CACHEDIR` at real disk and **do not use
  `--containall`** — a tmpfs `/tmp` OOMs heavy VEP/sort (documented failure in the group's
  original annotate script; `common.sh` already sets a disk-backed workdir).
- **Reproducibility hardening:** `environment.yml` pins versions; for byte-identical rebuilds run
  `conda-lock` and point the Dockerfile at the lockfile. Pull the image by `@sha256:` digest.

## Running

- **Resources (one-time, small):** under the VEP-only contract you need exactly two things on the
  host — a **VEP 115 GRCh38 cache** (~24 GB; it carries gnomAD v4.1 frequencies + ClinVar) and the
  **CADD** SNV+indel files. Nothing else: no gnomAD, ClinVar, dbNSFP, SpliceAI or LOFTEE download.
  `prepare_resources.sh --dir DIR fetch|verify|emit-env` still fetches these, runs INSIDE the image
  (baked in at `/opt/hprv/scripts`, on PATH, with its pinned `/opt/hprv/resources/manifest.env`;
  `HPRV_RESOURCE_MANIFEST` re-pins without a rebuild), and `emit-env` writes the `${ENV}` exports
  the config expects. Never bake resource DATA into the image (`.dockerignore` keeps `resources/*`
  except the manifest out of the build context). See [docs/resources.md](docs/resources.md).
  Already have a VEP VCF? Skip all of it: set `resources.vep.annotated_vcf` and Step 2 ingests it.
- **HPC (primary):** `apptainer exec --cleanenv --bind ... hprv.sif run_pipeline.sh --config
  config/config.yaml` (add `--from N --to M` for a subset).
- **Dev/host:** run individual step scripts; python steps need `PYTHONPATH=src` and the container
  env (cyvcf2/pysam/scipy/pyyaml) — easiest is to exec them inside the image.

## Testing

- **Host, no heavy deps:** `python3 -m py_compile` all scripts; `bash -n` all shell;
  `python3 tests/test_pure.py` (pure-logic: config, ped, trios-file parsing, annotation getters,
  genotype QC, selection funnel, Step-6 helpers).
- **End-to-end integration** (`tests/integration/`): `run_integration.sh` generates a tiny
  self-consistent mock genome + trios (`make_mock_data.py`) engineered to exercise every mode
  and filter path, runs resolve + Steps 0,1,3,4,5,6 with REAL bcftools + the python steps (only
  Step 2's VEP call is mocked via `mock_annotate.py`), and asserts the resolution, funnel, and
  calls (`assert_integration.py`). Runs in CI on host bcftools — no image build needed. To run
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
