# CLAUDE.md — implementer guide

Orientation for anyone (human or Claude) extending this repository. Read this first, then
[docs/pipeline_design.md](docs/pipeline_design.md) and
[docs/README.md#canonical-defaults](docs/README.md#canonical-defaults).

## What this is

A container-based pipeline that screens **GMKF Kids First per-trio VCFs** (GRCh38, GATK
genotype-refinement output, **not** jointly genotyped across the cohort) for high-priority rare
variants and for genes carrying more rare variants than expected, targeting **rare disease and
germline pediatric cancer**. Runs under Apptainer on HPC.

## Golden rules (do not violate)

1. **Public repo, controlled-access data.** Never commit VCF/BAM/CRAM/PED, real filesystem
   paths, sample/subject IDs, or PHI. All paths are `${ENV}` placeholders resolved at runtime.
   `.gitignore` enforces this; before every commit, sanity-check with
   `grep -rE '/Users/|/scratch|/home/[a-z]|BS_[A-Z0-9]{8}'`.
2. **gnomAD v4.1 `faf95` is the ONLY population-frequency oracle.** These trios are not jointly
   genotyped, so internal cohort AC/AN is meaningless (absent ≠ hom-ref). Never `bcftools merge`
   the trios into a genotype matrix and never compute population frequency from internal counts;
   internal recurrence is valid only as an artifact/blocklist signal.
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
  whatshap/python. The conda env is deliberately perl-free so it never shadows VEP's Perl.
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
  `bcftools annotate` — it never re-runs VEP. Keep it that way.
- **Step 2 INFO fields** (the contract `src/hprv/annotations.py` owns): VEP/plugin fields are
  lifted by `bcftools +split-vep` with a `vep_` prefix (`vep_Consequence`, `vep_IMPACT`,
  `vep_SYMBOL`, `vep_REVEL_score`, `vep_AlphaMissense_score`, `vep_MPC_score`, `vep_CADD_PHRED`,
  `vep_SpliceAI_pred_DS_*`, `vep_LoF`, ...). External transfers are `hprv_gnomad_af`,
  `hprv_gnomad_grpmax_af`, `hprv_gnomad_faf95`, `hprv_gnomad_nhomalt`, `hprv_clnsig`,
  `hprv_clnrevstat`, `hprv_clnsigconf`. Add a new annotation by wiring it through Step 2's
  split-vep/annotate list AND `annotations.F`.
- **Step 4 output**: per-trio `*.candidates.annotated.vcf.gz` + `trios.candidates.tsv`
  (`trio_id  candidates_vcf  ped`). Per-trio VCFs are the authoritative unit — no cohort genotype
  matrix is ever built.
- **Step 5 output**: `candidates.calls.tsv` (one row per candidate; comp-het pairs share a
  `pair_id`). **Step 6 output**: `genes.ranked.tsv`.

## Gotchas that WILL bite you

- **gnomAD INFO tag names.** Defaults (`AF_joint`, `AF_grpmax_joint`, `fafmax_faf95_max_joint`,
  `nhomalt_joint`) match gnomAD v4.1 *joint* conventions but **vary by how the file was
  subset/downloaded**. Verify with `bcftools view -h $GNOMAD_SITES | grep INFO` and set the
  `resources.gnomad.*_tag` config keys. Step 2 fails loudly if a tag is absent.
- **LOFTEE data filenames vary** across builds/forks. If the default LoF plugin string is wrong,
  set `HPRV_LOF_PLUGIN` to the full plugin argument. LOFTEE must be the **GRCh38** fork.
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

- **HPC (primary):** `apptainer exec --cleanenv --bind ... hprv.sif run_pipeline.sh --config
  config/config.yaml` (add `--from N --to M` for a subset).
- **Dev/host:** run individual step scripts; python steps need `PYTHONPATH=src` and the container
  env (cyvcf2/pysam/scipy/pyyaml) — easiest is to exec them inside the image.

## Testing

- **Host, no heavy deps:** `python3 -m py_compile` all scripts; `bash -n` all shell; pure-logic
  checks for `config`, `ped`, `annotations._STAR`, and the Step-6 helpers (`bh_fdr`, `classify`).
- **In-container integration** (TODO): add a tiny synthetic trio fixture under `tests/data/` and
  run the flow end-to-end inside the built image in CI.
- **Validation (TODO):** GIAB/CMRG truth sets + a positive-control variant panel to measure
  sensitivity/precision of de novo and prioritization logic; extend the synonymous-λ≈1
  calibration pipeline-wide.

## Open TODOs

- **CNV/SV module** (GATK-gCNV / Manta / ExomeDepth) — the biggest coverage gap.
- **Pseudogene/seg-dup handling** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) — flag/annotate.
- **Proband mosaicism tier** (low-VAF calls outside the het AB band).
- **Phenotype layer** — Exomiser/LIRICAL HPO ranking as a *prior* (never-drop), plus HPO ingestion
  that degrades gracefully when phenotype is sparse.
- **Burden calibration** — supply an unfiltered synonymous de novo count so Step 6 can verify
  λ ≈ 1; add TRAPD as the corroborative case-vs-gnomAD signal.
- **Pediatric-cancer overlay** — implement the ACMG SF v3.3 / PanelApp-green tiering and
  second-hit boost as a reporting overlay.
- **conda-lock + digest-pinned base image**; **real-data validation**.

## Commit conventions

Small, logical commits (`chore:` / `docs:` / `build:` / `feat:`). End messages with
`Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`. Work on a branch; open a PR into
`main`.
