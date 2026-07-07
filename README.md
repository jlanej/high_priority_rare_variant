# high_priority_rare_variant

Screen **GMKF Kids First per-trio VCFs** for *high-priority rare variants* — and for
**genes carrying more rare variants than expected** — in rare disease and germline
pediatric cancer.

The inputs are per-trio VCFs (GRCh38) produced by GATK's genotype-refinement workflow,
**not** jointly genotyped across the cohort. Everything runs from **one container** under
Apptainer on HPC, driven by a single config file.

> ⚠️ **Public repository — controlled-access data.** GMKF/Kids First data is dbGaP
> controlled-access. **Never** commit VCF/BAM/CRAM/PED files, real filesystem paths,
> sample/subject identifiers, or any PHI. All inputs, resources, and paths are supplied at
> runtime via `${ENV}` placeholders in the config. The `.gitignore` enforces this — keep it
> that way.

## What it does

A resolve preflight + seven-step flow (see **[docs/pipeline_design.md](docs/pipeline_design.md)**
for the vetted design and the artifact each step produces):

| Step | What | Output |
|------|------|--------|
| resolve | Map each `kid/dad/mom` trio to the VCF containing all three members (exact sample-ID match; extras OK); generate PEDs | `trios.resolved.tsv`, `trio_resolution.tsv`, `peds/` |
| 0 | Per-trio QC gate (Mendelian error + sex) | `qc_report.tsv` |
| 1 | Subset to trio members, normalize, build a **site-only union** of loci (never a genotype merge) | `cohort.sites.vcf.gz` |
| 2 | Annotate the union **once** (VEP + LOFTEE + dbNSFP + SpliceAI + gnomAD faf95 + ClinVar) — **VEP is never run per trio** | `cohort.sites.annotated.vcf.gz` |
| 3 | Select biologically-plausible sites (rarity + function; ClinVar P/LP override); tag each with *why* it was kept | `plausible.sites.vcf.gz` |
| 4 | Recover **real per-trio genotypes** at plausible sites + transfer annotations | per-trio `*.candidates.annotated.vcf.gz` |
| 5 | Pedigree-aware inheritance screen + genotype QC (de novo / recessive / comp-het / X) | `candidates.calls.tsv` |
| 6 | Cross-pedigree gene burden (de novo enrichment + constraint weighting) | `genes.ranked.tsv` |

Every step records input/output counts and funnel tallies to `audit/counts.tsv`, assembled into
`audit/summary.md` — a global + per-trio "what went where and why" (see [Auditing](#auditing)).

The methodology — thresholds, tool choices, and the evidence behind them — is documented and
source-cited in **[docs/](docs/README.md)**. Every default lives in one place:
**[Canonical defaults](docs/README.md#canonical-defaults)**.

## Design principles

- **gnomAD v4.1 `faf95` is the only population-frequency oracle.** Because the trios are not
  jointly genotyped, internal cohort AC/AN is meaningless (absent ≠ hom-ref) and is used only
  as an artifact/blocklist signal.
- **Gene lists and constraint are priors/tiers, never hard filters** ("never-drop rule") — so
  novel-gene discovery survives.
- **One container, config-driven, no hard paths.** The same image runs on a laptop (Docker) and
  on HPC (Apptainer); scripts only ever call tools inside it.
- **Fail loudly, verify before claiming, be idempotent** (`.done` files; integrity checks).

## Quickstart

```bash
# 1. Get the image (built + published to GHCR on every commit). Pull by digest on HPC.
apptainer pull hprv.sif docker://ghcr.io/<owner>/high_priority_rare_variant:latest

# 2. Configure. Copy the example and point the ${ENV} placeholders at YOUR resources.
cp config/config.example.yaml config/config.yaml     # config.yaml is git-ignored
export REF_FASTA=/path/to/GRCh38.fa VEP_CACHE=/path/to/vep GNOMAD_SITES=/path/to/gnomad.vcf.gz
export CLINVAR_VCF=... HPRV_WORK=/path/to/work ...
#    (see config.example.yaml for the full list; verify gnomAD INFO tag names with
#     `bcftools view -h $GNOMAD_SITES | grep INFO`)

# 3. Provide inputs (git-ignored):
#    - a trios file: TSV with a header naming kid/dad/mom (any order); IDs match the VCFs:
#          #kid   dad    mom
#          CH1    FA1    MO1
#      export TRIOS_FILE=/path/to/trios.tsv
#    - the VCF source (a directory and/or a list file):
#          export VCF_DIR=/path/to/vcfs        # globbed for *.vcf.gz/*.vcf/*.bcf
#    The pipeline finds the VCF containing all three members for each trio and generates
#    the internal manifest + PEDs automatically — sample order within a VCF does not matter,
#    and a VCF may contain additional members.

# 4. Run the whole pipeline inside the container.
apptainer exec --cleanenv \
    --bind "$(dirname "$REF_FASTA")" --bind "$VEP_CACHE" --bind "$HPRV_WORK" --bind "$VCF_DIR" \
    hprv.sif run_pipeline.sh --config config/config.yaml
```

## Auditing

The run is fully answerable — "what went where and why":
- **`trio_resolution.tsv`** — for every kid: resolved / unresolved (and which member was
  missing), the chosen VCF, and whether multiple VCFs matched.
- **`audit/counts.tsv`** — every step's input/output counts and funnel tallies (global and
  per-trio), including Step 3's keep/drop reasons.
- **`audit/summary.md`** — assembled global variant funnel (union → annotated → plausible) and a
  per-trio table (candidate genotypes → candidate calls by inheritance mode).
- Each retained variant carries an **`hprv_keep_reason`** INFO tag; each Step-5 call row carries
  its inheritance mode, the evidence annotations, and (for de novo) a `review_prior_crosscheck` flag.

Run a subset with `--from N --to M`. Every step is idempotent, so re-running resumes where it
stopped. On the laptop/dev path the same step scripts work through Docker automatically (the
container-exec layer auto-detects the runtime).

## Repository layout

```
docs/            source-cited methods reference + vetted pipeline design
config/          config.example.yaml (the contract; every tunable, no real paths)
env/             environment.yml (pinned conda toolchain layered onto the VEP image)
Dockerfile       one image: Ensembl VEP 115 base + bcftools/slivar/somalier/python...
pipeline/        resolve_trios.py + step scripts (00..06) + run_pipeline.sh + lib/common.sh
src/hprv/        shared python: config, annotations, genotype QC, ped, selection, audit
.github/workflows build + publish to GHCR on every commit (provenance + SBOM)
```

## Known scope limitations

SNV/indel only initially — **CNV/SV are a real blind spot** (10–15% of pediatric-cancer/rare-
disease diagnoses); pseudogene/seg-dup regions (*PMS2*, *CYP21A2*, *SMN1*) are low-confidence
from short reads; proband mosaicism needs a dedicated tier; the phenotype (Exomiser/HPO) layer
and de novo enrichment **calibration** (synonymous λ ≈ 1) are planned. The scripts are
implemented and syntax-checked but **not yet validated on real data** — GIAB/CMRG truth sets and
a positive-control panel are the next step. See
[pipeline_design.md](docs/pipeline_design.md#known-scope-limitations-stated-honestly-not-hidden).

## License

MIT — see [LICENSE](LICENSE).
