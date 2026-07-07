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

A seven-step flow (see **[docs/pipeline_design.md](docs/pipeline_design.md)** for the vetted
design and the artifact each step produces):

| Step | What | Output |
|------|------|--------|
| 0 | Per-trio QC gate (Mendelian error + sex) | `qc_report.tsv` |
| 1 | Normalize + build a **site-only union** of loci (never a genotype merge) | `cohort.sites.vcf.gz` |
| 2 | Annotate the union **once** (VEP + LOFTEE + dbNSFP + SpliceAI + gnomAD faf95 + ClinVar) | `cohort.sites.annotated.vcf.gz` |
| 3 | Select biologically-plausible sites (rarity + function; ClinVar P/LP override) | `plausible.sites.vcf.gz` |
| 4 | Recover **real per-trio genotypes** at plausible sites + transfer annotations | per-trio `*.candidates.annotated.vcf.gz` |
| 5 | Pedigree-aware inheritance screen + genotype QC (de novo / recessive / comp-het / X) | `candidates.calls.tsv` |
| 6 | Cross-pedigree gene burden (de novo enrichment + constraint weighting) | `genes.ranked.tsv` |

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
export CLINVAR_VCF=... TRIO_MANIFEST=/path/to/trios.tsv HPRV_WORK=/scratch/run1 ...
#    (see config.example.yaml for the full list; verify gnomAD INFO tag names with
#     `bcftools view -h $GNOMAD_SITES | grep INFO`)

# 3. Provide a trio manifest (git-ignored): a TSV with header  trio_id  vcf  ped

# 4. Run the whole pipeline inside the container.
apptainer exec --cleanenv \
    --bind "$(dirname "$REF_FASTA")" --bind "$VEP_CACHE" --bind "$HPRV_WORK" \
    hprv.sif run_pipeline.sh --config config/config.yaml
```

Run a subset with `--from N --to M`. Every step is idempotent, so re-running resumes where it
stopped. On the laptop/dev path the same step scripts work through Docker automatically (the
container-exec layer auto-detects the runtime).

## Repository layout

```
docs/            source-cited methods reference + vetted pipeline design
config/          config.example.yaml (the contract; every tunable, no real paths)
env/             environment.yml (pinned conda toolchain layered onto the VEP image)
Dockerfile       one image: Ensembl VEP 115 base + bcftools/slivar/somalier/python...
pipeline/        step scripts (00..06) + run_pipeline.sh + lib/common.sh (container exec)
src/hprv/        shared python: config, annotations, genotype QC, ped
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
