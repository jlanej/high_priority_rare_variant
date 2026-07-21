# Containerization, CI & Reproducibility

How this pipeline packages its tools, builds and publishes a signed container in CI, and runs it byte-reproducibly under Apptainer on HPC while keeping the public repository free of paths and PHI.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- **One image**, not many: `FROM ensemblorg/ensembl-vep:release_115.0` (the group's validated VEP base image) with a micromamba layer (`env/environment.yml`) for the CLI/Python tools. **conda-lock** hash-pinning for byte-identical rebuilds is a **TARGET** — today `environment.yml` pins versions, not hashes (see [Base image and dependency pinning](#base-image-and-dependency-pinning)).
- **Pinned tool set:** bcftools / htslib / samtools **1.23** (not 1.22 — slivar 0.3.4 on bioconda needs htslib ≥ 1.23.1), bedtools **2.31.1**, vcfanno **0.3.3**, slivar **0.3.4**, whatshap **2.3**, somalier **0.2.19**, plus Python `cyvcf2` / `pysam` / `pandas` / `numpy` / `scipy` — all version-pinned in `environment.yml`. **ensembl-vep 115** comes from the **base image**, not conda.
- **VEP cache is external:** release-matched to the VEP binary, bind-mounted at runtime, **never baked into the image** (species/assembly/release-specific, multi-GB).
- **CI:** buildx → **GHCR** on each commit, tagged by git SHA + branch, `provenance: true` + `sbom: true` and `actions/attest-build-provenance` (SLSA Build L3). **amd64-only.** Consumers pull by `@sha256:` digest, never `:latest`.
- **Apptainer:** point `APPTAINER_TMPDIR`/`APPTAINER_CACHEDIR` at real node-local scratch; **avoid `--containall`** (tmpfs `/tmp` OOMs heavy VEP/bcftools sorts); `--cleanenv`; bind read-only refs.
- **Repo hygiene:** config-driven paths + env vars, `.gitignore` every `*.vcf*` / `*.bam` / `*.cram` / `*.ped` / results dir; Kids First data is dbGaP controlled-access and is never committed.
- **Version-pin annotation resources too** — the VEP cache release, gnomAD v4.1, and the ClinVar release date are as load-bearing for reproducibility as the tool versions.

## Why reproducibility is non-negotiable here

This pipeline makes **variant-tiering decisions** on controlled-access pediatric data. Two builds a year apart must produce the same call on the same VCF, and a reviewer must be able to trace any reported variant back to an exact tool + annotation-resource stack. That requires pinning three separate layers, each of which can silently drift:

1. **Tool binaries** — a bcftools or VEP minor bump can change normalization or annotation output.
2. **Annotation resources** — VEP cache release, gnomAD version, ClinVar release date. These change more often than the tools and are the usual source of "it moved and nobody noticed."
3. **The container itself** — a floating base-image tag or an unlocked `environment.yml` re-solves to different package files on every rebuild.

## Base image and dependency pinning

The tension: **VEP is heavy** (Perl + Ensembl API + BioPerl + external caches) while the CLI tools (bcftools/htslib/samtools/bedtools/vcfanno/slivar) and the Python environment are conda-friendly. Three patterns are viable.

| Pattern | Reproducibility | When to use |
| --- | --- | --- |
| **`FROM ensemblorg/ensembl-vep` + micromamba layer (this project's choice)** | Medium — VEP exactly as upstream ships + a version-pinned conda layer for the CLI/Python tools; two package managers to pin | **Chosen here:** VEP's Perl/plugin stack is finicky to reproduce on conda, so it is taken from the group's validated upstream image |
| Single conda/micromamba image | High — one lockfile, one package manager | A pure-conda toolchain with no heavy Perl/plugin dependency |
| biocontainers, one image per tool | High per-tool (each is one bioconda package, versioned + hashed) | If a **workflow manager** pulls one container per process |

**What this project builds: `FROM ensemblorg/ensembl-vep:release_115.0` + a micromamba layer.** VEP's Perl + Ensembl API + plugin stack is finicky to reproduce on conda (bioconda's `perl` even shadows VEP's — the Dockerfile deletes it after the solve so `env perl` resolves to the base image's Perl), so VEP comes from the group's validated upstream image and the CLI/Python tools are layered on from `env/environment.yml`. Two hardening steps are **TARGET**, not yet done: pin the base image **by `@sha256:` digest** (the Dockerfile currently pins the tag `release_115.0`, and CI records the resolved digest for consumers to pull by), and generate a **conda-lock** lockfile so the micromamba layer installs byte-identical package files — a plain version string does *not* guarantee the identical artifact; conda-lock pins each package with a hash. Both are tracked in `env/environment.yml` and CLAUDE.md's Open TODOs.

Adopt **biocontainers** only if you move to a workflow manager (see below), where one container per process is the natural granularity.

### Keep the VEP cache OUT of the image

The VEP cache is species/assembly/release-specific and multi-GB. Baking it into the image bloats every layer and couples the cache release to the image tag. Instead:

- Version the cache **explicitly** and keep the **cache release matched to the VEP binary release** (VEP 115 → VEP 115 cache).
- Store it on the reference filesystem and **bind-mount it at runtime** (see Apptainer section).
- The same discipline applies to gnomAD v4.1 (the [frequency oracle](allele_frequency.md)) and the dated ClinVar release used by [clinical classification](clinical_classification.md): pin them, don't float them.

### Exact tool versions to pin

| Tool | Pinned version | Note |
| --- | --- | --- |
| bcftools / htslib / samtools | **1.23** | `env/environment.yml` pins **1.23** (not 1.22): slivar 0.3.4 on bioconda is built against htslib ≥ 1.23.1, and bcftools/samtools 1.23 resolve htslib to the same 1.23.1. Keep all three htslib-linked tools on the same release. |
| bedtools | **2.31.1** | Latest stable (Nov 2023); fixes GCC-13 builds. |
| vcfanno | **0.3.3** | — |
| slivar | **0.3.4** | Ships prebuilt static binaries + a bioconda recipe; pin the exact tag. |
| ensembl-vep | **115** | From the **base image** (`ensemblorg/ensembl-vep:release_115.0`), **not** conda. Do **not** float VEP — annotation output changes between releases. (Ensembl 116 exists as of June 2026; 115 is the pinned default.) |
| Python: `cyvcf2`, `pysam`, `pandas`, `numpy`, `scipy` | version-pinned in `environment.yml` | Link `pysam`'s bundled htslib to the **same 1.23 line** as the CLI tools. |

> **Related tooling versions (pinned elsewhere in the reference):** phenotype ranking uses **Exomiser 15.1.0** (June 2026; 15.0.0 was Feb 2026), which requires **Java 21** — pin the binary together with its *matching* Exomiser data release (do not pair a 15.x binary with the older `2406` data bundle, which targets Exomiser 14). See [gene lists & phenotype](gene_lists_and_phenotype.md). VEP r113 (Oct 2024) already updated its built-in gnomAD annotations to v4.1; if you rely on VEP's bundled gnomAD, note the release alignment.

## Publishing to GHCR via GitHub Actions

Build and push on each commit with `docker/setup-buildx-action` + `docker/build-push-action`. Log into `ghcr.io` with the built-in `GITHUB_TOKEN`; grant the workflow `packages: write`, `id-token: write`, and `attestations: write`.

```yaml
# .github/workflows/build.yml (illustrative, parameterized — no real org/paths)
permissions:
  contents: read
  packages: write
  id-token: write        # OIDC token for keyless signing
  attestations: write
jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@<pinned-sha>
      - uses: docker/setup-buildx-action@<pinned-sha>
      - uses: docker/login-action@<pinned-sha>
        with:
          registry: ghcr.io
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}
      - id: build
        uses: docker/build-push-action@<pinned-sha>
        with:
          push: true
          platforms: linux/amd64          # amd64-only; see below
          tags: |
            ghcr.io/<org>/<pipeline>:${{ github.sha }}
            ghcr.io/<org>/<pipeline>:${{ github.ref_name }}
          provenance: true                 # SLSA provenance
          sbom: true                       # embed SBOM
      - uses: actions/attest-build-provenance@<pinned-sha>
        with:
          subject-name: ghcr.io/<org>/<pipeline>
          subject-digest: ${{ steps.build.outputs.digest }}
          push-to-registry: true
```

- **Provenance + SBOM:** `provenance: true` and `sbom: true` on build-push-action, plus `actions/attest-build-provenance` (GA since May 2024). The OIDC token is signed by Sigstore Fulcio and satisfies **SLSA Build L3**.
- **Multi-arch:** build **`linux/amd64` only** unless you genuinely target ARM HPC — slivar/VEP ARM support is uneven and multi-arch doubles build time. Most HPC is amd64; single-arch is the defensible default.
- **Digest pinning downstream:** tag with the immutable **git SHA** (and optionally a semver release), record the resulting `sha256:` digest, and have consumers (and the pipeline itself) pull `ghcr.io/<org>/<pipeline>@sha256:...` — **never `:latest`**.
- **CI as a reproducibility test:** run the pipeline in CI against a **de-identified toy trio** committed for that purpose (see repo hygiene). This doubles as the smoke test that the locked environment still resolves and the tools still run.

## Running under Apptainer / Singularity on HPC

Convert the OCI image to a `.sif` from the digest-pinned reference, then exec against bind-mounted refs.

```bash
# Convert once (per digest). Point temp/cache at real node-local scratch, NOT $HOME or tmpfs.
export APPTAINER_TMPDIR="$SCRATCH/apptainer/tmp"
export APPTAINER_CACHEDIR="$SCRATCH/apptainer/cache"
apptainer pull pipeline.sif \
  docker://ghcr.io/<org>/<pipeline>@sha256:<digest>

# Exec: bind real-disk /tmp and the read-only reference bundle; clean the host env.
apptainer exec \
  --cleanenv \
  --bind "$SCRATCH/tmp:/tmp" \
  --bind "$REF_DIR:/refs:ro" \
  --bind "$WORK:/work" \
  pipeline.sif \
  vep --offline --cache --dir_cache /refs/vep_cache_115 --assembly GRCh38 \
      --input_file /work/trio.vcf.gz --output_file /work/trio.vep.vcf.gz
```

- **TMPDIR/CACHEDIR to scratch.** `.sif` conversion is I/O- and temp-heavy. Set `APPTAINER_TMPDIR` and `APPTAINER_CACHEDIR` to large **node-local scratch** — not `$HOME` (quota) and **not a tmpfs `/tmp`** (tmpfs shares RAM and OOMs on large-layer conversion).
- **Avoid `--containall`.** It isolates `/tmp` into an in-memory tmpfs; VEP and bcftools sort spill large temp files there → OOM. Prefer default containment, or explicitly `--bind $SCRATCH/tmp:/tmp` to real disk.
- **Bind the read-only reference bundle** (FASTA, the **release-matched** VEP cache, gnomAD/annotation files) and the data dir. Apptainer auto-binds `$HOME`, `$PWD`, `/tmp`, `/proc`, `/sys` by default.
- **`--cleanenv`** stops host environment variables from leaking in and perturbing reproducibility.

## Config-driven, PHI-free public repository

Kids First VCFs are **controlled-access (dbGaP)**. Nothing patient-derived may enter git.

- **No hardcoded paths, no sample IDs, no PHI.** All reference/annotation/output locations come from a config file (YAML) + env vars (e.g. `$REF_DIR`, `$VEP_CACHE`, `$WORK`).
- **`.gitignore`** must exclude `*.vcf*`, `*.bam`, `*.cram`, `*.ped` (real IDs), results directories, and `.env`.
- **Commit only:** a `config/config.example.yaml` and a **de-identified toy trio** used exclusively for CI.
- Every threshold in this reference lives in `config/config.example.yaml` as an **overridable default** — a gene-specific ClinGen VCEP value overrides a generic cutoff.

## Do you need a workflow manager?

For a **single-purpose screen** (annotate + filter a handful of per-trio VCFs), a plain, well-pinned container + a driver script is defensible and simpler to audit. Adopt a manager when you need per-sample parallelism across many trios, resumability, and provenance.

| Option | Fit for this project | Notes |
| --- | --- | --- |
| **Nextflow / nf-core** | Best fit if scaling | Native Apptainer support; one container per process (pairs with biocontainers); Kids First / GA4GH ecosystems lean Nextflow/WDL; nf-core gives reproducibility scaffolding for free. |
| **WDL / Cromwell** | Choose for consortium sharing | GMKF/Kids First DRC and Broad ecosystem is WDL-native; pick if running on Cavatica/Terra. Call-caching gives resumability. |
| **Snakemake** | Fine if Python-first, single cluster | — |

**Recommendation:** start scriptable + containerized; graduate to **Nextflow (Apptainer profile)** or **WDL** if you scale to the full cohort or submit to Cavatica.

## The non-joint, genotype-refined input constrains what "reproducible" means

These are **per-trio, GATK genotype-refinement** VCFs (posteriors added; **not** cohort joint-genotyped). Two reproducibility-relevant consequences:

- **Cohort-internal allele frequencies are uninterpretable for filtering** — absent genotype ≠ hom-ref across a non-joint merge — so rarity filtering relies on **external** population AF (gnomAD v4.1). Pinning that external resource is therefore part of reproducibility, not an afterthought. See [cohort construction](cohort_construction.md) and [allele frequency](allele_frequency.md).
- **GATK refinement tags are inputs, not something the pipeline recomputes.** Use `hiConfDeNovo` as the primary de novo screen and re-verify per the QC gates. One knob is genuinely resource-dependent and must be documented per run: GATK's `--num-reference-samples-if-no-call` value depends on the `--supporting-callsets` resource shipped. The commonly used `af-only-gnomad.hg38.vcf.gz` is a gnomAD **exome** AF resource whose N is not "~76k genomes"; set N to the shipped resource's documented sample size and record it as version-dependent. Details in [inheritance & genotype QC](inheritance_and_genotype_qc.md).

## Known scope limitations affecting reproducibility

Reproducibility of a *decision* is only as good as the pipeline's coverage. State these honestly:

- **SNV/indel only initially.** CNV/SV are a real blind spot — 10–15% of pediatric-cancer and rare-disease diagnoses are CNV/SV (single-exon *RB1*/*SMARCB1*/*DICER1*/*NF1* deletions, *PMS2* rearrangements). A reproducible SNV/indel run still systematically misses these until a GATK-gCNV / Manta / ExomeDepth module is added.
- **Pseudogene / segmental-duplication genes** (*PMS2*/*PMS2CL*, *CYP21A2*, *SMN1/2*, *NEB*, *GBA*) are low-confidence from short reads regardless of how well the container is pinned; flag those regions.
- **Proband post-zygotic mosaicism** (*NF1*, overgrowth) — low-VAF calls fall outside the het AB band and need a dedicated mosaic tier.
- **No calibration/validation harness yet.** Reproducible tiering decisions should be measured, not assumed: add GIAB/CMRG truth sets, synthetic-diploid benchmarking, a synonymous-λ ≈ 1 check, and a positive-control variant panel to quantify de novo and prioritization sensitivity/precision. Pin the truth-set versions like any other resource.

## Recommended defaults (this pipeline)

| Item | Canonical default |
| --- | --- |
| Base image | `FROM mambaorg/micromamba@sha256:...` (digest-pinned) + **conda-lock** lockfile |
| bcftools / htslib / samtools | **1.22** (overridable) |
| bedtools | **2.31.1** |
| vcfanno / slivar | **0.3.3** / **0.3.4** |
| ensembl-vep | **115** (matches working annotate script) |
| Python | `cyvcf2` / `pysam` (htslib 1.22) / `pandas` / `numpy`, all lock-pinned |
| VEP cache | **external, bind-mounted, release-matched** to the VEP binary; annotate against **GRCh38** |
| Frequency oracle | **gnomAD v4.1** (external only), read from the **VEP cache**; filter field is the grpmax **proxy** (point estimate). `faf95` = TARGET — the cache carries no AC/AN ([limitations.md §2](limitations.md)) |
| CI / publish | buildx → **GHCR** per commit; tag by git SHA + branch; `provenance: true` + `sbom: true` + `attest-build-provenance`; **amd64-only**; consume by `@sha256:` digest |
| Apptainer runtime | `APPTAINER_TMPDIR`/`CACHEDIR` → node scratch; **no `--containall`**; `--bind $SCRATCH/tmp:/tmp` + refs; `--cleanenv` |
| Repo hygiene | config-driven paths + env vars; `.gitignore` all VCF/BAM/CRAM/PED/results; dbGaP data never committed |
| Resource pinning | VEP cache release, gnomAD v4.1, dated ClinVar release, Exomiser 15.1.0 + matching data release — all version-pinned |

All values are **configurable defaults** in `config/config.example.yaml`, not immutable law; a gene-specific ClinGen VCEP value overrides any generic cutoff.

## Sources

- Ensembl releases / VEP: <https://github.com/Ensembl/ensembl-vep/releases> ; Ensembl 115 (Sep 2025) <https://www.ensembl.info/2025/09/02/ensembl-115-has-been-released/> ; VEP r113 → gnomAD v4.1 <https://www.ensembl.info/2024/10/18/ensembl-113-has-been-released/>
- bcftools / htslib / samtools releases: <https://github.com/samtools/bcftools/releases> ; <http://www.htslib.org>
- bedtools 2.31.1: <https://github.com/arq5x/bedtools2/releases/tag/v2.31.1>
- slivar releases: <https://github.com/brentp/slivar/releases> ; vcfanno releases: <https://github.com/brentp/vcfanno/releases>
- Exomiser releases (15.1.0, June 2026; Java 21): <https://github.com/exomiser/Exomiser/releases> ; data-version compatibility: <https://github.com/exomiser/Exomiser/discussions/562>
- gnomAD v4.1 (Apr 2024; 730,947 exomes + 76,215 genomes, GRCh38): <https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/>
- GATK Genotype Refinement workflow: <https://gatk.broadinstitute.org/hc/en-us/articles/360035531432-Genotype-Refinement-workflow-for-germline-short-variants>
- GMKF / Kids First harmonization (BWA-MEM, GRCh38, GATK best practices): <https://aacrjournals.org/cancerres/article/79/13_Supplement/2465/634458>
- micromamba-docker (FROM by sha256 digest; lockfile create): <https://micromamba-docker.readthedocs.io/en/stable/advanced_usage.html> ; conda-lock (per-package checksums): <https://github.com/conda/conda-lock>
- GitHub Actions buildx / GHCR: <https://github.com/docker/build-push-action> ; SLSA provenance / SBOM attestations: <https://docs.docker.com/build/ci/github-actions/attestations/> ; SLSA provenance metadata: <https://docs.docker.com/build/metadata/attestations/slsa-provenance/>
- Apptainer TMPDIR/CACHEDIR, bind, tmpfs `/tmp`, `--containall`: <https://apptainer.org/docs/user/main/> ; USC CARC <https://www.carc.usc.edu/user-guides/hpc-systems/software/apptainer> ; UCL <https://www.rc.ucl.ac.uk/docs/Software_Guides/Singularity/>
- Workflow-manager comparison (Nextflow/nf-core vs Snakemake vs WDL/Cromwell): <https://www.nature.com/articles/s41598-021-99288-8> ; nf-core <https://nf-co.re>
