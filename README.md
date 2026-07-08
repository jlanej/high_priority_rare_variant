# high_priority_rare_variant

Screen **GMKF Kids First per-trio VCFs** for high-priority **inherited** rare variants, and
consolidate **genes where rare functional variants recur across multiple individuals** — in
rare disease and germline pediatric cancer.

The inputs are per-trio VCFs (GRCh38) produced by GATK's genotype-refinement workflow,
**not** jointly genotyped across the cohort. Everything runs from **one container** under
Apptainer on HPC, driven by a single config file.

**Scope.** This pipeline focuses on **inherited germline variation** — dominant (a rare
functional heterozygous variant that recurs across individuals), recessive (homozygous /
compound-het-in-trans), and X-linked. **De novo** filtering and review, and **mtDNA
heteroplasmy**, are handled by separate dedicated machinery; de novo variants are detected here
only as a lightweight cross-reference, and mtDNA is out of scope.

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
| 0 | Per-trio QC gate (Mendelian error + chrX sex + contamination: verifyBamID FREEMIX, else VCF-only CHARR) | `qc_report.tsv` |
| 1 | Subset to trio members, normalize, build a **site-only union** of loci (never a genotype merge) | `cohort.sites.vcf.gz` |
| 2 | Annotate the union **once** (VEP + LOFTEE + dbNSFP + SpliceAI + gnomAD faf95 + ClinVar) — **VEP is never run per trio** | `cohort.sites.annotated.vcf.gz` |
| 3 | Select biologically-plausible sites (rarity + function; ClinVar P/LP override); tag each with *why* it was kept | `plausible.sites.vcf.gz` |
| 4 | Recover **real per-trio genotypes** at plausible sites + transfer annotations | per-trio `*.candidates.annotated.vcf.gz` |
| 5 | Pedigree-aware inheritance screen + genotype QC: **dominant** (inherited het), recessive (hom / comp-het-in-trans), X-linked; de novo is secondary | `candidates.calls.tsv` |
| 6 | **Cross-pedigree gene consolidation**: tally distinct individuals per gene by model (dominant het / biallelic / X-linked), weighted by constraint | `genes.ranked.tsv` |
| 7 | Consolidated **.xlsx** supplemental-table summary (documented: gene consolidation, calls, resolution, QC, audit) | `hprv_summary.xlsx` |
| 8 | **igv.js** trio variant-review export: `variants.tsv` + mini-CRAM slices (child/mother/father) + per-trio VCF tracks | `igv/` |

Every step records input/output counts and funnel tallies to `audit/counts.tsv`, assembled into
`audit/summary.md` — a global + per-trio "what went where and why" (see [Auditing](#auditing)).

The methodology — thresholds, tool choices, and the evidence behind them — is documented and
source-cited in **[docs/](docs/README.md)**. Every default lives in one place:
**[Canonical defaults](docs/README.md#canonical-defaults)**.

## Methods summary

What the pipeline actually does, in enough detail to follow it. Every threshold named here is a
configurable default from the **[Canonical defaults](docs/README.md#canonical-defaults)** table;
the evidence behind each choice is in **[docs/](docs/README.md)**.

**1. Cohort site list (no internal frequencies).** Each trio VCF is subset to its 3 members,
kept to `PASS` sites, normalized (`bcftools norm -m- -f`, split multiallelics + left-align), and
reduced to variant *loci only* (`view -G`). The per-trio site files are unioned
(`concat -a -D` + sort + dedup) into one cohort site list. This is a *union of loci*, never a
genotype `merge`: because the trios are not jointly genotyped, an absent record ≠ hom-ref, so any
internal cohort AC/AN would be fiction. The trios' stale embedded annotations (old VEP / gnomAD)
are stripped here and re-computed fresh.

**2. Annotate once.** The cohort site list is annotated a single time (VEP is *never* run per
trio): VEP consequence/IMPACT + LOFTEE (HC/LC pLoF), dbNSFP predictors (REVEL, AlphaMissense,
MPC), SpliceAI, CADD, plus external **gnomAD v4.1** frequency and dated **ClinVar** transferred
by `bcftools annotate`.

**3. Frequency oracle.** Rarity is judged on **gnomAD v4.1 group-max `faf95`** (the filtering
allele frequency — a CI lower bound), not internal counts. Benign-common variants
(`faf95 ≥ 0.05`, ClinGen BA1) are dropped and never rescued.

**4. Plausible-variant selection.** An inheritance-agnostic filter keeps a site if it is rare
(permissive-union cutoff) **and** functionally credible — HIGH/MODERATE VEP impact, LOFTEE-HC
pLoF, calibrated missense (REVEL / AlphaMissense), or SpliceAI — with **ClinVar P/LP** (≥1★) kept
as an override. Each kept site is tagged with *why* (`hprv_keep_reason`). Gene lists and
constraint are **not** applied here (never-drop rule), so novel genes survive.

**5. Per-trio inheritance screen (inherited focus).** Real per-trio genotypes are recovered at
the plausible sites and classified with refined-`GQ` genotype QC (GQ ≥ 20, DP ≥ 10, allele
balance bands from `AD`):
- **Dominant** — a rare (`faf95 < 1e-4`), functional **heterozygous** variant transmitted from
  ≥ 1 parent (origin recorded). This is the signal Step 6 consolidates.
- **Recessive** — homozygous, or **compound het in trans** (parent-of-origin: maternal + paternal).
- **X-linked recessive** — male hemizygous with a carrier mother (sex-aware ploidy; kid sex
  inferred from chrX heterozygosity when the PED is unknown).
- **De novo** — detected via GATK `hiConfDeNovo` (child-membership checked) but treated as a
  *secondary cross-reference* only; dedicated de novo filtering/review lives in separate machinery.

**6. Cross-pedigree gene consolidation.** Candidate calls are aggregated per gene into a count of
**distinct individuals** carrying a qualifying variant under each model (dominant het / biallelic
/ X-linked; de novo counted separately as secondary). A gene is flagged **recurrent** at
≥ `min_carriers` (default 2) distinct individuals, and genes are ranked recurrent-first, weighted
by **gene constraint** (LOEUF / pLI / s\_het) — a recurrent het in a haploinsufficient gene is far
more compelling than one in a constraint-tolerant gene. An optional de novo Poisson enrichment vs
a Samocha mutation model is reported as a secondary column when a mutation-rate table is supplied.

**7. Outputs for review.** A single documented **`.xlsx`** workbook consolidates the run
(gene consolidation, candidate calls, trio resolution, QC, audit) as a supplemental table. An
**igv.js** export produces a `variants.tsv` (the fork's variant-review schema — `chrom/pos/ref/alt`
+ inheritance mode + genotypes + our annotations as filterable columns), plus **mini-CRAM** slices
(±1 kb) for child/mother/father around each candidate locus (from a `sample→CRAM` map) and per-trio
VCF tracks — ready to serve with the [jlanej/igv.js](https://github.com/jlanej/igv.js) trio
variant-review server.

## Design principles

- **Focus on inherited variation; dominant recurrence is a first-class signal.** Heterozygous
  variants become interesting when they *stack up across individuals* in the same gene. De novo
  and mtDNA are handled by separate dedicated pipelines.
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
pipeline/        resolve_trios.py + step scripts (00..08) + run_pipeline.sh + lib/common.sh
src/hprv/        shared python: config, annotations, genotype QC, ped, selection, audit, report, igv
.github/workflows build + publish to GHCR on every commit (provenance + SBOM)
```

## Scope boundaries & known limitations

**Handled by separate dedicated pipelines (out of scope here):** **de novo** variant filtering
and review (detected here only as a lightweight cross-reference), and **mtDNA heteroplasmy**.

**Known limitations of this pipeline:** SNV/indel only — **CNV/SV are a real blind spot**
(10–15% of pediatric-cancer/rare-disease diagnoses); pseudogene/seg-dup regions (*PMS2*,
*CYP21A2*, *SMN1*) are low-confidence from short reads; the phenotype (Exomiser/HPO) prior is
planned. The pipeline is exercised end-to-end by an **integration test** on generated mock data
(`tests/integration/`, real bcftools) but is **not yet validated on real data** — GIAB/CMRG truth
sets and a positive-control panel are the next step. See
[pipeline_design.md](docs/pipeline_design.md#known-scope-limitations-stated-honestly-not-hidden).

## License

MIT — see [LICENSE](LICENSE).
