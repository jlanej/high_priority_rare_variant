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

> 📋 **VEP-centric contract — know what the screen can and cannot see *before* you run it.**
> Almost every annotation comes from one source: a **VEP 115 GRCh38 cache + the CADD, SpliceAI,
> REVEL, AlphaMissense, NMD and SpliceVault plugins** (SpliceAI and REVEL/AlphaMissense are required
> by default; SpliceVault is optional review evidence — what the cell does when a splice site is lost).
> Exactly **two** things are `bcftools annotate`-transferred, because the cache cannot supply them
> at any price: the **ClinVar sites VCF** (review status ⇒ gold stars, which rank in Step 9 and
> never gate) and the **gnomAD v4.1 joint slim** (real **`faf95`** + **`nhomalt`**). The slim is
> **required by the default configuration**: `resources.gnomad.oracle: faf95` makes faf95 the ONE
> rarity oracle for the whole run and halts at preflight without it; `grpmax_proxy` is the
> deliberate opt-down to the VEP-cache point estimate, and the two arms never mix within a run.
> What remains missing: **no LOFTEE** (pLoF confidence), no exome/genome discordance flag. The full
> ledger, with the cost to close each, is **[docs/limitations.md](docs/limitations.md)**. Read it
> before you interpret a negative result.

## What it does

A resolve preflight + ten-step flow (Steps 0–9; see **[docs/pipeline_design.md](docs/pipeline_design.md)**
for the vetted design and the artifact each step produces):

| Step | What | Output |
|------|------|--------|
| resolve | Map each `kid/dad/mom` trio to the VCF containing all three members (exact sample-ID match; extras OK); generate PEDs | `trios.resolved.tsv`, `trio_resolution.tsv`, `peds/` |
| 0 | Per-trio QC gate (Mendelian error on the first `qc.max_sites` QC-passing autosomal biallelic sites + chrX sex + contamination: verifyBamID FREEMIX, else VCF-only CHARR) | `qc_report.tsv` |
| 1 | Subset to trio members, normalize, build a **site-only union** of loci (never a genotype merge) | `cohort.sites.vcf.gz` |
| 2 | Annotate the union **once** (VEP 115 cache + CADD/SpliceAI/REVEL/AlphaMissense plugins; gnomAD v4.1 AFs and ClinVar `CLIN_SIG` ride in the cache), then two `bcftools annotate` transfers: ClinVar review status, and the gnomAD joint slim for `faf95`/`nhomalt` (required under the default `faf95` oracle) — **VEP is never run per trio** | `cohort.sites.annotated.vcf.gz` |
| 3 | Select biologically-plausible sites (rarity + function; ClinVar P/LP override); tag each with *why* it was kept | `plausible.sites.vcf.gz` |
| 4 | Recover **real per-trio genotypes** at plausible sites + transfer annotations | per-trio `*.candidates.annotated.vcf.gz` |
| 5 | Pedigree-aware inheritance screen + genotype QC: **dominant** (inherited het), recessive (hom / comp-het-in-trans), X-linked; de novo is secondary. Each call carries HGVS, exon/CDS geometry, the NMD verdict, the decomposed SpliceAI event and every INFO field verbatim | `candidates.calls.tsv` |
| 5b | **Wide-window SpliceAI** over the called set (live `-D 4999`, one chunk per SLURM array task, scores cached per variant): `spliceai_wide_*` evidence columns beside the precomputed score the gate used | `candidates.calls.tsv` (+ `spliceai_rescore/scores.tsv`) |
| 6 | **Cross-pedigree gene consolidation**: tally distinct individuals per gene by model (dominant het / biallelic / X-linked), weighted by constraint; a case-only recurrence rank plus, with a mutational-target table, the expected-carrier excess (`p_carrier_excess`) | `genes.ranked.tsv` |
| 7 | Consolidated **.xlsx** supplemental-table summary (documented: gene consolidation, calls, resolution, QC, audit) | `hprv_summary.xlsx` |
| 8 | **igv.js** trio variant-review export: `variants.tsv` + mini-CRAM slices (child/mother/father) + per-trio VCF tracks | `igv/` |
| 9 | **Prioritization**: gene excess over its mutational target (NB2, trimmed fit, mid-p calibration) + a six-signal artifact panel + a graded gene down-weight, then per-variant tiering and an additive `priority_points` composite — **a re-rank, never a drop** | `variants.prioritized.tsv`, `genes.prioritized.tsv`, `igv/variants.prioritized.tsv` (the igv.js review table: Step 8's columns + every triage column, so filtering/sorting happens in the review tool) |

Every step records input/output counts and funnel tallies to `audit/counts.tsv`, assembled into
`audit/summary.md` — a global + per-trio "what went where and why" (see [Auditing](#auditing)).

The methodology — thresholds, tool choices, and the evidence behind them — is documented and
source-cited in **[docs/](docs/README.md)**. Every default lives in one place:
**[Canonical defaults](docs/README.md#canonical-defaults)**.

## Methods summary

> **Publication-quality Methods:** [docs/methods.md](docs/methods.md) (rendered as
> [docs/methods.pdf](docs/methods.pdf)) is the full, code-derived description of every
> filter, model and default — written as a Methods section, with the thresholds, resource
> versions and per-step outputs tabulated. The paragraphs below are the short form.

What the pipeline actually does, in enough detail to follow it. Every threshold named here is a
configurable default from the **[Canonical defaults](docs/README.md#canonical-defaults)** table;
the evidence behind each choice is in **[docs/](docs/README.md)**.

**1. Cohort site list (no internal frequencies).** Each trio VCF is subset to its 3 members,
kept to FILTER `PASS`/`.` sites, normalized (`bcftools norm -m- -f`, split multiallelics + left-align), and
reduced to variant *loci only* (`view -G`). The per-trio site files are unioned
(`concat -a -D` + sort + dedup) into one cohort site list. This is a *union of loci*, never a
genotype `merge`: because the trios are not jointly genotyped, an absent record ≠ hom-ref, so any
internal cohort AC/AN would be fiction. The trios' stale embedded annotations (old VEP / gnomAD)
are stripped here and re-computed fresh.

**2. Annotate once.** The cohort site list is annotated a single time (VEP is *never* run per
trio): VEP 115 (consequence/IMPACT/SYMBOL/HGVS/MANE) + the **CADD**, **SpliceAI**, **REVEL** and
**AlphaMissense** plugins, with gnomAD v4.1 per-population point AFs (`--af_gnomade`/`--af_gnomadg`)
and ClinVar `CLIN_SIG` (`--check_existing`) coming out of the cache itself. The CSQ fields are
lifted to INFO with a `vep_` prefix (`bcftools +split-vep`). Exactly **two** external transfers
follow, each under its own INFO namespace: the ClinVar sites VCF (`clinvar_CLNREVSTAT` ⇒
`clinvar_stars`) and the gnomAD v4.1 joint slim (`gnomad_faf95`, `gnomad_faf95_group`,
`gnomad_nhomalt`, `gnomad_AF_joint`, `gnomad_AF_grpmax`). Already have a VEP 115 VCF? Point
`resources.vep.annotated_vcf` at it and Step 2 skips the VEP call (the two transfers still run).

**3. Frequency oracle.** Rarity is judged on **gnomAD v4.1**, never on internal counts, through a
single chokepoint (`annotations.frequency()`) that reads **one oracle per run**, selected by
`resources.gnomad.oracle`:

- **`faf95` (the default)** — gnomAD's published filtering allele frequency
  (`fafmax_faf95_max_joint`), the *lower bound of the 95% CI* and the quantity ACMG/ClinGen specify
  for frequency filtering. Requires the gnomAD joint slim (`resources.gnomad.sites_slim`, ~10 GB,
  `--only gnomad_sites`); the run halts at preflight without it. An allele gnomAD has but published
  no faf95 for resolves to 0 (`rarity_basis=zero_ci`, rarest); an allele with no gnomAD record
  resolves to absent (rarest). The proxy is never consulted.
- **`grpmax_proxy`** — a deliberate opt-down: the max point-estimate AF over the **grpmax-eligible**
  ancestry groups (AFR/AMR/EAS/NFE/SAS), read from the VEP cache. No CI correction is possible, so
  it errs slightly toward *dropping* low-count alleles.

`rarity_oracle` is recorded once per run and `rarity_basis` (`measured`/`zero_ci`/`absent`) is the
per-variant provenance; every row's `rarity_af` is the value the gates actually applied. Because
faf95 ≤ the point estimate, the default **retains more** at the same cutoffs than the proxy would
([the ledger](docs/limitations.md#2-faf95--implemented-and-the-default-oracle-the-proxy-is-a-deliberate-opt-down)).
VEP's `MAX_AF` and the global AFs are **reporting only, never filter fields under either arm** —
they fail in opposite directions and neither is a safe substitute. Benign-common variants
(`≥ 0.05`, ClinGen BA1) are dropped and never rescued.

**4. Plausible-variant selection.** An inheritance-agnostic filter keeps a site if it is rare
(permissive-union cutoff) **and** functionally credible. The functional ladder is deliberately
**three rungs**: HIGH/MODERATE VEP impact, else **SpliceAI Δ ≥ 0.2**, else **CADD ≥ 25.3**. **ClinVar P/LP** (no conflicts) is
an override, and also rescues a variant from the rarity gate. Each kept site is tagged with *why*
(`hprv_keep_reason`). Gene lists and constraint are **not** applied here (never-drop rule), so
novel genes survive.

Two caveats a reader must hold onto. **CADD is the entire NON-SPLICE non-coding screen** (SpliceAI
covers the splice-disrupting class a rung above it) — every missense is
IMPACT=MODERATE and is kept at rung 1, so 25.3 (Pejaver-2022's PP3-supporting cutoff,
calibrated on *missense only*) is applied exclusively to non-coding variants it was never
calibrated for. Treat it as a discovery rank (≈ top 0.3% genome-wide), **not** as ACMG PP3
evidence. And the P/LP override is **star-blind by design** — review status is available
(`clinvar_stars`, from the ClinVar transfer) and ranks in Step 9, but gating the screen on it would
violate never-drop; this over-retains (more to curate) rather than over-drops.

**5. Per-trio inheritance screen (inherited focus).** Real per-trio genotypes are recovered at
the plausible sites and classified with refined-`GQ` genotype QC (GQ ≥ 20, DP ≥ 10, allele
balance bands from `AD`):
- **Dominant** — a rare (`dominant_max`, default `rarity_af` < 1e-4 on the run's oracle), functional **heterozygous** variant
  transmitted from ≥ 1 parent (origin recorded). This is the signal Step 6 consolidates.
- **Recessive** — homozygous, or **compound het in trans** (parent-of-origin: maternal + paternal).
- **X-linked recessive** — male hemizygous with a carrier mother (sex-aware ploidy; kid sex
  inferred from chrX heterozygosity when the PED is unknown).
- **De novo** — detected via GATK `hiConfDeNovo` (child-membership checked) but treated as a
  *secondary cross-reference* only; dedicated de novo filtering/review lives in separate machinery.
- **Splice events and transcript geometry** — every call carries the SpliceAI event behind its
  score (`spliceai_event`, the affected site's position, a paired event, the exon-boundary shift
  and its frame, `spliceai_effect`), VEP's exon/intron numbering and CDS position, MANE status and
  the NMD verdict (`nmd_status`: `escaping`/`triggering`/`not_assessed`, from the VEP NMD plugin);
  **Step 5b** re-scores every call live at a 4,999 bp window (`spliceai_wide_*`), the class the
  50 bp precomputed set cannot see. Evidence for review — the gate keeps reading the precomputed max.

**6. Cross-pedigree gene consolidation.** Candidate calls are aggregated per gene into a count of
**distinct individuals** carrying a qualifying variant under each model (dominant het / biallelic
/ X-linked; de novo counted separately as secondary). A gene is flagged **recurrent** at
≥ `min_carriers` (default 2) distinct individuals, and genes are ranked recurrent-first, weighted
by **gene constraint** (LOEUF / pLI / s\_het / pHaplo) — a recurrent het in a haploinsufficient
gene is far more compelling than one in a constraint-tolerant gene. The per-model recurrence
p-values (`p_recurrence`, BH `q`, exome-wide flag) are a **case-only rank**, never a calibrated
test — two carriers of private variants at N = 200 already give p ≈ 3e-7. When a mutational-target
table is supplied, Step 6 also reports each gene's expected carrier count from its Samocha
mutation rate (`exp_carriers_mu`, `carrier_excess_ratio`, `p_carrier_excess`) and, with
`burden.rank_by_mutational_target: true` (the default), orders recurrent genes by that excess so
long genes stop leading by size; an optional de novo Poisson enrichment is a secondary column.

**7. Outputs for review.** A single documented **`.xlsx`** workbook consolidates the run
(gene consolidation, candidate calls, trio resolution, QC, audit) as a supplemental table. An
**igv.js** export produces a `variants.tsv` (the fork's variant-review schema — `chrom/pos/ref/alt`
+ inheritance mode + genotypes + our annotations as filterable columns), plus **mini-CRAM** slices
(±1 kb) for child/mother/father around each candidate locus (from a `sample→CRAM` map) and per-trio
VCF tracks — ready to serve with the [jlanej/igv.js](https://github.com/jlanej/igv.js) trio
variant-review server.

**9. Prioritization — a re-rank, never a drop.** Steps 6–8 leave the reviewer a flat list. Step 9
orders it, in two layers ([docs/prioritization.md](docs/prioritization.md)):

- **Gene layer — excess over the gene's mutational target.** `E_g = C·μ_g` from gnomAD v2.1.1's
  per-gene Samocha targets, with a **negative-binomial** null whose `(C, α)` are re-fit per cohort
  on an iteratively **trimmed** bulk, so the artifact tail cannot calibrate its own null. The NB is
  not an assumption: dispersion measured **φ = 18.7** raw and **1.29** trimmed, and in a two-fold
  cross-fit the Poisson tail was **2.41× anti-conservative** at α = 1e-3 (with its own arm
  trimmed; 1.82× untrimmed) while the NB was
  conservative — the correct direction of error under the never-drop rule. A **mid-p calibration
  diagnostic** for both nulls lands in the audit on every run (the pre-implementation science
  audit's **A-3** calibration gap). `C` is fit over the **full** gene universe including zero-count genes, because
  the candidate list is a zero-truncated sample.
- **Six orthogonal artifact signals** → an integer corroboration count: cohort saturation (**195×**
  enriched), segdup overlap (10.1×), artifact-prone gene family (6.5×), synonymous o/e departure
  (5.2×), gnomAD's own constraint flag (4.0×), low cumulative allele frequency (3.9×). The
  statistic says a gene is *anomalous*; the panel is the **independent second witness** required
  before any penalty.
- **A four-tier graded down-weight**, with an **auditable established-gene ceiling** (a
  control-union gene never enters T2/T3, and is flagged `established_gene_high_excess` instead) and
  a CDS-fallback ceiling. Measured at the shipped default (`min_n_for_ratio_rule: 3`): **92 genes /
  2,369 variants (9.57%) triaged at 100% established-gene retention** (228 / 2,565 / 10.36% with
  the count floor at 1). The maximum penalty (−3.0) **cannot by itself** demote a variant
  carrying strong molecular evidence in a constrained gene — it is a re-rank, not a veto.
- **Variant layer** — a V0–V5 tier under ACMG/ClinGen **SVI mechanism gating** (gene constraint
  counts only when the variant has a credible molecular effect; a molecularly-benign prediction
  caps the total regardless of the gene), plus an additive Tavtigian-style **`priority_points`**
  composite in which **every term is its own reported column**, so a reviewer reads
  `spliceai=+4, rarity=+2, constraint=+1, gene_artifact=−3` rather than one opaque number.
- **Two rankings always ship**: fully agnostic, and prior-informed with an optional phenotype
  overlay that **defaults OFF** (so hprv stays phenotype-agnostic). `rank_delta` exposes exactly
  which calls a gene list promoted — the set to scrutinise for confirmation bias.

Three honest limits carried in the output rather than papered over: **V5 is reached only by a
nonsense or frameshift the VEP NMD plugin assessed and did not flag as escaping**
(`nmd_status=triggering`; canonical splice and every unassessed pLoF stay V4), the missense tier reads REVEL, then
AlphaMissense, then CADD in a **fixed precedence** (`missense_evidence_source` names which spoke;
the CADD fallback is labelled `cadd_offlabel` and claims no graded strength), and
**`priority_points` is not an ACMG score** — never read against Tavtigian's P/LP/VUS bands.

## Design principles

- **Focus on inherited variation; dominant recurrence is a first-class signal.** Heterozygous
  variants become interesting when they *stack up across individuals* in the same gene. De novo
  and mtDNA are handled by separate dedicated pipelines.
- **gnomAD v4.1 is the only population-frequency oracle**, read as ONE quantity per run: real
  `faf95` from the gnomAD joint slim by default, or the VEP-cache grpmax point-estimate proxy as a
  deliberate opt-down (`resources.gnomad.oracle`; see
  [limitations](docs/limitations.md#2-faf95--implemented-and-the-default-oracle-the-proxy-is-a-deliberate-opt-down)).
  Because the trios are not jointly genotyped, internal cohort AC/AN is meaningless
  (absent ≠ hom-ref) and is used only as an artifact/blocklist signal.
- **Gene lists and constraint are priors/tiers, never hard filters** ("never-drop rule") — so
  novel-gene discovery survives.
- **One container, config-driven, no hard paths.** The same image runs on a laptop (Docker) and
  on HPC (Apptainer); scripts only ever call tools inside it.
- **Fail loudly, verify before claiming, be idempotent** (`.done` files; integrity checks).

## Quickstart

```bash
# 1. Get the image (built + published to GHCR on every commit). Pull by digest on HPC.
apptainer pull hprv.sif docker://ghcr.io/<owner>/high_priority_rare_variant:latest

# 2. Prepare the annotation resources ONCE (the image ships software; this fetches data).
#    A bare `fetch` prepares the default set: a reference FASTA, a VEP 115 GRCh38 cache, CADD
#    (license-gated, ~82 GB), the per-gene constraint + mutational-target tables (Steps 6/9),
#    ClinVar (~0.18 GB -> gold stars), and REVEL + AlphaMissense (~1.3 GB -> Step 9's calibrated
#    missense tier, REQUIRED by default). Do NOT narrow it with --only: that silently skips the
#    last three. Two things a bare `fetch` does NOT pull, both needed by the default config:
#      - the gnomAD v4.1 JOINT slim (~10 GB lands; preparing it streams ~877 GB) — REQUIRED by the
#        default `resources.gnomad.oracle: faf95` (opt down with `grpmax_proxy`);
#      - SpliceAI (step 2b below), whose useful files are login-gated — REQUIRED by default
#        (`resources.vep.spliceai_required: true` HALTS the run at preflight when missing).
#    LOFTEE/dbNSFP data are genuinely unused and stay opt-in. prepare_resources.sh + its pinned
#    manifest ship IN the image (on PATH). See docs/resources.md.
apptainer exec --bind /data hprv.sif \
    prepare_resources.sh --dir /data/hprv_resources --accept-license fetch
apptainer exec --bind /data hprv.sif \
    prepare_resources.sh --dir /data/hprv_resources --only gnomad_sites fetch

# 2b. SpliceAI raw hg38 scores. `--only spliceai` cannot finish the job: the SNV mirror it can
#     reach is MANE-only and the indel file is login-gated with no no-login mirror. Use the
#     BaseSpace helper (needs an authenticated `bs` CLI on the host). The `spliceai/` subdir is
#     load-bearing — emit-env exports $DIR/spliceai/<file>, which is where these land.
scripts/download_spliceai.sh --dir /data/hprv_resources/spliceai --ref "$REF_FASTA"

#    `verify` mirrors the default config: gnomad_sites, revel and alphamissense are reported as
#    REQUIRED (non-zero exit when missing) unless the config opts down.
apptainer exec --bind /data hprv.sif \
    prepare_resources.sh --dir /data/hprv_resources verify
apptainer exec --bind /data hprv.sif \
    prepare_resources.sh --dir /data/hprv_resources emit-env --out /data/hprv_resources/resources.env
#    (emit-env writes GNOMAD_SITES uncommented once the slim is present.)

#    Already have a VEP 115 GRCh38 VCF of your sites? Set resources.vep.annotated_vcf and Step 2
#    skips the VEP call entirely — no cache/CADD/plugin fetch needed. The two bcftools transfers
#    still run, so the gnomAD slim is still required under the default oracle.

# 3. Configure. Copy the example; point the ${ENV} placeholders at your prepared resources.
cp config/config.example.yaml config/config.yaml     # config.yaml is git-ignored
source /data/hprv_resources/resources.env            # exports REF_FASTA / VEP_CACHE / CADD_SNV / ...
export HPRV_WORK=/path/to/work
#    (resources.env also exports SPLICEAI_SNV / SPLICEAI_INDEL — these ARE read: Step 2 wires them
#     as the VEP SpliceAI plugin, and spliceai_required (default true) HALTS preflight if missing.)

# 4. Provide inputs (git-ignored):
#    - a trios file: TSV with a header naming kid/dad/mom (any order); IDs match the VCFs:
#          #kid   dad    mom
#          CH1    FA1    MO1
#      export TRIOS_FILE=/path/to/trios.tsv
#    - the VCF source (a directory and/or a list file):
#          export VCF_DIR=/path/to/vcfs        # globbed for *.vcf.gz/*.vcf/*.bcf
#    The pipeline finds the VCF containing all three members for each trio and generates
#    the internal manifest + PEDs automatically — sample order within a VCF does not matter,
#    and a VCF may contain additional members.

# 5. Run the whole pipeline inside the container.
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
pipeline/        resolve_trios.py + step scripts (00..09) + run_pipeline.sh + lib/common.sh
src/hprv/        shared python: config, annotations, genotype QC, ped, selection, prioritize, audit, report, igv
.github/workflows build + publish to GHCR on every commit (provenance + SBOM)
```

## Scope boundaries & known limitations

**Handled by separate dedicated pipelines (out of scope here):** **de novo** variant filtering
and review (detected here only as a lightweight cross-reference), and **mtDNA heteroplasmy**.

**Known limitations of this pipeline:** the full ledger — what the screen cannot see, why, and
what each item costs to fix — is **[docs/limitations.md](docs/limitations.md)**. It is the anchor;
the headlines are:

- **From the VEP-centric contract:** no LOFTEE (pLoF confidence / PVS1 grading) and no
  exome/genome discordance flag. `faf95`, `nhomalt`, ClinVar stars, REVEL/AlphaMissense and
  SpliceAI are all wired (the first three via the two `bcftools annotate` transfers, the rest as
  VEP plugins).
- **Structural, independent of the contract:** SNV/indel only — **CNV/SV are a real blind spot**
  (10–15% of pediatric-cancer/rare-disease diagnoses); pseudogene/seg-dup regions (*PMS2*,
  *CYP21A2*, *SMN1*) are low-confidence from short reads; the phenotype (Exomiser/HPO) prior is
  planned.
- **Not yet validated on real data.** The pipeline is exercised end-to-end by an **integration
  test** on generated mock data (`tests/integration/`, real bcftools), but sensitivity/precision
  are **unmeasured**, not measured-and-acceptable — GIAB/CMRG truth sets and a positive-control
  panel are the next step.

**Reading a negative result:** "no candidate" means no *coding*, SpliceAI-high or CADD-high
non-coding variant passing the run's rarity gate (gnomAD `faf95` by default; a point-estimate
proxy only if opted down) in a *SNV/indel* callset, without phenotype weighting or CNV calling, and
with ClinVar review status used only to rank, never to gate. That is a useful screen; it is not an
exclusion. See also
[pipeline_design.md](docs/pipeline_design.md#known-scope-limitations-stated-honestly-not-hidden).

## License

MIT — see [LICENSE](LICENSE).
