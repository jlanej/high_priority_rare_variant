# Functional Consequence & In-Silico Impact Prediction

Assigns molecular consequence and calibrated deleteriousness evidence to rare SNV/indels so downstream tiering can weigh loss-of-function, missense, and splice impact defensibly.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

> ### ⚠ Status: most of this document is reference science, not running code
>
> The pipeline runs a **VEP-only contract** — a VEP 115 GRCh38 cache plus the CADD plugin, and
> nothing else. No LOFTEE, dbNSFP, SpliceAI or ClinVar file is downloaded or bind-mounted.
>
> **What is IMPLEMENTED at the Step-3 screen is a two-rung ladder and nothing more:**
> 1. VEP `IMPACT` ∈ `keep_impacts` (`[HIGH, MODERATE]`) → keep;
> 2. else `CADD_PHRED ≥ 25.3` → keep; else drop `not_functional`.
> Plus a ClinVar `CLIN_SIG` P/LP override (**unstarred** — the cache carries no review status) and
> a BA1 drop at AF ≥ 0.05. See [`src/hprv/selection.py`](../src/hprv/selection.py).
>
> **Not available:** LOFTEE (HC/LC, flags), REVEL, AlphaMissense, MPC, MetaRNN, BayesDel, VEST4,
> SpliceAI, ClinVar review status. Every threshold below that names one of them is a **TARGET**
> for the tiering step in [ROADMAP.md](ROADMAP.md), not a filter that fires today. The science is
> retained deliberately: it is the justification for that roadmap.
>
> Full accounting of what the screen therefore cannot see, and what each gap costs to close:
> **[limitations.md](limitations.md)**. Canonical, marked IMPLEMENTED-vs-TARGET thresholds:
> [Canonical defaults](README.md#canonical-defaults). If either disagrees with this file, they win.

## TL;DR

- **IMPLEMENTED — the whole functional ladder:** **VEP IMPACT** (HIGH/MODERATE kept) → **CADD PHRED ≥ 25.3** → drop. One annotation source (VEP 115 cache), one plugin (CADD). Version-pin both.
- **IMPLEMENTED — ClinVar `CLIN_SIG` P/LP keeps a variant** and rescues it from the rarity gate, at **any** review status; no star gate exists (the cache has no `CLNREVSTAT`). This over-retains rather than over-drops — safe for a screen, costly in curation.
- **CADD is on an off-label threshold.** 25.3 is Pejaver-2022's PP3-supporting cutoff, calibrated on **missense only**; missense never reaches the CADD rung (it is MODERATE, kept a rung earlier), so 25.3 is applied *exclusively* to the non-coding variants it was never calibrated for. It is a **discovery rank (≈ top 0.3% genome-wide), not PP3 evidence**. See [CADD as the only predictor](#cadd-as-the-only-predictor--implemented-and-off-label).
- **VEP IMPACT** (HIGH/MODERATE/LOW/MODIFIER) is a coarse convenience tier only — never rely on it alone for splice-adjacent variants. With no SpliceAI to layer on, **that is exactly what the screen currently does**, and deep-intronic / exonic-synonymous splice variants are invisible to it.
- **TARGET — pLoF confidence = LOFTEE HC with no flags**; grade PVS1 via the Abou Tayoun 2018 decision tree (NMD-escape / last-exon / 3′-terminal-50 bp / single-exon downgrade), gated on ClinGen gene–disease validity ≥ Moderate and a known LoF mechanism.
- **TARGET — missense primary = REVEL** (Pejaver-2022 calibrated): PP3 supporting **≥ 0.644**, moderate **≥ 0.773**, strong **≥ 0.932**; BP4 supporting **≤ 0.290**, moderate **≤ 0.183**. **AlphaMissense** (likely_pathogenic **≥ 0.564**) as orthogonal support only. These are **reporting/tiering** values: as a *screen* filter they are provably inert (see [the unreachability note](#why-the-missense-predictors-never-filtered-anything)).
- **TARGET — splicing = Walker-2023 calibrated SpliceAI, on RAW scores**: PP3 Δ **≥ 0.2**, BP4 Δ **≤ 0.1**, 0.1–0.2 uninformative; canonical ±1,2 with Δ **≥ 0.5** = high tier. If PVS1(splice) applies, do not also apply PP3.
- **Report ONE predictor per variant — never stack correlated missense tools** as independent evidence; cap each tool at its Pejaver-calibrated strength.
- A gene-specific **ClinGen VCEP** score cutoff always overrides these generic defaults.
- **Scope limits:** SNV/indel only (no CNV/SV); pseudogene / segmental-duplication regions (e.g. PMS2/PMS2CL) are low-confidence from short reads and are **neither flagged nor masked today** (a TARGET). Full ledger: [limitations.md](limitations.md).

## VEP consequence ontology → IMPACT tiers

Ensembl VEP assigns Sequence Ontology consequence terms and bins them into four IMPACT tiers. The tier is a convenience filter, not a pathogenicity call.

| IMPACT | Representative SO consequences | Interpretation |
|---|---|---|
| **HIGH** | `transcript_ablation`, `splice_acceptor_variant`, `splice_donor_variant`, `stop_gained`, `frameshift_variant`, `stop_lost`, `start_lost`, `transcript_amplification` | Assumed protein-truncating / LoF / NMD-triggering; PVS1-eligible pending LOFTEE + decision tree |
| **MODERATE** | `missense_variant`, `inframe_insertion`, `inframe_deletion`, `protein_altering_variant` | Needs a calibrated missense/impact predictor |
| **LOW** | `splice_region_variant`, `synonymous_variant`, `stop_retained_variant`, `start_retained_variant`, `incomplete_terminal_codon_variant` | Usually low priority — but `splice_region_variant` can be highly deleterious |
| **MODIFIER** | non-coding / intronic / UTR / intergenic / regulatory | Predictions difficult; never sole evidence |

Caveat: `splice_region_variant` is binned LOW yet can abolish splicing. Never use IMPACT alone for any splice-adjacent variant — layer SpliceAI on top (see [Splicing predictors](#splicing-predictors--target-not-implemented), which is a TARGET: no SpliceAI layer exists today). Prioritize MANE Select / canonical transcripts, but retain per-transcript annotations so a consequence on a clinically relevant non-canonical transcript is not lost.

**As implemented**, IMPACT is not a convenience tier — it is the *first and dominant rung* of the screen. `keep_impacts: [HIGH, MODERATE]` keeps every pLoF and every missense outright, and `functional_reason()` **returns** there without consulting any score. LOW and MODIFIER variants get exactly one further chance: CADD. Two consequences worth stating plainly:

- A `splice_region_variant` (LOW) survives **only** if CADD carries it — the SpliceAI layer this section tells you to add does not exist yet.
- Step 2 runs VEP with `--flag_pick` (not `--pick`), so **every** consequence block is retained and merely marked `PICK=1`; the per-transcript annotations are preserved as advised, and `bcftools +split-vep -s` chooses the selection rule downstream.

## Predicted loss-of-function: LOFTEE and NMD-escape — TARGET, not implemented

> **Not running.** No LOFTEE data files are fetched, so there is no `LoF` / `LoF_filter` /
> `LoF_flags` annotation and no HC/LC label: a `stop_gained` in the last exon (likely
> NMD-escaping, plausibly benign) is **indistinguishable from a true null allele** here.
> The effect on *selection* is near zero — `keep_impacts` already keeps every HIGH-impact pLoF,
> and the retired `loftee_hc` branch could only fire on LoF calls VEP had *not* rated
> HIGH/MODERATE, close to an empty set. The cost is **tiering**: PVS1 grading needs it. The
> plugin code is already baked into the image at `/plugins`, so re-enabling is data acquisition
> (~13 GB, mostly the GERP bigwig) and config — not a rebuild. See
> [limitations.md §5](limitations.md).

LOFTEE is a VEP plugin that flags predicted LoF (stop-gained, frameshift, essential splice-site) as **HC** (high-confidence, passes all filters) or **LC** (low-confidence, fails ≥ 1 filter).

- **LC filters** include the **50-bp rule** (variant in the last exon or within 50 bp of the last exon–exon junction → likely NMD-escape), non-canonical splice sites of the affected intron, in-frame rescue by nearby splice sites, and small-intron GT/AG issues.
- **HC flags** (single-exon gene, weak PhyloCSF conservation, NAGNAG acceptor, non-canonical splice) warrant caution but are overridable by gene knowledge — roughly 14% of gnomAD HC pLoF carry a flag.
- Use `konradjk/loftee`'s **`grch38` branch** — the same repo, not a third-party fork; only its `master` branch is GRCh37-only. The Dockerfile clones that branch into `/plugins` (with `loftee_path:/plugins`), because the `ensembl-vep` base ships every other VEP plugin but builds with `--skip_plugins LoF`. Pin a commit SHA: the branch tip moves.
- Only the plugin **code** is baked into the image; the **data** (human_ancestor FASTA, GERP bigwig, `loftee.sql`) is host-fetched — which is the part that is not currently fetched.

For rare-disease and germline pediatric-cancer screening, treat **HC, no-flag pLoF** in a haploinsufficient / known disease gene as PVS1-eligible. Grade the PVS1 strength with the **Abou Tayoun 2018 / ClinGen SVI decision tree**: NMD-escape (last exon or 3′-terminal 50 bp) and single-exon context downgrade PVS1 from Very Strong. PVS1 application is gated on **ClinGen gene–disease validity ≥ Moderate** and a **known LoF disease mechanism**. See [gene_constraint.md](gene_constraint.md) for the haploinsufficiency weighting (pHaplo, ClinGen HI) that informs "is LoF the mechanism here."

## Missense / pathogenicity predictors — calibrated cutoffs — TARGET, not implemented

> **Not running.** No dbNSFP, REVEL or AlphaMissense file is present, so none of the scores in
> this section are annotated. Everything below is the calibration reference for the planned
> ACMG PP3/BP4 tiering step — **but read the next subsection before concluding the screen lost
> discriminative power, because it did not.**

### Why the missense predictors never filtered anything

This document previously credited REVEL/AlphaMissense/MPC with doing work they **provably never did**, and the correction matters more than the resource note:

- Every one of these is a **missense-only** score — by construction, a variant carrying a REVEL, AlphaMissense, MPC or MetaRNN value is a missense variant.
- **Every missense variant is `IMPACT=MODERATE`.**
- `MODERATE ∈ keep_impacts`, and `functional_reason()` **returns at the impact rung** before any predictor is consulted.

Therefore a variant could only reach a missense-predictor branch if it were simultaneously missense (to *have* a score) and not-MODERATE (to *get past* rung 1) — an empty set. **These branches were unreachable even back when the code contained them and dbNSFP was configured.** Their calibrated cutoffs did none of the discrimination this document advertised; removing dbNSFP cost the screen **exactly zero** selection power. `tests/integration/assert_integration.py` now asserts the `revel` / `alphamissense` / `mpc` / `spliceai` / `loftee_hc` keep-reasons never fire.

The genuine loss is **reporting and tiering**, not screening: a curator no longer sees a REVEL score beside a missense candidate, and PP3/BP4 will need one. If that step is built, ClinGen SVI requires committing to **one** predictor chosen *before* seeing results — REVEL is the ClinGen-calibrated option (AlphaMissense postdates the 2022 calibration). The re-add is ~1.3 GB of *dedicated* files, **not** dbNSFP (see [dbNSFP](#dbnsfp-as-the-aggregation-resource--dead-url-and-not-the-right-vehicle)).

### The calibration reference

The reference standard is **Pejaver 2022 (ClinGen SVI)**, which replaced the historic "count concordant tools" approach with single-tool, strength-tiered thresholds derived from local positive predictive value. **Report one predictor per variant** — correlated predictors do not constitute independent lines of evidence, and none may exceed its Pejaver-calibrated strength.

| Predictor | PP3 strong | PP3 moderate | PP3 supporting | BP4 supporting | BP4 moderate | BP4 strong |
|---|---|---|---|---|---|---|
| **REVEL** (primary) | ≥ 0.932 | ≥ 0.773 | ≥ 0.644 | ≤ 0.290 | ≤ 0.183 | ≤ 0.016 |
| **BayesDel** (noAF) | ≥ 0.50 | ≥ 0.27 | ≥ 0.13 | ≤ −0.18 | ≤ −0.36 | — |
| **VEST4** | ≥ 0.965 | ≥ 0.861 | ≥ 0.764 | < 0.764 | — | — |
| **CADD** (PHRED) | — | ≥ 28.1 | ≥ 25.3 | ≤ 22.7 | ≤ 17.3 | — |
| **MutPred2** | ≥ 0.932 | ≥ 0.829 | ≥ 0.737 | ≤ 0.391 | ≤ 0.197 | ≤ 0.010 |

REVEL is the only tool in this set reaching both PP3-strong and BP4-very-strong, and is this pipeline's **primary** missense predictor.

**AlphaMissense** (Cheng 2023) is used as **orthogonal** support, not stacked with REVEL. Its three bins are calibrated to ~90% ClinVar precision:

| AlphaMissense bin | Score |
|---|---|
| likely_benign | ≤ 0.34 |
| ambiguous | 0.34 – 0.564 |
| likely_pathogenic | ≥ 0.564 |

Independent ClinGen-style calibration supports AlphaMissense up to roughly PP3-moderate / BP4-strong.

**CADD v1.7** (Jan 2024; adds protein language models + regulatory CNNs; PHRED-scaled; GRCh38 chr1–22,X,Y; splice-aware since CADD-Splice / v1.6) is a genome-wide predictor covering non-missense sites, but is a weaker missense discriminator than the meta-predictors above — as tiering evidence, use it to fill gaps, not to override REVEL/AlphaMissense. It is the one predictor this pipeline actually has; see the next section for how that changes what it means.

Other predictors in dbNSFP — **PrimateAI-3D, EVE, MPC, MetaRNN, ClinPred** — would be available as orthogonal support (**TARGET**: none is annotated today). **MPC ≥ 2** strongly up-weights a missense variant in a constrained region (missense Z > 3.09 gives gene-level support; see [gene_constraint.md](gene_constraint.md)) — note MPC is missense-only and so falls under the unreachability argument above: it can inform *tiering*, never the screen. Do not sum multiple correlated missense predictors as independent evidence.

## CADD as the only predictor — IMPLEMENTED, and off-label

CADD v1.7 (via the dedicated plugin: `whole_genome_SNVs.tsv.gz` + the gnomAD indel table, so SNVs **and** indels are scored genome-wide) is the pipeline's **only** functional predictor, and the **only** keep-path for anything VEP rates below MODERATE. `cadd_phred_supporting: 25.3` is therefore not one threshold among many — **it is the entire non-coding screen**. Two problems, stated plainly rather than buried:

- **The number is named after a calibration that never applies to it.** 25.3 is Pejaver-2022's PP3-*supporting* cutoff for CADD, derived on **missense variants only**. But missense is MODERATE and is kept a rung earlier, so no missense variant ever reaches this comparison. In practice 25.3 is applied **exclusively to the non-coding variants it was not calibrated for**. Treat a `cadd` keep-reason as a **discovery rank** (≈ top 0.3% of genome-wide CADD scores), **not** as ACMG PP3 evidence, and never quote it as PP3 in a report.
- **There is no ClinGen-endorsed non-coding CADD threshold** to replace it with. This is not a lookup we skipped; the calibration does not exist. Lowering the cutoff widens discovery at a steep review cost (PHRED ≈ 20 → top 1%; ≈ 15 → top 3%), and the right value is region-dependent (promoter vs deep intron vs UTR).

Fixing this is **calibration work, not a download** — region-stratified thresholds or a purpose-built non-coding predictor. Until then, the honest reading of a non-coding candidate is "ranked high by a genome-wide scalar", not "predicted pathogenic". If CADD is not configured, Step 2 warns and the screen degrades to an impact-only filter with **no** evidence for any sub-MODERATE variant. See [limitations.md §4](limitations.md).

## Splicing predictors — TARGET, not implemented

> **Not running, and this is the single largest loss in the screen.** No SpliceAI file is present,
> so splice detection is whatever VEP's **positional** consequence terms reach:
> `splice_donor_variant` / `splice_acceptor_variant` (the ±1,2 dinucleotides, HIGH) and
> `splice_region_variant` (≈ ±3–8 intronic / 1–3 exonic, LOW). A variant creating a **cryptic
> splice site deep in an intron** (`intron_variant`, MODIFIER) or an **exonic synonymous change
> that disrupts splicing** (`synonymous_variant`, LOW) is **invisible to this screen** unless CADD
> happens to rescue it. CADD v1.6+ ingests SpliceAI and MMSplice as *input features*, so it is a
> **lossy proxy** — a strong splice signal tends to raise CADD — **not a substitute**: it is a
> single genome-wide scalar with no calibrated non-coding threshold. Cost to fix ≈ 0.6 GB; see
> [limitations.md §1](limitations.md) for the mirror, the SNV-only trap and the contig-naming trap.

- **SpliceAI** (delta score 0–1). The original heuristics were 0.2 (high recall), 0.5 (recommended), 0.8 (high precision). The **ClinGen SVI Splicing Subgroup (Walker 2023)** calibrated cutoffs for variants outside the essential ±1,2 dinucleotides:

  | SpliceAI Δ (**raw**) | Evidence | Note |
  |---|---|---|
  | ≥ 0.2 | **PP3** (supporting) | ~78% sensitivity |
  | 0.1 – 0.2 | uninformative | no evidence assigned |
  | ≤ 0.1 | **BP4** | ~87% specificity |

  Canonical ±1,2 splice variants with Δ ≥ 0.5 are treated as a high tier. **If PVS1(splice) already applies at any strength, do NOT also apply PP3** — that would double-count the splice evidence.

  **Raw vs masked — an earlier version of this document had this backwards.** Walker 2023 calibrated these thresholds on **RAW** scores: verbatim, *"We used the maximum raw SpliceAI delta score"*; the word "masked" does not appear in the paper. Any claim that 0.2 is masked-calibrated is wrong. This matters directionally: **masked ≤ raw always** (masking zeroes deltas that strengthen an existing annotated site), so applying a raw-calibrated 0.2 to masked scores is **more stringent than calibrated** and silently loses sensitivity.

  There is a **genuine tension between primary sources here, and it should be treated as such rather than resolved by fiat**: Illumina recommends *masked* scores for interpretation (they suppress deltas at sites already annotated as splice sites, which are usually not the interesting signal), while ClinGen calibrated on *raw*. Pick one deliberately, state which, and use its matching cutoff — do not mix a masked score with a raw-derived threshold.

  **The under-appreciated caveat is the window, not the mask.** Walker scored with a ±4,999 nt window, whereas SpliceAI's **precomputed** files use the default `-D 50`. Cryptic sites beyond ±50 nt of the variant simply are not in the precomputed table, so a precomputed lookup is not the thing Walker calibrated regardless of mask choice. A free **masked MANE mirror** does exist on the Ensembl FTP (SNV-only).
- **Pangolin** (Genome Biol 2023) shows comparable-or-superior sensitivity to SpliceAI in massively-parallel-assay benchmarks; good orthogonal confirmation at an analogous ~0.2 cutoff.
- **MaxEntScan** — classic motif model, useful for quantifying donor/acceptor strength change at canonical sites; supplementary to SpliceAI.

## dbNSFP as the aggregation resource — dead URL, and not the right vehicle

> **Not running, and not the recommended way back.** dbNSFP is no longer fetched, and its pinned
> download URL is **dead** — the S3 bucket returns `NoSuchBucket`; distribution moved to
> registration-gated downloads. More to the point, 30 GB delivered exactly five columns this
> pipeline ever read, none of which reached the screen (see
> [the unreachability note](#why-the-missense-predictors-never-filtered-anything)).
>
> **If these scores are wanted for tiering, use the dedicated files instead:**
> `AlphaMissense_hg38.tsv.gz` (643 MB) and `revel-v1.3_all_chromosomes.zip` (667 MB) — ~1.3 GB
> total versus 30 GB. Trap: Ensembl's `AlphaMissense.pm` emits `am_pathogenicity` / `am_class`,
> **not** `AlphaMissense_score`, so the INFO field name in `annotations.F` must match the plugin,
> not this document's old VEP command.

**dbNSFP** precomputes transcript-specific scores for nonsynonymous and splice-site SNVs, so a single annotation join delivers all missense predictors without per-tool installs: REVEL, AlphaMissense, CADD, BayesDel, VEST4, MetaRNN, ClinPred, MPC, PrimateAI, EVE, ESM1b, MutPred2, plus conservation (GERP, phyloP). Two structural limits are worth knowing even if it is never re-added: it is **coding/splice-site SNVs only** (no genome-wide non-coding coverage, no indels), and its CADD column is coding-only — which is why this pipeline takes CADD from the dedicated plugin instead.

- Current is **v5.3** (Oct 2025); **v4.9a** (Aug 2024) is the stable v4 academic branch. VEP 115's dbNSFP plugin cannot parse v5 filenames — pin **4.9a** if you use it at all.
- Use the **`a` (academic) branch** — it retains REVEL/CADD/VEST4/PolyPhen2/ClinPred, which the `c` (commercial) branch strips. AlphaMissense (CC-BY) was added at v4.7 and is present in both branches thereafter.
- **Pin the exact dbNSFP build** in the container manifest; a silent build change shifts every downstream score.

## Example annotation commands

Generic, parameterized invocations — substitute release, cache, and plugin data paths via config; never hard-code real paths.

**IMPLEMENTED** — what `pipeline/02_annotate_sites.sh` actually runs, once, on the cohort union (never per trio; Step 4 transfers annotations with `bcftools annotate`):

```bash
# VEP-only contract: release-matched cache + the CADD plugin. The cache supplies the
# gnomAD v4.1 per-population AFs (the rarity oracle) and ClinVar CLIN_SIG — --af_gnomade/
# --af_gnomadg and --check_existing are load-bearing, not extras.
vep \
  --cache --offline --dir_cache "${VEP_CACHE_DIR}" --cache_version "${VEP_RELEASE}" \
  --species homo_sapiens --assembly GRCh38 --fasta "${REF}" \
  --vcf --compress_output bgzip --force_overwrite --no_stats \
  --symbol --biotype --numbers --hgvs --canonical --mane \
  --af_gnomade --af_gnomadg --max_af --check_existing \
  --flag_pick --pick_order mane_select,mane_plus_clinical,canonical,rank \
  --plugin CADD,snv="${CADD_SNV}",indels="${CADD_INDEL}" \
  --fork "${THREADS}" -i "${IN_VCF}" -o "${OUT_VCF}"
```

`--flag_pick` (**not** `--pick`) keeps every consequence block and merely marks the chosen one `PICK=1`, so per-transcript annotation is retained as this document advises and the selection rule is applied downstream.

```bash
# Lift CSQ -> INFO with a vep_ prefix. NB: `bcftools +split-vep` is a PLUGIN and does NOT
# accept --threads ("unrecognized option") — passing it aborts the step.
bcftools +split-vep -c "${FIELDS}" -s "${SEL}" -p vep_ -Oz -o "${SPLIT_VCF}" "${OUT_VCF}"
```

**TARGET** — the full stack this document describes, for when the tiering step is built. None of these plugins is configured today; note `AlphaMissense.pm` emits `am_pathogenicity`/`am_class`, and prefer the dedicated REVEL/AlphaMissense files over dbNSFP:

```bash
vep ... \
  --plugin LoF,loftee_path:/plugins,human_ancestor_fa:"${ANCESTOR_FA}",gerp_bigwig:"${GERP_BW}" \
  --plugin AlphaMissense,file="${ALPHAMISSENSE_TSV}" \
  --plugin REVEL,"${REVEL_TSV}" \
  --plugin SpliceAI,snv="${SPLICEAI_SNV}",indel="${SPLICEAI_INDEL}"
```

## Recommended defaults (this pipeline)

For per-trio GATK genotype-refinement VCFs (GRCh38, not cohort-joint), rare-disease + germline pediatric cancer. All values are configurable defaults in `config/config.example.yaml`; a gene-specific ClinGen VCEP cutoff overrides any of them.

| Layer | Status | Default | Notes |
|---|---|---|---|
| Annotation stack | **IMPLEMENTED** | VEP 115 (GRCh38 cache, MANE prioritized) + **CADD v1.7 plugin only** | The whole contract. No LOFTEE/dbNSFP/SpliceAI/ClinVar file |
| Functional ladder | **IMPLEMENTED** | `keep_impacts: [HIGH, MODERATE]` → keep; **else** `cadd_phred_supporting: 25.3` → keep; else drop `not_functional` | Two rungs; nothing else. `src/hprv/selection.py` |
| CADD threshold | **IMPLEMENTED** | PHRED ≥ **25.3** | **Off-label**: missense-calibrated, applied only to non-coding. A discovery rank (≈ top 0.3%), *not* PP3 |
| ClinVar override | **IMPLEMENTED** | `CLIN_SIG` P/LP → keep; also rescues from `too_common`; `conflicting` excluded | **No star gate** — the cache has no `CLNREVSTAT`. Unstarred 1★ P/LP is honored. ClinVar pinned at 2025-02 |
| Hard benign | **IMPLEMENTED** | drop if AF ≥ **0.05** (BA1) | Never rescued, not even by ClinVar P/LP |
| Annotation stack (full) | TARGET | + LOFTEE (`grch38` branch) + SpliceAI + REVEL/AlphaMissense | Version-pin all; ~15 GB total. [limitations.md](limitations.md) |
| pLoF confidence | TARGET | LOFTEE **HC, no flags** | NMD-escape / last-exon / 3′-50 bp / single-exon → downgrade PVS1 (Abou Tayoun tree) |
| PVS1 gating | TARGET | ClinGen gene–disease validity ≥ Moderate + known LoF mechanism | See [gene_constraint.md](gene_constraint.md) |
| Missense (primary) | TARGET | **REVEL** PP3 supporting ≥ 0.644 / moderate ≥ 0.773 / strong ≥ 0.932; BP4 ≤ 0.290 / ≤ 0.183 | **Tiering only** — inert as a screen filter (missense is MODERATE, kept a rung earlier) |
| Missense (orthogonal) | TARGET | **AlphaMissense** likely_pathogenic ≥ 0.564; ambiguous 0.34–0.564; likely_benign ≤ 0.34 | Do not stack with REVEL as independent evidence |
| Regional missense | TARGET | **MPC ≥ 2** up-weights; missense Z > 3.09 gene support | Missense-only ⇒ tiering only |
| Splicing | TARGET | **SpliceAI** PP3 Δ ≥ 0.2; BP4 Δ ≤ 0.1; 0.1–0.2 uninformative; canonical ±1,2 Δ ≥ 0.5 high tier | Walker calibrated on **RAW** Δ, window ±4,999 nt. No PP3 if PVS1(splice) applies |
| Benign deprioritize (BP4) | TARGET | REVEL ≤ 0.183 **AND** AlphaMissense ≤ 0.34 **AND** SpliceAI Δ ≤ 0.1 | Deprioritize, not discard, if in a critical gene |
| Evidence hygiene | TARGET | One calibrated predictor per variant; cap at Pejaver strength | Never sum correlated predictors |

**Ordering:** apply the rarity gate first and genotype/QC gates (see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)) before functional tiering, so functional scoring runs only on variants that survive frequency and quality screening. Note the rarity gate is **not** `faf95` as implemented: the VEP cache carries no AC/AN, so faf95 cannot be computed at any price and the gate uses a grpmax **point-estimate proxy** over the grpmax-eligible groups (AFR/AMR/EAS/NFE/SAS) — see [allele_frequency.md](allele_frequency.md). Calibrated functional evidence would then feed ACMG/AMP classification (see [clinical_classification.md](clinical_classification.md)); that step is not built, so today's output is a screen, not a classification.

## Scope limitations (state honestly)

The complete ledger — every gap, why it exists, and what each costs to close — is
**[limitations.md](limitations.md)**. Read it before interpreting a negative result. Specific to
this layer:

- **The predictor stack is one predictor.** CADD, on a threshold calibrated for a variant class that never reaches it. There is no splice prediction, no pLoF confidence, no missense score, and no star-gated clinical evidence. A negative result from this layer means "no HIGH/MODERATE-impact variant and no CADD-high non-coding variant", **not** "nothing functional here".
- **SNV/indel only.** This layer does not detect CNV/SV, which account for ~10–15% of pediatric-cancer and rare-disease diagnoses (single-exon RB1/SMARCB1/DICER1/NF1 deletions, PMS2 rearrangements). LOFTEE and the missense predictors cannot see these either; a future GATK-gCNV / Manta / ExomeDepth module is required.
- **Pseudogene / segmental-duplication regions** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) yield low-confidence short-read calls and functional annotation on paralog-mapping variants is unreliable. **These regions are currently neither flagged nor masked** — that is a TARGET, not present behavior.
- **Non-coding regulatory variants** fall in VEP MODIFIER, where prediction is weak; they are not scored to reportable tiers by this layer alone. They are nonetheless the *only* class the CADD rung acts on, which is precisely the off-label problem above.
- **Predictor calibration** is anchored to ClinVar-derived truth sets (Pejaver, Walker, AlphaMissense); performance on genes/regions under-represented in those sets is not guaranteed. Pipeline-wide sensitivity/precision **are unmeasured, not measured-and-acceptable** — GIAB/CMRG truth sets and a positive-control panel remain TODO.

Every one of these is **additive to fix**: the VEP-only contract is a single seam — one `bcftools annotate`/plugin addition in Step 2 plus the INFO field in `annotations.F`. Nothing in the architecture forecloses any of it.

## Sources

- Pejaver et al. 2022, ClinGen SVI PP3/BP4 calibration, *Am J Hum Genet* — https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/ (DOI 10.1016/j.ajhg.2022.10.013); ClinGen summary — https://clinicalgenome.org/docs/calibration-of-computational-tools-for-missense-variant-pathogenicity-classification-and-clingen-recommendations-for-pp3-bp4-cri/
- Walker et al. 2023, ClinGen SVI Splicing Subgroup, *Am J Hum Genet* — https://pmc.ncbi.nlm.nih.gov/articles/PMC10357475/
- Cheng et al. 2023, AlphaMissense, *Science* — https://www.science.org/doi/10.1126/science.adg7492; thresholds — https://www.ebi.ac.uk/training/online/courses/alphafold/classifying-the-effects-of-missense-variants-using-alphamissense/understanding-pathogenicity-scores-from-alphamissense/
- Abou Tayoun et al. 2018, ClinGen SVI PVS1 decision tree, *Hum Mutat* — https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/
- Ensembl VEP calculated consequences & IMPACT — https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
- LOFTEE README — https://github.com/konradjk/loftee/blob/master/README.md (GRCh38 use requires the repo's `grch38` branch, not `master`); gnomAD flags context — https://gnomad.broadinstitute.org
- SpliceAI (Illumina) — https://github.com/Illumina/SpliceAI (masked-vs-raw and the `-D` distance default); lookup — https://spliceailookup.broadinstitute.org/
- Pangolin / splicing benchmark, *Genome Biol* 2023 — https://link.springer.com/article/10.1186/s13059-023-03144-z
- CADD v1.7 release notes (Jan 2024) — https://cadd.gs.washington.edu/static/ReleaseNotes_CADD_v1.7.pdf; CADD-Splice, *Genome Med* 2021 — https://link.springer.com/article/10.1186/s13073-021-00835-9
- dbNSFP releases/changelog — https://www.dbnsfp.org/releases ; https://sites.google.com/site/jpopgen/dbNSFP/changelog
