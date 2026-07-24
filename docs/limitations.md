# Limitations of the first pass

What this pipeline **cannot currently see**, why, and what each would cost to fix.

> Part of the high_priority_rare_variant methods reference. This is the honest counterpart to
> [Canonical defaults](README.md#canonical-defaults): that table says what the code does, this
> one says what it *doesn't*. If you are deciding whether a negative result means "not there"
> or "not looked for", start here.

## Why the first pass looks like this

The pipeline runs on a **VEP-centric contract**: VEP 115 GRCh38 — its cache plus its plugins (CADD,
and optionally SpliceAI). No gnomAD, ClinVar, dbNSFP or LOFTEE file is bcftools-transferred in.

That was a deliberate trade. The alternative was ~1.4 TB of resource acquisition (gnomAD joint
sites alone are 877 GB), each piece with its own version pinning, license gate, index, contig-naming
hazard and failure mode — before a single trio could be screened. The VEP cache already carries
gnomAD v4.1 frequencies and ClinVar, and the group already runs VEP. So the first pass buys a
**simple, sound, reproducible spine** — one annotation source, one frequency chokepoint, one
functional ladder — and pays for it in the coverage documented below.

Most items here are **additive to fix** — and §1 (SpliceAI) already has been. The contract is a
narrow seam: each is re-enabled by **either a VEP plugin (as SpliceAI and CADD are) or one `bcftools
annotate` transfer** in `02_annotate_sites.sh`, plus its INFO field in `annotations.F`. Nothing in
the architecture forecloses any of it.

## The ledger

Ordered by what they cost a real diagnosis, worst first.

### 1. SpliceAI — now WIRED (was the largest loss); residual caveats below

**Resolved.** SpliceAI is now a first-class keep-path. Step 2 runs the SpliceAI VEP plugin over
the precomputed raw delta scores (`resources.vep.spliceai_snv`/`spliceai_indel`), split-vep lifts
`vep_SpliceAI_pred_DS_*`, and `selection.py` keeps a variant whose max delta score ≥
`filters.functional.spliceai_ds_min` (default **0.2**, the ClinGen SVI PP3-supporting cutoff). So a
variant that **creates a cryptic splice site 200 bp into an intron**, or a **synonymous exonic
change that disrupts splicing** — annotated `intron_variant` (MODIFIER) or `synonymous_variant`
(LOW), both below `keep_impacts` — is now nominated on its splice signal instead of being dropped.
The raw score rides through to `candidates.calls.tsv` / the IGV export / the xlsx for reviewer
tiering. It is **keep-only** (a missing/None score never drops a variant, it only fails to rescue).

**Residual caveats (why this is not a clean "solved"):**
- **The PRECOMPUTED set is not exhaustive — but there is now an optional live backfill.** Illumina's
  tables score genome-wide SNVs and a large indel set, but not every indel/context, at a narrow
  window; a missing score is **not** "no splice effect" (see *Using SpliceAI to triage splice-altering
  variants in 7,220 individuals*, medRxiv 2025). **Step 2b** (`resources.vep.spliceai_backfill.enabled`,
  off by default) closes most of this gap: it runs the stock Illumina model live over just the
  cohort variants that carry **no** precomputed score (default: indels only — SNVs are complete), at
  a **wider `-D` window (500)** to reach deep-intronic cryptic sites the precomputed `-D 50` set
  misses, and folds the result into the same `vep_SpliceAI_pred_DS_*` fields BEFORE Step 3. The model
  + GENCODE annotation are bundled in the image's isolated `spliceai` env (no extra download); it is
  TensorFlow inference (~1 var/s/CPU) but the unscored set is small, so it stays cheap. See
  [resources.md#spliceai](resources.md#spliceai).
- **The FULL raw genome-wide set needs a one-time manual download.** It lives on Illumina BaseSpace
  (login-gated). The free no-login Ensembl mirror is **MANE-select SNV only** (no non-MANE
  transcripts, no indel file) — a real subset. Use the full raw for "don't miss anything". See
  [resources.md](resources.md#spliceai).
- **Contig-naming trap.** The Ensembl mirror uses `1`/`X`; GMKF is `chr`-prefixed. A mismatched
  `tabix` query returns empty with **exit code 0** — Step 2's presence guard (`no
  vep_SpliceAI_pred_DS_* lifted`) catches a silently-dead plugin, but confirm the score VCF's
  contigs match your reference.
- CADD v1.6+ also ingests SpliceAI as an input feature, so it remains a weak backstop for splice
  signal below the SpliceAI keep threshold.

### 2. No faf95 — the rarity gate is a point estimate

`frequency()` returns a **grpmax proxy**: the max AF across the grpmax-eligible ancestry groups
(AFR/AMR/EAS/NFE/SAS). gnomAD's published `faf95` is the *lower bound of the 95% CI* of that
frequency. Computing it requires AC/AN. **The VEP cache carries neither**, so faf95 is not
approximated here — it is unrecoverable at any price.

**Consequence:** since a point estimate is always ≥ its own CI lower bound, every rarity gate
fires slightly *more* often than a faf95 gate would. The pipeline therefore **errs toward
dropping** on low-count alleles — a false-negative direction. The error shrinks as the group's
AN grows, so excluding the small bottlenecked groups (which the proxy does by construction, §2a)
removes the large half of it; the residual is bounded by AC and is worst for singletons in the
smaller eligible groups.

**Cost to fix:** ~10 GB and one long prep job. Stream-slim the 24 gnomAD v4.1 joint chromosome
VCFs down to 5 of their 664 INFO fields, keeping `fafmax_faf95_max_joint` + `nhomalt_joint`;
nothing but the slim output lands on disk, and GCS egress is free. This also restores §3.

#### 2a. A trap that is NOT a limitation — do not "fix" it

VEP's `MAX_AF` is right there in the CSQ and looks like a better rarity field. **It is not, and
using it would be a regression.** It maximises over the bottlenecked founder groups gnomAD's own
grpmax *deliberately excludes* (`ami` AN≈900, `asj`, `fin`, `mid`) **and** over the tiny 1000
Genomes phase-3 populations. A single allele in `ami` reads as AF ≈ 1.1e-3 — ten-fold over
`dominant_max` — so MAX_AF silently kills real ultra-rare dominant candidates. The global
`gnomADe_AF` / `gnomADg_AF` fail the opposite way: they dilute an ancestry-enriched benign
polymorphism across the whole cohort and retain it. **The two wrong substitutions err in opposite
directions; there is no single safe fallback.** Guarded by
`tests/test_pure.py::test_frequency_excludes_bottlenecked_pops` and the `GENEFND` integration case.

### 3. No nhomalt — no homozygote sanity check

gnomAD's homozygote count is the classic tell for a false recessive call: an allele with many
homozygotes in a population reference is unlikely to cause severe recessive disease. It is
unavailable, so `hom_recessive` / `compound_het` / `x_linked_recessive` calls carry no such check.

Note this gate was **never** applied to the biallelic modes anyway — the retired
`filters.denovo.require_gnomad_absent_or_singleton` only touched de novo (secondary here), and
was implemented as `nhomalt > 1`: a homozygote-count test, never the allele-count test its name
promised. So the practical loss is smaller than it looks, and the pre-existing gap is the more
interesting one.

**Cost to fix:** free, alongside §2 — `nhomalt_joint` is one of the fields the gnomAD slim keeps.

### 4. CADD is the only functional predictor, on an off-label threshold

CADD is the sole keep-path for anything VEP rates below MODERATE, which makes
`cadd_phred_supporting: 25.3` the **entire non-coding screen**. Two honest problems:

- **The number is named after a calibration that never applies to it.** 25.3 is Pejaver-2022's
  PP3-*supporting* cutoff, derived on **missense only**. Missense never reaches the CADD rung —
  it is MODERATE, kept a rung earlier — so in practice 25.3 is applied *exclusively* to the
  non-coding variants it was not calibrated for. Treat it as a discovery rank (≈ top 0.3%
  genome-wide), **not** as ACMG PP3 evidence.
- **There is no ClinGen-endorsed non-coding CADD threshold** to replace it with. Lowering it
  widens discovery at a steep review cost (PHRED 20 ≈ top 1%; 15 ≈ top 3%).

**Cost to fix:** not a resource problem — a calibration problem. Region-stratified thresholds, or
a purpose-built non-coding predictor, are research work rather than a download.

### 5. No LOFTEE — no pLoF confidence

No HC/LC label, so a `stop_gained` in the last exon (likely NMD-escaping and benign) is
indistinguishable from a true null allele.

**Effect on *selection*: near zero** — `keep_impacts` already keeps every HIGH-impact pLoF, and
the old `loftee_hc` branch only fired on LoF calls VEP had *not* rated HIGH/MODERATE, which is
close to an empty set. The real cost is **tiering**: PVS1 strength grading (Abou-Tayoun) needs it.
The plugin code is still baked into the image, so re-enabling is config, not a rebuild.

**Cost to fix:** ~13 GB (mostly the GERP bigwig).

### 6. No ClinVar review status — the ≥2★ gate is retired

The cache exposes `CLIN_SIG` but **no `CLNREVSTAT`**, so star ratings do not exist. A 1★
single-submitter P/LP assertion is now indistinguishable from an expert-panel one, and the
override honors any unstarred P/LP. This **over-retains** (more to review) rather than
over-dropping — the safe direction for a screen, but it does admit known-noisy assertions.

ClinVar is also **as stale as the cache**: VEP 115 pins ClinVar 2025-02, where a ClinVar VCF
ships monthly. Reclassification is real; treat P/LP as a triage prior, never an answer.

**Cost to fix:** ~0.18 GB — the cheapest item on this list by a wide margin, and it also
un-stales ClinVar.

### 7. No REVEL / AlphaMissense / MPC — and this costs the screen nothing

Listed for completeness, because it looks alarming and isn't.

These are **missense-only** scores. Every missense is `IMPACT=MODERATE`. `selection.py` keeps
MODERATE at the impact rung and **returns before any predictor is consulted**. So these branches
were **unreachable even when the code contained them and dbNSFP was configured** — their
calibrated cutoffs did none of the discriminative work the docs advertised. Removing dbNSFP cost
exactly zero selection power, and CI now asserts these keep-reasons never fire.

The genuine loss is **reporting/tiering**: a curator no longer sees a REVEL score next to a
missense candidate, and the planned ACMG PP3/BP4 step will need one. If that step is built,
ClinGen SVI says commit to **one** predictor chosen before seeing results — REVEL is the
ClinGen-calibrated option (AlphaMissense postdates the 2022 calibration).

**Cost to fix:** ~1.3 GB via the *dedicated* files (`AlphaMissense_hg38.tsv.gz` 643 MB,
`revel-v1.3_all_chromosomes.zip` 667 MB) — **not** dbNSFP, whose 30 GB delivered 5 columns we
read and whose pinned URL is dead anyway (the S3 bucket returns `NoSuchBucket`; it moved to
registration-gated downloads). Trap: Ensembl's `AlphaMissense.pm` emits `am_pathogenicity` /
`am_class`, **not** `AlphaMissense_score`.

## Scale

**Step 2 (VEP) is distributed.** It was the tall pole — one un-resumable single-node `--fork` run
that a WGS union (~57M sites) could not finish inside a 24 h walltime. It now shards by contig
(one VEP run per contig, each with its own `.done`), so it resumes across walltime kills and runs
one-contig-per-node as a coherent SLURM job graph. See
[pipeline/slurm/](../pipeline/slurm/README.md) for the `prep → plan → scatter[array] → gather →
downstream` graph, and `resources.vep.shard_by_contig` for the single-job in-process version. The
output is byte-identical to a single pass (`tests/integration/assert_shard_equivalence.sh`).

**Steps 0/1/4 still run per-trio serially — now the tall pole on WGS.** Not a correctness issue
and not a blocker: the run completes and resumes (`.done` guards). But each loops over trios one at
a time, and every iteration is independent (one trio VCF in, one per-trio file out, no cross-trio
state until Step 1's `concat`), so the work is embarrassingly parallel and simply is not dispatched
yet. `runtime.threads` does not help: bcftools `--threads` only adds BGZF (de)compression workers,
while `norm`'s reference lookups and left-alignment (the actual cost) are single-threaded.

Measured shape at ~200 trios:

| Input | Per trio (Step 1) | Serial total (Steps 0+1+4) |
|---|---|---|
| **Exome** (~150k variants/trio) | ~20–40 s | **~1.5 h — nothing to fix** |
| **WGS** (~4.5M variants/trio) | ~2.5–5 min | **~12–26 h** |

So: on exome, ignore this. On WGS, the natural next increment is a trio-array for Steps 0/1/4,
mirroring the Step-2 scatter — the `.done` idempotency needed to do it safely already exists; only
the dispatch is missing. [docs/tooling_and_reproducibility.md](tooling_and_reproducibility.md)
names this exact trigger — "adopt a manager when you need per-sample parallelism across many trios".

Two findings worth recording so nobody re-derives them:

- **Do not shrink `HPRV_PLAUSIBLE_PAD`.** It looks like an obvious win (a 1000 bp pad covers a lot
  of genome) and it is backwards: Step 4's region-restrict cost scales with the *number of
  regions*, not the coverage, because each costs a BGZF block re-decompression. Measured: pad=50 →
  15,949 regions → 3.46 s; pad=1000 → 11,633 regions → **2.63 s**. Shrinking the pad makes Step 4
  *slower*.
- **Step 8's serial CRAM slicing is deliberate**, not an oversight — it protects a shared
  network/FUSE mount, and `outputs.igv.extract_jobs` drives `samtools -@` for real intra-slice
  parallelism.

## Structural gaps (independent of the VEP-only contract)

These predate the contract and are tracked in [ROADMAP.md](ROADMAP.md):

- **No CNV/SV calling** — the largest coverage gap overall. A deletion removing an exon is invisible.
- **No phenotype layer** — no HPO/Exomiser prior, so a gene is ranked without regard to whether
  it fits the patient.
- **No co-segregation (PP1/BS4)**, no UPD rescue, no ROH.
- **Pseudogene/seg-dup regions** (PMS2/PMS2CL, SMN1/2, CYP21A2) are neither flagged nor masked.
- **No real-data validation** — GIAB/CMRG truth sets and a positive-control panel are still TODO,
  so sensitivity/precision are currently unmeasured rather than measured-and-acceptable.

### Inheritance-model residuals (Step 5)

Known, bounded gaps left open after the Step-5 model review. Each is a *recall* or *visibility*
limit, not a wrong call. See [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md) §3.

- **Compound-het pairing keys on the single VEP-PICK'd gene.** Each variant is indexed under one
  gene (the `--flag_pick` block `split-vep -s` selected), so two hits that are both damaging in the
  *same* gene fail to pair when the picked block names a different overlapping gene for one of them.
  This bites only where a variant is damaging at **equal** consequence rank in two overlapping MANE
  genes (shared exon, readthrough/bicistronic loci). Because hets are gathered at the permissive
  1e-2 gate but the dominant fall-through needs 1e-4, an unpaired second hit between those bounds is
  emitted under **no mode at all**. Fix: lift `Gene` with `-s all` for the pairing key only, and
  index each variant under every gene it hits. (slivar's `comphet` evaluates all transcripts for
  exactly this reason.)
- **A male non-PAR chrX hemizygote is dropped by every mode when the mother's genotype is
  uninformative.** For a male child on chrX the only reachable branches split the maternal genotype
  into `{HOM_REF}` (de novo) and `{HET, HOM_ALT}` (X-linked recessive); a maternal **no-call** falls
  through both with no row and **no audit counter**. A hemizygous LoF in an affected boy is causally
  self-sufficient — the maternal genotype separates inherited from de novo/germline-mosaic, it does
  not establish causality — so the exclusion is most costly exactly here. Rate-limited (the mother is
  diploid on X at ~2× the son's coverage). Fix: emit with `flags=maternal_gt_uninformative`, or at
  minimum add a dropped-count audit metric so a negative is distinguishable from "not looked for".
- **A permissive comp-het partner can still suppress a dominant call.** Hets are pooled for pairing
  at `recessive_max` (1e-2) but the dominant gate is `dominant_max` (1e-4). A *phase-confirmed*
  `mat × pat` pair consumes both legs, so a genuinely dominant-grade variant is re-labelled
  `compound_het` whenever the child happens to carry any other sub-1e-2 functional het in the same
  gene — near-certain in long genes (TTN, NEB, RYR1, DMD). The variant is **not** lost (Step 6 unions
  the per-model trio sets, so `n_carriers`/`recurrent` are unaffected), but it leaves `n_dominant`
  and therefore the headline `p_recurrence`. The *unconfirmed* (de-novo-partner) half of this was
  fixed — such pairs no longer consume. Fix for the rest: require the partner to be biallelic-credible
  before consuming, or emit the dominant row too with an `also_comphet_partner` flag.
- **No Y-linked inheritance model.** chrY is deliberately routed away from the mother-keyed
  hemizygous models, so male non-PAR chrY produces no rows at all. Clinically near-empty (Y-linked
  Mendelian SNV disease is essentially confined to spermatogenic failure, whose lesions are CNVs).

## Reading a negative result

Given the above, "no candidate found" for a trio means: no **coding** variant (or CADD-high
non-coding variant) passing a **point-estimate** rarity gate, in a **SNV/indel** callset, without
phenotype weighting, splice prediction, CNV calling, or star-gated clinical evidence. That is a
useful screen. It is not an exclusion.
