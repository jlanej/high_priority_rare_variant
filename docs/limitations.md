# Limitations of the first pass

What this pipeline **cannot currently see**, why, and what each would cost to fix.

> Part of the high_priority_rare_variant methods reference. This is the honest counterpart to
> [Canonical defaults](README.md#canonical-defaults): that table says what the code does, this
> one says what it *doesn't*. If you are deciding whether a negative result means "not there"
> or "not looked for", start here.

## Why the first pass looks like this

The pipeline runs on a **VEP-only contract**: a VEP 115 GRCh38 cache plus the CADD plugin, and
nothing else. No gnomAD, ClinVar, dbNSFP, SpliceAI or LOFTEE file is downloaded or bind-mounted.

That was a deliberate trade. The alternative was ~1.4 TB of resource acquisition (gnomAD joint
sites alone are 877 GB), each piece with its own version pinning, license gate, index, contig-naming
hazard and failure mode — before a single trio could be screened. The VEP cache already carries
gnomAD v4.1 frequencies and ClinVar, and the group already runs VEP. So the first pass buys a
**simple, sound, reproducible spine** — one annotation source, one frequency chokepoint, one
functional ladder — and pays for it in the coverage documented below.

Everything here is **additive to fix**. The contract is a single seam: each item below is
re-enabled by one `bcftools annotate` transfer in `02_annotate_sites.sh` plus its INFO field in
`annotations.F`. Nothing in the architecture forecloses any of it.

## The ledger

Ordered by what they cost a real diagnosis, worst first.

### 1. No SpliceAI — deep-intronic and synonymous splice variants are invisible

**The largest loss, and the only one that removes a whole variant class.**

Without SpliceAI, splice detection is whatever VEP's *positional* consequence terms can reach:
`splice_donor_variant` / `splice_acceptor_variant` (the ±1,2 dinucleotides, HIGH) and
`splice_region_variant` (roughly ±3–8 intronic / 1–3 exonic, LOW). That is it.

So a variant that **creates a cryptic splice site 200 bp into an intron**, or a **synonymous
exonic change that disrupts splicing**, is annotated `intron_variant` (MODIFIER) or
`synonymous_variant` (LOW) — both below `keep_impacts` — and is dropped at Step 3 unless CADD
happens to rescue it. This is a recognised rare-disease diagnostic category, and the screen
cannot nominate it.

**Partial mitigation, not a substitute:** CADD v1.6+ ingests SpliceAI and MMSplice as *input
features*, so a strong splice signal tends to raise CADD. It is a lossy re-encoding through a
single genome-wide scalar, with no calibrated non-coding threshold (see §4).

**Cost to fix:** ~0.6 GB. Filter the free Ensembl MANE mirror to Δ ≥ 0.1 (only ~1% of records
reach the 0.2 supporting cutoff, so the full 27 GB table is ~99% dead weight for this pipeline).
Because SpliceAI is a *keep-only* rule — absence yields `None`, which never drops anything — a
filtered subset is lossless for current behaviour. Two traps: the free mirror is **SNV-only**
(no indel file exists on it at all; indels are BaseSpace/login-gated), and it uses Ensembl-style
contigs (`1`, `X`) while GMKF is `chr`-prefixed — a mismatched `tabix` query returns empty with
**exit code 0**, so assert a non-zero record count.

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

## Scale: the per-trio steps run serially

Not a correctness issue and not a blocker — the run completes, and it resumes (`.done` guards) if
a walltime kill interrupts it. But it is the one place where a large **WGS** cohort costs an order
of magnitude of wall-clock it does not need to.

Steps 0, 1 and 4 each loop over trios one at a time, and `run_pipeline.sh` is a single serial
process. Every iteration is independent — one trio VCF in, one per-trio file out, no cross-trio
state until Step 1's `concat` — so the work is embarrassingly parallel and simply is not
dispatched. `runtime.threads` does not help: bcftools `--threads` only adds BGZF (de)compression
workers, while `norm`'s reference lookups and left-alignment (the actual cost) are single-threaded.

Measured shape at ~200 trios:

| Input | Per trio (Step 1) | Serial total (Steps 0+1+4) |
|---|---|---|
| **Exome** (~150k variants/trio) | ~20–40 s | **~1.5 h — nothing to fix** |
| **WGS** (~4.5M variants/trio) | ~2.5–5 min | **~12–26 h** (≈2–3 h at 10-way) |

So: on exome, ignore this. On WGS, either accept a long resumable run, or dispatch the per-trio
loops (`xargs -P`, or a scheduler array over manifest shards). The `.done` idempotency needed to
do it safely already exists; only the dispatch is missing.
[docs/tooling_and_reproducibility.md](tooling_and_reproducibility.md) names this exact trigger —
"adopt a manager when you need per-sample parallelism across many trios".

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

## Reading a negative result

Given the above, "no candidate found" for a trio means: no **coding** variant (or CADD-high
non-coding variant) passing a **point-estimate** rarity gate, in a **SNV/indel** callset, without
phenotype weighting, splice prediction, CNV calling, or star-gated clinical evidence. That is a
useful screen. It is not an exclusion.
