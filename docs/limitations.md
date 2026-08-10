# Limitations of the first pass

What this pipeline **cannot currently see**, why, and what each would cost to fix.

> Part of the high_priority_rare_variant methods reference. This is the honest counterpart to
> [Canonical defaults](README.md#canonical-defaults): that table says what the code does, this
> one says what it *doesn't*. If you are deciding whether a negative result means "not there"
> or "not looked for", start here.

## Why the first pass looks like this

The pipeline runs on a **VEP-centric contract**: VEP 115 GRCh38 — its cache plus its plugins (CADD,
and SpliceAI, required by default). No gnomAD, ClinVar, dbNSFP or LOFTEE file is bcftools-transferred in.

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
- **The PRECOMPUTED set is not exhaustive, and the screen currently accepts that.** Illumina's
  tables cover all theoretical SNVs but only **1 nt insertions and deletions of up to 4 nt**, at a
  narrow ±50 nt window, computed against GENCODE V24lift37 transcripts; a missing score is **not**
  "no splice effect" (see *Using SpliceAI to triage splice-altering variants in 7,220 individuals
  with rare conditions highlights limitations of the precomputed scores*, medRxiv 2025,
  doi:10.1101/2025.08.27.25334471 — re-running with updated annotations and `-D 500` recovered
  18.2% more predicted splice-altering variants and an 11.7% diagnostic increase, and larger indels
  were ~4x enriched for splice effects: 4.7% vs 1.1%). **Step 2b**
  (`resources.vep.spliceai_backfill.enabled`, **OFF by default**; if enabled and the image's
  isolated `spliceai` env is absent the run HALTS at preflight) closes most of this gap: it runs the stock Illumina model live over just the
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

### 2. faf95 — RESOLVED by an opt-in resource; the proxy remains the fallback

`frequency()` now prefers gnomAD's published **faf95** (the lower bound of the 95% Poisson CI),
transferred in Step 2 from the gnomAD v4.1 **joint** slim (`resources.gnomad.sites_slim`,
`prepare_resources.sh --only gnomad_sites fetch`, ~10 GB). It falls back **per variant** to the
grpmax **point-estimate proxy** — the max AF across AFR/AMR/EAS/NFE/SAS — wherever gnomAD
published no faf95. `rarity_oracle` reports which one fired on every row. The VEP cache itself
still carries no AC/AN, so without the slim the proxy is all there is; that path is unchanged.

Three properties worth stating, each verified against the real v4.1 data rather than assumed:

- **The FAF group set is afr/amr/eas/mid/nfe/sas** — `GRPMAX_POPS` plus `mid`, and it EXCLUDES the
  bottlenecked ami/asj/fin. So faf95 does not reintroduce the `MAX_AF` trap of §2a. `mid` is the
  one deviation, so `faf95_group` reports the producing group on every row.
- **Absent faf95 ≠ AF 0.** gnomAD emits `fafmax` only where a group's CI lower bound exceeds zero
  (74% of a chr22 sample carried none). The proxy fallback there is the *more stringent* of the
  two, so it can only ever filter more — never silently retain what faf95 would have caught.
- **Supplying the slim RETAINS MORE.** faf95 ≤ the point estimate, so the same cutoffs stop
  discarding low-count alleles the interval never justified discarding. A smaller candidate list
  after enabling it means a broken join, not a better filter — check Step 2's match count.

**Consequence of the fallback path** (no slim configured): since a point estimate is always ≥ its
own CI lower bound, every rarity gate fires slightly *more* often than a faf95 gate would. That
path **errs toward dropping** on low-count alleles — a false-negative direction. The error shrinks
as the group's
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

### 3. nhomalt — AVAILABLE with the gnomAD slim, reported rather than scored

gnomAD's homozygote count is the classic tell for a false recessive call: an allele with many
homozygotes in a population reference is unlikely to cause severe recessive disease. It rides in
with the faf95 slim (§2), so a `hom_recessive` / `compound_het` / `x_linked_recessive` call whose
allele gnomAD already carries homozygotes for now raises **`nhomalt_recessive_conflict`**.

The flag charges **0 points by default**, deliberately: how many homozygotes should disqualify a
recessive candidate depends on the condition's penetrance, age of onset and prevalence, none of
which hprv knows. Shipping a number would be exactly the uncalibrated value this document
criticises elsewhere. The flag is reported and filterable in the review table, and
`prioritization.composite.weights.quality.nhomalt_conflict` charges it if you decide on one.
Without the slim `nhomalt` is absent — which is NOT 0, and never raises the flag.

Note this gate was **never** applied to the biallelic modes anyway — the retired
`filters.denovo.require_gnomad_absent_or_singleton` only touched de novo (secondary here), and
was implemented as `nhomalt > 1`: a homozygote-count test, never the allele-count test its name
promised. So the practical loss is smaller than it looks, and the pre-existing gap is the more
interesting one.

**Cost to fix:** free, alongside §2 — `nhomalt_joint` is one of the fields the gnomAD slim keeps.

### 4. CADD is the general-purpose predictor, on an off-label threshold

Below MODERATE impact there are two keep-paths: SpliceAI (checked first, and only reaching
splice-disrupting variants) and CADD. So `cadd_phred_supporting: 25.3` is the **entire non-splice
non-coding screen** — everything intronic / synonymous / UTR / regulatory without a splice signal
rests on this one number. Two honest problems:

- **The number is named after a calibration that never applies to it.** 25.3 is Pejaver-2022's
  PP3-*supporting* cutoff, derived on **missense only**. Missense never reaches the CADD rung —
  it is MODERATE, kept at rung 1 — so in practice 25.3 is applied *exclusively* to the
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

### 6. ClinVar review status — RESOLVED, via the one bcftools transfer

The VEP cache exposes `CLIN_SIG` but **no `CLNREVSTAT`** at any price, so stars cannot come from
it. They now come from the ClinVar sites VCF itself, transferred in Step 2
(`resources.clinvar.vcf`) — the single non-CSQ annotation in the pipeline, and the reason the
transfer machinery exists at all. `clinvar_stars` (0-4) rides on every candidate row and into the
igv.js review table. Transferring the VCF also **un-stales** ClinVar: VEP 115 pins ClinVar
2025-02, where the VCF ships weekly, and the transferred release is version-pinned and recorded.

Three properties worth stating, because each is a place this could have gone wrong:

- **Stars RANK, they never gate.** The screen (Step 3) is deliberately star-blind: a 1★ assertion
  still reaches review, it is merely ranked below a 3★ one. Reinstating the old ≥2★ keep/drop gate
  would violate never-drop.
- **Blank is not zero.** Absent = the transfer did not run (nobody looked); `0` = ClinVar has a
  record whose submitter provided no assertion criteria. Conflating them would silently damp every
  P/LP assertion in a run with no ClinVar resource, so an absent star count leaves the Step-9
  clinical term at **full** weight.
- **Only the positive limb is damped.** A low-star *benign* assertion is not shrunk toward zero —
  that would promote a poorly-reviewed benign call, the opposite of the intent.

Reclassification is still real; treat P/LP as a triage prior, never an answer.

**Cost:** ~0.18 GB.

### 7. REVEL / AlphaMissense change the SCREEN by nothing — and that is structural

Listed because it looks alarming and isn't, and because the conclusion survives having wired them.

These are **missense-only** scores. Every missense is `IMPACT=MODERATE`. `selection.py` keeps
MODERATE at the impact rung and **returns before any predictor is consulted**. So these branches
are **unreachable whether or not the resource is configured** — a property of the ladder, not of
the annotation contract. CI asserts these keep-reasons never fire.

**Both are now wired** (VEP plugins over the dedicated files, not dbNSFP — see
[resources.md](resources.md)), and the sentence above is unchanged by that: they add and remove no
candidates. What they buy is the thing this section used to call the genuine loss — **reporting
and tiering**. A curator now sees REVEL/AlphaMissense next to a missense candidate, and Step 9's
missense tier is consulted in a fixed precedence (REVEL -> AlphaMissense -> off-label CADD ->
none) that always reports which predictor spoke, in `missense_evidence_source`. Fixed order, not a
max: ClinGen SVI says commit to **one** predictor chosen before seeing results, so best-of-N would
be an uncalibrated cherry-pick. REVEL leads because its cut points are the ones Pejaver 2022
calibrated; AlphaMissense follows (SVI now endorses it on par, and it postdates that calibration).

hprv still assigns **no ACMG weight** — the cut points ORDER candidates. MPC remains unwired.

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
