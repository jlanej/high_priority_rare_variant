# Annotation resources — acquisition & preparation

The container image ships the **software** (VEP 115, bcftools/tabix, LOFTEE plugin code, the conda
env). It deliberately does **not** ship the large reference **data** — those are prepared once on
your system and **bind-mounted** at runtime. This doc explains why, what you need, and how the
`scripts/prepare_resources.sh` helper fetches/verifies/wires it all.

> Every resource path in [`config/config.example.yaml`](../config/config.example.yaml) is already a
> `${ENV}` placeholder, and Step 2 fails loudly if a configured resource is missing. So the only
> gap between "image pulled" and "pipeline runs" is preparing this data and pointing the env vars at
> it — which `prepare_resources.sh --dir DIR emit-env` does for you.

> **STATUS — the VEP-only contract (read this before you download anything).**
> The pipeline's annotation source is **VEP 115 GRCh38 — its cache + its plugins: CADD, and the
> optional SpliceAI plugin** (see [## SpliceAI](#spliceai)). Step 2 performs **no external `bcftools
> annotate` transfers**. **gnomAD, ClinVar, dbNSFP and LOFTEE data are no longer fetched, bind-mounted
> or read** — the config keys that pointed at them are gone. The **required** acquisition therefore
> collapses to: **VEP cache (~24 GB) + CADD SNV+indel (~82 GB)** — plus, to turn on the splice
> keep-path, the **SpliceAI** raw score files — plus two small optional tables for Step-6 ranking.
> Everything else on this page is retained as the **shopping list for building on top** and is
> clearly marked *not currently used*. What the reduced set costs the screen — and what each
> re-addition buys — is the ledger in **[limitations.md](limitations.md)**; the declared source of
> truth for thresholds is [Canonical defaults](README.md#canonical-defaults). Neither is restated here.
>
> `prepare_resources.sh` follows the contract: bare `fetch` prepares only the required set
> (reference, VEP cache, CADD, constraint) — it will **not** start the ~877 GB gnomAD download —
> `verify` requires only that set and reports the rest as *not required*, and `emit-env` exports
> only the `${ENV}` vars the config still has keys for. The retired resources stay reachable
> behind an explicit `--only gnomad_sites,clinvar,…` so the [roadmap
> restorations](ROADMAP.md) are one flag away. One upstream fact to know regardless: the
> **pinned dbNSFP URL is dead** (see below).

## Why a prepare script, not a bundled image

1. **Size.** VEP cache ≈ 24 GB and CADD SNV ≈ 81 GB — even the *reduced* required set is ~107 GB.
   Baking that into a public image makes it un-pullable and violates the golden rule (the VEP cache
   is *never* baked; resources bind-mount).
2. **License.** CADD (not redistributable) legally **cannot** ship inside a public image. The VEP
   cache, reference FASTA, and the constraint tables are free to fetch and prepare automatically.
3. **Freshness / reproducibility.** The cache and CADD are versioned. You pin a version + checksum
   per resource (recorded in [`resources/manifest.env`](../resources/manifest.env)) so a run is
   reproducible and re-poolable, exactly as the reproducibility hardening prescribes.
   Note the corollary the contract inherits: **the cache freezes its own gnomAD and ClinVar
   snapshots** (VEP 115 pins ClinVar 2025-02), so "re-pin ClinVar" now means "move the cache".

So: **image = software; `prepare_resources.sh` = data.** Free resources auto-download, verify
(sha256), and index; license-gated resources are validated *if you provide them*, with precise
acquisition instructions printed otherwise. Nothing is installed — only downloaded.

## The resources you actually need

| Resource | Feeds | Config key (`${ENV}`) | Acquisition | ~Size |
|---|---|---|---|---|
| GRCh38 reference FASTA | VEP + `bcftools norm` | `reference.fasta` (`REF_FASTA`) | free | ~1 GB (gz) |
| VEP indexed cache (r115) | **everything**: consequence/IMPACT, gnomAD v4.1 AFs, ClinVar `CLIN_SIG` | `resources.vep.cache_dir` (`VEP_CACHE`) | free | ~24 GB |
| CADD SNV + indel | CADD plugin — the sole functional predictor, genome-wide, SNV+indel | `resources.vep.cadd_snv` / `cadd_indel` (`CADD_SNV`/`CADD_INDEL`) | **license-gated**, huge | ~82 GB |
| VEP plugin **code** (`.pm`) | CADD (LOFTEE code is baked but unused) | `resources.vep.plugins_dir` (`VEP_PLUGINS`) | **in the image** at `/plugins` (not fetched) | — |
| Constraint per-gene TSV | Step-6 ranking (LOEUF/pLI/s_het/pHaplo) — **optional**, skipped if unset | `resources.constraint.*` (`GNOMAD_V2_CONSTRAINT`) | free | small |
| Samocha mutation-rate table | Step-6 de-novo Poisson (secondary) — **optional**, skipped if unset | `resources.mutation_rate_table` (`MUTRATE_TABLE`) | free | small |

That is the whole required set. The constraint/mutation-rate tables are genuinely optional:
`run_pipeline.sh` passes them to Step 6 only when the path is set *and* exists, so Step 6 degrades
rather than fails without them.

Gene-list / phenotype resources (`ACMG_SF_LIST`, `PANELAPP_*`, `RECESSIVE_CPG_LIST`, `HPO_TERMS`)
are **reserved** (the overlay is not yet wired) and are not required to run the pipeline.

### What the cache supplies — and the one place it silently under-delivers

The VEP 115 GRCh38 cache carries gnomAD **v4.1** frequencies (since VEP r113) for **both** exomes
and genomes, plus a ClinVar snapshot, so `--af_gnomade --af_gnomadg --max_af --check_existing` is
the entire annotation acquisition story. Two properties of that source are load-bearing and are
**not** fixable by downloading something else:

- **Point AFs only.** The cache has **no `faf95`/`fafmax` field and no AC/AN**, so faf95's CI
  correction cannot be recomputed downstream at any price. `frequency()` is a *grpmax proxy* —
  the max AF over the grpmax-eligible groups (AFR/AMR/EAS/NFE/SAS) — not faf95. See
  [limitations.md §2](limitations.md).
- **Cache frequencies exist only for alleles accessioned into dbSNP.** An un-accessioned gnomAD
  variant silently returns *no* frequency and reads as "absent ⇒ rarest". Ensembl itself
  recommends `--custom` with the gnomAD VCF over `--af_gnomad*` for exactly this reason. The bias
  is toward **retention** (extra review), not toward missed calls — which is why the first pass
  accepts it.

**Do not "improve" the rarity field with `MAX_AF` or global `AF`.** Both are present in the CSQ,
both look better, and both are regressions in *opposite* directions. See
[limitations.md §2a](limitations.md); it is guarded by tests.

## Skipping acquisition entirely — bring your own VEP VCF

If you **already have a VEP 115 GRCh38 VCF** (e.g. from the group's existing DNM annotate run),
point `resources.vep.annotated_vcf` (`${VEP_ANNOTATED_VCF}`) at it, or pass `--vep-vcf` to Step 2.
Step 2 then **skips the VEP call entirely** and ingests your file as `cohort.sites.annotated.vcf.gz`
— **no cache, no CADD, no downloads at all**. It is not a blind trust: Step 2 **verifies the
assembly/build** rather than trusting the filename (a GRCh37 or older-release VCF would otherwise
annotate "successfully") and **asserts that frequencies actually landed**, failing loudly if they
did not — because with no gnomAD fields anywhere, rarity filtering would silently pass everything.

This is the cheapest path to a first run by a wide margin, and for many sites it is the *only*
step needed. Your VCF must carry the CSQ fields listed in `src/hprv/annotations.py` (`F`).

**Plugin code vs. plugin data.** VEP plugin *code* (the `.pm` scripts) is **software → baked into the
image**, not fetched. The `ensemblorg/ensembl-vep:release_115.0` base bundles the `Ensembl/VEP_plugins`
`.pm` at `/plugins` (`CADD.pm`, `dbNSFP.pm`, `SpliceAI.pm`, …) but builds with `--skip_plugins LoF`, and
**LOFTEE is a separate repo** — so our Dockerfile additionally bakes in the **`konradjk/loftee` grch38
branch** (master is GRCh37-only) at `/plugins` and installs its one missing Perl dep (`DBD::SQLite`;
`Bio::DB::BigFile`/Kent lib is already compiled into the base). `VEP_PLUGINS` therefore defaults to
`/plugins` via the image and needs no setup. Under the current contract **only `CADD.pm` is invoked**;
the rest of the plugin code sits inert in the image, which is what makes re-enabling any of the
optional resources below a *config* change rather than an image rebuild.

Exact URLs, versions, and checksums are pinned in [`resources/manifest.env`](../resources/manifest.env)
(re-pin there). Verified specifics for the **required** set:

- **VEP cache** — Ensembl **r115 indexed** cache `homo_sapiens_vep_115_GRCh38.tar.gz` (Ensembl FTP; use
  the plain Ensembl build, **not** refseq/merged). The cache version **must equal** the VEP binary (115).
  Reference FASTA = Ensembl **primary_assembly** (not toplevel). Open license. Note the cache
  precomputes **only SIFT and PolyPhen-2** as predictors — no REVEL/AlphaMissense/CADD/SpliceAI/MPC —
  which is why CADD is a separate download and why the others are absent from the screen.
- **CADD v1.7 — the complete CADD source.** `whole_genome_SNVs.tsv.gz` scores every possible SNV
  genome-wide (coding **and** non-coding) plus the precomputed indel set. Step 2 sources CADD only
  from the plugin (`vep_CADD_PHRED`). Fetched with `--accept-license`; point `CADD_SNV`/`CADD_INDEL`
  at your existing files to reuse them (**the group's existing DNM VEP call already runs this exact
  plugin**, so those files can very likely be reused as-is). Free-academic, non-redistributable.
  Read [limitations.md §4](limitations.md) before trusting `cadd_phred_supporting: 25.3`: because
  every missense is MODERATE and kept a rung earlier, that missense-calibrated cutoff is in practice
  applied *exclusively* to non-coding variants it was never calibrated for. It is a discovery rank,
  not ACMG PP3 evidence.
- **Constraint** — gnomAD v2.1.1 `lof_metrics.by_gene` (LOEUF `oe_lof_upper` + pLI) + Zeng-2024 s_het
  (Zenodo) + Collins-2022 pHaplo (Zenodo), left-joined by `scripts/join_constraint.py` into one
  per-gene TSV. All free (CC-BY/CC0). gnomAD v4 has no equivalent by-gene LOEUF flatfile — v2.1.1 is canonical.

> **Disk budget.** The required set is roughly: VEP cache ~24 GB (+ ~24 GB transient for the
> tarball) + FASTA ~1 GB + CADD ~82 GB + constraint <10 MB ≈ **110 GB**. Point
> `TMPDIR`/`APPTAINER_TMPDIR` at real disk with headroom. If you bring your own VEP VCF (above),
> the budget is **zero**.

## kraken2 database (optional, Step-8b NHF)

**Outside** the VEP-only annotation contract — this is a review-layer aid, not an annotation
source, and it is entirely optional. If you want the **non-human-fraction (NHF)** columns in the
igv.js review export (Step 8b: what fraction of a candidate's ALT-supporting reads classify as
non-human — a contamination / mis-mapping down-rank signal), point `resources.kraken2_db`
(`${KRAKEN2_DB}`) at a kraken2 database **directory**. Leave it unset and Step 8b simply warns and
skips, leaving the NHF columns blank (`outputs.igv.nonhuman_screen.enabled` defaults to `true`, but
"enabled" only means "screen *if* a DB is present").

| Resource | Role | Config key | Size / notes |
|---|---|---|---|
| kraken2 DB | taxonomic classifier index for ALT-read NHF | `resources.kraken2_db` (`KRAKEN2_DB`) | **bind-mounted DATA, never baked**; ~tens of GB |

- **Recommended DB: PrackenDB `k2_NCBI_reference`** (one NCBI reference assembly per species + human
  + RefSeq viral + UniVec-Core, built at kraken2's default k=35). The **nonhuman-screen package**
  (not this repo) ships a fetch helper `scripts/download_kraken2_db.sh --db DIR` that downloads,
  extracts, and validates it.
- **The dir MUST contain** the index files `hash.k2d`, `opts.k2d`, `taxo.k2d` **AND** the taxonomy
  dumps `taxonomy/nodes.dmp` + `taxonomy/names.dmp` (either under `taxonomy/` or at the DB root).
  Without the dumps kraken2 falls back to exact-taxid matching and the NHF signal is **unreliable in
  both directions** — Step 8b preflights for them and warns loudly, but there is no per-variant flag
  for it, so use a DB that ships them.
- **Software is in the image, DATA is not.** The `kraken2` binary (source-built at a pinned version)
  and `nonhuman-screen` (pip-pinned to a commit) are baked in; the DB is host-fetched and
  bind-mounted, exactly like the VEP cache and CADD. `.dockerignore` keeps resource DATA out of the
  build context.
- **Put the DB on local NVMe/tmpfs, not a FUSE mount.** kraken2's cost is the per-process DB load;
  Step 8b runs invocations **serially** with `--memory-mapping` (`nonhuman_screen.memory_mapping:
  true`) so the OS page cache stays warm across trios and the DB is paged in ~once per node instead
  of once per invocation. On a slow/network mount that win evaporates. At cohort scale, scatter like
  Step 2 (one kraken2 DB-load per array task).

## Can the VEP cache replace these? (Deliberately, yes — with a documented price)

An earlier version of this page argued "no, never skip the downloads." **The pipeline now does
exactly that, on purpose.** The trade was ~1.4 TB of acquisition — each piece with its own version
pin, license gate, index and contig-naming hazard — against a simple, sound, reproducible spine:
one annotation source, one frequency chokepoint, one functional ladder. The first pass
takes the spine. The price, per resource, is honest and bounded:

| Dropped | What the cache gives instead | The real price |
|---|---|---|
| gnomAD sites VCF | point AFs per population, v4.1, exomes + genomes | no **faf95** (unrecoverable — no AC/AN), no **nhomalt**. Rarity is a point estimate, erring toward *dropping* |
| ClinVar VCF | cache-frozen `CLIN_SIG` (2025-02) | no **CLNREVSTAT** ⇒ no star gate; stale. Over-*retains* |
| LOFTEE | VEP `IMPACT` | near-zero for *selection*; costs PVS1 tiering |

(SpliceAI was on this list as "the biggest loss" — it is now **wired** as the third functional
rung; see [## SpliceAI](#spliceai) and [limitations.md §1](limitations.md).)
| dbNSFP (REVEL/AM/MPC/MetaRNN) | SIFT/PolyPhen (unused) | **zero selection power** — see below |

**The dbNSFP row is the one that surprises people, so it is worth stating plainly:** REVEL,
AlphaMissense, MPC and MetaRNN are **missense-only** scores; every missense is `IMPACT=MODERATE`;
`selection.py` keeps MODERATE at the impact rung and **returns before any predictor is consulted**.
Those branches were therefore **unreachable even back when the code contained them and dbNSFP was
configured**. Removing a 30 GB download cost the screen **exactly zero** discrimination — CI now
asserts those keep-reasons never fire. The genuine loss is *reporting/tiering*, not detection.

Full ledger, with the cost-to-fix for each: **[limitations.md](limitations.md)**. Do not read the
table above as a to-do list — read §2a there first.

## Usage

Run the helper **inside the image** so bcftools/tabix/vep are on PATH (no host installs needed).
The script ships in the image at `/opt/hprv/scripts/` and is on `PATH`, so call it by name — you do
not need a checkout of this repo on the host:

> **A bare `fetch` already prepares only the required set** (`reference`, `vep_cache`, `cadd`,
> `constraint`) — it will **not** start the ~877 GB gnomAD download. The retired resources
> (`gnomad_sites`, `clinvar`, `loftee`, `dbnsfp`, `spliceai`) run **only** when you name them with
> `--only`, so the roadmap restorations are one flag away. Passing `--only reference,vep_cache,cadd,constraint`
> (as below) is therefore explicit-but-equivalent to a bare `fetch`.

```bash
# 1. fetch + prepare what the VEP-only contract reads (bare `fetch` does the same set;
#    the explicit --only just documents it)
apptainer exec --bind /data hprv.sif \
    prepare_resources.sh --dir /data/hprv_resources fetch \
    --only reference,vep_cache,cadd,constraint --accept-license

# 2. emit the export lines your config's ${ENV} placeholders expect
apptainer exec --bind /data hprv.sif \
    prepare_resources.sh --dir /data/hprv_resources emit-env --out /data/hprv_resources/resources.env
source /data/hprv_resources/resources.env      # then run_pipeline.sh --config ...
```

Valid `--only` ids: `reference`, `vep_cache`, `gnomad_sites`, `clinvar`, `loftee`, `constraint`,
`dbnsfp`, `spliceai`, `cadd`.

The pinned manifest ships alongside the script at `/opt/hprv/resources/manifest.env`. To re-pin a
version without rebuilding the image, bind-mount an edited copy and point `HPRV_RESOURCE_MANIFEST`
at it.

`fetch` is idempotent (skips a target that already exists and matches its checksum) and resumable.
`--accept-license` acknowledges the non-commercial terms of the gated resources (CADD, here).

### What each mode covers

- **`fetch`** (default) prepares `reference`, `vep_cache`, `cadd`, `constraint` — the required set,
  nothing more. Pass `--only gnomad_sites,clinvar,loftee,dbnsfp,spliceai` (any subset) to also
  prepare a retired resource for a [roadmap restoration](ROADMAP.md); they are never fetched by
  default, because an ~877 GB gnomAD download for data no step reads is not a sane default.
- **`verify`** requires exactly the set `fetch` prepares, and reports the retired resources as
  *not required* rather than failing on them. `run_pipeline.sh` additionally preflights what it
  actually needs before doing work.
- **`emit-env`** exports only the `${ENV}` placeholders `config/config.example.yaml` still has
  keys for. It also prints a commented `VEP_ANNOTATED_VCF` line — uncomment it to skip the VEP
  call entirely (see above).

## License-gated resources — what you must provide

**CADD** is free for academic/non-commercial use but **cannot be redistributed**; the script fetches
it from the UW host under `--accept-license`, or you can point `CADD_SNV`/`CADD_INDEL` at copies your
institution already has and skip it (`--only reference,vep_cache,constraint`).

**dbNSFP's pinned URL is dead** — `https://dbnsfp.s3.amazonaws.com/dbNSFP4.9a.zip` returns
`NoSuchBucket`; distribution moved to registration-gated downloads. The manifest entry and
`prep_dbnsfp` are both retained but **cannot succeed**. This is moot under the current contract
(dbNSFP is not used, and per the section above it never contributed selection power), and if you
ever do want those scores, the replacement is the **dedicated files, not dbNSFP** — see below.

## Optional resources — NOT currently fetched or used

**None of the following is downloaded, bind-mounted or read by the pipeline today.** This is the
shopping list for building on top: each is re-enabled by **one `bcftools annotate` transfer in
`02_annotate_sites.sh` plus its INFO field in `annotations.F`** — the contract is a single seam, and
the plugin code is already in the image. Ordered by value-per-GB. Sizes and rationale come from
[limitations.md](limitations.md); the diagnostic cost of each absence lives there, not here.

| Resource | ~Size | Why you'd add it | Acquisition trap |
|---|---|---|---|
| **ClinVar** GRCh38 VCF | ~0.18 GB | restores `CLNREVSTAT` ⇒ the ≥2★ gate, `CLNSIGCONF`, and un-stales ClinVar (monthly vs. the cache's 2025-02) | none — free, public domain, dated monthly release from NCBI |
| **REVEL + AlphaMissense** (dedicated) | ~1.3 GB | missense scores for *reporting/tiering* and a future PP3/BP4 step — **not** selection power | use `AlphaMissense_hg38.tsv.gz` (643 MB) + `revel-v1.3_all_chromosomes.zip` (667 MB), **never** the 30 GB dbNSFP (dead URL, 5 useful columns). Ensembl's `AlphaMissense.pm` emits `am_pathogenicity`/`am_class`, **not** `AlphaMissense_score` |
| **gnomAD v4.1 joint** slim | ~10 GB | restores real **faf95** + **nhomalt** (the homozygote sanity check) | stream-slim the 24 chromosome VCFs to ~5 of their 664 INFO fields; nothing but the slim output lands on disk and GCS egress is free. Confirmed v4.1 joint tags: `AF_joint`, `AF_grpmax_joint`, `fafmax_faf95_max_joint`, `nhomalt_joint` |
| **LOFTEE** GRCh38 data | ~13 GB | HC/LC pLoF confidence ⇒ PVS1 strength grading | mostly the GERP bigwig. Use the plugin's **`grch38` branch** (already baked in; master is GRCh37-only) |

## SpliceAI

**WIRED** (not optional-unused): the SpliceAI VEP plugin runs in Step 2 over precomputed raw delta
scores, and `selection.py` keeps a variant whose max delta score ≥ `filters.functional.spliceai_ds_min`
(default 0.2). Config keys `resources.vep.spliceai_snv` / `spliceai_indel` (`${SPLICEAI_SNV}` /
`${SPLICEAI_INDEL}`). Optional + graceful: unset ⇒ the splice keep-path is inactive (Step 2 warns).
Fetch/verify with `prepare_resources.sh --only spliceai`.

| File | Config key (`${ENV}`) | Source | ~Size |
|---|---|---|---|
| `spliceai_scores.raw.snv.hg38.vcf.gz` (+ `.tbi`) | `resources.vep.spliceai_snv` (`SPLICEAI_SNV`) | **full**: Illumina BaseSpace (login); no-login **MANE-only** mirror on Ensembl FTP | ~27 GB |
| `spliceai_scores.raw.indel.hg38.vcf.gz` (+ `.tbi`) | `resources.vep.spliceai_indel` (`SPLICEAI_INDEL`) | Illumina BaseSpace (login) — **no** free mirror | ~1 GB |

- **Use the FULL, RAW genome-wide files for "don't miss anything".** The **raw** (not masked) set
  keeps scores everywhere, so deep-intronic and non-canonical sites are not zeroed out. The free
  no-login Ensembl mirror is **MANE-select SNV only** (no non-MANE transcripts, no indel file) — a
  genuine subset; the full genome-wide raw SNV **+** indel VCFs are on Illumina BaseSpace (login),
  a one-time manual download. Both must be **bgzipped + `tabix -p vcf`**-indexed.
- **Threshold.** The default `spliceai_ds_min: 0.2` is the ClinGen SVI PP3-supporting cutoff.
  Walker-2023 calibrated on **raw** (not masked) scores, so a raw-calibrated 0.2 is correct for the
  raw files here. Lower it toward 0.1/0.05 for higher deep-intronic recall (the raw score rides
  through to the outputs regardless, so a below-cutoff score is still visible to a reviewer).
- **Two traps.** (1) The precomputed files default to a `-D 50` window — distant cryptic-site
  effects are outside the precomputed set regardless of build (running SpliceAI live catches more).
  (2) Contig naming: the Ensembl mirror uses `1`/`X`, GMKF is `chr`-prefixed — a mismatched `tabix`
  query returns empty with **exit code 0**. Step 2's presence guard ("no `vep_SpliceAI_pred_DS_*`
  lifted") catches a silently-dead plugin, but confirm the score VCF's contigs match your reference.
