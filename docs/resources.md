# Annotation resources — acquisition & preparation

The container image ships the **software** (VEP 115, bcftools/tabix, LOFTEE plugin code, the conda
env). It deliberately does **not** ship the large reference **data** — those are prepared once on
your system and **bind-mounted** at runtime. This doc explains why, what you need, and how the
`scripts/prepare_resources.sh` helper fetches/verifies/wires it all.

> Every resource path in [`config/config.example.yaml`](../config/config.example.yaml) is already a
> `${ENV}` placeholder, and Step 2 fails loudly if a configured resource is missing. So the only
> gap between "image pulled" and "pipeline runs" is preparing this data and pointing the env vars at
> it — which `prepare_resources.sh --dir DIR emit-env` does for you.

## Why a prepare script, not a bundled image

1. **Size.** VEP cache ≈ 25 GB, gnomAD v4.1 joint sites are hundreds of GB (we slim them to tens),
   dbNSFP ≈ 30 GB, CADD SNV ≈ 80 GB, SpliceAI ≈ 40 GB. Baking that into a public image makes it
   un-pullable and violates the golden rule (the VEP cache is *never* baked; resources bind-mount).
2. **License.** dbNSFP (redistribution-restricted), SpliceAI/Illumina (login-gated, non-commercial),
   and CADD (not redistributable) legally **cannot** ship inside a public image. ClinVar, gnomAD,
   the LOFTEE data, and the constraint tables are free to fetch and prepare automatically.
3. **Freshness / reproducibility.** ClinVar is weekly, gnomAD/dbNSFP are versioned. You pin a
   version + checksum per resource (recorded in [`resources/manifest.env`](../resources/manifest.env))
   so a run is reproducible and re-poolable, exactly as the reproducibility hardening prescribes.

So: **image = software; `prepare_resources.sh` = data.** Free resources auto-download, verify
(sha256), and index; license-gated resources are validated *if you provide them*, with precise
acquisition instructions printed otherwise. Nothing is installed — only downloaded.

## The resources

| Resource | Feeds | Config key (`${ENV}`) | Acquisition | ~Size |
|---|---|---|---|---|
| GRCh38 reference FASTA | VEP + `bcftools norm` | `reference.fasta` (`REF_FASTA`) | free | ~3 GB |
| VEP indexed cache (r115) | VEP core | `resources.vep.cache_dir` (`VEP_CACHE`) | free | ~25 GB |
| VEP plugin **code** (`.pm`) | CADD/dbNSFP/SpliceAI + LOFTEE grch38 | `resources.vep.plugins_dir` (`VEP_PLUGINS`) | **in the image** at `/plugins` (not fetched) | — |
| gnomAD v4.1 **joint** sites | rarity oracle (**faf95**) | `resources.gnomad.sites_vcf` (`GNOMAD_SITES`) | free, large (slimmed) | ~tens GB |
| ClinVar GRCh38 VCF | clinical evidence | `resources.clinvar.vcf` (`CLINVAR_VCF`) | free | ~0.2 GB |
| dbNSFP | REVEL/AlphaMissense/MPC/MetaRNN (its CADD_phred is **not** requested — CADD comes from the plugin) | `resources.vep.dbnsfp` (`DBNSFP`) | **license-gated** | ~30 GB |
| CADD SNV + indel | CADD plugin — **the complete, genome-wide CADD** (dbNSFP CADD dropped) | `CADD_SNV` / `CADD_INDEL` | **license-gated**, huge | ~80 GB |
| SpliceAI precomputed | splicing (SpliceAI plugin) | `SPLICEAI_SNV` / `SPLICEAI_INDEL` | **license-gated** | ~40 GB |
| LOFTEE GRCh38 data | pLoF confidence (LoF plugin) | `resources.vep.loftee_data` (`LOFTEE_DATA`) | free | ~several GB |
| Constraint per-gene TSV | Step-6 ranking (LOEUF/pLI/s_het/pHaplo) | `resources.constraint.*` (`GNOMAD_V2_CONSTRAINT`) | free | small |
| Samocha mutation-rate table | Step-6 de-novo Poisson (secondary) | `resources.mutation_rate_table` (`MUTRATE_TABLE`) | free | small |

Gene-list / phenotype resources (`ACMG_SF_LIST`, `PANELAPP_*`, `RECESSIVE_CPG_LIST`, `HPO_TERMS`)
are **reserved** (the overlay is not yet wired) and are not required to run the pipeline.

**Plugin code vs. plugin data.** VEP plugin *code* (the `.pm` scripts) is **software → baked into the
image**, not fetched. The `ensemblorg/ensembl-vep:release_115.0` base bundles the `Ensembl/VEP_plugins`
`.pm` at `/plugins` (`CADD.pm`, `dbNSFP.pm`, `SpliceAI.pm`, …) but builds with `--skip_plugins LoF`, and
**LOFTEE is a separate repo** — so our Dockerfile additionally bakes in the **`konradjk/loftee` grch38
branch** (master is GRCh37-only) at `/plugins` and installs its one missing Perl dep (`DBD::SQLite`;
`Bio::DB::BigFile`/Kent lib is already compiled into the base). `VEP_PLUGINS` therefore defaults to
`/plugins` via the image and needs no setup. Only the plugin **data** (CADD tsv, dbNSFP tsv, SpliceAI
VCF, and LOFTEE's `human_ancestor.fa.gz` / GERP bigwig / `loftee.sql`) is host-fetched by this script.

Exact URLs, versions, and checksums are pinned in [`resources/manifest.env`](../resources/manifest.env)
(re-pin there). Verified specifics:

- **VEP cache** — Ensembl **r115 indexed** cache `homo_sapiens_vep_115_GRCh38.tar.gz` (Ensembl FTP; use
  the plain Ensembl build, **not** refseq/merged). The cache version **must equal** the VEP binary (115).
  Reference FASTA = Ensembl **primary_assembly** (not toplevel). Open license.
- **gnomAD v4.1 joint** — `gs://gcp-public-data--gnomad/release/4.1/vcf/joint/…` (also S3 / HTTPS).
  Raw is **~877 GB across 24 files**; the script downloads each chromosome, keeps **only** the 4 INFO
  tags below, deletes the raw, and concatenates → tens of GB persistent. The v4.1 joint tag spellings
  are **confirmed to match the config defaults exactly**: `AF_joint`, `AF_grpmax_joint`,
  `fafmax_faf95_max_joint`, `nhomalt_joint`. Public domain.
- **ClinVar** — dated monthly GRCh38 VCF from NCBI (pin the date in `manifest.env`). Public domain.
- **dbNSFP `4.9a`** — the **academic** build (S3, direct URL, `--accept-license`). Pin **4.9a**, not
  v5.x (the VEP-115 plugin can't parse v5 filenames) and the **`a`** build, not `c` (which omits
  REVEL/CADD/AlphaMissense). Columns confirmed present: `REVEL_score, AlphaMissense_score,
  AlphaMissense_pred, MPC_score, MetaRNN_score, CADD_phred`. Heavy prep (a genome-wide sort — point
  scratch at real disk with tens of GB free). Non-commercial + citation.
- **CADD v1.7 — the complete CADD source.** `whole_genome_SNVs.tsv.gz` scores every possible SNV
  genome-wide (coding **and** non-coding) plus the precomputed indel set — strictly more complete than
  dbNSFP's coding-only `CADD_phred`, which is therefore **dropped** from the dbNSFP column list. Step 2
  sources CADD only from the plugin (`vep_CADD_PHRED`). Fetched by default (with `--accept-license`);
  point `CADD_SNV`/`CADD_INDEL` at your existing files to reuse them. Free-academic, non-redistributable.
- **SpliceAI** — the full genome-wide set is Illumina **BaseSpace (login-gated)**; the script instead
  auto-fetches the **no-login Ensembl MANE mirror** for SNVs (`--accept-license`), and instructs you to
  supply the indel file from BaseSpace. The MANE mirror covers MANE-select transcripts only. Non-commercial.
- **LOFTEE** — GRCh38 data from the Broad host (`human_ancestor.fa.gz` + `.fai`/`.gzi`, `loftee.sql.gz`
  → gunzip, 12 GB GERP bigwig). Use the plugin **`grch38` branch** (master is GRCh37-only). Free.
- **Constraint** — gnomAD v2.1.1 `lof_metrics.by_gene` (LOEUF `oe_lof_upper` + pLI) + Zeng-2024 s_het
  (Zenodo) + Collins-2022 pHaplo (Zenodo), left-joined by `scripts/join_constraint.py` into one
  per-gene TSV. All free (CC-BY/CC0). gnomAD v4 has no equivalent by-gene LOEUF flatfile — v2.1.1 is canonical.

> **Disk budget.** Plan for roughly: VEP cache ~25 GB + FASTA ~1 GB + gnomAD slimmed ~tens of GB
> (transient ~900 GB streamed if you don't `--only` a subset) + ClinVar ~0.2 GB + dbNSFP ~30 GB
> (+ ~50 GB transient sort) + LOFTEE ~13 GB + constraint <10 MB. Point `TMPDIR`/`APPTAINER_TMPDIR`
> at real disk with headroom.

## Can the VEP cache replace any of these? (No — a documented trap)

A `--everything` VEP run emits gnomAD **AF** (`gnomADe_AF`/`gnomADg_AF`/`MAX_AF`), a coarse
cache-frozen `CLIN_SIG`, and `SIFT`/`PolyPhen` straight from the cache — so it is tempting to skip
the downloads below. **Do not.** None of them replace a load-bearing resource for this pipeline:

- **gnomAD** — the cache gives point **AF** (exome & genome *separately*), not **`faf95`** (the sole
  rarity oracle, golden rule #2), not the **v4.1 *joint*** combine, and not **`nhomalt`** (needed for
  the recessive / de-novo `nhomalt ≤ 1` checks). `frequency()` intentionally has **no AF fallback**;
  wiring cache AF in would silently re-break the faf95 contract. → keep the sites VCF.
- **ClinVar** — the cache's `CLIN_SIG` has **no `CLNREVSTAT`** (star rating → the ≥2★ auto-promote
  gate can't run), no `CLNSIGCONF`, and is frozen at the cache's ClinVar version. → keep the dated VCF.
- **dbNSFP / SpliceAI** — the cache only carries the legacy `SIFT`/`PolyPhen`, not
  **REVEL / AlphaMissense / MPC / MetaRNN** or **SpliceAI** (the Step-3 functional/splice gate). → keep both.
- **LOFTEE** — VEP's own `IMPACT` is not LOFTEE's **`LoF`/`LoF_flags`** HC/LC confidence. → keep it.
- **CADD** — the one overlap. We keep the **dedicated plugin** (genome-wide `CADD_PHRED`, the most
  complete) and **drop dbNSFP's coding-only `CADD_phred`** so there is a single authoritative CADD
  field. dbNSFP is still required for REVEL/AlphaMissense/MPC. (Your existing DNM VEP call already runs
  this exact CADD plugin, so those files can be reused.)

Where the cache *does* help: `gnomADe_AF`/`MAX_AF` are a good **cheap pre-filter + QC cross-check**
against faf95, and since AF is free from the cache, the gnomAD slim strictly only needs
`fafmax_faf95_max_joint` + `nhomalt_joint` (we also keep `AF_joint`/`AF_grpmax_joint` for reporting).

## Usage

Run the helper **inside the image** so bcftools/tabix/vep are on PATH (no host installs needed):

```bash
# 1. fetch + prepare everything free; validate any gated files you've placed under --dir
apptainer exec --bind /data hprv.sif \
    scripts/prepare_resources.sh --dir /data/hprv_resources fetch

# 2. check every expected file is present + indexed (gnomAD INFO tags verified against config)
apptainer exec --bind /data hprv.sif \
    scripts/prepare_resources.sh --dir /data/hprv_resources verify

# 3. emit the export lines your config's ${ENV} placeholders expect
apptainer exec --bind /data hprv.sif \
    scripts/prepare_resources.sh --dir /data/hprv_resources emit-env --out /data/hprv_resources/resources.env
source /data/hprv_resources/resources.env      # then run_pipeline.sh --config ...
```

`fetch` is idempotent (skips a target that already exists and matches its checksum) and resumable.
Use `--only gnomad_sites,clinvar` to prepare a subset, and `--accept-license` to acknowledge the
non-commercial terms of the gated resources you supply.

## License-gated resources — what you must provide

dbNSFP, CADD, and SpliceAI are free for academic/non-commercial use but **cannot be redistributed**,
so the script will not download them for you. Obtain them through your institution and drop them
under `--dir` (the script prints the exact expected filenames); `fetch`/`verify` then validate and
index them. See the per-resource instructions the script prints, and the pinned sources in
`resources/manifest.env`.
