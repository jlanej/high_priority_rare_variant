# Roadmap — high-priority gaps vs state of the art

From a 10-domain SOTA literature review (2023–2026) against the current repo. Items are ranked by
**impact × defensibility-per-unit-effort** for the pipeline's headline claim: *recurrent +
constrained inherited variants = candidates*. Effort is "to add here." This is a living document —
update it as items land. (Out of scope by design: de novo review, mtDNA — separate pipelines.)

## Restoring the VEP-only contract's gaps (start here — best ROI in the repo)

The first pass runs on a VEP 115 cache + CADD alone; see **[limitations.md](limitations.md)** for
the full ledger of what that costs and why each was an acceptable trade. Each row below is
**additive**: one `bcftools annotate` transfer in `02_annotate_sites.sh` plus its INFO field in
`annotations.F`. Nothing in the architecture blocks any of them. Ranked by clinical value per GB.

| # | Restore | Size | Buys | Notes |
|---|---------|:----:|------|-------|
| R1 | **ClinVar VCF** | ~0.18 GB | `CLNREVSTAT` ⇒ the ≥2★ auto-promote gate, + a monthly release instead of the cache's ClinVar 2025-02 | Cheapest win by an order of magnitude. |
| R2 | **SpliceAI slim (Δ≥0.1)** | ~0.6 GB | Deep-intronic + exonic-synonymous splice — **a whole variant class the screen currently cannot nominate** | Highest clinical value. Filter the free Ensembl MANE mirror; lossless (keep-only rule). SNV-only — no indel file exists on the free mirror. Guard contig naming (`1` vs `chr1`): a mismatched tabix query returns empty at **exit 0**. |
| R3 | **gnomAD slim (5 of 664 INFO fields)** | ~10 GB | True `faf95` (restores the CI correction) + `nhomalt` (the recessive false-positive tell) | Stream-slim the 24 joint chrom VCFs; nothing but the slim lands. GCS egress is free. Requires htslib built with libcurl — check `samtools --version` (not `bcftools --version`, which prints no feature line). |
| R4 | **Dedicated REVEL (+ AlphaMissense)** | ~1.3 GB | Reporting/tiering columns; the ClinGen-calibrated predictor for a PP3/BP4 step | **Buys the *screen* nothing** — see limitations.md §7. Use the dedicated files, **not** dbNSFP (30 GB for 5 columns, and its URL is dead: S3 `NoSuchBucket`, now registration-gated). Trap: Ensembl's `AlphaMissense.pm` emits `am_pathogenicity`, not `AlphaMissense_score`. |
| R5 | **LOFTEE data** | ~13 GB | HC/LC pLoF confidence ⇒ PVS1 strength grading | Plugin code already in the image; near-inert for *selection* (HIGH impact already keeps every pLoF), so this is a tiering prerequisite. |

Not a download, but on this list because it gates the same reasoning: **CADD's threshold is
off-label.** 25.3 is Pejaver-2022's *missense* PP3-supporting cutoff, and missense never reaches
the CADD rung — so it is applied only to the non-coding variants it was never calibrated for, with
no ClinGen-endorsed alternative to swap in. Region-stratified calibration is research, not a fetch.

## Dependency spine (build in this order)

1. **Per-sample ancestry + relatedness/swap QC (somalier)** — prerequisite for *any* honest
   gnomAD-as-control burden and for de-duplicating "recurrence." Already in the container image.
2. **Calibrated per-gene recurrence null + FDR** — converts the headline signal from a heuristic
   sort into a falsifiable statistic. Reuses machinery already present (`scipy.poisson.sf`, `bh_fdr`).
3. **CoCoRV external-control engine** — depends on #1; supersedes the documented-but-unbuilt TRAPD plan.

## Quick wins (high impact, low–med effort — do these first)

| # | Gap | Impact | Effort | Why now |
|---|-----|:------:|:------:|---------|
| 1 | **✅ DONE.** **Calibrated recurrence null + cross-gene FDR.** From each qualifying variant's gnomAD `faf95` + N_trios, compute expected carriers → binomial/Poisson tail → BH-q **per gene**. | High | Low | *The* defensibility gap: today 2 carriers at faf95≈1e-4 rank the same as 2 truly-private carriers, and long/mutable genes float up uncorrected. Reuses the de-novo arm's `poisson.sf`/`bh_fdr`. *(Step 6: `p_recurrence`/`q_recurrence`/`recurrence_exome_wide_sig`.)* |
| 2 | **somalier: per-sample ancestry (1KG/HGDP PCs) + cross-cohort relatedness/dup/swap + joint sex.** | High | Low | A swapped/dup proband fabricates recurrence; ancestry-mismatched faf95 mis-estimates rarity. Already imaged, unused. **Unblocks CoCoRV**; also fixes the fragile chrX-only sex check. |
| 3 | **✅ DONE.** **CHARR contamination screen (VCF-only freemix).** | High | Low | 1–3% contamination turns hom-ref→apparent-het, manufacturing false inherited hets / comp-het second hits. Computable from AD/DP/GQ already parsed; no BAM needed. *(Step 0: verifyBamID `FREEMIX` if `resources.selfsm_dir` set, else CHARR proxy; `contam_flag` folds into `overall_pass`.)* |
| 4 | **PP1/BS4 co-segregation points** (ingest parent affected status from PED col 6) **+ a variant-keyed meiosis ledger** to sum segregations across families. | High | Low | The one informative meiosis per trio is discarded today; the ledger turns the cohort into the extended pedigree a single trio lacks. |
| 5 | **UTRannotator** (5′UTR/uORF) VEP plugin. | High | Low | One-line VEP fix: uAUG-creating/uORF-disrupting variants in haploinsufficient CPS genes (NF1, RB1) are currently dropped as `not_functional`. Ships with VEP 115. |
| 6 | **UPD screen (UPDhmm/UPDio)** to *rescue* apparent-Mendelian-error homozygous recessives. | High | Low | The recessive logic currently deletes the UPD case (1/1 child + 0/0 parent → "Mendelian error"); paternal UPD(11p15) → ~20% of Beckwith-Wiedemann. |
| 7 | **AlphaMissense at calibrated Strong/Moderate** (2025 SVI); wire or drop the dead MetaRNN field. | Med | Low | AlphaMissense is now SVI-endorsed on par with REVEL — reaches Strong on constrained genes where REVEL sits at Supporting. |
| 8 | **conda-lock lockfile + `@sha256`-pinned base image + CI drift gate.** | Med | Low | The Dockerfile pins a mutable tag + unlocked `>=` specs; a silent htslib/numpy bump can change normalization/tiering silently. Docs already prescribe this. |

*Ride-free with the above:* F_ROH consanguinity prior (on #12), phase-confidence field (on #11),
robust sex-check (on #2).

## Strategic (high impact, high effort — the big bets)

| # | Gap | Impact | Effort | Why now |
|---|-----|:------:|:------:|---------|
| 9 | **Germline CNV calling (GATK-gCNV) + ACMG/ClinGen dosage annotation (AnnotSV/ClassifyCNV) + CNV-in-trans into the comp-het resolver** (ExomeDepth as a concordance second caller). | High | High | Largest true blind spot: 10–15% of pediatric-CPS diagnoses are CNV/SV (single-exon RB1/SMARCB1/DICER1/NF1/PMS2 deletions). **Needs re-accessing Kids First CRAMs + a ≥100–150-sample batch** — a new data-model dependency. |
| 10 | **Phenotype ranker (Exomiser hiPHIVE) as an additive prior**, tuned per the 2025 optimization (human-only associations; REVEL+AlphaMissense+SpliceAI blend we already compute) + graceful sparse/absent-HPO degradation. | High | High | Phenotype-blind ranking buries the true diagnosis when every proband carries many rare functional variants. **Ship tuned** — mis-tuned Exomiser underperforms its own baseline. Needs HPO ingestion. |
| 11 | **Read-backed + population phasing** (WhatsHap where reads span; **gnomAD variant co-occurrence** for distant pairs) with a **three-class phase output** (trans / cis-reject / unknown-review). | High | Med–High | Comp-het is only called for the mat×pat case today; parent-of-origin-only cis pairs inflate carrier counts. gnomAD co-occurrence is a public lookup (no data cost); read-backed needs BAMs (couples with #9). |
| 12 | **Runs-of-homozygosity / homozygosity mapping (AutoMap / `bcftools roh`)** to prioritize homozygous-recessive candidates inside ROH tracts. | High | Med | ~50% reduction in candidate hom variants with 92.5% of causal variants inside ROH; a per-proband signal currently discarded. Unlocks F_ROH priors; complements UPD (#6). |

*Strategic follow-ons (gated):* **CoCoRV** external-control burden (ancestry-stratified CMH,
empirical-null λ, discrete-aware FDR) — gated on #2, natural phase-2 after #1. **Somatic
second-hit / LOH overlay** from matched Kids First tumor data (PBTA/OpenPedCan) — ~1/3 of true
carriers show a second hit; high impact but a new matched-tumor pipeline.

> **Shared dependency.** #2 (somalier), #9 (gCNV), #11 (read-backed phasing), and the true
> verifyBamID upgrade to #3 all gate on one thing: re-accessing the source **CRAMs**. They are
> grouped — with ordering, cost, and the CoCoRV unlock — in
> **[docs/cram_access_phase.md](cram_access_phase.md)**. Pay the access cost once.

## Nice-to-have

- **ACAT-O / Cauchy multi-mask omnibus** (Med/Low) — hard-depends on #1's per-mask p-values.
- **GIAB HG002/CMRG benchmarking harness** (hap.py/vcfeval) — impacts *credibility* not yield; run
  *after* the caller set stabilizes (post-CNV).
- **Gene-list tier priors** (OMIM/PanelApp/ClinGen/ACMG-SF v3.3/CGC) + MOI-consistency cross-check.
- **Constitutional-mosaic tier** (VAF 0.03–0.30, beta-binomial vs DP) — the rigid 0.25 AB floor drops
  mosaic TP53/NF1; shares VAF machinery with a **CHIP confounder flag**.
- **Extended-window splicing** (SpliceVault / Pangolin) for deep-intronic/cryptic pseudoexons.
- **Age-dependent penetrance model** (gated on #4); **per-gene predictor calibration override table**
  (ClinGen VCEP cutoffs); **PMS2/PMS2CL paralog resolution** (short-read partial rescue);
  **STR/ExpansionHunter** (low CPS relevance — neuro/general-rare-disease only).

## Sequencing traps
- ACAT-O needs #1's p-values; CoCoRV needs #2's ancestry; F_ROH and robust sex-check ride free on
  #12 and #2 (don't build standalone); GIAB benchmarking must come *after* the caller set stabilizes.

## Key sources
- Recurrence null / external control: TRAPD (Guo AJHG 2018, [PMC6174288](https://pmc.ncbi.nlm.nih.gov/articles/PMC6174288/)); CoCoRV (Chen Nat Commun 2022, [PMC9095601](https://pmc.ncbi.nlm.nih.gov/articles/PMC9095601/)).
- somalier (Pedersen Genome Med 2020, [PMC7362544](https://pmc.ncbi.nlm.nih.gov/articles/PMC7362544/)); CHARR (Lu AJHG 2023, [PMC10716339](https://pmc.ncbi.nlm.nih.gov/articles/PMC10716339/)).
- PP1/BS4 (ClinGen SVI, AJHG 2024, [PMC10806742](https://pmc.ncbi.nlm.nih.gov/articles/PMC10806742/)); UTRannotator (Zhang Bioinformatics 2021; Whiffin Nat Commun 2020).
- UPDhmm ([Bioinformatics 2026](https://academic.oup.com/bioinformatics/article/42/3/btag062/8529595)); ROH/AutoMap (Quinodoz Nat Commun 2021; [PMC7477492](https://pmc.ncbi.nlm.nih.gov/articles/PMC7477492/)).
- Exomiser optimization ([Genome Med 2025](https://link.springer.com/article/10.1186/s13073-025-01546-1)); GATK-gCNV (Babadi Nat Genet 2023, [PMC10904014](https://pmc.ncbi.nlm.nih.gov/articles/PMC10904014/); exome-CNV yield 2.6% AJHG 2024).
- WhatsHap ([docs](https://whatshap.readthedocs.io/)); CMRG/GIAB (Wagner Nat Biotechnol 2022; GA4GH stratifications); AlphaMissense (Cheng Science 2023).
