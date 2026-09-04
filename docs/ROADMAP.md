# Roadmap — high-priority gaps vs state of the art

From a 10-domain SOTA literature review (2023–2026) against the current repo. Items are ranked by
**impact × defensibility-per-unit-effort** for the pipeline's headline claim: *recurrent +
constrained inherited variants = candidates*. Effort is "to add here." This is a living document —
update it as items land. (Out of scope by design: de novo review, mtDNA — separate pipelines.)

## Restoring the VEP-centric contract's gaps (mostly done)

The pipeline runs on a VEP 115 cache + the CADD/SpliceAI/REVEL/AlphaMissense plugins and two
`bcftools annotate` transfers (the ClinVar sites VCF and the gnomAD joint slim); see
**[limitations.md](limitations.md)** for
the full ledger of what that costs and why each was an acceptable trade. Each row below is
**additive**: one `bcftools annotate` transfer in `02_annotate_sites.sh` plus its INFO field in
`annotations.F`. Nothing in the architecture blocks any of them. Ranked by clinical value per GB.

| # | Restore | Size | Buys | Notes |
|---|---------|:----:|------|-------|
| R1 | **✅ DONE — ClinVar VCF** | ~0.18 GB | `CLNREVSTAT` ⇒ `clinvar_stars`, + an independently pinned release beside the cache's ClinVar 2025-02 | Transferred in Step 2. Stars RANK in Step 9 (positive limb damped below `min_review_stars`); the ≥2★ auto-promote *gate* this row once proposed is deliberately NOT built — it would violate never-drop. |
| R2 | **✅ DONE — SpliceAI** | ~28 GB | Deep-intronic + exonic-synonymous splice, as rung 2 of the functional ladder (Δ ≥ `spliceai_ds_min`, default 0.2) | Shipped BETTER than this row proposed: the FULL raw SNV+indel set as a VEP **plugin** (not a bcftools transfer, not the Δ≥0.1 MANE-only slim), **required by default** (`spliceai_required`). Remaining gap: the precomputed set covers only 1 nt insertions and deletions ≤ 4 nt, at ±50 nt — see limitations.md. Contig-naming guard still applies. |
| R3 | **✅ DONE — gnomAD slim (5 of 664 INFO fields)** | ~10 GB | True `faf95` (restores the CI correction) + `nhomalt` | Implemented: `prepare_resources.sh --only gnomad_sites fetch` stream-slims the 24 joint chrom VCFs (nothing raw lands; needs htslib+libcurl, checked via `samtools --version`), Step 2 transfers them, and `annotations.frequency()` reads faf95 as the ONE oracle for the run (`resources.gnomad.oracle: faf95`, the default — required, halts without the slim); `grpmax_proxy` is the deliberate opt-down and the arms never cross. `rarity_oracle` is recorded once per run, `rarity_basis` (`measured`/`zero_ci`/`absent`) per variant. FAF groups verified as afr/amr/eas/mid/nfe/sas. |
| R4 | **✅ DONE — Dedicated REVEL + AlphaMissense** | ~1.3 GB | Step 9's missense tier (REVEL → AlphaMissense → CADD, a fixed precedence); required by default | **Buys the *screen* nothing** — see limitations.md §7. Wired as VEP plugins over the dedicated files, **not** dbNSFP (30 GB for 5 columns, and its URL is dead: S3 `NoSuchBucket`, now registration-gated). Trap handled: Ensembl's `AlphaMissense.pm` emits `am_pathogenicity`, not `AlphaMissense_score`. |
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
| 1 | **✅ DONE (as a rank).** **Case-only recurrence null + cross-gene FDR, plus the mutational-target carrier expectation.** From each qualifying variant's `rarity_af` + N_trios, compute expected carriers → binomial tail → BH-q **per gene**; and, with a mutational-target table, `exp_carriers_mu`/`p_carrier_excess` (`burden.rank_by_mutational_target`). | High | Low | *The* defensibility gap: 2 carriers at faf95≈1e-4 once ranked the same as 2 truly-private carriers, and long/mutable genes floated up uncorrected. The binomial null is a case-only RANK — it saturates for private variants, so never call it calibrated; the mutational-target expectation is what stops long genes leading by size. *(Step 6: `p_recurrence`/`q_recurrence`/`recurrence_exome_wide_sig`, `p_carrier_excess`.)* |
| 2 | **somalier: per-sample ancestry (1KG/HGDP PCs) + cross-cohort relatedness/dup/swap + joint sex.** | High | Low | A swapped/dup proband fabricates recurrence; ancestry-mismatched faf95 mis-estimates rarity. Already imaged, unused. **Unblocks CoCoRV**; also fixes the fragile chrX-only sex check. |
| 3 | **◐ PARTIAL.** **Contamination screen (verifyBamID FREEMIX + a VCF-only raw proxy).** | High | Low | 1–3% contamination turns hom-ref→apparent-het, manufacturing false inherited hets / comp-het second hits. **FREEMIX path is production-ready**; the VCF-only fallback is a raw (uncorrected) ref-read fraction — a CHARR-*like* proxy, **not** the calibrated Lu-2023 statistic — that flags only gross (≳5–8%) contamination, so it does **not** yet catch the 1–3% band. *(Step 0: verifyBamID `FREEMIX` if `resources.selfsm_dir` set, else the raw proxy; `contam_flag` folds into the **advisory** `overall_pass`, which no step auto-excludes yet.)* **TODO:** a corrected CHARR (per-genotype mean `/mean(1−AF)`, baseline-subtracted, threshold re-derived from spike-ins) post-annotation, and/or config-gated exclusion of flagged trios from the recurrence tally. |
| 3b | **◐ PARTIAL.** **Null calibration diagnostic + a graded artifact down-weight (Step 9).** | High | Med | *The* other half of the A-3 gap. **DONE:** a **mid-p calibration diagnostic** is computed for the fitted null and for the Poisson alternative on every run (`audit/counts.tsv` → `calibration.*`), so the NB-vs-Poisson choice is auditable rather than asserted — measured 2.41× Poisson anti-conservatism at α=1e-3 (arm trimmed; 1.82× untrimmed) vs 0.31× for the NB. Plus a per-gene **excess-over-mutational-target** statistic (NB2, trimmed fit, BH-FDR), a six-signal artifact panel, and a four-tier graded down-weight that triaged **9.57% of the candidate list at 100% established-gene retention** (default count floor; 10.36% at floor 1) — never a drop. *(Step 9: `variants.prioritized.tsv` / `genes.prioritized.tsv`; see [prioritization.md](prioritization.md).)* **TODO:** the **synonymous-λ** check and a positive-control **recovery** measurement on real data — neither needs a new resource, and λ ≫ 1 would mean every rank in the scheme inherits a filter bias. |
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
| 10 | **Phenotype ranker (Exomiser hiPHIVE) as an additive prior**, tuned per the 2025 optimization (human-only associations; of its REVEL+AlphaMissense+SpliceAI blend we currently compute only SpliceAI) + graceful sparse/absent-HPO degradation. | High | High | Phenotype-blind ranking buries the true diagnosis when every proband carries many rare functional variants. **Ship tuned** — mis-tuned Exomiser underperforms its own baseline. Needs HPO ingestion. |
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
- **◐ PARTIAL — Gene-list tier priors** (OMIM/PanelApp/ClinGen/ACMG-SF v3.3/CGC) + MOI-consistency
  cross-check. **DONE (Step 9):** the *mechanism* is built and it is deliberately Class-A —
  a gene-list prior enters as **one additive term (+2)** in a **second, separately reported
  ranking** (`rank_prior` beside `rank_agnostic`, with `rank_delta` exposing exactly which calls
  list membership promoted), it is **mechanism-gated** like constraint (zero at V0, so a list can
  never rescue a molecularly-benign prediction — that is how gene-list priors turn into
  confirmation bias), it is a **prior and never a filter**, and it **defaults OFF** with the list
  living in a config file path rather than in code, so hprv itself names no gene. The
  **MOI-consistency cross-check is also live** (`moi_coherence` ∈ `coherent`/`discordant`/`unknown`, with
  `unknown` scoring *exactly* 0 so novel genes are never punished, and the audit-A-6 long-gene
  comp-het drift suppressing the discordance penalty). **Done since:** the version-pinned list
  *content* — PanelApp 243 v5.12 / 259 v1.30 green and ACMG SF v3.3 were fetched and added to the
  established-gene union (2,218 genes); the 100%-retention figure did **not** move. Residual: ACMG
  SF v3.3 membership is recalled rather than machine-read — replace it with the published table
  before publication (see prioritization.md §6).
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
