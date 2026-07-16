# Clinical Pathogenicity: ClinVar & ACMG/AMP

How this pipeline uses ClinVar as clinical evidence today (an unstarred P/LP keep-override), and the review-gated, points-based ACMG/AMP classification it is designed to grow into.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

> ### ⚠ Status: almost all of this document is TARGET, not what runs
>
> Under the **VEP-only contract** (a VEP 115 GRCh38 cache + the CADD plugin; no ClinVar VCF is
> downloaded or bind-mounted), the *only* clinical evidence the pipeline reads is the cache's
> `CLIN_SIG` string. That has two hard consequences for this document:
>
> - **There is no `CLNREVSTAT`, therefore no stars.** The ≥ 2★ auto-promote gate that this
>   document treats as central **cannot currently be applied at all** — it is not relaxed, it is
>   unimplementable. Unstarred P/LP (excluding `conflicting`) is honored instead, which
>   **over-retains** rather than over-drops. Restoring stars costs **~0.18 GB** — the cheapest
>   item in [limitations.md](limitations.md).
> - **No ACMG classification is computed.** AutoGVP, the Tavtigian points backbone, graded PVS1,
>   PM2/PP3/BP4 strengths, and the ClinGen gene–disease validity gate are **not implemented**.
>   The pipeline is a *screen*: it assigns no ACMG weight and emits no P/LP/VUS call of its own.
>
> The science below is retained deliberately — it is the tiering roadmap and the justification
> for the resource spend. Sections that describe intent rather than behaviour are marked
> **TARGET**. For the full ledger of what the first pass cannot see, see
> [limitations.md](limitations.md); the authority on what runs is
> [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- **IMPLEMENTED (all of it):** ClinVar `CLIN_SIG` is read from the **VEP cache** (`--check_existing`; lowercase, `&`-joined). A **P/LP** assertion — `conflicting` excluded — is a **keep override** at Step 3: it rescues a variant that fails the rarity or functional screen (`hprv_keep_reason=clinvar_plp`). **BA1-common (AF ≥ 0.05) is never rescued**, P/LP or not. The raw `clnsig` string is carried through to `candidates.calls.tsv` and the IGV export for the curator. That is the entire clinical layer.
- **NOT implemented — no stars.** The cache carries no `CLNREVSTAT`, so a 1★ single-submitter P/LP is indistinguishable from a 3★ expert-panel call and both are honored. **Release is pinned by the cache** (VEP 115 ⇒ **ClinVar 2025-02**), not independently — ClinVar itself ships monthly, so the assertions are stale by construction.
- *TARGET* — ClinVar star mapping: **4★** practice guideline, **3★** expert panel (a ClinGen VCEP call overrides all submitters), **2★** multiple submitters/no conflicts, **1★** single submitter or conflicting, **0★** no assertion criteria. Needs a ClinVar VCF.
- *TARGET* — auto-promote **P/LP at ≥ 2★** (no conflicts); 1★ P/LP → prioritize + human review; `Conflicting_classifications` / `Uncertain_significance` (VUS) → flag, never auto-promote; **0★ excluded** from auto-logic.
- *TARGET* — **ACMG classifier backbone = AutoGVP** (CHOP/Kids First; ClinVar + modified InterVar with graded PVS1 and PP5/BP6 removed) — purpose-built for this GMKF pediatric-cancer / rare-disease GRCh38 use case. Nothing of it is wired in.
- *TARGET* — **Combining = Tavtigian/ClinGen Bayesian points**: **P ≥ 10, LP 6–9, VUS 0–5, LB −1…−5, B ≤ −6** (Supporting ±1, Moderate ±2, Strong ±4, Very Strong ±8).
- *TARGET* — **PM2 at Supporting strength only** — it is ACMG *evidence*, distinct from the pipeline's rarity screening gate; do not conflate "passed the rarity filter" with "PM2 met." (The screening gate exists; PM2 does not.)
- *TARGET* — **PVS1** at **graded strength** via the Abou Tayoun 2018 decision tree, gated on ClinGen gene–disease validity **≥ Moderate** and a known loss-of-function mechanism — never naive full-strength. Blocked twice over: no LOFTEE and no gene–disease validity table.
- *TARGET* — **Missense in-silico:** report **one** Pejaver-2022–calibrated predictor per variant (primary **REVEL**), never stack correlated predictors. No missense predictor is available, and note it would change **no keep decision** — see [functional_annotation.md](functional_annotation.md) and [limitations.md](limitations.md) §7.
- *TARGET* — **Gene-level gate:** restrict high-priority auto-calls to ClinGen **Definitive / Strong / Moderate** gene–disease validity; a gene-specific VCEP threshold overrides any generic cutoff. No such gate runs; nothing is currently auto-promoted for it to gate.

## ClinVar as an evidence source

### Record model and review status

ClinVar aggregates submitted records (SCVs) into variant-level (VCV) and variant/condition (RCV) records, and assigns each a **gold-star review status (0–4)** reflecting the strength of the submitting evidence. Since January 2024 ClinVar splits classification into three axes — **germline**, **somatic clinical impact**, and **oncogenicity** — each with its own VCF INFO tags (`CLNSIG`/`CLNREVSTAT`, `ONC*`/`ONCREVSTAT`, `SCI*`/`SCIREVSTAT`). This pipeline uses the **germline** axis for rare-disease and germline pediatric-cancer screening.

> **What the pipeline actually sees.** All of the above describes the **ClinVar VCF**, which this
> pipeline does not read. It reads VEP's cached `CLIN_SIG`, which is a **classification string
> only**: no `CLNREVSTAT`, no `CLNSIGCONF`, no star, no submitter breakdown, no axis split. The
> star table immediately below is therefore reference material for the TARGET design — **none of
> it can be evaluated at runtime today.** The one distinction the code can and does make is
> `conflicting` (matched as a substring of either the old
> `conflicting_interpretations_of_pathogenicity` or the current
> `conflicting_classifications_of_pathogenicity`), which is excluded from the P/LP override.

| Stars | `CLNREVSTAT` token | Meaning |
|------:|--------------------|---------|
| 4★ | `practice_guideline` | Professional practice guideline |
| 3★ | `reviewed_by_expert_panel` | ClinGen VCEP / expert-panel call — **overrides all submitters** |
| 2★ | `criteria_provided,_multiple_submitters,_no_conflicts` | Concordant multi-lab |
| 1★ | `criteria_provided,_single_submitter` | Single lab, criteria provided |
| 1★ | `criteria_provided,_conflicting_classifications` | Conflicting (inspect `CLNSIGCONF`) |
| 0★ | `no_assertion_criteria_provided` / `no_classification_provided` | No criteria — excluded from auto-logic |

*(TARGET table — requires `CLNREVSTAT` from a ClinVar VCF; unavailable under the VEP-only contract.)*

`CLNSIG` VCF tokens use underscores: `Pathogenic`, `Likely_pathogenic`, `Pathogenic/Likely_pathogenic`, `Uncertain_significance`, `Likely_benign`, `Benign`, and `Conflicting_classifications_of_pathogenicity` (formerly `Conflicting_interpretations_of_pathogenicity`), plus low-penetrance / risk-allele and non-standard terms (`drug_response`, `association`, `protective`, `Affects`).

### Consumption rules — IMPLEMENTED

What Step 3 does today (`src/hprv/selection.py`, `src/hprv/annotations.py:clnsig_is_plp`), in order:

- **BA1 first, and it wins.** AF ≥ 0.05 → dropped, **before** ClinVar is consulted. A P/LP assertion on a common allele does **not** rescue it. This is the one place a ClinVar P/LP label is deliberately overruled, and it is the correct direction: BA1 is stand-alone benign evidence.
- **Extract P/LP** by matching `pathogenic` in `CLIN_SIG`, excluding any string containing `conflicting`, `likely_benign`, or `benign/likely`. Both spellings are matched — VEP's lowercase `&`-joined form (`pathogenic&likely_pathogenic`) and the ClinVar VCF's capitalised `/`-joined form (`Pathogenic/Likely_pathogenic`) — so the predicate survives a future switch to a real ClinVar VCF.
- **P/LP is a keep override, not a promotion.** It rescues a variant from the rarity gate (`too_common`) and from the functional ladder, tagging it `hprv_keep_reason=clinvar_plp`. It does not assign a class, a tier, or ACMG weight.
- **No star gate — it is unimplementable, not relaxed.** Every P/LP is honored regardless of review status, including 0★/1★ single-submitter assertions that the TARGET design would exclude or route to review. Effect: **over-retention** (curation load), not missed calls. gnomAD's own guidance — filter to ≥ 1★ with a specified classification — cannot be followed here.
- **VUS / Conflicting are never promoted** (they do not match the predicate) and are **never dropped on ClinVar grounds** — they simply face the normal rarity + functional screen, and the `clnsig` string rides along into `candidates.calls.tsv` and the IGV export. Note the difference from the TARGET: the string is *visible to a curator*, but nothing **routes** it to review, and `CLNSIGCONF` (the submitter breakdown) does not exist to inspect.
- **The release is pinned by the cache, not by us.** VEP 115's GRCh38 cache carries **ClinVar 2025-02**; the ClinVar VCF ships monthly. So provenance is recorded (via the cache version) but staleness is not controllable independently of a VEP upgrade. ClinVar reclassifies continuously and can carry outdated or single-lab calls: treat P/LP as a **triage prior, never an answer**, and plan periodic re-classification of previously reported VUS.

### Consumption rules — TARGET (needs a ClinVar VCF; ~0.18 GB)

- **Auto-promote at ≥ 2★** with no conflicts. **1★ P/LP** is prioritized but routed to human review.
- **Exclude 0★** from auto-logic (gnomAD guidance: keep ≥ 1★ with a specified classification).
- Inspect the `CLNSIGCONF` submitter breakdown for conflicting records.
- Pin a **dated release** independent of the VEP cache and record the release date in run provenance.

The recipe below is the TARGET wiring — **it is not in `02_annotate_sites.sh`**, which performs no
`bcftools annotate` transfers at all. It is the shape re-enabling would take: one transfer, plus
the new INFO fields in `annotations.F`.

```bash
# TARGET (not implemented): annotate with a version-pinned ClinVar release (GRCh38).
# clinvar_YYYYMMDD.vcf.gz is the dated release; record the date in run metadata.
bcftools annotate \
  -a "${CLINVAR_VCF}" \
  -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNSIGCONF \
  -Oz -o "${OUT_VCF}" "${IN_VCF}"

# High-confidence auto-promotable P/LP: P or LP in CLNSIG AND >=2-star, no conflicts.
# Today only the first clause is expressible, against vep_CLIN_SIG rather than CLNSIG.
bcftools view -i \
  'INFO/CLNSIG ~ "Pathogenic" &&
   INFO/CLNREVSTAT ~ "multiple_submitters" &&
   INFO/CLNREVSTAT !~ "conflicting"' \
  "${OUT_VCF}"
```

## ACMG/AMP framework and ClinGen SVI refinements

> **TARGET — none of this section is implemented.** The pipeline computes no ACMG criteria, no
> points total, and no class. It emits a *keep reason* (`impact_high` | `impact_moderate` |
> `cadd` | `clinvar_plp`), which is screening provenance, not evidence weight. The section is
> retained as the specification for the tiering step, and because it is the reason the
> resource-restoration items in [limitations.md](limitations.md) are worth their cost.

### The 2015 framework, points-based

The 2015 ACMG/AMP framework (Richards et al.) defines 28 criteria (16 pathogenic PVS1–PP5, 12 benign) combined by verbal rules into P / LP / VUS / LB / B. ClinGen's Sequence Variant Interpretation (SVI) working group has progressively replaced the verbal combining rules with the **Tavtigian Bayesian points system**, in which each criterion contributes exponentially-scaled points and the sum determines the class. This lets a criterion be applied at a **tunable strength** rather than a fixed weight.

| Strength | Pathogenic points | Benign points |
|----------|------------------:|--------------:|
| Supporting | +1 | −1 |
| Moderate | +2 | −2 |
| Strong | +4 | −4 |
| Very Strong | +8 | — |

| Total points | Classification |
|-------------:|----------------|
| ≥ 10 | Pathogenic |
| 6 – 9 | Likely Pathogenic |
| 0 – 5 | Uncertain Significance (VUS) |
| −1 … −5 | Likely Benign |
| ≤ −6 | Benign |

### Key SVI refinements the tiering step would apply (TARGET)

- **PM2 → Supporting by default** (SVI Recommendation v1.0). Absence/rarity is weak evidence; do not apply PM2 at Moderate. PM2 is ACMG *evidence* and must be kept distinct from the upstream rarity **screening** gate (see [allele_frequency.md](allele_frequency.md)) — passing the rarity filter is not the same as "PM2 met." Two extra cautions specific to this contract: the screening gate reads a **point estimate, not `faf95`**, and the cache reports frequencies only for **dbSNP-accessioned** alleles, so "absent" is weaker evidence here than it looks. PM2 built naively on this field would inherit both flaws.
- **PP3 / BP4 calibration (Pejaver 2022)**: continuous predictors receive strength-stratified thresholds. Use **one** predictor per variant; do not sum correlated predictors as independent evidence. The **REVEL** primary defaults are restated below and detailed in [functional_annotation.md](functional_annotation.md). AlphaMissense is *not* part of the Pejaver 2022 calibration and lacks an SVI PP3 stratification — treat it only as orthogonal support. **No calibrated predictor is available today** (no REVEL/AlphaMissense/MPC); CADD is present but its 25.3 cutoff is a missense-derived number applied only to non-coding variants, so it is a **discovery rank, not PP3 evidence**.
- **PVS1 decision tree (Abou Tayoun 2018)**: loss-of-function variants receive graded strength (PVS1 / _Strong / _Moderate / _Supporting) by consequence, NMD escape, exon/region context (last exon, 3′-terminal 50 bp, single-exon), and only when LoF is the disease mechanism and gene–disease validity is **≥ Moderate**. Naive full-strength PVS1 is a major over-calling source. This is **blocked twice over**: no LOFTEE HC/LC confidence and no ClinGen gene–disease validity table. It is also the single strongest argument for restoring LOFTEE, whose value to *selection* is otherwise near zero ([limitations.md](limitations.md) §5). See [functional_annotation.md](functional_annotation.md) and [gene_constraint.md](gene_constraint.md).
- **PS3 / BS3 functional (Brnich 2020)**: assay-based strength is set via OddsPath from the number of validated controls (≥ 11 controls → Moderate; more → Strong).

No single consolidated 2023–2025 replacement guideline is published yet; refinements continue to arrive as individual SVI recommendations and as gene-specific VCEP specifications, which take precedence over the generic rules whenever they exist.

## Automated ACMG classifiers

Rule-based ACMG automation varies widely in rigor; the common failure is applying **PVS1 at full strength** and using the reputable-source criteria PP5/BP6, both of which inflate P/LP calls.

| Tool | Basis | Caveat |
|------|-------|--------|
| **AutoGVP** (CHOP / Kids First / NCI) | ClinVar + **modified InterVar** (graded PVS1, PP5/BP6 removed), dockerized R workflow | Purpose-built for this GMKF pediatric-cancer / rare-disease GRCh38 use case — the reference backbone adopted here |
| InterVar (Wang lab) | Rule-based on ANNOVAR | Applies **PVS1 at full strength**, uses PP5/BP6 → P/LP over-calling; no VCEP specs |
| TAPES, GeneBe | Rule-based, no phenotype integration | Lower causal-variant prioritization in benchmarks |
| Franklin (Genoox) | Proprietary | Strong benchmarks but black-box; ToS constraints for bulk/container use |

**AutoGVP is the intended backbone — TARGET, not wired in.** It integrates a dated ClinVar release with a modified InterVar (PVS1-strength adjustment, PP5/BP6 removed) and is built for exactly this consortium use case, which is why it is the chosen target. Adopting it requires a dated ClinVar VCF at minimum, and its InterVar half needs the annotation set the VEP-only contract does not carry. Until then the pipeline runs **no** automated classifier — which also means it inherits none of the PVS1-full-strength / PP5-BP6 over-calling above. Universal caveat that applies to the target: automated calls are **screening aids, not diagnostic** — VUS and conflicts require human review, and none replace applicable VCEP rules.

## ClinGen gene–disease validity gate (TARGET)

ClinGen classifies each gene→disease relationship as **Definitive / Strong / Moderate / Limited / Disputed / Refuted / No Known Disease Relationship**. ACMG recommends diagnostic panels include only **Definitive / Strong / Moderate** genes.

**Not implemented.** No ClinGen gene–disease validity table is ingested; no gene-level gate runs, and nothing is auto-promoted for such a gate to act on. In the target design it would gate *auto-promotion* and PVS1 applicability (≥ Moderate validity plus a known LoF mechanism) — never variant retention, consistent with the never-drop principle: curated lists act as priors, not hard filters (see [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md)). Its absence costs the screen nothing, since a gate on auto-promotion is vacuous without auto-promotion; it is a prerequisite for the tiering step, not for the first pass.

## Integration pattern (per-trio VCF)

**IMPLEMENTED — one annotation source, four keep reasons, no classification:**

```text
Per-trio VCF (GRCh38, GATK genotype-refined)
  └─ VEP 115 cache + CADD plugin  -> split-vep lifts CSQ to INFO   (functional_annotation.md)
     └─ Step 3 classify (src/hprv/selection.py):
        ├─ AF >= 0.05 (BA1)                        -> DROP (never rescued, P/LP included)
        ├─ CLIN_SIG P/LP, not conflicting          -> KEEP  'clinvar_plp'   (unstarred!)
        ├─ AF >= 1e-2 and not P/LP                 -> DROP 'too_common'
        ├─ IMPACT in {HIGH, MODERATE}              -> KEEP  'impact_high'|'impact_moderate'
        ├─ CADD_PHRED >= 25.3                      -> KEEP  'cadd'  (non-coding only, in practice)
        └─ otherwise                               -> DROP 'not_functional'
```

*TARGET — the tiering step this document specifies. Every join below is currently absent:*

```text
Per-trio VCF (GRCh38, GATK genotype-refined)
  └─ VEP: canonical / MANE consequence            (functional_annotation.md)
     └─ join dated ClinVar VCF: CLNSIG, CLNREVSTAT, CLNSIGCONF     [absent: cache CLIN_SIG only]
        └─ join gnomAD v4.1 joint grpmax faf95     (allele_frequency.md)   [absent: point-estimate proxy]
           └─ REVEL (single predictor)             (functional_annotation.md)  [absent]
              └─ AutoGVP: ClinVar + modified InterVar                     [absent]
                 └─ Tavtigian points  (PVS1 graded per Abou Tayoun; PM2 → Supporting)  [absent]
                    └─ gate on ClinGen gene–disease validity ≥ Moderate   [absent]
                       └─ germline VCEP specs override generic cutoffs    [absent]
```

Because Kids First trios are **GATK genotype-refinement output (posterior/PP-refined), not cohort joint-genotyped**, this pipeline does **not** derive internal-cohort allele frequency. Population rarity comes from gnomAD v4.1 — but as **cached point-estimate AFs read through VEP**, not `faf95` from the gnomAD sites VCF (see [allele_frequency.md](allele_frequency.md) and [limitations.md](limitations.md) §2); de-novo confidence comes from the trio PP/GQ, not cohort frequency (see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)). Germline pediatric-cancer specifics — dominant vs recessive CPS handling, second-hit boosts, PMS2/PMS2CL — are in [pediatric_cancer.md](pediatric_cancer.md).

## Recommended defaults (this pipeline)

| Parameter | Default | Status | Notes |
|-----------|---------|--------|-------|
| ClinVar source | VEP 115 cache `CLIN_SIG` (`--check_existing`) | **IMPLEMENTED** | Release pinned by the cache ⇒ **ClinVar 2025-02**; not independently datable |
| ClinVar P/LP | Keep override: `pathogenic` matched, `conflicting`/benign excluded | **IMPLEMENTED** | Rescues from `too_common` + the functional ladder; tags `clinvar_plp`. **No class assigned** |
| Star gate | — | **RETIRED (unimplementable)** | No `CLNREVSTAT` in the cache ⇒ unstarred P/LP is honored. Over-retains; ~0.18 GB to restore |
| BA1 vs P/LP | AF ≥ 0.05 drops **even if P/LP** | **IMPLEMENTED** | BA1 is evaluated first and is never rescued |
| Conflicting / VUS | Not promoted, not dropped on ClinVar grounds; `clnsig` emitted to output | **IMPLEMENTED (partial)** | Never-drop holds; but no review **routing**, and no `CLNSIGCONF` to inspect |
| ClinVar release | Dated, version-pinned; recorded in provenance | *TARGET* | Needs a ClinVar VCF. Re-classify prior VUS periodically |
| Auto-promote P/LP | `CLNSIG` P or LP **and ≥ 2★**, no conflicts | *TARGET* | 1★ P/LP → prioritize + human review |
| 0★ records | Excluded from auto-logic | *TARGET* | gnomAD guidance: keep ≥ 1★ with a specified classification |
| ACMG classifier | **AutoGVP** (ClinVar + modified InterVar, PP5/BP6 removed) | *TARGET* | Graded PVS1. Nothing wired in |
| Combining rule | Tavtigian points: **P ≥ 10, LP 6–9, VUS 0–5, LB −1…−5, B ≤ −6** | *TARGET* | Supporting ±1 / Moderate ±2 / Strong ±4 / Very Strong ±8 |
| PM2 | **Supporting** strength only | *TARGET* | ACMG evidence, distinct from the rarity screening gate (which exists; PM2 does not) |
| PVS1 | Graded (Abou Tayoun 2018); gene validity ≥ Moderate + known LoF | *TARGET* | Blocked twice: no LOFTEE, no validity table |
| PP3/BP4 missense | **REVEL** (single predictor): PP3 supporting ≥ 0.644, moderate ≥ 0.773, strong ≥ 0.932; BP4 supporting ≤ 0.290, moderate ≤ 0.183 | *TARGET* | No missense predictor available. Do not stack; ranges from Pejaver 2022 |
| Gene gate | ClinGen **Definitive / Strong / Moderate** validity | *TARGET* | VCEP threshold overrides generic cutoff |

Values marked IMPLEMENTED are **configurable defaults** (`config/config.example.yaml`), not immutable law; values marked *TARGET* are specifications with no code behind them yet. A gene-specific ClinGen VCEP specification (its own PVS1 strength, PM2/BA1/BS1 thresholds, or PP3 calibration) **overrides** any generic default here — a rule the target tiering step must honor, and which nothing in the first pass currently implements. If this table ever disagrees with [Canonical defaults](README.md#canonical-defaults), that table wins.

## Scope limitations (stated honestly)

The contract-level gaps are catalogued once in **[limitations.md](limitations.md)** — read it before
interpreting a negative result. The two that bear directly on this document:

- **No review status ⇒ no star gate** ([limitations.md](limitations.md) §6). The gate this document
  treats as the centrepiece of ClinVar consumption cannot be applied at all. Every P/LP is honored
  at face value, so **a `clinvar_plp` keep is not evidence of curated confidence** — it may be one
  submitter with no criteria. Curators must check the star level in ClinVar's web UI by hand. The
  direction of error is over-retention, and the fix is the cheapest on the ledger (~0.18 GB).
- **ClinVar is stale by construction.** Pinned to the cache (VEP 115 ⇒ **2025-02**) rather than
  ClinVar's monthly release. A variant reclassified since is read at its old assertion — in both
  directions: a since-downgraded P/LP still overrides the screen, and a since-upgraded VUS gets no
  override.

And the pre-existing scope limits:

- **SNV/indel only initially.** The ClinVar override here acts on SNVs and indels (no ACMG classification is computed for any variant class). CNV/SV pathogenic calls (e.g. single-exon *RB1*/*SMARCB1*/*DICER1*/*NF1* deletions, *PMS2* rearrangements — 10–15% of pediatric-cancer and rare-disease diagnoses) are a known blind spot pending a future GATK-gCNV / Manta / ExomeDepth module.
- **Pseudogene / segmental-duplication genes** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) are low-confidence from short reads; ClinVar/ACMG calls in those regions must be treated as suspect regardless of star level (and star level is unavailable here anyway). Note the pipeline does **not** currently flag or mask these regions — [limitations.md](limitations.md), structural gaps.
- **Proband post-zygotic mosaicism** (e.g. *NF1*, overgrowth): low-VAF calls fall outside the standard heterozygous allele-balance band and can be missed before classification is even reached (see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)).
- **Genotype-refinement caveat:** gnomAD priors in GATK `CalculateGenotypePosteriors` can push a genuine ultra-rare pathogenic call toward hom-ref — for top ClinVar/ACMG candidates, cross-check the pre-refinement PL/GT.
- **Calibration/validation:** classification sensitivity/precision should be measured against GIAB/CMRG truth sets and a positive-control variant panel; no truth-set benchmark is wired in yet. Screening sensitivity is therefore **unmeasured, not measured-and-acceptable**.

## Sources

- ClinVar review status & stars: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
- ClinVar classification / `CLNSIG`: https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
- gnomAD ClinVar review-status filter guidance: https://gnomad.broadinstitute.org/news/2023-09-clinvar-variants-filter-by-review-status/
- gnomAD v4.1: https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/ ; ClinGen guidance on gnomAD v4 for VCEPs (Mar 2024): https://clinicalgenome.org/site/assets/files/9445/clingen_guidance_to_vceps_regarding_the_use_of_gnomad_v4_march_2024.pdf
- Tavtigian points (Bayesian ACMG): https://onlinelibrary.wiley.com/doi/10.1002/humu.24088 (DOI 10.1002/humu.24088)
- ACGS 2023 UK variant-classification guidelines (points / PM2): https://www.acgs.uk.com/media/12443/uk-practice-guidelines-for-variant-classification-v1-2023.pdf
- Pejaver 2022 PP3/BP4 calibration: https://www.cell.com/ajhg/pdfExtended/S0002-9297(22)00461-X (DOI 10.1016/j.ajhg.2022.10.013)
- PVS1 decision tree, Abou Tayoun 2018: https://onlinelibrary.wiley.com/doi/abs/10.1002/humu.23626 (DOI 10.1002/humu.23626)
- PS3/BS3 functional-assay strength, Brnich 2020: https://link.springer.com/article/10.1186/s13073-019-0690-2 (DOI 10.1186/s13073-019-0690-2)
- ClinGen gene–disease validity framework: https://clinicalgenome.org/docs/evaluating-the-clinical-validity-of-gene-disease-associations-an-evidence-based-framework-developed-by-the-clinical-genome/
- AutoGVP (CHOP / Kids First / NCI): https://academic.oup.com/bioinformatics/article/40/3/btae114/7616989 (DOI 10.1093/bioinformatics/btae114)
- Automated ACMG classifier benchmark: https://academic.oup.com/bioinformatics/article/42/2/btaf623/8483023
