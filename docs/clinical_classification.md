# Clinical Pathogenicity: ClinVar & ACMG/AMP

How this pipeline uses ClinVar as clinical evidence today (a star-blind P/LP keep-override at the screen, plus a star-damped ranking term in Step 9), and the review-gated, points-based ACMG/AMP classification it is designed to grow into.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

> ### ⚠ Status: the ACMG classification in this document is TARGET; the ClinVar layer is live
>
> Two clinical inputs run today. **ClinVar `CLIN_SIG`** comes from the VEP 115 cache
> (`--check_existing`) and drives the Step-3 P/LP keep override. **ClinVar review status**
> (`CLNREVSTAT` ⇒ `clinvar_stars` 0–4) arrives via the **ClinVar sites VCF**, the first of Step 2's
> two `bcftools annotate` transfers (`resources.clinvar.vcf`). That changes two things this
> document used to say:
>
> - **Stars exist — and they RANK, they do not gate.** The ≥ 2★ auto-promote/exclude gate this
>   document treats as central is **deliberately not reinstated at the screen**: a keep/drop gate
>   on review status would violate never-drop, so a 1★ single-submitter assertion is still kept
>   and reviewed. Stars are consumed in **Step 9**, where they damp the *positive* limb of the
>   clinical term below `resources.clinvar.min_review_stars` (2) by `low_star_scale` (0.5). Absent
>   stars (no transfer) read `UNAVAILABLE` = **full weight**, never 0★.
> - **No ACMG classification is computed.** AutoGVP, graded PVS1, PM2/PP3/BP4 *strengths* and the
>   ClinGen gene–disease validity gate are **not implemented**, and no P/LP/VUS label is ever
>   emitted. Step 9 does compute an additive **`priority_points`** composite in the *style* of
>   Tavtigian — but it is a triage rank, not an ACMG total, and [must never be read against the
>   P ≥ 10 / LP 6–9 / VUS 0–5 bands](#why-priority_points-is-not-an-acmg-total).
>
> The science below is retained deliberately — it is the tiering roadmap and the justification
> for the resource spend. Sections that describe intent rather than behaviour are marked
> **TARGET**. For the full ledger of what the first pass cannot see, see
> [limitations.md](limitations.md); the authority on what runs is
> [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- **IMPLEMENTED (all of it):** ClinVar `CLIN_SIG` is read from the **VEP cache** (`--check_existing`; lowercase, `&`-joined). A **P/LP** assertion — `conflicting` excluded — is a **keep override** at Step 3: it rescues a variant that fails the rarity or functional screen (`hprv_keep_reason=clinvar_plp`). **BA1-common (AF ≥ 0.05) is never rescued**, P/LP or not. The raw `clnsig` string is carried through to `candidates.calls.tsv` and the IGV export for the curator. That is the entire clinical layer.
- **IMPLEMENTED — stars, as a RANKING input.** The ClinVar sites VCF is transferred in Step 2 (`resources.clinvar.vcf` ⇒ `clinvar_CLNREVSTAT` ⇒ `clinvar_stars` 0–4). The **screen stays star-blind**; Step 9 damps the *positive* limb of its clinical term (+4) by `low_star_scale` (**0.5**) below `min_review_stars` (**2**). Only the positive limb: shrinking a low-star *benign* call toward zero would promote a poorly-reviewed benign assertion. **Blank ≠ 0★** — blank means the transfer did not run (`clinvar_review_status=UNAVAILABLE`, full weight); 0★ means ClinVar has a record whose submitter provided no assertion criteria.
- **ClinVar star mapping (IMPLEMENTED, `annotations.clinvar_stars`):** **4★** practice guideline, **3★** expert panel (a ClinGen VCEP call overrides all submitters), **2★** multiple submitters/no conflicts, **1★** single submitter or conflicting, **0★** no assertion criteria. An *unrecognised* status string returns `None` (unknown), never 0★.
- **Two releases, two dates.** `CLIN_SIG` is pinned by the **cache** (VEP 115 ⇒ **ClinVar 2025-02**); the transferred VCF is **independently version-pinned** and should be dated in run provenance. ClinVar ships monthly, so the cached assertions are stale by construction — which is also why the transfer carries `clinvar_CLNSIG` alongside, so a run can compare the two.
- *TARGET* — auto-promote **P/LP at ≥ 2★** (no conflicts) and **exclude 0★** from auto-logic. Deliberately **not** adopted at the screen (it would violate never-drop); it belongs to the tiering step, if one is ever built.
- *TARGET* — **ACMG classifier backbone = AutoGVP** (CHOP/Kids First; ClinVar + modified InterVar with graded PVS1 and PP5/BP6 removed) — purpose-built for this GMKF pediatric-cancer / rare-disease GRCh38 use case. Nothing of it is wired in.
- *TARGET* — **Combining = Tavtigian/ClinGen Bayesian points**: **P ≥ 10, LP 6–9, VUS 0–5, LB −1…−5, B ≤ −6** (Supporting ±1, Moderate ±2, Strong ±4, Very Strong ±8).
- *TARGET* — **PM2 at Supporting strength only** — it is ACMG *evidence*, distinct from the pipeline's rarity screening gate; do not conflate "passed the rarity filter" with "PM2 met." (The screening gate exists; PM2 does not.)
- *TARGET* — **PVS1** at **graded strength** via the Abou Tayoun 2018 decision tree, gated on ClinGen gene–disease validity **≥ Moderate** and a known loss-of-function mechanism — never naive full-strength. Blocked twice over: no LOFTEE and no gene–disease validity table.
- *TARGET* — **Missense in-silico:** report **one** Pejaver-2022–calibrated predictor per variant (primary **REVEL**), never stack correlated predictors. No missense predictor is available, and note it would change **no keep decision** — see [functional_annotation.md](functional_annotation.md) and [limitations.md](limitations.md) §7.
- *TARGET* — **Gene-level gate:** restrict high-priority auto-calls to ClinGen **Definitive / Strong / Moderate** gene–disease validity; a gene-specific VCEP threshold overrides any generic cutoff. No such gate runs; nothing is currently auto-promoted for it to gate.

## ClinVar as an evidence source

### Record model and review status

ClinVar aggregates submitted records (SCVs) into variant-level (VCV) and variant/condition (RCV) records, and assigns each a **gold-star review status (0–4)** reflecting the strength of the submitting evidence. Since January 2024 ClinVar splits classification into three axes — **germline**, **somatic clinical impact**, and **oncogenicity** — each with its own VCF INFO tags (`CLNSIG`/`CLNREVSTAT`, `ONC*`/`ONCREVSTAT`, `SCI*`/`SCIREVSTAT`). This pipeline uses the **germline** axis for rare-disease and germline pediatric-cancer screening.

> **What the pipeline actually sees — from TWO sources.** The **VEP cache** supplies `CLIN_SIG`, a
> classification string only (no review status, no `CLNSIGCONF`, no submitter breakdown, no axis
> split); the **ClinVar sites VCF transfer** supplies `CLNREVSTAT` ⇒ `clinvar_stars` and a second
> copy of `CLNSIG` (as `clinvar_CLNSIG`, reporting-only, so a run can see the pinned release
> disagree with the cache's older one). `CLNSIGCONF` is **not** transferred, so the submitter
> breakdown behind a conflicting record still cannot be inspected. The star table below is
> therefore live for `clinvar_stars`; the ≥ 2★ *gating* rules beside it remain TARGET by choice.
> The distinction the screen itself makes is `conflicting` (matched as a substring of either the
> old `conflicting_interpretations_of_pathogenicity` or the current
> `conflicting_classifications_of_pathogenicity`), which is excluded from the P/LP override.

| Stars | `CLNREVSTAT` token | Meaning |
|------:|--------------------|---------|
| 4★ | `practice_guideline` | Professional practice guideline |
| 3★ | `reviewed_by_expert_panel` | ClinGen VCEP / expert-panel call — **overrides all submitters** |
| 2★ | `criteria_provided,_multiple_submitters,_no_conflicts` | Concordant multi-lab |
| 1★ | `criteria_provided,_single_submitter` | Single lab, criteria provided |
| 1★ | `criteria_provided,_conflicting_classifications` | Conflicting (inspect `CLNSIGCONF`) |
| 0★ | `no_assertion_criteria_provided` / `no_classification_provided` | No criteria — excluded from auto-logic |

*(IMPLEMENTED — `annotations.clinvar_stars` maps exactly these tokens, canonicalised to
lowercase-alphanumeric so the comma-separated array form, the string form and ClinVar's periodic
renames all collapse to one key. The pre-2024 spellings are carried too, because a pinned older
release is a legitimate input. Stars RANK in Step 9; they never gate the screen.)*

`CLNSIG` VCF tokens use underscores: `Pathogenic`, `Likely_pathogenic`, `Pathogenic/Likely_pathogenic`, `Uncertain_significance`, `Likely_benign`, `Benign`, and `Conflicting_classifications_of_pathogenicity` (formerly `Conflicting_interpretations_of_pathogenicity`), plus low-penetrance / risk-allele and non-standard terms (`drug_response`, `association`, `protective`, `Affects`).

### Consumption rules — IMPLEMENTED

What Step 3 does today (`src/hprv/selection.py`, `src/hprv/annotations.py:clnsig_is_plp`), in order:

- **BA1 first, and it wins.** AF ≥ 0.05 → dropped, **before** ClinVar is consulted. A P/LP assertion on a common allele does **not** rescue it. This is the one place a ClinVar P/LP label is deliberately overruled, and it is the correct direction: BA1 is stand-alone benign evidence.
- **Extract P/LP** by matching `pathogenic` in `CLIN_SIG`, excluding any string containing `conflicting`, `likely_benign`, or `benign/likely`. Both spellings are matched — VEP's lowercase `&`-joined form (`pathogenic&likely_pathogenic`) and the ClinVar VCF's capitalised `/`-joined form (`Pathogenic/Likely_pathogenic`) — so the predicate survives a future switch to a real ClinVar VCF.
- **P/LP is a keep override, not a promotion.** It rescues a variant from the rarity gate (`too_common`) and from the functional ladder at **Step 3**, tagging it `hprv_keep_reason=clinvar_plp`. It does not assign a class, a tier, or ACMG weight. **Caveat — the frequency-band rescue is inert past `recessive_max`:** Step 5's inheritance modes have no ClinVar bypass and every mode gates on `rarity_af` `< recessive_max` (1e-2) or `< dominant_max` (1e-4). So a P/LP variant with `rarity_af` in `[recessive_max, benign_ba1)` (e.g. a recessive founder allele like *GJB2* c.35delG ~1%) survives Step 3 but is called by **no** Step-5 mode. That drop is a deliberate consequence of the recessive frequency ceiling, but it is now **audited** (`05_inheritance` metric `clinvar_plp_dropped_ge_recessive_max`) rather than silent; the functional-ladder bypass for genuinely rare (`< recessive_max`) P/LP variants is unaffected.
- **No star gate at the screen — by CHOICE, not by absence.** Stars are available (`clinvar_stars`, from the Step-2 ClinVar transfer) and are deliberately not consulted here: a keep/drop gate on review status would violate never-drop, so every P/LP is honored regardless, including the 0★/1★ single-submitter assertions gnomAD's own guidance would filter out. Effect: **over-retention** (curation load), not missed calls. The star level is instead a **Step-9 ranking** input, so a 1★ assertion is kept, reviewed, and ranked below a 3★ one.
- **VUS / Conflicting are never promoted** (they do not match the predicate) and are **never dropped on ClinVar grounds** — they simply face the normal rarity + functional screen, and the `clnsig` string rides along into `candidates.calls.tsv` and the IGV export. Note the difference from the TARGET: the string is *visible to a curator*, but nothing **routes** it to review, and `CLNSIGCONF` (the submitter breakdown) does not exist to inspect.
- **The release is pinned by the cache, not by us.** VEP 115's GRCh38 cache carries **ClinVar 2025-02**; the ClinVar VCF ships monthly. So provenance is recorded (via the cache version) but staleness is not controllable independently of a VEP upgrade. ClinVar reclassifies continuously and can carry outdated or single-lab calls: treat P/LP as a **triage prior, never an answer**, and plan periodic re-classification of previously reported VUS.

### The ClinVar transfer — IMPLEMENTED

`02_annotate_sites.sh` transfers review status from the version-pinned ClinVar sites VCF. This is
one of **exactly two** `bcftools annotate` transfers in the pipeline (the other is the gnomAD joint
slim), and it exists because the cache cannot supply `CLNREVSTAT` at any price. It is prefixed
`clinvar_` rather than `vep_` on purpose, so which oracle a field came from is readable at a glance.

```bash
# IMPLEMENTED (02_annotate_sites.sh): transfer review status from a version-pinned ClinVar release.
# `:=` renames on transfer. CLNREVSTAT is Number=. and its VALUES contain commas
# ("criteria_provided,_multiple_submitters,_no_conflicts"), so it arrives as an array —
# annotations._canon_revstat canonicalises every rendering to one comparable key.
bcftools annotate \
  -a "${CLINVAR_VCF}" \
  -c "INFO/clinvar_CLNREVSTAT:=INFO/CLNREVSTAT,INFO/clinvar_CLNSIG:=INFO/CLNSIG" \
  -Oz -o "${OUT_VCF}" "${IN_VCF}"
```

Two acquisition traps, both handled in `scripts/prepare_resources.sh` and both **silent** if they
are not:

- **ClinVar ships BARE contig names** (`1`, not `chr1`) while these callsets are chr-prefixed.
  `bcftools annotate` matches on the contig STRING, so an unrenamed ClinVar transfers **zero**
  records and exits 0. `prep_clinvar` renames and then proves the rename landed; Step 2
  additionally **dies** on a 0-match transfer, because a cohort union always overlaps ClinVar
  somewhere.
- **Absent ≠ 0★.** If the transfer does not run at all, `clinvar_stars` is blank and Step 9 reads
  `UNAVAILABLE` = full weight. Conflating that with 0★ would silently damp every P/LP assertion in
  a run that simply had no ClinVar resource.

### Consumption rules — TARGET (a tiering step that does not exist)

- **Auto-promote at ≥ 2★** with no conflicts; **1★ P/LP** prioritized but routed to human review.
- **Exclude 0★** from auto-logic (gnomAD guidance: keep ≥ 1★ with a specified classification).
- Inspect the `CLNSIGCONF` submitter breakdown for conflicting records — **not transferred today**;
  it would be a third field on the existing transfer plus its entry in `annotations.F`.
- Record the transferred release date in run provenance (the file is pinned; the date is not yet
  written into the audit).

## ACMG/AMP framework and ClinGen SVI refinements

> **TARGET — none of this section is implemented.** The pipeline computes no ACMG criteria, no
> ACMG points total, and no class. Step 3 emits a *keep reason* (`impact_high` | `impact_moderate`
> | `spliceai` | `cadd` | `clinvar_plp`) — screening provenance, not evidence weight — and Step 9
> emits `priority_points`, which is **not** an ACMG total (see immediately below). The section is
> retained as the specification for the tiering step, and because it is the reason the
> resource-restoration items in [limitations.md](limitations.md) are worth their cost.

### Why `priority_points` is not an ACMG total

Step 9 computes an additive, per-term composite in the *style* of Tavtigian 2020, and the
resemblance is deliberate — auditability, and honest degradation when a term is unavailable. **The
analogy stops at the arithmetic.** The total must never be read against the P ≥ 10 / LP 6–9 /
VUS 0–5 bands below, and no P/LP/VUS label may be emitted from it, because:

- the criteria are **not ACMG criteria** — a CADD-based term is not PP3, and the molecular tier is
  a screening rank, not a strength assignment;
- there is **no phenotype, segregation or functional evidence** in the pipeline at all, so PP1/BS4,
  PS3/BS3, PS4 and PM3/BP2 contribute nothing and their absence is not scored;
- the ClinVar term applies an **uncalibrated review-status damp**, which is not ACMG PP5/BP6
  (SVI removed those); and
- the artifact-penalty terms (`pts_gene_artifact`, the NHF and genotype-QC penalties) have **no
  ACMG analogue whatsoever**.

The column is named `priority_points`, never `acmg_points`, and `09_prioritize.py` prints this
caveat on every run. See [prioritization.md](prioritization.md).

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

- **PM2 → Supporting by default** (SVI Recommendation v1.0). Absence/rarity is weak evidence; do not apply PM2 at Moderate. PM2 is ACMG *evidence* and must be kept distinct from the upstream rarity **screening** gate (see [allele_frequency.md](allele_frequency.md)) — passing the rarity filter is not the same as "PM2 met." One caution specific to this contract: the screening gate reads **`rarity_af`**, which by default *is* faf95 (the quantity PM2's own guidance is written around) but on the `grpmax_proxy` opt-down is a point estimate — and on that arm the cache reports frequencies only for **dbSNP-accessioned** alleles, so "absent" is weaker evidence than it looks. PM2 built naively on the proxy arm would inherit that flaw; on the faf95 arm the gnomAD joint slim carries every allele regardless of accession.
- **PP3 / BP4 calibration (Pejaver 2022)**: continuous predictors receive strength-stratified thresholds. Use **one** predictor per variant; do not sum correlated predictors as independent evidence. The **REVEL** primary defaults are restated below and detailed in [functional_annotation.md](functional_annotation.md). AlphaMissense is *not* part of the Pejaver 2022 calibration and lacks an SVI PP3 stratification — treat it only as orthogonal support. **Both calibrated missense predictors ARE available** — REVEL and AlphaMissense are VEP plugins, required by default (`resources.vep.missense_predictors_required`) — but they feed **Step 9's missense tier**, not a scored PP3, and the ladder is a **fixed precedence** (REVEL → AlphaMissense → off-label CADD → none, reported in `missense_evidence_source`), never a max over whatever scored: best-of-N is exactly the uncalibrated stacking SVI warns against. Note the tier stops at REVEL's moderate cut (0.773 ⇒ V4); no 0.932 "strong" rung is coded, because V5 is unreachable without NMD annotation. For SPLICING, SpliceAI is wired and its `spliceai_ds_min: 0.2` IS the ClinGen SVI PP3-supporting threshold — used as a keep-path and a tier, not as scored PP3 evidence, because no ACMG tiering step exists. CADD's 25.3 cutoff is a missense-derived number applied in the screen only to non-coding variants, so it is a **discovery rank, not PP3 evidence**.
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

**AutoGVP is the intended backbone — TARGET, not wired in.** It integrates a dated ClinVar release with a modified InterVar (PVS1-strength adjustment, PP5/BP6 removed) and is built for exactly this consortium use case, which is why it is the chosen target. Adopting it needs no new ClinVar resource — the dated VCF is already transferred — but its InterVar half needs an annotation set this contract still does not carry (LOFTEE HC/LC, a ClinGen gene–disease validity table, exon/CDS position for graded PVS1). Until then the pipeline runs **no** automated classifier — which also means it inherits none of the PVS1-full-strength / PP5-BP6 over-calling above. Universal caveat that applies to the target: automated calls are **screening aids, not diagnostic** — VUS and conflicts require human review, and none replace applicable VCEP rules.

## ClinGen gene–disease validity gate (TARGET)

ClinGen classifies each gene→disease relationship as **Definitive / Strong / Moderate / Limited / Disputed / Refuted / No Known Disease Relationship**. ACMG recommends diagnostic panels include only **Definitive / Strong / Moderate** genes.

**Not implemented.** No ClinGen gene–disease validity table is ingested; no gene-level gate runs, and nothing is auto-promoted for such a gate to act on. In the target design it would gate *auto-promotion* and PVS1 applicability (≥ Moderate validity plus a known LoF mechanism) — never variant retention, consistent with the never-drop principle: curated lists act as priors, not hard filters (see [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md)). Its absence costs the screen nothing, since a gate on auto-promotion is vacuous without auto-promotion; it is a prerequisite for the tiering step, not for the first pass.

## Integration pattern (per-trio VCF)

**IMPLEMENTED — five keep reasons at the screen, a star-damped ranking term at Step 9, no
classification:**

```text
Per-trio VCF (GRCh38, GATK genotype-refined)
  └─ VEP 115 cache + CADD/SpliceAI/REVEL/AlphaMissense plugins   (functional_annotation.md)
     └─ + 2 bcftools transfers: ClinVar (CLNREVSTAT -> clinvar_stars), gnomAD joint slim (faf95)
        └─ Step 3 classify (src/hprv/selection.py) — STAR-BLIND on purpose:
           ├─ rarity_af >= 0.05 (BA1)               -> DROP (never rescued, P/LP included)
           ├─ CLIN_SIG P/LP, not conflicting        -> KEEP  'clinvar_plp'  (any star level)
           ├─ rarity_af >= 1e-2 and not P/LP        -> DROP 'too_common'
           ├─ IMPACT in {HIGH, MODERATE}            -> KEEP  'impact_high'|'impact_moderate'
           ├─ SpliceAI max delta >= 0.2             -> KEEP  'spliceai'  (checked BEFORE cadd)
           ├─ CADD_PHRED >= 25.3                    -> KEEP  'cadd'  (non-coding only, in practice)
           └─ otherwise                             -> DROP 'not_functional'
           └─ Step 9 RANK (src/hprv/prioritize.py):
              clinical term +4 (p_lp) / -4 (benign), positive limb x0.5 below 2 stars;
              blank stars = UNAVAILABLE = full weight (never 0 stars).
```

*TARGET — the ACMG tiering step this document specifies. What is now present vs still absent:*

```text
Per-trio VCF (GRCh38, GATK genotype-refined)
  └─ VEP: canonical / MANE consequence            (functional_annotation.md)      [present]
     └─ join dated ClinVar VCF: CLNSIG, CLNREVSTAT            [present]  CLNSIGCONF [absent]
        └─ join gnomAD v4.1 joint grpmax faf95     (allele_frequency.md)          [present]
           └─ REVEL / AlphaMissense (fixed precedence)  (functional_annotation.md) [present]
              └─ AutoGVP: ClinVar + modified InterVar                             [absent]
                 └─ Tavtigian points  (PVS1 graded per Abou Tayoun; PM2 → Supporting)  [absent]
                    └─ gate on ClinGen gene–disease validity ≥ Moderate            [absent]
                       └─ germline VCEP specs override generic cutoffs             [absent]
```

What remains missing is the **ACMG layer itself** — criterion assignment, graded strengths, and a
class — not the annotations under it. The evidence a tiering step would read is now largely on the
table; `priority_points` deliberately does not pretend to be that layer.

Because Kids First trios are **GATK genotype-refinement output (posterior/PP-refined), not cohort joint-genotyped**, this pipeline does **not** derive internal-cohort allele frequency. Population rarity comes from gnomAD v4.1 — by default as real **`faf95`** from the gnomAD joint sites slim, the quantity ACMG/ClinGen specify (see [allele_frequency.md](allele_frequency.md) and [limitations.md](limitations.md) §2); de-novo confidence comes from the trio PP/GQ, not cohort frequency (see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md)). Germline pediatric-cancer specifics — dominant vs recessive CPS handling, second-hit boosts, PMS2/PMS2CL — are in [pediatric_cancer.md](pediatric_cancer.md).

## Recommended defaults (this pipeline)

| Parameter | Default | Status | Notes |
|-----------|---------|--------|-------|
| ClinVar significance | VEP 115 cache `CLIN_SIG` (`--check_existing`) | **IMPLEMENTED** | Release pinned by the cache ⇒ **ClinVar 2025-02** |
| ClinVar review status | ClinVar sites VCF transfer ⇒ `clinvar_CLNREVSTAT` ⇒ `clinvar_stars` 0–4 | **IMPLEMENTED** | `resources.clinvar.vcf`, independently version-pinned. One of Step 2's two transfers; 0-match dies |
| ClinVar P/LP | Keep override: `pathogenic` matched, `conflicting`/benign excluded | **IMPLEMENTED** | Rescues from `too_common` + the functional ladder; tags `clinvar_plp`. **No class assigned** |
| Star gate at the SCREEN | — | **RETIRED BY CHOICE** | Stars exist; gating on them would violate never-drop, so every P/LP is kept and reviewed. Over-retains |
| Star damp at Step 9 | positive limb × `low_star_scale` (**0.5**) below `min_review_stars` (**2**) | **IMPLEMENTED** | Positive limb only — damping a low-star *benign* call would promote it. Blank = `UNAVAILABLE` = full weight ≠ 0★ |
| BA1 vs P/LP | AF ≥ 0.05 drops **even if P/LP** | **IMPLEMENTED** | BA1 is evaluated first and is never rescued |
| Conflicting / VUS | Not promoted, not dropped on ClinVar grounds; `clnsig` emitted to output | **IMPLEMENTED (partial)** | Never-drop holds; but no review **routing**, and no `CLNSIGCONF` to inspect |
| ClinVar release | Dated, version-pinned; recorded in provenance | *TARGET* | Needs a ClinVar VCF. Re-classify prior VUS periodically |
| Auto-promote P/LP | `CLNSIG` P or LP **and ≥ 2★**, no conflicts | *TARGET* | 1★ P/LP → prioritize + human review |
| 0★ records | Excluded from auto-logic | *TARGET* | gnomAD guidance: keep ≥ 1★ with a specified classification |
| ACMG classifier | **AutoGVP** (ClinVar + modified InterVar, PP5/BP6 removed) | *TARGET* | Graded PVS1. Nothing wired in |
| Combining rule | Tavtigian points: **P ≥ 10, LP 6–9, VUS 0–5, LB −1…−5, B ≤ −6** | *TARGET* | Supporting ±1 / Moderate ±2 / Strong ±4 / Very Strong ±8 |
| PM2 | **Supporting** strength only | *TARGET* | ACMG evidence, distinct from the rarity screening gate (which exists; PM2 does not) |
| PVS1 | Graded (Abou Tayoun 2018); gene validity ≥ Moderate + known LoF | *TARGET* | Blocked twice: no LOFTEE, no validity table |
| PP3/BP4 missense (as ACMG evidence) | **REVEL** (single predictor): PP3 supporting ≥ 0.644, moderate ≥ 0.773, strong ≥ 0.932; BP4 supporting ≤ 0.290, moderate ≤ 0.183 | *TARGET* | Ranges from Pejaver 2022. The predictors ARE wired (Step-9 tier, fixed precedence, no 0.932 rung); what is absent is the ACMG *criterion assignment* |
| Gene gate | ClinGen **Definitive / Strong / Moderate** validity | *TARGET* | VCEP threshold overrides generic cutoff |

Values marked IMPLEMENTED are **configurable defaults** (`config/config.example.yaml`), not immutable law; values marked *TARGET* are specifications with no code behind them yet. A gene-specific ClinGen VCEP specification (its own PVS1 strength, PM2/BA1/BS1 thresholds, or PP3 calibration) **overrides** any generic default here — a rule the target tiering step must honor, and which nothing in the first pass currently implements. If this table ever disagrees with [Canonical defaults](README.md#canonical-defaults), that table wins.

## Scope limitations (stated honestly)

The contract-level gaps are catalogued once in **[limitations.md](limitations.md)** — read it before
interpreting a negative result. The two that bear directly on this document:

- **Stars RANK rather than gate** ([limitations.md](limitations.md) §6): the ClinVar VCF is
  transferred in Step 2 ⇒ `clinvar_stars` (0–4), which scales the Step-9 clinical term's positive
  limb only. The SCREEN remains star-blind on purpose — gating it would violate never-drop. So
  **a `clinvar_plp` keep is still not evidence of curated confidence**: it may be one submitter
  with no criteria. What changed is that the reviewer can now *see* which, in the `clinvar_stars`
  column, instead of checking ClinVar's web UI by hand. The direction of error at the screen
  remains over-retention.
- **The cached `CLIN_SIG` is stale by construction.** The P/LP *override* reads the cache
  (VEP 115 ⇒ **2025-02**), not the transferred release, so a variant reclassified since is
  overridden at its old assertion in both directions: a since-downgraded P/LP still rescues, and a
  since-upgraded VUS gets no rescue. The transferred `clinvar_CLNSIG` rides along precisely so the
  disagreement is visible — comparing them is a curator step, not an automated one.

And the pre-existing scope limits:

- **SNV/indel only initially.** The ClinVar override here acts on SNVs and indels (no ACMG classification is computed for any variant class). CNV/SV pathogenic calls (e.g. single-exon *RB1*/*SMARCB1*/*DICER1*/*NF1* deletions, *PMS2* rearrangements — 10–15% of pediatric-cancer and rare-disease diagnoses) are a known blind spot pending a future GATK-gCNV / Manta / ExomeDepth module.
- **Pseudogene / segmental-duplication genes** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) are low-confidence from short reads; ClinVar/ACMG calls in those regions must be treated as suspect regardless of star level — a 3★ assertion in a paralogous region is a statement about the variant, not about your ability to call it here. Note the pipeline does **not** currently flag or mask these regions — [limitations.md](limitations.md), structural gaps.
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
