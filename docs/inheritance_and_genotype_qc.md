# Inheritance Models, Genotype Refinement & QC

How this pipeline uses GATK genotype-refinement annotations, per-genotype QC gates, and Mendelian inheritance logic to nominate high-priority rare **inherited germline** variants from per-trio GRCh38 VCFs. The pipeline's focus is inherited germline variation; de novo filtering/review and mtDNA heteroplasmy are handled by separate dedicated pipelines.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- Inputs are **per-trio** (mother–father–child) VCFs on GRCh38 from the GATK Genotype-Refinement workflow — **not** jointly genotyped across the cohort, so population frequency comes from external gnomAD v4.1, never from callset AC/AN (see [allele_frequency.md](allele_frequency.md)).
- **Focus is inherited germline variation.** The four first-class inheritance modes are: **dominant** (rare inherited het), **autosomal recessive homozygous**, **compound het (in trans)**, and **X-linked recessive**. **De novo** is a *secondary* cross-reference only (its dedicated filtering/review lives in separate bespoke machinery). **mtDNA heteroplasmy** is out of scope here — handled by a separate dedicated pipeline.
- **Dominant (new, first-class):** a rare (grpmax **faf95 < 1e-4**), functional, **inherited heterozygous** variant transmitted from ≥ 1 parent (parent-of-origin recorded: maternal / paternal / both) and **not** part of a compound-het pair. This is the key new signal — heterozygous variants become interesting when they **recur across multiple individuals in the same gene** (see gene consolidation, [gene_burden.md](gene_burden.md)).
- Trust the **refined PP-derived GQ**, not raw PL/GQ. Default per-genotype gates: **GQ ≥ 20**, **DP ≥ 10** (DP ≥ 20 for de novo), het **AB 0.25–0.75**, hom-alt **AB ≥ 0.90**, hom-ref **AB ≤ 0.10**, **FILTER = PASS** only. AB is derived from AD (it is not a native FORMAT field).
- **Compound het** = two rare hets in the *same gene* in **trans**, established by parent-of-origin (mat-only + pat-only) or WhatsHap read-backed phasing; a de novo second hit is a valid partner.
- **De novo (secondary cross-reference):** detected via GATK **`hiConfDeNovo`** (child-membership checked via `annotations.is_hiconf_denovo_for`), then **re-verified** with DP/AB plus **parental cleanliness** (each parent alt AD ≤ 1, DP ≥ 10) and gnomAD absent/singleton. It is **not** the driver here — dedicated bespoke machinery handles de novo filtering and review.
- **X/hemizygous:** apply sex-aware ploidy, **drop male non-PAR chrX/chrY het calls**, separate PAR from non-PAR. **mtDNA heteroplasmy** is handled by a dedicated pipeline (out of scope).
- Sample/pedigree QC via **Peddy** (parent–child IBS0 ≈ 0, relatedness ≈ 0.5), genome-wide **Mendelian-error rate < 2%**, and **UPDhmm** per-chromosome UPD screening.
- **Known failure mode:** the gnomAD prior in `CalculateGenotypePosteriors` can push a genuine ultra-rare pathogenic call toward hom-ref — cross-check pre-refinement PL/GT for top candidates.

---

## 1. What the Genotype-Refinement workflow gives us

The GMKF Kids First inputs are the output of the GATK **Genotype-Refinement workflow**, which post-processes trio genotypes to *sharpen genotype accuracy* — it is not variant discovery. It runs three tools in sequence and adds the annotations this pipeline keys on.

### 1.1 CalculateGenotypePosteriors (PL → posterior)

Recomputes each genotype probability by combining the HaplotypeCaller likelihoods (PL) with two priors:

- **Population prior** from a `--supporting-callsets` VCF (the standard Kids First choice is `af-only-gnomad.hg38.vcf.gz`). At matching sites the tool derives priors from the supporting AC; away from those sites it would use AC discovered among input samples unless `--discovered-ac-priors-off`. **Discovered-AC priors are only meaningful with ≥10 called samples** — a single per-trio VCF has 3 samples and is not jointly genotyped, so this condition is not met. Population priors therefore come almost entirely from gnomAD plus the pedigree, **not** from cohort AC.
- **Family/Mendelian prior** from the trio `-ped` pedigree.

Outputs (all Phred-scaled): **PP** (posterior per genotype; supersedes PL for the called GT), **GQ overwritten** from PP, **PG** (prior vector HomRef/Het/HomVar), and for trios **JL**/**JP** (joint likelihood/posterior that the trio genotypes are wrong). The called GT can flip when a genotype is judged highly likely wrong.

> **`--num-reference-samples-if-no-call` is resource-dependent.** Set it to the documented sample size of the *specific* supporting resource you ship, and treat the exact N as version-dependent. Note that the commonly-shipped `af-only-gnomad.hg38.vcf.gz` is a gnomAD **exome** AF resource (N on the order of 10^5–10^6 depending on release), so a "76,156 genomes" value is likely wrong for it. Do not leave the default at 0 for a 3-sample callset.

**Trust model:** after this step, **PP / refined GQ is the field to trust**, not the original PL/GQ.

### 1.2 VariantFiltration — low-GQ genotype flag

Flags genotypes with **refined GQ < 20** (`--genotype-filter-expression "GQ < 20" --genotype-filter-name lowGQ`) at the FORMAT/genotype level. GQ 20 ≈ 99% confidence in the genotype. `PossibleDeNovo` ignores genotypes already tagged filtered, so run filtration **before** annotation when you want `lowGQ` genotypes excluded from de novo consideration (the usual choice).

### 1.3 VariantAnnotator PossibleDeNovo — hiConf/loConf tags

Adds INFO `hiConfDeNovo` / `loConfDeNovo` from the pedigree. Exact criteria (biallelic sites, requires a Mendelian violation with parents ref/ref → child het/hom-alt):

| Tag | Child GQ | Parent GQ | AC/AF cutoff |
|-----|----------|-----------|--------------|
| `hiConfDeNovo` | ≥ 20 | both ≥ 20 | AC below `max(4, 0.001 × N_samples)` |
| `loConfDeNovo` | ≥ 10 | both > 0 | same cutoff |

The tool applies a de novo prior of **1e-6** per site. It is a **genotype-configuration + GQ classifier only** — no allele-balance, depth (default `depthThreshold=0`), or PL-ratio filter is applied. With a single trio the AC/AF cutoff (≥ 4 samples) is essentially inert, so classification reduces to the GQ thresholds. Treat `PossibleDeNovo` as a **screen, not a validated de novo caller**.

### 1.4 Site filtering upstream of genotypes

A per-trio callset cannot support VQSR (which needs cohort scale). Verify the FILTER column and build only on **PASS** sites; if hard-filtered, the GATK reference thresholds are:

| Class | Hard-filter reference thresholds |
|-------|----------------------------------|
| SNP | `QD < 2`, `FS > 60`, `MQ < 40`, `MQRankSum < -12.5`, `ReadPosRankSum < -8`, `SOR > 3` |
| INDEL | `QD < 2`, `FS > 200`, `ReadPosRankSum < -20` |

Deep-learning site filtering has moved from the deprecated CNNScoreVariants to **NVScoreVariants** (PyTorch) in recent GATK.

---

## 2. Cross-cutting genotype & sample QC (apply before mode logic)

These gates apply to every inheritance mode. AB (allele balance) is computed from **AD** as `alt / (ref + alt)`; it is not a native GATK FORMAT field.

| Gate | Default | Notes |
|------|---------|-------|
| FILTER | `PASS` only | site-level |
| Genotype quality | **GQ ≥ 20** (from PP) | flag `lowGQ` at GQ < 20 |
| Depth | **DP ≥ 10** per sample | **DP ≥ 20** preferred for de novo |
| Het allele balance | **0.25 ≤ AB ≤ 0.75** | common practical band 0.2–0.8 |
| Hom-alt allele balance | **AB ≥ 0.90** | |
| Hom-ref allele balance | **AB ≤ 0.10** | |

**Sample / pedigree QC:**

- **Peddy** confirms reported vs. genotype-inferred sex (chrX heterozygosity) and pedigree relatedness via IBS0 / Rel; parent–child pairs should show **IBS0 ≈ 0** and **relatedness ≈ 0.5**. Elevated IBS0 signals a sample swap or mislabeled trio → block downstream inheritance logic. Kids First runs Peddy explicitly.
- **Genome-wide Mendelian-error rate < 2%** is the concordance backstop: higher signals a bad trio/swap; localized MIE clusters flag CNV/UPD.
- **Contamination gate.** kid/dad/mom roles are already peddy-verified upstream, so Step 0's job is to catch the *less-curated* trios. It flags a trio if any member is contaminated: **verifyBamID2 `FREEMIX` > 0.05** when a directory of `*.selfSM` files is supplied (`resources.selfsm_dir`; this mirrors the group's DNM freemix QC), otherwise a **VCF-only CHARR** proxy — the reference-allele read fraction at high-quality homozygous-ALT SNV sites, thresholded at **> 0.02**. An uncontaminated sample sits near 0 (sequencing error + reference bias); cross-sample contamination injects reference reads at hom-alt sites and raises it. This matters because 1–3% contamination turns hom-ref → apparent-het, manufacturing false inherited hets / comp-het second hits. (CHARR: Lu et al., *Am J Hum Genet* 2023, [PMC10716339](https://pmc.ncbi.nlm.nih.gov/articles/PMC10716339/).)

> **Known failure mode (gnomAD-prior suppression).** Because `CalculateGenotypePosteriors` folds a population prior into every genotype, a genuinely rare/private pathogenic variant can be *pushed toward hom-ref* by a low-AF prior, suppressing a real finding. For candidate high-priority variants, **cross-check the pre-refinement PL/GT** to ensure refinement did not down-weight a true rare call.

---

## 3. Inheritance models

Genotype-level rules assume a mother–father–child trio VCF on GRCh38 with GT, GQ, DP, AD, and PL/PP available. Rarity gates below are stated as grpmax `faf95` (the pipeline's filter field); see [allele_frequency.md](allele_frequency.md) for the frequency oracle and [inheritance engine tooling](#5-inheritance-engines--tooling).

The four first-class inheritance modes are **dominant (inherited het)**, **autosomal recessive homozygous**, **compound het (in trans)**, and **X-linked recessive**. **De novo** is retained only as a lightweight *secondary* cross-reference ([§3.5](#35-de-novo-secondary-cross-reference)); **mtDNA heteroplasmy** is out of scope ([§3.6](#36-mitochondrial-chrm--out-of-scope)).

### 3.1 Dominant — rare inherited heterozygous

- **Genotypes:** child `0/1`, with the alt allele **transmitted from ≥ 1 parent** — mother `0/1` (maternal), father `0/1` (paternal), or both parents `0/1` (recorded as parent-of-origin: maternal / paternal / both). A child het with both parents `0/0` is de novo, not this mode (see [§3.5](#35-de-novo-secondary-cross-reference)).
- **Requirement:** the variant must be **functional** and **not** part of a compound-het pair in the same gene ([§3.3](#33-compound-heterozygous-two-hets-in-trans)) — a single inherited het, not one half of a biallelic hit.
- **QC:** child confident het (GQ ≥ 20, DP ≥ 10, AB 0.25–0.75) *and* the transmitting parent(s) confident het on the same criteria.
- **Rarity:** external gnomAD grpmax **faf95 < 1e-4** (dominant rarity gate; same as de novo).
- **Why it matters:** heterozygous inherited variants are individually low-specificity, but become compelling when they **recur across multiple distinct individuals in the same gene**. Gene consolidation ([§5](#5-inheritance-engines--tooling), [gene_burden.md](gene_burden.md)) tallies dominant-het carriers per gene and ranks recurrent genes first, **weighted by gene constraint** (a recurrent het in a haploinsufficient gene is the most compelling). Do **not** reject solely because an unaffected parent carries it — reduced/age-dependent penetrance is expected ([§6](#6-co-segregation--penetrance-modifiers-never-hard-filters)).

### 3.2 Autosomal recessive — homozygous

- **Genotypes:** child `1/1`, mother `0/1`, father `0/1`.
- **QC:** both parents confident hets (GQ ≥ 20, AB 0.25–0.75), child confident hom-alt (AB ≥ 0.90, DP ≥ 10).
- A `1/1` child with a `0/0` parent is a Mendelian error → suspect a hemizygous "false hom" (deletion on the other allele) or **UPD** (see [§4](#4-upd-imprinting--mosaicism-flags-not-primary-calls)).
- **Rarity:** per-allele grpmax **faf95 < 1e-2** (permissive discovery default), with a stricter **1e-3** high-confidence tier. Applied per variant, not per gene. **Do not** down-weight biallelic candidates using pLoF constraint ([gene_constraint.md](gene_constraint.md)).

### 3.3 Compound heterozygous (two hets in trans)

- **Requirement:** two rare het variants in the same gene on **opposite** haplotypes (trans). Cis pairs are non-causal.
- **Trio phasing (preferred, "free"):** variant A from mother (mat `0/1`, pat `0/0`) and variant B from father (pat `0/1`, mat `0/0`) ⇒ trans by descent. Both variants tracing to the same parent ⇒ cis ⇒ reject. A **de novo second hit** is a legitimate partner (one inherited + one de novo).
- **Read-backed phasing (WhatsHap):** resolves phase directly when both variants lie within one read/fragment, and uniquely combines read-based *and* pedigree phasing. Caveat: WhatsHap drops variants with missing or Mendelian-inconsistent parental genotypes, lowering the phasing rate when parental data is incomplete. For exon-spanning pairs too far apart for short reads, trio phasing is the fallback.
- **QC:** both variants must independently pass het QC (GQ ≥ 20, DP ≥ 10, AB 0.25–0.75) in the child *and* the transmitting parent.
- **Rarity:** per-allele grpmax **faf95 < 1e-2** (1e-3 high-confidence).

### 3.4 X-linked / hemizygous

- **X-recessive (affected male):** child hemizygous alt (`1/1` from a diploid caller, or `1` if ploidy-aware), mother `0/1` carrier, father `0/0` (in region). Affected female = `1/1` with both parents carriers.
- **X-dominant:** child het/hemi affected (de novo or transmitted); watch for male-lethal patterns.
- **Hemizygous chrX / chrY (male, outside PAR):** males are haploid, so a diploid caller emits `0/0` or `1/1` — **a het call in a male non-PAR region is a QC red flag** (mapping artifact, PAR misplacement, or XXY). Enforce **sex-aware ploidy** and **drop** male non-PAR chrX/chrY het calls. Distinguish PAR1/PAR2 (diploid) from non-PAR. QC the hemizygous alt call at DP ≥ 10, GQ ≥ 20, AB ≥ 0.90. Peddy sex-check must pass first.

### 3.5 De novo (secondary cross-reference)

De novo is **not** the driver of this pipeline: dedicated bespoke machinery handles de novo filtering **and** review. Here it is retained only as a lightweight secondary cross-reference against the inherited-variation signal, and counted **separately** from the inherited modes in gene consolidation ([§5](#5-inheritance-engines--tooling)).

- **Genotypes:** child `0/1`, mother `0/0`, father `0/0` (autosomal); hemizygous de novo on male chrX is `0 → 1` (see [X-linked](#34-x-linked--hemizygous)).
- **Detection:** GATK `hiConfDeNovo` present (all three trio-member GQ ≥ 20), with **child membership** verified via `annotations.is_hiconf_denovo_for` (confirming the tag applies to *this* child). Use `loConfDeNovo` (child GQ ≥ 10) only as a lower-sensitivity tier.
- **Re-verify (the tool does not):** child **DP ≥ 20**, het **AB 0.25–0.75**, all three GQ ≥ 20, and **parental cleanliness** — each parent alt AD ≤ 1 with DP ≥ 10 (a parental alt fraction of a few percent suggests inherited or parental mosaicism, not de novo).
- **Rarity:** external gnomAD grpmax **faf95 < 1e-4** and low `nhomalt` (≤ 1, enforced). *(An AC-based **absent-or-singleton** test is a target refinement, not yet enforced — de novo is secondary here.)*
- **ACMG evidence:** PS2 (confirmed de novo) / PM6 (assumed de novo) are scored by the **ClinGen SVI point system**. See [clinical_classification.md](clinical_classification.md).
- **Enrichment (optional secondary):** de novo Poisson enrichment vs the Samocha mutation model (denovolyzeR-style; exome-wide P < 2.5e-6, BH q < 0.05) is reported only when a mutation-rate table is supplied ([gene_burden.md](gene_burden.md)).

### 3.6 Mitochondrial (chrM) — out of scope

- **mtDNA heteroplasmy is handled by a separate dedicated pipeline** and is out of scope here; chrM is not among the active inheritance modes in this pipeline.

---

## 4. UPD, imprinting & mosaicism (flags, not primary calls)

- **Uniparental disomy (UPD):** suspect when a chromosome shows excess Mendelian errors from one parent or long isodisomic ROH; it can unmask a "homozygous recessive" from a single carrier parent or cause imprinting disorders. Screen per-chromosome with **UPDhmm** (HMM on the trio VCF; avoids ROH confounding by consanguinity) or UPDio.
- **Parental mosaicism:** a "de novo" with low-level parental alt reads (AB ~1–10%, below the het cutoff) is transmitted from a mosaic parent → alters recurrence risk. Detect by inspecting parental AB on de novo candidates rather than hard `0/0` gating.
- **Proband (post-zygotic) mosaicism** — *known scope limitation.* Low-VAF somatic/post-zygotic calls in the proband (relevant for NF1, overgrowth, and some cancer-predisposition phenotypes) fall **outside** the het AB 0.25–0.75 band and will be filtered. A dedicated mosaic tier is future work.

---

## 5. Inheritance engines & tooling

| Tool | Role | Status / notes |
|------|------|----------------|
| **GATK** CalculateGenotypePosteriors / VariantFiltration / PossibleDeNovo | Source of PP/GQ and `hiConfDeNovo`/`loConfDeNovo` | Already upstream; reuse these annotations — provenance-clean. GATK 4.6.2.0 (2025-04-13). |
| **slivar** (v0.3.4) | Core segregation engine | JS-expression filtering with built-in dominant/inherited-het, recessive/hom-alt, `comphet`, `x_denovo`, and `denovo` helpers driven by a PED file. Parent-of-origin for inherited hets comes from the trio genotypes. Single static binary → trivially containerizable. Pinned in [tooling_and_reproducibility.md](tooling_and_reproducibility.md). |
| **WhatsHap** | Read-based + pedigree phasing | For compound-het trans resolution when reads span both variants. |
| **Peddy** | Sex + relatedness + Mendelian-error QC | IBS0/Rel checks; run before inheritance logic. Upstream source of the kid/dad/mom role assignments. |
| **verifyBamID2 / CHARR** | Per-sample contamination | Step 0 gate: ingest verifyBamID `FREEMIX` (`*.selfSM`) if available, else a VCF-only CHARR estimate from AD at hom-alt sites. |
| **UPDhmm** | Per-chromosome UPD detection | HMM on the trio VCF. |
| **cyvcf2 / pysam / bcftools** | Custom QC gates & edge cases | Hemizygous handling, parental-mosaicism thresholds, de novo QC (AB/DP), gene-list joins — trivial to tune where slivar's canned expressions don't expose the knob. |
| **TrioDeNovo** (v0.06) | Optional orthogonal DNM caller | Stale (C++, ~2015); use only for a second-caller consensus on high-stakes DNMs. |
| **GEMINI** | — | **Deprecated** (Python 2.7 + SQLite); the author points users to slivar. Do not adopt. |

Per-mode calls feed **recurrence-based gene consolidation** (Step 6, [gene_burden.md](gene_burden.md)): tally the number of **distinct individuals** per gene carrying a qualifying variant under each inherited model (**dominant het** / **biallelic** [hom + comp-het] / **X-linked**), with **de novo counted separately** as a secondary signal. A gene is **recurrent** at ≥ **min_carriers** (default **2**) distinct individuals; rank recurrent-first, **weighted by gene constraint** (LOEUF / pLI / s_het — a recurrent het in a haploinsufficient gene is the most compelling; see [gene_constraint.md](gene_constraint.md)). De novo Poisson enrichment vs the Samocha model (denovolyzeR-style; exome-wide P < 2.5e-6, BH q < 0.05) is an optional secondary signal, reported only when a mutation-rate table is supplied; TRAPD vs matched gnomAD is an optional corroboration (not yet implemented).

Phenotype-driven prioritization (Exomiser / LIRICAL) and curated gene lists are covered in [gene_lists_and_phenotype.md](gene_lists_and_phenotype.md); they act as ranking priors under the never-drop rule, never as hard reporting gates.

### 5.1 Normalize before any AF match or segregation

```bash
# Left-align and split multiallelics FIRST, against the GRCh38 reference,
# before any allele-frequency join or segregation logic.
bcftools norm -f GRCh38.fa -m -both input.vcf.gz -Oz -o normalized.vcf.gz
bcftools index -t normalized.vcf.gz
```

### 5.2 Genotype-QC gating (illustrative, parameterized)

```bash
# Keep PASS sites, then flag low-confidence genotypes on refined GQ.
# (GQ here is the PP-derived, refined GQ after CalculateGenotypePosteriors.)
bcftools view -f PASS normalized.vcf.gz \
  | bcftools +setGT -- -t q -n . -i 'FMT/GQ<20 | FMT/DP<10' \
  > gt_filtered.vcf
```

Allele balance is derived from `FMT/AD` in downstream cyvcf2/slivar logic (e.g. het band `0.25 ≤ AD[alt]/sum(AD) ≤ 0.75`), since AB is not a native FORMAT field.

---

## 6. Co-segregation & penetrance (modifiers, never hard filters)

Do **not** reject an inherited candidate solely because an unaffected parent carries it — reduced or age-dependent penetrance is expected in dominant rare-disease and pediatric-cancer genes. Conversely, a phenocopy (variant not tracking with phenotype) weakens co-segregation. Treat co-segregation as ACMG/ClinGen **supporting** evidence (PP1 / BS4), not a hard filter. For germline pediatric cancer, retain heterozygous LP/P in known predisposition genes regardless of parental affected status (see [pediatric_cancer.md](pediatric_cancer.md)).

---

## 7. Known scope limitations

- **SNV/indel only initially.** CNV/SV are a real blind spot — ~10–15% of pediatric-cancer and rare-disease diagnoses are CNV/SV (single-exon RB1/SMARCB1/DICER1/NF1 deletions, PMS2 rearrangements). A het `1/1` child with a `0/0` parent, or a "hemizygous false hom," can be the *shadow* of an undetected deletion. Future module: GATK-gCNV / Manta / ExomeDepth.
- **Pseudogene / segmental-duplication regions** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) are low-confidence from short reads — flag those regions rather than trust naïve genotypes.
- **Proband post-zygotic mosaicism** falls outside the het AB band (see [§4](#4-upd-imprinting--mosaicism-flags-not-primary-calls)).
- **Calibration/validation** of the de novo and segregation logic should be measured against GIAB/CMRG truth sets and a positive-control variant panel; this is not yet wired in.
- **Phenotype (HPO) dependency:** the phenotype-scoring layer needs per-proband HPO terms, which are not part of genotype QC but gate downstream prioritization; handle sparse/absent phenotype gracefully.

---

## Recommended defaults (this pipeline)

| Setting | Default | Override / note |
|---------|---------|-----------------|
| Genotype quality | **GQ ≥ 20** (refined PP-derived) | flag `lowGQ` at GQ < 20 |
| Depth | **DP ≥ 10** per sample | **DP ≥ 20** for de novo |
| Het allele balance | **AB 0.25–0.75** | from AD |
| Hom-alt allele balance | **AB ≥ 0.90** | from AD |
| Hom-ref allele balance | **AB ≤ 0.10** | from AD |
| Site filter | **FILTER = PASS** | hard-filter thresholds in [§1.4](#14-site-filtering-upstream-of-genotypes) |
| Dominant (inherited het) | rare functional **het** transmitted from ≥ 1 parent; parent-of-origin recorded (mat / pat / both); **not** a comp-het partner | grpmax **faf95 < 1e-4**; recurrence-consolidated |
| Recessive / comp-het rarity | grpmax **faf95 < 1e-2** per allele | **1e-3** high-confidence tier |
| Recessive hom | child `1/1` + both parents `0/1`, all GQ ≥ 20, DP ≥ 10 | 1/1 vs 0/0 parent → suspect deletion/UPD |
| Compound het | two rare hets, same gene, **trans** (parent-of-origin or WhatsHap) | de novo second hit is valid |
| X-linked recessive | male hemizygous + carrier mother; sex-aware ploidy; **drop** male non-PAR chrX/chrY het calls; separate PAR/non-PAR | |
| De novo (secondary) | GATK **`hiConfDeNovo`** (child-membership via `is_hiconf_denovo_for`); re-verify child DP ≥ 20, AB 0.25–0.75, each parent alt AD ≤ 1 / DP ≥ 10, gnomAD absent/singleton; grpmax **faf95 < 1e-4** | cross-reference only; filtering/review in separate machinery. Optional: Poisson enrichment vs Samocha (P < 2.5e-6, BH q < 0.05) when a mutation-rate table is supplied |
| Gene consolidation | recurrent at ≥ **min_carriers** (default **2**) distinct individuals; rank recurrent-first, weighted by constraint (LOEUF / pLI / s_het) | dominant / biallelic / X-linked; de novo counted separately |
| mtDNA heteroplasmy | out of scope — handled by a separate dedicated pipeline | chrM not an active mode here |
| Sample/pedigree QC | **Peddy** (parent–child IBS0 ≈ 0, rel ≈ 0.5); genome-wide MIE **< 2%**; **UPDhmm** per chromosome | |
| `--num-reference-samples-if-no-call` | set to the shipped supporting-resource's documented N | resource-dependent; version-pin |
| Rare-variant safeguard | cross-check pre-refinement PL/GT for top candidates | gnomAD-prior suppression |

All values are **configurable defaults** in `config/config.example.yaml`, not immutable law. A gene-specific ClinGen VCEP threshold overrides any generic cutoff.

---

## Sources

- GATK Genotype-Refinement workflow (CGP → VariantFiltration → PossibleDeNovo; 1e-6 de novo prior; lo/hi-conf GQ 10/20): https://gatk.broadinstitute.org/hc/en-us/articles/360035531432-Genotype-Refinement-workflow-for-germline-short-variants ; archived mirror: https://github.com/broadgsa/gatk/blob/master/doc_archive/methods/Genotype_Refinement_workflow.md
- GATK CalculateGenotypePosteriors (supporting-callsets, discovered-AC, ≥10 samples): https://gatk.broadinstitute.org/hc/en-us/articles/360037226592-CalculateGenotypePosteriors
- CalculateGenotypePosteriors gnomAD command example (NIH Biowulf): https://hpc.nih.gov/training/gatk_tutorial/gt-posteriors.html
- PossibleDeNovo source (exact hiConf/loConf constants): https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/PossibleDeNovo.java ; doc: https://gatk.broadinstitute.org/hc/en-us/articles/13832750995739-PossibleDeNovo
- Hard-filtering germline short variants: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
- Filter with VQSR or hard-filtering: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
- NVScoreVariants (replaces CNNScoreVariants): https://gatk.broadinstitute.org/hc/en-us/articles/10064202674971-Introducing-NVIDIA-s-NVScoreVariants-a-new-deep-learning-tool-for-filtering-variants
- GATK 4.6.2.0 release (2025-04-13): https://github.com/broadinstitute/gatk/discussions/9148
- Kids First genotype-refinement workflow (Peddy, VEP, per-trio): https://github.com/kids-first/kf-genotype-refinement-workflow
- De novo consensus / force-call precision, parental-AB filtering (Life Science Alliance 2025): https://www.life-science-alliance.org/content/8/6/e202403039
- De novo filter thresholds vs coverage (Heredity 2025): https://www.nature.com/articles/s41437-025-00754-0
- ClinGen SVI PS2/PM6 recommendation v1.1: https://www.clinicalgenome.org/site/assets/files/3461/svi_proposal_for_de_novo_criteria_v1_1.pdf
- slivar releases & rare-disease wiki: https://github.com/brentp/slivar/releases · https://github.com/brentp/slivar/wiki/rare-disease
- GEMINI (deprecated): https://github.com/arq5x/gemini
- TrioDeNovo: https://genome.sph.umich.edu/wiki/Triodenovo · https://pmc.ncbi.nlm.nih.gov/articles/PMC4410659/
- WhatsHap read-based + pedigree phasing: https://www.biorxiv.org/content/10.1101/085050v2.full
- SmartPhase (WhatsHap parental-genotype caveat; trans/cis for diagnosis): https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007613
- Peddy (sex/relatedness/IBS0/Mendelian error, AJHG 2017): https://github.com/brentp/peddy · https://www.cell.com/ajhg/pdf/S0002-9297(17)30017-4.pdf
- mtDNA caller benchmark, Mutserve/Haplocheck, heteroplasmy thresholds (Frontiers Genetics 2022): https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.692257/full · Mutserve: https://github.com/seppinho/mutserve
- gnomAD mtDNA (Laricchia et al., *Genome Research* 2022;32(3):569–582; gnomAD v3.1, 56,434 individuals): https://genome.cshlp.org/content/32/3/569 · https://pmc.ncbi.nlm.nih.gov/articles/PMC8896463/
- UPDhmm (trio UPD detection, Bioinformatics 2026): https://academic.oup.com/bioinformatics/article/42/3/btag062/8529595
- UPD MIE values / trio genotype method (King et al., Genome Research 2014): https://genome.cshlp.org/content/24/4/673.full.html
- ACMG/AMP variant classification (Richards et al. 2015) and ACGS 2024 UK practice guidelines: https://www.genomicseducation.hee.nhs.uk/wp-content/uploads/2024/08/ACGS-2024_UK-practice-guidelines-for-variant-classification.pdf
