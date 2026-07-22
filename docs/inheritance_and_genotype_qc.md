# Inheritance Models, Genotype Refinement & QC

How this pipeline uses GATK genotype-refinement annotations, per-genotype QC gates, and Mendelian inheritance logic to nominate high-priority rare **inherited germline** variants from per-trio GRCh38 VCFs. The pipeline's focus is inherited germline variation; de novo filtering/review and mtDNA heteroplasmy are handled by separate dedicated pipelines.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

> ### ⚠ The rarity field in this document is **not** `faf95`
>
> The pipeline runs a **VEP-only contract**: annotations come from a VEP 115 GRCh38 cache plus the
> CADD plugin, and nothing else. The cache carries **no AC/AN**, so `faf95` (a 95% CI lower bound)
> **cannot be computed at any price** — it is a **TARGET**, not what runs. Every rarity gate below
> is applied to the **grpmax proxy**: the max gnomAD v4.1 point-estimate AF over the
> grpmax-eligible groups (`AFR/AMR/EAS/NFE/SAS`), `annotations.frequency()`. The **numbers** are
> unchanged; the **field** is a point estimate, so the gates run slightly stringent on low-count
> alleles (a false-negative direction). Read `faf95` below as "the target field", `grpmax proxy` as
> "what Step 5 actually reads". Never substitute VEP's `MAX_AF` (it maxes over bottlenecked founder
> groups ⇒ drops real candidates) or a global AF (dilutes ⇒ over-retains).
> Also **not available**: `nhomalt`, so the de novo homozygote condition is **retired**.
> Full ledger: **[limitations.md](limitations.md)** (§2 faf95, §2a the MAX_AF trap, §3 nhomalt).

## TL;DR

- Inputs are **per-trio** (mother–father–child) VCFs on GRCh38 from the GATK Genotype-Refinement workflow — **not** jointly genotyped across the cohort, so population frequency comes from external gnomAD v4.1, never from callset AC/AN (see [allele_frequency.md](allele_frequency.md)).
- **Focus is inherited germline variation.** The four first-class inheritance modes are: **dominant** (rare inherited het), **autosomal recessive homozygous**, **compound het (in trans)**, and **X-linked recessive**. **De novo** is a *secondary* cross-reference only (its dedicated filtering/review lives in separate bespoke machinery). **mtDNA heteroplasmy** is out of scope here — handled by a separate dedicated pipeline.
- **Dominant (new, first-class):** a rare (grpmax **proxy AF < 1e-4** — not `faf95`, see the banner above), functional, **inherited heterozygous** variant transmitted from ≥ 1 parent (parent-of-origin recorded: maternal / paternal / both) and **not** part of a compound-het pair. This is the key new signal — heterozygous variants become interesting when they **recur across multiple individuals in the same gene** (see gene consolidation, [gene_burden.md](gene_burden.md)).
- Trust the **refined PP-derived GQ**, not raw PL/GQ. Default per-genotype gates: **GQ ≥ 20**, **DP ≥ 10** (DP ≥ 20 for de novo), het **AB 0.25–0.75**, hom-alt **AB ≥ 0.90**, hom-ref **AB ≤ 0.10**, **FILTER = PASS** only. AB is derived from AD (it is not a native FORMAT field).
- **Compound het** = two rare hets in the *same gene* in **trans**, established by parent-of-origin (mat-only + pat-only) from the trio genotypes. A de novo second hit is a valid partner, but it **cannot be phased** from trio genotypes (~50/50 cis/trans), so such a pair is emitted with a `unphased_denovo_partner` flag rather than as confirmed trans. Read-backed **WhatsHap** phasing — which would resolve it — is a **TARGET**, not wired (the binary ships in the image but no step invokes it).
- **De novo (secondary cross-reference):** detected via GATK **`hiConfDeNovo`** (child-membership checked via `annotations.is_hiconf_denovo_for`), then **re-verified** with DP/AB plus **parental cleanliness** (each parent alt AD ≤ 1, DP ≥ 10). The former gnomAD absent/singleton + `nhomalt` condition is **retired** — the cache carries no `nhomalt` ([limitations.md](limitations.md) §3). It is **not** the driver here — dedicated bespoke machinery handles de novo filtering and review.
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

> **Multiallelic caveat — why Step 4 passes `--keep-sum AD`.** That formula is only correct if
> `AD[0]` really counts "reads not supporting this ALT". `bcftools norm -m-` subsets the
> `Number=R` AD array per allele and **discards the other ALT's reads**, so a non-ref/non-ref
> (`1/2`) genotype splits into two legs with `ref_ad ≈ 0` and AB reads **≈ 1.0 on both** — the het
> band then rejects a genuine trans compound het, silently. `--keep-sum AD`
> (`04_subset_and_annotate_trios.sh`) folds the other ALT's reads back into `AD[0]`, restoring the
> intended semantic: the same legs read AB 0.487 / 0.513, and a biallelic het is unchanged.
> Regression-tested by the `GENECH2` case in `tests/integration/`.

| Gate | Default | Notes |
|------|---------|-------|
| FILTER | `PASS` only | site-level |
| Genotype quality | **GQ ≥ 20** (from PP) | flag `lowGQ` at GQ < 20 |
| Depth | **DP ≥ 10** per sample | **DP ≥ 20** preferred for de novo |
| Het allele balance | **0.25 ≤ AB ≤ 0.75** | common practical band 0.2–0.8 |
| Hom-alt allele balance | **AB ≥ 0.90** | |
| Hom-ref allele balance | **AB ≤ 0.10** | |

> **Missing-genotype semantics — `strict_gt=True` is load-bearing.** Every "is this parent a
> no-call?" test in Step 5 assumes a no-call is distinguishable from a confident `0/0`. cyvcf2's
> **default does not give you that**: with `strict_gt=False` a *half*-called genotype such as `0/.`
> is reported as **hom-ref**, not unknown (its `as_gts` helper: "if a single allele is missing e.g
> `0/.` it's still encoded as hom ref because it has no alts"). That would silently satisfy the de
> novo requirement that both parents be confidently hom-ref, and would make a half-called parent
> look like an affirmative non-carrier when establishing compound-het trans. `05_inheritance_screen.py`
> therefore opens every trio VCF with `strict_gt=True`, so any `.` in a genotype reads as a no-call.
> (Step 0's QC pass does not set it: for Mendelian-error counting, treating `0/.` as hom-ref is the
> conservative direction.)

**Sample / pedigree QC:**

- **Peddy** confirms reported vs. genotype-inferred sex (chrX heterozygosity) and pedigree relatedness via IBS0 / Rel; parent–child pairs should show **IBS0 ≈ 0** and **relatedness ≈ 0.5**. Elevated IBS0 signals a sample swap or mislabeled trio → block downstream inheritance logic. Kids First runs Peddy explicitly.
- **Genome-wide Mendelian-error rate < 2%** is the concordance backstop: higher signals a bad trio/swap; localized MIE clusters flag CNV/UPD.
- **Contamination gate.** kid/dad/mom roles are already peddy-verified upstream, so Step 0's job is to catch the *less-curated* trios. It flags a trio if any member is contaminated: **verifyBamID2 `FREEMIX` > 0.05** when a directory of `*.selfSM` files is supplied (`resources.selfsm_dir`; this mirrors the group's DNM freemix QC), otherwise a **VCF-only CHARR** proxy — the reference-allele read fraction at high-quality homozygous-ALT SNV sites, thresholded at **> 0.02**. An uncontaminated sample sits near 0 (sequencing error + reference bias); cross-sample contamination injects reference reads at hom-alt sites and raises it. This matters because 1–3% contamination turns hom-ref → apparent-het, manufacturing false inherited hets / comp-het second hits. (CHARR: Lu et al., *Am J Hum Genet* 2023, [PMC10716339](https://pmc.ncbi.nlm.nih.gov/articles/PMC10716339/).)

> **Known failure mode (gnomAD-prior suppression).** Because `CalculateGenotypePosteriors` folds a population prior into every genotype, a genuinely rare/private pathogenic variant can be *pushed toward hom-ref* by a low-AF prior, suppressing a real finding. For candidate high-priority variants, **cross-check the pre-refinement PL/GT** to ensure refinement did not down-weight a true rare call.

---

## 3. Inheritance models

Genotype-level rules assume a mother–father–child trio VCF on GRCh38 with GT, GQ, DP, AD, and PL/PP available. Rarity gates below are stated as grpmax `faf95` because that is the **target** field the literature cites; the pipeline's **actual** filter field is the grpmax **proxy** (a point estimate — see the banner at the top and [allele_frequency.md](allele_frequency.md)). The cutoff values are identical either way. See also [inheritance engine tooling](#5-inheritance-engines--tooling).

The four first-class inheritance modes are **dominant (inherited het)**, **autosomal recessive homozygous**, **compound het (in trans)**, and **X-linked recessive**. **De novo** is retained only as a lightweight *secondary* cross-reference ([§3.5](#35-de-novo-secondary-cross-reference)); **mtDNA heteroplasmy** is out of scope ([§3.6](#36-mitochondrial-chrm--out-of-scope)).

### 3.1 Dominant — rare inherited heterozygous

- **Genotypes:** child `0/1`, with the alt allele **transmitted from ≥ 1 parent**. A transmitting parent may be `0/1` **or `1/1`** — a hom-alt parent still carries a transmissible alt (consanguinity, a common-ish recessive allele, an affected parent, or — on chrX — a hemizygous father rendered `1/1` by a diploid caller). Parent-of-origin is recorded as `flags=origin=mat|pat|both`:
  - exactly one parent carrying ⇒ `mat` / `pat`;
  - one parent `1/1` and the other `0/1` ⇒ **deterministic**, because a `1/1` parent transmits the alt obligately, so a het child took the alt from it and the ref from the other (`mom 1/1 × dad 0/1` ⇒ `mat`; `mom 0/1 × dad 1/1` ⇒ `pat`);
  - both parents `0/1` ⇒ a genuine 50/50 ⇒ `both` (and, being unphaseable, never used for compound-het pairing, [§3.3](#33-compound-heterozygous-two-hets-in-trans));
  - both parents `1/1` with a het child is a Mendelian error ⇒ no call.
  A child het with both parents `0/0` is de novo, not this mode (see [§3.5](#35-de-novo-secondary-cross-reference)).
- **Requirement:** the variant must be **functional** and **not** part of a *phase-confirmed* compound-het pair in the same gene ([§3.3](#33-compound-heterozygous-two-hets-in-trans)) — a single inherited het, not one half of a biallelic hit. A pair that is only *inferred* (a de-novo partner) does **not** suppress the dominant call.
- **QC:** child confident het (GQ ≥ 20, DP ≥ 10, AB 0.25–0.75) *and* each carrying parent confident **on its own zygosity band** — the het band for `0/1`, the hom-alt band (AB ≥ 0.90) for `1/1`.
- **`origin_unverified` flag:** when only one parent carries and the *other* parent is not an affirmative, QC-passing `0/0` (a no-call, or a `0/0` failing GQ/DP/AB), the call is still emitted but flagged — that parent might silently carry the allele too, so the recorded origin is inferred rather than established.
- **Rarity:** external gnomAD grpmax **proxy AF < 1e-4** (dominant rarity gate; same as de novo). Target field is `faf95`; see the banner.
- **Why it matters:** heterozygous inherited variants are individually low-specificity, but become compelling when they **recur across multiple distinct individuals in the same gene**. Gene consolidation ([§5](#5-inheritance-engines--tooling), [gene_burden.md](gene_burden.md)) tallies dominant-het carriers per gene and ranks recurrent genes first, **weighted by gene constraint** (a recurrent het in a haploinsufficient gene is the most compelling). Do **not** reject solely because an unaffected parent carries it — reduced/age-dependent penetrance is expected ([§6](#6-co-segregation--penetrance-modifiers-never-hard-filters)).

### 3.2 Autosomal recessive — homozygous

- **Genotypes:** child `1/1`; **both** parents carrying — each `0/1` **or `1/1`** (a hom-alt parent is accepted for consanguinity, a common-ish recessive allele, or an affected parent). A parent no-call means inheritance is unestablished and the call is not made.
- **QC:** child confident hom-alt (AB ≥ 0.90, DP ≥ 10); each carrying parent confident **on its own zygosity band** — het band for `0/1`, hom-alt band for `1/1` (GQ ≥ 20, DP ≥ 10).
- A `1/1` child with a `0/0` parent is a Mendelian error → suspect a hemizygous "false hom" (deletion on the other allele) or **UPD** (see [§4](#4-upd-imprinting--mosaicism-flags-not-primary-calls)).
- **Rarity:** per-allele grpmax **proxy AF < 1e-2** (permissive discovery default; target field is `faf95`), with a stricter **1e-3** high-confidence tier. Applied per variant, not per gene. **Do not** down-weight biallelic candidates using pLoF constraint ([gene_constraint.md](gene_constraint.md)).

### 3.3 Compound heterozygous (two hets in trans)

- **Requirement:** two rare het variants in the same gene on **opposite** haplotypes (trans). Cis pairs are non-causal.
- **Trio phasing (the ONLY method wired):** variant A from mother (mat `0/1`, pat `0/0`) and variant B from father (pat `0/1`, mat `0/0`) ⇒ trans by descent. Both variants tracing to the same parent ⇒ cis ⇒ reject. Only `mat` and `pat` legs pair; a `both`-origin leg (both parents `0/1`) is unphaseable and is never paired.
- **The non-transmitting parent must be *affirmatively* hom-ref.** Trans-by-descent rests on that parent *not* carrying the allele, so a no-call or an unqualified `0/0` (allele dropout) there is not evidence of anything — if that parent silently carries it too, the "trans" pair may be **cis**. Such a pair is still emitted (never-drop) but flagged **`origin_unverified`**. Note this needs `strict_gt=True` on the VCF reader: with cyvcf2's default a half-called `0/.` is reported as hom-ref rather than as a no-call ([§2](#2-cross-cutting-genotype--sample-qc-apply-before-mode-logic)).
- **A `1/1` transmitting parent is deterministic, not ambiguous.** `mom 1/1 × dad 0/1` with a het child is unambiguously maternal (see [§3.1](#31-dominant--rare-inherited-heterozygous)), so it **does** pair in trans with a paternal leg. This matters most on chrX, where a diploid caller renders every hemizygous carrier father as `1/1`.
- **De novo second hit:** a legitimate partner biologically, but **unphaseable** from trio genotypes (it may sit cis or trans with the inherited hit at ~50/50), so that pair is emitted with a `unphased_denovo_partner` flag — a candidate to confirm, not a confirmed biallelic hit. Because it is unconfirmed it does **not** suppress the dominant call on the inherited leg.
- **Read-backed phasing (WhatsHap) — TARGET, not implemented.** `whatshap` is pinned in the image but **no step invokes it**; there is no read-backed phasing today. When wired it would resolve phase directly where both variants lie within one read/fragment (and uniquely combine read-based *and* pedigree phasing), which is exactly what would settle the de-novo-partner case above. Caveat for that future work: WhatsHap drops variants with missing or Mendelian-inconsistent parental genotypes, lowering the phasing rate when parental data is incomplete.
- **QC:** both variants must independently pass het QC (GQ ≥ 20, DP ≥ 10, AB 0.25–0.75) in the child *and* the transmitting parent.
- **Rarity:** per-allele grpmax **proxy AF < 1e-2** (1e-3 high-confidence; target field is `faf95`).

### 3.4 X-linked / hemizygous

- **X-recessive (affected male):** child hemizygous alt (`1/1` from a diploid caller, or `1` if ploidy-aware — cyvcf2 maps both to hom-alt), mother `0/1` **or `1/1`** carrier. **The father's chrX is NOT required and is not examined**: he transmits his Y to a son, so his chrX is never transmitted. An affected or carrier father, or a paternal chrX no-call, must not drop the call — if he does carry, it is recorded as `flags=father_carries_x_allele`, never as a filter.
- **X-recessive (affected female):** child `1/1`, mother `0/1`/`1/1` carrier, **father `1/1`** — he *does* transmit his single X to every daughter, so a homozygous daughter requires a hemizygous-affected father. All three are QC'd on their own zygosity band.
- **X-dominant is not a separate emitted mode.** It is genotypically indistinguishable from the modes above and is covered by them: a **female** proband's non-PAR chrX het flows through the ordinary het collector and is emitted as `dominant` (with parent-of-origin), while a **male** hemizygote is emitted as `x_linked_recessive`. Filter on `chrom` to separate X calls from autosomal ones; the mode label does not encode it. Watch for male-lethal patterns when interpreting.
- **Hemizygous chrX / chrY (male, outside PAR):** males are haploid, so a diploid caller emits `0/0` or `1/1` — **a het call in a male non-PAR region is a QC red flag** (mapping artifact, PAR misplacement, or XXY) and is **dropped** from the het collector on both chrX and chrY. PAR1/PAR2 are diploid in both sexes and are correctly routed through the autosomal models. QC the hemizygous alt call at DP ≥ 10, GQ ≥ 20, AB ≥ 0.90. Peddy/Step-0 sex inference gates this: a trio whose child sex is unresolved has its sex chromosomes skipped entirely rather than assumed female.
- **chrY yields no inherited call.** The hemizygous models above are keyed on the **mother's** genotype, which is meaningful only on chrX — on chrY the father is the sole transmitter and the mother has no Y at all. chrY is therefore excluded from both the hemizygous de novo path and `x_linked_recessive` (`male_x_chrx` in `05_inheritance_screen.py`), so male non-PAR chrY records produce no rows. Routing them through the X logic would have reported father-to-son transmission as *de novo* and let recurrent Yq / X-transposed mismapping artifacts — which carry no gnomAD chrY AF, so the rarity gate passes unconditionally — accumulate in Step 6's X-linked tier. **There is no Y-linked inheritance model**; the clinical cost is negligible (Y-linked Mendelian SNV disease is essentially confined to spermatogenic failure, and the lesions there are CNVs, already a documented blind spot).

### 3.5 De novo (secondary cross-reference)

De novo is **not** the driver of this pipeline: dedicated bespoke machinery handles de novo filtering **and** review. Here it is retained only as a lightweight secondary cross-reference against the inherited-variation signal, and counted **separately** from the inherited modes in gene consolidation ([§5](#5-inheritance-engines--tooling)).

- **Genotypes:** child `0/1`, mother `0/0`, father `0/0` (autosomal); hemizygous de novo on male **chrX** is `0 → 1`, requiring only the transmitting **mother** to be hom-ref (see [X-linked](#34-x-linked--hemizygous)). **chrY is excluded from this path** — the mother-keyed rule is inverted there, so a father-transmitted Y allele would otherwise be reported as de novo. Both parents must be *confidently* hom-ref: a no-call does not qualify, and `strict_gt=True` is what makes a half-called `0/.` register as a no-call rather than as hom-ref ([§2](#2-cross-cutting-genotype--sample-qc-apply-before-mode-logic)).
- **Detection:** GATK `hiConfDeNovo` present (all three trio-member GQ ≥ 20), with **child membership** verified via `annotations.is_hiconf_denovo_for` (confirming the tag applies to *this* child). Use `loConfDeNovo` (child GQ ≥ 10) only as a lower-sensitivity tier.
- **Re-verify (the tool does not):** child **DP ≥ 20**, het **AB 0.25–0.75**, all three GQ ≥ 20, and **parental cleanliness** — each parent alt AD ≤ 1 with DP ≥ 10 (a parental alt fraction of a few percent suggests inherited or parental mosaicism, not de novo).
- **Rarity:** external gnomAD grpmax **proxy AF < 1e-4** (target: `faf95`). The `nhomalt` condition is **RETIRED, not enforced** — no `nhomalt` field exists in the VEP cache, and the retired `filters.denovo.require_gnomad_absent_or_singleton` key only ever implemented `nhomalt > 1`, a homozygote-count test rather than the allele-count test its name promised ([limitations.md](limitations.md) §3). *(An AC-based **absent-or-singleton** test remains a target refinement — de novo is secondary here.)*
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
| **slivar** (v0.3.4) | Reference segregation engine — **pinned but NOT invoked** | The actual Step-5 engine is custom cyvcf2 (`05_inheritance_screen.py`); slivar ships in the image for ad-hoc work but no step calls it. Its PED-driven JS-expression helpers (dominant/inherited-het, recessive/hom-alt, `comphet`, `x_denovo`, `denovo`) are the reference design. Single static binary. |
| **WhatsHap** | Read-based + pedigree phasing | **TARGET — pinned in the image but not invoked.** Trans is resolved by trio parent-of-origin only; a de-novo-partner pair is flagged `unphased_denovo_partner` rather than phased. |
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
| Dominant (inherited het) | rare functional **het** transmitted from ≥ 1 parent; parent-of-origin recorded (mat / pat / both); **not** a comp-het partner | grpmax **proxy AF < 1e-4** (not `faf95`); recurrence-consolidated |
| Recessive / comp-het rarity | grpmax **proxy AF < 1e-2** per allele (not `faf95`) | **1e-3** high-confidence tier |
| Recessive hom | child `1/1` + both parents `0/1`, all GQ ≥ 20, DP ≥ 10 | 1/1 vs 0/0 parent → suspect deletion/UPD |
| Compound het | two rare hets, same gene, **trans** (parent-of-origin from trio GTs; WhatsHap is a TARGET, not wired); the non-transmitting parent must be an affirmative QC-passing `0/0`; a `1/1` transmitting parent is deterministic and **does** pair | de novo second hit is valid but **unphaseable** → `unphased_denovo_partner` (and does not suppress the dominant call); unobserved other-parent → `origin_unverified` |
| X-linked recessive | affected male = hemizygous + carrier mother (**father's chrX not required**, flagged if he carries); affected female = `1/1` + carrier mother + `1/1` father; sex-aware ploidy; **drop** male non-PAR chrX/chrY het calls; separate PAR/non-PAR | X-dominant is not a separate mode (female X het → `dominant`, male hemizygote → `x_linked_recessive`); **chrY yields no inherited call** |
| De novo (secondary) | GATK **`hiConfDeNovo`** (child-membership via `is_hiconf_denovo_for`); re-verify child DP ≥ 20, AB 0.25–0.75, each parent alt AD ≤ 1 / DP ≥ 10; grpmax **proxy AF < 1e-4** (not `faf95`; the gnomAD absent/singleton + `nhomalt` condition is **retired** — no `nhomalt` field) | cross-reference only; filtering/review in separate machinery. Optional: Poisson enrichment vs Samocha (P < 2.5e-6, BH q < 0.05) when a mutation-rate table is supplied |
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
