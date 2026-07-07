# Allele-Frequency Filtering & the Frequency Oracle

How this pipeline decides whether a variant is rare enough to be a high-priority candidate, using an external population reference rather than the untrustworthy internal cohort.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- **Oracle = gnomAD v4.1 (GRCh38)**, using the **joint (exome + genome) AF/AN**; heed the v4.1 exome/genome **discordance flag**.
- **Filter on grpmax `faf95`** (group-max filtering allele frequency, 95% CI lower bound) — *not* the point-estimate popmax/grpmax AF. Fall back to grpmax AF only where faf95 is missing.
- **Dominant / de novo:** keep if `faf95 < 1e-4`; for de novo additionally require absent-or-singleton in gnomAD and low `nhomalt`.
- **Recessive / compound-het:** keep if `faf95 < 1e-2` (permissive discovery default), with a stricter **`< 1e-3`** high-confidence tier; applied **per variant**, not per gene.
- **Hard benign (all modes):** drop if `faf95 ≥ 0.05` (ClinGen BA1) — never rescue.
- **Never** use internal cohort AC/AN as a population frequency: the non-joint per-trio design makes AN uninterpretable (absent genotype ≠ hom-ref). Internal recurrence is valid **only** as an artifact/blocklist signal.
- A **gene-specific ClinGen VCEP** BA1/BS1 or a **Whiffin/Ware maximum credible AF** overrides any generic cutoff.
- The rarity gate is a **screening filter**, distinct from the ACMG **PM2** criterion (applied at *Supporting* strength only). Passing the gate is not the same as "PM2 met."

## Why an external frequency oracle

This pipeline screens **GMKF Kids First per-trio VCFs** that are GATK Genotype-Refinement output but are **not jointly genotyped across the cohort**. That design makes internal allele counts unusable as a population frequency:

- **No consistent cohort-wide AN.** Each trio is called independently, so a variant's internal frequency reflects only 2–6 chromosomes; there is no shared denominator across trios.
- **Absence ≠ reference.** In a non-joint merge, a variant missing from another trio may be a no-call or low-depth site, not a confident hom-ref. Internal AC/AN therefore mis-estimates both numerator and denominator.
- **gnomAD provides the defensible denominator.** It is large, uniformly joint-genotyped, ancestry-resolved, and ships proper filtering-allele-frequency confidence intervals.

Internal data still has one legitimate frequency-adjacent use: **artifact detection**. A variant recurring across many unrelated trios is more likely a systematic sequencing/mapping artifact than a truly common allele. Use that as a panel-of-normals-style **blocklist** signal (tune the recurrence count `N` empirically), never as a population AF. See [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md) and [cohort_construction.md](cohort_construction.md).

## The reference dataset: gnomAD v4.1 (GRCh38)

- **Composition:** 730,947 exomes (416,555 UK Biobank + 314,392 non-UKB) plus 76,215 genomes, all unrelated, aligned to GRCh38/hg38. v4.1 is the current release (Apr 2024). The union callset is ~807k samples, but exome vs genome N differs per site.
- **v4.1 key fixes:** corrects the v4.0 allele-number (AN) bug; adds a **joint (combined exome + genome) AN and AF** at every site called in either data type; adds a **discordant-frequency flag** where a contingency/CMH test between exomes and genomes gives p < 1e-4 (~2.5% of variants).
- **Practical rule:** prefer the joint AF/AN and heed the discordance flag before trusting a single subset.

> VEP release r113 (Oct 2024) updated its built-in gnomAD annotation to **v4.1** (genomes + exomes), so a release-matched VEP cache can supply these fields directly; a standalone gnomAD VCF / Hail table is the alternative. See [tooling_and_reproducibility.md](tooling_and_reproducibility.md).

## Global AF vs grpmax vs FAF — use the right number

| Metric | What it is | Why we do / don't use it |
| --- | --- | --- |
| **Global AF** | AF across all samples | Dilutes an ancestry-enriched variant; a variant common in one group looks rare globally. **Do not filter on this.** |
| **grpmax AF** (formerly popmax) | Highest point-estimate AF across genetic-ancestry groups | Better than global, but a point estimate is noisy when a group's AN is small. |
| **FAF (faf95 / faf99)** | Lower bound of the 95% (or 99%) Poisson CI on the AF | The frequency you can be ≥95% confident the true AF is *at least*. Conservative for *filtering out* benign variants — you only exclude a variant as "too common" when confident it really is common. **This is our filter field.** |
| **grpmax FAF** | faf95 from the ancestry group with the highest FAF | The value ClinGen VCEPs use for BA1/BS1. |

**Founder-group exclusion.** gnomAD excludes bottlenecked/founder groups (Amish, Ashkenazi Jewish, Finnish, Middle Eastern) from grpmax FAF, because pathogenic founder alleles legitimately reach high frequency there and would wrongly inflate the filter. Rely on this exclusion; do **not** re-introduce those groups' frequencies into the gate.

Using the CI lower bound is deliberately conservative: it protects against false exclusion of true pathogenic alleles that happen to appear by chance in a small sample.

## Inheritance-mode–dependent rarity gates

The maximum tolerated frequency depends on inheritance mode. These are **screening defaults**, overridable in `config/config.example.yaml`, and are superseded by any gene-specific value.

| Candidate class | Filter field | Keep if | Notes |
| --- | --- | --- | --- |
| **Dominant / de novo** | grpmax faf95 | `< 1e-4` | For de novo, additionally require absent-or-singleton in gnomAD and low `nhomalt`. Treat any appreciable gnomAD frequency as strong evidence against a de novo call. |
| **Recessive / compound-het** | grpmax faf95 | `< 1e-2` (discovery) | Stricter high-confidence tier at `< 1e-3`. Applied **per variant**, not per gene, and compound-het/biallelic aware. |
| **Hard benign (all modes)** | grpmax faf95 | drop if `≥ 0.05` | ClinGen general-purpose **BA1**; never rescue. |

The literature range for the generic recessive cutoff spans roughly 1e-3 to 1e-2 (5e-3 is a commonly cited midpoint); this pipeline uses the permissive 1e-2 discovery default with a 1e-3 high-confidence tier so that biallelic candidates are not lost early. For dominant conditions, ClinGen general practice sits near grpmax faf95 `< 1e-4` absent a gene-specific value.

**Gene-specific override.** Where a ClinGen VCEP publishes calibrated BA1/BS1 values (e.g. cardiomyopathy, RASopathy), use those instead of the generic cutoff. The RASopathy VCEP, for instance, calibrates BA1 grpmax faf well below the generic 0.05, illustrating why gene-specific tables win when they exist. See [clinical_classification.md](clinical_classification.md).

## Maximum credible population allele frequency (Whiffin/Ware)

Whiffin et al. (Genet Med 2017) define the maximum AF a variant can have and still plausibly cause the disease:

```
maxAF = (prevalence × max allelic heterogeneity × 1/penetrance-factor)
        / (inheritance-adjusted allele count)
```

Inputs: **disease prevalence**, **allelic heterogeneity** (max fraction of cases from one variant), **genetic heterogeneity** (max fraction from one gene), **penetrance**, and **inheritance mode** (monoallelic vs biallelic). Lower penetrance → higher tolerated AF. Compare the variant's observed **grpmax faf95** to this maxAF; if `faf95 > maxAF`, filter out. In cardiomyopathy this removed ~two-thirds of candidates without losing true positives. Compute per gene with the CardioDB calculator. This is the preferred override for [pediatric_cancer.md](pediatric_cancer.md) predisposition genes where a defensible prevalence/penetrance model exists.

## Rarity gate vs the ACMG PM2 criterion

These are related but **not the same**, and the pipeline keeps them separate:

- The **rarity gate** above is a screening filter that decides whether a variant continues through prioritization.
- **PM2** is an ACMG evidence criterion, met when a variant is **absent from, or at extremely low frequency in, population controls** (a low count is tolerated for recessive). ClinGen SVI (2020) **downgraded PM2 to *Supporting* by default**. It is evidence toward a classification, not a hard filter — do not let it alone drive a call, and do not equate "passed the rarity filter" with "PM2 met." PM2 handling lives in [clinical_classification.md](clinical_classification.md).

## Annotation & filtering commands (generic, parameterized)

Annotate with the release-matched VEP cache (gnomAD v4.1 fields) — no baked-in paths:

```bash
vep \
  --offline --cache --dir_cache "${VEP_CACHE_DIR}" \
  --assembly GRCh38 --fasta "${REF_FASTA}" \
  --vcf --compress_output bgzip \
  --input_file "${IN_VCF}" --output_file "${OUT_VCF}" \
  --af_gnomade --af_gnomadg \
  --plugin gnomADc  # or annotate grpmax faf95 from a bind-mounted gnomAD VCF/Hail table
```

Filter on grpmax faf95 with a fallback to grpmax AF when faf95 is absent (dominant example, `1e-4`):

```bash
bcftools view -i \
  '(INFO/gnomad_grpmax_faf95 != "." && INFO/gnomad_grpmax_faf95 < 1e-4) ||
   (INFO/gnomad_grpmax_faf95 == "." && INFO/gnomad_grpmax_af   < 1e-4)' \
  "${ANNOTATED_VCF}" -Oz -o "${RARE_DOMINANT_VCF}"
```

Never-rescue benign cutoff (BA1) and discordance-flag awareness applied together:

```bash
# Drop common variants (BA1); optionally down-rank sites carrying the exome/genome discordance flag
bcftools view -e \
  'INFO/gnomad_grpmax_faf95 >= 0.05' \
  "${RARE_VCF}" -Oz -o "${RARE_NOT_BENIGN_VCF}"
```

Field names above are placeholders — align them to whatever your annotation source emits, and pin the gnomAD/VEP release. See [functional_annotation.md](functional_annotation.md) for the downstream consequence layer.

## Known limitations

- **SNV/indel only.** This frequency logic applies to short-variant calls. CNV/SV are a real blind spot (10–15% of pediatric-cancer and rare-disease diagnoses), and gnomAD SNV FAF does not address them; a future GATK-gCNV / Manta / ExomeDepth module is required. See [pipeline_design.md](pipeline_design.md).
- **Pseudogene / segmental-duplication genes** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) are low-confidence from short reads; gnomAD frequencies in those paralogous regions can themselves be unreliable, so flag those regions rather than trusting the AF.
- **Proband mosaicism.** Low-VAF post-zygotic calls are handled by the genotype-QC layer, not here, but note that a genuinely rare pathogenic mosaic call must survive both the AB band and this rarity gate.
- **Founder/bottlenecked populations** are excluded from grpmax FAF by gnomAD; a pathogenic founder allele common in an excluded group is intentionally not counted against the filter.

## Recommended defaults (this pipeline)

| Parameter | Default | Notes |
| --- | --- | --- |
| Frequency oracle | gnomAD **v4.1** (GRCh38), **joint** AF/AN | Heed exome/genome **discordance flag**. |
| Filter field | grpmax **faf95** | Fall back to grpmax AF only where faf95 absent. |
| Dominant / de novo keep | faf95 `< 1e-4` | De novo also: absent-or-singleton + low `nhomalt`. |
| Recessive / comp-het keep | faf95 `< 1e-2` (discovery); `< 1e-3` high-confidence | Per variant, not per gene. |
| Hard benign (all modes) | faf95 `≥ 0.05` → drop | ClinGen BA1; never rescue. |
| Internal cohort AC/AN | **Never** as population AF | Use only as recurrence/artifact blocklist signal. |
| PM2 | Supporting strength only | Evidence, not the screening gate. |
| Gene-specific override | ClinGen VCEP BA1/BS1 or Whiffin/Ware maxAF | Overrides the generic cutoffs whenever available. |

All values are configurable defaults, not immutable law. Prefer a gene-specific ClinGen VCEP value over any generic cutoff.

## Sources

- gnomAD v4.1 release notes: https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/
- gnomAD FAF help: https://gnomad.broadinstitute.org/help/faf
- ClinGen guidance to VCEPs on using gnomAD v4 (BA1/BS1, grpmax FAF), March 2024: https://clinicalgenome.org/site/assets/files/9445/clingen_guidance_to_vceps_regarding_the_use_of_gnomad_v4_march_2024.pdf
- ClinGen SVI PM2 recommendation v1.0 (2020; Supporting default): https://clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf
- Whiffin et al., "Using high-resolution variant frequencies to empower clinical genome interpretation," Genet Med 2017: https://www.nature.com/articles/gim201726.pdf (PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5563454)
- Maximum credible allele frequency calculator (CardioDB): https://www.cardiodb.org/allelefrequencyapp/
- ClinGen RASopathy ACMG specification (calibrated grpmax faf thresholds): https://www.sciencedirect.com/science/article/pii/S2949774425014694
- ACGS 2024 UK variant-classification best-practice guidelines: https://www.genomicseducation.hee.nhs.uk/wp-content/uploads/2024/08/ACGS-2024_UK-practice-guidelines-for-variant-classification.pdf
