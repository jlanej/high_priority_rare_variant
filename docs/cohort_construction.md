# Building a Cohort Site List from Non-Jointly-Genotyped Trios

How to combine independently-called GMKF Kids First per-trio VCFs into a defensible cohort site list without fabricating genotypes or inventing a population frequency.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

## TL;DR

- The trios are **variant-only, per-trio** VCFs (GATK Genotype-Refinement output), **not jointly genotyped** across the cohort. Any merge that fills absent cells invents genotypes.
- **Absent record ≠ hom-ref.** "No record for this sample at this site" conflates true hom-ref with a genuine no-call (uncovered / unassessed). The GVCF reference-block info that would disambiguate is already gone.
- Therefore an internal cohort **AC/AN is uninterpretable and is never a population frequency.** Do **not** `bcftools merge` these files into a genotype matrix, and never use `--missing-to-ref`.
- Build a **site-only union** instead: normalize each trio, drop genotypes, then union with dedup. This is legitimate for "which loci appear anywhere" and a coarse trio-**recurrence** count — nothing more.
- The frequency oracle is **external gnomAD v4.1**, read through the **VEP cache** as a **grpmax point-estimate proxy** — **not** `faf95`, which the cache cannot supply (no AC/AN) and which remains the TARGET (see [allele_frequency.md](allele_frequency.md), [limitations.md §2](limitations.md)). The argument below is unaffected: the point is that the denominator is *external*, never internal. Internal recurrence is valid **only** as an artifact/blocklist signal.
- Normalize **before** any union: `norm -m-` (split multiallelics), left-align/`-c` against the **exact GRCh38 build** the trios were called on, then `view -G` to drop genotypes.
- If a real genotype matrix is ever required, re-run true **joint genotyping from the per-trio gVCFs** (GLnexus or GATK `GenomicsDBImport`→`GenotypeGVCFs`) — you cannot recover it from variant-only VCFs.
- Pin **bcftools/htslib 1.23** (project default, overridable; see [tooling_and_reproducibility.md](tooling_and_reproducibility.md)).

## The core hazard: a naive merge fabricates genotypes and corrupts AC/AN

Each trio VCF is a **variant-only** record set for that trio. When `bcftools merge` combines them, any sample that had no record at a site gets a genotype the caller never emitted. The absent record is treated as reference, but "absent" collapses two genuinely different states:

1. the site was truly hom-ref in that trio, and
2. the site was never assessed / had insufficient coverage (a real no-call).

Because each trio was called **independently**, a site absent in trio B may be a no-call there, yet the merge can still report it as `0/0`. Whether a filled cell becomes `0/0` or `./.` is not reliable when inputs are not jointly genotyped — a long-standing, documented ambiguity in bcftools ([#1891](https://github.com/samtools/bcftools/issues/1891), [#2236](https://github.com/samtools/bcftools/issues/2236), [#402](https://github.com/samtools/bcftools/issues/402)). The GVCF `<NON_REF>`/reference-block information that would disambiguate hom-ref from no-call is already gone once you only have per-sample variant VCFs.

Consequences for any cohort-derived frequency:

- **AN is unknowable.** You cannot distinguish N confidently-hom-ref samples from N uncalled samples, so any denominator is fiction.
- **AC is inflated or deflated** depending on how merge resolves each fill. `--missing-to-ref` forces `0/0`, which *guarantees* AN inflation and hides genuine missingness.
- **Missingness metrics are meaningless**, defeating one of their main QC/batch-filter uses.

**Bottom line:** an internal cohort AC/AN computed from a non-joint merge is **not a valid population frequency** and must never be used as a frequency oracle. This is the single most important rule on this page.

## Two legitimate constructs — keep them strictly separate

### 1. True joint genotyping (only possible from gVCFs)

If you retain per-trio **gVCFs** (with reference blocks), you can produce a real project-level VCF with a genotype for every sample at every variant site:

- **GATK**: `GenomicsDBImport` → `GenotypeGVCFs` yields a complete genotype matrix; each step is single-threaded and needs external scatter/parallelization ([GATK how-to](https://gatk.broadinstitute.org/hc/en-us/articles/360035889971)).
- **GLnexus**: internally multithreaded, benchmarked ~8× faster than the GATK pair on a 32-vCPU node (0.84 h vs 6.83 h) and scaling better; diploid-only and computes fewer QC metrics ([Yun et al., *Bioinformatics* 2020, PMC8023681](https://pmc.ncbi.nlm.nih.gov/articles/PMC8023681/)).

Only a joint callset yields trustworthy cohort AC/AN/missingness. **This pipeline does not assume gVCFs are available**, so joint genotyping is an optional escalation path, not the default flow.

### 2. Site-only union list (the default here — variant VCFs only)

With variant-only VCFs you cannot recover cohort genotypes, so build a **sites list**: the set of distinct variant loci across trios, genotypes dropped. This is legitimate for:

- defining "which variants appear anywhere in the cohort" (annotation targets, candidate loci), and
- a coarse **recurrence** count (in how many trios a site was seen).

It is **never** a frequency denominator. Recurrence feeds the artifact/blocklist logic described below, and the per-pedigree screening in [gene_burden.md](gene_burden.md) works from trios directly rather than from a synthetic genotype matrix.

## Building a defensible site-only union

### Normalization is non-negotiable — do it before any union

Otherwise the "same" variant carries multiple representations and dedup/annotation silently fail.

1. **Split multiallelics to biallelic**: `bcftools norm -m-` (default splits all types).
2. **Left-align and normalize indels against the reference**: requires `-f <GRCh38.fa>`; left-alignment/normalization only fires when `--fasta-ref` is supplied. Use the **exact GRCh38 build the trios were called against** — matching contig naming (e.g. `chr1`) and including the alt/decoy contigs per Kids First.
3. **Check REF against the reference** for mismatches. The pipeline uses `-c w` (**warn**, never silently rewrite) so a build/contig mismatch stays visible rather than being masked by an automatic `-c s` rewrite.
4. **Drop genotypes** to make the file sites-only: `bcftools view -G`.
5. **Sort, bgzip, index** each file (`bcftools sort`; `bcftools index -t`).
6. **Union with dedup**: `bcftools concat -a -D` on the normalized sites files removes duplicate identical records; follow with `sort` + `norm -d exact` to collapse residual duplicates that differ only in representation.

See [bcftools norm/view/concat docs](https://samtools.github.io/bcftools/bcftools.html).

### Reconciling FILTER / INFO across trios

FILTER and INFO were computed **per trio** and are not comparable across the cohort. For a sites list, either strip them or compute any recurrence tag yourself with semantics you control:

- Strip the incomparable per-trio INFO: `bcftools annotate -x INFO`. (The pipeline keeps the FILTER
  column — records are already reduced to PASS/`.` by an upstream `bcftools view -f`, so there is
  nothing incomparable left there to strip.)
- If you keep a recurrence count, compute it with `bcftools +fill-tags` (or your own `INFO/NS`-style counter over the input files) — **do not** carry per-trio AC/AN into the cohort file as though it were population frequency.

Per-trio genotype/QC gates (GQ, DP, allele balance, FILTER=PASS, refined-PP handling, de novo re-verification) belong upstream, at the trio level — see [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md). Reconciling FILTER here is only about not laundering incomparable per-trio calls into a cohort-level field.

### Example commands (generic, parameterized — no real paths)

Per-trio prep (run on each trio file):

```bash
bcftools norm -m- -f "${REF_FASTA}" -c w "${TRIO_VCF}" -Ou \
  | bcftools view -G -Oz -o "${TRIO_SITES}" \
  && bcftools index -t "${TRIO_SITES}"
```

Cohort site-only union with dedup:

```bash
bcftools concat -a -D -Ou ${TRIO_SITES_GLOB} \
  | bcftools sort -Ou \
  | bcftools norm -d exact -f "${REF_FASTA}" -Oz -o "${COHORT_SITES}" \
  && bcftools index -t "${COHORT_SITES}"
```

Strip incomparable per-trio annotations (retain only a self-computed recurrence tag if desired):

```bash
bcftools annotate -x INFO "${COHORT_SITES}" -Oz -o "${COHORT_SITES_CLEAN}"
```

Escalation only, and only if gVCFs exist — a real genotype matrix via joint genotyping:

```bash
# GLnexus (fast, diploid) — from per-trio gVCFs, NOT variant-only VCFs
glnexus_cli --config gatk ${TRIO_GVCF_GLOB} > cohort.joint.bcf
# or GATK: GenomicsDBImport -> GenotypeGVCFs (scatter per-interval externally)
```

## Internal frequency: an artifact filter, never a population frequency

Internal cohort recurrence has exactly one valid use: a **batch/artifact filter**. A site seen in an implausibly large fraction of your trios — given that gnomAD says it is rare — flags a recurrent sequencing/mapping/annotation artifact or a batch effect; remove or down-weight it. It must **never** substitute for gnomAD as the population-rarity denominator, precisely because the non-joint merge makes AN uninterpretable.

The population-rarity decision lives entirely with the external oracle: **gnomAD v4.1** (GRCh38-native; 730,947 exomes + 76,215 genomes). **As implemented**, those frequencies are read from the **VEP cache** and reduced to a **grpmax proxy** — the max *point-estimate* AF over the grpmax-eligible groups (AFR/AMR/EAS/NFE/SAS). Filtering on the 95%-CI-based **grpmax `faf95`** (and heeding the v4.1 **exome/genome discordance flag**) is the **TARGET**: both need the gnomAD sites VCF, since the cache carries no AC/AN and no FAF ([limitations.md §2](limitations.md)). gnomAD v4.1 (April 2024) fixed the v4.0 allele-number bug and added the harmonized joint AF/AN/FAF. Full frequency logic, gene-specific ClinGen VCEP BA1/BS1 overrides, and the FAF rationale are in [allele_frequency.md](allele_frequency.md); ClinGen issued guidance on using gnomAD v4 for BA1/BS1 evidence in March 2024.

## Known scope limitations (stated honestly)

- **No CNV/SV from this construct.** A site-only SNV/indel union cannot represent copy-number or structural events. ~10–15% of pediatric-cancer and rare-disease diagnoses are CNV/SV (single-exon *RB1*/*SMARCB1*/*DICER1*/*NF1* deletions, *PMS2* rearrangements). That is a future module (GATK-gCNV / Manta / ExomeDepth), not covered here.
- **Pseudogene / segmental-duplication regions** (*PMS2*/*PMS2CL*, *CYP21A2*, *SMN1/2*, *NEB*, *GBA*) are low-confidence from short reads. Normalization and union do nothing to fix mismapped paralog calls — flag those regions rather than trusting a recurrence count there.
- **No cohort genotype matrix means no cohort-level missingness/QC** (e.g. site-level call-rate filtering, cohort Hardy-Weinberg). Those require true joint genotyping.
- **Recurrence is coarse.** It counts trios in which a site was seen, gated by per-trio calling sensitivity; it is not, and cannot be corrected into, an allele count.

## Recommended defaults (this pipeline)

| Decision | Default | Notes |
|---|---|---|
| Cohort construct | **Site-only union** (genotypes dropped) | Variant-only trios; no genotype matrix synthesized |
| `bcftools merge` of trios | **Prohibited** for cohort AC/AN | Fabricates genotypes; AN uninterpretable |
| `--missing-to-ref` | **Never** | Guarantees AN inflation, hides missingness |
| Multiallelic handling | `bcftools norm -m-` | Split to biallelic before union |
| Left-align / normalize | `norm -f <GRCh38.fa>` on the exact build the trios used | Contig-name matched; alt/decoy per Kids First |
| REF check | `norm -c w` (warn on REF/ALT mismatch; never rewrite) | Surface build/contig problems, don't mask them |
| Drop genotypes | `bcftools view -G` | Makes the file sites-only |
| Union + dedup | `concat -a -D` → `sort` → `norm -d exact` | Collapse residual representation duplicates |
| Per-trio INFO | Strip (`annotate -x INFO`) | Incomparable across independent calls (FILTER kept — already PASS-filtered upstream) |
| Internal recurrence | Self-computed; **artifact/blocklist signal only** | Never a population frequency |
| Frequency oracle | **External gnomAD v4.1** via the VEP cache, **grpmax proxy** (point estimate). `faf95` = TARGET, not implemented | See [allele_frequency.md](allele_frequency.md), [limitations.md §2](limitations.md) |
| Genotype matrix (if ever needed) | Joint-genotype from **gVCFs** (GLnexus / GATK) | Cannot be recovered from variant-only VCFs |
| Tool versions | **bcftools/htslib 1.23** (default, overridable) | See [tooling_and_reproducibility.md](tooling_and_reproducibility.md) |

All values are **defaults**, overridable in `config/config.example.yaml`. A gene-specific ClinGen VCEP threshold (applied downstream in the frequency/classification layers) overrides any generic cutoff.

## Sources

- gnomAD v4.0 release notes — https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/
- gnomAD v4.1 release notes — https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/
- gnomAD filtering allele frequency (FAF) help — https://gnomad.broadinstitute.org/help/faf
- ClinGen guidance on gnomAD v4 (March 2024) — https://clinicalgenome.org/site/assets/files/9445/clingen_guidance_to_vceps_regarding_the_use_of_gnomad_v4_march_2024.pdf
- bcftools manual (norm, view -G, concat, merge, annotate) — https://samtools.github.io/bcftools/bcftools.html
- bcftools merge missing-vs-hom-ref ambiguity — https://github.com/samtools/bcftools/issues/1891 · https://github.com/samtools/bcftools/issues/2236 · https://github.com/samtools/bcftools/issues/402
- GATK Consolidate GVCFs / GenotypeGVCFs how-to — https://gatk.broadinstitute.org/hc/en-us/articles/360035889971
- GLnexus + DeepVariant cohort calling (Yun et al., *Bioinformatics* 2020) — https://pmc.ncbi.nlm.nih.gov/articles/PMC8023681/
- Broad WARP JointGenotyping pipeline — https://broadinstitute.github.io/warp/docs/Pipelines/JointGenotyping_Pipeline/README
