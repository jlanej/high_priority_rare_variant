# Allele-Frequency Filtering & the Frequency Oracle

How this pipeline decides whether a variant is rare enough to be a high-priority candidate, using an external population reference rather than the untrustworthy internal cohort.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

> ### ⚠ Status: what of this document runs
>
> The pipeline runs a **VEP-only contract** — a VEP 115 GRCh38 cache plus the CADD plugin, and
> no other resource file. That determines which half of this document is code and which half is
> reference science:
>
> | Section | Status |
> |---|---|
> | Why an external oracle, not internal AC/AN | **IMPLEMENTED** — and structural; it is why gnomAD is read at all |
> | The gnomAD v4.1 reference dataset | **IMPLEMENTED**, but read from the **VEP cache**, not a sites VCF — so only the *point* AFs of it exist here |
> | Rarity field = **grpmax proxy** (max AF over AFR/AMR/EAS/NFE/SAS) | **IMPLEMENTED** (`src/hprv/annotations.py:frequency`) |
> | Rarity field = **grpmax `faf95`** | **TARGET, not implemented.** The cache carries no AC/AN, so the CI lower bound is *unrecoverable* — see [§ faf95 is unavailable](#faf95-is-unavailable-and-why-a-cache-cannot-supply-it) |
> | `nhomalt` conditions (de novo absent-or-singleton, homozygote sanity check) | **TARGET, not implemented** — no `nhomalt` field exists |
> | Whiffin/Ware maximum credible AF; ClinGen VCEP gene-specific BA1/BS1 | **TARGET, not implemented** — no gene-specific override table is wired |
> | Exome/genome discordance flag | **TARGET, not implemented** — not a cache field |
> | Why `MAX_AF` and global AF must never be substituted | **IMPLEMENTED and test-enforced** — [read this before "simplifying" the rarity field](#the-max_af-trap-the-most-dangerous-simplification-in-this-pipeline) |
>
> The reference science below is retained deliberately: it is the justification for the roadmap,
> not decoration. The full ledger of what the first pass cannot see lives in
> **[limitations.md](limitations.md)** (§2 covers faf95, §2a the MAX_AF trap, §3 nhomalt).

## TL;DR

- **Oracle = gnomAD v4.1 (GRCh38)**, read from the **VEP cache** (`--af_gnomade` / `--af_gnomadg`). External, never internal.
- **The rarity field is a grpmax _proxy_:** the max **point-estimate AF** across the grpmax-**eligible** ancestry groups only — `AFR, AMR, EAS, NFE, SAS`. One function owns it for the whole pipeline: `annotations.frequency()`.
- **It is not `faf95`.** faf95 (the 95% CI lower bound) is the *target* and remains the right field; it needs AC/AN, which the cache does not carry, so it is **absent — not approximated**. The proxy therefore runs ~one CI-width **high** on low-count alleles, i.e. the gates err toward **dropping**.
- **It is not `MAX_AF` and must never become `MAX_AF`.** MAX_AF maxes over the bottlenecked founder groups grpmax deliberately excludes; global AF fails in the *opposite* direction. [Neither is a safe fallback.](#the-max_af-trap-the-most-dangerous-simplification-in-this-pipeline)
- **Dominant / de novo:** keep if proxy AF `< 1e-4` (applied at Step 5, per mode).
- **Recessive / compound-het:** keep if proxy AF `< 1e-2` (permissive discovery default), with a `< 1e-3` **high-confidence tier** that *flags* (`high_conf_rarity`) rather than drops; applied **per variant**, not per gene.
- **Hard benign (all modes):** drop if proxy AF `≥ 0.05` (ClinGen BA1) — never rescue, not even by ClinVar P/LP.
- **Never** use internal cohort AC/AN as a population frequency: the non-joint per-trio design makes AN uninterpretable (absent genotype ≠ hom-ref). Internal recurrence is valid **only** as an artifact/blocklist signal.
- *Target:* a **gene-specific ClinGen VCEP** BA1/BS1 or a **Whiffin/Ware maximum credible AF** should override any generic cutoff. Neither is wired yet.
- The rarity gate is a **screening filter**, distinct from the ACMG **PM2** criterion (applied at *Supporting* strength only). Passing the gate is not the same as "PM2 met."

## Why an external frequency oracle

This pipeline screens **GMKF Kids First per-trio VCFs** that are GATK Genotype-Refinement output but are **not jointly genotyped across the cohort**. That design makes internal allele counts unusable as a population frequency:

- **No consistent cohort-wide AN.** Each trio is called independently, so a variant's internal frequency reflects only 2–6 chromosomes; there is no shared denominator across trios.
- **Absence ≠ reference.** In a non-joint merge, a variant missing from another trio may be a no-call or low-depth site, not a confident hom-ref. Internal AC/AN therefore mis-estimates both numerator and denominator.
- **gnomAD provides the defensible denominator.** It is large, uniformly joint-genotyped, ancestry-resolved, and ships proper filtering-allele-frequency confidence intervals. (This pipeline reaches gnomAD through the VEP cache, which relays the ancestry-resolved frequencies but *not* the confidence intervals — see [below](#faf95-is-unavailable-and-why-a-cache-cannot-supply-it). The argument for an external oracle is unaffected: a point estimate over ~807k uniformly genotyped samples is still categorically better than an AN of 6.)

Internal data still has one legitimate frequency-adjacent use: **artifact detection**. A variant recurring across many unrelated trios is more likely a systematic sequencing/mapping artifact than a truly common allele. Use that as a panel-of-normals-style **blocklist** signal (tune the recurrence count `N` empirically), never as a population AF. See [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md) and [cohort_construction.md](cohort_construction.md).

## The reference dataset: gnomAD v4.1 (GRCh38)

- **Composition:** 730,947 exomes (416,555 UK Biobank + 314,392 non-UKB) plus 76,215 genomes, all unrelated, aligned to GRCh38/hg38. v4.1 is the current release (Apr 2024). The union callset is ~807k samples, but exome vs genome N differs per site.
- **v4.1 key fixes:** corrects the v4.0 allele-number (AN) bug; adds a **joint (combined exome + genome) AN and AF** at every site called in either data type; adds a **discordant-frequency flag** where a contingency/CMH test between exomes and genomes gives p < 1e-4 (~2.5% of variants).
- **Practical rule (TARGET):** prefer the joint AF/AN and heed the discordance flag before trusting a single subset. Not implemented — joint AF/AN and the discordance flag are sites-VCF fields. The cache exposes exome and genome AFs separately, and `frequency()` simply takes the max across both, which is the conservative reading (it cannot under-call a frequency) but is **not** the joint estimate and knows nothing about exome/genome discordance.

### faf95 is unavailable, and why a cache cannot supply it

VEP release r113 (Oct 2024) updated its built-in gnomAD annotation to **v4.1** for both genomes
and exomes, and this pipeline reads it (`--af_gnomade` / `--af_gnomadg`). But be precise about
*what* it supplies, because the obvious assumption is wrong and expensive:

- The cache carries **point-estimate AFs only** — per-population and global. It has **no `faf95`,
  no `fafmax`, and no `AC`/`AN`.**
- faf95 is a **Poisson CI lower bound computed from AC and AN**. With neither numerator nor
  denominator, it cannot be recomputed downstream. **faf95 is not approximated here; it is
  unrecoverable at any price** — no amount of post-processing recovers a confidence interval from
  a point estimate that arrived with no counts attached.
- Restoring it means the **actual gnomAD data**: slim the 24 v4.1 joint chromosome VCFs down to
  `fafmax_faf95_max_joint` + `nhomalt_joint` (~10 GB kept, GCS egress free). See
  [limitations.md §2](limitations.md).
- Second cache caveat: cache frequencies exist only for alleles **accessioned into dbSNP**. An
  un-accessioned gnomAD variant silently returns *no* frequency and reads as "absent ⇒ rarest".
  Ensembl itself recommends `--custom` with the gnomAD VCF over `--af_gnomad*` for this reason.
  This biases toward **retention** (extra review), not toward missed calls.

The **exome/genome discordance flag** and the **joint AN** described above are likewise sites-VCF
fields, not cache fields — the "prefer joint AF/AN, heed the discordance flag" rule is a target
here, not a behaviour. See [tooling_and_reproducibility.md](tooling_and_reproducibility.md).

## Global AF vs grpmax vs FAF — use the right number

| Metric | What it is | Why we do / don't use it |
| --- | --- | --- |
| **Global AF** | AF across all samples | Dilutes an ancestry-enriched variant; a variant common in one group looks rare globally. **Do not filter on this.** Carried as `vep_gnomAD{e,g}_AF` for **reporting only**. |
| **VEP `MAX_AF`** | Max AF over *all* gnomAD groups **and** the 1000 Genomes phase-3 populations | **Never a filter field — it is a trap.** [See below.](#the-max_af-trap-the-most-dangerous-simplification-in-this-pipeline) Carried as `vep_MAX_AF` / `vep_MAX_AF_POPS` for reporting, so a reviewer can spot a call whose founder-group frequency is high. |
| **grpmax AF** (formerly popmax) | Highest point-estimate AF across the **grpmax-eligible** genetic-ancestry groups | Better than global, but a point estimate is noisy when a group's AN is small. **This — reconstructed as a proxy from the per-population cache AFs — is our current filter field.** |
| **FAF (faf95 / faf99)** | Lower bound of the 95% (or 99%) Poisson CI on the AF | The frequency you can be ≥95% confident the true AF is *at least*. Conservative for *filtering out* benign variants — you only exclude a variant as "too common" when confident it really is common. **The right filter field, and the TARGET — [unavailable under the VEP-only contract](#faf95-is-unavailable-and-why-a-cache-cannot-supply-it).** |
| **grpmax FAF** | faf95 from the ancestry group with the highest FAF | The value ClinGen VCEPs use for BA1/BS1. **TARGET.** |

**The proxy, precisely.** `annotations.frequency()` returns the max cache AF over
`GRPMAX_POPS = (AFR, AMR, EAS, NFE, SAS)` across both exome and genome fields. Mirroring gnomAD's
own grpmax *inclusion set* is exactly what makes it a defensible stand-in. Its one honest error:
a point estimate is always ≥ its own CI lower bound, so every gate fires slightly **more** often
than a faf95 gate would — the pipeline **errs toward dropping** low-count alleles (a
false-negative direction, bounded by AC, worst for singletons in the smaller eligible groups).

**Founder-group exclusion.** gnomAD excludes bottlenecked/founder groups (Amish, Ashkenazi Jewish,
Finnish, Middle Eastern, and "remaining") from grpmax FAF, because pathogenic founder alleles
legitimately reach high frequency there and would wrongly inflate the filter. The proxy reproduces
that exclusion by construction. Rely on it; do **not** re-introduce those groups' frequencies into
the gate. This also removes the *large* half of the point-estimate error above — a CI correction
matters most exactly where AN is small.

Using the CI lower bound is deliberately conservative: it protects against false exclusion of true pathogenic alleles that happen to appear by chance in a small sample. That protection is what the proxy currently lacks.

### The `MAX_AF` trap: the most dangerous simplification in this pipeline

`vep_MAX_AF` is right there in the CSQ, it is a single field instead of ten, and it is *labelled*
as the maximum population frequency. Substituting it for the grpmax proxy looks like an obvious
cleanup. **It is a regression that silently destroys real candidates**, and it is the single most
likely wrong "fix" a future maintainer will make. It is guarded by
`tests/test_pure.py::test_frequency_excludes_bottlenecked_pops` and by the `GENEFND` integration
case — if you find yourself deleting either, stop.

**Why it fails — the worked example.** MAX_AF maximises over the bottlenecked founder groups
gnomAD's own grpmax *deliberately excludes* (`ami`, `asj`, `fin`, `mid`, `remaining`) **and** over
the tiny 1000 Genomes phase-3 populations. The Amish subset has **AN ≈ 900**. So a **single
observed allele** in `ami` — one chromosome, in a founder population, carrying no information
whatsoever about the general population — reads as:

```
AF_ami = 1 / 900 ≈ 1.1e-3
```

That is **ten-fold over `dominant_max = 1e-4`**. A genuinely ultra-rare, absent-everywhere-else
dominant candidate is thrown away on the strength of one chromosome in a group gnomAD explicitly
tells you not to filter on. The failure is **silent**: no warning, no flag, just a variant that
never appears in the output. This is precisely the false-negative mode grpmax was invented to
prevent, and MAX_AF re-introduces all of it.

**Global AF fails in the opposite direction.** `gnomADe_AF` / `gnomADg_AF` average an
ancestry-enriched benign polymorphism across the whole cohort: a variant at 3% in AFR and absent
elsewhere dilutes to well under the gate and is **retained** as a false positive.

**So there is no single safe fallback.** The two available shortcuts err in **opposite
directions** — MAX_AF over-drops (kills true positives), global AF over-retains (floods review).
Neither can be swapped in "just to simplify"; the per-population max over the eligible groups is
the only field that is wrong in neither direction. Both are kept in the output as **reporting**
columns, next to `grpmax_af`, and are never read by a gate. See [limitations.md §2a](limitations.md).

## Inheritance-mode–dependent rarity gates

The maximum tolerated frequency depends on inheritance mode. These are **screening defaults**, overridable in `config/config.example.yaml`. The filter field is the **grpmax proxy** described above (`annotations.frequency()`) — read `faf95` in the literature citations, `grpmax proxy AF` in the code.

| Candidate class | Keep if (proxy AF) | Applied | Notes |
| --- | --- | --- | --- |
| **Dominant / de novo** | `< 1e-4` (`rarity.dominant_max`) | **Step 5**, per mode | *Target, absent:* the de novo "absent-or-singleton + low `nhomalt`" condition. No `nhomalt` field exists, and the retired config key implemented it as `nhomalt > 1` — a homozygote-count test, never the allele-count test its name promised. The `< 1e-4` gate still applies to de novo calls and does the bulk of the work. |
| **Recessive / compound-het** | `< 1e-2` (`rarity.recessive_max`) | **Step 3** (permissive union) **+ Step 5** | The `< 1e-3` tier (`rarity.recessive_strict`) **flags** the call `high_conf_rarity` — it does **not** drop. Applied per variant, not per gene. |
| **Hard benign (all modes)** | drop if `≥ 0.05` (`rarity.benign_ba1`) | **Step 3** | ClinGen general-purpose **BA1**; never rescued — not even by a ClinVar P/LP assertion. |

**Where the gates actually fire.** Step 3 (`selection.py`) is inheritance-agnostic, so it applies
a **permissive union**: BA1 drops, then anything at or above `recessive_max` drops as `too_common`
*unless* ClinVar P/LP rescues it. The mode-specific `dominant_max` is applied later, in Step 5,
once inheritance is known. A consequence worth knowing: a ClinVar P/LP variant at AF 2e-4 survives
Step 3 via the P/LP rescue but is **not** called dominant at Step 5 — Step 5's rarity check has no
ClinVar override. It reaches the candidate VCF, not the dominant call set.

The literature range for the generic recessive cutoff spans roughly 1e-3 to 1e-2 (5e-3 is a commonly cited midpoint); this pipeline uses the permissive 1e-2 discovery default with a 1e-3 high-confidence tier so that biallelic candidates are not lost early. For dominant conditions, ClinGen general practice sits near grpmax faf95 `< 1e-4` absent a gene-specific value — the pipeline adopts that number, applied to the proxy rather than to faf95.

**Gene-specific override — TARGET, not implemented.** No gene-specific BA1/BS1 table is wired into
the config; the generic cutoffs above apply uniformly to every gene today. Where a ClinGen VCEP
publishes calibrated BA1/BS1 values (e.g. cardiomyopathy, RASopathy), those *should* win. The
RASopathy VCEP, for instance, calibrates BA1 grpmax faf well below the generic 0.05 — so the
generic 0.05 is, for those genes, knowably too permissive. See [clinical_classification.md](clinical_classification.md).

## Maximum credible population allele frequency (Whiffin/Ware) — TARGET

> **Not implemented.** No per-gene maxAF table is wired into the config; nothing in the pipeline
> computes or consults one. This section is the reference science and the design for it — it is
> the most defensible available replacement for a generic cutoff, and it is retained because it
> is what the roadmap should build, not because it runs.

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

## How this is wired (IMPLEMENTED)

**Step 2** annotates the cohort union once with the cache's gnomAD v4.1 fields and lifts them to
INFO with `bcftools +split-vep -p vep_`:

```bash
vep \
  --offline --cache --dir_cache "${VEP_CACHE_DIR}" \
  --assembly GRCh38 --fasta "${REF_FASTA}" \
  --vcf --compress_output bgzip \
  --input_file "${IN_VCF}" --output_file "${OUT_VCF}" \
  --af_gnomade --af_gnomadg          # per-population + global point AFs. No faf95 exists.
```

There is **no** `bcftools annotate` transfer from a gnomAD sites VCF — that is the whole VEP-only
contract. Step 2 fails loudly if *none* of the ten grpmax-eligible AF fields
(`gnomAD{e,g}_{AFR,AMR,EAS,NFE,SAS}_AF`) survives into the CSQ, because a silently-absent rarity
oracle reads as "everything is rare" and the screen would keep every common polymorphism. It also
asserts on the **values**, not just the header: a cache built without frequency data yields a
fully-populated column of empty strings, which fails the same way.

**Filtering is not a `bcftools` expression.** Every rarity decision goes through one Python
chokepoint — `annotations.frequency()` — which Steps 3, 5 and 6 all read:

```python
GRPMAX_POPS = ("AFR", "AMR", "EAS", "NFE", "SAS")   # gnomAD's own grpmax inclusion set

def grpmax_af(variant):        # max over vep_gnomAD{e,g}_{AFR,AMR,EAS,NFE,SAS}_AF
    ...

def frequency(variant):        # THE rarity field. None => no eligible group reports it => rarest.
    return grpmax_af(variant)
```

Keeping it a single function is deliberate: it is the one place a maintainer could quietly swap in
`MAX_AF` and break the screen invisibly, so it is also the one place the tests watch. If you add a
frequency source, add it here — nothing else in the codebase reaches around this contract. See
[functional_annotation.md](functional_annotation.md) for the downstream consequence layer.

## Known limitations

The frequency-specific gaps — **no faf95** (§2), **the MAX_AF trap** (§2a), **no nhomalt** (§3) —
are documented once, in **[limitations.md](limitations.md)**, with the cost to fix each. Summarised
above rather than restated here. What is specific to this layer:

- **The rarity gate is a point estimate, so it errs toward dropping.** Every gate fires slightly
  more often than a faf95 gate would. Direction matters: this is a **false-negative** bias on
  low-count alleles, not a false-positive one. Bounded by AC, worst for singletons in the smaller
  eligible groups.
- **"Absent" is weaker evidence than it looks.** The cache only carries frequencies for alleles
  accessioned into dbSNP, so an un-accessioned gnomAD variant reads as absent ⇒ rarest. Biases
  toward retention (extra review), not toward missed calls — the opposite direction to the above.
- **SNV/indel only.** This frequency logic applies to short-variant calls. CNV/SV are a real blind spot (10–15% of pediatric-cancer and rare-disease diagnoses), and gnomAD SNV FAF does not address them; a future GATK-gCNV / Manta / ExomeDepth module is required. See [pipeline_design.md](pipeline_design.md).
- **Pseudogene / segmental-duplication genes** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) are low-confidence from short reads; gnomAD frequencies in those paralogous regions can themselves be unreliable, so flag those regions rather than trusting the AF. **Not currently flagged or masked.**
- **Proband mosaicism.** Low-VAF post-zygotic calls are handled by the genotype-QC layer, not here, but note that a genuinely rare pathogenic mosaic call must survive both the AB band and this rarity gate.
- **Founder/bottlenecked populations** are excluded from the rarity field by construction, mirroring gnomAD's own grpmax; a pathogenic founder allele common in an excluded group is intentionally not counted against the filter. This is a **feature, not a gap** — see the [MAX_AF trap](#the-max_af-trap-the-most-dangerous-simplification-in-this-pipeline).
- **No gene-specific overrides, no maximum-credible-AF model, no discordance flag.** All three are targets; a generic cutoff is applied to every gene.

## Recommended defaults (this pipeline)

| Parameter | Default | Status | Notes |
| --- | --- | --- | --- |
| Frequency oracle | gnomAD **v4.1** (GRCh38), from the **VEP cache** | **IMPLEMENTED** | Point AFs only. Joint AF/AN and the exome/genome discordance flag are sites-VCF fields — *target*. |
| Filter field | **grpmax proxy** = max AF over `AFR/AMR/EAS/NFE/SAS` | **IMPLEMENTED** | `annotations.frequency()`. A point estimate. **Never** `MAX_AF` (over-drops) and **never** global AF (over-retains). |
| Filter field | grpmax **faf95** | **TARGET** | The right field. Needs AC/AN ⇒ unrecoverable from the cache; needs the gnomAD sites VCF (~10 GB slim). |
| Dominant / de novo keep | proxy AF `< 1e-4` | **IMPLEMENTED** (Step 5) | De novo "absent-or-singleton + low `nhomalt`" is **removed** — no `nhomalt` field exists. |
| Recessive / comp-het keep | proxy AF `< 1e-2` (discovery) | **IMPLEMENTED** (Steps 3 + 5) | Per variant, not per gene. |
| High-confidence recessive tier | proxy AF `< 1e-3` | **IMPLEMENTED** (Step 5) | **Flags** `high_conf_rarity`; does not drop. |
| Hard benign (all modes) | proxy AF `≥ 0.05` → drop | **IMPLEMENTED** (Step 3) | ClinGen BA1; never rescued, not even by ClinVar P/LP. |
| Internal cohort AC/AN | **Never** as population AF | **IMPLEMENTED** (structural) | Use only as recurrence/artifact blocklist signal. |
| PM2 | Supporting strength only | **TARGET** | No ACMG criterion assignment step exists yet. Evidence, not the screening gate. |
| Gene-specific override | ClinGen VCEP BA1/BS1 or Whiffin/Ware maxAF | **TARGET** | Nothing is wired; the generic cutoff applies to every gene. |

All values are configurable defaults, not immutable law. Prefer a gene-specific ClinGen VCEP value over any generic cutoff — once there is machinery to express one.

## Sources

- gnomAD v4.1 release notes: https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/
- gnomAD FAF help: https://gnomad.broadinstitute.org/help/faf
- ClinGen guidance to VCEPs on using gnomAD v4 (BA1/BS1, grpmax FAF), March 2024: https://clinicalgenome.org/site/assets/files/9445/clingen_guidance_to_vceps_regarding_the_use_of_gnomad_v4_march_2024.pdf
- ClinGen SVI PM2 recommendation v1.0 (2020; Supporting default): https://clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf
- Whiffin et al., "Using high-resolution variant frequencies to empower clinical genome interpretation," Genet Med 2017: https://www.nature.com/articles/gim201726.pdf (PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5563454)
- Maximum credible allele frequency calculator (CardioDB): https://www.cardiodb.org/allelefrequencyapp/
- ClinGen RASopathy ACMG specification (calibrated grpmax faf thresholds): https://www.sciencedirect.com/science/article/pii/S2949774425014694
- ACGS 2024 UK variant-classification best-practice guidelines: https://www.genomicseducation.hee.nhs.uk/wp-content/uploads/2024/08/ACGS-2024_UK-practice-guidelines-for-variant-classification.pdf
