# Allele-Frequency Filtering & the Frequency Oracle

How this pipeline decides whether a variant is rare enough to be a high-priority candidate, using an external population reference rather than the untrustworthy internal cohort.

> Part of the high_priority_rare_variant methods reference. Thresholds here are the
> configurable defaults defined in [Canonical defaults](README.md#canonical-defaults).

> ### ⚠ Status: what of this document runs
>
> The pipeline runs a **VEP-centric contract** — a VEP 115 GRCh38 cache plus the CADD / SpliceAI /
> REVEL / AlphaMissense plugins and exactly two `bcftools annotate` transfers (the ClinVar sites VCF
> and the gnomAD v4.1 joint slim). That determines which half of this document is code and which
> half is reference science:
>
> | Section | Status |
> |---|---|
> | Why an external oracle, not internal AC/AN | **IMPLEMENTED** — and structural; it is why gnomAD is read at all |
> | The gnomAD v4.1 reference dataset | **IMPLEMENTED** — read from the gnomAD v4.1 **joint slim** (faf95, the default oracle) and from the **VEP cache** (per-population point AFs, the opt-down arm) |
> | Rarity field = **grpmax `faf95`** | **IMPLEMENTED and the DEFAULT** (`resources.gnomad.oracle: faf95`), via the gnomAD joint slim transferred in Step 2 (`resources.gnomad.sites_slim`, required; the run halts at preflight without it) — see [§ one oracle per run](#one-oracle-per-run-faf95-default-or-the-grpmax-proxy) |
> | Rarity field = **grpmax proxy** (max point AF over AFR/AMR/EAS/NFE/SAS) | **IMPLEMENTED as a deliberate opt-down** (`oracle: grpmax_proxy`; `src/hprv/annotations.py:grpmax_af`). ONE oracle per run — the arms never cross |
> | `nhomalt` conditions | **PARTIAL.** `nhomalt` is transferred with the slim and reported per variant, and a biallelic call in a gene where gnomAD carries homozygotes raises `nhomalt_recessive_conflict`. The flag costs **0 points by default** — there is no calibration for how many homozygotes should disqualify a recessive candidate, so hprv reports it rather than inventing a penalty. |
> | Whiffin/Ware maximum credible AF; ClinGen VCEP gene-specific BA1/BS1 | **TARGET, not implemented** — no gene-specific override table is wired |
> | Exome/genome discordance flag | **TARGET, not implemented** — not a cache field |
> | Why `MAX_AF` and global AF must never be substituted | **IMPLEMENTED and test-enforced** — [read this before "simplifying" the rarity field](#the-max_af-trap-the-most-dangerous-simplification-in-this-pipeline) |
>
> The reference science below is retained deliberately: it is the justification for the roadmap,
> not decoration. The full ledger of what the first pass cannot see lives in
> **[limitations.md](limitations.md)** (§2 covers faf95, §2a the MAX_AF trap, §3 nhomalt).

## TL;DR

- **Oracle = gnomAD v4.1 (GRCh38).** External, never internal. Reached through one function for the whole pipeline, `annotations.frequency()`, which reads **ONE quantity per run** (`resources.gnomad.oracle`).
- **`faf95` is the DEFAULT.** gnomAD's published filtering allele frequency (`fafmax_faf95_max_joint`, the 95% CI lower bound), transferred in Step 2 from the gnomAD v4.1 joint slim as `gnomad_faf95`. Requires `resources.gnomad.sites_slim`; the run halts at preflight without it. An allele gnomAD has but published no faf95 for resolves to **0** (`rarity_basis=zero_ci`, rarest); an allele with no gnomAD record resolves to **absent** (rarest). The proxy is never consulted on this arm.
- **`grpmax_proxy` is a deliberate opt-down:** the max **point-estimate AF** across the grpmax-**eligible** ancestry groups only — `AFR, AMR, EAS, NFE, SAS` — read from the VEP cache. It runs ~one CI-width **high** on low-count alleles, i.e. it errs toward **dropping**; faf95 removes that error, so the default **retains more** at the same cutoffs. `rarity_oracle` is recorded once per run; `rarity_basis` (`measured`/`zero_ci`/`absent`) is the per-variant provenance within that oracle.
- **It is not `MAX_AF` and must never become `MAX_AF`, under either arm.** MAX_AF maxes over the bottlenecked founder groups grpmax deliberately excludes; global AF fails in the *opposite* direction. [Neither is a safe fallback.](#the-max_af-trap-the-most-dangerous-simplification-in-this-pipeline)
- **Dominant / de novo:** keep if `rarity_af` `< 1e-4` (applied at Step 5, per mode).
- **Recessive / compound-het:** keep if `rarity_af` `< 1e-2` (permissive discovery default), with a `< 1e-3` **high-confidence tier** that *flags* (`high_conf_rarity`) rather than drops; applied **per variant**, not per gene.
- **Hard benign (all modes):** drop if `rarity_af` `≥ 0.05` (ClinGen BA1) — never rescue, not even by ClinVar P/LP.
- **Never** use internal cohort AC/AN as a population frequency: the non-joint per-trio design makes AN uninterpretable (absent genotype ≠ hom-ref). Internal recurrence is valid **only** as an artifact/blocklist signal.
- *Target:* a **gene-specific ClinGen VCEP** BA1/BS1 or a **Whiffin/Ware maximum credible AF** should override any generic cutoff. Neither is wired yet.
- The rarity gate is a **screening filter**, distinct from the ACMG **PM2** criterion (applied at *Supporting* strength only). Passing the gate is not the same as "PM2 met."

## Why an external frequency oracle

This pipeline screens **GMKF Kids First per-trio VCFs** that are GATK Genotype-Refinement output but are **not jointly genotyped across the cohort**. That design makes internal allele counts unusable as a population frequency:

- **No consistent cohort-wide AN.** Each trio is called independently, so a variant's internal frequency reflects only 2–6 chromosomes; there is no shared denominator across trios.
- **Absence ≠ reference.** In a non-joint merge, a variant missing from another trio may be a no-call or low-depth site, not a confident hom-ref. Internal AC/AN therefore mis-estimates both numerator and denominator.
- **gnomAD provides the defensible denominator.** It is large, uniformly joint-genotyped, ancestry-resolved, and ships proper filtering-allele-frequency confidence intervals. (By default this pipeline reads those intervals directly — `faf95` from the gnomAD joint slim; the opt-down arm reaches gnomAD through the VEP cache, which relays the ancestry-resolved point frequencies but *not* the confidence intervals — see [below](#one-oracle-per-run-faf95-default-or-the-grpmax-proxy). The argument for an external oracle is unaffected either way: even a point estimate over ~807k uniformly genotyped samples is categorically better than an AN of 6.)

Internal data still has one legitimate frequency-adjacent use: **artifact detection**. A variant recurring across many unrelated trios is more likely a systematic sequencing/mapping artifact than a truly common allele. Use that as a panel-of-normals-style **blocklist** signal (tune the recurrence count `N` empirically), never as a population AF. See [inheritance_and_genotype_qc.md](inheritance_and_genotype_qc.md) and [cohort_construction.md](cohort_construction.md).

## The reference dataset: gnomAD v4.1 (GRCh38)

- **Composition:** 730,947 exomes (416,555 UK Biobank + 314,392 non-UKB) plus 76,215 genomes, all unrelated, aligned to GRCh38/hg38. v4.1 is the current release (Apr 2024). The union callset is ~807k samples, but exome vs genome N differs per site.
- **v4.1 key fixes:** corrects the v4.0 allele-number (AN) bug; adds a **joint (combined exome + genome) AN and AF** at every site called in either data type; adds a **discordant-frequency flag** where a contingency/CMH test between exomes and genomes gives p < 1e-4 (~2.5% of variants).
- **Practical rule:** prefer the joint AF/AN. **Satisfied by the default oracle**: `faf95` is read from the gnomAD **joint** release (`fafmax_faf95_max_joint`), so the default oracle *is* the joint estimate. The opt-down proxy is a max over the cache's separate exome/genome per-population AFs — conservative (it cannot under-call) but not the joint estimate. The exome/genome **discordance flag** remains unwired (it is a sites-VCF field the slim does not keep).

### One oracle per run: faf95 (default) or the grpmax proxy

VEP release r113 (Oct 2024) updated its built-in gnomAD annotation to **v4.1** for both genomes
and exomes, and this pipeline reads it (`--af_gnomade` / `--af_gnomadg`). But be precise about
*what* it supplies, because the obvious assumption is wrong and expensive:

- The cache carries **point-estimate AFs only** — per-population and global. It has **no `faf95`,
  no `fafmax`, and no `AC`/`AN`.**
- faf95 is a **Poisson CI lower bound computed from AC and AN**. With neither numerator nor
  denominator, it cannot be recomputed downstream. **faf95 is not approximated here; it is
  unrecoverable at any price** — no amount of post-processing recovers a confidence interval from
  a point estimate that arrived with no counts attached.
- **So faf95 comes from the actual gnomAD data**, and now does: `prepare_resources.sh --only
  gnomad_sites fetch` stream-slims the 24 v4.1 **joint** chromosome VCFs to 5 of their ~664 INFO
  fields (~10 GB lands, GCS egress is free, nothing raw is staged), and Step 2 `bcftools
  annotate`-transfers them as `gnomad_faf95` / `gnomad_faf95_group` / `gnomad_nhomalt` (+ two
  reporting-only AFs). This is the second and last external transfer in the pipeline.

**Which groups faf95 covers — verified against the real v4.1 joint data, not assumed.** The FAF
group set is **afr, amr, eas, mid, nfe, sas**: `GRPMAX_POPS` **plus `mid`**, and critically it
EXCLUDES the bottlenecked ami/asj/fin that [the MAX_AF trap](#the-max_af-trap-the-most-dangerous-simplification-in-this-pipeline)
is about. So faf95 does **not** reintroduce that failure mode. `mid` is the single deviation from
hprv's own proxy, which is why the producing group rides along as `faf95_group` on every row —
the same reason `max_af_pops` rides beside `max_af`.

**ONE oracle per run — the arms never cross.** `resources.gnomad.oracle` picks the quantity for
the whole run — **`faf95` by default**, because it is the quantity ACMG/ClinGen specify and the
one you can cite without caveat. It requires the gnomAD joint slim and **halts at preflight**
without it rather than quietly running on the point estimate; `grpmax_proxy` is the deliberate
opt-down. The choice is recorded once in `audit/counts.tsv`, and `annotations.frequency()` never
consults the other arm. An earlier
design preferred faf95 per variant and fell back to the proxy, which made two rows in one run
comparable on different quantities — undescribable in a methods section — and got the fallback
direction wrong (below). `rarity_basis` records the per-variant provenance *within* the chosen
oracle: `measured` | `zero_ci` | `absent`.

**Both arms are gnomAD v4.1.** The proxy is the VEP cache's own gnomAD AFs, so this is a choice of
QUANTITY (CI lower bound vs point estimate), never of database. There is no such thing as "absent
from gnomAD but has a proxy".

**`zero_ci`: an absent faf95 on an allele gnomAD HAS is 0, not unknown.** gnomAD emits `fafmax` as
**missing**, never as `0`, wherever no group's CI lower bound clears zero — 80% of a chr22 sample.
Of the records with no faf95 but a proxy ≥ 1e-4, **96.5% are AC ≤ 2**: singletons whose point
estimate is inflated by a small group's AN (AC=1 / AN≈4,500 reads as 2.2e-4). Resolving those to 0
(rarest) is the whole point of a filtering allele frequency. `gnomad_AF_joint` is the witness that
distinguishes "gnomAD looked and could not bound it" from "gnomAD has no record".

**The two arms differ on `mid`, and this is the one behavioural difference to state in a methods
section.** The grpmax proxy excludes `mid` (mirroring gnomAD's own grpmax); the FAF group set
**includes** it. So an allele enriched only in the Middle Eastern group is invisible to the proxy
and *is* seen by faf95 — the arms can reach opposite conclusions on exactly that class.
`faf95_group` reports the producing group so those rows are identifiable.

**Expect a LARGER candidate list.** faf95 ≤ the point estimate, so the same cutoffs stop
discarding low-count alleles whose confidence interval never justified the call. If supplying the
slim made your list *smaller*, something is wrong — check the Step-2 join count.
- Cache caveat (proxy arm only): cache frequencies exist only for alleles **accessioned into
  dbSNP**. An un-accessioned gnomAD variant silently returns *no* frequency and reads as "absent ⇒
  rarest". Ensembl itself recommends `--custom` with the gnomAD VCF over `--af_gnomad*` for this
  reason. This biases toward **retention** (extra review), not toward missed calls. The joint slim
  carries every gnomAD allele, so the default arm does not share this gap.

The **exome/genome discordance flag** and the **joint AN** described above are likewise sites-VCF
fields, not cache fields — the "prefer joint AF/AN, heed the discordance flag" rule is a target
here, not a behaviour. See [tooling_and_reproducibility.md](tooling_and_reproducibility.md).

## Global AF vs grpmax vs FAF — use the right number

| Metric | What it is | Why we do / don't use it |
| --- | --- | --- |
| **Global AF** | AF across all samples | Dilutes an ancestry-enriched variant; a variant common in one group looks rare globally. **Do not filter on this.** Carried as `vep_gnomAD{e,g}_AF` for **reporting only**. |
| **VEP `MAX_AF`** | Max AF over *all* gnomAD groups **and** the 1000 Genomes phase-3 populations | **Never a filter field — it is a trap.** [See below.](#the-max_af-trap-the-most-dangerous-simplification-in-this-pipeline) Carried as `vep_MAX_AF` / `vep_MAX_AF_POPS` for reporting, so a reviewer can spot a call whose founder-group frequency is high. |
| **grpmax AF** (formerly popmax) | Highest point-estimate AF across the **grpmax-eligible** genetic-ancestry groups | Better than global, but a point estimate is noisy when a group's AN is small. **Reconstructed as a proxy from the per-population cache AFs, this is the `grpmax_proxy` arm — the filter field only when a run opts down from faf95.** |
| **FAF (faf95 / faf99)** | Lower bound of the 95% (or 99%) Poisson CI on the AF | The frequency you can be ≥95% confident the true AF is *at least*. Conservative for *filtering out* benign variants — you only exclude a variant as "too common" when confident it really is common. **The right filter field, and the DEFAULT** (`gnomad_faf95` from the joint slim). |
| **grpmax FAF** | faf95 from the ancestry group with the highest FAF | The value ClinGen VCEPs use for BA1/BS1. This is exactly `fafmax_faf95_max_joint`, i.e. what hprv reads. **IMPLEMENTED.** |

**The proxy, precisely.** Under `oracle: grpmax_proxy`, `annotations.frequency()` returns
`grpmax_af()`: the max cache AF over `GRPMAX_POPS = (AFR, AMR, EAS, NFE, SAS)` across both exome
and genome fields. Mirroring gnomAD's own grpmax *inclusion set* is exactly what makes it a
defensible stand-in. Its one honest error: a point estimate is always ≥ its own CI lower bound, so
every gate fires slightly **more** often than a faf95 gate would — that arm **errs toward
dropping** low-count alleles (a false-negative direction, bounded by AC, worst for singletons in
the smaller eligible groups).

**Founder-group exclusion.** gnomAD excludes bottlenecked/founder groups (Amish, Ashkenazi Jewish,
Finnish, Middle Eastern, and "remaining") from grpmax FAF, because pathogenic founder alleles
legitimately reach high frequency there and would wrongly inflate the filter. The proxy reproduces
that exclusion by construction. Rely on it; do **not** re-introduce those groups' frequencies into
the gate. This also removes the *large* half of the point-estimate error above — a CI correction
matters most exactly where AN is small.

Using the CI lower bound is deliberately conservative: it protects against false exclusion of true pathogenic alleles that happen to appear by chance in a small sample. That protection is what the opt-down proxy lacks, and why faf95 is the default.

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

The maximum tolerated frequency depends on inheritance mode. These are **screening defaults**, overridable in `config/config.example.yaml`. The filter field is `rarity_af` — whatever `annotations.frequency()` returned under the run's oracle: real `faf95` by default, the grpmax proxy only if a run opts down.

| Candidate class | Keep if (`rarity_af`) | Applied | Notes |
| --- | --- | --- | --- |
| **Dominant / de novo** | `< 1e-4` (`rarity.dominant_max`) | **Step 5**, per mode | The de novo "absent-or-singleton + low `nhomalt`" condition stays **retired by choice**: the retired config key implemented it as `nhomalt > 1` — a homozygote-count test, never the allele-count test its name promised. `nhomalt` IS transferred with the slim and is reported, not gated (Step 9's `nhomalt_recessive_conflict`, 0 points by default). The `< 1e-4` gate does the work. |
| **Recessive / compound-het** | `< 1e-2` (`rarity.recessive_max`) | **Step 3** (permissive union) **+ Step 5** | The `< 1e-3` tier (`rarity.recessive_strict`) **flags** the call `high_conf_rarity` — it does **not** drop. Applied per variant, not per gene. |
| **Hard benign (all modes)** | drop if `≥ 0.05` (`rarity.benign_ba1`) | **Step 3** | ClinGen general-purpose **BA1**; never rescued — not even by a ClinVar P/LP assertion. |

**Where the gates actually fire.** Step 3 (`selection.py`) is inheritance-agnostic, so it applies
a **permissive union**: BA1 drops, then anything at or above `recessive_max` drops as `too_common`
*unless* ClinVar P/LP rescues it. The mode-specific `dominant_max` is applied later, in Step 5,
once inheritance is known. A consequence worth knowing: a ClinVar P/LP variant at AF 2e-4 survives
Step 3 via the P/LP rescue but is **not** called dominant at Step 5 — Step 5's rarity check has no
ClinVar override. It reaches the candidate VCF, not the dominant call set.

The literature range for the generic recessive cutoff spans roughly 1e-3 to 1e-2 (5e-3 is a commonly cited midpoint); this pipeline uses the permissive 1e-2 discovery default with a 1e-3 high-confidence tier so that biallelic candidates are not lost early. For dominant conditions, ClinGen general practice sits near grpmax faf95 `< 1e-4` absent a gene-specific value — the pipeline adopts that number, applied to faf95 by default (to the proxy only under `oracle: grpmax_proxy`).

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
  --af_gnomade --af_gnomadg          # per-population + global point AFs: the grpmax_proxy arm
                                     # (reporting columns under faf95). faf95 does not come from
                                     # here — see the joint-slim transfer below.
```

Then the second of Step 2's two `bcftools annotate` transfers brings in the joint slim (the first
is the ClinVar sites VCF):

```bash
bcftools annotate -a "${GNOMAD_SITES_SLIM}" \
  -c INFO/gnomad_faf95:=INFO/fafmax_faf95_max_joint,INFO/gnomad_faf95_group:=INFO/fafmax_faf95_max_gen_anc_joint,INFO/gnomad_nhomalt:=INFO/nhomalt_joint,INFO/gnomad_AF_joint:=INFO/AF_joint,INFO/gnomad_AF_grpmax:=INFO/AF_grpmax_joint \
  -Oz -o "${ANNOTATED_VCF}" "${SPLIT_VCF}"
```

Step 2 logs `gnomAD joint matched N / M sites (faf95 present on K)` — the denominator is sites with
ANY gnomAD INFO, not sites with faf95, because an absent faf95 is legitimate (`zero_ci`) — and dies
on a 0-match join. Under the default `oracle: faf95` a missing or failed transfer is a hard stop,
including on the `resources.vep.annotated_vcf` ingest path; Step 3 additionally asserts
`gnomad_AF_joint` is declared in the header. Step 2 also fails loudly if *none* of the ten
grpmax-eligible AF fields (`gnomAD{e,g}_{AFR,AMR,EAS,NFE,SAS}_AF`) survives into the CSQ, because a
silently-absent proxy arm reads as "everything is rare" and an opted-down run would keep every
common polymorphism. It asserts on the **values**, not just the header: a cache built without
frequency data yields a fully-populated column of empty strings, which fails the same way.

**Filtering is not a `bcftools` expression.** Every rarity decision goes through one Python
chokepoint — `annotations.frequency()` — which Steps 3, 5 and 6 all read:

```python
GRPMAX_POPS = ("AFR", "AMR", "EAS", "NFE", "SAS")   # gnomAD's own grpmax inclusion set

def grpmax_af(variant):        # max over vep_gnomAD{e,g}_{AFR,AMR,EAS,NFE,SAS}_AF (point estimate)
    ...

def faf95(variant):            # gnomad_faf95 from the joint slim (95% CI lower bound), or None
    ...

def gnomad_observed(variant):  # True when the slim carries ANY record (gnomad_AF_joint present)
    ...

def rarity_oracle(cfg):        # "faf95" (default) | "grpmax_proxy" — ONE per run, from config
    ...

def frequency(variant, cfg):   # THE rarity field. None => the chosen oracle has no value => rarest.
    if rarity_oracle(cfg) == "faf95":
        f = faf95(variant)
        if f is not None:
            return f                                      # rarity_basis = measured
        return 0.0 if gnomad_observed(variant) else None  # zero_ci (rarest) | absent (rarest)
    return grpmax_af(variant)                             # the opt-down arm; never mixed with faf95
```

Keeping it a single function is deliberate: it is the one place a maintainer could quietly swap in
`MAX_AF`, or let the arms cross, and break the screen invisibly — so it is also the one place the
tests watch. If you add a frequency source, add it here — nothing else in the codebase reaches
around this contract. See [functional_annotation.md](functional_annotation.md) for the downstream
consequence layer.

## Known limitations

The frequency-specific entries — **faf95** (§2, the default oracle via the joint slim), **the MAX_AF
trap** (§2a, still live under either arm), **nhomalt** (§3, reported via the same slim) —
are documented once, in **[limitations.md](limitations.md)**, with the cost to fix each. Summarised
above rather than restated here. What is specific to this layer:

- **Under `oracle: grpmax_proxy` the rarity gate is a point estimate, so it errs toward
  dropping.** Every gate fires slightly more often than a faf95 gate would. Direction matters: this
  is a **false-negative** bias on low-count alleles, not a false-positive one. Bounded by AC, worst
  for singletons in the smaller eligible groups. The default `faf95` oracle does not have it.
- **"Absent" is weaker evidence than it looks — on the proxy arm.** The cache only carries
  frequencies for alleles accessioned into dbSNP, so an un-accessioned gnomAD variant reads as
  absent ⇒ rarest. Biases toward retention (extra review), not toward missed calls — the opposite
  direction to the above. The joint slim carries every gnomAD allele, so under `faf95` "absent"
  really does mean "not in gnomAD".
- **`mid`.** The FAF group set includes the Middle Eastern group and the proxy excludes it, so an
  allele enriched only in `mid` can be gated by one arm and not the other; `faf95_group` names the
  producing group on every row.
- **SNV/indel only.** This frequency logic applies to short-variant calls. CNV/SV are a real blind spot (10–15% of pediatric-cancer and rare-disease diagnoses), and gnomAD SNV FAF does not address them; a future GATK-gCNV / Manta / ExomeDepth module is required. See [pipeline_design.md](pipeline_design.md).
- **Pseudogene / segmental-duplication genes** (PMS2/PMS2CL, CYP21A2, SMN1/2, NEB, GBA) are low-confidence from short reads; gnomAD frequencies in those paralogous regions can themselves be unreliable, so flag those regions rather than trusting the AF. **Not currently flagged or masked.**
- **Proband mosaicism.** Low-VAF post-zygotic calls are handled by the genotype-QC layer, not here, but note that a genuinely rare pathogenic mosaic call must survive both the AB band and this rarity gate.
- **Founder/bottlenecked populations** are excluded from the rarity field by construction, mirroring gnomAD's own grpmax; a pathogenic founder allele common in an excluded group is intentionally not counted against the filter. This is a **feature, not a gap** — see the [MAX_AF trap](#the-max_af-trap-the-most-dangerous-simplification-in-this-pipeline).
- **No gene-specific overrides, no maximum-credible-AF model, no discordance flag.** All three are targets; a generic cutoff is applied to every gene.

## Recommended defaults (this pipeline)

| Parameter | Default | Status | Notes |
| --- | --- | --- | --- |
| Frequency oracle | gnomAD **v4.1** (GRCh38), ONE quantity per run (`resources.gnomad.oracle`) | **IMPLEMENTED** | `faf95` (default) from the joint slim; `grpmax_proxy` (opt-down) from the VEP cache. The exome/genome discordance flag is a sites-VCF field the slim does not keep — *target*. |
| Filter field (default) | grpmax **faf95** = `fafmax_faf95_max_joint` | **IMPLEMENTED, DEFAULT** | `gnomad_faf95` via the ~10 GB joint slim, transferred in Step 2; the run halts at preflight without it. FAF groups = afr/amr/eas/**mid**/nfe/sas — verified, and excluding ami/asj/fin. Absent-but-observed ⇒ 0 (`zero_ci`); no record ⇒ absent. The proxy is never consulted. |
| Filter field (opt-down) | **grpmax proxy** = max point AF over `AFR/AMR/EAS/NFE/SAS` | **IMPLEMENTED** (`oracle: grpmax_proxy`) | `annotations.grpmax_af()`. A point estimate; errs toward dropping. **Never** `MAX_AF` (over-drops) and **never** global AF (over-retains) — under either arm. |
| Dominant / de novo keep | `rarity_af` `< 1e-4` | **IMPLEMENTED** (Step 5) | De novo "absent-or-singleton + low `nhomalt`" stays **retired**; `nhomalt` is reported (Step 9 flag), not gated. |
| Recessive / comp-het keep | `rarity_af` `< 1e-2` (discovery) | **IMPLEMENTED** (Steps 3 + 5) | Per variant, not per gene. |
| High-confidence recessive tier | `rarity_af` `< 1e-3` | **IMPLEMENTED** (Step 5) | **Flags** `high_conf_rarity`; does not drop. |
| Hard benign (all modes) | `rarity_af` `≥ 0.05` → drop | **IMPLEMENTED** (Step 3) | ClinGen BA1; never rescued, not even by ClinVar P/LP. |
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
