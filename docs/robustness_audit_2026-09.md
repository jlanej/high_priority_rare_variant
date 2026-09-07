# Robustness audit (2026-09) — the plausibility filter and the genotype-to-meaning assignment

An audit of the two function families that decide what this pipeline can find at all, commissioned
because everything downstream of them is never-drop. Conventions follow
`docs/pipeline_review_2026-07.md` and `docs/pipeline_review_2026-09.md`: the code is ground truth,
every finding carries a file:line, and a finding demonstrated by running code is marked apart from
one argued from the code. Nothing in the repository was modified; this is a consolidation pass.

**How to read the severities.** Every finding carries two: the one the finder assigned, and the
median of three independent verifiers who re-examined it. They disagree a lot — the finders called
39 findings major, the verifiers' median calls **4** major — and the
verifiers are the better estimate, because each one had to either reproduce the failure or argue it
away. The document is ordered by the verifiers' severity; the finder's is shown in parentheses where
they differ. The verifiers overwhelmingly confirmed the *mechanisms* and downgraded the *stakes*:
"real, reproduced, and rare in GMKF-shaped data" is the modal verdict.

Method: three mapping agents established the drop inventory, the funnel arithmetic and the
structural premise; fourteen independent lenses swept the two families; a completeness critic named
five angles no lens had covered, and six follow-up finders investigated them; every finding was then
put to three verifiers — one that tried to reproduce it by execution, one that tried to overturn it
from the code, one that checked whether it was already documented or already fixed on this branch.

## The one-paragraph answer

The premise holds, with one correction that matters. Step 9 asserts row-count conservation twice
and Step 6 only aggregates, so **no step after the inheritance screen can lose a call** — verified
by reconciliation on a real run (41 = 41 = 41) and by a deliberately broken manifest that still
produced 41 rows. But never-drop conserves *the file Step 9 was handed*, not *the run's call set*,
and the audit found three ways those diverge. Meanwhile the two families that DO decide sensitivity
are the two least instrumented steps in the pipeline: Steps 1-4 report a reason for every drop
(Step 3 reports one for every keep as well), and Steps 5 and 6 have **no input counter at all**.
51 of the 106 findings are silent losses — a variant removed with no counter, no
flag and no reason string, which a reviewer cannot distinguish from "we looked and there was nothing
there" — and the four findings the verifiers left at major are all of that kind or worse: two ways
a trios file silently swaps the parents for a whole trio, the transmitting-parent QC asymmetry, and
the missing Step-5 instrumentation itself.

## The structural finding: never-drop is satisfied on the wrong denominator

What never-drop guarantees: given file *F*, Step 9 emits one row per row of *F*, asserted at
`pipeline/09_prioritize.py:942` and `:992`. (Worth knowing: the Step-9 scoring loop contains no
`continue` or `break`, so those assertions are structurally unreachable today. They defend future
edits, not present behaviour.)

What it does not guarantee, in descending order of practical danger:

1. **That *F* is the run's call set.** `pipeline/run_pipeline.sh:368` selects Step 9's input on
   non-emptiness alone, with no freshness check against `candidates.calls.tsv`. Demonstrated: a
   stale 10-row `igv/variants.tsv` beside a fresh 41-row `candidates.calls.tsv` made Step 9 process
   10 calls, pass both never-drop assertions, and exit 0 with no warning. Reachable via
   `outputs.igv.enabled: false`, a `--from 9` re-run, or SLURM array mode, which exits before
   assembling the table. `hprv_summary.xlsx` reads `candidates.calls.tsv` and would show 41.
   Nothing compares the two.
2. **That the call reaches the gene layer.** `pipeline/06_gene_burden.py:196` and
   `pipeline/09_prioritize.py:525` both `continue` on a call with no gene and no symbol, neither
   counted. Reachable: an intergenic variant kept by Step 3's CADD rung has no gene, so it
   vanishes from the gene layer while surviving in the variant layer with every gene column blank
   and no marker.
3. **That the call is reviewable.** A `trio_id` absent from the resolved manifest blanks every
   CRAM and VCF track path for every affected row. Demonstrated: 35 of 41 calls opened in the
   review server with no reads for any family member, while never-drop held at 41 and the audit
   read `08_igv variants 41`.
4. **That anything reached Step 5's output at all.** Step 5 has no input counter, and compound-het
   row inflation makes the audit funnel read as a *gain* where there was a loss. On the shipped
   integration fixture the audit says 31 candidate genotypes in and 35 calls out for CH_A — so no
   arithmetic a reviewer can do reveals that **4 of those 31 produced no row**, a 13% loss in the
   repo's own regression harness.

## The drop ledger

Every place a variant can be removed, and whether the removal leaves a trace. Measured on a real
run unless marked otherwise.

| Step | Guard | Removes | Trace? |
|---|---|---|---|
| resolve | no VCF holds all three members | whole trio | **yes** — `trios_unresolved` + `trio_resolution.tsv` |
| 1 | `view -f 'PASS,.'` | non-PASS records | **NO** — `input_sites` counts the POST-filter file |
| 1 | `-t "^chrM,chrMT,M,MT"` | all chrM | **NO** — `excluded_contigs_remaining` counts SURVIVORS, always 0 by construction |
| 1 | `view --min-ac 1` | AC=0 alleles; also a genuine carrier when INFO/AC is stale on the no-samples branch | **NO** |
| 1 | `norm -c e` | REF mismatch | dies loudly (correct) |
| 2 | `vep` without `--dont_skip` | records on contigs the cache lacks (`chrUn_*`, `HLA-*`, decoys) | partial — the delta is computable, nothing computes it |
| 2 | `gather_shards` | a shard where VEP dropped every record passes as a valid 0-record bgzip | **NO** |
| 3 | `fr >= benign_ba1` | BA1-common, incl. ClinVar P/LP | **yes** — `reason.ba1` |
| 3 | `>= recessive_max` and not P/LP | the 1e-2..0.05 band | **yes** — `reason.too_common` (but a zero-count reason is absent from the file, so "0 dropped" and "never emitted" look identical) |
| 3 | no impact / SpliceAI / CADD rung | the bulk drop | **yes** — `reason.not_functional` (but see the mis-filing findings) |
| 4 | `bcftools annotate -c INFO` re-attachment | nothing — but a transfer that lands partially or not at all is unchecked, and Step 5 then reports every allele as `rarity_basis=absent` | **NO** (critic INV-1) |
| 4 | `isec -c none -n=2` | trio records not allele-exactly in the plausible set | **NO** — only the post-isec count exists |
| 4 | `view --min-ac 1` | loci the child does not carry AND loci the child was not callable at — indistinguishable | **NO** |
| **5** | `require_pass and v.FILTER` | records non-PASS in this trio but PASS in another (verifiers: no loss vs policy — the record IS non-PASS here) | **NO** |
| **5** | sex unresolved | all non-PAR chrX/chrY for that trio | one stderr WARN per trio, no metric |
| **5** | chrY | every chrY call incl. father-to-son transmission | **NO** (documented) |
| **5** | `if gene:` | every gene-less inherited het — mode-asymmetric, since de novo / hom-rec / X-linked need no gene | **NO** |
| **5** | matched no mode | the entire no-row space (see the findings) | **NO — and no input counter exists** |
| **6** | `if not gene: continue` | calls with no gene attribution | **NO — and no input counter exists** |

**The instrumentation stops exactly where the sensitivity is decided.** Steps 1-4 are honestly
instrumented, Step 3 exemplarily so. Steps 5 and 6 record only what they emitted.

## Summary table

106 distinct findings: 92 from the first sweep and 14 from the critic's follow-up round (20 were
raised there; overlapping reports of one mechanism are merged, marked ×N). **Sev** is the verifiers'
median, with the finder's original in parentheses where it differs. **Ev** = reproduced by
execution / argued from the code. **R?** marks the 26 findings where one of the three
verifiers voted to refute — read those with the corrections in mind. FAM: FILTER = plausibility
filter, GENOTYPE = genotype-to-meaning, UPSTREAM = what feeds them.

| # | Sev | Class | Fam | Location | Ev | R? | Doc? | Finding |
|---|---|---|---|---|---|---|---|---|
| 1 | maj | silent loss | GENOTYPE | 05_inheritance_screen.py:305 | exec |  | **NEW** | A QC failure in the TRANSMITTING parent deletes the call silently; the identical failure in the NON-transmitting parent emits it w |
| 2 | maj | silent loss | GENOTYPE | 05_inheritance_screen.py:442 | exec |  | yes | Step 5 records only what it emitted: the entire no-row space is invisible, and comp-het row inflation makes the audit funnel read  |
| 3 | maj | wrong call | GENOTYPE | ped.py:46 | exec |  | **NEW** | A trios-file header spelling of `mother_id`/`father_id` silently swaps the parent roles for the whole trio, inverting every origin |
| 4 | maj | wrong call | UPSTREAM | ped.py:48 | exec |  | yes | A trios file whose parent columns are not exactly alias-named silently transposes mother and father, and the MIE gate is provably  |
| 5 | mod (maj) | silent loss | BOTH | 04_subset_and_annotate_trios.sh:264 | exec | R | **NEW** | Step 4's `bcftools annotate -c INFO` transfer — the one that populates the VCFs Step 5 reads — has no match count, no audit metric (×2) |
| 6 | mod (maj) | silent loss | BOTH | 05_inheritance_screen.py:469 | exec |  | **NEW** | Nothing verifies that Step 4's annotation re-attachment landed, and Step 5's join-coverage guard is structurally blind to it becau (×3) |
| 7 | mod (maj) | doc mismatch | GENOTYPE | 00_qc.py:44 | exec |  | yes | Step 0's Mendelian-error rate is mathematically symmetric under a father/mother swap and parental sex is never inferred, so the do |
| 8 | mod (maj) | silent loss | GENOTYPE | 05_inheritance_screen.py:165 | exec | R | yes | The ClinVar P/LP inert-band counter is scoped to recessive_max while the gate that actually kills an unpaired inherited het is dom |
| 9 | mod (maj) | silent loss | GENOTYPE | 05_inheritance_screen.py:210 | exec |  | yes | `has_hiconf` keys a per-record de novo requirement on a CALLSET-WIDE header line, so a declared-but-unpopulated (or wrong-child) h |
| 10 | mod (maj) | silent loss | GENOTYPE | 05_inheritance_screen.py:230 | exec |  | yes | A HOM_ALT child with fewer than two carrier parents produces no row under ANY mode at ANY frequency, with no fall-through and no c |
| 11 | mod (maj) | silent loss | GENOTYPE | 05_inheritance_screen.py:230 | exec |  | yes | A hom-alt child with one non-carrying or no-call parent produces zero rows and zero counters — the deletion-in-trans / UPD shape t |
| 12 | mod (maj) | wrong call | GENOTYPE | 05_inheritance_screen.py:274 | exec | R | yes | A spanning-deletion `*` in the child's genotype is invisible to Step 5: the co-located real allele is emitted as a plain inherited |
| 13 | mod (maj) | unverifiable | GENOTYPE | 05_inheritance_screen.py:274 | exec |  | yes | Compound-het trans evidence rests on the loose `hom_ref` band (AB <= 0.10), not the `clean_parent` test the de novo path uses (alt |
| 14 | mod (maj) | unverifiable | GENOTYPE | 05_inheritance_screen.py:288 | exec |  | **NEW** | The `origin="denovo"` classification rests on a `clean_parent` pass that can be vacuous, and it is the ONE such site in Step 5 tha |
| 15 | mod (maj) | silent loss | GENOTYPE | 05_inheritance_screen.py:296 | exec |  | yes | A deterministic parent-of-origin is discarded when the OTHER carrying parent's genotype QC fails, even though the determination do |
| 16 | mod (maj) | wrong call | GENOTYPE | 05_inheritance_screen.py:297 | exec |  | yes | The obligate-transmission arms (a 1/1 parent) require the OTHER, non-determining parent's genotype QC — so a genetically certain p |
| 17 | mod (maj) | silent loss | GENOTYPE | 05_inheritance_screen.py:305 | exec |  | yes | The transmitting parent's QC failure drops the call silently, while the non-transmitting parent's failure is fail-open and flagged |
| 18 | mod (maj) | silent loss | GENOTYPE | 05_inheritance_screen.py:313 | exec |  | yes | Never-drop-with-a-flag is applied when one parent is uninformative and the other CARRIES, but not when the other is a confident 0/ |
| 19 | mod (maj) | silent loss | GENOTYPE | 05_inheritance_screen.py:331 | exec |  | yes | A genuine biallelic hit loses BOTH legs when one leg fails any single collector gate — the pairing model has no partial credit and |
| 20 | mod (maj) | silent loss | GENOTYPE | 05_inheritance_screen.py:331 | exec | R | yes | A `both`-origin het (both parents 0/1) in the [dominant_max, recessive_max) band is emitted under no mode at all — the same unphas |
| 21 | mod (maj) | wrong call | GENOTYPE | 05_inheritance_screen.py:341 | exec |  | **NEW** | A compound-het pair whose trans evidence FAILED or was VACUOUS still consumes both legs, deleting the dominant call — violating th |
| 22 | mod (maj) | silent loss | GENOTYPE | genotype.py:40 | exec |  | **NEW** | The genotype-QC allele-balance knobs have no unit or range check, and no step guards a zero-call run: percent-scale values reduce  |
| 23 | mod (maj) | silent loss | GENOTYPE | genotype.py:129 | exec |  | yes | A trio VCF missing FORMAT/AD or FORMAT/GQ yields ZERO calls for that trio, silently; Step 4 warns only for AD and says calls 'may  |
| 24 | mod (maj) | fragility | GENOTYPE | genotype.py:129 | exec |  | yes | A callset with no FORMAT/AD (or no FORMAT/GQ) produces zero calls in every mode, exit 0, with a Step-4 audit identical to a health |
| 25 | mod (maj) | silent loss | GENOTYPE | genotype.py:149 | exec |  | yes | `clean_parent`'s alt-AD cap is an ABSOLUTE count, so its stringency inverts with depth — and a failure removes the child's variant |
| 26 | mod (maj) | unverifiable | GENOTYPE | igv.py:193 | exec |  | yes | Every Step-5 evidence flag except `origin=` is dropped from the reviewer-facing tables — an unphased compound-het reads as a confi (×2) |
| 27 | mod (maj) | silent loss | FILTER | 02_annotate_sites.sh:516 | exec |  | **NEW** | The SpliceAI keep-path reads only the selected CSQ block, so a scored splice variant on a non-picked overlapping gene is dropped a |
| 28 | mod (maj) | silent loss | FILTER | 03_select_plausible.py:46 | exec |  | yes | Step 3 hard-fails when the rarity field is undeclared but performs no equivalent check for the two functional-rung fields; on the  |
| 29 | mod (maj) | wrong call | FILTER | 05_inheritance_screen.py:427 | exec |  | **NEW** | Step 5 accepts a per-trio VCF carrying no oracle fields at all and reports the run as faf95; Step 3 refuses the identical input |
| 30 | mod (maj) | silent loss | FILTER | selection.py:39 | exec |  | **NEW** | A YAML scalar (or lower-case) `keep_impacts` silently deletes the entire impact rung — every stop_gained and missense is then drop |
| 31 | mod (maj) | silent loss | FILTER | selection.py:39 | exec |  | **NEW** | `filters.functional.keep_impacts` written as a YAML scalar silently disables the entire impact keep-path and files every HIGH-impa |
| 32 | mod (maj) | silent loss | FILTER | selection.py:39 | exec |  | **NEW** | `filters.functional.keep_impacts` is unvalidated: four plausible YAML shapes silently disable the entire impact rung, and the loss |
| 33 | mod (maj) | silent loss | FILTER | selection.py:42 | exec |  | yes | The impact rung reads only the MANE-Select transcript's IMPACT: a stop_gained on an alternative transcript of the SAME gene is fil |
| 34 | mod (maj) | silent loss | FILTER | selection.py:78 | exec |  | yes | A NOVEL non-coding indel (and any MNV) can reach no keep-path at all, and the drop is filed as `not_functional` — the same string  |
| 35 | mod (maj) | fragility | UPSTREAM | 04_subset_and_annotate_trios.sh:257 | exec |  | **NEW** | Step 4 rewrites the per-trio intermediates in the persistent tmpdir but index_vcf refuses to re-index them, so a stale index silen (×2) |
| 36 | mod (maj) | silent loss | UPSTREAM | resolve_trios.py:126 | exec |  | **NEW** | resolve_trios re-runs on every invocation and can re-point a trio onto the pipeline's OWN derived output, replacing the source cal |
| 37 | mod | silent loss | BOTH | 04_subset_and_annotate_trios.sh:207 | exec |  | yes | A trio VCF with no FORMAT/AD yields ZERO calls in every inheritance mode, and the only trace is a per-trio WARN that scopes the lo |
| 38 | mod | silent loss | GENOTYPE | 05_inheritance_screen.py:194 | exec |  | yes | The parental-mosaicism de novo the docs instruct you to detect is hard-gated away, and parental allele balance is never written to |
| 39 | mod | silent loss | GENOTYPE | 05_inheritance_screen.py:230 | exec |  | yes | hom_recessive truth table: a 1/1 child with one 0/0 parent (deletion-in-trans / UPD) or one no-call parent produces no row, no fla |
| 40 | mod | silent loss | GENOTYPE | 05_inheritance_screen.py:249 | exec |  | yes | A homozygous chrX call in a daughter is dropped whenever the father is not a pristine `1/1` — the exact discipline the affected-ma |
| 41 | mod | silent loss | GENOTYPE | 05_inheritance_screen.py:288 | exec |  | yes | `filters.denovo.parent_max_alt_ad` silently gates the RECESSIVE model: at the shipped default of 1 absolute alt read, a compound-h |
| 42 | mod | silent loss | GENOTYPE | 05_inheritance_screen.py:313 | exec |  | yes | A het child with one affirmatively clean parent and one no-call parent is emitted under no mode at all — 2 of the 64 truth-table c |
| 43 | mod | silent loss | GENOTYPE | genotype.py:140 | exec |  | **NEW** | The AB bands are pure fractions with no read-count floor and no depth-dependent widening, so at the configured `min_dp: 10` the he |
| 44 | mod | fragility | GENOTYPE | genotype.py:141 | exec |  | yes | A trio VCF with no FORMAT/AD yields ZERO inherited calls of any mode — every QC limb fails closed on allele balance, and the only  |
| 45 | mod | fragility | GENOTYPE | genotype.py:142 | exec |  | yes | `--keep-sum AD` fabricates reference reads under a HOM-ALT genotype, so a true homozygous hit at a multiallelic site fails the hom |
| 46 | mod | doc mismatch | FILTER | functional_annotation.md:139 | code |  | yes | docs claim CADD scores indels 'genome-wide'; it does not, and a private indel below MODERATE impact therefore has neither rescue r |
| 47 | mod | wrong call | UPSTREAM | 00_qc.py:107 | exec |  | yes | Step 0 opens the trio VCFs without `strict_gt`, so a half-called `1/.` counts as a chrX HET in sex inference — a male mis-inferred |
| 48 | mod | silent loss | UPSTREAM | 01_make_cohort_sites.sh:77 | exec |  | yes | Neither Step 1's nor Step 4's cache key includes the reference FASTA, and Step 4's omits HPRV_PLAUSIBLE_PAD — so a changed referen |
| 49 | mod | fragility | UPSTREAM | 02_annotate_sites.sh:684 | exec |  | yes | Step 2's rarity-oracle guards still police the retired oracle: the value-level guard is satisfiable by `vep_MAX_AF` alone, and the |
| 50 | mod | silent loss | UPSTREAM | ped.py:36 | exec |  | yes | A headerless trios file silently consumes its first trio as a header row |
| 51 | min (maj) | wrong call | GENOTYPE | 05_inheritance_screen.py:210 | exec | R | yes | The hiConfDeNovo gate re-imposes on the male-X hemizygous de novo path the exact paternal requirement that path exists to remove |
| 52 | min (maj) | unverifiable | FILTER | 02_annotate_sites.sh:643 | exec |  | yes | A PARTIALLY-covered gnomAD slim passes Step 2's 0-match die, Step 3's header guard and `prepare_resources verify`, and then report |
| 53 | min (maj) | silent loss | UPSTREAM | 04_subset_and_annotate_trios.sh:232 | exec | R | yes | Step 1 filters each trio to FILTER PASS/. but Step 4 re-derives genotypes with no FILTER at all, so a record filtered in trio A bu |
| 54 | min (mod) | unverifiable | BOTH | 05_inheritance_screen.py:31 | exec | R | **NEW** | `hprv_keep_reason` is carried into every per-trio candidate record and then read by nothing — the free per-record witness of the h |
| 55 | min (mod) | wrong call | BOTH | 05_inheritance_screen.py:101 | exec | R | **NEW** | `bcftools norm -m-` copies the site-level hiConfDeNovo tag to every split ALT, so the reported `hiConfDeNovo` column reads 1 on in |
| 56 | min (mod) | silent loss | BOTH | 05_inheritance_screen.py:123 | exec |  | **NEW** | Every Python boolean knob is `bool(value)`, so a quoted or `${ENV}`-templated `false` reads TRUE — the gates a user tried to relax |
| 57 | min (mod) | wrong call | GENOTYPE | 05_inheritance_screen.py:138 | exec |  | yes | `high_conf_rarity` never fires on the rarest class: a variant ABSENT from gnomAD is tiered below one gnomAD carries (×2) |
| 58 | min (mod) | wrong call | GENOTYPE | 05_inheritance_screen.py:169 | exec | R | yes | The male non-PAR hemizygous-het red flag is enforced on the proband only, never on the father — so a het-called father on chrX col |
| 59 | min (mod) | silent loss | GENOTYPE | 05_inheritance_screen.py:204 | exec |  | **NEW** | `denovo_min_dp` is applied unadjusted to the HAPLOID male chrX, so a hemizygous de novo at typical X depth is dropped with no fall |
| 60 | min (mod) | silent loss | GENOTYPE | 05_inheritance_screen.py:204 | exec | R | yes | The male-X hemizygous de novo carries a ploidy-blind DP >= 20 floor on a chromosome with half the autosomal depth, so an identical |
| 61 | min (mod) | unverifiable | GENOTYPE | 05_inheritance_screen.py:210 | exec | R | **NEW** | `use_hiconf_tag` gates the de novo MODE row but not the `origin="denovo"` label that admits a leg to trans pairing — the pipeline  |
| 62 | min (mod) | silent loss | GENOTYPE | 05_inheritance_screen.py:248 | exec |  | yes | A daughter's homozygous chrX call is deleted without trace when the father's chrX is a no-call or is rendered het — the opposite o |
| 63 | min (mod) | silent loss | GENOTYPE | 05_inheritance_screen.py:248 | exec |  | yes | chrX affected-female homozygote is lost unless the father is affirmatively HOM_ALT, and line 230 explicitly denies it the autosoma |
| 64 | min (mod) | silent loss | GENOTYPE | 05_inheritance_screen.py:249 | exec |  | yes | The female X-linked-recessive model hard-requires the father's chrX genotype, so the same record is emitted for a son and silently |
| 65 | min (mod) | unverifiable | GENOTYPE | 05_inheritance_screen.py:288 | exec |  | yes | The `parent_ad_unmeasured` witness is attached to the de novo ROW but never to the de novo ORIGIN that a compound-het pair rests o |
| 66 | min (mod) | fragility | GENOTYPE | 05_inheritance_screen.py:331 | exec | R | **NEW** | Compound-het pair enumeration is quadratic (k x m) with no cap or de-duplication, and is the mechanism that makes candidate_calls  |
| 67 | min (mod) | silent loss | GENOTYPE | 05_inheritance_screen.py:415 | exec | R | **NEW** | A PED sex code outside {0,'',None,1,2} silently discards Step 0's inference and blanks every non-PAR chrX/chrY call, with no warni |
| 68 | min (mod) | wrong call | GENOTYPE | 05_inheritance_screen.py:416 | exec | R | yes | A declared PED sex silently overrides a contradicting Step-0 chrX inference, and no output column records which sex the ploidy mod |
| 69 | min (mod) | fragility | GENOTYPE | run_pipeline.sh:212 | exec |  | yes | Step 0's cache key ignores VCF content, so a re-delivered trio VCF re-runs Steps 1 and 4 but serves a STALE inferred sex to Step 5 |
| 70 | min (mod) | fragility | GENOTYPE | genotype.py:32 | exec |  | **NEW** | `GtThresholds.from_config` validates nothing: an inverted het AB band silently deletes every dominant and compound-het call, and ` |
| 71 | min (mod) | doc mismatch | GENOTYPE | genotype.py:74 | exec |  | yes | `dp()` returns sum(AD) when AD is present and FORMAT/DP when it is absent — one threshold, two different depth quantities, and the |
| 72 | min (mod) | doc mismatch | GENOTYPE | genotype.py:132 | code | R | yes | There is no `denovo_min_gq`: the de novo arm is filtered at the same GQ >= 20 as inherited calls, on a posterior that CalculateGen |
| 73 | min (mod) | silent loss | GENOTYPE | genotype.py:146 | exec |  | yes | `parent_max_alt_ad = 1` is an absolute read count with no depth-relative alternative: two alt reads at DP 100 (2%) removes the var |
| 74 | min (mod) | silent loss | FILTER | prepare_resources.sh:217 | exec |  | yes | The faf95 oracle transfers frequencies from gnomAD records that FAILED gnomAD's own QC filters (AS_VQSR/AC0), and the gnomAD FILTE |
| 75 | min (mod) | fragility | FILTER | annotations.py:339 | exec |  | yes | Under `oracle: faf95` an absent faf95 resolves to 0.0 (rarest) regardless of how large the transferred `gnomad_AF_joint` is, so BA |
| 76 | min (mod) | silent loss | FILTER | annotations.py:434 | exec |  | **NEW** | `clnsig_is_plp` is vetoed by a co-asserted `likely_benign` but not by `benign`, so the ClinVar P/LP override silently fails on rea |
| 77 | min (mod) | silent loss | FILTER | genotype.py:140 | exec | R | yes | The fail-CLOSED half of the AD asymmetry deletes a carrier parent and its call with no witness; sample_qc_ad_measured only covers  |
| 78 | min (mod) | silent loss | FILTER | selection.py:70 | exec |  | **NEW** | No coherence check on the four `filters.rarity.*` cutoffs: `dominant_max > recessive_max` is silently inert, and `benign_ba1 < rec |
| 79 | min (mod) | unverifiable | FILTER | selection.py:78 | exec |  | **NEW** | `not_functional` is one bucket for 'the scores said no' and 'no score existed' — the distinction the rest of the pipeline treats a |
| 80 | min (mod) | fragility | UPSTREAM | 00_qc.py:90 | exec |  | yes | Step 0's chrX sex inference — the sole determinant of X ploidy in the shipped flow — ignores FILTER and is boundary-inclusive on t |
| 81 | min (mod) | fragility | UPSTREAM | 00_qc.py:107 | exec |  | yes | Step 0 computes the Mendelian-error rate WITHOUT `strict_gt`, so half-called parents Step 5 refuses to interpret are charged as Me |
| 82 | min (mod) | unverifiable | UPSTREAM | 00_qc.py:203 | exec | R | yes | Step 0 writes contam_flag=0 and contam_source=charr when CHARR could not be computed at all, and a real 0.15 contamination flag di |
| 83 | min (mod) | doc mismatch | UPSTREAM | 00_qc.py:221 | exec |  | yes | Step 0's sex-swap gate is unreachable in the shipped pipeline: `sex_match` is a constant 1, so `overall_pass` can never fail on se |
| 84 | min (mod) | unverifiable | UPSTREAM | 00_qc.py:241 | exec |  | yes | `overall_pass=1` is emitted for a trio on which NOTHING was measured — all three Step-0 gates fail open, and Step 0 writes no audi |
| 85 | min (mod) | doc mismatch | UPSTREAM | 01_make_cohort_sites.sh:30 | exec |  | **NEW** | `filters.genotype_qc.require_pass: false` cannot re-admit a non-PASS record — Step 1's `-f 'PASS,.'` is a hardcoded script default |
| 86 | min (mod) | silent loss | UPSTREAM | 01_make_cohort_sites.sh:148 | exec |  | yes | Step 1's only per-trio metric is named `input_sites` but counts the OUTPUT: the FILTER drop, the chrM exclusion and the multiallel |
| 87 | min (mod) | silent loss | UPSTREAM | 01_make_cohort_sites.sh:165 | exec |  | yes | Step 1 strips INFO before the union `norm`, so a symbolic <DEL> loses its INFO/END and bcftools left-aligns it to position 1 of th |
| 88 | min (mod) | silent loss | UPSTREAM | 01_make_cohort_sites.sh:169 | exec |  | **NEW** | `bcftools view --min-ac 1` trusts a stale INFO/AC instead of recomputing from genotypes, and it is applied on branches that never  |
| 89 | min (mod) | unverifiable | UPSTREAM | 02_annotate_sites.sh:642 | exec | R | yes | Step 2's gnomAD and ClinVar transfer coverage is a stderr log line only, never an audit metric, and the faf95 count has no guard a |
| 90 | min (mod) | silent loss | UPSTREAM | 04_subset_and_annotate_trios.sh:107 | exec |  | yes | Step 4's region-BED restriction is not the pure optimization its comment claims: an indel whose raw representation sits more than  |
| 91 | min (mod) | fragility | UPSTREAM | prepare_resources.sh:182 | exec | R | **NEW** | `prepare_resources.sh` keeps a gnomAD slim that failed its own integrity check and adopts it as `cached` on the next run; `verify` |
| 92 | min | doc mismatch | BOTH | 05_inheritance_screen.py:138 | exec |  | yes | `tag_strict` inverts on the rarest class: a variant absent from gnomAD never earns `high_conf_rarity`, while a `zero_ci` variant a |
| 93 | min | unverifiable | BOTH | 05_inheritance_screen.py:152 | code | R | yes | FORMAT/FT is never read anywhere in the pipeline: a genotype the caller itself marked as failing is treated as a confident call, i |
| 94 | min | fragility | GENOTYPE | 05_inheritance_screen.py:162 | exec |  | yes | Step 5 asserts no biallelic precondition; an un-normalised multiallelic record silently yields zero calls because ref_AD is 0 and  |
| 95 | min | doc mismatch | GENOTYPE | 05_inheritance_screen.py:184 | exec | R | **NEW** | `inheritance.emit_denovo: false` does not omit de novo rows: the de novo variant is still emitted, as a compound_het leg |
| 96 | min | unverifiable | GENOTYPE | 05_inheritance_screen.py:341 | exec |  | yes | A phase-confirmed pair consumes both legs, so an ultra-rare dominant-grade variant is relabelled compound_het by a partner up to 1 |
| 97 | min | doc mismatch | GENOTYPE | 05_inheritance_screen.py:353 | exec |  | yes | The methods reference promises `origin_unverified` for an unqualified hom-ref parent, but the code emits two flags that appear in  |
| 98 | min | doc mismatch | GENOTYPE | annotations.py:159 | code |  | yes | `loConfDeNovo` is transferred and declared but read by nothing, while two docs mark the lower-sensitivity tier IMPLEMENTED |
| 99 | min | fragility | GENOTYPE | genotype.py:113 | exec |  | **NEW** | chrX/chrY ALT and unlocalized contigs are not recognised as sex chromosomes, so a male's het calls there are emitted as autosomal  |
| 100 | min | fragility | FILTER | annotations.py:295 | exec | R | **NEW** | `rarity_oracle()` silently coerces any unrecognised `resources.gnomad.oracle` value to `faf95`, so a config typo selects the oppos |
| 101 | min | fragility | FILTER | annotations.py:295 | exec | R | **NEW** | `resources.gnomad.oracle` is normalised two incompatible ways: Python maps every unrecognised spelling to `faf95`, the shell compa |
| 102 | min | doc mismatch | FILTER | selection.py:36 | code |  | **NEW** | Step 3's classifier documents a permissive union over 'the looser of the dominant/recessive cutoffs' but hardcodes `recessive_max` |
| 103 | min | unverifiable | FILTER | selection.py:41 | exec | R | **NEW** | Rung order makes `reason.spliceai` / `reason.cadd` count only the rescues that had no impact rung, so the audit under-reports the  |
| 104 | min | doc mismatch | UPSTREAM | CLAUDE.md:312 | exec | R | **NEW** | CLAUDE.md:312-314 and docs/pipeline_design.md:228 promise every `bcftools annotate` transfer dies on a 0-match join; the third tra |
| 105 | min | fragility | UPSTREAM | 01_make_cohort_sites.sh:165 | exec | R | yes | Step 1 strips INFO/END, so a symbolic ALT record enters the annotated union with no span |
| 106 | min | silent loss | UPSTREAM | 04_subset_and_annotate_trios.sh:232 | exec |  | **NEW** | Step 4's -R fast path and its whole-genome fallback are not equivalent: an indel whose left-alignment shift exceeds HPRV_PLAUSIBLE |

## Findings — major (verifiers' median)

The four findings three verifiers could neither overturn nor talk down. Two corrupt every call in a trio rather than one call.

### M1. A QC failure in the TRANSMITTING parent deletes the call silently; the identical failure in the NON-transmitting parent emits it with a flag

**`pipeline/05_inheritance_screen.py:305`** · silent loss · GENOTYPE · verifiers: confirmed/partly/confirmed · reproduced by execution · **not documented anywhere**

`mom_ok`/`dad_ok` (:262-265) gate `origin` at :297/:299/:303/:305/:309 — a transmitting parent that fails its own zygosity band yields `origin = None`, the variant never enters `hets` (:314), and no counter records it. The mirror-image evidence — the NON-transmitting parent failing the same gate — sets `unverified`/`vacuous` (:306-311) and the call is emitted with `origin_unverified` / `parent_ad_unmeasured`. The never-drop machinery therefore exists and is applied to the weaker claim (non-transmission) but not to the stronger one (transmission), while the child's own genotype is perfect in both.

**Failure scenario.** One VCF, one gate (GQ), two role assignments. chr4:1000 ASYM_NONTX — child 0/1 AB 0.50 GQ 99, mother 0/1 GQ 99 (transmitting), father 0/0 GQ 12 (non-transmitting): EMITTED as `dominant flags=origin=mat;origin_unverified`. chr4:1010 ASYM_TX — child 0/1 AB 0.50 GQ 99, mother 0/1 GQ 12 (transmitting), father 0/0 GQ 99: DELETED, no row, no counter. Same with FORMAT/AD absent instead of GQ: AD_NONTX (father 0/0, AD absent) emits `dominant flags=origin=mat;parent_ad_unmeasured`; AD_TX (mother 0/1, AD absent) is deleted. Correct behaviour: emit with a `transmitting_parent_qc_failed` flag, as the code already does for the other parent, or at minimum audit the count.

**Who it hits.** Every inherited het (dominant + both comp-het legs) whose transmitting parent has allelic imbalance, low GQ/DP, or a missing AD field. Parental AB drift outside 0.25-0.75 is routine in repeat-adjacent and GC-extreme regions; the child's call is unaffected. docs/inheritance_and_genotype_qc.md:170 states the transmitting-parent QC requirement but nothing states that failing it produces silence rather than a flag.

**What the verifiers corrected.**

- Two corrections, both making the finding stronger, plus one documentation nuance.

1) The asymmetry is not limited to GQ and missing AD — it holds for EVERY sample_qc gate. Executed and confirmed for four: GQ (12 vs 99), DP (8 vs 40, below min_dp 10), het AB band (0.80 vs 0.50), and FORMAT/AD absence. In all four, the…
- The executed facts hold: a transmitting parent failing het/hom-alt QC (GQ, DP, AB, or absent AD) yields origin=None at 05_inheritance_screen.py:297/299/303/305/309, the variant never enters `hets` (:314-316), and no audit metric records it — so it is invisible downstream (Step 8 reads candidates.calls.tsv, 08_igv_expor…
- Two refinements, neither weakening the verdict.

(a) The finding's "CLAIMED DOCUMENTATION STATUS: NO" is right for the GQ/DP/AB gates but slightly overstated for the AD gate. The AD-absent drop IS documented at CLAUDE.md:487-488 ("`het`/`hom_alt`/`denovo_child` require `ab is not None` and so DROP a carrier when FORMAT…

**Fix.** When `mom_ok`/`dad_ok` fails but the parent's GT carries the alt, still set `origin` and emit with a `transmitting_parent_qc_failed` flag; otherwise audit a per-trio `origin_unresolved` counter.

### M2. Step 5 records only what it emitted: the entire no-row space is invisible, and comp-het row inflation makes the audit funnel read as a GAIN where 13% of one trio's candidates were lost

**`pipeline/05_inheritance_screen.py:442`** · silent loss · GENOTYPE · verifiers: confirmed/partly/confirmed · reproduced by execution · documented at PARTIAL — docs/limitations.md:303-304 asks for 'at minimum a dropped-count audit metric' for ONE class (the male chrX maternal-no-call residual). The general absence of any Step-5 input counter, and the fact that row inflation makes the funnel read as a gain, are documented nowhere.

All ten `audit.record` calls in Step 5 (:442-486) record OUTPUTS — `candidate_calls`, `mode.*`, `trios_screened`, `rarity_basis.*`. There is no `variants_examined` metric and none of the three `continue`s (:152, :156, :159) nor any of the mode-condition failures is counted. Because a comp-het leg emits one row per pair and one variant can appear under up to three modes, `candidate_calls` routinely EXCEEDS Step 4's `candidate_genotypes`, so no arithmetic a reviewer can perform on audit/counts.tsv reveals a loss.

**Failure scenario.** On the shipped integration run: `04_subset CH_A candidate_genotypes 31` and `05_inheritance CH_A candidate_calls 35`. 35 > 31 makes loss look arithmetically impossible, yet only 26 distinct input records produced a row — 4 of 31 vanished (child GQ 12; half-called father; a ClinVar Pathogenic at faf95 1.6e-4; a missense at faf95 1.6e-3). Exhaustively: over the full clean-QC grid of 5 sequence contexts x 3 rarity bands x 64 (child,father,mother) genotype triples, restricting to configurations in which the child CARRIES the alt at the rarest band, 102 of 160 produce no row (autosome 16/32, chrX-PAR male 16/32, chrX female 18/32, chrX male 20/32, chrY male 32/32). Exactly ONE of those 102 is a true Mendelian error (het child of two `1/1` parents). ZERO of the 102 is counted, flagged, or given a reason string. A reviewer reading a negative result cannot distinguish any of them from 'no such variant existed'.

**Who it hits.** Every no-row class in the taxonomy: rarity band-gap, parental uninformativeness, hom-alt single-carrier, sex-chromosome routing, and every genotype-QC failure. It is the mechanism by which all the other findings here are invisible rather than merely conservative.

**What the verifiers corrected.**

- Three corrections to supporting detail; none touches the mechanism, the severity, or the fix.

(a) "All ten `audit.record` calls in Step 5 record OUTPUTS" is inaccurate. One IS a drop counter: `clinvar_plp_dropped_ge_recessive_max` (:444 per-trio, :484 global). But its existence SHARPENS the finding rather than weakeni…
- Four corrections to the finding as written:

(1) "no arithmetic a reviewer can perform on audit/counts.tsv reveals a loss" is FALSE. Σ(`04_subset candidate_genotypes`) − `09_prioritize observations_distinct` = 37 − 33 = 4, exactly the no-row count. Step 4's counter is `bcftools index -n` (pipeline/lib/common.sh:186-192…
- Four corrections, none of which touch the core claim:

1. "All ten `audit.record` calls record OUTPUTS" is inaccurate. `clinvar_plp_dropped_ge_recessive_max` (05_inheritance_screen.py:444, :484) IS a drop counter — the P4 fix from the 2026-07 review. It is bounded to `fr >= rec_max` (:170-172), and that bound provably…

**Fix.** Add per-trio `variants_examined`, `candidates_no_mode`, and per-`continue` counters (`skipped_filter`, `skipped_sex_unresolved`, `skipped_y_nonpar`) at pipeline/05_inheritance_screen.py:442, and report `distinct_observations` beside `candidate_calls` so the funnel is monotone.

### M3. A trios-file header spelling of `mother_id`/`father_id` silently swaps the parent roles for the whole trio, inverting every origin and fabricating false X-linked calls

**`src/hprv/ped.py:46`** · wrong call · GENOTYPE · verifiers: confirmed/confirmed/confirmed · reproduced by execution · **not documented anywhere**

`read_trios_file` matches header names by EXACT membership in `_DAD_ALIASES`/`_MOM_ALIASES` (ped.py:22-24, loop at :38-45), then falls back to a positional index PER COLUMN and INDEPENDENTLY (:46-51). Any suffixed spelling (`father_id`, `mother_id`, `paternal_id`, `dad_id`) matches nothing, so that column silently reverts to a fixed position while a recognised `sample_id`/`kid` column is still resolved by name — mixing name resolution with a positional assumption, with no warning and no record of which columns were used.

**Failure scenario.** A collaborator hands over `sample_id<TAB>mother_id<TAB>father_id` / `CH<TAB>MO<TAB>FA` (mother listed first, the alphabetical order of a spreadsheet export). kid is matched by name to col 0; dad falls back to col 1 (= the MOTHER) and mom to col 2 (= the FATHER). resolve_trios reports "1/1 trios resolved", writes a PED with MO in the father column, and every Step-5 call for that trio is mis-assigned: both autosomal `dominant` rows have `flags=origin=` inverted, the genuine `x_linked_recessive` call (mother the carrier) DISAPPEARS with no counter, and a FABRICATED `x_linked_recessive` is emitted on the father's hemizygous X allele — an allele a son cannot inherit — carrying `flags=high_conf_rarity`. Row COUNT is unchanged (3 vs 3), so no count-based check sees it. It should either resolve all three columns positionally or all three by name, and hard-stop (or at minimum WARN naming the resolved column indices) when only some header names are recognised.

**Who it hits.** EVERY call in the affected trio. Autosomal dominant/comp-het rows keep their trans logic (mat x pat pairing is symmetric) but carry an inverted `origin=`; all hemizygous X models, which key on the MOTHER's genotype (05_inheritance_screen.py:180 `male_x_chrx`), are inverted — true maternal-carrier X-linked recessive calls are lost and paternal X alleles are reported as maternally transmitted. Affects 100% of trios whose row was transcribed under a non-alias header with mother-before-father ordering.

**What the verifiers corrected.**

- Two refinements, neither weakening the finding.

1. The X-linked outcome is a MISLABEL first and a loss second. Under the swapped PED the genuine maternally-inherited X-linked recessive at ordinary depth is emitted as `denovo_x_hemi`, not dropped — the true father (0/0) sits in the mother slot and trips the male-X de n…
- Three refinements, none of which rescue the code:
1. "Row COUNT is unchanged (3 vs 3), so no count-based check sees it" is an artifact of the finder's engineered mock. On the real integration trio the total changed (6 -> 5) and the per-mode audit changed (x_linked_recessive 2 -> 1). The loss is still silent in the oper…
- The finding is correct in mechanism and location (src/hprv/ped.py:22-24, 38-51), is genuinely undocumented, and is untouched by the ten fix commits on this branch. Four corrections/refinements to its write-up:

1. WRONG RATIONALE in the failure scenario: "mother listed first, the alphabetical order of a spreadsheet exp…

**Fix.** Make header resolution all-or-nothing: if ANY of the three columns is matched by name, require all three by name (else `return`/raise), and log the resolved column indices; accept `*_id` suffixes in the alias sets.

### M4. A trios file whose parent columns are not exactly alias-named silently transposes mother and father, and the MIE gate is provably blind to it

**`src/hprv/ped.py:48`** · wrong call · UPSTREAM · verifiers: confirmed/confirmed/confirmed · reproduced by execution · documented at NO (adjacent only: docs/ROADMAP.md:43 requests somalier for CROSS-COHORT duplicate/swap QC, a different failure, and it is an open item)

read_trios_file matches parent columns against an EXACT alias list (_DAD_ALIASES=('dad','father','paternal'), _MOM_ALIASES=('mom','mother','maternal')) and, when neither matches, falls back to POSITIONAL indices dad=1, mom=2 (ped.py:46-51) with no warning. Nothing downstream can detect the transposition: mendelian_violation (00_qc.py:44-52) is exactly symmetric in its dad/mom arguments, and scan_sex (00_qc.py:116) only ever inspects the CHILD, so parent sex is asserted and never verified.

**Failure scenario.** A collaborator supplies `proband_id / maternal_id / paternal_id` (or `sample_id / mother_id / father_id`, or `kid / parent_1 / parent_2`) — none of those parent names is in the alias list, and the maternal column comes first. read_trios_file returns ('CH1','MO1','FA1') as (kid,dad,mom); resolve_trios writes a PED with the mother in the father slot. Every dominant call then reports the wrong parent-of-origin, and a genuine maternally-inherited hemizygous X-linked recessive (son 1/1, mother 0/1, father 0/0) is re-read as a paternal-side event: with the shipped `use_hiconf_tag: true` it is DELETED (0 rows, no counter), and with the tag off it is emitted as `denovo_x_hemi`. It should have been `x_linked_recessive` with origin from the mother.

**What the verifiers corrected.**

- Two refinements, neither of which weakens the finding.

(1) The mechanism is slightly worse than stated. The fallbacks at ped.py:46-51 fire INDEPENDENTLY per role, not as an all-or-nothing "neither parent matched -> positional". So a header that matches the KID alias but not the parents assigns dad to column 1 — which…
- Two scope corrections; neither touches the core mechanism.

1. THE TITLE OVERSTATES THE TRIGGER. The fallback is POSITIONAL (kid=0, dad=1, mom=2), not "transposing". A non-alias header transposes only when the parent columns are in mother-then-father order. `proband_id / paternal_id / maternal_id` — the conventional PE…
- Two refinements, neither weakening the finding. (1) The finding says this is undocumented; it is stronger than that — the docs affirmatively claim the opposite. CLAUDE.md:129 ("header names kid/dad/mom in any order") and docs/pipeline_design.md:196-197 ("any column order") assert unconditional order-independence, which…

**Fix.** Accept `*_id`-suffixed aliases, and make the positional fallback a hard stop (or a WARN plus an audit metric) instead of assuming column order; verify parent sex genotypically before trusting the roles.

## Findings — moderate

Confirmed mechanisms whose verifiers judged the triggering configuration real but uncommon in GMKF-shaped trio data, or the loss bounded. Full write-ups for the ones the finders had called major; one-line entries for the rest.

### D1. Step 4's `bcftools annotate -c INFO` transfer — the one that populates the VCFs Step 5 reads — has no match count, no audit metric and no 0-match die, and an under-match deletes every dominant and compound_het call while denovo/hom_recessive survive

**`pipeline/04_subset_and_annotate_trios.sh:264`** · silent loss · BOTH · verifiers: partly/refuted/confirmed · reproduced by execution · **not documented anywhere**

04:264 `bcftools annotate -a "$PLAUSIBLE_TX" -c INFO -o "$out" "$cand"` is the third and most consequential annotation transfer in the pipeline, and unlike Step 2's two external transfers (which each count matched records and `die` at 0 — 02:580-588 ClinVar, 02:638-651 gnomAD) it is followed only by `require_intact_bgzip` (04:268) and a record COUNT (`audit 04_subset candidate_genotypes`, 04:273). A partial or total under-match leaves every record present with a fully-populated INFO header and empty values, so the audit number is unchanged and `.done` is stamped. Step 5 then reads `gene = A._str(v,"gene") or A.symbol(v)` at 05:258 and the `if gene:` at 05:259 gates the ENTIRE `hets` collector, which is the sole source of both `dominant` and `compound_het` rows; `denovo`, `hom_recessive` and `x_linked_recessive` need no gene and still fire, so the loss is mode-selective rather than total.

**Failure scenario.** EXECUTED on the real integration run's CH_A candidate VCF (31 records, healthy: 35 rows = 2 denovo + 2 hom_recessive + 14 compound_het + 17 dominant). (a) ZERO-match transfer (INFO values stripped, header lines preserved exactly as a 0-match `annotate -c INFO` leaves them): 35 rows -> 4 rows; all 17 dominant and all 14 compound_het calls vanish, only denovo 2 + hom_recessive 2 remain. (b) PARTIAL transfer (annotations survive on chr1, lost from chr2 on — the 'high-coordinate' shape 04's own comment at :104-109 says happened in the field): 35 rows -> 8 rows, and all four modes are still present so the output shape looks like a small healthy run. In both cases `bcftools view -H` still shows 31 records, so `audit 04_subset CH_A candidate_genotypes` records 31 exactly as in the healthy run, `require_intact_bgzip` passes, `.done` is stamped, and the run exits 0. Should instead: verify per-record transfer coverage and die. `INFO/hprv_keep_reason` is the exact witness — Step 3 stamps it on 100% of plausible sites (03_select_plausible.py:71) and `isec -c none` guarantees every $cand record has an allele-exact plausible match, so coverage must be exactly 100% with no 'rare cohort' caveat (verified 31/31 healthy, 10/31 partial, 0/31 zero-match). Nothing in the repo reads `hprv_keep_reason` today (grep over pipeline/ src/ tests/ returns nothing), so it is a free, exact detector.

**Who it hits.** Every inherited het — the pipeline's headline model ('the dominant model — recurrent inherited rare functional hets — is the new emphasis', CLAUDE.md) — plus every compound-het leg, for any trio whose transfer under-matched. On the integration fixture that is 31 of 35 rows (89%) for a total loss and 27 of 35 (77%) for a partial one.

**What the verifiers corrected.**

- The mechanism, the consequence and every number are correct and were reproduced with a real `bcftools annotate -c INFO` under-match (35 rows -> 4 on a 0-match, losing all 17 dominant and all 14 compound_het; 35 -> 8 on a partial), with exit 0, `bgzip -t` passing, `bcftools index -n` = 31 and `audit 04_subset candidate_…
- 04:264's transfer cannot under-match: `bcftools isec -c none -n=2 -w1` at 04:259 makes every `$cand` record allele-exactly present in `$PLAUSIBLE`, and `$PLAUSIBLE_TX` (04:89-91) is that same record set minus two INFO tags — so 100% transfer coverage is structural, not incidental. That is why Step 2's ClinVar/gnomAD tr…
- Three corrections, none of which refute the finding.

1. SEVERITY: major -> moderate, on trigger frequency. The mechanism and blast radius are exactly as claimed and I reproduced both, but the triggering configuration cannot arise from GMKF trio data shape. `bcftools isec -c none -n=2 -w1` at 04:259 guarantees every `$…

**Fix.** After 04:264, count `bcftools query -i 'INFO/hprv_keep_reason!="."'` in $out, audit it as `candidate_annotations_transferred`, and die unless it equals `count_variants "$out"`.

### D2. Nothing verifies that Step 4's annotation re-attachment landed, and Step 5's join-coverage guard is structurally blind to it because its witness rides in the same transfer

**`pipeline/05_inheritance_screen.py:469`** · silent loss · BOTH · verifiers: confirmed/partly/confirmed · reproduced by execution · **not documented anywhere**

Step 4 re-attaches every annotation with one blanket `bcftools annotate -a "$PLAUSIBLE_TX" -c INFO` (04_subset_and_annotate_trios.sh:264) that has no match count and no 0-match die, unlike Step 2's two transfers; Step 5 then opens that copy (05:427) with no header check and no per-record check before stamping rarity_af/rarity_oracle/rarity_basis for the whole downstream. The one guard that looks like it covers this, `rarity_faf95_absent_but_cache_has_af` (05:469), tests `rarity_basis == "absent" and grpmax_af` -- but grpmax_af is a `vep_*` field carried by the SAME `-c INFO` transfer, so when the transfer fails both sides vanish together and the counter reads 0. Step 3 has the equivalent header guard (03_select_plausible.py:46-56) and run_pipeline.sh's oracle preflight (:104,:113) is gated on `run_step 2`, so a `--from 4`/`--from 5` invocation (documented in CLAUDE.md) reaches Step 5 with none of the three checks having run.

**Failure scenario.** A per-trio candidate VCF whose `-c INFO` transfer did not land -- a hand-built or resumed Step-5 manifest, a Step-4 regression, or the stale-index class of bug the code already fixed once for PLAUSIBLE_TX (04:92-99). Executed on the real integration run's CH_A candidate VCF: 31 candidate genotypes that legitimately produce 35 calls produce 4. All 17 `dominant` and all 14 `compound_het` calls disappear, because `05:259 if gene:` never collects a het whose vep_Gene/vep_SYMBOL is blank; only the gene-less-tolerant `denovo` and `hom_recessive` modes survive, so the loss is silently biased against precisely the inherited focus of the pipeline. The 4 survivors are written with rarity_oracle=faf95, rarity_basis=absent and blank rarity_af -- a positive claim that gnomAD holds no record for the allele, i.e. the RAREST verdict -- and the same blank becomes igv/variants.tsv's headline `frequency` column. The guard records 0, no WARN is printed, exit 0. The inverse control proves the blind spot: strip ONLY the gnomad_* fields (a Step-2 partial-join failure, which is what the guard was written for) and it correctly fires with 29 orphans. It should die, exactly as Step 2 dies on a 0-match transfer and Step 3 dies on a missing witness header.

**Who it hits.** Every candidate in an affected trio. The gene-keyed inherited modes (dominant, compound_het) lose 100 percent -- 31 of 35 calls in the demonstration -- while denovo and hom_recessive survive un-annotated, so a reviewer sees a short but non-empty call set rather than an obviously empty one.

**What the verifiers corrected.**

- Three refinements. (1) The finding understates the blindness: a header guard of the Step-3 kind would NOT catch the realistic in-run failure, because a 0-match `bcftools annotate` still writes every source INFO header line — verified, `zeromatch_A.vcf.gz` declares `ID=gnomad_AF_joint` while zero records carry a value.…
- The mechanism is accurately described and I reproduced it, but the triggering state is not producible by the pipeline. `bcftools isec -c none -n=2 -w1 "$norm" "$PLAUSIBLE"` (04:257) guarantees every record handed to the `-c INFO` transfer has an exact allele match in PLAUSIBLE_TX, which is PLAUSIBLE with two never-pres…
- Two corrections, both to severity framing rather than to the mechanism.

SEVERITY: moderate, not major. The blind spot and the total absence of verification are real, demonstrated and unfixed, but the finding's own triggering configuration does not arise on its own in an ordinary GMKF run. `$cand` is produced by `bcfto…

**Fix.** Assert INFO/hprv_keep_reason is present on every record at 05_inheritance_screen.py:427 (its absence is unambiguous proof of transfer failure, because every candidate is by construction a plausible-set member, whereas an absent gnomad_faf95 is a legitimate biological state), and add a transferred-record count + die to 04_subset_and_annotate_trios.sh:264.

### D3. Step 0's Mendelian-error rate is mathematically symmetric under a father/mother swap and parental sex is never inferred, so the documented swap backstop has ZERO power against a parent-label swap

**`pipeline/00_qc.py:44`** · doc mismatch · GENOTYPE · verifiers: confirmed/partly/confirmed · reproduced by execution · documented at docs/ROADMAP.md:43 (somalier, listed as an open roadmap item, calls the chrX-only sex check "fragile"); docs/cram_access_phase.md:28 scopes the undetectable case to CROSS-COHORT duplicate/swap. Both UNDERSTATE: neither says the within-trio parent swap is undetectable, and the docstring at 00_qc.py:10 plus docs/inheritance_and_genotype_qc.md:128 actively claim MIE covers "mislabeled parents".

`mendelian_violation(gc, gd, gm)` tests only whether the parental allele POOL can produce the child's genotype (:44-52); every branch is symmetric in `gd`/`gm`. `scan_sex` is called with `ped["child"]` only (:116), so the parents' chrX heterozygosity is never computed and a male sample sitting in the mother slot is invisible. The module docstring (:10) calls the MIE rate "a sensitive proxy for sample swaps / mislabeled parents" and docs/inheritance_and_genotype_qc.md:128 says a high rate "signals a bad trio/swap"; neither is true for the parent-label swap, the one swap this pipeline's own hand-transcribed trios file can introduce (see the ped.py finding).

**Failure scenario.** A trio whose dad/mom labels are exchanged (bad transcription, or the ped.py header fallback above) is scanned by Step 0 and receives a byte-identical qc_report row to the correctly-labelled trio: same `mie_rate`, same `mie_flag`, `sex_match=1`, `overall_pass=1`, and "0 trio(s) flagged" on stderr. The reviewer reads the QC sheet and the IGV `sample_qc.tsv` as confirmation that the roles are right. Step 5 then emits inverted origins and false/missing X-linked calls for every variant in that trio, and Step 6 aggregates them as ordinary carriers. Step 0 should infer sex for all three members (a male in the mother slot is a one-line check) and/or state plainly that MIE cannot detect a parent swap.

**Who it hits.** Every call in a parent-swapped trio (see the ped.py finding for the concrete downstream damage). Frequency depends on transcription hygiene, but the pipeline offers no detection at any rate.

**What the verifiers corrected.**

- Two refinements, neither reducing the finding. (a) "IDENTICAL QC REPORT" is exact for all QC verdict columns, but byte-identity holds only when the two parents have equal contamination — on the real integration trios the dad_contam/mom_contam VALUES swap slots (0.15 moves from dad_contam to mom_contam) while mie_rate,…
- Three corrections. (1) The QC row is not byte-identical: dad_contam/mom_contam VALUES exchange columns under the swap (0.15 moved dad->mom in trio CH_B). Every QC DECISION column — mie_rate, mie_flag, inferred_sex, sex_match, contam_flag, overall_pass — is identical, which is the load-bearing claim and is exact. (2) Th…
- The finding stands on substance — not fixed, not documented, and every executed claim reproduces — but three framing points need correcting.

1. NOT A BLANKET OVER-CLAIM. `mendelian_violation` is blind ONLY to the father<->mother exchange, where both members remain true parents. It retains full power against every othe…

**Fix.** Run `scan_sex` for the father and mother as well and flag any member whose inferred sex contradicts its PED role; document that MIE is blind to a dad/mom exchange.

### D4. The ClinVar P/LP inert-band counter is scoped to recessive_max while the gate that actually kills an unpaired inherited het is dominant_max — so the docs' "now audited rather than silent" claim is false for the band that fires

**`pipeline/05_inheritance_screen.py:165`** · silent loss · GENOTYPE · verifiers: partly/refuted/partly · reproduced by execution · documented at docs/limitations.md:292-294 (PARTIAL — documents the [1e-4,1e-2) no-mode drop but scopes it to the comp-het gene-key mismatch; docs/clinical_classification.md:87 asserts the opposite of this finding, that the P/LP inert-band drop "is now audited ... rather than silent", while scoping the audit to [recessive_max, benign_ba1) only)

Hets are pooled into `hets` at `recessive_max` (1e-2, line 257) but the dominant emission gates at `dominant_max` (1e-4, line 366), so an inherited het in [1e-4, 1e-2) with no trans partner is emitted under no mode. `n_plp_inert` — the one counter written to make that drop auditable (its comment at :146-148 even says "all gate at rec_max/dom_max") — tests only `_fr >= rec_max`, a band the variant never reaches. Its code comment also says "grpmax AF" while the code reads `A.frequency()` (faf95 under the default oracle).

**Failure scenario.** Child het for a ClinVar Pathogenic stop_gained, gnomAD faf95 = 2e-4, mother 0/1 AB 0.50 GQ 99 (transmitting), father affirmative 0/0 AD 40,0 GQ 99, no second hit in the gene. Step 3 rescues it (`hprv_keep_reason=clinvar_plp`), Step 4 genotypes it, Step 5 collects it as origin=mat and then discards it at the 1e-4 dominant gate. Output: zero rows, and `audit/counts.tsv` records `clinvar_plp_dropped_ge_recessive_max 0`. It should either emit the row or count it; a reviewer reads the 0 as "no ClinVar P/LP variant was lost".

**Who it hits.** Every functional inherited het with rarity_af in [1e-4, 1e-2) and no trans partner in its gene — on real WGS the default state of most surviving hets — and specifically every ClinVar P/LP het in that band, which is the class the counter was written for.

**What the verifiers corrected.**

- The MECHANISM reproduces exactly and is real: pipeline/05_inheritance_screen.py:165 counts only `_fr >= rec_max` (1e-2) while the dominant emission at :366 gates at `dom_max` (1e-4), so a ClinVar P/LP inherited het with faf95 in [1e-4, 1e-2) and no trans partner is emitted under no mode AND not counted — executed: GENE…
- The counter is correctly scoped and the documentation is accurate. `n_plp_inert` covers `[recessive_max, benign_ba1)` — a reachable band I demonstrated fires (n_plp_inert=1 at faf95=2e-2, both het and hom-alt children) — and docs/clinical_classification.md:87 explicitly scopes its "now audited" claim to that same band.…
- The core mechanism is real and I reproduced it twice, including in the repo's own shipped integration run (chr2:8000 GENE7, P/LP, faf95=1.6e-4, clean paternal transmission -> zero rows in candidates.calls.tsv while audit/counts.tsv reports `clinvar_plp_dropped_ge_recessive_max 0` alongside `reason.clinvar_plp 2`). It i…

**Fix.** Test `_fr >= dom_max` (or emit a second metric for the [dom_max, rec_max) band) and add a per-trio `candidates_no_mode` counter.

### D5. `has_hiconf` keys a per-record de novo requirement on a CALLSET-WIDE header line, so a declared-but-unpopulated (or wrong-child) hiConfDeNovo silently deletes every de novo call in the run

**`pipeline/05_inheritance_screen.py:210`** · silent loss · GENOTYPE · verifiers: confirmed/partly/partly · reproduced by execution · documented at NO — and docs/gene_burden.md:137 states the opposite reassurance ("IMPLEMENTED (only required when the tag is present in the callset header)"), presenting the header check as the safety valve when the header is precisely the wrong thing to key on

Line 74 sets `has_hiconf = "ID=hiConfDeNovo" in vcf.raw_header` — a HEADER fact — and line 210 then requires the per-RECORD tag to name this child (`require_hiconf and gt.has_hiconf and not A.is_hiconf_denovo_for(...)` -> `ok=False`). The header survives every downstream transform (Step 4 keeps it explicitly at 04:214 `-x '^INFO/hiConfDeNovo,INFO/loConfDeNovo'`) whether or not a single record carries the tag, so "GATK declared the annotation" is silently taken to mean "GATK evaluated this variant and said no".

**Failure scenario.** EXECUTED. Two VCFs identical in every record; one declares `##INFO=<ID=hiConfDeNovo,...>` in the header, the other does not. A textbook clean autosomal de novo (GDN1: child 0/1 AD=20,20 DP=40 GQ=99; both parents 0/0 AD=40,0 DP=40 GQ=99; stop_gained, HIGH, no gnomAD record) is emitted as `mode=denovo` from the header-less file and produces NO ROW AND NO COUNTER from the header-bearing one. Run totals: 9 rows -> 5 rows. Three concrete real-data triggers, none of which errors: (a) the trio VCF was re-headered after GATK (`bcftools reheader -s`, the routine BS_*->subject-ID remap on Kids First data) — reheader renames the sample columns but never rewrites INFO strings, so hiConfDeNovo still lists the OLD IDs while the PED/VCF carry the new ones, and `is_hiconf_denovo_for` is False for every record in the callset; (b) a multi-family VCF where GATK's PED covered a different family; (c) per this repo's OWN docs (docs/inheritance_and_genotype_qc.md:64, "Exact criteria (biallelic sites, ...)") PossibleDeNovo evaluates only biallelic sites, while Steps 1 and 4 `bcftools norm -m-` split multiallelics into perfectly ordinary biallelic de novo records that were never eligible for a tag. Undetectable: Step 5 audits only emitted modes (`tmodes` is keyed on observed events), so `mode.denovo` is simply ABSENT from audit/counts.tsv at zero — verified on the shipped run, where `mode.x_linked_recessive` is likewise absent for CH_A and `mode.denovo_x_hemi` never appears at all. Also inconsistent within one function: the same variant that the gate refuses to call `denovo` still enters `hets` as `origin="denovo"` at :288 and is emitted as a `compound_het` leg flagged `unphased_denovo_partner` — verified in the header-bearing run — so the output asserts an unphased de novo partner for a variant the same run declined to call de novo.

**Who it hits.** Every de novo in the run when the callset was re-headered / GATK's PED disagreed; otherwise every de novo at a site that was multiallelic in the raw trio VCF (~2-5% of WGS sites, enriched in repeats where the DNM rate is highest).

**What the verifiers corrected.**

- The mechanism, the silence, the audit blindness, the compound_het inconsistency and the absent regression coverage all reproduce exactly as described. Three corrections/qualifications:

1. THE FINDING OMITS THAT THIS IS DELIBERATE AND UNIT-TESTED. tests/test_pure.py:2902-2903 asserts precisely this behaviour as CORRECT…
- The core mechanism is real and I reproduced it independently: `has_hiconf` (05:74) is a header fact, :210 turns it into a per-record requirement, the drop is uncounted, and Step 4's exact `bcftools annotate -x '^INFO/hiConfDeNovo,INFO/loConfDeNovo'` provably preserves the header line with zero tagged records — so the d…
- Not fixed and not documented — both gate answers are NO, so the finding stands. Corrections: (a) severity is moderate, not major — the gate touches only `denovo`/`denovo_x_hemi`, a mode this repo explicitly designates secondary and out of scope (CLAUDE.md scope paragraph; docs/gene_burden.md:93), and executed checks co…

**Fix.** Require the tag only when this callset actually populates it (probe for >=1 tagged record per trio, or count tagged records) and, when the gate fires, emit the call with `flags=hiconf_tag_absent` plus a per-trio `denovo_gated_by_hiconf` audit metric instead of dropping.

### D6. A HOM_ALT child with fewer than two carrier parents produces no row under ANY mode at ANY frequency, with no fall-through and no counter — including for a ClinVar Pathogenic allele

**`pipeline/05_inheritance_screen.py:230`** · silent loss · GENOTYPE · verifiers: confirmed/partly/partly · reproduced by execution · documented at PARTIAL, and the documentation is actively misleading — docs/inheritance_and_genotype_qc.md:159 tells the reader 'A 1/1 child with a 0/0 parent is a Mendelian error -> suspect a hemizygous false hom (deletion on the other allele) or UPD (see §4)', and §4 (docs/inheritance_and_genotype_qc.md:198-200) calls UPD a flag rather than a primary call. No such configuration is ever emitted or counted, so the reader is instructed to suspect something the output can never show them. docs/limitations.md's 'Inheritance-model residuals' list (:282-317) omits it entirely.

`hom_recessive` (:230-231) requires BOTH parents in (HET, HOM_ALT). A hom-alt child is also excluded from the het collector (:257 requires `gc == G.HET`), so there is no fall-through to `dominant`. Every configuration with one carrier parent plus a `0/0` or no-call parent therefore produces nothing. The one counter that exists for lost clinically-significant alleles, `n_plp_inert` (:163-166), only fires when `frequency >= rec_max`, so it stays 0 for a rare pathogenic allele lost this way.

**Failure scenario.** Child `1/1` AD 0,40 DP 40 GQ 99, `vep_CLIN_SIG=pathogenic`, gnomad_faf95 1e-6; mother `0/1` carrier; father `0/0` clean (or `./.`). Step 5 emits no row and `n_plp_inert = 0`. This is the canonical signature of (a) a point mutation in trans with a CNV deletion — CLAUDE.md names CNV as the pipeline's biggest coverage gap, so the SNV half of a CNV+SNV biallelic hit is also invisible; (b) uniparental disomy; (c) parental allele dropout. Making both parents carriers emits `hom_recessive` — the only difference is the second parent's call.

**Who it hits.** Every homozygous-alt candidate in the child where one parent is hom-ref, no-call, or QC-failing: the hemizygous 'false hom' (deletion in trans), UPD, and parental allele dropout. On real WGS these are exactly the biallelic candidates a rare-disease screen most wants to see.

**What the verifiers corrected.**

- Three corrections, none of which changes the verdict.

1. SCOPE IS BROADER THAN CLAIMED (finding understates itself). I probed a case the finding does not make: with BOTH parents genuine `0/1` carriers, a single parent failing genotype QC produces the identical zero-row, zero-counter outcome. Executed (qcfail.py), chil…
- The code mechanism is exactly as claimed and I reproduced it: a HOM_ALT child with fewer than two carrier parents yields no row under any mode at any frequency, and n_plp_inert stays 0. Two corrections to the finding's framing.

FIRST, the documentation status. The finder calls the docs "actively misleading"; they are…
- Three corrections. (1) The title's "under ANY mode" is wrong: on chrX non-PAR with a male child, a hom-alt (hemizygous) child with only a carrier MOTHER does emit x_linked_recessive (executed: dad 0/0 or ./. with mom 0/1 both emit). The loss is confined to autosomes and to female chrX (which requires father 1/1). (2) T…

**Fix.** Emit the call with `mode=hom_alt_single_carrier_parent` + `flags=mendelian_error_suspect_deletion_or_upd` (never-drop), or at minimum audit `hom_alt_child_uncarried_parent` per trio.

### D7. A hom-alt child with one non-carrying or no-call parent produces zero rows and zero counters — the deletion-in-trans / UPD shape the docs tell the reader to "suspect"

**`pipeline/05_inheritance_screen.py:230`** · silent loss · GENOTYPE · verifiers: confirmed/partly/partly · reproduced by execution · documented at docs/inheritance_and_genotype_qc.md:159 and :256 (PARTIAL — they describe the shape and instruct the reader to "suspect a hemizygous false hom (deletion on the other allele) or UPD", which is unachievable because nothing is emitted; docs/limitations.md:283-316 "Inheritance-model residuals" omits this case). The documentation UNDERSTATES the impact by framing it as a CNV-calling blind spot rather than a Step-5 silent drop.

`hom_recessive` requires BOTH parents in (HET, HOM_ALT) at :230-231, and the het collector at :257 requires `gc == G.HET`, so a HOM_ALT child has no fall-through path. There is no counter for the rejection and no flag anywhere in the output.

**Failure scenario.** Child 1/1 for a frameshift at faf95 1e-5, AD 0,40 GQ 99; mother 0/1 carrier; father an affirmative 0/0 (AD 40,0 GQ 99) — the classic point mutation in trans with an undetected exon deletion, or UPD. Zero rows, nothing in audit/counts.tsv. Same result with the father a `./.` no-call. The pipeline should emit it with a `mendelian_error_hemizygous_candidate` flag (never-drop) or at minimum count it, so the reviewer the docs address has something to inspect.

**Who it hits.** Every recessive candidate where the second allele is a CNV, a UPD segment, or a parental allele dropout — the documented 10-15% CNV/SV diagnostic fraction reaches Step 5 exactly in this genotype shape.

**What the verifiers corrected.**

- Three corrections to the finding as written, none of which change the verdict.

1. OVERSTATED IMPACT SCOPE. "the documented 10-15% CNV/SV diagnostic fraction reaches Step 5 exactly in this genotype shape" is not true. Most CNV/SV diagnoses are CNV-only (whole-gene or multi-exon deletions/duplications with no SNV partne…
- The code behaviour is exactly as described and I reproduced it, but three of the finding's claims are wrong or overstated.

1. SCOPE. The loss is autosomal (and female chrX) only. A male non-PAR chrX hemizygote with a 0/0 or `./.` father IS emitted as `x_linked_recessive` (05_inheritance_screen.py:236-242 deliberately…
- The behaviour is real, current, and reproducible by execution, and it is NOT fixed on this branch. But it is already documented — most precisely at docs/ROADMAP.md:48 (open item #6, High impact / Low effort): "The recessive logic currently deletes the UPD case (1/1 child + 0/0 parent -> 'Mendelian error')". The reporte…

**Fix.** Emit with `flags=mendelian_error_possible_deletion_or_upd` under never-drop, or add a per-trio counter for hom-alt children with a non-carrying parent.

### D8. A spanning-deletion `*` in the child's genotype is invisible to Step 5: the co-located real allele is emitted as a plain inherited het with a `child_gt` naming a reference allele the child does not carry and a falsely-clean transmitting parent

**`pipeline/05_inheritance_screen.py:274`** · wrong call · GENOTYPE · verifiers: confirmed/refuted/partly · reproduced by execution · documented at NO for the plausibility-filter / genotype-meaning families. src/hprv/igv.py:112 is the only place in the repo that acknowledges `*` at all, and it is Step 8b's NHF join. docs/inheritance_and_genotype_qc.md and docs/limitations.md's "Inheritance-model residuals" do not mention spanning deletions.

`bcftools norm -m-` splits `A  T,*` into an `A>T` record (child rendered `1/0`) and an `A>*` record; Step 5 reads only the record in front of it. `gt_types` for the `A>T` record reports HET, and `dad_clear`/`mom_clear` (05:274-275) test `sample_qc(..., "hom_ref")` on the parent who transmitted the deletion — who after the split reads `0/0` with AD `[38,0]`, i.e. an AFFIRMATIVE, QC-passing, allele-depth-MEASURED hom-ref, so neither `origin_unverified` nor `parent_ad_unmeasured` fires. Nothing in selection.py, genotype.py or 05_inheritance_screen.py ever inspects the child's own `*` allele, and `base_row` (05:96,100) writes `gt_bases` verbatim.

**Failure scenario.** KID `1/2` at chr1:202 (maternal T in trans with a paternal 4 bp frameshift deletion at chr1:200), DAD `0/2`, MOM `0/1`. TWO wrong outcomes, one per branch. (i) The deletion record does not itself survive Step 3 (a novel intronic/UTR deletion — see the indel finding): Step 5 emits exactly one row, `mode=dominant flags=origin=mat child_gt=T/A father_gt=A/A`. The child has no A at 202 and the father transmitted a deletion covering it; a true biallelic hit is reported as a single dominant het with phase affirmed. (ii) The `*` record IS annotated and kept: Step 5 emits a SECOND compound_het pair `T1:CH2` pairing the same maternal SNV with the `*` SHADOW of the deletion already paired in `T1:CH1` — one biological compound het becomes two pair_ids and four rows, and chr1:202 A>* counts as a distinct (trio,chrom,pos,ref,alt) observation for the gene. A related third effect: `--keep-sum AD` (04:206) folds the deletion-supporting reads into AD[0], so the child's allele balance for the real allele is alt/(alt + deletion reads) and is judged against the diploid het band 0.25-0.75 (genotype.py:141) even though the child is functionally hemizygous there — with AD `0,15,45` the AB is exactly 0.250, one read from silent rejection. Should instead: detect `*` in the child's genotype, flag the co-located allele `child_hemizygous_spanning_deletion`, never pair the `*` record itself, and never write `father_gt=A/A` for a parent whose transmitted haplotype is deleted.

**Who it hits.** Every site where a deletion called in one trio member overlaps a variant site in another — GATK emits `*` there by construction, so this is routine in per-trio jointly-genotyped GATK output rather than an edge case. Branch (i) hits SNV-in-trans-with-a-deletion, the classic recessive second hit; branch (ii) hits any gene carrying both.

**What the verifiers corrected.**

- Three refinements, none of which weaken the finding materially:

(a) In branch (ii) the PHASE the pipeline infers is actually CORRECT — the maternal T genuinely is in trans with the paternal deletion. The harm there is duplication (two pair_ids, four rows, a spurious `A>*` observation), not a wrong trans call. The find…
- What is actually true: (1) `bcftools norm -m-` does render the co-located `A>T` leg with the deletion-carrying parent as an affirmative, AD-measured `0/0` — I confirmed this — but that is correct for the T allele and is the very evidence the trans test needs; with the deletion partner present (the finder's own frameshi…
- Three corrections. (1) `src/hprv/igv.py:112` is not the only acknowledgement of `*` in the repo — the Step-2b SpliceAI-backfill `_is_indel` helper also excludes it, guarded by `tests/test_pure.py:476`. Both are outside the plausibility-filter and genotype-meaning families, so "undocumented for these two families" still…

**Fix.** Treat `*` as a first-class allele in Step 5: skip `*` records for pairing/mode assignment, and when the child's raw genotype at a locus contains `*`, emit the co-located call with a `spanning_deletion_in_trans` flag instead of an unqualified dominant/hom-ref assertion.

### D9. Compound-het trans evidence rests on the loose `hom_ref` band (AB <= 0.10), not the `clean_parent` test the de novo path uses (alt AD <= 1) — a parent showing 4 alt reads yields an unflagged "confirmed trans" pair

**`pipeline/05_inheritance_screen.py:274`** · unverifiable · GENOTYPE · verifiers: confirmed/partly/partly · reproduced by execution · documented at NO. docs/inheritance_and_genotype_qc.md:166 and :277 say the non-transmitting parent must be an "affirmative QC-passing 0/0" and that "an unqualified 0/0 (allele dropout) ... is not evidence of anything", without saying which of the two available predicates is used; the code uses the looser one.

`mom_clear`/`dad_clear` (:274-275) call `sample_qc(..., "hom_ref")`, which at genotype.py:145 tests only `ab <= homref_ab_max` (0.10). The de novo path calls `sample_qc(..., "clean_parent")` (genotype.py:146-149), which additionally requires `alt_ad <= parent_max_alt_ad` (1). So the ONLY evidence a comp-het pair is in trans tolerates up to 10% alt reads in the "non-transmitting" parent — the signature of an under-called parental het, which is the single most likely way a real CIS pair is reported as TRANS.

**Failure scenario.** Non-transmitting parent called 0/0 with AD 40,4 (AB 0.091), DP 44, GQ 99 — precisely what the pipeline's own documented failure mode produces (inheritance_and_genotype_qc.md:131: CalculateGenotypePosteriors' gnomAD prior pushes a genuine rare het toward hom-ref). The pair is emitted as `compound_het` with NO flags: not `origin_unverified` (that fires only when the test FAILS) and not `trans_evidence_unmeasured` (that fires only when AD is ABSENT). Parent AD/AB/GQ are not among `COLS` (:31-54), so the trans claim cannot be audited from candidates.calls.tsv at all.

**Who it hits.** Every compound_het pair in a callset with any parental allelic imbalance or refinement-driven het->hom-ref suppression; the false-trans direction manufactures biallelic diagnoses.

**What the verifiers corrected.**

- Two immaterial precision fixes. (1) The dirty pair's flags are literally empty only above the 1e-3 high-confidence rarity tier; at the report's own faf95=2e-3 I reproduced flags=[] exactly, while at faf95=2e-5 the column reads "high_conf_rarity" — a rarity tag, not a phase tag — so the claim (no phase flag, output indi…
- The comp-het trans clearance does use the looser `hom_ref` band (AB ≤ 0.10, no alt-count cap) rather than `clean_parent` (alt AD ≤ 1) — executed and confirmed at pipeline/05_inheritance_screen.py:274-275 vs genotype.py:145-149, with an unflagged `compound_het` pair emitted for a non-transmitting parent at AD 40,4. But:…
- The core mechanism is real and I reproduced it independently: 05_inheritance_screen.py:274-275 gates comp-het trans evidence on `hom_ref` (AB <= 0.10, genotype.py:145), not on `clean_parent` (alt AD <= 1, genotype.py:146-149) which the de novo branch uses at :288-289. A non-transmitting parent at 0/0 AD=40,4 (AB 0.0909…

**Fix.** Use `clean_parent` for trans evidence, or emit the non-transmitting parent's alt AD in the call row and flag `trans_parent_alt_reads=N`.

### D10. The `origin="denovo"` classification rests on a `clean_parent` pass that can be vacuous, and it is the ONE such site in Step 5 that consults no witness and sets no flag

**`pipeline/05_inheritance_screen.py:288`** · unverifiable · GENOTYPE · verifiers: confirmed/partly/confirmed · reproduced by execution · **not documented anywhere**

Line 288-289 calls `G.sample_qc(v,m,thr,'clean_parent')` and the same for the father to decide that an unexplained het is de novo — the label that makes the leg eligible for trans pairing. Unlike the de novo MODE branch (:199,:202) and the one-carrier branch (:307,:311), it never calls `G.sample_qc_ad_measured`, and the both-hom-ref branch leaves `vacuous` at its `False` initialisation (:283). `clean_parent` returns True when AD is absent (genotype.py:148-149: `a is None or a <= parent_max_alt_ad`), so two GATK ref-block `0/0` parents affirm the de novo on zero allele evidence, and the pair emitted at :353-358 can never carry `trans_evidence_unmeasured`.

**Failure scenario.** Trio, gene GVAC, two functional hets at faf95 5e-4 (inside the recessive band). Leg A chr2:3000 child 0/1 AD=20,20 DP=40 GQ=99, mother 0/1, father 0/0 AD=40,0 -> origin=mat, measured. Leg B chr2:3500 child 0/1 AD=20,20 DP=40 GQ=99, BOTH parents `0/0:.:40:99` (ref-block shape, AD absent) -> `clean_parent` passes vacuously for both -> origin=denovo. The two legs are emitted as `compound_het` pair T1:CH1 with flags `high_conf_rarity;unphased_denovo_partner`. What SHOULD happen: the pair carries `trans_evidence_unmeasured` (or leg B carries `parent_ad_unmeasured`) exactly as the de novo row and the dominant row do for the identical parental shape. The code comment at :285-287 says this branch exists to stop 'a dropped-out parental het masquerading as de novo and getting paired cis' — and the predicate it uses passes vacuously in precisely the case where dropout is invisible. `base_row` (:80-104) emits no parental DP/AD/AB columns, so a reviewer cannot recover the vacuity from `candidates.calls.tsv` either.

**Who it hits.** Any compound-het whose second leg is (or looks) de novo in a callset where `0/0` genotypes come from GVCF ref blocks and therefore carry no AD — the dominant `0/0` shape in GATK per-trio VCFs, so potentially most comp-het pairs with a de novo leg. The pipeline's own docs (docs/inheritance_and_genotype_qc.md:166) and CLAUDE.md:486 build the whole trans-evidence story on this witness; this branch is the hole in it.

**What the verifiers corrected.**

- The finding is accurate as written, with three refinements. (1) The witness is not merely unconsulted — `mom_meas`/`dad_meas` are computed at :280-281 and then discarded unread inside the both-hom-ref branch, so the fix costs nothing. (2) Documentation status is stronger than "NO": the behaviour is actively counter-doc…
- The code behaviour is exactly as described and I reproduced it, but this is an observability/consistency gap, not a sensitivity or correctness defect. Nothing is dropped: the vacuous `clean_parent` pass RETAINS a leg that would otherwise vanish with no row (`origin=None`), and unphased pairs do not consume their legs (…
- Three corrections/refinements, none of which overturn the finding.

1. SCOPE. The claim "the pair emitted at :353-358 can never carry trans_evidence_unmeasured" is exactly right, but total invisibility of the vacuity in the OUTPUT FILE is narrower than the finding implies. When leg B's frequency is below dominant_max (…

**Fix.** In the `gmm == HOM_REF and gd == HOM_REF` branch, set `vacuous = not (G.sample_qc_ad_measured(v,m,'clean_parent') and G.sample_qc_ad_measured(v,d,'clean_parent'))` so the existing `wa or wb` -> `trans_evidence_unmeasured` path at :353 fires.

### D11. A deterministic parent-of-origin is discarded when the OTHER carrying parent's genotype QC fails, even though the determination does not use that parent — worst on chrX, the case the doc wrote the branch for

**`pipeline/05_inheritance_screen.py:296`** · silent loss · GENOTYPE · verifiers: partly/partly/confirmed · reproduced by execution · documented at docs/inheritance_and_genotype_qc.md:150 (PARTIAL and self-contradictory — :145 justifies the branch as "deterministic, because a 1/1 parent transmits the alt obligately", which needs only the 1/1 parent's call, while :150 requires "each carrying parent confident on its own zygosity band"; the resulting silent drop is documented nowhere)

At :296-299 the deterministic branch sets `origin = "mat"`/`"pat"` only `if (mom_ok and dad_ok)`. When one parent is HOM_ALT the origin is obligate from that parent alone (the child is het, so the other parent transmitted ref regardless of its zygosity), yet a QC failure on the non-transmitting carrier collapses origin to None and the variant leaves `hets` entirely — no dominant row, no pairing, no counter.

**Failure scenario.** chrX non-PAR, female proband 0/1 AB 0.50 GQ 99, father 1/1 (a hemizygous carrier rendered 1/1 by a diploid caller — the configuration docs/inheritance_and_genotype_qc.md:167 says "matters most on chrX"), mother 0/1 with AB 0.182 (36,8 at DP 44, just under the 0.25 het floor). A father transmits his single X to every daughter, so the alt is obligately paternal whatever the mother's call is; the code emits nothing. The autosomal twin (mother 1/1, father 0/1 at AB 0.182) fails the same way.

**Who it hits.** Any inherited het where one parent is hom-alt and the other carrying parent has an allele-imbalanced or low-GQ call — systematically the chrX carrier-father case, since a diploid caller renders every hemizygous father 1/1 and the mother's het is then the QC bottleneck.

**What the verifiers corrected.**

- Two corrections. (1) PREVALENCE: the deterministic branch at pipeline/05_inheritance_screen.py:290-299 is guarded by "elif mom_carries and dad_carries", so it fires only when BOTH parents carry the allele. The plain chrX case the finding calls systematic — hemizygous carrier father rendered 1/1, mother hom-ref or no-ca…
- The defect is real and reproduces, but its scope is narrower than stated. The deterministic branch at 05_inheritance_screen.py:290-299 is entered ONLY when both parents carry the allele (`mom_carries and dad_carries`). The ordinary chrX hemizygous-carrier-father configuration the finding headlines - father 1/1, mother…
- Three corrections. (1) Scope: the deterministic branch fires only when BOTH parents carry, so the chrX emphasis is overstated — executed control, chrX dad 1/1 + mom 0/0 with the identical AB-0.18 call, emits `dominant origin=pat;origin_unverified`. A lone hemizygous carrier father is NOT affected; the mother must also…

**Fix.** In the HOM_ALT x HET branch require only the hom-alt parent's QC, and record the other parent's QC failure as a flag rather than a drop.

### D12. The obligate-transmission arms (a 1/1 parent) require the OTHER, non-determining parent's genotype QC — so a genetically certain parent-of-origin is discarded on the other parent's allele balance

**`pipeline/05_inheritance_screen.py:297`** · wrong call · GENOTYPE · verifiers: confirmed/partly/partly · reproduced by execution · documented at NO — and docs/inheritance_and_genotype_qc.md:145 and :167 assert the opposite ("a 1/1 parent transmits the alt obligately … deterministic … so it DOES pair in trans"). docs/inheritance_and_genotype_qc.md:150 states the conjunction as a QC rule ("each carrying parent confident on its own zygosity band") without noting that it defeats the determinism claimed two lines earlier, and docs/limitations.md:282-316 does not list it.

Lines :297 and :299 read `origin = "mat" if (mom_ok and dad_ok) else None` and `origin = "pat" if (mom_ok and dad_ok) else None`. When one parent is HOM_ALT its transmission is obligate and the child's single alt is that parent's, full stop — the other parent's GQ/DP/AB has no bearing on the inference. The `and` conjunction nonetheless makes the non-determining parent's QC a veto; `origin` becomes None, `if origin:` at :314 never appends to `hets`, and the variant is emitted under no mode and counted by nothing.

**Failure scenario.** Autosome, chr1:1100, child 0/1 AD 20,20 DP 40 GQ 99 (pristine); mother 1/1 AD 0,40 DP 40 GQ 99 (obligate transmitter, passes hom_alt band); father 0/1 AD 32,8 (AB 0.20, fails the het band). Correct call: dominant, origin=mat, certain. Actual: nothing emitted, no flag, no counter. The identical record with the father at AD 20,20 emits dominant[origin=mat]. Same result when the father instead has DP 8 (P1c) or GQ 12 (P1d), and symmetrically when the father is 1/1 and the MOTHER's het fails (P1e). The chrX case is worse because it also destroys the phase-CONFIRMED compound het the block was written to rescue: gene XCH on chrX non-PAR, daughter het at two loci, father 1/1 hemizygous carrier at leg 1, mother a clean het — emits compound_het pair T1:CH1. Gene XCH2, identical except the mother's het at leg 1 has AB 0.20 — leg 1 emits NOTHING, the pair is gone, and leg 2 is downgraded to a lone dominant[origin=mat].

**Who it hits.** Every inherited het where one parent is HOM_ALT and the other is HET. Autosomally this is founder alleles, consanguinity, and affected parents (the exact cases the code comment at :291-295 lists). On chrX for a DAUGHTER it is the general case, because a diploid caller renders every hemizygous carrier father as 1/1 — the code's own comment says the arm exists for that. Loss rate is the other parent's het-QC failure rate: ~11% at DP 10, 3.5% at DP 15, 1.2% at DP 20 (binomial floor; real AB is heavier-tailed), plus every parent below DP 10 or GQ 20 outright.

**What the verifiers corrected.**

- The finding is accurate as written; two refinements. (1) The claim that the docs "assert the opposite" applies to inheritance_and_genotype_qc.md:145 and :167 (the determinism claims) but not to :150 — the QC rule there ("each carrying parent confident on its own zygosity band") literally describes what the code does, s…
- The core defect is confirmed: pipeline/05_inheritance_screen.py:297,299 gate a genetically certain parent-of-origin on the NON-determining parent's genotype QC, and when that QC fails the variant is emitted under no mode with no flag and no audit counter (Step 5 records only emitted rows, :442-446/:484). Reachable — St…
- Three corrections. (1) chrX: "for a DAUGHTER it is the general case" is FALSE — the arm at :290 fires only when BOTH parents carry at the same site. Executed control: dad 1/1 hemizygous + mom 0/0 -> dominant{origin=pat}, and a failing mother yields only origin_unverified. The X trigger is "carrier mother AND hemizygous…

**Fix.** On the HOM_ALT arms gate only on the determining parent (`origin = "mat" if mom_ok else None`), and record the other parent's QC failure as a flag (e.g. `other_parent_qc_fail`) rather than a veto.

### D13. The transmitting parent's QC failure drops the call silently, while the non-transmitting parent's failure is fail-open and flagged — a parent with NO data preserves the call, a parent with slightly imperfect data destroys it

**`pipeline/05_inheritance_screen.py:305`** · silent loss · GENOTYPE · verifiers: confirmed/partly/confirmed · reproduced by execution · documented at docs/inheritance_and_genotype_qc.md:150 and :170 state the QC requirement, but nothing documents the consequence (a silent, uncounted drop of a pristine child call) or the asymmetry against the `origin_unverified` treatment two lines away in the code. docs/limitations.md:282-316 omits it. The docs UNDERSTATE: a reader of :150 would expect the call to be flagged, as the adjacent :151 rule flags the non-transmitting parent.

Lines :305 / :309 set `origin = "mat" if mom_ok else None` / `"pat" if dad_ok else None`; `if origin:` at :314 then silently omits the variant from `hets`, so it gets neither a dominant row (:366) nor a compound-het leg (:331), and none of Step 5's ten audit metrics (:442-486) counts an examined-but-unemitted variant. Two lines further on, :306 / :310 handle the NON-transmitting parent's identical QC failure the opposite way — `unverified = not dad_clear` — which emits the call and marks it `origin_unverified` (never-drop). The same class of missing parental evidence therefore fails closed on one parent and fails open on the other.

**Failure scenario.** chr1:2200, a stop_gained at faf95 1e-6: child 0/1 AD 20,20 DP 40 GQ 99; father 0/1 AD 5,4 DP 9 (one read below min_dp); mother a clean 0/0 at DP 40. Nothing is emitted and nothing is counted. Now chr1:2400, same variant: father 0/1 AD 20,20 DP 40 (clean transmitter); mother 0/0 AD 2,0 DP 2 GQ 20 — a parent with essentially no data at all. That one IS emitted, as dominant[origin=pat;origin_unverified]. The configuration with MORE parental information (a father called het at DP 9) is discarded; the configuration with LESS (a mother with 2 reads) survives with a flag. In both the child's own call is identical and pristine.

**Who it hits.** Every inherited het with exactly one carrying parent — i.e. the pipeline's headline dominant mode and every compound-het leg. Under a binomial(DP,0.5) read model a TRUE het parent fails the 0.25-0.75 AB band with probability 0.109 at DP 10, 0.035 at DP 15, 0.012 at DP 20, 0.005 at DP 30 (executed with scipy). Real het AB is heavier-tailed than binomial (reference bias at indels, mapping bias), and on top of that every parent below DP 10 or GQ 20 fails outright. So a low-single-digit percentage of all dominant candidates genome-wide, concentrated in exactly the low-coverage regions where recall already suffers.

**What the verifiers corrected.**

- The finding is correct as written. Two consequences it UNDERSTATES, both reproduced by execution:

(a) IT IS NOT ONLY A LOSS -- IT SILENTLY CONVERTS A RECESSIVE DIAGNOSIS INTO A DOMINANT ONE. Same gene CHGENE, two rare stop_gained hets in the child (both AD 20,20 DP 40 GQ 99), leg 1 maternal (mom 0/1 AD 20,20; dad clea…
- The mechanism is real and I reproduced it verbatim: a transmitting parent that fails its own zygosity-band QC yields `origin = None` (05_inheritance_screen.py:305/:309), the variant never enters `hets` (:314), no row is emitted under any mode, and none of Step 5's ten audit metrics counts it. `igv/variants.tsv` is buil…

**Fix.** Apply never-drop symmetrically: emit the call with `flags=transmitting_parent_qc_fail` (and the failing limb) instead of returning None, or at minimum add a per-trio `origin_none_transmitting_parent_qc` audit metric.

### D14. Never-drop-with-a-flag is applied when one parent is uninformative and the other CARRIES, but not when the other is a confident 0/0 — the second case emits nothing and counts nothing

**`pipeline/05_inheritance_screen.py:313`** · silent loss · GENOTYPE · verifiers: confirmed/partly/partly · reproduced by execution · documented at PARTIAL — tests/integration/make_mock_data.py:126-128 states the drop is intended for the half-call case; docs/limitations.md:297-304 documents only the chrX-male analogue and asks there for 'at minimum a dropped-count audit metric'. The autosomal case, the QC-failing-0/0 variants of it, and the inconsistency with the adjacent origin_unverified branch are documented NOWHERE. The documentation understates the impact by scoping it to chrX.

At :304-311 (`elif mom_carries` / `elif dad_carries`) a het child with ONE carrying parent and one uninformative parent is emitted with `origin_unverified`. At :312-313 the neighbouring configuration — one parent uninformative (no-call, half-call, or a called 0/0 that fails GQ/DP QC) and the other an affirmative clean 0/0 — falls to `origin = None`, so the variant never enters `hets` at :314-316 and no mode can fire; the de novo branch at :184 also fails because it requires both parents `HOM_REF` and an UNKNOWN parent is not HOM_REF. No counter, no flag, no reason string.

**Failure scenario.** Child `0/1` at a `stop_gained`, AD 20,20 / DP 40 / GQ 99, faf95 absent (rarest); father `0/.` (a GATK half-call, DP 18 GQ 45); mother `0/0` AD 40,0 DP 40 GQ 99. Step 5 emits NOTHING. Change the father to `0/1` and the identical variant is emitted as `dominant flags=origin=pat;origin_unverified`. Same for father `0/0` at DP 8 (below `parent_min_dp` 10) or GQ 15 (below `min_gq` 20) — both configurations are the ordinary low-coverage parent, not a pathology. It should be emitted with a `parent_uninformative` flag exactly as the adjacent branch does, or at minimum counted.

**Who it hits.** Every rare functional het in the child where exactly one parent's genotype is uninformative and the other is a confident hom-ref. Parental no-calls / low-depth calls are ubiquitous on real WGS, so this is a systematic recall loss on the pipeline's HEADLINE mode (recurrent inherited dominant hets). One of the 4 uncounted losses in the repo's own integration fixture is this exact shape.

**What the verifiers corrected.**

- Two refinements, neither changing the verdict. (1) Line attribution: the two "low-coverage 0/0 parent" examples in the report (father 0/0 at DP 8, and at GQ 15, with a clean 0/0 mother) do NOT reach :312-313 — both parents read HOM_REF, so they fall through the de-novo-origin branch at :288-292 where `origin = "denovo"…
- The drop is the specified, tested contract, not a silent loss of a callable variant. A het child with one parent uninformative and the other an affirmative clean 0/0 has NO observed transmitting parent, so it satisfies neither the dominant definition (docs/inheritance_and_genotype_qc.md:135, "transmitted from >= 1 pare…
- The behaviour is real and NOT fixed on this branch — I reproduced all six claimed cases (plus an AB-failing-0/0 case) by driving the current `screen_trio` on my own synthetic trio, and confirmed Step 5's audit has no drop metric at all. Three corrections to the finding's framing:

1. LINE ATTRIBUTION. Only the UNKNOWN-…

**Fix.** In the `else` at :312, emit the het with `origin='unknown'` + `flags=parent_uninformative` (never-drop, as the sibling branch already does), or at minimum increment a per-trio `candidates_parent_uninformative` audit metric.

### D15. A genuine biallelic hit loses BOTH legs when one leg fails any single collector gate — the pairing model has no partial credit and no counter

**`pipeline/05_inheritance_screen.py:331`** · silent loss · GENOTYPE · verifiers: confirmed/partly/confirmed · reproduced by execution · documented at docs/limitations.md:288

A leg that fails child-het QC (:257), transmitting-parent QC (:305/:309) or the gene-key match (:258) never enters `hets`, so `pairs` at :331-332 is empty. The surviving leg then has only one exit: the dominant emission at :366, which re-gates at `dominant_max` (1e-4) although `hets` was pooled at `recessive_max` (1e-2). In the recessive band both legs therefore vanish with no row, no flag and no counter.

**Failure scenario.** Gene GATEA, both legs missense at faf95 5e-4: leg 1 child 0/1 AB 0.50 from mother; leg 2 child 0/1 AB 0.22 (39,11 — below het_ab_min 0.25) from father. Executed: ZERO rows for both. Gene GATEB, identical but the child is clean on both legs and the transmitting FATHER is AB 0.18 on leg 2: ZERO rows for both. Control GATEOK (same coordinates, same 5e-4, both legs clean): the `compound_het` pair is emitted — so the single failing gate is the whole cause. Control DOMFT (identical to GATEA but both legs at 1e-5) emits a `dominant` row for the survivor, proving the rarity band is the second, binding constraint. Correct behaviour: the child unambiguously carries two rare functional hets in one gene; at minimum the surviving leg should fall through to a reported row, and the dropped leg should be counted.

**Who it hits.** Any proband with two rare functional hets in one gene in the 1e-4..1e-2 band where either leg carries a QC blemish, a `both`-origin parent pair, or a split gene key. This is the pipeline's ONLY two-hit recessive discovery model. The documentation UNDERSTATES it: limitations.md:288 scopes the no-row outcome to the pick-key case ('damaging at EQUAL consequence rank in two overlapping MANE genes'), and the executed evidence shows child AB, transmitting-parent AB, and HET x HET parents each produce the identical total loss with no overlapping gene involved.

**What the verifiers corrected.**

- The finding is correct as written; its SCOPE is narrower than reality. Two refinements:

1. The trigger list is larger than the three cases reported. Executed and confirmed, each causing both legs to vanish with no row/flag/counter: child het-AB blemish (:257), child GQ below min_gq (:257), transmitting-parent AB blemi…
- The reproduced behaviour is real but misattributed. (1) Compound-het pairing causes none of the loss: the QC-failing leg is dropped by the uniform genotype-QC gate that applies to every mode — proven by chr4:1010 in the finder's own DOMFT control, which gets NO CALL at dominant-grade 1e-5 rarity — and the surviving cle…
- Three corrections to the finding as written. (1) It omits the highest-frequency limb: a single parent NO-CALL on one leg (origin=None, :313) produces the identical total loss of both legs — I reproduced it (chr3:8000/8010, zero rows). That limb matters more in real GMKF data than the origin=both case it does cite. (2)…

**Fix.** Emit an unpaired recessive-band leg as a reported row (e.g. mode `het_unpaired` with the failing-gate reason), or at minimum audit a per-trio `recessive_band_unpaired` counter.

### D16. A `both`-origin het (both parents 0/1) in the [dominant_max, recessive_max) band is emitted under no mode at all — the same unphaseable situation a de-novo partner is emitted-and-flagged for

**`pipeline/05_inheritance_screen.py:331`** · silent loss · GENOTYPE · verifiers: confirmed/refuted/partly · reproduced by execution · documented at docs/inheritance_and_genotype_qc.md:146,165 (documents that a `both` leg is never paired, but not that it then yields no row at all; docs/limitations.md:283-316 "Inheritance-model residuals" omits this case entirely while listing the exactly-parallel pick-key case in those words). The documentation UNDERSTATES the impact.

`by["both"]` is populated at :303/:316 and never read: `pairs` is built only from mat x (pat + denovo) and pat x denovo. A `both` leg is in the dominant allowlist at :366 but that emission also requires `rare(v, dominant_max)` (1e-4), so a `both` leg pooled at 1e-2 and sitting above 1e-4 is paired by nothing and emitted by nothing, with no counter.

**Failure scenario.** Consanguineous/founder trio, two HIGH-impact hets in one gene, faf95 3e-3 each, both parents 0/1 at AB 0.50 GQ 99 (the canonical carrier-couple configuration). Both legs get origin=both, neither pairs, neither clears the dominant gate: zero rows. The identical pair with one parent hom-ref instead emits a `compound_het` pair. It should at minimum be emitted with a `phase_unresolved_both_origin` flag, exactly as an equally-unphaseable de-novo-partner pair is emitted with `unphased_denovo_partner` (docs/inheritance_and_genotype_qc.md:168).

**Who it hits.** Any biallelic-candidate gene where the child's het is carried by both parents and rarity_af is in [1e-4, 1e-2) — enriched by consanguinity and founder alleles, i.e. the target cohorts.

**What the verifiers corrected.**

- Two refinements, neither of which weakens the finding.

(1) The precise defect is narrower and sharper than "any both-origin het in the band is lost". A LONE `both` het in that band is symmetric with a lone `mat`/`pat` het in that band — I confirmed by execution that GENEM's solitary `mat` leg at faf95=3e-3 also emits…
- `origin=both` is irrelevant to the outcome. ANY unpaired het in [dominant_max, recessive_max) emits no row — executed: a single `mat`-origin het at faf95=3e-3 with an affirmatively QC-passing hom-ref father (KID 0/1, DAD 0/0, MOM 0/1) produces NO ROW, identically to the `both`/`both` case. This is the behaviour docs/li…
- Three corrections to the finding as written. (1) The triggering condition is a `both`-origin leg WITH a same-gene partner, not "origin=both AND the frequency band": an unpaired band het yields no row regardless of origin (executed: single mat@3e-3 -> NO ROWS), and that drop is asserted as intended at tests/test_pure.py…

**Fix.** Pair `both` legs and flag `phase_unresolved_both_origin`, or at minimum audit the count.

### D17. A compound-het pair whose trans evidence FAILED or was VACUOUS still consumes both legs, deleting the dominant call — violating the rule stated in the code three lines above

**`pipeline/05_inheritance_screen.py:341`** · wrong call · GENOTYPE · verifiers: partly/partly/partly · reproduced by execution · **not documented anywhere**

`unphased` at :336 is defined ONLY as `ka in denovo_keys or kb in denovo_keys`. The unverified flags `ua`/`ub` (the non-transmitting parent was not an affirmative QC-passing 0/0) and the vacuous-pass flags `wa`/`wb` (that parent had no AD at all) are computed and used for flags at :347-354 but never enter the consume test at :341. So a pair the code itself marks as possibly-CIS suppresses the `dominant` row on both legs — contradicting the comment at :337-340, 'ONLY a phase-confirmed pair may veto the dominant model'.

**Failure scenario.** Gene GVETO, both legs stop_gained at faf95 1e-5 (dominant grade). Leg 1: child 0/1, mother 0/1, father ./. (no-call, so trans is unverified). Leg 2: child 0/1, father 0/1, mother affirmative 0/0. Executed: two `compound_het` rows flagged `origin_unverified`, `consumed` populated, and ZERO `dominant` rows. The identical file's de-novo-partner control (GDN) correctly does NOT consume and yields its `dominant` row; the vacuous-pass case (GVAC, father 0/0 with FORMAT/AD absent) also consumes and yields no dominant row. Correct behaviour: a pair carrying `origin_unverified` or `trans_evidence_unmeasured` is epistemically identical to the de-novo-partner pair — 50/50 cis/trans — so it must be emitted AND leave the dominant call standing. Since `n_dominant` is the PRIMARY headline recurrence signal (06_gene_burden.py:352, sorted on at :495), every parental no-call at a second hit silently removes a carrier from the headline tally.

**Who it hits.** Any comp-het pair in which the non-transmitting parent is a no-call, a low-GQ/low-DP 0/0, or a 0/0 with no FORMAT/AD (GATK ref-block shape). Parent no-calls run a few percent of sites in real trio WGS and ref-block 0/0 without AD is the NORMAL shape of a non-carrier parent under a child het, so `trans_evidence_unmeasured` in particular should fire on a large fraction of pairs.

**What the verifiers corrected.**

- The defect reproduces exactly as described — executed, not reasoned. Three corrections to the finding's framing:

(a) DOCUMENTATION STATUS IS NOT "NO". docs/limitations.md:305-313 ("A permissive comp-het partner can still suppress a dominant call") documents the suppression mechanism and names the precise downstream co…
- The mechanism is real and reproducible, but it is already documented as an open residual at docs/limitations.md:305-313, it is flagged in the output (`origin_unverified` / `trans_evidence_unmeasured`), and it removes NO carrier: 06_gene_burden.py:346 unions dom|bi|x into `n_carriers` and :443 flags `recurrent` on that…
- The mechanism is real and reproduces on current HEAD; it is neither documented nor fixed on this branch. But the impact claims need correcting: (1) the variant is NOT silently lost — both legs are emitted as compound_het rows carrying `origin_unverified` / `trans_evidence_unmeasured`, so the deficiency is visible in th…

**Fix.** Change :341 to `if not (unphased or ua or ub or wa or wb):` so only a genuinely phase-confirmed, measured pair consumes its legs.

### D18. The genotype-QC allele-balance knobs have no unit or range check, and no step guards a zero-call run: percent-scale values reduce Step 5 to 0 calls with exit 0 everywhere

**`src/hprv/genotype.py:40`** · silent loss · GENOTYPE · verifiers: partly/partly/partly · reproduced by execution · **not documented anywhere**

`GtThresholds.from_config` coerces `het_ab_min/het_ab_max/homalt_ab_min/homref_ab_max` with bare `float()` and never checks they lie in [0,1]. A value on a percent scale (25/75/90/10 instead of 0.25/0.75/0.90/0.10) makes `thr.het_ab_min <= ab <= thr.het_ab_max` unsatisfiable for every het and `ab >= homalt_ab_min` unsatisfiable for every hom-alt, while `ab <= homref_ab_max` becomes vacuously true. Step 5 then emits nothing and `main()` returns 0 (05_inheritance_screen.py:491) — there is no analogue of the guard `04_subset_and_annotate_trios.sh:123` already applies to the identical failure one step earlier.

**Failure scenario.** A user expresses allele balance as a percentage — one coherent, self-consistent edit of the whole `filters.genotype_qc` block. Executed end-to-end with the real Step-5 CLI on the shipped integration manifest: `Step 5 complete: 0 candidate calls across 2 trios`, EXIT=0, candidates.calls.tsv has 0 rows, audit records `candidate_calls 0` per trio and `candidate_calls_total 0`. Step 6 then reports `0 genes, 2 trios` and also exits 0; the workbook, igv/variants.tsv and Step 9 all complete empty. Nothing distinguishes this from a cohort in which no trio carried a qualifying variant. Partial versions are worse because they are not empty: `homalt_ab_min: 90` alone silently removes every hom_recessive call (2 -> 0) and two compound_het legs (14 -> 12) while the run still looks populated.

**Who it hits.** Everything. With all four knobs on a percent scale the entire call set is lost; with `homalt_ab_min` alone every homozygous-recessive and hemizygous X-linked call, plus every compound-het leg whose carrier parent is 1/1.

**What the verifiers corrected.**

- Two corrections. (1) Step 9 is NOT silent on an empty call set: pipeline/09_prioritize.py:449-454 refuses and returns 1 ("has no data rows — refusing to write an empty prioritization"), executed and confirmed (EXIT=1, no files written). Since run_pipeline.sh is `set -euo pipefail` with default TO=9, a default full run…
- Three corrections. (a) CLASS: this is not silent_loss — Step 5 emits "Step 5 complete: 0 candidate calls across 2 trios / by mode: {}" to stderr and records audit metrics candidate_calls=0 per trio and candidate_calls_total=0. The loss IS counted; what is absent is a guard distinguishing "zero because misconfigured" fr…
- Three specifics in the finding are off, all in the direction of understating or mislabelling:

1. The baseline is wrong because the reporter ran Step 5 WITHOUT `--qc-report`. Their baseline is 39 rows with no x_linked_recessive (CH_A's X modes were skipped: "child sex unresolved"). Run the way pipeline/run_pipeline.sh:…

**Fix.** Range-check the four AB knobs to [0,1] at `GtThresholds.from_config` and add a zero-call guard in Step 5 mirroring 04_subset_and_annotate_trios.sh:123.

### D19. A trio VCF missing FORMAT/AD or FORMAT/GQ yields ZERO calls for that trio, silently; Step 4 warns only for AD and says calls 'may be dropped' when in fact every call is dropped

**`src/hprv/genotype.py:129`** · silent loss · GENOTYPE · verifiers: confirmed/confirmed/confirmed · reproduced by execution · documented at PARTIAL — pipeline/04_subset_and_annotate_trios.sh:207 warns 'no FORMAT/AD in the trio VCF — allele-balance QC is unavailable ... so such calls may be dropped'. That understates it (every call is dropped, not 'such calls'), it is a WARN with no audit metric, and there is no equivalent check or warning for FORMAT/GQ anywhere.

`sample_qc` returns False immediately when `gq()` is None (:131-133), and the `het` / `hom_alt` / `denovo_child` limbs require `ab is not None` (:140-143). Every Step-5 mode requires the CHILD to pass one of those three limbs, so an absent FORMAT/GQ or FORMAT/AD makes every mode unreachable for the whole trio. Step 5 then records `candidate_calls 0` — indistinguishable from 'screened and found nothing'.

**Failure scenario.** A trio VCF carrying `GT:DP:GQ` (no AD) — e.g. a re-delivered/legacy callset, an externally supplied parent VCF, or any file that has been through `bcftools annotate -x FORMAT/AD`. Three textbook calls (a dominant het, a hom-recessive, a de novo) all produce 0 rows. Same total loss for `GT:AD:DP` (no GQ), for which NOTHING in the pipeline warns at all — `grep -rn 'FORMAT/GQ' pipeline/ src/` returns nothing. The run exits 0 and audit/counts.tsv shows `05_inheritance <trio> candidate_calls 0`.

**Who it hits.** ALL variants in an affected trio — a 100% loss for that trio, not a class-specific one. Reachable for any trio whose VCF provenance differs from the GATK genotype-refinement default.

**What the verifiers corrected.**

- Three corrections to the finding as written, none of which change the verdict. (1) Documentation status: the requirement is stated at docs/inheritance_and_genotype_qc.md:137 — "Genotype-level rules assume a mother-father-child trio VCF on GRCh38 with GT, GQ, DP, AD, and PL/PP available." The finding cites only 04_subse…
- The finding stands as written; two nuances soften the `silent_loss` class label. (1) It is not literally counter-free: Step 5 records `candidate_calls 0` per trio and Step 4 records `candidate_genotypes N > 0`, so audit/summary.md does show the funnel collapsing — it is indistinguishable in KIND from a genuinely empty…
- Two refinements, neither of which overturns the finding.

(a) The AD half of the MECHANISM is already documented — CLAUDE.md:486-489 states it verbatim ("`het`/`hom_alt`/`denovo_child` require `ab is not None` and so DROP a carrier when FORMAT/AD is absent"), repeated in the `sample_qc_ad_measured` docstring at src/hpr…

**Fix.** Add a Step-5 (or Step-4 preflight) hard check that the trio VCF header declares FORMAT/GQ and FORMAT/AD, and `die` — or at minimum audit `sample_qc_unavailable` per trio so a zero call count is distinguishable from a genuine negative.

### D20. A callset with no FORMAT/AD (or no FORMAT/GQ) produces zero calls in every mode, exit 0, with a Step-4 audit identical to a healthy run

**`src/hprv/genotype.py:129`** · fragility · GENOTYPE · verifiers: confirmed/partly/partly · reproduced by execution · documented at PARTIAL - CLAUDE.md:488 documents only the per-SAMPLE missing-AD asymmetry (a ref-block parent); pipeline/04_subset_and_annotate_trios.sh:207 warns but understates the scope to multiallelic hets; GQ is documented nowhere. The whole-callset zeroing is UNDOCUMENTED and the existing warning text materially understates it.

sample_qc (genotype.py:129-150) tests GQ first and allele balance last, and every limb that can AFFIRM a carrier -- het / hom_alt / denovo_child -- requires `ab is not None`. allele_balance() returns None whenever FORMAT/AD is absent and gq() returns None whenever FORMAT/GQ is absent, so a callset lacking either field fails every carrier limb for every sample at every site and Step 5 emits nothing under any mode. The only guard in the whole pipeline is a single WARN at pipeline/04_subset_and_annotate_trios.sh:207 whose text scopes the damage to 'a multiallelic (1/2) het ... may be dropped'; there is no guard of any kind for GQ.

**Failure scenario.** A trio callset delivered without FORMAT/AD (a non-GATK caller such as freebayes, which emits RO/AO; or any VCF that has been through `bcftools annotate -x FORMAT/...` upstream). Executed on the shipped integration cohort: baseline is 41 calls (dominant 20, compound_het 14, denovo 3, hom_recessive 2, x_linked_recessive 2). After stripping ONLY FORMAT/AD from the two raw trio VCFs -- leaving GT:DP:GQ = 0/1:40:99, i.e. perfect depth and quality -- Step 4 still audits `candidate_genotypes 31` and `6`, byte-identical to the healthy run, emits one WARN per trio, and Step 5 then reports `0 candidate calls across 2 trios`, `by mode: {}`, EXIT=0. Stripping ONLY FORMAT/GQ instead gives the identical 41 -> 0 with no warning emitted anywhere. What it should do: halt at preflight when the trio VCF header declares neither field, or at minimum audit a per-trio `genotypes_qc_unmeasurable` count so that a zero is distinguishable from a clean negative.

**Who it hits.** Every variant in every trio, all five modes. Not triggered by GMKF GATK output as delivered; triggered by any callset lacking AD or GQ cohort-wide (non-GATK caller, re-headered or FORMAT-stripped delivery, collaborator-supplied VCF).

**What the verifiers corrected.**

- Two corrections, neither changing the verdict.

1. "GQ is documented nowhere" is overstated. docs/inheritance_and_genotype_qc.md:137 explicitly states
   the assumption: "Genotype-level rules assume a mother-father-child trio VCF on GRCh38 with GT, GQ,
   DP, AD, and PL/PP available." What is undocumented and unenforce…
- Three corrections. (a) The freebayes trigger is refuted by execution: cyvcf2 0.31.4 falls back to FORMAT/RO+AO, so a freebayes record (GT:DP:GQ:RO:AO = 0/1:40:99:20:20) yields ref_ad 20 / alt_ad 20 / ab 0.5 and sample_qc(...,'het') -> True. Real triggers are a FORMAT-stripped or re-headered delivery, or a plain `bcftoo…
- The finding is mechanically accurate and reproduces on current HEAD; three adjustments. (1) Documentation citation is incomplete: besides CLAUDE.md:487-489 (which states the drop-on-missing-AD mechanism verbatim), docs/inheritance_and_genotype_qc.md:137 explicitly declares "GT, GQ, DP, AD, and PL/PP available" as an in…

**Fix.** Hard-stop at Step 4 / preflight when the trio VCF header declares no FORMAT/AD or no FORMAT/GQ, and audit a per-trio `genotypes_qc_unmeasurable` counter.

### D21. `clean_parent`'s alt-AD cap is an ABSOLUTE count, so its stringency inverts with depth — and a failure removes the child's variant under EVERY mode, not just de novo

**`src/hprv/genotype.py:149`** · silent loss · GENOTYPE · verifiers: partly/partly/partly · reproduced by execution · documented at PARTIAL — docs/inheritance_and_genotype_qc.md:187 states the rule ("each parent alt AD <= 1") but justifies it with a FRACTION rationale it does not implement ("a parental alt fraction of a few percent suggests inherited or parental mosaicism, not de novo") and implies re-classification ("not de novo") where the code deletes the variant entirely; the depth inversion, the mosaic-parent class and the comp-het erasure are undocumented, and docs/limitations.md:283-317 does not list this residual

`sample_qc(..., 'clean_parent')` returns `(a is None or a <= thr.parent_max_alt_ad) and (ab is None or ab <= homref_ab_max)` with `parent_max_alt_ad = 1` (genotype.py:29, config.example.yaml:342). Because the limb counts READS rather than a fraction, a deeper parent is judged more harshly. And the consequence is not confined to the de novo row: at 05_inheritance_screen.py:288-289, when both parents read HOM_REF the het collector sets `origin = "denovo" if both clean_parent else None`, and `origin=None` means the variant is never appended to `hets` (:314-316) — so it also cannot become a `dominant` call or a compound-het leg. One noise read deletes the variant from the run.

**Failure scenario.** EXECUTED, two ways. (1) Depth sweep over a hom-ref father at a site where the child is a clean 0/1: DP 30 / 1 alt (3.33% alt fraction) PASSES clean_parent, while DP 1000 / 3 alt (0.30%) FAILS — a ten-fold cleaner parent is rejected. Every failing row in the sweep also passes the sibling predicate `sample_qc(...,'hom_ref')`, which is what the code uses when the OTHER parent carries, so the same parent is 'clean enough' in one branch and disqualifying in the other. (2) End-to-end through the real `screen_trio`: GDN7 (stop_gained, HIGH, no gnomAD record; child 0/1 AD=150,150 GQ=99; father 0/0 AD=298,2 — 0.67% alt) and GDN8 (a 8%-VAF mosaic father, AD=92,8 DP=100) each produce NO ROW UNDER ANY MODE and no counter, while GDN9 (father AD=29,1, a FIVE-FOLD higher alt fraction) is emitted as `denovo`. The comp-het layer is corrupted by the same mechanism: gene GCH2 with a maternal het at chr2:12000 and a de novo at chr2:12500, both faf95=5e-4, yields a `compound_het` pair when the father is clean; giving the father 2 alt reads out of 200 at the de novo site deletes BOTH legs — no pair (the de novo leg never enters `hets`) and no dominant row for the maternal leg (5e-4 sits in the [1e-4, 1e-2) inert band). A single 1%-noise site erases an entire biallelic candidate with zero trace.

**Who it hits.** Every de novo whose parent carries >=2 alt reads: parental gonosomal mosaics (high clinical value — recurrence risk), index-hopped/contaminated libraries, and systematically noisy loci. Rate scales with sequencing depth, so it worsens on deeper or PCR-free callsets.

**What the verifiers corrected.**

- The mechanism, the depth inversion, the total absence of any counter, and the hom_ref/clean_parent asymmetry all reproduce exactly as claimed, on the current branch, via the real functions. Three corrections to the write-up: (1) "cannot become a dominant call" is inert — in the both-parents-hom-ref branch that produces…
- Three of the finding's load-bearing claims are wrong.

(1) "Removes the child's variant under EVERY mode, not just de novo" is FALSE. `clean_parent` has exactly two call sites (`05_inheritance_screen.py:194/201` and `:288-289`), both in the de-novo path. Executed: both parents HOM_REF, no gnomAD record (rarest), clean…
- Three corrections.

1. NOT NEW — already documented AND explicitly declared a non-defect on this branch. docs/pipeline_review_2026-09.md:262 ("`parent_max_alt_ad = 1` is an absolute count; 1 alt read in 10 and 1 in 200 pass identically. A fraction-based cleanliness test is the usual refinement.") sits under "Observatio…

**Fix.** Make the limb depth-aware — `alt_ad <= max(parent_max_alt_ad, ceil(frac * dp))` with a configurable `parent_max_alt_frac` (~0.02) — and never let a clean_parent failure suppress the het collector: keep `origin=None` for the de novo row but still admit the het with `flags=parental_alt_reads=N` so it can reach the dominant/comp-het layers.

### D22. Every Step-5 evidence flag except `origin=` is dropped from the reviewer-facing tables — an unphased compound-het reads as a confirmed biallelic hit

**`src/hprv/igv.py:193`** · unverifiable · GENOTYPE · verifiers: confirmed/partly/partly · reproduced by execution · documented at NO — docs/artifact_gene_triage_review.md:288 obliquely notes its phase check 'needs candidates.calls.tsv', but nothing states that igv/variants.tsv and both prioritized tables carry no flags column; docs/prioritization.md:712 documents a different gap (parental GQ/DP/AB absent)

`build_variants_tsv` reads `r.get("flags")` solely to feed `_origin()` (:67-72), which extracts the `origin=` token and DISCARDS every other token; `COLUMNS` (:31-64) has no `flags` entry, so `igv/variants.tsv` has no flags column. Step 9 then loses them on both output paths: `variants.prioritized.tsv` is written from the curated `VARIANT_COLUMNS` (pipeline/09_prioritize.py:84, :953) which omits `flags` even when the input WAS `candidates.calls.tsv`, and `igv/variants.prioritized.tsv` is the input header plus its complement (:982-984), so it can only carry what Step 8 wrote. Result: `unphased_denovo_partner`, `origin_unverified`, `trans_evidence_unmeasured`, `parent_ad_unmeasured`, `father_carries_x_allele` and `high_conf_rarity` exist only in `candidates.calls.tsv` and the workbook sheet that copies it — never in the table CLAUDE.md:241 calls "the table a reviewer actually opens".

**Failure scenario.** A compound-het pair whose second hit is a de novo. Step 5 emits both legs with `flags=unphased_denovo_partner` precisely so, per docs/inheritance_and_genotype_qc.md:168 and CLAUDE.md:187-189, the pair 'is not read as a confirmed biallelic hit' (trio genotypes cannot phase a de novo against an inherited allele; it is ~50/50 cis/trans, and WhatsHap is not wired). In `igv/variants.prioritized.tsv` that pair renders as `inheritance=compound_het`, `pair_id=CH_A:CH1`, `origin=` (blank), with no phase qualifier anywhere in the 100+ columns — indistinguishable from a mat x pat pair that IS trans by descent. Step 9 adds `partner_leg_quality_unknown` (parental GQ/DP/AB absent) but nothing about phase. The same erasure hides `trans_evidence_unmeasured`, the witness added specifically so a vacuous `clean_parent` pass (a GATK ref-block 0/0 parent with no FORMAT/AD) is distinguishable from a measured one, and `origin_unverified`, the marker that the non-transmitting parent was never affirmatively observed hom-ref — i.e. exactly the three markers that separate an inferred phase from a confirmed one, all invisible in the file the phase claim is reviewed from.

**Who it hits.** Every compound_het pair with a de novo leg or an unverified non-transmitting parent, every call resting on a vacuous AD-unmeasured parental pass, and every recessive/X-linked call otherwise tiered by high_conf_rarity. 17/41 (41%) of the shipped run's calls; on real data the comp-het share is quadratic in a gene's het count, so the erasure concentrates in exactly the long genes where phase is least trustworthy.

**What the verifiers corrected.**

- Two refinements, both narrowing rather than contradicting the finding:

(1) The 17/41 figure is right but conflates two kinds of erasure. Only 10/41 (24%) carry an EVIDENCE-critical erased flag (`unphased_denovo_partner` 8, `parent_ad_unmeasured` 1, `father_carries_x_allele` 1). The other 7 are `high_conf_rarity` alone…
- The mechanism (no `flags` column in igv/variants.tsv, variants.prioritized.tsv or igv/variants.prioritized.tsv) is confirmed by execution, but two claims need correcting. (1) An unphased de-novo comp-het pair is NOT strictly indistinguishable from a trans-by-descent pair in the review table: `mother_gt` and `father_gt`…
- The defect is REAL, CURRENT and UNDOCUMENTED — not fixed on this branch and absent from the resolution table — but two load-bearing claims are overstated, and both prop up the "major" rating.

1. "Indistinguishable from a mat x pat pair that IS trans by descent" is wrong. It is distinguishable, just not FILTERABLE. The…

**Fix.** Add `flags` to `igv.COLUMNS` (keeping the derived `origin` column) and to `09_prioritize.VARIANT_COLUMNS`, so the phase/QC qualifiers travel with the call into every reviewer-facing table.

### D23. The SpliceAI keep-path reads only the selected CSQ block, so a scored splice variant on a non-picked overlapping gene is dropped as `not_functional`

**`pipeline/02_annotate_sites.sh:516`** · silent loss · FILTER · verifiers: confirmed/confirmed/partly · reproduced by execution · **not documented anywhere**

The SpliceAI VEP plugin emits per-GENE delta scores (which is why `SpliceAI_pred_SYMBOL` is in Step 2's want list at :418), but `bcftools +split-vep -c "$have_fields" -s "$sel"` lifts ONE block's values for every field alike. Under the auto-resolved default `-s pick` (:471, always taken on a `--flag_pick` VCF), the delta score belonging to any non-picked gene is discarded before `annotations.spliceai_ds()` (annotations.py:392) ever sees it, and `selection.py:50` then has nothing to rescue with. The counter-selectors do not fix it either: `-s mane`/`-s all` recover the score but emit `vep_IMPACT` as `HIGH,MODIFIER` at a two-block locus, which is not in `keep_impacts` — so no configured value of `resources.vep.csq_select` keeps both the SpliceAI rung and the IMPACT rung intact at an overlapping-gene locus, and both failures are filed under the same reason string.

**Failure scenario.** A deep-intronic variant in an intron of GENEA (SpliceAI DS_DL 0.87, scored by the plugin under symbol GENEA) that also lies in the 3'UTR of the overlapping GENEB, whose transcript VEP PICKs (the shipped --pick_order at docs/functional_annotation.md:220 leads with mane_select, so the MANE gene wins regardless of which gene the splice site belongs to). EXECUTED on bcftools 1.22 with exactly that two-block CSQ, then through the real getters and the real Step-3 classifier: `-s pick` and `-s worst` both yield `vep_SpliceAI_pred_DS_DL=.`, `spliceai_ds()` -> None, `classify()` -> (False, 'not_functional'); `-s all` and `-s mane` yield `.,0.87`, `spliceai_ds()` -> 0.87, `classify()` -> (True, 'spliceai'). What should happen: the score is a property of the variant/gene pair the pipeline already lifts a SYMBOL for, so the variant should be kept as `spliceai` under every selector. What actually happens: it is dropped under the reason that means 'we looked and found no functional signal', and downstream `prioritize.spliceai_status` (:913-921) would label it `not_covered`, asserting the precomputed set has no score for it — which is false, the score exists in the SpliceAI file.

**Who it hits.** Variants VEP rates below MODERATE (deep-intronic, exonic-synonymous, UTR) at overlapping / nested / readthrough gene loci, where SpliceAI is the ONLY keep-path. I could not run real VEP here, so I have not measured how often PICK lands on the gene without the SpliceAI record; the mechanism and both endpoints are demonstrated, the frequency is not. Note docs/pipeline_design.md:221 frames `csq_select: mane` as costing recall, which is backwards for this rung, and docs/functional_annotation.md:60 states the opposite requirement ('retain per-transcript annotations so a consequence on a clinically relevant non-canonical transcript is not lost').

**What the verifiers corrected.**

- Four corrections; none defeats the finding, and the first strengthens it.

1. `-s mane` does NOT reliably recover the score — the finding's claim that mane yields `.,0.87` is construction-dependent. When only the PICKed gene's transcript carries MANE_SELECT (the realistic readthrough/nested case, and the one the shippe…
- Three corrections, none of which rescues the code.

(a) "-s mane recovers the score" is an artifact of the finder's mock, which set MANE_SELECT on BOTH blocks. Executed with MANE_SELECT on the picked gene only — the realistic shape where pick and mane agree — `-s mane` behaves exactly like `pick`: SYM=GENEB, DS_DL=. So…
- Two detail corrections, both verified by execution; the core mechanism stands.

(a) The `-s mane` result is CONFIGURATION-DEPENDENT, not a general property. Reproduced all three cases: when BOTH overlapping genes carry their own MANE Select transcript (the realistic two-overlapping-MANE-genes case, and the one the repo…

**Fix.** Lift the four SpliceAI DS fields (and any other per-transcript plugin score) in a second `+split-vep -s all` pass so `_max_float` takes the max over blocks, and keep `$sel` only for the identity fields (SYMBOL/Gene/Feature/IMPACT/Consequence).

### D24. Step 3 hard-fails when the rarity field is undeclared but performs no equivalent check for the two functional-rung fields; on the ingest and `--from 3` paths a dead SpliceAI produces no output at all

**`pipeline/03_select_plausible.py:46`** · silent loss · FILTER · verifiers: confirmed/partly/partly · reproduced by execution · documented at CLAUDE.md:648 (the preflight halt is documented as not enforced on the annotated_vcf path) — but it understates: nothing else detects it either, and the SpliceAI warning at pipeline/02_annotate_sites.sh:454 is conditioned on the resource being configured, so that path is silent end to end.

`03_select_plausible.py:46-55` returns 1 when `INFO/gnomad_AF_joint` is absent from the header, on the stated reasoning that a silently-absent field would make every gate read None while the run exits 0. The two fields that decide 29 of VEP's 41 consequence terms — `vep_CADD_PHRED` and `vep_SpliceAI_pred_DS_*` — get no such check. Upstream, `run_pipeline.sh:119-121` skips the whole CADD/SpliceAI preflight when `resources.vep.annotated_vcf` is set (and `run_step 2` at :104 skips it entirely on `--from 3`), and `02_annotate_sites.sh:454` gates its SpliceAI warning on `is_set HPRV_SPLICEAI_SNV` — which is exactly what is NOT set on an ingest run. CADD at least gets an unconditional warn at :447; SpliceAI gets none.

**Failure scenario.** A group annotates centrally with a VEP cache only (no CADD, no SpliceAI plugin) and the run sets `resources.vep.annotated_vcf`. Executed with the real Step-3 CLI and the shipped `config/config.example.yaml` on a VCF whose header declares neither field: `Step 3 complete: 2/6 sites plausible`, EXIT=0, audit `reason.not_functional 4`. The dropped four are `synonymous_variant`, `intron_variant`, `5_prime_UTR_variant`, `splice_region_variant` — the whole below-MODERATE class. `spliceai_required: true` is the shipped default and CLAUDE.md:316 calls it a preflight HALT; here it is inert and emits nothing. Detection is doubly blocked: the only audit observable that a live SpliceAI leaves is `reason.spliceai`, and `03_select_plausible.py:69,83` only emits reason rows for reasons that occurred, so 'SpliceAI is dead' and 'SpliceAI rescued nothing' produce byte-identical `counts.tsv`. Expected: the same header guard the rarity oracle gets, or at minimum an unconditional warn for SpliceAI to match CADD's.

**Who it hits.** All 29 of VEP's 41 consequence terms that sit below MODERATE — intronic, synonymous, UTR, splice-region, regulatory, non-coding-transcript. Reachable on any `resources.vep.annotated_vcf` run and on any `--from 3` re-run against an older annotated union.

**What the verifiers corrected.**

- Three refinements, none of which weakens the core claim:

1. The dropped variants ARE counted — as `reason.not_functional` (4 in my run). So this is not "no counter"; it is a MISATTRIBUTED reason: the sites are recorded as "not functional" when the truth is "not scored". The finding's own text states this correctly, bu…
- The residual true statement is narrow: on the `resources.vep.annotated_vcf` ingest path, SpliceAI's absence produces no warning, because 02_annotate_sites.sh:454 conditions the warn on `is_set HPRV_SPLICEAI_SNV` while CADD's warn at :447 is unconditional. The one configuration that is genuinely unannounced is CADD PRES…
- Three corrections to the finding as written. (1) It is not undocumented: the ingest-path exemption is stated in config/config.example.yaml:103, docs/README.md:117 and :125-126, docs/functional_annotation.md:157 and CLAUDE.md:648. What is undocumented is that the path emits no signal whatsoever for SpliceAI — no warn, n…

**Fix.** Add a Step-3 header/value guard for `vep_CADD_PHRED` and `vep_SpliceAI_pred_DS_AG` mirroring the `gnomad_AF_joint` guard, and always emit every reason key (including zeros) to `audit/counts.tsv`.

### D25. Step 5 accepts a per-trio VCF carrying no oracle fields at all and reports the run as faf95; Step 3 refuses the identical input

**`pipeline/05_inheritance_screen.py:427`** · wrong call · FILTER · verifiers: confirmed/partly/partly · reproduced by execution · **not documented anywhere**

03_select_plausible.py:41-55 hard-stops (`return 1`) when `oracle: faf95` and the input header declares no INFO/gnomad_AF_joint, because "every variant would read as ABSENT from gnomAD (= rarest) and no rarity gate would fire". Step 5 opens its input at :427 with no equivalent check — `Trio.__init__` inspects the header only for `ID=hiConfDeNovo` (:74) — so the same file is screened, every `A.frequency()` returns None, `rare(v, limit)` is True for every allele and every gate passes, while `rarity_oracle` is stamped `faf95` on all 41+ rows and audited as `rarity_oracle.faf95 1`. Reachable without any file corruption: run_pipeline.sh's preflight (:110-116) only checks that the slim FILE exists, and Step 4's per-trio cache key is `plausible-cksum|vcf|stat|samples` — nothing about the oracle — so `--from 5` (or `--from 4`, which cache-hits) after switching `resources.gnomad.oracle` from grpmax_proxy to faf95 screens candidate VCFs built without the slim.

**Failure scenario.** Take the real CH_A candidates VCF, strip the five gnomad_* INFO fields (values and header, as a run without the slim would have) and leave every vep_* cache field intact; screen it under the shipped default `oracle: faf95`. Step 3 on that exact file exits 1 with the ERROR text. Step 5 exits 0 and emits 37 calls versus the healthy 35 — the list grew. The two extra rows are `dominant chr2:8000 GENE7` (true gnomad_faf95=0.00016, ClinVar pathogenic, admitted at Step 3 by the clinvar_plp override) and `dominant chr2:17500 GENEMID` (true gnomad_faf95=0.0016), both above `dominant_max`=1e-4 and both promoted because the lost oracle made them read "absent = rarest". On real data this is the whole BA1/common tail entering the candidate list as dominant calls while candidates.calls.tsv, igv/variants.tsv and the workbook all assert the run was gated on faf95. Step 5 should apply Step 3's header guard to every per-trio VCF it opens.

**What the verifiers corrected.**

- Two corrections. (1) The run is not silent in the finding's own scenario: Step 5's existing join-coverage guard (pipeline/05_inheritance_screen.py:468-478) emitted "WARN: 29 of 37 calls have NO gnomAD joint record yet DO carry a VEP-cache gnomAD AF" and recorded `rarity_faf95_absent_but_cache_has_af 29` in audit/counts…
- Step 5 does have a guard for this exact failure mode, one line below the code the finding cites as absent: `pipeline/05_inheritance_screen.py:461-478` counts calls with `rarity_basis == "absent"` that still carry a cache `grpmax_af`, records `rarity_faf95_absent_but_cache_has_af` in the audit, and writes a loud WARN. I…
- Step 5 DOES have a faf95 guard — pipeline/05_inheritance_screen.py:468-477 — and it fires in this exact scenario: my run printed "WARN: 29 of 37 calls have NO gnomAD joint record yet DO carry a VEP-cache gnomAD AF ... Those calls read as RAREST and may be common alleles" and recorded `rarity_faf95_absent_but_cache_has_…

**Fix.** Move Step 3's gnomad_AF_joint header guard into Step 5's per-trio VCF open (05:427), and halt when rarity_basis is 'absent' for ~100% of rows under oracle faf95.

### D26. A YAML scalar (or lower-case) `keep_impacts` silently deletes the entire impact rung — every stop_gained and missense is then dropped as `not_functional`

**`src/hprv/selection.py:39`** · silent loss · FILTER · verifiers: confirmed/partly/partly · reproduced by execution · **not documented anywhere**

`keep_impacts = set(get(cfg, "filters.functional.keep_impacts", ["HIGH","MODERATE"]))` applies `set()` to whatever YAML parsed. A scalar string becomes a set of CHARACTERS, so the membership test at :42 (`A.impact(v) in keep_impacts`) can never match a VEP IMPACT value; a lower-case list fails the same test on case. There is no schema validation anywhere in `src/hprv/config.py` (no `validate`, no allowed-key list), no default-shape assertion, and `tests/test_pure.py:411` exercises only the correct list shape.

**Failure scenario.** A user narrows the screen with `filters.functional.keep_impacts: HIGH` (a natural YAML slip for `[HIGH]`). Executed with the real `build_classifier`: `keep_impacts` parses to `{'G','H','I'}`; a rare `stop_gained` with no CADD score returns `(False,'not_functional')`; a rare `missense_variant` with CADD 22.0 returns `(False,'not_functional')`; only a `frameshift_variant` with CADD 33.0 survives, and then under `reason.cadd`. The impact rung is gone: the screen becomes CADD/SpliceAI-only, every HIGH-impact pLoF below CADD 25.3 is lost, the audit shows a large, plausible `reason.not_functional`, and the run exits 0. `keep_impacts: [high, moderate]` behaves identically. Expected: an unparseable/empty/unknown `keep_impacts` should be a hard stop, as `03_select_plausible.py:46-55` already does for the rarity witness field.

**Who it hits.** Every HIGH- and MODERATE-impact variant in the run whose CADD is below 25.3 and whose SpliceAI delta is below 0.2 — i.e. the large majority of pLoF and essentially all missense. Fires only on a mis-shaped config, but with no error and no distinguishing counter.

**What the verifiers corrected.**

- No correction needed — every claim reproduced. Two refinements: (1) the finding's evidence prints "keep_impacts-as-parsed=['G','H','I']" for the scalar case, which is the sorted SET; the raw parse is the string 'HIGH' (cosmetic, the mechanism is stated correctly). (2) The finding stops at the pure function; I additiona…
- The mechanism is exactly as described and I reproduced every variant of it, including two the finder did not list (`keep_impacts: HIGH, MODERATE` without brackets, and `keep_impacts: []`), plus the safe cases (shipped list, key omitted). No guard exists: no schema validation in config.py, only one `set(get(...))` in th…
- The mechanism, the executed evidence, the "not documented" claim and the "not fixed" status are all correct — src/hprv/selection.py:39 is unchanged from main and no schema validation exists in src/hprv/config.py. Two corrections. (a) Severity major -> moderate: the trigger is a user config-authoring slip only. config/c…

**Fix.** In `build_classifier`, reject a non-list `keep_impacts` and any member not in {HIGH,MODERATE,LOW,MODIFIER} with a `die`, mirroring Step 3's existing header guard.

### D27. `filters.functional.keep_impacts` written as a YAML scalar silently disables the entire impact keep-path and files every HIGH-impact loss under `reason.not_functional`

**`src/hprv/selection.py:39`** · silent loss · FILTER · verifiers: confirmed/partly/confirmed · reproduced by execution · **not documented anywhere**

Line 39 is `keep_impacts = set(get(cfg, "filters.functional.keep_impacts", ["HIGH","MODERATE"]))`. `set()` over a string yields a set of single CHARACTERS, so the membership test at line 42 (`A.impact(v) in keep_impacts`) can never match 'HIGH' or 'MODERATE'. There is no type check and no warning; the variant falls through the SpliceAI and CADD rungs and is returned as `(False, 'not_functional')` at line 78.

**Failure scenario.** A user narrows the screen by editing config.yaml from `keep_impacts: [HIGH, MODERATE]` to `keep_impacts: HIGH` (or to the quoted string `"HIGH,MODERATE"`). Every `stop_gained` / `frameshift_variant` / `splice_donor_variant` without a SpliceAI score >= 0.2 or a CADD >= 25.3 is dropped at Step 3. The run exits 0, `audit/counts.tsv` reports a plausible-looking `reason.not_functional` count, and `summary.md` shows a normal funnel — a reviewer reading a negative result cannot tell the primary keep-path was off. This is the worst finding class in the brief: an unrecoverable loss counted under a reason that means the opposite of what happened.

**Who it hits.** All HIGH and MODERATE impact variants lacking a rescuing SpliceAI/CADD score — i.e. the bulk of the screen's yield. Requires a config edit away from config/config.example.yaml:299, but that edit ('only keep HIGH') is a natural one and produces no error.

**What the verifiers corrected.**

- The finding is correct as written but UNDERSTATES the trigger surface in one way and slightly overstates it in another.

Understates: the bug is not confined to YAML scalars. A genuine YAML LIST in lowercase — `keep_impacts: [high, moderate]` — disables the impact rung identically, because selection.py:42 compares the…
- The mechanism is real and I reproduced it, but the classification and severity are wrong. It is `fragility` (a config combination), not `silent_loss`: Step 3 records the reason histogram, and in the broken state `reason.impact_high` and `reason.impact_moderate` vanish from audit/counts.tsv entirely while `not_functiona…
- Two adjustments, neither of which undermines the finding. (a) Severity: moderate, not major — the trigger is a config edit only; the shipped default config/config.example.yaml:299 is a proper YAML list and no data shape can trigger it, so a default GMKF-style run is unaffected. (b) The finding scopes the trigger too na…

**Fix.** Coerce and validate: reject a non-list `keep_impacts` (or split a string on commas) and `die` on any token outside VEP's four IMPACT values.

### D28. `filters.functional.keep_impacts` is unvalidated: four plausible YAML shapes silently disable the entire impact rung, and the losses are filed as `not_functional`

**`src/hprv/selection.py:39`** · silent loss · FILTER · verifiers: confirmed/partly/partly · reproduced by execution · **not documented anywhere**

`keep_impacts = set(get(cfg, "filters.functional.keep_impacts", ["HIGH","MODERATE"]))` then `A.impact(v) in keep_impacts`. A YAML scalar (`keep_impacts: HIGH`) becomes `set("HIGH") == {'H','I','G'}`; a title-case or lowercase list never matches VEP's uppercase IMPACT; an empty list matches nothing. Nothing anywhere validates the tokens, the type, or the case, and the drop is recorded under `reason.not_functional` — a reason that means the opposite of what happened.

**Failure scenario.** A user hand-edits config.yaml to `filters.functional.keep_impacts: HIGH` (or `[High, Moderate]`, the spelling used in prose). Step 3's first and dominant rung is now dead: every stop_gained, frameshift, canonical-splice and missense falls through to SpliceAI and CADD. Executed on the real integration annotated union (43 sites): plausible sites fall 36 -> 7 and `reason.not_functional` rises 2 -> 31, with exit 0 and no warning. The only net is `04_subset_and_annotate_trios.sh:123`, which dies solely when Step 3 keeps ZERO sites — SpliceAI/ClinVar keeps (7 here) keep it silent. Amplifier: CADD absence is only a WARN (`run_pipeline.sh:127`), so keep_impacts malformed + no CADD leaves only SpliceAI hits and ClinVar P/LP. Expected behaviour: reject a non-list or a value outside {HIGH,MODERATE,LOW,MODIFIER} at load, or at minimum log the resolved set and audit an `impact_rung_matched` count.

**Who it hits.** Every coding variant the screen exists to find — all pLoF and all missense — except those rescued by CADD >= 25.3 or SpliceAI >= 0.2. On the mock that is 29 of 43 union sites (67%). On real WGS most nonsense is rescued by CADD (relabelled `cadd`), so the net loss concentrates on missense, in-frame indels and MODERATE splice-region calls, i.e. the majority of the coding candidate list.

**What the verifiers corrected.**

- Two small refinements to the finding as written, neither of which changes the verdict:

1. "the drop is recorded under reason.not_functional" is right, but it is worth being precise that counts.tsv is not entirely trace-free: `reason.impact_high` and `reason.impact_moderate` DISAPPEAR from the audit rather than reading…
- The mechanism is real and I reproduced it end to end: scalar, title-case, lowercase and empty keep_impacts all silently disable the impact rung (36/43 -> 7/43 plausible on the real integration union, exit 0, no warning), there is no config validation anywhere in the repo, and selection.py:39 is the only config-driven s…
- The mechanism, the executed numbers (36/43 -> 7/43, not_functional 2 -> 31, exit 0), the "not documented" status and the "not fixed on this branch" status all hold exactly - I reproduced them against the current code. Two claims do not hold. (1) The CLASS is wrong: this is not silent_loss as the brief defines it. 03_se…

**Fix.** Validate at load: require a list, upper-case each token, reject unknown IMPACT tokens, and audit the count of variants kept at the impact rung so a dead rung is visible in counts.tsv.

### D29. The impact rung reads only the MANE-Select transcript's IMPACT: a stop_gained on an alternative transcript of the SAME gene is filed as `not_functional`

**`src/hprv/selection.py:42`** · silent loss · FILTER · verifiers: confirmed/partly/confirmed · reproduced by execution · documented at NO — docs/pipeline_design.md:216-221 and CLAUDE.md:167 discuss `--pick_order` only as a GENE-ATTRIBUTION hazard (SYMBOL naming a readthrough gene), and :221 says 'set csq_select: mane if gene attribution matters more than recall', implying the default `pick` is the recall-preferring option. For the impact rung it is the opposite: `worst` is strictly more sensitive.

Step 2 runs VEP with `--flag_pick --pick_order mane_select,mane_plus_clinical,canonical,rank` (`02_annotate_sites.sh:313`) and split-vep selects `-s pick`, so within one gene the MANE transcript wins at criterion 1 and `rank` (worst consequence) is never consulted. `A.impact(v)` therefore returns the MANE block's IMPACT, and the ladder's first rung judges the variant on that alone; a more severe consequence on any alternative transcript of the same gene is discarded before `classify()` is called.

**Failure scenario.** CDKN2A: a variant that is `stop_gained` (HIGH) on ENST00000579755 (p14ARF, MANE Plus Clinical) and `intron_variant` (MODIFIER) on ENST00000304494 (p16INK4a, MANE Select). Executed with real bcftools 1.22 `+split-vep`: `-s pick` yields `vep_IMPACT=MODIFIER`, and the real classifier returns `(False,'not_functional')`; `-s worst` on the identical record yields `vep_IMPACT=HIGH` and `(True,'impact_high')`. A second record, `5_prime_UTR_variant`/MODIFIER on MANE and `missense_variant`/MODERATE on the alternative transcript with CADD 22.1, behaves the same way. The variant is counted as `reason.not_functional`, which asserts the functional evidence was consulted and was insufficient. Expected: either the worst consequence within a gene should feed the impact rung, or the MANE-only recall cost should be stated where the pick hazard is documented.

**Who it hits.** Any variant whose damaging consequence falls on a non-MANE transcript of the same gene: alternative first/last exons, and the ~60 MANE-Plus-Clinical genes where the second isoform is clinically required (CDKN2A p14ARF/p16INK4a dual reading frames is the in-scope pediatric-cancer example). CADD often rescues true LoF here (nonsense CADD is typically >30) but does not rescue missense or inframe changes, and CADD is optional and warn-only.

**What the verifiers corrected.**

- Three corrections/additions; the core claim stands.

1. THE CDKN2A EXAMPLE AS WRITTEN DOES NOT TRIGGER — use the shared-exon case instead. The finding's scenario is `stop_gained` on p14ARF (ENST00000579755) / `intron_variant` on p16INK4a (ENST00000304494). CDKN2A exon 1beta sits ~19 kb outside the p16 transcript (p16 s…
- The impact rung does read only the PICK'd (normally MANE Select) block — confirmed by execution. But the claimed consequence is overstated on four counts.

(a) The residual lost class is much narrower than "any variant whose damaging consequence falls on a non-MANE transcript". It is: MODIFIER/LOW on the MANE Select tr…
- Two mechanism imprecisions, neither refuting the finding. (1) `--flag_pick` marks ONE block per RECORD, not per gene, so the loss also spans overlapping genes, slightly broadening the affected class beyond "the same gene". (2) The claim that `rank` "is never consulted" holds only when a single transcript is MANE Select…

**Fix.** Lift a second `vep_*_worst` IMPACT column (split-vep can be run twice, or `-s worst` into a parallel prefix) and let the impact rung take the max of picked and worst, keeping the picked block for SYMBOL/gene attribution.

### D30. A NOVEL non-coding indel (and any MNV) can reach no keep-path at all, and the drop is filed as `not_functional` — the same string a scored-and-benign variant gets

**`src/hprv/selection.py:78`** · silent loss · FILTER · verifiers: confirmed/partly/partly · reproduced by execution · documented at NO — and three places assert the opposite: docs/functional_annotation.md:139 ("whole_genome_SNVs.tsv.gz + the gnomAD indel table, so SNVs **and** indels are scored genome-wide"), pipeline/02_annotate_sites.sh:318-320, and resources/manifest.env:64 ("the COMPLETE CADD source"). docs/limitations.md:51-57 DOES document the SpliceAI precomputed indel gap and then names CADD as the backstop — which is exactly the claim that fails — and docs/limitations.md:318-325 ("Reading a negative result") presents a "SNV/indel callset" as uniformly screened. The documentation therefore understates this badly.

Below MODERATE impact there are exactly two keep-paths: SpliceAI (selection.py:52-53) and CADD (selection.py:56-57). The SpliceAI precomputed set covers only 1 nt insertions and deletions <= 4 nt — the repo states this itself at src/hprv/prioritize.py:917-919 and docs/limitations.md:51-57 — and the CADD resource for indels is `gnomad.genomes.r4.0.indel.tsv.gz` (resources/manifest.env:68, 1.2 GB), CADD's precomputed scores for the indels OBSERVED in gnomAD r4.0, NOT all possible indels (contrast whole_genome_SNVs.tsv.gz, 81 GB, which genuinely is every possible SNV). A novel indel is by construction absent from gnomAD, so both getters return None and selection.py:78 emits `not_functional`. Neither resource contains MNVs at all. Step 2's VALUE-level 0-lift guards (02_annotate_sites.sh:527-545) count scored SITES cohort-wide, so an indel score arm that is 100% dead is masked by the SNV majority and the run logs a healthy "CADD scores present on N / M sites".

**Failure scenario.** A private 12 bp deep-intronic deletion in a recessive gene: VEP IMPACT=MODIFIER, absent from gnomAD so no CADD score, >4 nt so no precomputed SpliceAI score. Executed through the real `build_classifier`: DROP (not_functional). Identical outcome for a novel 6 bp 5'UTR insertion and for an `AC>GT` synonymous MNV. The control — the same 12 bp deletion rated MODERATE (inframe_deletion) — is KEPT, and a SNV at the same locus with CADD 28.0 is KEPT via the cadd rung. So the loss is representation-specific, not biology-specific, and it lands squarely on the ultra-rare class the screen exists to find; docs/limitations.md:53-57 notes larger indels are ~4x enriched for splice effects (4.7% vs 1.1%). Should instead be retained (or dropped under a distinct, audited reason such as `unscored_representation`) so a reviewer can tell "scored and benign" from "no predictor covers this representation".

**Who it hits.** Every indel not present in gnomAD r4.0 genomes, and every MNV, whose VEP IMPACT is below MODERATE. In WGS that is on the order of a thousand novel non-coding indels per proband; the coding subset is rescued by the impact rung, the non-coding subset (~98% of indels) is not.

**What the verifiers corrected.**

- Four sharpenings, none of which undercut the finding:

1. THE MNV LIMB IS WEAK and should be demoted to a footnote. It reproduces in the classifier (pos 1000, AC>GT, LOW -> DROP not_functional) and Step 1's `norm -m-` does not atomize it. But the input here is GMKF Kids First GATK output (CLAUDE.md: "GATK genotype-refi…
- The drop behaviour is real and reproduces exactly, but the finding is misclassified and its documentation claim is wrong.

(a) NOT `silent_loss`. pipeline/03_select_plausible.py:82-83 audits `reason.{r}` for drops too, so every such variant is counted as `reason.not_functional` in audit/counts.tsv, with a reason string…
- The mechanism is CONFIRMED by execution and is NOT fixed on this branch (selection.py is absent from `git diff main..HEAD --stat`; the review log has no indel finding). Four corrections:

1. HALF OF IT IS ALREADY DOCUMENTED and must be labelled so. The SpliceAI precomputed-coverage gap (all SNVs, but only 1 nt insertio…

**Fix.** Split the Step-3 drop reason into `not_functional` (scored, below cutoff) vs `unscored_representation` (indel/MNV with neither a CADD nor a SpliceAI value), audit it, and stratify Step 2's 0-lift guards by variant class (SNV / indel) so a dead indel arm dies loudly.

### D31. Step 4 rewrites the per-trio intermediates in the persistent tmpdir but index_vcf refuses to re-index them, so a stale index silently truncates the isec

**`pipeline/04_subset_and_annotate_trios.sh:257`** · fragility · UPSTREAM · verifiers: confirmed/confirmed/confirmed · reproduced by execution · **not documented anywhere**

`$norm` and `$cand` live in $HPRV_TMPDIR, which run_pipeline.sh:163 defaults to the PERSISTENT $W/tmp, and are deleted only on the success path (04:269). `index_vcf` returns early whenever any .tbi/.csi exists (lib/common.sh:176). 04:162 removes the stale index of `$out` before rebuilding but not those of `$norm`/`$cand`, so a rewritten `$norm` is read by `bcftools isec` (04:259) through the previous run's index. The identical defence -- `rm -f "$PLAUSIBLE_TX".tbi "$PLAUSIBLE_TX".csi` -- sits 160 lines above at 04:98 under a CRITICAL comment naming exactly this failure ('bcftools annotate then reads offsets that no longer match the rewritten data and SILENTLY drops the transfer ... exit 0, .done still stamped').

**Failure scenario.** Step 4 is killed by walltime/OOM during trio T (after 04:257 stamped the index, before 04:269 cleaned up). Step 3 is then re-run with different thresholds, or Step 2b backfills SpliceAI, so plausible.sites.vcf.gz changes. On resume, `_tkey` correctly invalidates trio T, `$norm` is rebuilt with a different record set from the changed region BED, `index_vcf` no-ops on the surviving .tbi, and `bcftools isec` returns only the records the old index covers. `.done` is stamped with the CORRECT key, `require_intact_bgzip "$out"` passes (the file is valid, just short), and `audit 04_subset candidate_genotypes` reports the truncated count as though it were the truth -- so the loss is invisible in audit/counts.tsv and in summary.md. It should re-index, as line 98 already does for the annotation source.

**Who it hits.** An arbitrary coordinate-contiguous suffix of one trio's candidate set -- in the demonstration an entire contig. Hits only a trio rebuilt after an interrupted Step 4 whose plausible set or source VCF changed, but when it hits it removes whole chromosomes of candidates for that trio.

**What the verifiers corrected.**

- Two refinements to the finding as written, neither of which weakens the core claim:

1. `$cand`'s stale index is NOT load-bearing — only `$norm`'s is. The single consumer of `$cand` is `bcftools annotate -a "$PLAUSIBLE_TX" ... "$cand"` (04:264), which reads its main input as a stream and never opens `$cand.tbi`; only t…
- Three refinements to the finding as written, none of which weaken it materially:

1. The loss is NOT "a coordinate-contiguous suffix". Executed: a same-contig mid-file change kept 457 of 1000 records in a scattered pattern. The whole-contig case the finder demonstrated is the narrow regime; the realistic one (Step 3 re…
- Not documented anywhere (limitations.md / CLAUDE.md / inheritance_and_genotype_qc.md / pipeline_review_2026-09.md all searched) and not fixed on this branch — the only Step-4 change in main..HEAD is the S3 `_tkey` cache key, and 04:257/260 are still bare `index_vcf` calls with the cleanup at 04:269 on the success path…

**Fix.** Add `rm -f "$norm" "$norm".{tbi,csi} "$cand" "$cand".{tbi,csi}` beside the existing `rm -f "$out" ...` at 04:162, mirroring the PLAUSIBLE_TX defence at 04:98.

### D32. resolve_trios re-runs on every invocation and can re-point a trio onto the pipeline's OWN derived output, replacing the source callset with an already-filtered candidate VCF

**`pipeline/resolve_trios.py:126`** · silent loss · UPSTREAM · verifiers: confirmed/partly/confirmed · reproduced by execution · **not documented anywhere**

`enumerate_vcfs` globs `--vcf-dir` RECURSIVELY (`**/*.vcf.gz`, :37-38) and `chosen = sorted(candidates, key=lambda v: (len(v2s[v]), v))[0]` (:126) picks the candidate with the FEWEST SAMPLES, then lexical path. Nothing inspects record count, contigs, or provenance. run_pipeline.sh:191 re-runs resolve whenever `FROM <= 1` or the manifest is absent, so a resume re-resolves. The pipeline writes 3-sample per-trio VCFs into its own work dir (`trios/*.candidates.annotated.vcf.gz`, `igv/vcfs/*.vcf.gz`), which are strictly MORE trio-specific than a source callset carrying extra members — so the tie-break actively prefers them.

**Failure scenario.** `vcf_dir` is set to a delivery/project root and `work_dir` sits beneath it (the natural layout). The first run is correct. Any re-run — a resume after a walltime kill, or `--from 0/1` — re-resolves and, for any trio whose source VCF carries an extra member (a sibling), selects the pipeline's own Step-4/Step-8 output: an already rarity- and impact-filtered handful of records. Step 0 then measures MIE/sex on those few sites, Step 1's union for that trio becomes those sites, and every variant the trio could have contributed is gone. Exit code 0; the only trace is `status=resolved_multi` / `n_candidate_vcfs=3` in trio_resolution.tsv and `trios_multi_vcf` in the audit — neither says the chosen file is derived, and there is no per-trio stderr WARN.

**Who it hits.** All variants of every affected trio — the trio is analysed against a pre-filtered subset, so its entire contribution to the cohort union, to its own candidate list and to Step 6's carrier counts collapses. Triggered by directory layout, not by variant class, so it is all-or-nothing per trio.

**What the verifiers corrected.**

- The finding is correct and, in one respect, UNDERSTATED; it also omits one precondition it should state explicitly.

UNDERSTATED: the finding attributes the preference for derived files to the fewest-samples rule and therefore scopes the trigger to "any trio whose source VCF carries an extra member (a sibling)". My syn…
- Three corrections, none of which kills the finding.

(a) CLASS: not `silent_loss` as defined ("NO counter, NO flag, NO reason string"). `trio_resolution.tsv` writes the chosen path (which names the work dir), `status=resolved_multi` and `n_candidate_vcfs`; `audit/counts.tsv` gets `resolve/trios_multi_vcf`; the stderr s…
- Two qualifications, neither refuting: (a) the finding understates the trigger — the fewest-samples preference is only one route; the LEXICAL tie-break alone selects the derived output for a plain 3-sample source trio VCF whenever the work dir's path sorts before the source dir's (demonstrated: proj/analysis/... chosen…

**Fix.** Exclude the work dir from the glob, break ties on record count / contig coverage rather than sample count, and emit a per-trio stderr WARN listing every candidate whenever `n_candidate_vcfs > 1`.

### The remaining moderate findings

- **A trio VCF with no FORMAT/AD yields ZERO calls in every inheritance mode, and the only trace is a per-trio WARN that scopes the loss to multiallelic 1/2 hets**  
  `pipeline/04_subset_and_annotate_trios.sh:207` · silent loss · verifiers confirmed/partly/partly · documented at NO. The warning text at pipeline/04_subset_and_annotate_trios.sh:207 is itself the doc, an  
  `sample_qc` requires `ab is not None` for `het`, `hom_alt` and `denovo_child` (genotype.py:141,143), so every CARRIER limb fails closed when allele depths are absent, while `hom_ref`/`clean_parent` fail open. Step 4 gates `--keep-sum AD` on the presence of `##FORMAT=<ID=AD,` in the header (04:205-208) and, when it is missing, warns only that "a multiallelic (1/2) het cannot have its AB corrected, so such calls may be dropped" — but the consequence is total, not multiallelic-specific. Step 5 has no `variants_examined` counter and none of its `continue`s is counted, and neither Step 5, Step 6 no…  
  *Fix:* Hard-fail Step 4 (or Step 5) when a trio VCF declares no FORMAT/AD, and add a per-trio `variants_examined` / `candidates_no_mode` audit metric so any total loss is arithmetically visible.
- **The parental-mosaicism de novo the docs instruct you to detect is hard-gated away, and parental allele balance is never written to the output**  
  `pipeline/05_inheritance_screen.py:194` · silent loss · verifiers confirmed/partly/partly · documented at NO — docs/inheritance_and_genotype_qc.md:201 asserts the opposite approach, and docs/limit  
  The de novo path requires `sample_qc(..., "clean_parent")` on both parents (:194/:201), which fails at genotype.py:148-149 on `alt_ad > parent_max_alt_ad` (1, absolute). The het collector then also yields origin=None for that variant (:288-289 requires both parents to pass `clean_parent`), so there is no dominant fall-through either. `COLS` (:31-54) carries child_gq/child_dp/child_ab but no parental GQ/DP/AD/AB, so the inspection the doc prescribes is impossible even for the calls that survive.  
  *Fix:* Emit the call with `flags=parental_alt_reads=N` under never-drop, and add parent GQ/DP/AB columns to COLS.
- **hom_recessive truth table: a 1/1 child with one 0/0 parent (deletion-in-trans / UPD) or one no-call parent produces no row, no flag and no counter, while the docs instruct the reader to 'suspect' exactly those sites**  
  `pipeline/05_inheritance_screen.py:230` · silent loss · verifiers confirmed/partly/partly · documented at docs/ROADMAP.md:48  
  The branch at :230-234 requires BOTH parents in {HET, HOM_ALT}. A hom-alt child fails `gc == G.HET` at :257 so there is no fall-through into the het collector either. The configuration is discarded with no trace; Step 0's only Mendelian-error output is a genome-wide `mie_rate` over the first `qc.max_sites` (200000) sites (00_qc.py:151,246), so a single-locus event contributes 1/200000 and can never raise `mie_flag`.  
  *Fix:* Emit the 1/1-child / 0/0-parent and 1/1-child / no-call-parent configurations with a `mendelian_error_suspect_deletion_or_upd` flag, or audit a per-trio counter for each.
- **A homozygous chrX call in a daughter is dropped whenever the father is not a pristine `1/1` — the exact discipline the affected-male branch was written to avoid**  
  `pipeline/05_inheritance_screen.py:249` · silent loss · verifiers confirmed/partly/partly · documented at docs/inheritance_and_genotype_qc.md:176  
  The affected-female X-linked branch hard-requires `gd == G.HOM_ALT` (:249) and then QCs the father at the hom-alt allele-balance band (>= homalt_ab_min 0.90) at :251, while chrX non-PAR is excluded from the autosomal `hom_recessive` branch at :230 (`not G.is_x_nonpar(v)`) and the het collector is not reached for a HOM_ALT child — so there is no fall-through under any mode. The parallel affected-male branch at :239-241 deliberately does not consult the father at all, its own comment stating that 'an affected/carrier father or a father chrX no-call must not drop the call'. The same protective re…  
  *Fix:* Accept `gd in (HET, HOM_ALT, UNKNOWN)` for the daughter's call, record the paternal state as a flag (`father_gt_uninformative` / reuse `father_carries_x_allele`), and add a per-trio audit counter for the records this branch rejects.
- **`filters.denovo.parent_max_alt_ad` silently gates the RECESSIVE model: at the shipped default of 1 absolute alt read, a compound-het pair with an apparent de-novo second hit yields zero rows under any mode**  
  `pipeline/05_inheritance_screen.py:288` · silent loss · verifiers confirmed/partly/confirmed · documented at docs/pipeline_review_2026-09.md:51 (the absolute count is listed under 'Deliberately NOT c  
  At :288 a het whose parents are both HOM_REF gets `origin="denovo"` only if BOTH parents pass `sample_qc(..., "clean_parent")`, which requires `alt_ad <= thr.parent_max_alt_ad` (genotype.py:149) — an ABSOLUTE count defaulting to 1, read from `filters.denovo.parent_max_alt_ad`. `by["denovo"]` is the only source of de-novo-partner pairs (:331-332). Fail it and the leg is never collected, so there is no compound_het pair, no denovo row, and — for a frequency in the recessive band — no dominant row either, since the dominant fall-through gates at `dominant_max` (:366). No counter records any of it…  
  *Fix:* Make `clean_parent` depth-relative (a fraction plus a small absolute floor), move/duplicate the knob under `filters.genotype_qc` with documentation naming its recessive-model consumer, and audit the count of hets rejected for parental alt reads.
- **A het child with one affirmatively clean parent and one no-call parent is emitted under no mode at all — 2 of the 64 truth-table cells, and the autosomal case is undocumented**  
  `pipeline/05_inheritance_screen.py:313` · silent loss · verifiers confirmed/partly/partly · documented at docs/limitations.md:297-304 documents the analogous hole for a MALE non-PAR chrX hemizygot  
  At :312-313 the final `else: origin = None` catches every configuration where neither parent is a called carrier and the pair is not (HOM_REF, HOM_REF). A het child with father HOM_REF and mother UNKNOWN falls there: the de novo branch at :184 requires BOTH parents HOM_REF so it does not fire, and `if origin:` at :314 keeps the variant out of `hets`, so there is no dominant row and no comp-het leg. Genetically the allele is either maternally transmitted or de novo — both are reportable candidates — and the pipeline instead reports neither, with no flag and no counter.  
  *Fix:* In the `else` branch, when the OTHER parent is an affirmative clean 0/0, emit with `origin=<other>_or_denovo` and `flags=parent_gt_uninformative` rather than returning None; failing that, add a per-trio `origin_none_parent_nocall` audit metric.
- **The AB bands are pure fractions with no read-count floor and no depth-dependent widening, so at the configured `min_dp: 10` the het band alone rejects 10.9% of true hets per sample (20.7% of dominant calls, which need it twice) — silently**  
  `src/hprv/genotype.py:140` · silent loss · verifiers confirmed/partly/partly · **NEW**  
  `sample_qc` compares a point estimate `alt/(ref+alt)` against fixed constants (0.25-0.75 het, >=0.90 hom-alt, <=0.10 hom-ref) with no binomial tolerance and no minimum informative-read count. The band width is therefore constant while the sampling noise scales as 1/sqrt(DP), so the rejection rate is entirely determined by depth — and it is worst exactly at the `min_dp` floor, i.e. in the low-coverage GC-rich first exons and segdup-adjacent regions where disease genes are hardest to call. Step 5 requires the band to pass independently in the child AND in the transmitting parent (05_inheritance_…  
  *Fix:* Replace the fixed AB fractions with a binomial-tail test at the configured depth, or add a per-trio `ab_band_rejected` counter so the loss is at least visible.
- **A trio VCF with no FORMAT/AD yields ZERO inherited calls of any mode — every QC limb fails closed on allele balance, and the only signal is a Step 4 WARN that describes a far narrower consequence**  
  `src/hprv/genotype.py:141` · fragility · verifiers confirmed/confirmed/partly · documented at PARTIAL and it UNDERSTATES. CLAUDE.md:486-495 documents the per-sample AD fail-open/fail-c  
  `sample_qc` requires `ab is not None` for the `het`, `hom_alt` and `denovo_child` kinds (:141, :143). With no AD, `allele_balance` returns None for every sample, so the child fails het QC at :257, hom-alt QC at :231/:249, and denovo_child QC at :202 — the het collector, the hom-recessive branch, both X-linked branches and the de novo branch all return nothing. Step 5 then records `candidate_calls 0` and no other diagnostic. Step 4's guard at 04_subset_and_annotate_trios.sh:207 warns only that "a multiallelic (1/2) het cannot have its AB corrected, so such calls may be dropped", which describes…  
  *Fix:* Assert `##FORMAT=<ID=AD` in the Step-5 preflight (dying like the spliceai_required / faf95-slim gates), or degrade explicitly: skip only the AB limb and emit with `flags=ab_unmeasured`, recording a per-trio `ab_unmeasured` count.
- **`--keep-sum AD` fabricates reference reads under a HOM-ALT genotype, so a true homozygous hit at a multiallelic site fails the hom-alt AB band and is dropped with no counter**  
  `src/hprv/genotype.py:142` · fragility · verifiers confirmed/confirmed/confirmed · documented at NO - docs/inheritance_and_genotype_qc.md:90 documents this exact interaction for the HET d  
  Step 4 passes `norm -m- --keep-sum AD` (pipeline/04_subset_and_annotate_trios.sh:206) so that the other ALT's reads are folded back into AD[0] -- the correct per-allele semantic for a 1/2 HET, which is the case it was added for. Under a HOM-ALT genotype the identical fold turns third-allele reads into apparent REFERENCE reads, so allele_balance() = alt/(alt+other_alt) drops below `homalt_ab_min` (0.90) as soon as more than 10% of that member's reads support any other ALT. sample_qc(..., 'hom_alt') (genotype.py:142-143) then returns False and every mode gated on it emits nothing -- no row, no f…  
  *Fix:* Evaluate the hom-alt AB band against total site depth (or exclude other-ALT reads from the denominator) and emit a `multiallelic_ab_corrected` flag instead of dropping the call.
- **docs claim CADD scores indels 'genome-wide'; it does not, and a private indel below MODERATE impact therefore has neither rescue rung by construction**  
  `docs/functional_annotation.md:139` · doc mismatch · verifiers confirmed/partly/confirmed · documented at docs/limitations.md:44-52 documents the SpliceAI half accurately (1 nt insertions, deletio  
  The doc states CADD v1.7 is used 'via the dedicated plugin: whole_genome_SNVs.tsv.gz + the gnomAD indel table, so SNVs **and** indels are scored genome-wide'. Only the SNV file is exhaustive; the indel table is a gnomAD-derived lookup, as `pipeline/02_annotate_sites.sh:318-319` correctly calls it ('+ the precomputed indel set'). SpliceAI's indel file has the same property (limitations.md:44-52). Both rescue rungs at `selection.py:50,56` are therefore lookups keyed on previously-observed variants, so their coverage is anti-correlated with the rarity the pipeline selects for.  
  *Fix:* Change functional_annotation.md:139 to say CADD scores all SNVs genome-wide plus a precomputed (non-exhaustive) indel set, and note that a private indel below MODERATE impact has no keep-path.
- **Step 0 opens the trio VCFs without `strict_gt`, so a half-called `1/.` counts as a chrX HET in sex inference — a male mis-inferred as female loses every chrX hemizygous call with no counter and `sex_match=1`**  
  `pipeline/00_qc.py:107` · wrong call · verifiers confirmed/confirmed/confirmed · documented at docs/inheritance_and_genotype_qc.md:116 — PARTIAL, and it UNDERSTATES: it excuses Step 0's  
  `VCF(vcf_path)` at :69 and :107 uses cyvcf2's default `strict_gt=False`. Executed: under that default `0/.` reads HOM_REF and `1/.`/`./1` read **HET**. `scan_sex` (:88-102) increments `x_het` on any HET, so half-calls push `x_het_ratio` above `qc.x_het_male_max` (0.10) and `inferred_sex` flips from '1' to '2' (:156). Step 5 consumes that value and routes on it at :156/:169/:180. The MIE guard at :137-138 (`if G.UNKNOWN in (gc,gd,gm): continue`) is likewise blind to half-calls, so it cannot skip what it was written to skip.  
  *Fix:* Open with `strict_gt=True` in `scan_sex` and `qc_trio` (00_qc.py:69,107), and emit `sex_match` as empty rather than '1' when `inferred_sex` is None.
- **Neither Step 1's nor Step 4's cache key includes the reference FASTA, and Step 4's omits HPRV_PLAUSIBLE_PAD — so a changed reference silently reuses site files normalized against the old one (bypassing the `-c e` "halt loudly on a build mismatch" guard), and widening the pad to recover a missed record is a no-op**  
  `pipeline/01_make_cohort_sites.sh:77` · silent loss · verifiers confirmed/partly/confirmed · documented at NO (docs/pipeline_review_2026-09.md S3 fixed the source-VCF half of these keys; the refere  
  `trio_content_key` (01:76-79) hashes vcf path + size/mtime + FILTER + EXCLUDE_CONTIGS + samples, and the union key (01:90) is the manifest plus those — the reference FASTA appears in neither. Step 4's key (04:152-153) is plausible-set cksum + vcf + size/mtime + samples — no reference and no pad. But the reference is exactly what `bcftools norm -f` left-aligns against and what `-c e` checks; 01:154-158 states the `-c e` choice exists so "a build/contig mismatch must halt the run loudly". On any resume that halt cannot fire, and Step 4's `-c w` (04:203) only writes REF_MISMATCH lines to stderr.  
  *Fix:* Add the reference path + `hprv_stat_key` (and PLAUSIBLE_PAD for Step 4) to both content keys, and audit Step 4's `-c w` mismatch count.
- **Step 2's rarity-oracle guards still police the retired oracle: the value-level guard is satisfiable by `vep_MAX_AF` alone, and the header guard hard-dies on the proxy fields even under the faf95 arm that never reads them**  
  `pipeline/02_annotate_sites.sh:684` · fragility · verifiers confirmed/partly/confirmed · documented at NO — commit 55c892f ('collapse ~40 passages describing an oracle design the code retired')  
  The value guard at :684-693 ORs `vep_gnomADe_AF`, `vep_gnomADg_AF` and `vep_MAX_AF` into `freq_expr` alongside the ten grpmax-eligible fields, so a single populated MAX_AF satisfies a check whose stated purpose is to prove the RARITY ORACLE has values. MAX_AF and the global AFs are REPORTING-ONLY under both arms (golden rule 2), so the guard does not test the oracle at all. Symmetrically, the header guard at :442-443 (`_n_grpmax > 0`) and its comment 'The grpmax-eligible AF fields ARE the rarity oracle now' predate the faf95 default and are applied unconditionally, with no reference to `_oracl…  
  *Fix:* Build `freq_expr` from `$GRPMAX_AF_FIELDS` only (never MAX_AF/global AF), and gate both the header die and the value check on `_oracle == grpmax_proxy`, moving the `_oracle` assignment above :442.
- **A headerless trios file silently consumes its first trio as a header row**  
  `src/hprv/ped.py:36` · silent loss · verifiers confirmed/confirmed/confirmed · documented at NO — and CLAUDE.md documents the OPPOSITE precedent for the `--gene-prior` reader (headerl  
  `read_trios_file` unconditionally consumes the first line as the header (:36). When the file has no header, that line's kid/dad/mom IDs become column NAMES, none match an alias, the per-column positional fallback fires (:46-51), and parsing proceeds normally on the remaining rows. The first trio is dropped with no warning, no exception and no counter — `audit.record("resolve", "trios_input", ...)` reports the post-loss count.  
  *Fix:* Hard-stop (or WARN loudly) when NO header alias is recognised AND the first line's fields look like sample IDs; at minimum log which columns were used and the raw line count vs parsed trio count.

## Findings — minor

Real, reproduced, and either rare, retention-biased (the pipeline keeps too much rather than too little), documentation-only, or contradicted in part by a verifier. Listed so they are not re-found.

- **The hiConfDeNovo gate re-imposes on the male-X hemizygous de novo path the exact paternal requirement that path exists to remove** — `pipeline/05_inheritance_screen.py:210` · wrong call · verifiers confirmed/refuted/partly · finder said major  
  *Verifier:* The mechanism, the file:line, the reproduction and the "not already documented" status are all correct as stated. Three framing corrections:

(1) IMPACT IS NARROWER THAN CLAIMED. The finding argues from "a hemizygous LoF in an affected boy is causally self-suf…
- **A PARTIALLY-covered gnomAD slim passes Step 2's 0-match die, Step 3's header guard and `prepare_resources verify`, and then reports BA1-common alleles** — `pipeline/02_annotate_sites.sh:643` · unverifiable · verifiers confirmed/partly/partly · finder said major  
  *Verifier:* Three refinements, none of which lowers the severity.

1) DIRECTION. The finding is filed under the audit's false-negative hunt, but the failure is purely retention-biased: absent reads as rarest, so nothing is lost — common alleles are ADDED and every emitted…
- **Step 1 filters each trio to FILTER PASS/. but Step 4 re-derives genotypes with no FILTER at all, so a record filtered in trio A but PASS in trio B re-** — `pipeline/04_subset_and_annotate_trios.sh:232` · silent loss · verifiers partly/refuted/partly · finder said major  
  *Verifier:* Two corrections to the finding as written.

1. The gene-layer sentence "GENEA has 1 carrier instead of 2" implies the correct answer is 2. It is not. Under the pipeline's own policy the record IS non-PASS in TA, so 1 carrier is the correct count; applying `-f …
- **`hprv_keep_reason` is carried into every per-trio candidate record and then read by nothing — the free per-record witness of the handoff, and Step 3's** — `pipeline/05_inheritance_screen.py:31` · unverifiable · verifiers partly/refuted/partly · finder said moderate  
  *Verifier:* Three corrections to the finding as written.

1. The finding's worked example is wrong for the stock integration run. chr2:8000 GENE7 is NOT emitted as a `dominant` call — it appears in candidates.calls.tsv under NO mode at all. Its faf95 is 0.00016, above dom…
- **`bcftools norm -m-` copies the site-level hiConfDeNovo tag to every split ALT, so the reported `hiConfDeNovo` column reads 1 on inherited calls** — `pipeline/05_inheritance_screen.py:101` · wrong call · verifiers partly/refuted/partly · finder said moderate  
  *Verifier:* The per-allele stamping at pipeline/05_inheritance_screen.py:101 is genuinely allele-blind and I reproduced a `mode=dominant flags=origin=pat` row carrying hiConfDeNovo='1' using the real Step 5 CLI. But the input required to trigger it — a multiallelic record…
- **Every Python boolean knob is `bool(value)`, so a quoted or `${ENV}`-templated `false` reads TRUE — the gates a user tried to relax stay on, unwarned** — `pipeline/05_inheritance_screen.py:123` · silent loss · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* Two corrections, one of which makes the finding WORSE than stated.

1. UNDERSTATED — the failure runs in BOTH directions, and the other direction is a much larger silent loss. The finding only reports "string is truthy so the knob stays ON". But os.path.expand…
- **`high_conf_rarity` never fires on the rarest class: a variant ABSENT from gnomAD is tiered below one gnomAD carries** — `pipeline/05_inheritance_screen.py:138` · wrong call · verifiers confirmed/partly/confirmed · finder said moderate  
  *Verifier:* Two corrections, neither touching the core claim.

(a) THE AFFECTED COUNT IS OVERSTATED: 4 of 41, not "7 of 41 calls (17%)". tag_strict is applied at only four call sites -- :234 (hom_recessive), :244 (x_linked male), :252 (x_linked female) and :356 (compound_…
- **The male non-PAR hemizygous-het red flag is enforced on the proband only, never on the father — so a het-called father on chrX collapses a daughter's ** — `pipeline/05_inheritance_screen.py:169` · wrong call · verifiers confirmed/refuted/confirmed · finder said moderate  
  *Verifier:* Two small refinements, neither of which changes the verdict. (1) The affected-female X-recessive half (:248-252) is not actually caused by the `male_x` scoping at :169 — that branch never consults `male_x` for the father; it independently demands `gd == G.HOM_…
- **`denovo_min_dp` is applied unadjusted to the HAPLOID male chrX, so a hemizygous de novo at typical X depth is dropped with no fall-through and no coun** — `pipeline/05_inheritance_screen.py:204` · silent loss · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* Two factual corrections to the submitted text, neither affecting the verdict:

1. Wrong config key path. The threshold is `filters.genotype_qc.denovo_min_dp`, not `filters.denovo.denovo_min_dp`. src/hprv/genotype.py:39 reads it as `g + "denovo_min_dp"` with g …
- **The male-X hemizygous de novo carries a ploidy-blind DP >= 20 floor on a chromosome with half the autosomal depth, so an identical hemizygous LoF is k** — `pipeline/05_inheritance_screen.py:204` · silent loss · verifiers confirmed/refuted/partly · finder said moderate  
  *Verifier:* Two refinements, neither of which changes the verdict. (1) The mechanism is more precisely "two modes apply different child-depth floors to identical child evidence" rather than "chrX sensitivity depends on the mother's carrier status" — the mother's genotype …
- **`use_hiconf_tag` gates the de novo MODE row but not the `origin="denovo"` label that admits a leg to trans pairing — the pipeline refuses to report th** — `pipeline/05_inheritance_screen.py:210` · unverifiable · verifiers confirmed/refuted/partly · finder said moderate  
  *Verifier:* Two parts of the submitted finding are wrong and should be corrected before it is filed.

1. THE STATED CAUSE IS NOT REALIZABLE. The scenario says leg B is untagged because "GATK judged it low-confidence — typically because of parental alt reads". Per this rep…
- **A daughter's homozygous chrX call is deleted without trace when the father's chrX is a no-call or is rendered het — the opposite of the policy the cod** — `pipeline/05_inheritance_screen.py:248` · silent loss · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* Two supporting claims need adjusting; the finding itself stands.

1. "The identical genotype configuration on chr2 emits hom_recessive" is true ONLY for the father-0/1 shape. Executed: with father ./. the chr2 record ALSO produces no row, because the autosomal…
- **chrX affected-female homozygote is lost unless the father is affirmatively HOM_ALT, and line 230 explicitly denies it the autosomal hom-recessive fall** — `pipeline/05_inheritance_screen.py:248` · silent loss · verifiers partly/partly/partly · finder said moderate  
  *Verifier:* Two sub-claims are wrong or imprecise; the core behaviour and its silence are correct.

1. MECHANISM ATTRIBUTION IS PARTLY WRONG. The claim that ":230 explicitly denies it the autosomal hom-recessive fall-through" is a true code fact, but it is the OPERATIVE c…
- **The female X-linked-recessive model hard-requires the father's chrX genotype, so the same record is emitted for a son and silently deleted for a daugh** — `pipeline/05_inheritance_screen.py:249` · silent loss · verifiers confirmed/partly/confirmed · finder said moderate  
  *Verifier:* Two details are wrong or understated, neither affecting the verdict.

1) Doc attribution error. The report calls father `0/1` "the diploid-caller rendering of a hemizygous male the docs themselves describe". The docs describe the opposite: docs/inheritance_and…
- **The `parent_ad_unmeasured` witness is attached to the de novo ROW but never to the de novo ORIGIN that a compound-het pair rests on** — `pipeline/05_inheritance_screen.py:288` · unverifiable · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* One sentence of the write-up is slightly overstated for its own fixture: in that fixture the solo `denovo` row for chr1:2500 IS emitted carrying `parent_ad_unmeasured`, so a reviewer joining on (trio_id, chrom, pos, ref, alt) could in principle recover the wit…
- **Compound-het pair enumeration is quadratic (k x m) with no cap or de-duplication, and is the mechanism that makes candidate_calls exceed candidate_gen** — `pipeline/05_inheritance_screen.py:331` · fragility · verifiers partly/refuted/partly · finder said moderate  
  *Verifier:* The mechanism is exactly as described and fully reproduces, but the finding is mis-labelled as undocumented and its downstream claim is stale.

1. NOT undocumented. The uncapped k x m cross product is finding R3 (severity **major**) in docs/pipeline_review_202…
- **A PED sex code outside {0,'',None,1,2} silently discards Step 0's inference and blanks every non-PAR chrX/chrY call, with no warning at all** — `pipeline/05_inheritance_screen.py:415` · silent loss · verifiers confirmed/refuted/confirmed · finder said moderate  
  *Verifier:* The technical claim is correct in every particular — mechanism, all four line numbers, both executed observations, and the "NO documentation" status all reproduce exactly. Two corrections to its framing:

(1) SEVERITY: moderate -> minor, because the default fl…
- **A declared PED sex silently overrides a contradicting Step-0 chrX inference, and no output column records which sex the ploidy model used** — `pipeline/05_inheritance_screen.py:416` · wrong call · verifiers partly/refuted/partly · finder said moderate  
  *Verifier:* The mechanism, the code lines and every executed consequence are correct. What is overstated is reachability and novelty. (1) The pipeline NEVER emits a PED with a declared kid sex: resolve_trios.py:132 hardcodes kid_sex="0" and run_pipeline.sh:180 regenerates…
- **Step 0's cache key ignores VCF content, so a re-delivered trio VCF re-runs Steps 1 and 4 but serves a STALE inferred sex to Step 5** — `pipeline/run_pipeline.sh:212` · fragility · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* Two small factual adjustments, neither affecting the verdict. (a) The row counts in the finding's evidence block ("3 rows -> 2 rows") come from a reduced fixture; on the full integration candidates VCF the real result is 6 CH_B rows -> 5 (41 -> 40 calls total,…
- **`GtThresholds.from_config` validates nothing: an inverted het AB band silently deletes every dominant and compound-het call, and `homref_ab_max: 0.9` ** — `src/hprv/genotype.py:32` · fragility · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* The finding is correct but UNDERSTATES its blast radius in two ways. (a) An inverted het band zeroes ALL FIVE modes, not just dominant and compound_het: hom_recessive and x_linked_recessive both route parent carrier status through carrier_ok -> sample_qc(..., …
- **`dp()` returns sum(AD) when AD is present and FORMAT/DP when it is absent — one threshold, two different depth quantities, and the AD-less (fail-open)** — `src/hprv/genotype.py:74` · doc mismatch · verifiers partly/partly/confirmed · finder said moderate  
  *Verifier:* The mechanism is real and reproduces, but the finding's stated output consequence is wrong. `child_dp` in `candidates.calls.tsv` is NOT a mixture of two quantities — it is always sum(AD). Every child QC limb (`het`/`hom_alt`/`denovo_child`, 05_inheritance_scre…
- **There is no `denovo_min_gq`: the de novo arm is filtered at the same GQ >= 20 as inherited calls, on a posterior that CalculateGenotypePosteriors move** — `src/hprv/genotype.py:132` · doc mismatch · verifiers partly/refuted/partly · finder said moderate  
  *Verifier:* Three corrections to the finding as written. (1) The "~80 phred" span is the WORST case, specific to the obligate-transmission configuration the finding chose (a 1/1 parent, where the child's hom-ref becomes the Mendelian violation and gains +60). Against the …
- **`parent_max_alt_ad = 1` is an absolute read count with no depth-relative alternative: two alt reads at DP 100 (2%) removes the variant from every path** — `src/hprv/genotype.py:146` · silent loss · verifiers partly/partly/partly · finder said moderate  
  *Verifier:* Three corrections to the submitted finding.

1. "removes the variant from every path, not just the de novo one ... also barred from the compound-het and dominant paths" — the DOMINANT half is false. `pipeline/05_inheritance_screen.py:365` emits dominant only f…
- **The faf95 oracle transfers frequencies from gnomAD records that FAILED gnomAD's own QC filters (AS_VQSR/AC0), and the gnomAD FILTER status is carried ** — `scripts/prepare_resources.sh:217` · silent loss · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* The mechanism, the silence and the failure scenario all reproduce exactly as described; two details are wrong. (1) FILTER STRINGS: the gnomAD v4.1 JOINT release that prep_gnomad slims does not use AS_VQSR/AC0 in its FILTER column — the declared values are PASS…
- **Under `oracle: faf95` an absent faf95 resolves to 0.0 (rarest) regardless of how large the transferred `gnomad_AF_joint` is, so BA1 and every rarity g** — `src/hprv/annotations.py:339` · fragility · verifiers partly/partly/partly · finder said moderate  
  *Verifier:* The mechanism and the end-to-end consequence reproduce against real gnomAD v4.1 data, but three specifics in the submission are wrong.

1. The stated trigger (AN_joint=120, AC_joint=6) cannot occur: the Clopper-Pearson 95% lower bound there is ~0.022 > 0, so g…
- **`clnsig_is_plp` is vetoed by a co-asserted `likely_benign` but not by `benign`, so the ClinVar P/LP override silently fails on real `&`-joined CLIN_SI** — `src/hprv/annotations.py:434` · silent loss · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* Two refinements, neither of which weakens the finding.

(1) The third clause of the predicate is DEAD CODE, which sharpens the diagnosis. `"benign/likely" not in s` can never fire independently: `"Benign/Likely_benign".lower()` = `"benign/likely_benign"`, whic…
- **The fail-CLOSED half of the AD asymmetry deletes a carrier parent and its call with no witness; sample_qc_ad_measured only covers the fail-OPEN half** — `src/hprv/genotype.py:140` · silent loss · verifiers confirmed/refuted/partly · finder said moderate  
  *Verifier:* Three refinements to the finding as written.

1) The finding scopes the loss to "the carrying parent". Executed evidence shows it is broader: the child's limbs fail closed identically, so a record with no FORMAT/AD anywhere loses the site under EVERY mode at o…
- **No coherence check on the four `filters.rarity.*` cutoffs: `dominant_max > recessive_max` is silently inert, and `benign_ba1 < recessive_max` silently** — `src/hprv/selection.py:70` · silent loss · verifiers confirmed/partly/confirmed · finder said moderate  
  *Verifier:* One detail in the finding is imprecise, and the reality is worse than stated. The claim 'nothing logs the effective cutoffs' is wrong: src/hprv/report.py:230-232 DOES echo all three configured values into the Step 7 About sheet — but as independently-effective…
- **`not_functional` is one bucket for 'the scores said no' and 'no score existed' — the distinction the rest of the pipeline treats as load-bearing is er** — `src/hprv/selection.py:78` · unverifiable · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* Three refinements, none of which overturn the finding.

1. The failure scenario's part 1 — "cannot tell whether the non-coding screen ran or merely appeared to" — is OVERSTATED on a DEFAULT run and EXACTLY RIGHT on the ingest path. `02_annotate_sites.sh:519-54…
- **Step 0's chrX sex inference — the sole determinant of X ploidy in the shipped flow — ignores FILTER and is boundary-inclusive on the female side** — `pipeline/00_qc.py:90` · fragility · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* Two corrections, neither of which touches the mechanism.

(a) The documentation claim is wrong on the boundary. docs/README.md:218 states verbatim "chrX-inferred sex vs. PED (het-ratio **< 0.10 → male**, `qc.x_het_male_max`; needs **≥ 20** informative chrX cal…
- **Step 0 computes the Mendelian-error rate WITHOUT `strict_gt`, so half-called parents Step 5 refuses to interpret are charged as Mendelian violations** — `pipeline/00_qc.py:107` · fragility · verifiers partly/partly/partly · finder said moderate  
  *Verifier:* Three corrections. (a) NOT undocumented: docs/inheritance_and_genotype_qc.md:116-117 documents the divergence and calls it deliberate — "(Step 0's QC pass does not set it: for Mendelian-error counting, treating `0/.` as hom-ref is the conservative direction.)"…
- **Step 0 writes contam_flag=0 and contam_source=charr when CHARR could not be computed at all, and a real 0.15 contamination flag disappears** — `pipeline/00_qc.py:203` · unverifiable · verifiers confirmed/refuted/confirmed · finder said moderate  
  *Verifier:* Three refinements, none of which weaken the finding.

1. The trigger is BROADER than "FORMAT/AD is absent". The CHARR accumulator (pipeline/00_qc.py:135-141) only sums hom-alt SNV sites passing `gq >= 20 and dp >= 10`, so `cdp` also stays 0 when AD is present …
- **Step 0's sex-swap gate is unreachable in the shipped pipeline: `sex_match` is a constant 1, so `overall_pass` can never fail on sex** — `pipeline/00_qc.py:221` · doc mismatch · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* Two corrections, both narrowing/redirecting rather than refuting.

(1) IMPACT IS MORE BOUNDED THAN THE FINDING IMPLIES. grep -rn 'sex_match|overall_pass' over src/ and pipeline/ outside 00_qc.py returns NOTHING — overall_pass has no consumer anywhere, and both…
- **`overall_pass=1` is emitted for a trio on which NOTHING was measured — all three Step-0 gates fail open, and Step 0 writes no audit metric at all** — `pipeline/00_qc.py:241` · unverifiable · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* Two surfaces are described slightly wrong. (a) `igv/sample_qc.tsv` carries NO `overall_pass` column — `src/hprv/igv.py:252-253` fixes the fieldnames as `trio_id, role, sample_id, mie_rate, inferred_sex, contam, contam_flag`. I generated it from the failing rep…
- **`filters.genotype_qc.require_pass: false` cannot re-admit a non-PASS record — Step 1's `-f 'PASS,.'` is a hardcoded script default that no config key ** — `pipeline/01_make_cohort_sites.sh:30` · doc mismatch · verifiers partly/partly/partly · finder said moderate  
  *Verifier:* The finding's mechanism, failure scenario and every line reference are correct and I reproduced all of them by execution. Two claims are wrong and change the filing.

1. CLASS is not doc_mismatch, and "documented: NO" is false. Every doc line that mentions req…
- **Step 1's only per-trio metric is named `input_sites` but counts the OUTPUT: the FILTER drop, the chrM exclusion and the multiallelic-split gain are ne** — `pipeline/01_make_cohort_sites.sh:148` · silent loss · verifiers confirmed/partly/confirmed · finder said moderate  
  *Verifier:* Two refinements, neither weakening the finding. (a) The finding's own fixture (a 2-record TA, one non-PASS, `input_sites` 1) was not available to me; I reproduced the identical mechanism at larger scale (5 raw -> `input_sites` 3, four netted events) and indepe…
- **Step 1 strips INFO before the union `norm`, so a symbolic <DEL> loses its INFO/END and bcftools left-aligns it to position 1 of the contig; Step 4 str** — `pipeline/01_make_cohort_sites.sh:165` · silent loss · verifiers confirmed/confirmed/confirmed · finder said moderate  
  *Verifier:* The mechanism needs no correction — every step reproduced as written. Two additions:

AMPLIFICATION (stronger than the finding states, executed): the loss is not merely a relocation, it is an N->1 COLLAPSE per contig. Two trios carrying DIFFERENT symbolic dele…
- **`bcftools view --min-ac 1` trusts a stale INFO/AC instead of recomputing from genotypes, and it is applied on branches that never ran `-s` to refresh ** — `pipeline/01_make_cohort_sites.sh:169` · silent loss · verifiers confirmed/partly/partly · finder said moderate  
  *Verifier:* Three corrections, none fatal to the finding.

1. THE DEFAULT ORCHESTRATED PATH IS IMMUNE, and for a stronger reason than the finding gives. pipeline/resolve_trios.py:132 is the only line that emits a resolved row and it hardcodes the samples field: `mf.write(…
- **Step 2's gnomAD and ClinVar transfer coverage is a stderr log line only, never an audit metric, and the faf95 count has no guard at all** — `pipeline/02_annotate_sites.sh:642` · unverifiable · verifiers confirmed/refuted/partly · finder said moderate  
  *Verifier:* Two refinements, neither of which undermines the finding.

(a) "no distinguishable signature in audit/counts.tsv" is too strong. Executed Step 3 on the broken input writes: sites_in 2 / sites_plausible 2 / reason.cadd 2 — the `reason.ba1` and `reason.too_commo…
- **Step 4's region-BED restriction is not the pure optimization its comment claims: an indel whose raw representation sits more than HPRV_PLAUSIBLE_PAD t** — `pipeline/04_subset_and_annotate_trios.sh:107` · silent loss · verifiers confirmed/partly/confirmed · finder said moderate  
  *Verifier:* The finding is correct in mechanism, in failure scenario, and in documentation status. Three refinements:

(a) EXACT THRESHOLD. Loss requires the left-alignment shift to exceed pad + len(REF_plausible) - 1, not merely pad. At defaults that is 1001 bp: a 1002 b…
- **`prepare_resources.sh` keeps a gnomAD slim that failed its own integrity check and adopts it as `cached` on the next run; `verify` reports it `ok` on ** — `scripts/prepare_resources.sh:182` · fragility · verifiers partly/refuted/partly · finder said moderate  
  *Verifier:* The mechanism is confirmed exactly; the failure SCENARIO is wrong at its last link. A `bcftools concat` interrupted by a walltime kill leaves NO index (index_vcf at :227 has not run), and `bcftools annotate -a` on an unindexed or truncated-unindexed file exits…
- **`tag_strict` inverts on the rarest class: a variant absent from gnomAD never earns `high_conf_rarity`, while a `zero_ci` variant at frequency 0 does** — `pipeline/05_inheritance_screen.py:138` · doc mismatch · verifiers confirmed/partly/confirmed  
  *Verifier:* The finding is accurate as written. Two refinements worth recording: (1) the affected modes are precisely those where `tag_strict` is applied — hom_recessive (line 234), the line-244 recessive row, x_linked_recessive (252) and compound_het (356) — I confirmed …
- **FORMAT/FT is never read anywhere in the pipeline: a genotype the caller itself marked as failing is treated as a confident call, including when a pare** — `pipeline/05_inheritance_screen.py:152` · unverifiable · verifiers confirmed/refuted/confirmed  
  *Verifier:* Two refinements to the finding as written, neither of which weakens it:

(1) The finding's fixture check only demonstrated survival through Step 1. I additionally confirmed the whole path: FT also survives Step 4 (which strips only INFO) and is readable at Ste…
- **Step 5 asserts no biallelic precondition; an un-normalised multiallelic record silently yields zero calls because ref_AD is 0 and the AB reads 1.0** — `pipeline/05_inheritance_screen.py:162` · fragility · verifiers confirmed/partly/partly  
  *Verifier:* Two corrections, neither of which refutes the finding.

1. "CLAIMED DOCUMENTATION STATUS: NO" is not quite right. The precondition IS documented, twice, in
   prose — docs/inheritance_and_genotype_qc.md:90-96 ("Multiallelic caveat — why Step 4 passes
   `--kee…
- **`inheritance.emit_denovo: false` does not omit de novo rows: the de novo variant is still emitted, as a compound_het leg** — `pipeline/05_inheritance_screen.py:184` · doc mismatch · verifiers confirmed/refuted/partly  
  *Verifier:* Two corrections to the finding's narration; the mechanism, line numbers and executed counts are all exactly right.

1. THE "INVISIBILITY" CLAIM IS OVERSTATED AND IS THE PART THAT ARGUES FOR HARM. The finding says the de novo leg is emitted "with no row anywher…
- **A phase-confirmed pair consumes both legs, so an ultra-rare dominant-grade variant is relabelled compound_het by a partner up to 100x too common to be** — `pipeline/05_inheritance_screen.py:341` · unverifiable · verifiers confirmed/partly/confirmed  
  *Verifier:* The finding's claim that "from candidates.calls.tsv alone a reviewer cannot tell which genes lost dominant carriers" is OVERSTATED. The emitted row carries every ingredient needed to reconstruct the absorption: rarity_af=1e-06 (below dominant_max), pair_id=T1:…
- **The methods reference promises `origin_unverified` for an unqualified hom-ref parent, but the code emits two flags that appear in no doc a scientist r** — `pipeline/05_inheritance_screen.py:353` · doc mismatch · verifiers confirmed/partly/partly  
  *Verifier:* Two refinements to the finding's wording, neither of which changes the verdict:

1. "appear in no doc a scientist reads" is right for the four docs named (grep 0 in inheritance_and_genotype_qc.md, limitations.md, README.md, pipeline_design.md), but `trans_evid…
- **`loConfDeNovo` is transferred and declared but read by nothing, while two docs mark the lower-sensitivity tier IMPLEMENTED** — `src/hprv/annotations.py:159` · doc mismatch · verifiers confirmed/partly/confirmed  
  *Verifier:* Two refinements, neither changing the verdict.

(a) The mismatch is in THREE places, not two. Beyond docs/gene_burden.md:137 and docs/inheritance_and_genotype_qc.md:186, two CODE COMMENTS assert a consumer that does not exist: pipeline/04_subset_and_annotate_t…
- **chrX/chrY ALT and unlocalized contigs are not recognised as sex chromosomes, so a male's het calls there are emitted as autosomal `dominant` — includi** — `src/hprv/genotype.py:113` · fragility · verifiers confirmed/partly/confirmed  
  *Verifier:* The finding is correct as written; three corrections/extensions.

(1) UNDERSTATED CLASS. The report files this as `fragility` only. It is also `silent_loss` and `doc_mismatch`. Executed: a textbook X-linked-recessive configuration on chrX_KI270880v1_alt (child…
- **`rarity_oracle()` silently coerces any unrecognised `resources.gnomad.oracle` value to `faf95`, so a config typo selects the opposite arm with no vali** — `src/hprv/annotations.py:295` · fragility · verifiers confirmed/refuted/partly  
  *Verifier:* Two corrections/additions to the finding as written.

(a) The MECHANISM's claim that "shell and Python agree — consistently on the wrong arm" is only half true, and the disagreeing half is worse. Python lowercases and strips (annotations.py:294 `.strip().lower…
- **`resources.gnomad.oracle` is normalised two incompatible ways: Python maps every unrecognised spelling to `faf95`, the shell compares exactly** — `src/hprv/annotations.py:295` · fragility · verifiers confirmed/refuted/partly  
  *Verifier:* Two precision corrections to the finding as written; neither changes the verdict.

1. "The two normalisations disagree in both directions" is exact only for CASE and WHITESPACE variants (GRPMAX_PROXY, Grpmax_Proxy, "grpmax_proxy " — Python honours the proxy, t…
- **Step 3's classifier documents a permissive union over 'the looser of the dominant/recessive cutoffs' but hardcodes `recessive_max`, so raising `domina** — `src/hprv/selection.py:36` · doc mismatch · verifiers confirmed/partly/confirmed  
  *Verifier:* The headline, the location (`src/hprv/selection.py:36`) and both failure limbs reproduce exactly as written. One part of the stated failure scenario is incomplete and should be corrected before anyone acts on it:

The finding says Step 3's drop is what stops t…
- **Rung order makes `reason.spliceai` / `reason.cadd` count only the rescues that had no impact rung, so the audit under-reports the score evidence in th** — `src/hprv/selection.py:41` · unverifiable · verifiers partly/refuted/partly  
  *Verifier:* Two supporting claims are overstated; the central mechanism and all quoted counts are correct.

1. "reason.spliceai is the ONLY audit observable that the SpliceAI resource is live" is NOT accurate. pipeline/02_annotate_sites.sh:525-534 computes the VALUE-level…
- **CLAUDE.md:312-314 and docs/pipeline_design.md:228 promise every `bcftools annotate` transfer dies on a 0-match join; the third transfer, at 04:264, is** — `CLAUDE.md:312` · doc mismatch · verifiers confirmed/refuted/partly  
  *Verifier:* Two refinements to the finding as written.

(a) CITATION. `docs/pipeline_design.md:228` is NOT a repeat of the false claim. That bullet is explicitly scoped ("The same discipline covers the rest of **Step 2**") and its statement — "the ClinVar and gnomAD trans…
- **Step 1 strips INFO/END, so a symbolic ALT record enters the annotated union with no span** — `pipeline/01_make_cohort_sites.sh:165` · fragility · verifiers confirmed/refuted/confirmed  
  *Verifier:* The mechanism and the file:line (pipeline/01_make_cohort_sites.sh:165, `bcftools annotate -x INFO`) are correct, but the stated outcome understates it. The finding says the record reaches VEP as "a coordinate with no end" that VEP "treats as a 1 bp event". Wha…
- **Step 4's -R fast path and its whole-genome fallback are not equivalent: an indel whose left-alignment shift exceeds HPRV_PLAUSIBLE_PAD is dropped by o** — `pipeline/04_subset_and_annotate_trios.sh:232` · silent loss · verifiers confirmed/confirmed/confirmed  
  *Verifier:* Three precision refinements; none of them refute the finding.

1. `bcftools view -R` filters on the record's SPAN, not its start. Executed: with BED `chr1 5000 5100`, a record at chr1:4900 with a 200 bp REF IS returned, while chr1:4000 C>T is not. So the losin…

## What the completeness critic caught

After the first 92 findings were verified, a critic was asked what the fourteen lenses had *not*
looked at. It named five angles, each demonstrated by execution, and six follow-up finders turned
them into the round-2 findings above. Three of the five are larger than any single finding.

**INV-1 — The Step 3 → 4 → 5 annotation handoff is unguarded on both sides.** Every lens examined
annotations as parsers or as gates. Nobody checked the handoff: Step 3 classifies the annotated
union, Step 4 re-derives per-trio genotypes from the raw VCF and re-attaches annotations with a
blanket `bcftools annotate -c INFO`, and Step 5 resolves `rarity_af`/`rarity_basis` once for the
whole downstream from that re-attached copy. Step 3 carries a header guard for the faf95 witness
with a comment naming the exact catastrophe; Step 4 has no post-transfer coverage check and Step 5
has no guard at all. Executed: stripping the `gnomad_*` values from the real CH_A candidate VCF
while leaving the header declared (the shape of a transfer that matched nothing) made the run
**bigger** — 35 → 37 rows, every row `rarity_basis=absent`, exit 0, no warning. That column is the
one CLAUDE.md tells Step 9 to trust rather than re-derive, and here it asserts "gnomAD has no
record for this allele" for alleles nobody checked.

**INV-5 — Nothing verifies the trio VCF was jointly genotyped, and the one gate that could notice
is blind by construction.** Every lens assumed a non-carrier parent is an affirmative `0/0`. A
`bcftools merge` of single-sample callsets, or GenotypeGVCFs output with reference blocks dropped
in transit, gives `./.` instead. Executed: rewriting every parental `0/0` to `./.` in the real
CH_A VCF deleted 10 of 35 calls with no counter and took the de novo arm to zero; the survivors all
carry `origin_unverified`, which becomes the norm rather than the exception, and nothing measures
its rate. Worse, `00_qc.py:139` drops every site with a no-call parent from the MIE denominator,
so a merge-shaped trio reports a *cleaner* Mendelian-error rate. `qc_report.tsv` carries no
call-rate or missingness column at all.

**INV-4 — At a split multiallelic, the QC bands are applied to different quantities than at a
biallelic site.** `--keep-sum AD` was checked for arity, not meaning. Executed on a homopolymer
indel `AT > A,ATT,ATTT` with child `0/1:14,14,8,4`: after the exact Step-4 transformation the
stutter reads fold into REF, so allele balance reads 0.35 where the balance over alleles the child
carries is 0.50 (fail-closed: a perfect het drops below `het_ab_min` once the other alleles' reads
exceed twice the ALT depth), while `dp()` returns the site depth of 40, so `denovo_min_dp` is
satisfied by reads supporting neither of the child's alleles (fail-open).

**INV-2 — Read-backed phase is already in the input and is never read.** GATK emits physical phase
as `FORMAT/PGT` + `PID`; `bcftools norm -m-` carries it through Step 4 intact; nothing in `src/` or
`pipeline/` reads it. `docs/ROADMAP.md:61` lists read-backed phasing as "high cost, needs BAMs".
Three Step-5 decisions each saw only one side: `unphased_denovo_partner` declares a pair
unphaseable when the caller may have observed the phase; the `both` origin is never paired for the
same reason; and a pair whose trans evidence failed still consumes both legs where PID could
decide it. This is not a defect but the cheapest sensitivity gain in the audit.

**INV-3 — Adjacent same-codon variants (MNVs) are annotated wrongly in both directions.** VEP
annotates each SNV against the reference codon independently, so two missenses that jointly form a
stop never earn a HIGH rung (silent loss) and two "damaging" missenses whose joint codon is
synonymous are each kept (wrong call); both double-count as two observations, and a parental
miscall at one position can produce a spurious trans pair. The mock contains no adjacent-SNV pair,
so CI is structurally blind. `docs/limitations.md` says nothing about MNVs.

**Noted, not pursued:** `audit.py:37` is last-value-wins over an append-only `counts.tsv`, and Step 3
emits `reason.*` only for reasons that fired, so a `--from 3` re-run with looser thresholds leaves
the previous run's reason counts standing in `summary.md`. An audit-layer defect, not a filter one.

## Robustness assessment

The question was whether these two families are robust. They are not equally so, and the split is
informative: **the plausibility filter is sound on what it sees; the genotype-to-meaning assignment
is sound on what it matches.** Both qualifiers are load-bearing.

### The plausibility filter — structurally sound, with a single-transcript blind spot

The ladder logic holds up. Every boundary was probed and behaved as documented, the rung order is
deliberate and correct, the ClinVar override does what it claims, and BA1 genuinely wins over a
pathogenic assertion. It is also **the best-instrumented step in the pipeline**: it emits a reason
for every drop *and* every keep, and those reasons reconcile exactly against the file. If the rest
of the pipeline met this standard, this audit would be much shorter.

Three real weaknesses, none in the arithmetic, and the verifiers rated all three moderate:

- **It reads one transcript.** The impact rung consults only the selected CSQ block, so a
  `stop_gained` on an alternative transcript of the same gene, or a scored splice variant on a
  non-picked overlapping gene, is filed as `not_functional`. A verifier narrowed the lost class
  usefully: MODIFIER/LOW on the picked transcript *and* HIGH elsewhere, with neither score rung
  firing — real, but not "any non-MANE damaging consequence".
- **`not_functional` is one bucket for two facts** — "the scores said no" and "no score existed".
  A verifier pointed out the decomposition is recoverable from `cohort.sites.annotated.vcf.gz`, so
  the information is not lost, merely not reported.
- **Config can silently disable a rung.** `keep_impacts` as a YAML scalar removes the entire impact
  keep-path. Verifiers: real, reproduced in four YAML shapes, but a config-authoring slip only.

**Verdict: trustworthy for what it was handed, on the transcript it chose.** A negative result from
Step 3 means "no keep-path fired on the picked block", which is narrower than "not plausible".

### The genotype-to-meaning assignment — the weak family, and the reason is structural

The individual mode conditions are each defensible, and several are more careful than the field
norm (obligate transmission from a hom-alt parent, the father's chrX correctly ignored for a son,
`strict_gt` so a half-call is not read as hom-ref). The problem is not any single condition. It is
that **the conditions are conjunctive and total, with no fall-through.** A variant that fails one
conjunct of every mode matches nothing, and matching nothing produces no row, no reason and no
counter. There is no "examined but unassigned" state.

That architecture produces the same failure repeatedly, from unrelated causes — a het whose
*transmitting* parent fails its allele-balance band (the child's perfect call is deleted, while the
identical failure in the non-transmitting parent is emitted with a flag; one of the four surviving
majors); a `both`-origin het in the 1e-4..1e-2 band that pairs with nothing and clears no dominant
gate; a hom-alt child with one hom-ref or no-call parent — the deletion-in-trans and UPD shapes the
documentation tells the reader to look for; every gene-less inherited het; an entire trio whose VCF
lacks FORMAT/AD or GQ. Most of these the verifiers rated moderate: real, reproduced, and rarer than
the finders implied on GMKF-shaped data, or already partly documented.

Two findings are not about sensitivity but about **correctness of meaning**, and both survived at
major with three confirmations each. A trios file whose header spells the parent columns
`mother_id`/`father_id` — or any spelling outside the exact alias list — silently transposes the
parents, inverting every parent-of-origin call for that trio; and Step 0's Mendelian-error backstop
is *mathematically symmetric* under a father/mother swap, so the documented guard provably cannot
detect it.

**Verdict: individually careful, structurally fragile, and unaudited.** A negative result from
Step 5 today means "no mode's full conjunction was satisfied", and there is no way to distinguish
that from "no candidate variant was present".

### Three categories, kept separate

- **Real defects that should be fixed** — wrong or missing behaviour with no documentation: the
  parent-role transposition and its blind backstop, the transmitting-parent QC asymmetry, the
  ClinVar inert-band counter thresholding on the wrong gate, the unguarded Step-4 handoff, the
  unverified joint-genotyping assumption, `keep_impacts` and the boolean-knob coercions, the
  comp-het consumption rule not matching its own stated contract.
- **Documented accepted limitations** — 69 of 106 findings are already written down
  somewhere, which speaks well of the project's honesty. But several are documented *more narrowly
  than they bite*: the 1e-4..1e-2 band drop is scoped in the docs to a compound-het gene-key edge
  case when the real condition is simply "no trans partner in this gene", the default state of most
  surviving hets on real WGS.
- **Inherent limits of trio screening** — no code change fixes these; they belong in the methods:
  phase for a de-novo partner is unresolvable from trio genotypes alone (though PGT/PID may
  resolve some — INV-2); a 1/1 child with a 0/0 parent is genuinely ambiguous between error,
  deletion and UPD; proband mosaicism sits outside any allele-balance band that also excludes noise.

### What would have to be true to trust a negative result

1. Step 5 emits `variants_examined` and a per-reason no-row tally; Step 6 emits `calls_in` and
   `calls_no_gene`. Until then the funnel cannot be reconciled, and comp-het row inflation makes it
   read as a gain.
2. Step 5 asserts, on the file it actually reads, the same witness Step 3 asserts on the union —
   and Step 0 reports a per-member call rate.
3. The parent-role and config-coherence traps are closed, since both corrupt every call in a trio
   rather than one call.
4. `not_functional` is split into "scored and rejected" and "no score available".
5. Step 9's input is checked for freshness against `candidates.calls.tsv`.

## Recommended order of work

Grouped by what the change buys. Effort: S = a few lines, M = a function, L = a component.

**Tier 1 — make the existing losses visible.** The cheapest work in the audit; it converts the whole
silent-loss class into an auditable one without changing a single call.

| | Change | Effort |
|---|---|---|
| 1 | Step 5: `variants_examined` plus a per-reason no-row tally (`no_mode`, `qc_fail`, `rarity_band`, `no_gene`, `sex_unresolved`, `chrY`) | S |
| 2 | Step 6: `calls_in` and `calls_no_gene` counters | S |
| 3 | Widen `clinvar_plp_dropped_ge_recessive_max` to the dominant gate — the counter was written for exactly this and misses the half that fires | S |
| 4 | Step 4: audit per-field lift coverage after the `-c INFO` transfer; Step 5: assert the faf95 witness on the per-trio file (INV-1) | S |
| 5 | Step 0: per-member genotype call rate / `./.` fraction, with a refuse threshold (INV-5) | S |
| 6 | Step 1: record the pre-filter count; Step 4: record the pre-`isec` count | S |
| 7 | Split `not_functional` into scored-and-rejected vs no-score-available | S |

**Tier 2 — close a wrong call.** These change results, so each needs a test alongside.

| | Change | Effort |
|---|---|---|
| 8 | Parent-role transposition: match parent columns strictly and fail loudly on any other header; run `scan_sex` on all three members so a male in the mother slot is caught | S |
| 9 | Transmitting-parent QC failure: emit with a flag rather than delete, matching how the non-transmitting parent is already handled | S |
| 10 | Comp-het consumption: only a pair with *measured* trans evidence should consume its legs, matching the documented rule | S |
| 11 | Validate `keep_impacts`, the rarity-cutoff coherence, the AB band ordering, and the boolean knobs; halt on an incoherent set | M |
| 12 | Step 9: compare its input's row count against `candidates.calls.tsv` and refuse a stale table | S |

**Tier 3 — close a silent loss.** Each recovers real variants and needs a decision, not just code.

| | Change | Effort |
|---|---|---|
| 13 | Emit an unassigned-call row (never-drop applied to Step 5) for a variant that failed a *single* named conjunct | M |
| 14 | Read `PGT`/`PID`: phase `both`-origin and de-novo-partner pairs where the caller observed it (INV-2) | M |
| 15 | Index each variant under every gene it hits for the pairing key, not only the picked block | M |
| 16 | Read the impact and SpliceAI rungs across all CSQ blocks | M |
| 17 | Compute allele balance and depth over ref + this-alt only at split multiallelics (INV-4) | M |
| 18 | hom-recessive with one non-carrying parent: emit flagged as possible deletion-in-trans or UPD | M |

**Not recommended.** Do not add a fall-through that emits every examined variant — that inverts the
screen. Item 13 should emit only variants that reached Step 5 and failed a single named conjunct.

## Verified sound

Checked deliberately and found correct. This section is what makes the rest trustworthy: the sweep
was not looking only for problems.

- **The never-drop chain in the variant layer.** `candidates.calls.tsv`, `igv/variants.tsv` and
  `variants.prioritized.tsv` reconcile exactly at 41 rows. `igv.build_variants_tsv` is strictly 1:1
  and stayed so when handed a deliberately broken manifest. Both Step 9 assertions read an
  unfiltered row count, so they are not vacuous. `report.py` truncates nothing.
- **`bcftools +split-vep` conserves records** — 4 in, 4 out, including a record with no CSQ, an
  intergenic record, and a record with no PICK flag.
- **Haploid chrX genotypes.** `GT=1` reads as `HOM_ALT` under cyvcf2 and routes correctly to
  `x_linked_recessive`.
- **`--keep-sum AD`** is correctly arity'd on a multiallelic (what it *means* is INV-4).
- **Obligate transmission.** A hom-alt parent yields a deterministic parent-of-origin, and 1/1 × 1/1
  with a het child is correctly refused as a Mendelian error.
- **The father's chrX is correctly ignored for a son**, and an affected father is flagged rather
  than used to drop the call.
- **`strict_gt=True`** genuinely makes a half-called `0/.` parent a no-call rather than hom-ref.
- **The Step-4 FILTER asymmetry causes no loss against policy** — a verifier applied the proposed
  fix and got a byte-identical `candidates.calls.tsv`.
- **The hiConfDeNovo gate does not re-impose a paternal requirement on the male-X path** — a
  verifier emitted `denovo_x_hemi` with the tag present; the finder's claim was overstated.
- **The documented maternal-no-call chrX gap reproduces exactly as written.**
- **Step 3's reason accounting reconciles exactly** on the real run.
- **The rarity chokepoint's arms never cross**, under both oracle settings, including `zero_ci` and
  `absent`.
- **`PGT:PID` and the `|` separator survive `bcftools norm -m-`** end to end, and
  `is_hom_alt_call` splits on `[/|]` — the phase information is intact, just unread.

## Provenance and confidence

| | |
|---|---|
| Findings raised | 112 (92 first sweep + 20 critic follow-up), 106 distinct after merging overlapping reports |
| Survived verification | all 106 — the rule was two of three verifiers voting confirmed *or partly*, and fewer than two voting not-a-finding |
| Had one verifier vote to refute | 26 |
| Severity, finder vs verifiers' median | major 39 → **4**, moderate 52 → **46**, minor 15 → **56** |
| Reproduced by executing code | 101 of 106 |
| Already documented somewhere in the repo | 69 |
| Not documented anywhere | 37 |

**Read the survival rule for what it is.** "Partly confirmed" counts toward survival, so this
filter keeps a finding whose mechanism is real even when its stakes were overstated — which is why
nothing was refuted outright while 26 findings carry a dissent and 35 of
39 finder-majors were talked down. The dissents are quoted under each finding. If you
want the strict list, it is the four majors plus the 42
moderates with no dissent.

Every finding was demonstrated by constructing an input and running the real function or the real
pipeline step against the mock genome and its outputs, not by reading the code — a deliberate
departure from the 2026-07 review, which recorded that no code was executed. Two caveats: no VEP
was available to the agents, so the findings that depend on what VEP annotates onto a
spanning-deletion `*` allele are partly reasoned; and the mock genome cannot produce every real
input shape (no adjacent-SNV pair, no merge-shaped trio), so INV-3 and INV-5 rest on synthetic
edits to the real fixture rather than on a naturally occurring case.

## Scope

This audit deliberately excluded the ranking layer. Step 9's scoring, the artifact panel, the tier
assignment and the composite were reviewed in `docs/pipeline_review_2026-09.md` and are out of
scope here, except where they consume a Step-5 output incorrectly. Nothing was implemented; this is
a consolidation pass, as requested.
