# Methods: screening non-jointly-genotyped trio callsets for inherited rare variation with hprv

*A methods description derived from the pipeline code and shipped configuration
(`config/config.example.yaml`) as of 2026-09-08 (the git history records the exact commit). Every threshold quoted is the
shipped default and is configurable; the implementing script for each step is listed in Table 3.
A rendered copy is at `docs/methods.pdf`.*

## 1. Overview and design principles

hprv is a container-based pipeline that screens per-trio germline VCFs for high-priority
inherited rare variants and consolidates genes in which rare functional variants recur across
unrelated probands. It was written for the Gabriella Miller Kids First (GMKF) per-trio callsets:
GRCh38 VCFs produced by the GATK genotype-refinement workflow, in which each trio was genotyped
independently and the cohort was never jointly genotyped. That property shapes the design. Because
a locus absent from one trio's VCF cannot be read as homozygous reference, the pipeline never
constructs a cohort genotype matrix and never derives allele frequency from internal counts;
population frequency is taken exclusively from gnomAD v4.1, as a single quantity fixed for the whole
run (Section 7). Cross-trio information is used only as a recurrence count of distinct individuals
(Section 12) and, in the prioritization layer, as a quality signal (Section 13).

The analysis proceeds through a trio-resolution preflight and ten steps (Steps 0 to 9). Variant
loci are unioned across trios into a site-only list (Step 1), annotated once (Step 2), reduced to
biologically plausible sites (Step 3), and then the true per-trio genotypes are recovered at those
sites (Step 4) and classified by inheritance mode with genotype-level quality control (Step 5),
then re-scored for splicing at a wide window (Step 5b). Calls are consolidated per gene (Step 6), summarised for review (Steps 7 and 8), and re-ranked by a
gene-level excess statistic, an artifact panel, and a per-variant evidence tier (Step 9).

Three rules govern every stage. First, rarity, functional impact, and genotype quality are gated
before, and independently of, any gene list or constraint metric; gene lists and constraint act only
as ranking priors, so a variant in a gene with no prior annotation is never removed (the
"never-drop" rule). Second, one population-frequency quantity is used for every variant in a run;
the run-level choice and the per-variant provenance are written into the output. Third, every step
is idempotent with content-keyed completion markers, halts on an incoherent configuration or a
transfer that matched nothing, and writes its input and output counts to a run audit
(Section 15).

The scope is inherited germline single-nucleotide and small indel variation under dominant,
recessive (homozygous and compound heterozygous in trans), and X-linked models. De novo variants
are detected only as a secondary cross-reference; mitochondrial contigs are excluded at the union
step; copy-number and structural variants are not assessed.

## 2. Software environment and reproducibility

All tools run from one container image built on the Ensembl VEP release 115.0 image
(`ensemblorg/ensembl-vep:release_115.0`), onto which a version-pinned micromamba environment is
layered: bcftools, samtools and htslib 1.23, bedtools 2.31.1, and Python 3.11 with cyvcf2 0.31.1,
pysam, NumPy 2, pandas 2, SciPy, PyYAML and openpyxl. The conda Perl is removed after the solve so
that VEP and its plugins run under the base image's Perl. For the optional non-human-read screen
(Section 14.3), kraken2 v2.17.1 is compiled from source and the `nonhuman-screen` package is
installed at a pinned commit (`9e9eb6ce`). An isolated conda environment holding the Illumina
SpliceAI model (TensorFlow 2.13 to 2.15, Keras 2) supports the optional live splice-score backfill
(Section 6.4). The pipeline is intended to run inside the image under Apptainer on HPC; the same
scripts detect a host invocation and wrap each tool call in Apptainer or Docker.

A single YAML configuration file is the contract between the user and the pipeline. Filesystem
locations are `${ENV}` placeholders resolved at run time, so no path, sample identifier, or
controlled-access artifact is ever committed. The YAML is loaded and re-emitted as shell exports
for the Bash steps, with booleans normalised to `true`/`false`.

The expensive steps (trio QC, the site union, annotation, per-trio genotype recovery, CRAM slicing
and prioritization) write completion markers whose content is a key over the inputs that
determine the output, and are skipped on re-run when the key matches; the inexpensive Python steps
(selection, inheritance classification, gene consolidation and the workbook) recompute every run.
The keys are: for the site union, the manifest plus each trio VCF's path, size and modification time,
the FILTER expression, the excluded contigs and the sample subset; for annotation, a checksum of
the union plus the identity of the cache and each plugin score file, the two transfer sources, the
transcript selector, the frequency oracle, the VEP flag recipe and the list of lifted fields; for
per-trio candidate VCFs, a checksum of the
plausible-site file plus the source VCF identity and sample subset; for prioritization, a SHA-256
over every input table, the configuration and the scalar arguments. A marker whose key differs
from the current inputs invalidates the cached product, and every bgzip product is
integrity-checked (`bgzip -t`) before it is trusted.

At the start of Steps 3 and 5 the screen's parameters are validated and the run halts on any
setting that would silently disable a filter: an unrecognised or empty impact list, a rarity ladder
that is not ordered (`dominant_max` ≤ `recessive_max` ≤ `benign_ba1`, and `recessive_strict` ≤
`recessive_max`), an allele-balance band outside [0, 1] or with its lower bound at or above its
upper bound, a de novo depth floor below the general depth floor, or a boolean given as a quoted
string.

## 3. Inputs and trio resolution

The user supplies a tab-separated trios file and a VCF source (a directory searched recursively for
`*.vcf.gz`, `*.vcf` and `*.bcf`, a file listing VCF paths, or both). The trios file must carry a
header naming the proband, father and mother columns; names are matched case-insensitively against
the aliases kid/child/proband, dad/father/paternal and mom/mother/maternal, with `_id`,
`_sample_id` and `sample_` affixes tolerated. Column order is never assumed when a header is
present, and a header that names some roles but not all three is an error rather than a positional
guess, because a transposed mother and father would invert every parent-of-origin call while
leaving the Mendelian-error rate unchanged. A file with no header is read positionally with a
warning. Duplicate proband rows are skipped.

Sample membership is read from each VCF header only. For each trio, the candidate VCFs are those
containing all three sample identifiers by exact match; a VCF may contain additional samples. When
more than one VCF qualifies, the one with the fewest samples is chosen, with a lexical tie-break on
path, and the multiplicity is recorded. A trio with no qualifying VCF is reported as unresolved in
`trio_resolution.tsv` (naming the missing members) and skipped; the run fails only when no trio
resolves. For each resolved trio a PED file is generated with the proband's sex unknown (`0`; it is
inferred in Step 0), the father coded male, the mother female, and the proband marked affected. The
resolved manifest (`trios.resolved.tsv`: trio identifier, VCF, PED, comma-separated samples) is the
input to Steps 0, 1 and 4, and the number of resolved trios is the screened-cohort size *N* used in
Steps 6 and 9.

## 4. Per-trio quality control (Step 0)

Step 0 computes an advisory set of per-trio quality flags. No step excludes a flagged trio
automatically; the flags are surfaced in the workbook and the review export for the curator.

For each trio, the VCF is scanned over autosomal biallelic records until 200,000 sites have passed
the site-level QC below (`qc.max_sites`); on genome-wide data the statistics therefore describe a
prefix of the genome, adequate for detecting sample swaps and gross contamination but blind to a
localised Mendelian-error cluster beyond it. A site enters the Mendelian-error denominator when all
three members are called and each has genotype quality (GQ) ≥ 20 and depth ≥ 10. A violation is
scored when a homozygous-alternate child has a homozygous-reference parent, a homozygous-reference
child has a homozygous-alternate parent, or a heterozygous child has two homozygous-reference or two
homozygous-alternate parents. The trio is flagged when the violation rate exceeds 0.02. The
per-member no-call rate over the scanned biallelic autosomal records is also reported and flagged
above 0.10, because a merged single-sample callset renders every non-carrier parent as missing
rather than homozygous reference and would otherwise pass the Mendelian test with a cleaner-looking
rate while defeating every parental test in Step 5.

Sex is inferred for all three members from chrX outside the pseudoautosomal regions (GRCh38 PAR1
chrX:10,001–2,781,479; PAR2 chrX:155,701,383–156,030,895). Among that individual's biallelic chrX
calls with GQ ≥ 20 and depth ≥ 10, a heterozygous fraction het/(het + hom-alt) below 0.10 is called
male and otherwise female, provided at least 20 informative sites are available (else unknown).
A father inferred female or a mother inferred male raises `parent_sex_flag`, the one direct
detector of transposed parental columns. The proband's inferred sex is compared with the PED only
when the PED states one.

Contamination is assessed from verifyBamID `FREEMIX` when a directory of `.selfSM` files is
configured (flag above 0.05). Otherwise a VCF-only proxy is used: for each member, the summed
reference-allele depth divided by the summed total depth over that member's homozygous-alternate SNV
sites with GQ ≥ 20 and depth ≥ 10 (flag above 0.02). This proxy is uncorrected for allele frequency
and reference bias and detects only gross contamination; the calibrated CHARR estimator is not
implemented. `overall_pass` is set when no flag fired.

## 5. Cohort site union (Step 1)

For each resolved trio the source VCF was reduced to its three members (`bcftools view -s`),
restricted to records with FILTER `PASS` or `.`, and stripped of the mitochondrial contigs
(`chrM`, `chrMT`, `M`, `MT`). Multiallelic records were split and indels left-aligned against the
reference (`bcftools norm -m- -f`, with a reference mismatch treated as fatal), alleles no member
carries after splitting were removed (`--min-ac 1`), genotypes were dropped (`view -G`) and all INFO
fields were removed. The per-trio site files were then concatenated with duplicate removal
(`bcftools concat -a -D`), sorted, and collapsed to unique representations (`norm -d exact -f`),
yielding one site-only cohort union. The union is verified to contain no record on an excluded
contig. Mitochondrial exclusion is enforced here rather than left unmodelled because every
inheritance model downstream is diploid, and homoplasmic rCRS-referenced haplogroup variants, which
carry no gnomAD frequency, would otherwise be called as homozygous recessive in every trio.

## 6. Variant annotation (Step 2)

### 6.1 Ensembl VEP

The cohort union was annotated once, never per trio, with Ensembl VEP release 115 against the
matching offline GRCh38 cache (`--cache --offline --assembly GRCh38 --fasta`), sharded by contig for
resumability and re-assembled before any downstream operation, so the product is identical to a
single pass. VEP was run with `--symbol --biotype --numbers --total_length --hgvs --canonical --mane
--af_gnomade --af_gnomadg --max_af --check_existing --flag_pick --pick_order
mane_select,mane_plus_clinical,canonical,rank`. The gnomAD v4.1 per-population exome and genome
allele frequencies and the ClinVar clinical significance therefore come from the VEP cache itself.
Five plugins were applied. Ensembl's stock NMD plugin, which needs no data file, marks a stop-gain,
frameshift or canonical splice variant predicted to escape nonsense-mediated decay (last exon,
within 50 nt of the penultimate exon's end, first 100 coding bases, or an intronless transcript).
The four score plugins were CADD v1.7 (the whole-genome SNV table and the gnomAD-genomes indel
table), SpliceAI (the precomputed raw genome-wide SNV and indel delta-score files, Illumina
`genome_scores_v1.3`; the plugin's own cutoff is deliberately not set so raw scores are carried),
REVEL v1.3 (the GRCh38-sorted distribution), and AlphaMissense (`AlphaMissense_hg38.tsv.gz`,
with `transcript_match=1` so a score is taken only from the transcript VEP annotated).

Because `--flag_pick` retains every consequence block and marks VEP's choice with `PICK=1`, the
consequence used downstream is selected when the CSQ fields are lifted to INFO with
`bcftools +split-vep`: the picked block when a PICK field exists, otherwise the worst consequence
(`resources.vep.csq_select` overrides this). The lifted fields, each prefixed `vep_`, are
Consequence, IMPACT, SYMBOL, Gene, Feature, BIOTYPE, EXON, INTRON, HGVSc, HGVSp, cDNA_position,
CDS_position, Protein_position (each with its total), MANE_SELECT, NMD, CADD_PHRED, CLIN_SIG, gnomADe_AF, gnomADg_AF, MAX_AF, MAX_AF_POPS, the ten grpmax-eligible population
frequencies (gnomADe and gnomADg for AFR, AMR, EAS, NFE and SAS), the SpliceAI delta scores and
positions (DS_AG, DS_AL, DS_DG, DS_DL, DP_AG, DP_AL, DP_DG, DP_DL, SYMBOL), REVEL,
am_pathogenicity and am_class.

The step halts if Consequence, IMPACT or SYMBOL is missing, if none of the ten grpmax-eligible
frequency fields is present, if SpliceAI or CADD is configured but not one site received a score,
or if no site carries any gnomAD frequency value; each of these would otherwise leave a run that
exits successfully with a silently disabled gate.

### 6.2 The two `bcftools annotate` transfers

Exactly two annotations are transferred from external VCFs, both because the VEP cache cannot supply
them. From the ClinVar GRCh38 sites VCF (pinned release 2026-07-06, contigs renamed to the `chr`
prefix during preparation), `CLNREVSTAT` and `CLNSIG` are transferred as `clinvar_CLNREVSTAT` and
`clinvar_CLNSIG`; the review status is later mapped to gold stars (practice guideline 4, expert
panel 3, multiple submitters without conflict 2, single submitter or conflicting 1, no assertion
criteria 0). From a slim of the gnomAD v4.1 joint (exomes plus genomes) sites VCF, prepared by
retaining five INFO fields (`AF_joint`, `AF_grpmax_joint`, `fafmax_faf95_max_joint`,
`fafmax_faf95_max_gen_anc_joint`, `nhomalt_joint`), those fields are transferred as
`gnomad_AF_joint`, `gnomad_AF_grpmax`, `gnomad_faf95`, `gnomad_faf95_group` and `gnomad_nhomalt`.
Both transfers are allele-exact (chromosome, position, reference and alternate) and each is guarded:
a transfer that matched zero sites is fatal, because a cohort union always overlaps both resources
and zero matches indicate a build or contig-naming mismatch. The gnomAD transfer is required under
the default frequency oracle (Section 7); the ClinVar transfer is optional and its absence leaves
the review status "unavailable", which is distinguished from zero stars everywhere downstream.

### 6.3 Ingesting an external VEP VCF

A pre-existing VEP output may be supplied in place of running VEP. It must carry a `##VEP` header
and must not declare a non-GRCh38 assembly; it must contain no multiallelic records, and it must
cover every site of the cohort union (uncovered sites would be silently lost at Step 4), unless a
partial ingest is explicitly permitted and the shortfall recorded. The same field lifting, guards
and transfers then apply.

### 6.4 Optional live SpliceAI backfill (Step 2b)

Off by default. When enabled, variants in the annotated union with no precomputed SpliceAI score
(indels only, by default, because the precomputed SNV set is genome-wide) are scored with the
Illumina SpliceAI model (`-A grch38 -D 500`) in the isolated environment. For each variant the
maximum of each delta score over the genes reported is written into the same
`vep_SpliceAI_pred_DS_*` fields, so Step 3 reads precomputed and backfilled scores identically. An
unavailable environment halts the run; a transient scoring failure leaves the union on precomputed
scores.

## 7. Population allele frequency: one oracle per run

All rarity decisions read one function, `annotations.frequency()`, which returns a single quantity
selected for the whole run by `resources.gnomad.oracle`.

Under the default, `faf95`, the value is gnomAD's filtering allele frequency
`fafmax_faf95_max_joint`: the lower bound of the 95% Poisson confidence interval on the allele
frequency, maximised over the FAF-eligible genetic-ancestry groups (afr, amr, eas, mid, nfe, sas),
which is the quantity the ACMG/AMP and ClinGen frequency-filtering recommendations specify. gnomAD
emits this field as missing, not zero, wherever no group's interval lower bound exceeds zero. The
pipeline therefore resolves an allele in three ways and records which applied in `rarity_basis`:
a published faf95 is used as is (`measured`); an allele present in the gnomAD record but with no
faf95 resolves to 0.0 (`zero_ci`, the rarest value); and an allele with no gnomAD record at all
resolves to missing (`absent`), which every gate also treats as rarest. The witness separating the
last two cases is the transferred `gnomad_AF_joint` field. The VEP-cache point estimates are never
consulted on this arm.

Under the explicit opt-down, `grpmax_proxy`, the value is the maximum of the VEP-cache point
estimates over the grpmax-eligible groups only (AFR, AMR, EAS, NFE, SAS, exome and genome), which
mirrors gnomAD's own exclusion of the bottlenecked ami, asj, fin and mid groups. The cache carries no
allele counts, so no interval correction is possible; the proxy sits about one interval width above
faf95 for low-count alleles and errs toward removing them. Cache frequencies exist only for
dbSNP-accessioned alleles, so an unaccessioned gnomAD allele reads as absent under this arm.

VEP's `MAX_AF` (which maximises over founder and 1000 Genomes populations) and the global gnomAD
frequencies are carried for reporting only and are not filter fields under either arm. Where a
field holds several values (multiple consequence blocks), the maximum is taken. Step 5 resolves the
value once per call and writes `rarity_af`, `rarity_oracle` and `rarity_basis`; Steps 6 and 9
consume those columns rather than re-deriving them. Under `faf95`, Step 5 also counts calls with no
gnomAD record that nevertheless carry a cache frequency; because the joint release is a superset of
the cache, a non-zero count indicates an under-matched transfer and is reported.

## 8. Selection of plausible sites (Step 3)

Step 3 is an inheritance-agnostic reduction of the annotated union applied in a fixed order. A
site with frequency ≥ 0.05 (the ClinGen BA1 threshold) is removed and never rescued, including by a
ClinVar assertion. A site is otherwise retained on rarity if its frequency is missing, is below the
permissive recessive cutoff 1 × 10⁻², or the site carries a ClinVar pathogenic or likely-pathogenic
assertion (`CLIN_SIG` containing "pathogenic" without "conflicting" or a benign qualifier; review
status is deliberately not consulted). A pathogenic or likely-pathogenic site is then kept outright.
Every other rare site must pass a three-rung functional ladder, any rung sufficing, evaluated in
order: VEP IMPACT of HIGH or MODERATE; maximum SpliceAI delta score (over acceptor and donor gain
and loss) ≥ 0.2; CADD PHRED ≥ 25.3. A missing SpliceAI or CADD score never removes a site; it merely
fails to retain it. Each retained site is tagged with its keep reason (`hprv_keep_reason`, one of
clinvar_plp, impact_high, impact_moderate, spliceai, cadd), and the audit records the drop reasons
(ba1, too_common, not_functional), splitting not_functional into sites both predictors scored and
rejected versus sites neither scored.

Two consequences of the ladder are worth stating. Every missense variant has IMPACT MODERATE and is
retained at the first rung, so REVEL and AlphaMissense cannot change which variants pass the
screen; they act only in Step 9. Conversely, CADD 25.3, a cutoff calibrated for missense, is in
practice applied only to variants below MODERATE impact (intronic, synonymous, UTR and regulatory),
where it functions as a discovery rank rather than as calibrated evidence. Gene lists and
constraint play no part in this step.

## 9. Per-trio genotype recovery (Step 4)

The true genotypes of each trio at the plausible sites were recovered from the source VCF rather
than from any merged product. The plausible loci were padded by 1,000 bp and merged into a BED, and
each trio VCF was restricted to those regions and its three members, normalised as in Step 1
(`norm -m- -f`, this time warning rather than failing on a reference mismatch) with `--keep-sum AD`
where an AD field exists, so that after splitting a multiallelic heterozygous genotype each allele's
reference depth includes the reads supporting the other alternate allele and allele balance remains
interpretable. Alleles the trio does not carry were removed, and all INFO was stripped except the
GATK `hiConfDeNovo` and `loConfDeNovo` tags. The result was intersected allele-exactly with the
plausible set (`bcftools isec -c none -n=2 -w1`) and the Step 2 annotations were transferred onto it
(`bcftools annotate -c INFO`) from a copy of the plausible set with the de novo tags removed, so the
transfer can never overwrite a trio's own de novo tag. The transfer is verified per trio: every
candidate record must carry `hprv_keep_reason`, otherwise the step halts.

## 10. Genotype quality control

Genotype-level QC is evaluated on the refined genotypes with cyvcf2 in strict mode, so that a
half-called genotype (`0/.`) is treated as missing rather than homozygous reference. GQ is read from
FORMAT/GQ, which in genotype-refined VCFs is derived from the posterior probabilities; depth is
the sum of allelic depths when AD is present, otherwise FORMAT/DP; allele balance (AB) is the
alternate-allele fraction of AD.

Five predicates are used (defaults from `filters.genotype_qc` and `filters.denovo`). A
heterozygous carrier requires GQ ≥ 20, depth ≥ 10 and AB within [0.25, 0.75]. A homozygous-alternate
carrier requires GQ ≥ 20, depth ≥ 10 and AB ≥ 0.90. A de novo child requires GQ ≥ 20, depth ≥ 20
and AB within [0.25, 0.75]. A homozygous-reference non-carrier requires GQ ≥ 20, depth ≥ 10 and
AB ≤ 0.10. A "clean" parent requires GQ ≥ 20, depth ≥ 10, at most one alternate-supporting read and
AB ≤ 0.10. The carrier predicates fail when AD is absent; the two non-carrier predicates pass when
AD is absent, because a GATK reference-block genotype (`GT:DP:GQ:MIN_DP:PL`, no AD) is the ordinary
shape of a non-carrier parent and failing it would discard every such trio. Such a pass is marked
(`parent_ad_unmeasured`, `trans_evidence_unmeasured`) so that a phase or de novo call resting on a
genotype call alone is distinguishable from one resting on read counts. Records with a FILTER value
other than PASS or `.` are skipped.

## 11. Inheritance-mode classification (Step 5)

Each trio's candidate VCF is screened with the proband's sex taken from the PED or, when unknown
there, from Step 0's inference; if the sex remains unknown, non-pseudoautosomal X and Y records are
skipped for that trio with a warning and autosomal modes still run. Non-PAR chrY records in a female
proband are skipped. A "rare" test at limit *L* passes when the oracle value is missing or below
*L*; the limits are 1 × 10⁻⁴ for the dominant and de novo models (`dominant_max`) and 1 × 10⁻² for
the recessive and X-linked models (`recessive_max`), with recessive and X-linked calls below
1 × 10⁻³ additionally flagged `high_conf_rarity`. Every emitted row carries the genotype of all three
members as allele strings, the child's GQ, depth and AB, the curated annotations (the HGVS coding
and protein descriptions, exon and intron numbering, CDS position and MANE status of the transcript
selected in Step 2, the NMD verdict of Section 13.4, and the SpliceAI event decomposition below) and
a `flags` field;
after the curated columns, every INFO field declared in the per-trio candidate VCF header is
written verbatim under an `info_` prefix (the union over trios, in header order), including the
complete VEP consequence string for every transcript, so no annotation computed in Step 2 is
lost in the projection from VCF to table.

**De novo (secondary).** An autosomal de novo requires a heterozygous child with both parents
homozygous reference, the child passing the de novo predicate, both parents passing the clean-parent
predicate, and rarity at 1 × 10⁻⁴. A male X-hemizygous de novo requires a homozygous-alternate son
with a homozygous-reference mother passing the clean-parent predicate (the father's X is not
transmitted to a son and is not consulted), the son passing the homozygous predicate with depth
≥ 20, and the same rarity. When the callset header declares `hiConfDeNovo`, the call additionally
requires that tag to list this child (`filters.denovo.use_hiconf_tag`); when the header lacks the
tag, detection rests on genotypes and QC alone. De novo rows carry `review_prior_crosscheck`,
a reminder that genotype refinement with population priors can suppress an ultra-rare call, and are
never the pipeline's primary result.

**Homozygous recessive.** An autosomal homozygous-alternate child with both parents carrying the
allele (heterozygous or homozygous), the child passing the homozygous predicate, each parent passing
the predicate matching its own genotype, and rarity at 1 × 10⁻². A homozygous child with a
homozygous-reference parent is not called (a deletion in trans or uniparental disomy is outside the
model) and is counted as Mendelian-inconsistent.

**X-linked recessive.** For a male proband, a homozygous (hemizygous) alternate call at a non-PAR
chrX site with a heterozygous or homozygous mother, both passing QC, and rarity at 1 × 10⁻²; a
father carrying the allele does not veto the call but is flagged (`father_carries_x_allele`). For a
female proband, a homozygous-alternate call with a carrier mother and a hemizygous
(homozygous-alternate) father, all three passing QC, and the same rarity. Non-PAR chrY is routed
to no inherited model.

**Heterozygous collection and parent of origin.** Every heterozygous child call outside male
hemizygous regions that passes the heterozygous predicate, is rare at 1 × 10⁻², and has a gene
annotation is pooled per gene with a parent-of-origin assignment. With both parents homozygous
reference the origin is "de novo" only if both pass the clean-parent predicate. When exactly one
parent carries the allele, the origin is that parent (maternal or paternal); the non-transmitting
parent must be an affirmative, QC-passing homozygous reference for the origin to count as verified,
otherwise the call is flagged `origin_unverified`; a verified non-transmitting parent with no allelic
depth is flagged as a vacuous pass. When both parents carry the allele, a homozygous parent
transmits obligately, so a heterozygous child of a homozygous and a heterozygous parent is assigned
to the homozygous parent; two heterozygous parents give the indeterminate origin "both"; two
homozygous parents cannot produce a heterozygous child and the record is counted as
Mendelian-inconsistent. If the transmitting parent fails its own genotype QC the call is retained
and flagged `transmitting_parent_qc_fail` rather than removed. A parent no-call leaves the origin
unestablished.

**Compound heterozygous.** Within each gene, trans pairs are formed between maternal and paternal
hets, and between a de novo het and a het of either parental origin; hets of origin "both" are
never paired. Each pair receives a `pair_id` and both legs are emitted as `compound_het`. A pair
with a de novo leg cannot be phased from trio genotypes and is flagged `unphased_denovo_partner`; a
pair whose trans evidence failed is flagged `origin_unverified`; a pair with a QC-failed
transmitting parent is flagged `transmitting_parent_qc_fail`; a pair whose trans evidence rested on a
parent without allelic depth is flagged `trans_evidence_unmeasured`. Only a pair carrying none of
the first three flags is treated as a confirmed biallelic hit and consumes its two legs from the
dominant model below.

**Dominant.** Every pooled het of maternal, paternal or indeterminate ("both") origin that was not
consumed by a confirmed compound-heterozygous pair and is rare at 1 × 10⁻⁴ is emitted as
`dominant` with `origin=mat|pat|both` and any of the flags above. This inherited heterozygous class
is the recurrence signal Step 6 consolidates. A variant may appear under more than one mode.

**Accounting.** Every record examined is classified as skipped (FILTER, unresolved sex, female
chrY), called, or assigned a no-row reason (child not a carrier, chrY, a male X het, child or parent
QC failure, rarity, no gene, Mendelian inconsistency, parental no-call, the inert
[1 × 10⁻⁴, 1 × 10⁻²) heterozygous band, a disabled mode, or the hiConf tag), so that per trio the
number examined equals skipped plus called plus no-row. ClinVar P/LP alleles the child carried that
yielded no row are counted separately.

**Splice event decomposition.** The screen gates on the maximum of SpliceAI's four delta scores;
each call also carries the event behind that maximum. The component with the largest delta score at
or above the screen's 0.2 floor is named (`spliceai_event`: acceptor or donor, gain or loss) with
the genomic position of the affected site (the variant position plus the model's offset), and a
second component at or above the floor is named alongside it. A gain and a loss of the same site
type are read as a cryptic-site shift whose exon-boundary displacement is the difference of the two
offsets and whose frame follows from that displacement modulo three; a lone loss is recorded as a
site loss (exon skipping or intron retention, frame undetermined without the exon length), a lone
gain as a site gain, two gains as a pseudoexon candidate, two losses as a whole-exon loss
candidate, and mixed site types as complex (`spliceai_effect`). A call whose SpliceAI gene symbol
differs from VEP's picked gene is flagged.

**Wide-window rescoring (Step 5b).** The precomputed SpliceAI files were generated with a 50 bp
window, so a cryptic site or pseudoexon partner further than 50 bp from the variant is invisible to
the screen. Every distinct variant in the calls table was therefore re-scored live with the SpliceAI
model at a 4,999 bp window (the window used by the ClinGen splicing recommendations), in chunks of
1,000 variants that run as independent scheduler array tasks, with results cached per variant so a
re-run of Step 5 scores nothing anew. The gene entry with the largest delta score is kept whole and
decomposed exactly as above into `spliceai_wide_*` columns, together with the wide-window maximum,
its difference from the precomputed maximum, and a flag marking an event more than 50 bp from the
variant. These columns are evidence for review; the selection gate and the Step 9 tier continue to
read the precomputed score so that results remain comparable between runs.

## 12. Gene-level consolidation (Step 6)

Calls were aggregated per gene symbol (falling back to the Ensembl gene identifier), counting the
number of distinct probands under each model: `n_dominant` (dominant rows), `n_biallelic`
(homozygous-recessive and compound-heterozygous rows), `n_xlinked` and, separately, `n_denovo`.
The headline carrier count `n_carriers` is the union of the three inherited sets, and a gene is
flagged `recurrent` at ≥ 2 distinct probands (`burden.min_carriers`). Recurrence is classified as
`same_variant` when every carrier's set of qualifying variant keys is identical (one shared allele,
or one shared compound-heterozygous pair) and `distinct_variant` otherwise; same-variant recurrence
is ranked below distinct-variant recurrence because a single recurrent site is as readily a founder
allele or a mapping artifact as a gene-level signal.

**Case-only recurrence null.** For each model with at least two carriers, the probability that a
random individual carries a qualifying genotype was computed under Hardy–Weinberg equilibrium
from the run-oracle frequencies *q*ᵥ of the gene's observed qualifying variants, with a missing or
zero frequency floored at 1 × 10⁻⁶ (`burden.absent_af_floor`): for the dominant model
*p* = 1 − Πᵥ(1 − *q*ᵥ)²; for the biallelic model *p* = (Σᵥ *q*ᵥ)²; for the X-linked model
*p* = 1 − Πᵥ(1 − *q*ᵥ). The recurrence *p*-value is the binomial upper tail P(X ≥ *n*) with
X ~ Binomial(*N*, *p*), where *N* is the number of resolved trios, or for the X-linked model the
number of male probands by Step 0 inference. Benjamini–Hochberg *q*-values are computed within each
model family over the nominated genes, and an exome-wide flag is set at *p* < 2.5 × 10⁻⁶. Because
the null is built only from variants observed in the cohort, it saturates on private variants and
is monotone in gene size; it is reported as a rank, not as a calibrated test.

**Size-normalised rank.** When a per-gene mutational-target table is supplied (the gnomAD v2.1.1
constraint file, unjoined), each gene's target is μ = μ_mis + μ_syn + μ_lof, with a missing μ_lof
imputed as 0.0516 × (μ_mis + μ_syn) and a missing μ_mis or μ_syn yielding no target. A cohort
scaling constant *C* = Σ *n*_carriers / Σ μ was fitted over every gene in the table, including genes
with zero carriers, and each gene received `exp_carriers_mu` = *C* μ, `carrier_excess_ratio` =
*n*_carriers / *C* μ and `p_carrier_excess`, the Poisson upper tail P(X ≥ *n*_carriers). With
`burden.rank_by_mutational_target` (default true) recurrent genes are ordered by this
size-normalised *p* where a target exists (`rank_basis` = mu_normalised), else by the smallest
case-only *p* across models.

**Constraint.** A gene is labelled constrained when any of LOEUF < 0.35, pLI ≥ 0.90,
s_het ≥ 0.10 or pHaplo ≥ 0.86 holds, read from a per-gene table joined during resource
preparation from gnomAD v2.1.1 (LOEUF, pLI), Zeng et al. 2024 GeneBayes s_het and Collins et al.
2022 pHaplo. Constraint is a ranking weight and never an exclusion. The column resolved for each
metric is recorded in the audit.

**Secondary de novo enrichment.** When a Samocha-style per-gene mutation-rate table is supplied,
the expected number of de novo loss-of-function plus missense variants is 2*N*(μ_lof + μ_mis)
(with μ_lof imputed as above when absent, and no test when μ_mis is absent), compared with the
observed count by a Poisson upper tail; *q*-values are computed over the full mutation-table
universe with untested genes at *p* = 1. This arm is uncalibrated against an observed synonymous
rate and is reported only.

Genes are written recurrent-first, distinct-variant before same-variant, then by the rank
*p*-value, then constrained before unconstrained, then by carrier counts and the de novo *p*.

## 13. Prioritization (Step 9)

Step 9 re-ranks the candidate calls without removing any: the number of output rows is asserted
to equal the number of input rows before anything is written. Its input is Step 8's
`igv/variants.tsv` when present (the only table carrying the non-human-read columns), otherwise
Step 5's `candidates.calls.tsv`; a row-count mismatch between the two is fatal because Step 8 is
one-to-one with Step 5.

### 13.1 Gene-level excess over mutational target

For each gene, `n_observed` is the number of distinct (proband, variant) observations, not rows,
because a compound-heterozygous leg appears once per pair and one variant can be emitted under
several modes. The expected count is E = *C* μ, with μ as in Section 12; a gene with no gnomAD
target but a CDS length receives the fallback log μ = −18.4484 + 1.0570 log(CDS length) (an
ordinary-least-squares fit over the gnomAD gene universe, roughly a ±30% band) and is labelled
`E_source = cds_fallback`.

The null is a negative-binomial (NB2) model with mean E and variance E(1 + αE), fitted over the
full gene universe (every gene in the mutational-target table with μ > 0, zero-count genes
included) by an iteratively trimmed procedure: *C* = Σ *n* / Σ μ over the retained genes; α by
golden-section maximum likelihood on [10⁻⁴, 20] (or by moments, `alpha_method`); genes with
P(X ≥ *n*) ≤ 10⁻³ under the current fit are removed from the next fit (never from the output); the
loop repeats to a fixed point (at most 20 iterations). The run halts if more than 5% of eligible
genes are trimmed, which indicates a null mis-specified for the cohort rather than an exome that is
5% artifact. The upper tail `p_nb` = P(X ≥ *n*) is evaluated as the regularised incomplete beta
function I_{E/(E+θ)}(*n*, θ) with θ = 1/α, not as a complement sum, so that tails below machine
epsilon remain resolvable. Benjamini–Hochberg `q_nb` is computed over the full universe, with
CDS-fallback genes corrected among themselves. A Poisson tail fitted with its own trimmed scaling
constant is reported alongside, and a mid-*p* calibration diagnostic (mean mid-*p* and the fraction
below 0.05 and 0.001 over genes with excess ratio < 5) is written to the audit for both nulls.

### 13.2 Six-signal artifact panel

Each gene is scored on six binary signals, summed without weights into `corroboration_count`:
cohort saturation, *n*_observed / *N* ≥ 0.10; segmental duplication, the fraction of the gene
covered by ≥ 98%-identity duplications ≥ 0.10 (from an optional table that must be on the same
coordinate build as the gnomAD v2.1.1 gene intervals, GRCh37; absent, the signal is off); membership
of an artifact-prone gene family by a configurable regular-expression list (olfactory receptors,
KRTAP, MUC, PRAME, NBPF, GOLGA, NPIP, HLA, TAS2R, defensins, LCE, SPRR and histone clusters);
synonymous observed/expected departure |oe_syn − 1| > 0.30; gnomAD's own constraint flag containing
`mis_too_many` or `syn_outlier`; and low cumulative pLoF allele frequency, classic_caf at or below
the 10th percentile of the candidate genes that carry a value, with a gene present in the table but
lacking a value also flagged. The last is a low-information-locus signal, not a frequency signal.

### 13.3 Graded gene down-weight

Tiers are assigned in ascending severity and the highest matching tier wins. T1 (watch, −0.5
points): excess ratio ≥ 3 and `q_nb` < 0.25. T2 (−1.5): ratio ≥ 5 with at least one corroborating
signal. T3 (−3.0): `q_nb` < 0.05 with at least one signal, or ratio ≥ 10 with at least two, or
ratio ≥ 5 with saturation and at least two signals. Ratio-based rules (T2 and the two ratio limbs of
T3) require `n_observed` ≥ 3 (`excess.min_n_for_ratio_rule`); the FDR limb does not. Two ceilings
are applied last: a gene in a supplied phenotype-agnostic established-gene union (for example
ClinGen definitive or strong validity, dosage haploinsufficiency 3, or actionability lists) is
capped at T1 and flagged `established_gene_high_excess`, which protects it from the score penalty
but does not vouch for its calls; a gene whose expectation came from the CDS-length fallback is
capped at T2. A ratio ≥ 10 with no corroborating signal is flagged `unexplained_excess` and
carries no penalty. Two guards halt the run: an established-gene file with fewer than 1,000 genes,
and T2 plus T3 covering more than 20% of distinct observations. Every tier decision is written as a
human-readable `downweight_reason` that reports values and asserts no biology; a high excess ratio
is a statement about calling and mapping quality, not about disease association.

### 13.4 Per-variant tier

Each call is assigned a molecular effect class from its consequence (canonical splice, pLoF,
missense, in-frame indel, otherwise synonymous/UTR/intronic; HIGH impact of any other term is
treated as pLoF) and a tier V0 to V5. V0 (molecularly benign prediction) requires impact LOW or
MODIFIER with both a SpliceAI score < 0.1 and a CADD score < 15 present; a missing score can never
satisfy V0. V5: a stop-gain or frameshift variant whose `nmd_status` is `triggering`, that is, one
the NMD plugin assessed and did not flag as escaping nonsense-mediated decay (`nmd_escape.enabled`,
default true). V4: SpliceAI ≥ 0.5, or a HIGH-impact pLoF or canonical splice variant not promoted to
V5. Canonical splice variants never reach V5, because the plugin judges the variant's own position
rather than the aberrant transcript, and a pLoF whose `nmd_status` is `escaping` or `not_assessed`
stays at V4: the absence of a verdict is never a promotion. pLoF confidence is `UNAVAILABLE`
because LOFTEE is not run. V3:
SpliceAI ≥ 0.2. For missense, a single predictor speaks in a fixed order and is named in
`missense_evidence_source`: REVEL (≥ 0.773 → V4, ≥ 0.644 → V3, ≤ 0.290 → V1, otherwise V2); if
REVEL is absent, AlphaMissense (≥ 0.564 → V3, ≤ 0.34 → V1, otherwise V2); if both are absent, CADD
≥ 25.3 → V3 labelled off-label, otherwise V2. A splice rung and a missense verdict on the same
variant are combined by taking the stronger tier and reporting both reasons. In-frame indels are
V2; LOW/MODIFIER and synonymous/UTR/intronic variants are V1 (a discovery rank); any remaining
MODERATE consequence is V2.

### 13.5 Additive composite and ranking

`priority_points` is an additive score in which every term is emitted as its own column (Table 2).
Rarity is banded on the run oracle's value: < 10⁻⁵ strong (+2), < 10⁻⁴ moderate (+1.5), < 10⁻³
supporting (+1), < 10⁻² permissive (+0.5), between 10⁻² and 0.05 zero, a missing value (absent from
gnomAD) +2, and ≥ 0.05 −8. A flag marks calls whose `MAX_AF` exceeds ten times the grpmax proxy,
indicating that a single excluded population drives the frequency. Gene constraint contributes +1
when the gene is constrained (Section 12), multiplied by a mechanism gate on the variant tier (V0
0, V1 and V2 0.5, V3 to V5 1) and set to zero for recessive and X-linked modes, because pLoF
constraint measures selection against heterozygotes. Recurrence contributes by carrier count, never
by *p*-value: two distinct-variant carriers +1, three or more +2, same-variant recurrence +0.5.
Quality penalties are −2 when the proband genotype fails GQ ≥ 20, depth ≥ 10 or the allele-balance
band appropriate to its zygosity (het [0.25, 0.75]; homozygous or hemizygous ≥ 0.90, with zygosity
derived from the allele-string genotype or the mode), −3 when any screened member's non-human read
fraction is flagged (Section 14.3), 0 when no member was screened, and −0.5 for every
compound-heterozygous leg because the partner leg's parental QC is not carried. A biallelic call
whose allele has a non-zero gnomAD homozygote count is flagged `nhomalt_recessive_conflict` and
charged 0 by default. The clinical term is +4 for a ClinVar pathogenic or likely-pathogenic
assertion and −4 for benign, with the positive term halved when the transferred review status is
below two stars; an unavailable review status leaves full weight. Mode-of-inheritance coherence
against an optional curated per-gene table contributes −1 when discordant and 0 when coherent or
unknown; the penalty is suppressed, and a caveat reported instead, for a heterozygous observation in
a recessive gene (the carrier-risk shape) and for a compound-heterozygous call in a gene with CDS
≥ 10,000 bp. The gene artifact penalty is the tier penalty of Section 13.3. An optional gene-list
prior (+2, scaled by an uncalibrated per-gene weight and by the same mechanism gate; off by default)
enters a second total only. Two hard caps follow: a V0 variant's total is capped at 0 and a BA1
frequency caps it at −4; neither removes the row.

Two rankings are always written: `rank_agnostic` on the total without the gene-list prior and
`rank_prior` with it (identical when no overlay is supplied), plus `rank_delta`. Ties are broken by
molecular points, rarity points, gene artifact penalty, then genomic coordinate, so the order is
reproducible. Every merged source table is also emitted verbatim under a per-source prefix
(`src_mutrate_`, `src_constraint_`, `src_segdup_`, `src_moi_`, `src_burden_`, `src_prior_`) so any
derived number can be re-computed from the output alone. `priority_points` is not an ACMG/AMP score
and is never mapped to a pathogenicity class.

## 14. Review outputs (Steps 7 and 8)

### 14.1 Workbook

Step 7 writes an `.xlsx` workbook with an About sheet (purpose, sheet legend, run summary, the
thresholds and frequency oracle read from the configuration, the VEP header line and the pipeline
commit), followed by the gene consolidation table, the candidate calls, the trio resolution table,
the QC report and the raw audit counts.

### 14.2 IGV review export

Step 8 writes a data directory for the igv.js trio variant-review server: `variants.tsv` with one
row per candidate call (coordinates, gene symbol, consequence, impact, HGVS coding and protein
descriptions, the oracle frequency as the headline `frequency` column, inheritance mode, origin,
pair identifier, the three genotypes, the proband's GQ, depth and AB, every frequency and
predictor column, ClinVar significance and stars, the non-human-read columns, and relative paths
to the alignment and VCF tracks), followed by every remaining column of the calls table verbatim
(the full `flags` field, the de novo tag, the Ensembl gene identifier as `gene_id`, and the
`info_` block), per-trio VCF
tracks, `sample_qc.tsv` from Step 0, a trio list and an empty curation state. When a
sample-to-CRAM map is supplied, mini-CRAMs are sliced for each trio member around the candidate
loci (± 1,000 bp, merged) with `samtools view -C --write-index`, retaining all reads by default,
against the reference the CRAMs were encoded with. Slices are content-keyed on the region set and
reused across runs.

### 14.3 Non-human read fraction (Step 8b, optional)

When a kraken2 database is supplied, the alternate-allele-supporting reads of each screened member
are classified against it with `nonhuman-screen` on a temporary copy of the mini-CRAM filtered
with SAM flag mask 1796 (unmapped, secondary, QC-fail and duplicate reads removed; supplementary
alignments and low mapping quality retained deliberately), at kraken2 confidence 0.05 with the
database memory-mapped. By default the proband is always screened and a parent only where it
carries the alternate allele. Results are joined to `variants.tsv` on the zero-based variant key,
yielding per-member `*_nhf` fractions with their read denominators and a derived `nhf_flag` set
when any screened member has a fraction ≥ 0.5 over ≥ 5 reads. A blank denotes an unscreened
member and is never read as a fraction of zero.

## 15. Audit and provenance

Every step appends (timestamp, step, scope, metric, value) rows to `audit/counts.tsv`, where scope
is `global` or a trio identifier, and the orchestrator assembles `audit/summary.md`. The record
covers trio resolution; the source, post-filter and union site counts; annotation coverage for each
plugin and transfer; Step 3 keep and drop reasons; per-trio candidate genotype and transfer counts;
Step 5's examined/skipped/called/no-row reconciliation and calls by mode; the run-level frequency
oracle and per-variant basis tallies; Step 6's resolved columns, scaling constant and significance
counts; and Step 9's null fit (C, α, θ, iterations, trimmed fraction, dispersion), calibration
diagnostics, tier and signal tallies, and never-drop counts. Provenance columns travel with the
data wherever two identical-looking values are different facts: `rarity_oracle` and
`rarity_basis`; `missense_evidence_source`; `E_source` and `mu_lof_src`; `constraint_source`;
`moi_coherence` and `moi_caveat`; `spliceai_status` (scored versus not covered); `nhf_status`
(clean, flagged, not screened); and `clinvar_review_status` (*N*-star or unavailable).

## 16. Limitations inherent in the implementation

The screen sees single-nucleotide and small indel variation only. Loss-of-function confidence
(LOFTEE) is not computed, and nonsense-mediated-decay escape is taken from the VEP NMD plugin's
positional rules rather than from a transcript model; for canonical splice variants it is a
position proxy and is not used for promotion. The precomputed SpliceAI set does not cover larger indels; a missing score is
recorded as not covered and cannot retain a variant. The case-only recurrence null is a rank, not
a calibrated association test, and the de novo enrichment is uncalibrated. Under the proxy
frequency arm the effective rarity stringency depends on the proband's ancestry, because the
maximum is taken over ancestry groups regardless of the proband's own. The Mendelian-error and
contamination statistics are computed on a capped prefix of the genome and the VCF-only
contamination proxy detects only gross contamination. Genotype QC operates on posterior-derived
qualities, which are not independent of the family prior the inheritance model then uses. The
variant-layer weights are reasoned rather than fitted, and the pipeline has been exercised
end-to-end on generated data (Section 17) but not yet validated against truth sets.

## 17. Availability and testing

The code is available under the MIT licence at
<https://github.com/jlanej/high_priority_rare_variant>; the container image is published to the
GitHub Container Registry on every commit and should be pulled by digest. External resources are
fetched and prepared by `scripts/prepare_resources.sh` from the pinned manifest
`resources/manifest.env` (Table 1). The repository carries a pure-logic test suite (configuration
validation, the annotation getters, genotype QC, the selection ladder, the whole Step 5
inheritance model against a stubbed VCF, Step 6 counting and ranking, and the Step 9 null,
tiering and never-drop invariant) and an end-to-end integration test that generates a small mock
genome and trios, runs the resolve step and Steps 0, 1, 3, 4, 5, 6, 8 and 9 with the real bcftools
toolchain and a mocked VEP call, and asserts the resolution, funnel and calls.

## Table 1. External resources and versions

| Resource | Version or release | Role |
|---|---|---|
| Ensembl VEP and offline cache | release 115, GRCh38 | consequence, IMPACT, symbol, HGVS, exon/CDS geometry, MANE, gnomAD v4.1 point AFs, ClinVar CLIN_SIG |
| Ensembl VEP NMD plugin | release 115 plugin set (no data file) | NMD-escape verdict behind `nmd_status` and the V5 rung |
| gnomAD joint sites (slim) | v4.1 | faf95, FAF group, nhomalt, joint and grpmax AF (default frequency oracle) |
| CADD | v1.7, GRCh38 (whole-genome SNVs; gnomAD-genomes r4.0 indels) | functional score, third selection rung |
| SpliceAI precomputed raw scores | Illumina genome_scores_v1.3, hg38 (SNV and indel) | splice delta scores, second selection rung |
| REVEL | v1.3 (GRCh38-sorted) | missense tier, first predictor |
| AlphaMissense | hg38 release | missense tier, second predictor |
| ClinVar sites VCF | 2026-07-06 (pinned; monthly) | review status → gold stars |
| gnomAD constraint | v2.1.1 lof_metrics.by_gene | LOEUF, pLI; and, unjoined, μ_mis/μ_syn/μ_lof, oe_syn, classic_caf, constraint_flag, CDS length |
| GeneBayes s_het | Zeng et al. 2024 | constraint predicate |
| pHaplo | Collins et al. 2022 | constraint predicate |
| kraken2 database (optional) | user-supplied, with taxonomy dumps | non-human read fraction |
| Reference FASTA | Ensembl release 115 GRCh38 primary assembly (or the calling reference) | normalisation, VEP, CRAM decoding |

## Table 2. Default thresholds and weights

| Parameter | Default | Applied in |
|---|---|---|
| FILTER retained | PASS or `.` | Steps 1, 5 |
| Excluded contigs | chrM, chrMT, M, MT | Step 1 |
| QC sites cap | 200,000 | Step 0 |
| Mendelian-error rate flag | > 0.02 | Step 0 |
| No-call rate flag | > 0.10 | Step 0 |
| chrX het fraction for male | < 0.10 over ≥ 20 sites | Step 0 |
| FREEMIX / VCF-proxy contamination flag | > 0.05 / > 0.02 | Step 0 |
| Frequency oracle | faf95 (gnomAD v4.1 joint) | Steps 3, 5, 6, 9 |
| BA1 (drop, never rescued) | ≥ 0.05 | Steps 3, 9 |
| Recessive / compound-het rarity | < 1 × 10⁻² (high-confidence tag < 1 × 10⁻³) | Steps 3, 5 |
| Dominant / de novo rarity | < 1 × 10⁻⁴ | Step 5 |
| Impact rung | HIGH, MODERATE | Step 3 |
| SpliceAI rung | max Δ ≥ 0.2 | Step 3 |
| CADD rung | PHRED ≥ 25.3 | Step 3 |
| GQ, depth (general) | ≥ 20, ≥ 10 | Steps 0, 5, 9 |
| De novo child depth | ≥ 20 | Step 5 |
| Het allele balance | 0.25–0.75 | Steps 5, 9 |
| Hom-alt allele balance | ≥ 0.90 | Steps 5, 9 |
| Hom-ref allele balance; clean-parent alt reads | ≤ 0.10; ≤ 1 | Step 5 |
| Recurrent gene | ≥ 2 distinct probands | Step 6 |
| Absent-allele floor in the recurrence null | 1 × 10⁻⁶ | Step 6 |
| Exome-wide *p*; FDR | 2.5 × 10⁻⁶; 0.05 | Step 6 |
| Constraint predicate | LOEUF < 0.35 or pLI ≥ 0.90 or s_het ≥ 0.10 or pHaplo ≥ 0.86 | Steps 6, 9 |
| μ_lof imputation factor | 0.0516 | Steps 6, 9 |
| Null model; trim *p*; max trim fraction | NB2 (MLE α); 10⁻³; 0.05 | Step 9 |
| CDS fallback offset | log μ = −18.4484 + 1.0570 log L | Step 9 |
| Saturation; segdup; oe_syn; caf percentile | ≥ 0.10 per trio; ≥ 0.10; \|oe_syn − 1\| > 0.30; 10th | Step 9 |
| Tier T1 / T2 / T3 | ratio ≥ 3 & q < 0.25 / ratio ≥ 5 & ≥ 1 signal / q < 0.05 & ≥ 1, or ratio ≥ 10 & ≥ 2, or ratio ≥ 5 & saturation & ≥ 2 | Step 9 |
| Count floor for ratio rules | *n* ≥ 3 | Step 9 |
| Tier penalties T0–T3 | 0, −0.5, −1.5, −3.0 | Step 9 |
| Established-gene ceiling; CDS ceiling | T1; T2 | Step 9 |
| Min control genes; max down-weight fraction | 1,000; 0.20 | Step 9 |
| SpliceAI strong / supporting / benign | ≥ 0.5 / ≥ 0.2 / < 0.1 | Step 9 |
| SpliceAI event floor | 0.2 (the SpliceAI rung) | Steps 5, 5b |
| Step 5b wide-window rescoring | enabled; 4,999 bp window; 1,000 variants per chunk | Step 5b |
| NMD escape grading (V5) | enabled; the VEP NMD plugin's fixed rules | Steps 2, 5, 9 |
| CADD benign; CADD missense (off-label) | < 15; ≥ 25.3 | Step 9 |
| REVEL moderate / supporting / benign | ≥ 0.773 / ≥ 0.644 / ≤ 0.290 | Step 9 |
| AlphaMissense supporting / benign | ≥ 0.564 / ≤ 0.34 | Step 9 |
| Molecular points V0–V5 | 0, 0.5, 1, 2, 4, 8 | Step 9 |
| Rarity points | strong 2, moderate 1.5, supporting 1, permissive 0.5, absent 2, BA1 −8 | Step 9 |
| Constraint point; gate V0 / V1–V2 / V3–V5 | +1; 0 / 0.5 / 1 (0 for recessive modes) | Step 9 |
| Recurrence points | 2 distinct +1; ≥ 3 distinct +2; same variant +0.5 | Step 9 |
| Quality points | GT fail −2; NHF flagged −3; not screened 0; comp-het partner −0.5; nhomalt conflict 0 | Step 9 |
| Clinical points | P/LP +4 (× 0.5 below 2 stars); benign −4 | Step 9 |
| MOI points | discordant −1; unknown/coherent 0 | Step 9 |
| Gene-list prior | +2 × weight × gate; disabled by default | Step 9 |
| Caps | V0 total ≤ 0; BA1 total ≤ −4 | Step 9 |
| Long-gene CDS for comp-het caveat | ≥ 10,000 bp | Step 9 |
| Mini-CRAM padding | ± 1,000 bp | Step 8 |
| NHF: classify-time mask; confidence; flag fraction; min reads | 1796; 0.05; ≥ 0.5; ≥ 5 | Steps 8b, 9 |

## Table 3. Steps, implementing code and principal outputs

All code paths are relative to the repository root: scripts under `pipeline/`, shared modules under
`src/hprv/`. Outputs are relative to the run's work directory.

| Step | Code (`pipeline/` unless noted) | Principal outputs |
|---|---|---|
| resolve | `resolve_trios.py`; `src/hprv/ped.py` | `trios.resolved.tsv`, `trio_resolution.tsv`, `peds/` |
| 0 | `00_qc.py`; `src/hprv/contamination.py`, `genotype.py` | `qc_report.tsv` |
| 1 | `01_make_cohort_sites.sh` | `cohort.sites.vcf.gz` |
| 2 (2b) | `02_annotate_sites.sh` (`02b_spliceai_backfill.sh`; `src/hprv/spliceai_backfill.py`) | `cohort.sites.annotated.vcf.gz` |
| 3 | `03_select_plausible.py`; `src/hprv/selection.py`, `annotations.py` | `plausible.sites.vcf.gz` |
| 4 | `04_subset_and_annotate_trios.sh` | `trios/<trio>.candidates.annotated.vcf.gz`; `trios.candidates.tsv` |
| 5 | `05_inheritance_screen.py`; `src/hprv/genotype.py`, `splice.py` | `candidates.calls.tsv` |
| 5b | `05b_spliceai_rescore.sh`; `src/hprv/spliceai_rescore.py`, `splice.py` | `candidates.calls.tsv` (`spliceai_wide_*` appended); `spliceai_rescore/scores.tsv` |
| 6 | `06_gene_burden.py` | `genes.ranked.tsv` |
| 7 | `07_report_xlsx.py`; `src/hprv/report.py` | `hprv_summary.xlsx` |
| 8 (8b) | `08_igv_export.sh`; `src/hprv/igv.py` | `igv/variants.tsv`, `igv/crams/`, `igv/vcfs/`, `igv/nhf/` |
| 9 | `09_prioritize.py`; `src/hprv/prioritize.py` | `variants.prioritized.tsv`, `genes.prioritized.tsv`, `igv/variants.prioritized.tsv` |
| all | `run_pipeline.sh`, `lib/common.sh`; `src/hprv/config.py`, `audit.py` | `audit/counts.tsv`, `audit/summary.md` |
