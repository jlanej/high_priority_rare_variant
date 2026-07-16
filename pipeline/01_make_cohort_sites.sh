#!/usr/bin/env bash
# =============================================================================
# 01_make_cohort_sites.sh  —  Pipeline Step 1: cohort site-only union VCF
#
# Build a normalized, SITE-ONLY union of variant loci across all per-trio VCFs.
# See docs/cohort_construction.md and docs/pipeline_design.md (Step 1).
#
# WHY site-only union, not `bcftools merge`:
#   These trios were NOT jointly genotyped. Merging them fabricates hom-ref
#   genotypes (absent != hom-ref), so any internal AC/AN is fiction. We therefore
#   drop genotypes and build a union of *loci only*. Population frequency comes
#   from external gnomAD downstream — never from this file.
#
# Per trio:  keep PASS -> norm (split multiallelic + left-align) -> drop GT ->
#            strip per-trio INFO (incomparable across trios) -> index
# Union:     concat -a -D (dedup identical) -> sort -> norm -d exact (collapse
#            residual duplicate representations) -> index
#
# Usage:
#   01_make_cohort_sites.sh --manifest M.tsv --ref GRCh38.fa --out cohort.sites.vcf.gz \
#       [--tmpdir DIR] [--threads N] [--filter 'PASS,.']
#
# Manifest: TSV with a header line; columns include at least `trio_id` and `vcf`.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"

MANIFEST="" REF="" OUT="" THREADS="${HPRV_THREADS:-4}" FILTER="PASS,."
while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest) MANIFEST="$2"; shift 2;;
        --ref)      REF="$2"; shift 2;;
        --out)      OUT="$2"; shift 2;;
        --tmpdir)   HPRV_TMPDIR="$2"; shift 2;;
        --threads)  THREADS="$2"; shift 2;;
        --filter)   FILTER="$2"; shift 2;;
        -h|--help)  grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done
[[ -n "$MANIFEST" && -n "$REF" && -n "$OUT" ]] || die "need --manifest, --ref, and --out"
[[ -f "$MANIFEST" ]] || die "manifest not found: $MANIFEST"
[[ -f "$REF" ]] || die "reference FASTA not found: $REF"

if is_done "$OUT"; then log "Step 1 already complete: $OUT (skipping)"; exit 0; fi

workdir="$(abspath_dir "$OUT")/sites_work"
mkdir -p "$workdir"

# --- resolve trio_id / vcf / samples columns from the manifest header ---
read -r header < "$MANIFEST"
idcol=0; vcfcol=0; scol=0; i=0
IFS=$'\t' read -ra cols <<< "$header"
for c in "${cols[@]}"; do
    i=$((i+1))
    case "$c" in trio_id) idcol=$i;; vcf) vcfcol=$i;; samples) scol=$i;; esac
done
[[ $idcol -gt 0 && $vcfcol -gt 0 ]] || die "manifest must have tab-separated 'trio_id' and 'vcf' header columns"

# --- collect inputs and build the bind set (dirs every tool call must see) ---
# HPRV_TMPDIR must be bound too: the union's `bcftools sort -T` scratch lives under it.
declare -a site_files=()
binds="$(abspath_dir "$OUT") $workdir $(abspath_dir "$REF") ${HPRV_TMPDIR:-}"
rows=()
while IFS= read -r _line || [[ -n "$_line" ]]; do rows+=("$_line"); done < <(tail -n +2 "$MANIFEST")
[[ ${#rows[@]} -gt 0 ]] || die "manifest has no data rows"

for row in "${rows[@]}"; do
    [[ -z "$row" || "$row" == \#* ]] && continue
    IFS=$'\t' read -ra f <<< "$row"
    trio="${f[$((idcol-1))]}"; vcf="${f[$((vcfcol-1))]}"
    [[ -n "$trio" && -n "$vcf" ]] || { warn "skipping malformed row: $row"; continue; }
    [[ -f "$vcf" ]] || die "trio VCF not found for $trio: $vcf"
    binds+=" $(abspath_dir "$vcf")"
done
# Deduplicate bind dirs and export for common.sh so all tool calls can see them.
HPRV_BIND="$(printf '%s\n' $binds | sort -u | tr '\n' ' ')"
export HPRV_BIND

# chrM is OUT OF SCOPE (see CLAUDE.md: mtDNA heteroplasmy has a dedicated pipeline). Excluded
# HERE, at the union, rather than with a guard in Step 5, so every downstream artifact agrees
# chrM was never analyzed — the annotation, the audit funnel, and the Step-7 workbook alike.
#
# It MUST be excluded, not merely left un-modelled: every inheritance mode here is DIPLOID, and
# genotype.py's is_x_nonpar/is_y_nonpar know only X and Y — so a chrM record routes through
# hom_recessive / dominant / compound_het as though it were an autosome. Against rCRS the
# near-fixed haplogroup variants (m.8860A>G in MT-ATP6, m.15326A>G in MT-CYB) are HOM_ALT in
# every member INCLUDING the father, so every trio fires hom_recessive at the same MT genes. The
# amplifier: the VEP cache carries no gnomAD mito AF, so frequency() returns None, rare() passes
# unconditionally, and Step 6 floors q to absent_af_floor -> p ~ 1e-12 -> those MT genes land in
# the recurrent, exome-wide-significant tier, ranked ABOVE every genuine non-recurrent nuclear
# candidate. Both chr-prefixed and Ensembl-style names are listed; naming an absent contig is a
# no-op in bcftools (verified), so this is safe on either convention.
EXCLUDE_CONTIGS="${HPRV_EXCLUDE_CONTIGS:-chrM,chrMT,M,MT}"

log "Step 1: building cohort site-only union from ${#rows[@]} trios"
log "  reference: $REF"
log "  output:    $OUT"
log "  excluding contigs (out of scope): $EXCLUDE_CONTIGS"

# --- per-trio: PASS -> norm -> site-only -> strip INFO -> index ---
for row in "${rows[@]}"; do
    [[ -z "$row" || "$row" == \#* ]] && continue
    IFS=$'\t' read -ra f <<< "$row"
    trio="${f[$((idcol-1))]}"; vcf="${f[$((vcfcol-1))]}"
    samples=""; [[ $scol -gt 0 ]] && samples="${f[$((scol-1))]}"
    [[ -n "$trio" && -n "$vcf" ]] || continue
    site="$workdir/${trio}.sites.norm.vcf.gz"

    if is_done "$site"; then
        log "  [$trio] cached"
        [[ -f "$site.tbi" || -f "$site.csi" ]] || index_vcf "$site"  # concat -a needs the index
        site_files+=("$site")
        audit 01_cohort_sites input_sites "$(count_variants "$site")" "$trio"
        continue
    fi
    # Subset to the trio's 3 members (dropping any extra members in a multi-sample VCF),
    # split multiallelics, then keep only alleles the trio actually carries (--min-ac 1
    # applied AFTER the split, so AC=0 alt alleles do not leak into the union). Then drop
    # genotypes and strip per-trio INFO. `norm -c w` warns (never silently rewrites) on a
    # REF mismatch — a build/contig problem should be visible, not masked.
    log "  [$trio] subsetting to trio + normalizing + dropping genotypes"
    if [[ -n "$samples" ]]; then
        bcftools view -s "$samples" -f "$FILTER" -t "^$EXCLUDE_CONTIGS" --threads "$THREADS" -Ou "$vcf" \
            | bcftools norm -m- -f "$REF" -c w -Ou - \
            | bcftools view --min-ac 1 -Ou - \
            | bcftools view -G -Ou - \
            | bcftools annotate -x INFO --threads "$THREADS" -Oz -o "$site" -
    else
        bcftools view -f "$FILTER" -t "^$EXCLUDE_CONTIGS" --threads "$THREADS" -Ou "$vcf" \
            | bcftools norm -m- -f "$REF" -c w -Ou - \
            | bcftools view --min-ac 1 -Ou - \
            | bcftools view -G -Ou - \
            | bcftools annotate -x INFO --threads "$THREADS" -Oz -o "$site" -
    fi
    index_vcf "$site"
    require_intact_bgzip "$site"
    mark_done "$site"
    site_files+=("$site")
    audit 01_cohort_sites input_sites "$(count_variants "$site")" "$trio"
done
[[ ${#site_files[@]} -gt 0 ]] || die "no per-trio site files were produced"

# --- union: concat (dedup) -> sort -> norm -d exact -> index ---
# Pass the inputs via a file list, not argv — thousands of trios would blow ARG_MAX / the
# open-file-descriptor limit when spread across the command line.
log "Step 1: unioning ${#site_files[@]} per-trio site files"
tmp_sort="$HPRV_TMPDIR/sort"; mkdir -p "$tmp_sort"
filelist="$workdir/union_filelist.txt"
printf '%s\n' "${site_files[@]}" > "$filelist"
bcftools concat -a -D -f "$filelist" --threads "$THREADS" -Ou \
    | bcftools sort -T "$tmp_sort" -Ou - \
    | bcftools norm -d exact -f "$REF" --threads "$THREADS" -Oz -o "$OUT" -
index_vcf "$OUT"
require_intact_bgzip "$OUT"
mark_done "$OUT"

n="$(count_variants "$OUT")"
audit 01_cohort_sites trios "${#site_files[@]}"
audit 01_cohort_sites union_sites "$n"
# Record the exclusion so it is an auditable decision rather than a silent disappearance: a
# reader of audit/counts.tsv can see chrM was dropped on purpose and how much was dropped.
# Counted on the union's OWN contigs (a $OUT that still contains them means the filter failed).
n_excluded=0
for _c in ${EXCLUDE_CONTIGS//,/ }; do
    _n="$(bcftools view -H -t "$_c" "$OUT" 2>/dev/null | wc -l | tr -d '[:space:]')"
    n_excluded=$((n_excluded + ${_n:-0}))
done
audit 01_cohort_sites excluded_contigs_remaining "$n_excluded"
[[ "$n_excluded" -eq 0 ]] || die "the cohort union still contains $n_excluded records on the \
out-of-scope contigs ($EXCLUDE_CONTIGS) — the exclusion did not take. Every inheritance mode \
downstream is diploid and would treat them as autosomal; chrM in particular would flood the \
recurrent tier of genes.ranked.tsv. Check the contig naming in your VCFs."
log "Step 1 complete: $OUT ($n unique sites; 0 on the excluded contigs)"
