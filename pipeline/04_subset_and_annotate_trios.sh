#!/usr/bin/env bash
# =============================================================================
# 04_subset_and_annotate_trios.sh  —  Pipeline Step 4: per-trio candidate VCFs
#
# For each trio: recover the REAL per-trio genotypes (PP/GQ/DP/AD/hiConfDeNovo)
# at the plausible sites, and attach the annotations computed once in Step 2.
# Keeps per-trio VCFs as the authoritative unit — we never build a genotype-merged
# cohort matrix (that would fabricate hom-ref). See docs/pipeline_design.md (Step 4).
#
# Per trio:  norm -m- -f ref  ->  isec (allele-aware) with plausible sites  ->
#            annotate -c INFO from the (annotated) plausible sites  -> index.
# Emits a manifest of per-trio candidate VCFs for Step 5.
#
# Usage:
#   04_subset_and_annotate_trios.sh --manifest M.tsv --plausible plausible.sites.vcf.gz \
#       --ref GRCh38.fa --outdir OUT [--threads N]
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"

MANIFEST="" PLAUSIBLE="" REF="${HPRV_REF_FASTA:-}" OUTDIR="" THREADS="${HPRV_THREADS:-4}"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)  MANIFEST="$2"; shift 2;;
        --plausible) PLAUSIBLE="$2"; shift 2;;
        --ref)       REF="$2"; shift 2;;
        --outdir)    OUTDIR="$2"; shift 2;;
        --threads)   THREADS="$2"; shift 2;;
        --tmpdir)    HPRV_TMPDIR="$2"; shift 2;;
        -h|--help)   grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done
[[ -n "$MANIFEST" && -n "$PLAUSIBLE" && -n "$REF" && -n "$OUTDIR" ]] \
    || die "need --manifest, --plausible, --ref, --outdir"
[[ -f "$MANIFEST" ]] || die "manifest not found: $MANIFEST"
[[ -f "$PLAUSIBLE" ]] || die "plausible sites not found: $PLAUSIBLE"
[[ -f "$REF" ]] || die "reference not found: $REF"

trio_dir="$OUTDIR/trios"; mkdir -p "$trio_dir"
out_manifest="$OUTDIR/trios.candidates.tsv"

# manifest columns
read -r header < "$MANIFEST"
idcol=0; vcfcol=0; pedcol=0; scol=0; i=0
IFS=$'\t' read -ra cols <<< "$header"
for c in "${cols[@]}"; do i=$((i+1)); case "$c" in trio_id) idcol=$i;; vcf) vcfcol=$i;; ped) pedcol=$i;; samples) scol=$i;; esac; done
[[ $idcol -gt 0 && $vcfcol -gt 0 ]] || die "manifest needs 'trio_id' and 'vcf' columns"

# bind set
binds="$OUTDIR $trio_dir $(abspath_dir "$PLAUSIBLE") $(abspath_dir "$REF")"
rows=()
while IFS= read -r _line || [[ -n "$_line" ]]; do rows+=("$_line"); done < <(tail -n +2 "$MANIFEST")
for row in "${rows[@]}"; do
    [[ -z "$row" || "$row" == \#* ]] && continue
    IFS=$'\t' read -ra f <<< "$row"
    v="${f[$((vcfcol-1))]}"; [[ -n "$v" && -f "$v" ]] && binds+=" $(abspath_dir "$v")"
done
mkdir -p "$HPRV_TMPDIR"; binds+=" $HPRV_TMPDIR"   # norm/cand + shared region BED + sort scratch live here; wrapped tools must see it
HPRV_BIND="$(printf '%s\n' $binds | sort -u | tr '\n' ' ')"; export HPRV_BIND

# ensure plausible sites indexed
index_vcf "$PLAUSIBLE"

# Build a padded, merged BED of plausible loci ONCE. Region-restricting each trio VCF to
# these windows BEFORE `norm` makes per-trio work scale with the candidate set, not the
# genome (norm's ref lookups + left-alignment are the expensive part). The pad absorbs
# left-alignment shifts — a trio indel's raw POS can sit a few bp from its normalized POS
# in the plausible set — so the exact-allele `isec` below stays the real gate and never
# drops a match; a generous window only adds a little norm work, never a spurious call.
PLAUSIBLE_PAD="${HPRV_PLAUSIBLE_PAD:-1000}"
region_bed="$HPRV_TMPDIR/plausible.regions.bed"
bcftools query -f '%CHROM\t%POS\t%REF\n' "$PLAUSIBLE" \
    | awk -v pad="$PLAUSIBLE_PAD" 'BEGIN{OFS="\t"}{s=$2-1-pad; if(s<0)s=0; print $1,s,($2-1+length($3)+pad)}' \
    | sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\t"}{if($1==c&&$2<=e){if($3>e)e=$3}else{if(c!="")print c,s,e;c=$1;s=$2;e=$3}}END{if(c!="")print c,s,e}' \
    > "$region_bed"
REGION_OK=0
if [[ -s "$region_bed" ]]; then REGION_OK=1; else warn "no plausible loci in $PLAUSIBLE; region-restrict disabled"; fi

printf 'trio_id\tcandidates_vcf\tped\n' > "$out_manifest"
log "Step 4: extracting candidate genotypes for ${#rows[@]} trios"

for row in "${rows[@]}"; do
    [[ -z "$row" || "$row" == \#* ]] && continue
    IFS=$'\t' read -ra f <<< "$row"
    trio="${f[$((idcol-1))]}"; vcf="${f[$((vcfcol-1))]}"
    ped=""; [[ $pedcol -gt 0 ]] && ped="${f[$((pedcol-1))]}"
    samples=""; [[ $scol -gt 0 ]] && samples="${f[$((scol-1))]}"
    [[ -n "$trio" && -f "$vcf" ]] || { warn "skipping $trio (missing VCF)"; continue; }

    out="$trio_dir/${trio}.candidates.annotated.vcf.gz"
    if is_done "$out"; then
        log "  [$trio] cached"
        printf '%s\t%s\t%s\n' "$trio" "$out" "$ped" >> "$out_manifest"
        audit 04_subset candidate_genotypes "$(count_variants "$out")" "$trio"
        continue
    fi

    norm="$HPRV_TMPDIR/${trio}.norm.vcf.gz"
    cand="$HPRV_TMPDIR/${trio}.cand.vcf.gz"
    # Region-restrict to plausible loci BEFORE norm when the trio VCF is indexed, so the
    # expensive per-trio norm scales with the candidate set, not the whole genome. `-R`
    # needs an index on the trio VCF; build one if missing, else fall back to the original
    # whole-genome norm (identical result, just slower) so unindexed inputs never regress.
    use_region=0
    if [[ "$REGION_OK" -eq 1 ]]; then
        if [[ -f "$vcf.tbi" || -f "$vcf.csi" ]]; then
            use_region=1
        elif index_vcf "$vcf" 2>/dev/null; then
            use_region=1
        else
            warn "  [$trio] trio VCF is unindexed and could not be indexed; using whole-genome norm"
        fi
    fi
    # Subset to the 3 trio members (dropping extras) so per-trio candidate VCFs carry exactly
    # the trio genotypes, split multiallelics, keep only alleles the trio carries (--min-ac 1,
    # post-split), and STRIP the source INFO (stale internal AC/AN/AF + GATK site fields must
    # not survive into the authoritative candidate VCF; the plausible-site INFO is transferred
    # below). `norm -c w` warns, never silently rewrites, on a REF mismatch.
    #
    # --keep-sum AD is LOAD-BEARING, not a nicety. `norm -m-` subsets the Number=R AD array per
    # allele and DISCARDS the other ALT's reads. For a non-ref/non-ref (1/2) genotype that leaves
    # ref_ad ~ 0 on BOTH legs, so genotype.allele_balance() — alt/(ref+alt) — reads ~1.0, fails the
    # het band, and a REAL trans compound het is dropped with no warning and no audit counter:
    #   child 1/2 AD=0,19,20 -> legs 1/0:0,19 and 0/1:0,20 -> AB=1.000, 1.000 -> het rejected.
    # --keep-sum folds the other ALT's reads back into AD[0], which is the CORRECT per-allele
    # semantic ("reads not supporting THIS alt"): the same legs then read AB 0.487/0.513, the
    # parents' genotypes stay right, and a biallelic het is byte-identical (verified). It
    # hard-errors when AD is absent, so gate on the header — and read the header ONCE into a
    # variable rather than piping into grep, which would SIGPIPE bcftools under `set -o pipefail`.
    nargs=(-m- -f "$REF" -c w)
    vcf_hdr="$(bcftools view -h "$vcf")"
    case "$vcf_hdr" in
        *'##FORMAT=<ID=AD,'*) nargs+=(--keep-sum AD) ;;
        *) warn "  [$trio] no FORMAT/AD in the trio VCF — allele-balance QC is unavailable, and a multiallelic (1/2) het cannot have its AB corrected, so such calls may be dropped" ;;
    esac
    # Strip the stale source INFO but KEEP GATK's de novo tags: they are FORMAT-less INFO fields
    # that Step 5 reads (annotations.F hiconf_denovo/loconf_denovo), and the plausible-site file
    # is site-only so nothing would restore them. A blanket `-x INFO` silently removed both,
    # which made Step 5's `has_hiconf` permanently False and filters.denovo.use_hiconf_tag a
    # no-op. `^` inverts: keep exactly these, drop the rest.
    xargs_info=(-x '^INFO/hiConfDeNovo,INFO/loConfDeNovo')
    if [[ "$use_region" -eq 1 ]]; then
        log "  [$trio] region-restrict + subset-to-trio + norm + intersect + annotate"
        # `-R` emits records in regions-file order (contigs grouped lexically); `bcftools sort`
        # (over the tiny candidate set, negligible) restores header/coordinate order so the
        # output is byte-identical to the whole-genome path, which relies on `norm` output
        # already being coordinate-sorted for its index step to succeed.
        if [[ -n "$samples" ]]; then
            bcftools view -s "$samples" -R "$region_bed" --threads "$THREADS" -Ou "$vcf" \
                | bcftools norm "${nargs[@]}" --threads "$THREADS" -Ou - \
                | bcftools view --min-ac 1 --threads "$THREADS" -Ou - \
                | bcftools annotate "${xargs_info[@]}" --threads "$THREADS" -Ou - \
                | bcftools sort -T "$HPRV_TMPDIR" -Oz -o "$norm" -
        else
            bcftools view -R "$region_bed" --threads "$THREADS" -Ou "$vcf" \
                | bcftools norm "${nargs[@]}" --threads "$THREADS" -Ou - \
                | bcftools view --min-ac 1 --threads "$THREADS" -Ou - \
                | bcftools annotate "${xargs_info[@]}" --threads "$THREADS" -Ou - \
                | bcftools sort -T "$HPRV_TMPDIR" -Oz -o "$norm" -
        fi
    else
        log "  [$trio] subset-to-trio + norm + intersect + annotate"
        if [[ -n "$samples" ]]; then
            bcftools view -s "$samples" --threads "$THREADS" -Ou "$vcf" \
                | bcftools norm "${nargs[@]}" --threads "$THREADS" -Ou - \
                | bcftools view --min-ac 1 --threads "$THREADS" -Ou - \
                | bcftools annotate "${xargs_info[@]}" --threads "$THREADS" -Oz -o "$norm" -
        else
            bcftools norm "${nargs[@]}" --threads "$THREADS" -Ou "$vcf" \
                | bcftools view --min-ac 1 --threads "$THREADS" -Ou - \
                | bcftools annotate "${xargs_info[@]}" --threads "$THREADS" -Oz -o "$norm" -
        fi
    fi
    index_vcf "$norm"
    # allele-aware intersection: trio records that match a plausible site exactly
    bcftools isec -c none -n=2 -w1 --threads "$THREADS" -Oz -o "$cand" "$norm" "$PLAUSIBLE"
    index_vcf "$cand"
    # transfer all INFO annotations from the (annotated) plausible sites
    bcftools annotate -a "$PLAUSIBLE" -c INFO --threads "$THREADS" -Oz -o "$out" "$cand"
    index_vcf "$out"
    require_intact_bgzip "$out"; mark_done "$out"
    rm -f "$norm" "$norm".{tbi,csi} "$cand" "$cand".{tbi,csi} 2>/dev/null || true

    printf '%s\t%s\t%s\n' "$trio" "$out" "$ped" >> "$out_manifest"
    nc="$(count_variants "$out")"
    audit 04_subset candidate_genotypes "$nc" "$trio"
    log "  [$trio] -> $out ($nc candidate genotypes)"
done

log "Step 4 complete. Per-trio candidate manifest: $out_manifest"
