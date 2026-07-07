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
HPRV_BIND="$(printf '%s\n' $binds | sort -u | tr '\n' ' ')"; export HPRV_BIND

# ensure plausible sites indexed
index_vcf "$PLAUSIBLE"

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
    log "  [$trio] subset-to-trio + norm + intersect + annotate"
    # Subset to the 3 trio members (dropping extras) so per-trio candidate VCFs carry
    # exactly the trio genotypes, then normalize to the same representation as sites.
    if [[ -n "$samples" ]]; then
        bcftools view -s "$samples" --threads "$THREADS" -Ou "$vcf" \
            | bcftools norm -m- -f "$REF" -c s --threads "$THREADS" -Oz -o "$norm" -
    else
        bcftools norm -m- -f "$REF" -c s --threads "$THREADS" -Oz -o "$norm" "$vcf"
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
