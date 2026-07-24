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
# The annotation SOURCE is a de-novo-free copy of the plausible sites (the GATK
# hiConfDeNovo/loConfDeNovo tags stripped once, up front): the blanket `-c INFO`
# transfer would otherwise overwrite a candidate's OWN de novo tags (which Step 5
# reads from the trio VCF) if the plausible file ever carried them. Today it does
# not (Step 1 strips all INFO), so this is a structural guard, not a live fix.
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

# Build the annotation SOURCE = plausible sites with the GATK de novo tags stripped, ONCE
# (shared across all trios). The per-trio transfer below is a blanket `bcftools annotate -c INFO`
# (allele-exact: keyed on CHROM+POS+REF+ALT); on any INFO field the candidate ALSO carries, the
# annotation source WINS. The only INFO a candidate carries into that step is its own
# hiConfDeNovo/loConfDeNovo (preserved from the trio VCF; Step 5 reads them). Stripping those two
# tags from the source makes it IMPOSSIBLE for the transfer to clobber the trio's real de novo
# child-list — regardless of what the plausible file grows to carry. Today the plausible file has
# no de novo tags (Step 1 does `annotate -x INFO`), so `-x` here just warns "tag not defined" and
# is a no-op; the guard is defense-in-depth, keeping the "annotations come from VEP, de novo comes
# from the trio" invariant structural rather than incidental. (The expected warnings are silenced;
# a real failure still aborts via the explicit `|| die`.)
PLAUSIBLE_TX="$HPRV_TMPDIR/plausible.tx.vcf.gz"
bcftools annotate -x 'INFO/hiConfDeNovo,INFO/loConfDeNovo' --threads "$THREADS" \
    -Oz -o "$PLAUSIBLE_TX" "$PLAUSIBLE" 2>/dev/null \
    || die "failed to build de-novo-free annotation source from $PLAUSIBLE"
# CRITICAL: PLAUSIBLE_TX is REWRITTEN every run, but index_vcf() is a no-op when ANY .tbi/.csi
# already exists — so a stale index from a prior run (the tmpdir defaults to the PERSISTENT $W/tmp,
# and PLAUSIBLE_TX is never cleaned) would be silently reused. `bcftools annotate` then reads
# offsets that no longer match the rewritten data and SILENTLY drops the vep_*/gnomAD-frequency
# transfer for high-coordinate variants (exit 0, .done still stamped) — corrupting the rarity
# oracle in the authoritative per-trio VCFs. Force a fresh index to match the fresh data.
rm -f "$PLAUSIBLE_TX".tbi "$PLAUSIBLE_TX".csi
index_vcf "$PLAUSIBLE_TX"

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
        #
        # `sort -T` is a PREFIX, NOT a directory — bcftools appends XXXXXX and mkdtemp()s that.
        # So `-T "$HPRV_TMPDIR"` would create the scratch dir as a SIBLING of the tmpdir, i.e. in
        # its PARENT (the work dir), which is how this failed in the field:
        #   mkdtemp(<work>/tmpSNIkLK) failed: Read-only file system
        # even though everything else written INTO $HPRV_TMPDIR (the region BED, the .norm.vcf.gz)
        # worked fine. Keep the prefix INSIDE the tmpdir — "$HPRV_TMPDIR/<trio>.sort" makes
        # bcftools create "$HPRV_TMPDIR/<trio>.sortXXXXXX" — which also keeps the scratch
        # per-trio, so this stays safe if the loop is ever parallelized. (Step 1 already does
        # this via "$HPRV_TMPDIR/sort", which is why it never hit the bug.)
        if [[ -n "$samples" ]]; then
            bcftools view -s "$samples" -R "$region_bed" --threads "$THREADS" -Ou "$vcf" \
                | bcftools norm "${nargs[@]}" --threads "$THREADS" -Ou - \
                | bcftools view --min-ac 1 --threads "$THREADS" -Ou - \
                | bcftools annotate "${xargs_info[@]}" --threads "$THREADS" -Ou - \
                | bcftools sort -T "$HPRV_TMPDIR/${trio}.sort" -Oz -o "$norm" -
        else
            bcftools view -R "$region_bed" --threads "$THREADS" -Ou "$vcf" \
                | bcftools norm "${nargs[@]}" --threads "$THREADS" -Ou - \
                | bcftools view --min-ac 1 --threads "$THREADS" -Ou - \
                | bcftools annotate "${xargs_info[@]}" --threads "$THREADS" -Ou - \
                | bcftools sort -T "$HPRV_TMPDIR/${trio}.sort" -Oz -o "$norm" -
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
    # transfer all INFO annotations from the (annotated) plausible sites. The source is the
    # de-novo-free copy built above, so this blanket `-c INFO` cannot overwrite the candidate's
    # own hiConfDeNovo/loConfDeNovo (see the PLAUSIBLE_TX note). Allele-exact: keyed CHROM+POS+REF+ALT.
    bcftools annotate -a "$PLAUSIBLE_TX" -c INFO --threads "$THREADS" -Oz -o "$out" "$cand"
    index_vcf "$out"
    require_intact_bgzip "$out"; mark_done "$out"
    rm -f "$norm" "$norm".{tbi,csi} "$cand" "$cand".{tbi,csi} 2>/dev/null || true

    printf '%s\t%s\t%s\n' "$trio" "$out" "$ped" >> "$out_manifest"
    nc="$(count_variants "$out")"
    audit 04_subset candidate_genotypes "$nc" "$trio"
    log "  [$trio] -> $out ($nc candidate genotypes)"
done

log "Step 4 complete. Per-trio candidate manifest: $out_manifest"
