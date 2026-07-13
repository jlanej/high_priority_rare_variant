#!/usr/bin/env bash
# =============================================================================
# 08_igv_export.sh  —  Pipeline Step 8: export for the igv.js trio review server
#
# Produces, under WORK/igv/ (the server's --data-dir):
#   variants.tsv     one row per candidate call (igv.js variant-review schema:
#                    chrom/pos/ref/alt + inheritance mode + genotypes + annotations
#                    + per-member track columns); any extra column is filterable.
#   crams/<trio>/<sample>.cram(+.crai)   mini-CRAMs sliced around candidate loci
#                    for child/mother/father (only if a sample->CRAM map is given).
#   vcfs/<trio>.vcf.gz(+.tbi)            per-trio candidate VCF track.
#   sample_qc.tsv    trio_id/role/sample_id + QC metrics (optional --sample-qc input).
#   trios.tsv        #kid mom dad.
#   curation.json    empty review state ({}).
#
# Serve with the jlanej/igv.js variant-review server:
#   node server.js --variants WORK/igv/variants.tsv --data-dir WORK/igv --genome hg38
#   (container: ghcr.io/jlanej/igv-variant-review)
#
# Usage:
#   08_igv_export.sh --work WORKDIR --ref GRCh38.fa [--cram-map map.tsv] [--padding 1000]
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"
export PYTHONPATH="${PYTHONPATH:-}:${HPRV_HOME:-$(cd "$HERE/.." && pwd)}/src"

# --ref is the CRAM ENCODING reference (needed to decode reference-compressed CRAMs),
# which may differ from the variant-calling reference used elsewhere; run_pipeline.sh
# passes resources.cram_ref (falling back to reference.fasta). --jobs bounds how many
# CRAM slices run concurrently: mini-slices are small index-based reads, but on a flaky
# FUSE/SBFS mount unbounded concurrency destabilizes it, so this is capped deliberately.
WORK="" REF="${HPRV_CRAM_REF:-${HPRV_REF_FASTA:-}}" CRAM_MAP="${HPRV_CRAM_MAP:-}" PAD=1000 GENOME=hg38 JOBS=1
while [[ $# -gt 0 ]]; do
    case "$1" in
        --work) WORK="$2"; shift 2;;
        --ref) REF="$2"; shift 2;;
        --cram-map) CRAM_MAP="$2"; shift 2;;
        --padding) PAD="$2"; shift 2;;
        --genome) GENOME="$2"; shift 2;;
        --jobs) JOBS="$2"; shift 2;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done
[[ "$JOBS" =~ ^[0-9]+$ && "$JOBS" -ge 1 ]] || die "--jobs must be a positive integer"
[[ -n "$WORK" ]] || die "need --work"
calls="$WORK/candidates.calls.tsv"; resolved="$WORK/trios.resolved.tsv"
cand="$WORK/trios.candidates.tsv"
[[ -f "$calls" && -f "$resolved" ]] || die "missing $calls or $resolved (run steps 4-6 first)"

is_set() { [[ -n "${1:-}" && "$1" != *'${'* ]]; }
DATA="$WORK/igv"; mkdir -p "$DATA/crams" "$DATA/vcfs"

# --- source-CRAM map (sample -> path); optional. Looked up on demand (no assoc
#     arrays, for bash 3.2 portability), matching the bespoke get_cram approach. ---
HAVE_MAP=0
binds="$DATA"; is_set "$REF" && [[ -f "$REF" ]] && binds+=" $(abspath_dir "$REF")"
cram_for() { [[ -f "${CRAM_MAP:-/nonexistent}" ]] || return 0; awk -F'\t' -v s="$1" '$1==s{print $2; exit}' "$CRAM_MAP"; }
if is_set "$CRAM_MAP" && [[ -f "$CRAM_MAP" ]]; then
    HAVE_MAP=1
    while IFS=$'\t' read -r s pth _; do
        [[ -z "$s" || "$s" == \#* ]] && continue
        [[ -f "$pth" ]] && binds+=" $(abspath_dir "$pth")"
    done < "$CRAM_MAP"
    log "Step 8: using CRAM map $CRAM_MAP"
else
    warn "no --cram-map; variants.tsv will be emitted without alignment tracks"
fi
HPRV_BIND="$(printf '%s\n' $binds | sort -u | tr '\n' ' ')"; export HPRV_BIND

# --- pre-group per-trio inputs in ONE pass each (was O(trios x calls) re-scanning) ---
mkdir -p "$HPRV_TMPDIR"
# 1) candidate loci -> one padded BED per trio. Tag each call with its trio, sort so a
#    trio's rows are contiguous, then split into per-trio BEDs keeping a single file
#    handle open at a time (scales to many trios without hitting awk's open-file limit).
beddir="$HPRV_TMPDIR/callbeds"; mkdir -p "$beddir"
awk -F'\t' -v pad="$PAD" '
    NR==1{for(i=1;i<=NF;i++){if($i=="trio_id")ti=i; if($i=="chrom")ci=i; if($i=="pos")pi=i}; next}
    ti&&ci&&pi{s=$pi-1-pad; if(s<0)s=0; print $ti"\t"$ci"\t"s"\t"($pi+pad)}' "$calls" \
    | sort -k1,1 -k2,2 -k3,3n \
    | awk -v dir="$beddir" '$1==""{next} {if($1!=prev){if(prev!="")close(pf); pf=dir"/"$1".bed"; prev=$1} print $2"\t"$3"\t"$4 > pf}'
# 2) trio -> candidate VCF path (trios.candidates.tsv has one row per trio). Write the
#    first path per trio to a tiny lookup file so the loop below is O(1) per trio.
vcfmapdir="$HPRV_TMPDIR/candvcf"; mkdir -p "$vcfmapdir"
[[ -f "$cand" ]] && awk -F'\t' -v dir="$vcfmapdir" 'NR>1 && $1!="" && !seen[$1]++{f=dir"/"$1; print $2 > f; close(f)}' "$cand"

# --- Pass 1 (serial, cheap): per trio, build the merged BED + copy the VCF track,
#     and emit one CRAM-slice TASK per member into a task list. No FUSE reads here. ---
tasks="$HPRV_TMPDIR/extract_tasks.tsv"; : > "$tasks"
while IFS=$'\t' read -r trio vcf ped samples; do
    [[ "$trio" == "trio_id" || -z "$trio" ]] && continue
    IFS=',' read -r kid dad mom <<< "$samples"

    bed="$beddir/${trio}.bed"
    if [[ -s "$bed" ]]; then
        merged="$HPRV_TMPDIR/${trio}.merged.bed"
        hprv_run -- bedtools merge -i "$bed" > "$merged" 2>/dev/null || cp "$bed" "$merged"
        if is_set "$REF" && [[ -f "$REF" && "$HAVE_MAP" -eq 1 ]]; then
            odir="$DATA/crams/$trio"; mkdir -p "$odir"
            for pair in "child:$kid" "mother:$mom" "father:$dad"; do
                role="${pair%%:*}"; sample="${pair#*:}"
                [[ -n "$sample" ]] || continue
                src="$(cram_for "$sample")"
                [[ -n "$src" && -f "$src" ]] || { warn "  [$trio] no CRAM for $role $sample"; continue; }
                # trio<TAB>role<TAB>sample<TAB>src<TAB>merged<TAB>ocram
                printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$trio" "$role" "$sample" "$src" "$merged" "$odir/${sample}.cram" >> "$tasks"
            done
        fi
    fi

    # per-trio VCF track (copy + index the candidate VCF under the data-dir)
    cvcf=""
    if [[ -f "$vcfmapdir/$trio" ]]; then IFS= read -r cvcf < "$vcfmapdir/$trio" || true; fi
    if [[ -n "$cvcf" && -f "$cvcf" ]]; then
        cp -f "$cvcf" "$DATA/vcfs/${trio}.vcf.gz"
        hprv_run --bind "$DATA" -- bcftools index -t -f "$DATA/vcfs/${trio}.vcf.gz" 2>/dev/null || true
    fi
done < <(tail -n +2 "$resolved")

# --- Pass 2 (serial): slice each candidate region set from its source CRAM, one CRAM
#     open at a time. Reads are index-based (small); serial access is the gentlest pattern
#     for a flaky FUSE/SBFS mount (the multi-core part is samtools' own -@ $JOBS per slice).
#     Each slice validates (quickcheck) and retries once on a transient read failure; a
#     persistently-bad source CRAM warns and is skipped (never aborts the whole export). ---
# Slice one candidate-region set from one source CRAM. Runs in a SUBSHELL (see the
# loop) so its `set +e` is local: we must handle samtools exit codes ourselves and
# retry/skip a transiently-bad read rather than let `set -e` abort the whole export
# on the first failure. Slices AND indexes in one pass (--write-index) because the
# igv.js server needs the .crai. Returns 0 if a valid slice landed, 1 to skip.
extract_one() {
    set +e
    local trio="$1" role="$2" sample="$3" src="$4" merged="$5" ocram="$6" attempt
    for attempt in 1 2; do
        if hprv_run -- samtools view -C -@ "$JOBS" -T "$REF" --regions-file "$merged" \
                --write-index -o "$ocram" "$src" 2>/dev/null \
           && [[ -s "$ocram" && -f "$ocram.crai" ]] \
           && hprv_run -- samtools quickcheck "$ocram" 2>/dev/null; then
            return 0
        fi
        rm -f "$ocram" "$ocram.crai"
        [[ "$attempt" -eq 1 ]] && { warn "  [$trio] slice transient-failed for $role $sample; retrying"; sleep 2; }
    done
    warn "  [$trio] CRAM slice failed for $role $sample ($src); skipping this alignment track"
    return 1
}

# Extract SERIALLY — one source CRAM open at a time. On a flaky FUSE/SBFS mount this
# is the gentlest access pattern (the user's stated requirement: "keep the mount well
# behaved"); the multi-core part is samtools' own --threads (-@ $JOBS) per slice.
# Run in a single interactive job where the CRAM mount is available.
# NB: `grep -c` prints 0 AND exits 1 on an all-blank/empty file, so `|| echo 0` would
# append a second "0" -> a "0\n0" that breaks the `-gt` test. Keep grep's own count.
n_tasks=$(grep -cve '^[[:space:]]*$' "$tasks" 2>/dev/null || true); n_tasks=${n_tasks:-0}
n_extracted=0
if [[ "$n_tasks" -gt 0 ]]; then
    log "Step 8: slicing $n_tasks mini-CRAM(s) serially (FUSE-safe; $JOBS thread(s)/slice)"
    while IFS=$'\t' read -r trio role sample src merged ocram; do
        [[ -z "$trio" ]] && continue
        if ( extract_one "$trio" "$role" "$sample" "$src" "$merged" "$ocram" ); then
            n_extracted=$((n_extracted + 1))
        fi
    done < "$tasks"
fi

# --- assemble variants.tsv + sample_qc.tsv + trios.tsv + curation.json ---
python3 - "$calls" "$resolved" "$DATA" "$WORK/qc_report.tsv" <<'PY'
import sys
from hprv import igv
calls, manifest, data, qc = sys.argv[1:5]
n = igv.build_variants_tsv(calls, manifest, data, f"{data}/variants.tsv")
m = igv.write_sample_qc(qc, manifest, f"{data}/sample_qc.tsv")
sys.stderr.write(f"  variants.tsv rows: {n}; sample_qc rows: {m}\n")
PY

# trios.tsv (#kid mom dad) + empty curation
{ printf '#kid\tmom\tdad\n'
  awk -F'\t' 'NR>1{split($4,s,","); if(s[1]!="") print s[1]"\t"s[3]"\t"s[2]}' "$resolved"; } > "$DATA/trios.tsv"
[[ -f "$DATA/curation.json" ]] || echo '{}' > "$DATA/curation.json"
printf '{"genome": "%s"}\n' "$GENOME" > "$DATA/config.json"   # server --genome hint

audit 08_igv variants "$(($(grep -cve '^[[:space:]]*$' "$DATA/variants.tsv") - 1))"
audit 08_igv minicrams "$n_extracted"
log "Step 8 complete: igv review export -> $DATA (variants.tsv, crams/, vcfs/, trios.tsv, curation.json)"
