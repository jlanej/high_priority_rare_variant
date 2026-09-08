#!/usr/bin/env bash
# =============================================================================
# 05b_spliceai_rescore.sh — Step 5b: wide-window SpliceAI over the CALLED set.
#
# The screen gates on Illumina's PRECOMPUTED SpliceAI scores (-D 50): a cryptic site or a pseudoexon
# partner further than 50 bp from the variant is invisible to them. This step re-scores every
# variant in candidates.calls.tsv LIVE with the stock model at a wide window (default -D 4999, the
# window Walker et al. 2023 used) and appends spliceai_wide_* columns beside the precomputed ones.
# The gate and the Step-9 tier keep reading the precomputed score; the wide columns ADD evidence.
#
# Scores are cached PER VARIANT in <outdir>/scores.tsv (Step 5 rewrites the calls table on every
# run, so a file-level marker would re-score everything every time); only unscored variants are
# planned into chunks, and the merge re-applies the cache to whatever Step 5 just wrote.
#
# Runs in the ISOLATED `spliceai` conda env baked into the image (TensorFlow; see Dockerfile). Exit
# 3 when that env is missing => the caller HALTS (the backfill's contract). Distributed sub-modes
# mirror Step 2's --annotate-* and Step 8b's --nhf-* (one chunk per SLURM array task):
#   --emit-manifest F : plan the pending chunks, copy the manifest to F, exit (no scoring).
#   --chunk N         : score chunk N of the manifest (0-based), exit.
#   --gather          : merge every scored chunk + the cache into the calls table, exit.
# With none of them: plan -> score every chunk serially -> gather (the laptop path).
#
# Usage:
#   05b_spliceai_rescore.sh --calls candidates.calls.tsv --ref GRCh38.fa --outdir W/spliceai_rescore
#       [--distance 4999] [--chunk-size 1000] [--floor 0.2] [--threads N]
#       [--emit-manifest F | --chunk N | --gather]
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"
export PYTHONPATH="${HPRV_HOME:-$(cd "$HERE/.." && pwd)}/src${PYTHONPATH:+:$PYTHONPATH}"

CALLS="" REF="${HPRV_REF_FASTA:-}" OUTDIR="" DISTANCE=4999 CHUNK_SIZE=1000 FLOOR=0.2
THREADS="${HPRV_THREADS:-4}"
SPLICEAI_ENV="${HPRV_SPLICEAI_ENV:-/opt/conda/envs/spliceai}"   # isolated env baked into the image
EMIT="" CHUNK_ID="" GATHER=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --calls) CALLS="$2"; shift 2;;
        --ref) REF="$2"; shift 2;;
        --outdir) OUTDIR="$2"; shift 2;;
        --distance) DISTANCE="$2"; shift 2;;
        --chunk-size) CHUNK_SIZE="$2"; shift 2;;
        --floor) FLOOR="$2"; shift 2;;
        --threads) THREADS="$2"; shift 2;;
        --tmpdir) HPRV_TMPDIR="$2"; shift 2;;
        --emit-manifest) EMIT="$2"; shift 2;;
        --chunk) CHUNK_ID="$2"; shift 2;;
        --gather) GATHER=1; shift;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done
[[ -n "$CALLS" && -f "$CALLS" ]] || die "need --calls <candidates.calls.tsv>"
[[ -n "$OUTDIR" ]] || die "need --outdir"
[[ "$DISTANCE" =~ ^[0-9]+$ && "$CHUNK_SIZE" =~ ^[0-9]+$ && "$CHUNK_SIZE" -ge 1 ]] \
    || die "--distance and --chunk-size must be positive integers"
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -ge 1 ]] || die "--threads must be a positive integer"
_nmode=0; [[ -n "$EMIT" ]] && _nmode=$((_nmode + 1)); [[ -n "$CHUNK_ID" ]] && _nmode=$((_nmode + 1))
[[ "$GATHER" -eq 1 ]] && _nmode=$((_nmode + 1))
[[ "$_nmode" -le 1 ]] || die "--emit-manifest / --chunk / --gather are mutually exclusive"
mkdir -p "$OUTDIR" "$HPRV_TMPDIR"
MANIFEST="$OUTDIR/manifest.txt"
FAI=""; [[ -n "$REF" && -f "$REF.fai" ]] && FAI="$REF.fai"

# Every dir a wrapped tool must see (a host invocation wraps each call in apptainer/docker).
HPRV_BIND="$(printf '%s\n' "$(abspath_dir "$CALLS")" "$OUTDIR" "$HPRV_TMPDIR" ${REF:+"$(abspath_dir "$REF")"} $HPRV_BIND | sort -u | tr '\n' ' ')"
export HPRV_BIND

# The isolated env must exist wherever scoring happens (the image). Exit 3 = UNAVAILABLE: the
# caller halts, as for the Step-2b backfill — a run that silently skipped this would look identical
# to one that scored everything and found nothing distal.
need_env() {
    if ! hprv_run -- test -x "$SPLICEAI_ENV/bin/spliceai"; then
        warn "Step 5b: the isolated 'spliceai' env ($SPLICEAI_ENV) is not present in the runtime. It ships only in the container image — run inside it, point HPRV_SPLICEAI_ENV at the env, or set resources.vep.spliceai_rescore.enabled: false."
        exit 3
    fi
    [[ -n "$REF" && -f "$REF" ]] || die "Step 5b needs --ref <GRCh38 FASTA> (the model reads the sequence around each variant)"
}

plan() {
    local n
    n="$(python3 -m hprv.spliceai_rescore plan --calls "$CALLS" --outdir "$OUTDIR" \
            --distance "$DISTANCE" --chunk-size "$CHUNK_SIZE" ${FAI:+--fai "$FAI"})"
    [[ -f "$MANIFEST" ]] || die "Step 5b: planning wrote no manifest at $MANIFEST"
    log "Step 5b: $n chunk(s) of up to $CHUNK_SIZE variant(s) pending at -D $DISTANCE"
}

# Score ONE chunk. Idempotent on (chunk content, window): a walltime-killed array task re-runs only
# the chunks without a matching .done. TensorFlow's thread pools are sized to the task's CPUs.
score_chunk() {  # $1 = chunk VCF
    local in="$1" out key
    out="${in%.vcf}.scored.vcf"
    key="$(cksum < "$in" | awk '{print $1"-"$2}')-D${DISTANCE}"
    if [[ -s "$out" && "$(cat "$out.done" 2>/dev/null)" == "$key" ]]; then
        log "  [$(basename "$in")] cached"; return 0
    fi
    rm -f "$out" "$out.done"
    log "  [$(basename "$in")] spliceai -D $DISTANCE ($(grep -cv '^#' "$in" || true) variant(s), $THREADS thread(s))"
    # `-p "$SPLICEAI_ENV"` (a prefix), not `-n spliceai`: the availability gate resolves the env
    # through $SPLICEAI_ENV, so a custom prefix must reach the same env here (see 02b).
    OMP_NUM_THREADS="$THREADS" TF_NUM_INTRAOP_THREADS="$THREADS" TF_NUM_INTEROP_THREADS=1 \
        hprv_run -- micromamba run -p "$SPLICEAI_ENV" spliceai -I "$in" -O "$out" -R "$REF" \
            -A grch38 -D "$DISTANCE" \
        || die "Step 5b: spliceai FAILED on $(basename "$in") — check the reference matches the cohort build/contigs and the spliceai env"
    [[ -s "$out" ]] || die "Step 5b: spliceai produced no output for $(basename "$in")"
    printf '%s\n' "$key" > "$out.done"
}

gather() {
    python3 -m hprv.spliceai_rescore gather --calls "$CALLS" --outdir "$OUTDIR" \
        --distance "$DISTANCE" --floor "$FLOOR" \
        || die "Step 5b: gather failed (an unscored chunk, or an unreadable calls table) — see above"
    log "Step 5b complete: spliceai_wide_* columns merged into $CALLS (cache: $OUTDIR/scores.tsv)"
}

if [[ -n "$EMIT" ]]; then
    plan
    [[ "$MANIFEST" -ef "$EMIT" ]] 2>/dev/null || cp -f "$MANIFEST" "$EMIT"
    log "Step 5b: chunk manifest -> $EMIT ($(grep -cve '^[[:space:]]*$' "$EMIT" || true) chunk(s))"
    exit 0
fi
if [[ -n "$CHUNK_ID" ]]; then
    [[ "$CHUNK_ID" =~ ^[0-9]+$ ]] || die "--chunk must be a non-negative integer"
    [[ -f "$MANIFEST" ]] || die "Step 5b: no manifest at $MANIFEST — run the plan first (--emit-manifest)"
    chunk="$(sed -n "$((CHUNK_ID + 1))p" "$MANIFEST")"
    [[ -n "$chunk" && -f "$chunk" ]] || die "Step 5b: no chunk at manifest line $((CHUNK_ID + 1))"
    need_env
    score_chunk "$chunk"
    exit 0
fi
if [[ "$GATHER" -eq 1 ]]; then
    gather
    exit 0
fi

# --- the self-contained serial path ------------------------------------------------------------
need_env
plan
n_pending="$(grep -cve '^[[:space:]]*$' "$MANIFEST" || true)"
if [[ "${n_pending:-0}" -gt 0 ]]; then
    log "Step 5b: scoring $n_pending chunk(s) serially (~1 variant/s/CPU — a SLURM array is minutes; see pipeline/slurm/README.md)"
    while IFS= read -r chunk; do
        [[ -n "$chunk" ]] && score_chunk "$chunk"
    done < "$MANIFEST"
fi
gather
