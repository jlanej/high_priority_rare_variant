#!/usr/bin/env bash
# =============================================================================
# phase.sbatch — one node's worth of the hprv pipeline. Invoked by submit_slurm.sh
# (and, for the scatter array + gather + downstream, by the `plan` phase) as:
#     sbatch [resources] phase.sbatch <phase>
# phase in: prep | plan | scatter | gather | downstream
#
# Every phase runs run_pipeline.sh INSIDE the container, so config->env resolution
# and the .done idempotency are identical to a normal run — this file is only the
# scheduler glue. Coherence is enforced by --dependency=afterok between phases (see
# submit_slurm.sh) AND by the gather step's own shard-completeness check.
# =============================================================================
set -euo pipefail
PHASE="${1:?phase.sbatch needs a phase: prep|plan|scatter|gather|downstream}"
# $HERE is RELIABLE only when this file is run directly (a test, or a manual invocation). Under
# sbatch it is NOT: SLURM copies the batch script into the job spool dir and runs it as
# `slurm_script`, so ${BASH_SOURCE[0]} points at /var/spool/slurmd/job<ID>/, not this file on
# shared storage. The `plan` phase therefore re-submits by $SLURM_DIR (below), which
# submit_slurm.sh resolved on the login node and carried through slurm_run.env — never by $HERE.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${HPRV_WORK:?HPRV_WORK must be exported (submit_slurm.sh does this via --export)}"
# shellcheck source=/dev/null
source "$HPRV_WORK/slurm_run.env"   # HPRV_CONFIG, container + resource settings, HPRV_SLURM_DIR
# The real shared-storage dir of this script (see comment above). Falls back to $HERE for a direct
# (non-sbatch) invocation, where $HERE is correct.
SLURM_DIR="${HPRV_SLURM_DIR:-$HERE}"

log() { printf '[phase:%s %s] %s\n' "$PHASE" "$(date +%H:%M:%S)" "$*" >&2; }

# Launch a step inside the container (or natively if HPRV_CONTAINER_BIN is empty).
run() {
    if [[ -n "${HPRV_CONTAINER_BIN:-}" ]]; then
        [[ -n "${APPTAINER_MODULE:-}" ]] && command -v module >/dev/null 2>&1 && module load "$APPTAINER_MODULE"
        "$HPRV_CONTAINER_BIN" exec --cleanenv --bind "$HPRV_BINDS" "$HPRV_SIF" "$@"
    else
        "$@"
    fi
}
rp() { run run_pipeline.sh --config "$HPRV_CONFIG" "$@"; }

MANIFEST="$HPRV_WORK/annotate_shards/manifest.txt"

case "$PHASE" in
    prep)
        # Resolve + Step 0 (QC) + Step 1 (cohort union). Serial for now; the union is the input the
        # contig scatter needs. (Making Steps 0/1/4 trio-arrays is the documented next step.)
        log "resolve + Steps 0-1 (QC + cohort union)"
        rp --from 0 --to 1
        ;;
    plan)
        # Enumerate the union's contigs into a manifest, then submit the scatter ARRAY + the
        # dependent gather + downstream. Job-submits-job is the standard SLURM pattern for an array
        # whose size (contig count) is only known after prep produced the union.
        log "enumerating contigs -> $MANIFEST"
        mkdir -p "$HPRV_WORK/annotate_shards"
        rp --from 2 --to 2 --annotate-emit-manifest "$MANIFEST"
        n="$(wc -l < "$MANIFEST" | tr -d '[:space:]')"
        [[ "$n" -gt 0 ]] || { echo "plan: empty contig manifest" >&2; exit 1; }
        # Re-submit by $SLURM_DIR, NOT $HERE — under sbatch $HERE is the job spool dir (see top).
        # Guard it: a wrong path here surfaces as sbatch's cryptic "Unable to open file", after the
        # whole graph is otherwise wired, so fail early and clearly instead.
        SELF="$SLURM_DIR/phase.sbatch"
        [[ -f "$SELF" ]] || { echo "plan: cannot re-submit — $SELF not found (HPRV_SLURM_DIR wrong? \
is $SLURM_DIR on shared storage the compute nodes can see?)" >&2; exit 1; }
        log "submitting scatter array of $n contigs (max ${SCATTER_CONCURRENCY} concurrent) + gather + downstream"
        arr=$(sbatch --parsable \
                --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" ${SLURM_QOS:+--qos="$SLURM_QOS"} \
                --cpus-per-task="$SCATTER_CPUS" --mem="$SCATTER_MEM" --time="$SCATTER_TIME" \
                --array=0-$((n - 1))%"${SCATTER_CONCURRENCY}" \
                --export=ALL,HPRV_WORK="$HPRV_WORK" \
                --job-name=hprv-scatter "$SELF" scatter)
        gath=$(sbatch --parsable --dependency=afterok:"$arr" \
                --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" ${SLURM_QOS:+--qos="$SLURM_QOS"} \
                --cpus-per-task="$GATHER_CPUS" --mem="$GATHER_MEM" --time="$GATHER_TIME" \
                --export=ALL,HPRV_WORK="$HPRV_WORK" \
                --job-name=hprv-gather "$SELF" gather)
        down=$(sbatch --parsable --dependency=afterok:"$gath" \
                --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" ${SLURM_QOS:+--qos="$SLURM_QOS"} \
                --cpus-per-task="$DOWN_CPUS" --mem="$DOWN_MEM" --time="$DOWN_TIME" \
                --export=ALL,HPRV_WORK="$HPRV_WORK" \
                --job-name=hprv-down "$SELF" downstream)
        log "submitted scatter=$arr gather=$gath downstream=$down"
        printf 'scatter %s\ngather %s\ndownstream %s\n' "$arr" "$gath" "$down" > "$HPRV_WORK/slurm_jobids.txt"
        ;;
    scatter)
        # One array task -> one contig. $SLURM_ARRAY_TASK_ID indexes the manifest (0-based).
        line=$((SLURM_ARRAY_TASK_ID + 1))
        contig="$(sed -n "${line}p" "$MANIFEST")"
        [[ -n "$contig" ]] || { echo "scatter: no contig at manifest line $line" >&2; exit 1; }
        log "VEP contig $contig (array task $SLURM_ARRAY_TASK_ID)"
        rp --from 2 --to 2 --annotate-shard-contig "$contig"
        ;;
    gather)
        # Depends afterok on the whole scatter array. Verifies every shard is complete (dies
        # otherwise — never a partial call set), concats, then finishes Step 2 (split-vep + guards).
        log "gather shards -> annotated union"
        rp --from 2 --to 2 --annotate-gather
        ;;
    downstream)
        # Range is configurable (DOWN_FROM..DOWN_TO, default 3..8). Some sites must stop the
        # distributed graph at Step 7 and run Step 8 (mini-CRAM/igv export) elsewhere — e.g. the
        # source CRAM store is a FUSE mount visible only on login/interactive nodes, not batch
        # compute. Set DOWN_TO=7 in cluster.env and run Step 8 separately (see slurm/README.md).
        df="${DOWN_FROM:-3}"; dt="${DOWN_TO:-8}"
        [[ "$df" =~ ^[0-9]+$ && "$dt" =~ ^[0-9]+$ && "$df" -ge 3 && "$dt" -le 8 && "$df" -le "$dt" ]] \
            || { echo "phase.sbatch: DOWN_FROM..DOWN_TO must satisfy 3<=FROM<=TO<=8 (got '$df'..'$dt')" >&2; exit 2; }
        log "Steps $df-$dt"
        rp --from "$df" --to "$dt"
        ;;
    *)
        echo "phase.sbatch: unknown phase '$PHASE'" >&2; exit 2;;
esac
