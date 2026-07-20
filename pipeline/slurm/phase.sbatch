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
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${HPRV_WORK:?HPRV_WORK must be exported (submit_slurm.sh does this via --export)}"
# shellcheck source=/dev/null
source "$HPRV_WORK/slurm_run.env"   # HPRV_CONFIG, container + resource settings (written at submit)

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
        log "submitting scatter array of $n contigs (max ${SCATTER_CONCURRENCY} concurrent) + gather + downstream"
        arr=$(sbatch --parsable \
                --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" ${SLURM_QOS:+--qos="$SLURM_QOS"} \
                --cpus-per-task="$SCATTER_CPUS" --mem="$SCATTER_MEM" --time="$SCATTER_TIME" \
                --array=0-$((n - 1))%"${SCATTER_CONCURRENCY}" \
                --export=ALL,HPRV_WORK="$HPRV_WORK" \
                --job-name=hprv-scatter "$HERE/phase.sbatch" scatter)
        gath=$(sbatch --parsable --dependency=afterok:"$arr" \
                --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" ${SLURM_QOS:+--qos="$SLURM_QOS"} \
                --cpus-per-task="$GATHER_CPUS" --mem="$GATHER_MEM" --time="$GATHER_TIME" \
                --export=ALL,HPRV_WORK="$HPRV_WORK" \
                --job-name=hprv-gather "$HERE/phase.sbatch" gather)
        down=$(sbatch --parsable --dependency=afterok:"$gath" \
                --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" ${SLURM_QOS:+--qos="$SLURM_QOS"} \
                --cpus-per-task="$DOWN_CPUS" --mem="$DOWN_MEM" --time="$DOWN_TIME" \
                --export=ALL,HPRV_WORK="$HPRV_WORK" \
                --job-name=hprv-down "$HERE/phase.sbatch" downstream)
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
        log "Steps 3-8 (select .. igv)"
        rp --from 3 --to 8
        ;;
    *)
        echo "phase.sbatch: unknown phase '$PHASE'" >&2; exit 2;;
esac
