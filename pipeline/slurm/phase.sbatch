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
PHASE="${1:?phase.sbatch needs a phase: prep|plan|scatter|gather|downstream|nhf-plan|nhf-scatter|nhf-gather}"
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
    nhf-plan)
        # Step 8b (NHF) as a PER-TRIO array. Requires a COMPLETED Step 8 — the array reads the
        # mini-CRAMs and per-trio VCFs that Step 8's Pass 1/2 produce, and deliberately does not
        # rebuild them (concurrent rebuilds would corrupt shared state). Unlike the Step-2 scatter,
        # 8b touches NO source-CRAM mount, so it is safe on batch compute even where slicing is not.
        MANIFEST="$HPRV_WORK/igv/nhf_trios.txt"
        mkdir -p "$(dirname "$MANIFEST")"
        rp --from 8 --to 8 --nhf-emit-manifest "$MANIFEST"
        # The manifest must EXIST. Step 8 now dies when a sub-task has no kraken2 DB, but an
        # absent file here must never read as "nothing to do" — that was a green no-op.
        [[ -f "$MANIFEST" ]] || { echo "nhf-plan: Step 8 wrote no NHF manifest at $MANIFEST (no kraken2 DB / no mini-CRAMs?)" >&2; exit 1; }
        n=$(grep -cve '^[[:space:]]*$' "$MANIFEST" 2>/dev/null || true); n=${n:-0}
        if [[ "$n" -eq 0 ]]; then
            log "no trios need NHF work — nothing to scatter; running gather to refresh variants.tsv"
            rp --from 8 --to 8 --nhf-gather
            exit 0
        fi
        SELF="$SLURM_DIR/phase.sbatch"
        [[ -f "$SELF" ]] || { echo "phase.sbatch: cannot find self at '$SELF'" >&2; exit 2; }
        log "submitting NHF array of $n trio(s) (max ${NHF_CONCURRENCY:-8} concurrent) + gather"
        arr=$(sbatch --parsable \
                --cpus-per-task="${NHF_CPUS:-8}" --mem="${NHF_MEM:-32G}" --time="${NHF_TIME:-04:00:00}" \
                --array=0-$((n - 1))%"${NHF_CONCURRENCY:-8}" \
                ${SLURM_ACCOUNT:+--account="$SLURM_ACCOUNT"} ${SLURM_PARTITION:+--partition="$SLURM_PARTITION"} \
                "$SELF" nhf-scatter)
        sbatch --parsable --dependency=afterany:"$arr" \
                --cpus-per-task="${NHF_GATHER_CPUS:-2}" --mem="${NHF_GATHER_MEM:-8G}" --time="${NHF_GATHER_TIME:-01:00:00}" \
                ${SLURM_ACCOUNT:+--account="$SLURM_ACCOUNT"} ${SLURM_PARTITION:+--partition="$SLURM_PARTITION"} \
                "$SELF" nhf-gather >/dev/null
        # afterANY, not afterok: a failed trio is a legitimate outcome here (blank NHF columns),
        # so the join must still run and report what is missing rather than the DAG stalling.
        ;;
    nhf-scatter)
        # One array task -> one trio. $SLURM_ARRAY_TASK_ID indexes the manifest (0-based).
        MANIFEST="$HPRV_WORK/igv/nhf_trios.txt"
        [[ -f "$MANIFEST" ]] || { echo "phase.sbatch: missing $MANIFEST (run nhf-plan first)" >&2; exit 2; }
        trio=$(awk -v i="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==i{print; exit}' "$MANIFEST")
        [[ -n "$trio" ]] || { echo "phase.sbatch: no trio at manifest line $((SLURM_ARRAY_TASK_ID + 1))" >&2; exit 2; }
        # Task-private TMPDIR: 8b writes a carriers list + per-member .err under $HPRV_TMPDIR, and
        # concurrent tasks sharing one path would read each other's half-written files.
        export HPRV_TMPDIR="${HPRV_TMPDIR:-/tmp}/nhf.$trio.$$"
        mkdir -p "$HPRV_TMPDIR"
        # `run()` execs the container with --cleanenv, which strips HPRV_TMPDIR; the *ENV_ prefixed
        # forms survive it, so the task-private scratch actually reaches the step.
        export APPTAINERENV_HPRV_TMPDIR="$HPRV_TMPDIR" SINGULARITYENV_HPRV_TMPDIR="$HPRV_TMPDIR"
        log "NHF trio $trio (array task $SLURM_ARRAY_TASK_ID)"
        rp --from 8 --to 8 --nhf-trio "$trio"
        ;;
    nhf-gather)
        # Join only: fold whatever NHF results exist into variants.tsv. Never dies on an
        # unscreened member — that is a legitimate state (see 08_igv_export.sh).
        log "NHF gather -> variants.tsv"
        rp --from 8 --to 8 --nhf-gather
        ;;
    downstream)
        # Range is configurable (DOWN_FROM..DOWN_TO, default 3..9 — Step 9 reads only $W and runs
        # here too). Some sites must stop the distributed graph at Step 7 and run Step 8
        # (mini-CRAM/igv export) elsewhere — e.g. the source CRAM store is a FUSE mount visible only
        # on login/interactive nodes, not batch compute. Set DOWN_TO=7 in cluster.env and run
        # Steps 8-9 separately (see slurm/README.md).
        df="${DOWN_FROM:-3}"; dt="${DOWN_TO:-9}"
        [[ "$df" =~ ^[0-9]+$ && "$dt" =~ ^[0-9]+$ && "$df" -ge 3 && "$dt" -le 9 && "$df" -le "$dt" ]] \
            || { echo "phase.sbatch: DOWN_FROM..DOWN_TO must satisfy 3<=FROM<=TO<=9 (got '$df'..'$dt')" >&2; exit 2; }
        log "Steps $df-$dt"
        rp --from "$df" --to "$dt"
        ;;
    *)
        echo "phase.sbatch: unknown phase '$PHASE'" >&2; exit 2;;
esac
