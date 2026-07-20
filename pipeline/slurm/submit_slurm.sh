#!/usr/bin/env bash
# =============================================================================
# submit_slurm.sh — submit the hprv pipeline as a coherent SLURM job graph.
#
#   prep ──afterok──> plan ──submits──> scatter[array] ──afterok──> gather ──afterok──> downstream
#
# Run ONCE from a login node:  ./submit_slurm.sh path/to/cluster.env
# It submits `prep` and `plan`; `plan` (after prep) enumerates the union's contigs
# and submits the scatter array + gather + downstream. Coherence: every edge is an
# --dependency=afterok, so any failure halts everything downstream (no partial call
# set), and gather independently re-verifies all shards.
#
# RESUME after a walltime kill or a failed contig: just re-run this script (or
# re-submit the failed scatter array elements). Completed pieces have .done files
# and are skipped; only unfinished work re-runs. To FORCE a full re-annotation,
# remove $HPRV_WORK/annotate_shards/ first.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ENVFILE="${1:?usage: submit_slurm.sh <cluster.env>}"
[[ -f "$ENVFILE" ]] || { echo "cluster.env not found: $ENVFILE" >&2; exit 1; }
# shellcheck source=/dev/null
source "$ENVFILE"

command -v sbatch >/dev/null 2>&1 || { echo "sbatch not on PATH — run this from a SLURM login node" >&2; exit 1; }

# Required settings — fail loudly rather than submitting a broken graph.
for v in HPRV_CONFIG HPRV_WORK SLURM_ACCOUNT SLURM_PARTITION; do
    [[ -n "${!v:-}" ]] || { echo "cluster.env: $v is not set" >&2; exit 1; }
done
[[ -f "$HPRV_CONFIG" ]] || { echo "HPRV_CONFIG not found: $HPRV_CONFIG" >&2; exit 1; }
if [[ -n "${HPRV_CONTAINER_BIN:-}" ]]; then
    [[ -f "${HPRV_SIF:-}" ]] || { echo "HPRV_SIF not found: ${HPRV_SIF:-<unset>} (needed for $HPRV_CONTAINER_BIN)" >&2; exit 1; }
fi

mkdir -p "$HPRV_WORK"
# Snapshot the resolved settings where every phase job can read them (shared storage).
cp "$ENVFILE" "$HPRV_WORK/slurm_run.env"

sb() { sbatch --parsable \
        --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" ${SLURM_QOS:+--qos="$SLURM_QOS"} \
        --export=ALL,HPRV_WORK="$HPRV_WORK" "$@"; }

prep=$(sb --cpus-per-task="${PREP_CPUS:-8}" --mem="${PREP_MEM:-32G}" --time="${PREP_TIME:-24:00:00}" \
        --job-name=hprv-prep "$HERE/phase.sbatch" prep)
plan=$(sb --dependency=afterok:"$prep" \
        --cpus-per-task="${PLAN_CPUS:-1}" --mem="${PLAN_MEM:-4G}" --time="${PLAN_TIME:-00:20:00}" \
        --job-name=hprv-plan "$HERE/phase.sbatch" plan)

cat <<EOF
Submitted the hprv SLURM graph:
  prep       = $prep   (resolve + QC + cohort union)
  plan       = $plan   (after prep: enumerates contigs, submits scatter[array] -> gather -> downstream)
The scatter/gather/downstream job IDs are chosen by 'plan' at runtime; find them with
  squeue -u \$USER --name=hprv-scatter,hprv-gather,hprv-down
or in $HPRV_WORK/slurm_jobids.txt after plan runs.
Re-run this script to resume; rm $HPRV_WORK/annotate_shards/ to force a full re-annotation.
EOF
