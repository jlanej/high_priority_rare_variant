#!/usr/bin/env bash
# =============================================================================
# prepare_resources_stage.sh — run the resource preparation ONCE as an idempotent
# upstream stage, gated by a .done sentinel, so downstream pipeline runs can
# depend on it. Wraps scripts/prepare_resources.sh (fetch -> verify -> emit-env)
# and only stamps success when `verify` passes.
#
# Idempotent: re-running with everything already cached is a fast no-op that just
# re-verifies. Safe to make a hard prerequisite of run_pipeline.sh.
#
#   scripts/prepare_resources_stage.sh --dir /path/to/hprv_resources [--accept-license] [--only id,id]
#
# INTERNET NOTE: the `fetch` step DOWNLOADS large resources, so it must run where
# the node has outbound internet. On many HPC clusters (e.g. MSI) batch/compute
# nodes have NO internet — run this on a login / data-transfer / interactive node
# (or a partition that permits egress). The heavy CPU/RAM work (gnomAD slimming,
# dbNSFP sort) benefits from many cores. Run it INSIDE the container so bcftools/
# tabix/samtools are on PATH, e.g.:
#
#   apptainer exec --cleanenv --bind "$RESDIR" "$HPRV_IMAGE" \
#       scripts/prepare_resources_stage.sh --dir "$RESDIR" --accept-license
#
# It can also be submitted with sbatch IF the chosen partition has internet.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DIR="" ; PASS_THRU=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dir) DIR="$2"; PASS_THRU+=(--dir "$2"); shift 2;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) PASS_THRU+=("$1"); shift;;
    esac
done
[[ -n "$DIR" ]] || { echo "ERROR: need --dir DIR" >&2; exit 2; }
mkdir -p "$DIR"
DONE="$DIR/.prepared.done"
ENVOUT="$DIR/resources.env"

# Idempotency: if the sentinel exists AND verify still passes, do nothing.
if [[ -f "$DONE" ]] && bash "$HERE/prepare_resources.sh" --dir "$DIR" verify >/dev/null 2>&1; then
    echo "[prep-stage] resources already prepared and verified: $DIR" >&2
    echo "[prep-stage] env: $ENVOUT" >&2
    exit 0
fi
rm -f "$DONE"

echo "[prep-stage] fetching/preparing resources into $DIR ..." >&2
bash "$HERE/prepare_resources.sh" "${PASS_THRU[@]}" fetch

echo "[prep-stage] verifying ..." >&2
bash "$HERE/prepare_resources.sh" --dir "$DIR" verify

echo "[prep-stage] emitting env -> $ENVOUT" >&2
bash "$HERE/prepare_resources.sh" --dir "$DIR" emit-env --out "$ENVOUT"

touch "$DONE"
echo "[prep-stage] DONE. Downstream runs may now depend on: $DONE" >&2
echo "[prep-stage] source the env before run_pipeline.sh:  source $ENVOUT" >&2
