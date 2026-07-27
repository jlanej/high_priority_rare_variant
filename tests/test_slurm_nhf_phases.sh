#!/usr/bin/env bash
# =============================================================================
# tests/test_slurm_nhf_phases.sh — smoke test for the SLURM glue around Step 8b.
#
# tests/test_nhf_scatter.sh drives 08_igv_export.sh DIRECTLY, so it never executes
# pipeline/slurm/phase.sbatch and cannot catch bugs in the orchestration layer. It missed
# exactly one: the nhf-plan / nhf-scatter phases referenced ${WORK:?}, a variable that is
# never assigned in pipeline/slurm/ and is not in cluster.env.example (only HPRV_WORK is).
# Both phases aborted immediately with "WORK: parameter null or not set" — and at sites that
# DO export $WORK as a scratch path (common on HPC) it silently resolved to the wrong
# directory instead, which is worse.
#
# So this test runs the phases with ONLY the documented variables set, and additionally with
# a DECOY $WORK exported, asserting the phases ignore it. `sbatch` and `run_pipeline.sh` are
# stubbed — we are testing the glue (paths, manifest handoff, array indexing), not the steps.
#
# Run: bash tests/test_slurm_nhf_phases.sh
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
PHASE_SH="$REPO/pipeline/slurm/phase.sbatch"

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
W="$T/run"; mkdir -p "$W" "$T/bin"

# The env submit_slurm.sh writes and phase.sbatch sources. HPRV_CONTAINER_BIN empty => run natively.
cat > "$W/slurm_run.env" <<EOF
HPRV_CONFIG=$T/config.yaml
HPRV_CONTAINER_BIN=
HPRV_SLURM_DIR=$REPO/pipeline/slurm
NHF_CONCURRENCY=4
EOF
: > "$T/config.yaml"

# Stub run_pipeline.sh: records its args, and for --nhf-emit-manifest writes a 2-trio manifest.
cat > "$T/bin/run_pipeline.sh" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> "$RP_LOG"
prev=""
for a in "$@"; do
    [[ "$prev" == "--nhf-emit-manifest" ]] && printf 'T1\nT2\n' > "$a"
    prev="$a"
done
STUB
chmod +x "$T/bin/run_pipeline.sh"

# Stub sbatch: records the submission and prints a job id (--parsable).
cat > "$T/bin/sbatch" <<'STUB'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$SBATCH_LOG"
echo 12345
STUB
chmod +x "$T/bin/sbatch"

export RP_LOG="$T/rp.log" SBATCH_LOG="$T/sbatch.log"
: > "$RP_LOG"; : > "$SBATCH_LOG"

fail=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fail=1; fi; }

# ---------------------------------------------------------------------------
# 1. nhf-plan with ONLY the documented variables (HPRV_WORK). WORK deliberately UNSET.
# ---------------------------------------------------------------------------
rc=0
env -u WORK PATH="$T/bin:$PATH" HPRV_WORK="$W" RP_LOG="$RP_LOG" SBATCH_LOG="$SBATCH_LOG" \
    bash "$PHASE_SH" nhf-plan > "$T/plan.out" 2>&1 || rc=$?
chk "nhf-plan runs with only HPRV_WORK set (WORK unset)" '[[ "$rc" -eq 0 ]]'
chk "nhf-plan did not abort on an unset variable" \
    '! grep -qi "parameter null or not set" "$T/plan.out"'
chk "nhf-plan wrote the manifest under \$HPRV_WORK/igv/" \
    '[[ -s "$W/igv/nhf_trios.txt" ]]'
chk "nhf-plan submitted a scatter array sized to the manifest (2 trios => 0-1)" \
    'grep -q -- "--array=0-1" "$SBATCH_LOG"'
chk "nhf-plan submitted a dependent gather" \
    'grep -q "nhf-gather" "$SBATCH_LOG"'
chk "gather depends afterany (a failed trio must not stall the join)" \
    'grep -q -- "--dependency=afterany" "$SBATCH_LOG"'

# ---------------------------------------------------------------------------
# 2. nhf-scatter maps $SLURM_ARRAY_TASK_ID -> the right manifest line
# ---------------------------------------------------------------------------
: > "$RP_LOG"
rc=0
env -u WORK PATH="$T/bin:$PATH" HPRV_WORK="$W" SLURM_ARRAY_TASK_ID=1 \
    RP_LOG="$RP_LOG" SBATCH_LOG="$SBATCH_LOG" \
    bash "$PHASE_SH" nhf-scatter > "$T/scatter.out" 2>&1 || rc=$?
chk "nhf-scatter runs with only HPRV_WORK set" '[[ "$rc" -eq 0 ]]'
chk "nhf-scatter task 1 selected the SECOND manifest trio (T2)" \
    'grep -q -- "--nhf-trio T2" "$RP_LOG"'
chk "nhf-scatter asked for Step 8 only" 'grep -q -- "--from 8 --to 8" "$RP_LOG"'

# ---------------------------------------------------------------------------
# 3. THE REGRESSION: a decoy $WORK must be ignored, not silently used.
#    On HPC $WORK is frequently exported as a scratch filesystem.
# ---------------------------------------------------------------------------
DECOY="$T/decoy"; mkdir -p "$DECOY"
: > "$RP_LOG"; : > "$SBATCH_LOG"
rc=0
PATH="$T/bin:$PATH" WORK="$DECOY" HPRV_WORK="$W" RP_LOG="$RP_LOG" SBATCH_LOG="$SBATCH_LOG" \
    bash "$PHASE_SH" nhf-plan > "$T/plan2.out" 2>&1 || rc=$?
chk "nhf-plan ignores a decoy \$WORK and still succeeds" '[[ "$rc" -eq 0 ]]'
chk "nhf-plan wrote NOTHING under the decoy \$WORK" \
    '[[ ! -e "$DECOY/igv" ]]'
chk "the manifest is still under \$HPRV_WORK" '[[ -s "$W/igv/nhf_trios.txt" ]]'

# ---------------------------------------------------------------------------
# 4. nhf-gather is join-only
# ---------------------------------------------------------------------------
: > "$RP_LOG"
rc=0
env -u WORK PATH="$T/bin:$PATH" HPRV_WORK="$W" RP_LOG="$RP_LOG" SBATCH_LOG="$SBATCH_LOG" \
    bash "$PHASE_SH" nhf-gather > "$T/gather.out" 2>&1 || rc=$?
chk "nhf-gather runs with only HPRV_WORK set" '[[ "$rc" -eq 0 ]]'
chk "nhf-gather passed --nhf-gather" 'grep -q -- "--nhf-gather" "$RP_LOG"'

# ---------------------------------------------------------------------------
# 5. empty manifest => plan short-circuits to a gather instead of a 0-task array
# ---------------------------------------------------------------------------
cat > "$T/bin/run_pipeline.sh" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> "$RP_LOG"
prev=""
for a in "$@"; do [[ "$prev" == "--nhf-emit-manifest" ]] && : > "$a"; prev="$a"; done
STUB
chmod +x "$T/bin/run_pipeline.sh"
: > "$RP_LOG"; : > "$SBATCH_LOG"
rc=0
env -u WORK PATH="$T/bin:$PATH" HPRV_WORK="$W" RP_LOG="$RP_LOG" SBATCH_LOG="$SBATCH_LOG" \
    bash "$PHASE_SH" nhf-plan > "$T/plan3.out" 2>&1 || rc=$?
chk "nhf-plan with nothing to do exits 0" '[[ "$rc" -eq 0 ]]'
chk "nhf-plan with nothing to do submits NO array" '! grep -q -- "--array" "$SBATCH_LOG"'
chk "nhf-plan with nothing to do still runs the gather" 'grep -q -- "--nhf-gather" "$RP_LOG"'

[[ "$fail" -eq 0 ]] && echo "All SLURM Step-8b phase tests passed." \
    || { echo "test_slurm_nhf_phases FAILED"; exit 1; }
