#!/usr/bin/env bash
# =============================================================================
# tests/test_slurm_rescore_phases.sh — smoke test for the SLURM glue around Step 5b.
#
# The `calls` phase runs Steps 3-5 with --rescore-emit-manifest, then must submit the Step-5b chunk
# array (sized from the manifest), a dependent rescore-gather, and the downstream tail — or, on an
# EMPTY manifest, downstream alone. rescore-scatter maps $SLURM_ARRAY_TASK_ID onto --rescore-chunk;
# downstream defaults to Steps 6..DOWN_TO because 3-5 ran in `calls`. `sbatch` and `run_pipeline.sh`
# are stubbed — this tests the glue (manifest handoff, array sizing, dependency edges, the empty and
# the missing-manifest branches), not the steps. Pure bash, no bio deps; a DECOY $WORK is exported
# in one run to assert the phases only ever read HPRV_WORK (the Step-8b regression).
#
# Run: bash tests/test_slurm_rescore_phases.sh
# =============================================================================
# shellcheck disable=SC2034  # `rc` is read inside chk's eval strings
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
PHASE_SH="$REPO/pipeline/slurm/phase.sbatch"

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
W="$T/run"; mkdir -p "$W" "$T/bin"

cat > "$W/slurm_run.env" <<EOF
HPRV_CONFIG=$T/config.yaml
HPRV_CONTAINER_BIN=
HPRV_SLURM_DIR=$REPO/pipeline/slurm
RESCORE_CONCURRENCY=5
DOWN_TO=8
EOF
: > "$T/config.yaml"

# Stub run_pipeline.sh: records its args; on --rescore-emit-manifest F writes RESCORE_N chunk paths
# (default 3) — or nothing at all when RESCORE_NO_MANIFEST is set, to exercise the loud failure.
cat > "$T/bin/run_pipeline.sh" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> "$RP_LOG"
prev=""
for a in "$@"; do
    if [[ "$prev" == "--rescore-emit-manifest" && -z "${RESCORE_NO_MANIFEST:-}" ]]; then
        : > "$a"
        i=0
        while [[ "$i" -lt "${RESCORE_N:-3}" ]]; do
            printf '%s/chunks/chunk_%04d.vcf\n' "$(dirname "$a")" "$i" >> "$a"; i=$((i + 1))
        done
    fi
    prev="$a"
done
STUB
# Stub sbatch: records the submission, prints an incrementing job id (--parsable).
cat > "$T/bin/sbatch" <<'STUB'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$SBATCH_LOG"
n=$(( $(cat "$SB_CTR" 2>/dev/null || echo 100) + 1 )); echo "$n" > "$SB_CTR"; echo "$n"
STUB
chmod +x "$T/bin/run_pipeline.sh" "$T/bin/sbatch"
export RP_LOG="$T/rp.log" SBATCH_LOG="$T/sbatch.log" SB_CTR="$T/ctr"
reset_logs() { : > "$RP_LOG"; : > "$SBATCH_LOG"; rm -f "$SB_CTR" "$W/slurm_jobids.txt" "$W/spliceai_rescore/manifest.txt"; }

fail=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fail=1; fi; }
run_phase() {  # $1 = phase, then extra env assignments
    local phase="$1"; shift
    env -u WORK "$@" PATH="$T/bin:$PATH" HPRV_WORK="$W" RP_LOG="$RP_LOG" SBATCH_LOG="$SBATCH_LOG" SB_CTR="$SB_CTR" \
        bash "$PHASE_SH" "$phase"
}

# ---------------------------------------------------------------------------
# 1. `calls` with 3 pending chunks: Steps 3-5 + plan, then array -> gather -> downstream.
#    A DECOY $WORK is exported: the phase must read HPRV_WORK only.
# ---------------------------------------------------------------------------
reset_logs
rc=0; run_phase calls WORK=/nonexistent/decoy RESCORE_N=3 > "$T/calls.out" 2>&1 || rc=$?
chk "calls exits 0 with 3 pending chunks" '[[ "$rc" -eq 0 ]]'
chk "calls ran Steps 3-5 once with --rescore-emit-manifest into HPRV_WORK (not \$WORK)" \
    'grep -q -- "--from 3 --to 5 --rescore-emit-manifest $W/spliceai_rescore/manifest.txt" "$RP_LOG" && ! grep -q nonexistent "$RP_LOG" "$SBATCH_LOG"'
chk "calls submitted a rescore-scatter array sized from the manifest (0-2) at RESCORE_CONCURRENCY" \
    'grep "phase.sbatch rescore-scatter" "$SBATCH_LOG" | grep -q -- "--array=0-2%5"'
chk "rescore-gather depends afterok on the array" \
    'grep "phase.sbatch rescore-gather" "$SBATCH_LOG" | grep -q -- "--dependency=afterok:101"'
chk "downstream depends afterok on rescore-gather" \
    'grep "phase.sbatch downstream" "$SBATCH_LOG" | grep -q -- "--dependency=afterok:102"'
chk "exactly three submissions (array, gather, downstream)" '[[ "$(grep -c . "$SBATCH_LOG")" -eq 3 ]]'
chk "job ids recorded for the operator" \
    'grep -q "^rescore-scatter 101" "$W/slurm_jobids.txt" && grep -q "^rescore-gather 102" "$W/slurm_jobids.txt" && grep -q "^downstream 103" "$W/slurm_jobids.txt"'
chk "every submission runs by the shared-storage phase.sbatch path" \
    '[[ "$(grep -c "$REPO/pipeline/slurm/phase.sbatch" "$SBATCH_LOG")" -eq 3 ]]'

# ---------------------------------------------------------------------------
# 2. rescore-scatter maps the array index onto --rescore-chunk; Step 5 proper is not re-run
#    (run_pipeline.sh handles that — here we only check the argument the phase hands it).
# ---------------------------------------------------------------------------
reset_logs
rc=0; run_phase rescore-scatter SLURM_ARRAY_TASK_ID=2 > /dev/null 2>&1 || rc=$?
chk "rescore-scatter task 2 -> run_pipeline --from 5 --to 5 --rescore-chunk 2" \
    '[[ "$rc" -eq 0 ]] && grep -q -- "--from 5 --to 5 --rescore-chunk 2" "$RP_LOG"'

# ---------------------------------------------------------------------------
# 3. rescore-gather -> --rescore-gather; downstream -> Steps 6..DOWN_TO by default.
# ---------------------------------------------------------------------------
reset_logs
rc=0; run_phase rescore-gather > /dev/null 2>&1 || rc=$?
chk "rescore-gather -> run_pipeline --from 5 --to 5 --rescore-gather" \
    '[[ "$rc" -eq 0 ]] && grep -q -- "--from 5 --to 5 --rescore-gather" "$RP_LOG"'
reset_logs
rc=0; run_phase downstream > /dev/null 2>&1 || rc=$?
chk "downstream defaults to Steps 6-8 (3-5 ran in calls; DOWN_TO=8)" \
    '[[ "$rc" -eq 0 ]] && grep -q -- "--from 6 --to 8" "$RP_LOG"'

# ---------------------------------------------------------------------------
# 4. `calls` with an EMPTY manifest (rescoring disabled / everything cached): downstream only.
# ---------------------------------------------------------------------------
reset_logs
rc=0; run_phase calls RESCORE_N=0 > "$T/calls_empty.out" 2>&1 || rc=$?
chk "empty manifest: calls exits 0" '[[ "$rc" -eq 0 ]]'
chk "empty manifest: no rescore array or gather submitted" '! grep -q "rescore-" "$SBATCH_LOG"'
chk "empty manifest: downstream submitted directly, once" \
    '[[ "$(grep -c "phase.sbatch downstream" "$SBATCH_LOG")" -eq 1 && "$(grep -c . "$SBATCH_LOG")" -eq 1 ]]'

# ---------------------------------------------------------------------------
# 5. `calls` when Step 5b wrote NO manifest at all: fail loudly (an absent file is never
#    "nothing to do" — that was the green no-op the nhf-plan phase once had). A STALE manifest
#    from an earlier attempt is planted first: the phase must remove it, not reuse it.
# ---------------------------------------------------------------------------
reset_logs
mkdir -p "$W/spliceai_rescore"; printf '%s/chunks/chunk_0000.vcf\n' "$W/spliceai_rescore" > "$W/spliceai_rescore/manifest.txt"
rc=0; run_phase calls RESCORE_NO_MANIFEST=1 > "$T/calls_none.out" 2>&1 || rc=$?
chk "missing manifest: calls fails" '[[ "$rc" -ne 0 ]]'
chk "missing manifest: the message names the missing file" 'grep -q "wrote no chunk manifest" "$T/calls_none.out"'
chk "missing manifest: nothing was submitted" '[[ ! -s "$SBATCH_LOG" ]]'

[[ "$fail" -eq 0 ]] && echo "ALL SLURM STEP-5B PHASE ASSERTIONS PASSED" || { echo "SLURM STEP-5B PHASE TEST FAILED"; exit 1; }
