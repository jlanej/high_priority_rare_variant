#!/usr/bin/env bash
# =============================================================================
# test_plan_resubmit_path.sh — regression test for the SLURM `plan` self-resubmit bug.
#
# Under sbatch, SLURM copies the batch script into the job spool dir and runs it as
# `slurm_script`, so ${BASH_SOURCE[0]} inside a running job points at the spool — NOT
# the script on shared storage. If `plan` re-submitted scatter/gather/downstream by a
# BASH_SOURCE-derived path, sbatch would fail with "Unable to open file .../phase.sbatch".
# The fix carries the real dir (HPRV_SLURM_DIR, resolved on the login node) through
# slurm_run.env. This test reproduces the spool scenario with a fake sbatch and asserts:
#   (1) with the fix, `plan` re-submits by the real shared-storage path (which exists);
#   (2) without it, the guard fails loudly instead of emitting a broken path.
#
# Pure bash — no bio deps — so it runs in CI's host-only test job.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SLURM_SRC="$REPO/pipeline/slurm"          # the REAL scripts under test
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
fail=0

# --- fakes on PATH: sbatch (logs args, emits ids), run_pipeline.sh (writes a manifest) ---
mkdir -p "$T/bin"
cat > "$T/bin/sbatch" <<'FAKE'
#!/usr/bin/env bash
echo "SBATCH $*" >> "$SB_LOG"
n=$(( $(cat "$SB_CTR" 2>/dev/null || echo 2000) + 1 )); echo "$n" > "$SB_CTR"; echo "$n"
FAKE
cat > "$T/bin/run_pipeline.sh" <<'FAKE'
#!/usr/bin/env bash
mf=""; while [[ $# -gt 0 ]]; do case "$1" in --annotate-emit-manifest) mf="$2"; shift 2;; *) shift;; esac; done
[[ -n "$mf" ]] && printf 'chr1\nchr2\nchrX\n' > "$mf"; exit 0
FAKE
chmod +x "$T/bin/sbatch" "$T/bin/run_pipeline.sh"
export PATH="$T/bin:$PATH" SB_LOG="$T/sbatch.log" SB_CTR="$T/ctr"

# --- a native (no-container) cluster.env, and a submit that writes slurm_run.env ---
mkdir -p "$T/work"; touch "$T/config.yaml"
cat > "$T/cluster.env" <<EOF
HPRV_CONFIG=$T/config.yaml
HPRV_WORK=$T/work
SLURM_ACCOUNT=acct
SLURM_PARTITION=part
HPRV_CONTAINER_BIN=
SCATTER_CPUS=16; SCATTER_MEM=32G; SCATTER_TIME=06:00:00
GATHER_CPUS=4; GATHER_MEM=16G; GATHER_TIME=04:00:00
DOWN_CPUS=8; DOWN_MEM=32G; DOWN_TIME=24:00:00
SCATTER_CONCURRENCY=12
EOF
: > "$T/sbatch.log"
bash "$SLURM_SRC/submit_slurm.sh" "$T/cluster.env" >/dev/null 2>&1 || true

check() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fail=1; fi; }

# (0) submit_slurm.sh must have recorded HPRV_SLURM_DIR = the real slurm dir.
check "submit wrote HPRV_SLURM_DIR" \
  'grep -q "^HPRV_SLURM_DIR=" "$T/work/slurm_run.env" && eval "$(grep ^HPRV_SLURM_DIR= "$T/work/slurm_run.env")"; [[ "$HPRV_SLURM_DIR" -ef "$SLURM_SRC" ]]'

# --- simulate SLURM running the plan job FROM THE SPOOL DIR (BASH_SOURCE points there) ---
spool="$T/spool/job12345"; mkdir -p "$spool"
cp "$SLURM_SRC/phase.sbatch" "$spool/slurm_script"    # SLURM's copy, renamed as it does

# (1) WITH the fix: plan resubmits by the REAL dir, and that path exists (not the spool).
: > "$T/sbatch.log"
HPRV_WORK="$T/work" bash "$spool/slurm_script" plan >/dev/null 2>&1 || echo "  (plan exit $?)"
echo "recorded re-submissions:"; grep -o 'phase.sbatch [a-z]*' "$T/sbatch.log" | sed 's/^/    /' || true
check "plan re-submits by the shared-storage path (not the spool)" \
  'grep -q "$SLURM_SRC/phase.sbatch scatter" "$T/sbatch.log" && ! grep -q "$spool/phase.sbatch" "$T/sbatch.log"'
check "the re-submitted path actually exists" '[[ -f "$SLURM_SRC/phase.sbatch" ]]'
check "all three phases re-submitted (scatter/gather/downstream)" \
  'grep -q "phase.sbatch scatter" "$T/sbatch.log" && grep -q "phase.sbatch gather" "$T/sbatch.log" && grep -q "phase.sbatch downstream" "$T/sbatch.log"'

# (2) WITHOUT the fix (HPRV_SLURM_DIR stripped): SLURM_DIR falls back to $HERE=spool, whose
#     phase.sbatch does not exist -> the guard must fail LOUDLY, not emit a broken sbatch path.
grep -v '^HPRV_SLURM_DIR=' "$T/work/slurm_run.env" > "$T/work/slurm_run.env.tmp" && mv "$T/work/slurm_run.env.tmp" "$T/work/slurm_run.env"
: > "$T/sbatch.log"
if HPRV_WORK="$T/work" bash "$spool/slurm_script" plan >"$T/neg.out" 2>&1; then
  echo "FAIL guard fires when the resubmit path is wrong (plan unexpectedly succeeded)"; fail=1
else
  grep -qi "cannot re-submit" "$T/neg.out" \
    && echo "PASS guard fires with an actionable message when the resubmit path is wrong" \
    || { echo "FAIL guard message missing"; cat "$T/neg.out"; fail=1; }
  ! grep -q "phase.sbatch scatter" "$T/sbatch.log" \
    && echo "PASS no phases submitted once the path is known-bad" \
    || { echo "FAIL submitted phases despite a bad path"; fail=1; }
fi

[[ "$fail" -eq 0 ]] && echo "ALL SLURM PLAN-RESUBMIT ASSERTIONS PASSED" || { echo "SLURM PLAN-RESUBMIT TEST FAILED"; exit 1; }
