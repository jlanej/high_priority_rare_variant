#!/usr/bin/env bash
# =============================================================================
# tests/test_spliceai_rescore.sh — Step 5b end to end against a STUB SpliceAI.
#
# Drives the REAL pipeline/05b_spliceai_rescore.sh + src/hprv/spliceai_rescore.py with a stub
# `micromamba run -p ENV spliceai` that writes deterministic SpliceAI= values, so the whole path —
# plan -> chunk -> score -> gather -> merge into candidates.calls.tsv — runs with nothing but bash and
# python3. Asserts:
#   * the spliceai_wide_* block lands after the curated spliceai_* columns, with the event decomposed
#     (a same-site gain+loss pair -> in-frame cryptic shift; a lone distal loss -> site_loss, distal=1)
#   * the per-variant cache: a re-run scores NOTHING new (Step 5 rewrites the calls table every run,
#     so the cache, not a file marker, is what makes 5b cheap), a re-merge never duplicates the block,
#     a changed window invalidates the cache
#   * serial == distributed (emit-manifest / --chunk N / --gather) byte for byte
#   * exit 3 when the isolated env is missing; gather refuses a PARTIAL merge
#
# Run: bash tests/test_spliceai_rescore.sh
# =============================================================================
# shellcheck disable=SC2034  # `rc` is read inside chk's eval strings
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
STEP="$REPO/pipeline/05b_spliceai_rescore.sh"

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
mkdir -p "$T/bin" "$T/env/bin" "$T/w" "$T/tmp" "$T/audit"
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"
export HPRV_RUNTIME=native HPRV_TMPDIR="$T/tmp" HPRV_AUDIT_DIR="$T/audit" HPRV_SPLICEAI_ENV="$T/env"
export SAI_LOG="$T/spliceai_calls.log"; : > "$SAI_LOG"

# --- stub `micromamba run -p PREFIX cmd args...` -> exec PREFIX/bin/cmd ---------------------------
cat > "$T/bin/micromamba" <<'STUB'
#!/usr/bin/env bash
[[ "${1:-}" == "run" ]] || { echo "stub micromamba: only 'run' is supported ($*)" >&2; exit 2; }
shift; prefix=""
while [[ $# -gt 0 ]]; do case "$1" in -p|--prefix) prefix="$2"; shift 2;; -n|--name) shift 2;; *) break;; esac; done
cmd="$1"; shift
exec "$prefix/bin/$cmd" "$@"
STUB
# --- stub spliceai: -I in -O out -D dist. Deterministic by position parity; logs every call. ----
#   even POS: donor loss 0.70 @+2 and donor gain 0.40 @+14 (+ a weaker second gene) -> the merge must
#             keep GENEA's entry whole and read a 12-nt in-frame cryptic shift
#   odd  POS: acceptor loss 0.60 @-120 -> site_loss, beyond the precomputed +/-50 bp window
cat > "$T/env/bin/spliceai" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
in="" out="" d=""
while [[ $# -gt 0 ]]; do case "$1" in -I) in="$2"; shift 2;; -O) out="$2"; shift 2;; -D) d="$2"; shift 2;; -R|-A) shift 2;; *) shift;; esac; done
printf '%s\t%s\t%s\n' "$in" "$d" "${OMP_NUM_THREADS:-unset}" >> "$SAI_LOG"
awk -F'\t' -v OFS='\t' '
  /^##/    { print; next }
  /^#CHROM/{ print "##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"stub\">"; print; next }
  { pos = $2 + 0
    if (pos % 2 == 0) info = "SpliceAI=" $5 "|GENEA|0.00|0.00|0.40|0.70|-3|-2|14|2," $5 "|GENEB|0.05|0.00|0.00|0.10|1|1|1|1"
    else              info = "SpliceAI=" $5 "|GENEA|0.00|0.60|0.00|0.00|10|-120|3|4"
    $8 = info; print }' "$in" > "$out"
STUB
chmod +x "$T/bin/micromamba" "$T/env/bin/spliceai"
export PATH="$T/bin:$PATH"

# --- a reference (only its .fai is read, for contig header lines) and a calls table ------------
printf '>chr1\nACGT\n>chr2\nACGT\n' > "$T/ref.fa"
printf 'chr1\t100000\t6\t4\t5\nchr2\t100000\t100013\t4\t5\n' > "$T/ref.fa.fai"
CALLS="$T/w/candidates.calls.tsv"
{ printf 'trio_id\tmode\tchrom\tpos\tref\talt\tcadd\tspliceai_ds\tspliceai_event\tinfo_vep_HGVSc\n'
  printf 'T1\tdominant\tchr1\t1000\tA\tT\t3\t0.55\tdonor_loss\tENST1:c.1A>T\n'      # even -> cryptic shift
  printf 'T1\tdominant\tchr1\t1001\tC\tG\t\t\t\tENST1:c.2C>G\n'                    # odd  -> distal site loss
  printf 'T2\tdominant\tchr1\t1000\tA\tT\t3\t0.55\tdonor_loss\tENST1:c.1A>T\n'      # same variant, 2nd trio
  printf 'T1\tcompound_het\tchr2\t2003\tG\t*\t\t\t\t\n'                             # symbolic: never scored
  printf 'T1\thom_recessive\tchr2\t2002\tG\tGA\t\t0.1\t\tENST2:c.9dup\n'            # even indel
} > "$CALLS"
cp "$CALLS" "$T/calls.orig.tsv"

fail=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fail=1; fi; }
cell() {  # cell FILE ROW(1-based data row) COLUMN -> value
    python3 - "$1" "$2" "$3" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
print(rows[int(sys.argv[2]) - 1].get(sys.argv[3], "<nocol>"))
PY
}
hdr() { head -1 "$1" | tr '\t' '\n'; }

# ---------------------------------------------------------------------------
# A. the serial path: plan -> score every chunk -> gather/merge
# ---------------------------------------------------------------------------
rc=0; bash "$STEP" --calls "$CALLS" --ref "$T/ref.fa" --outdir "$T/w/spliceai_rescore" \
        --distance 4999 --chunk-size 2 --threads 3 > "$T/A.log" 2>&1 || rc=$?
chk "A: serial run exits 0" '[[ "$rc" -eq 0 ]]'
chk "A: 3 distinct non-symbolic variants -> 2 chunks -> 2 model calls, each at -D 4999" \
    '[[ "$(wc -l < "$SAI_LOG" | tr -d " ")" -eq 2 ]] && [[ "$(cut -f2 "$SAI_LOG" | sort -u)" == "4999" ]]'
chk "A: TensorFlow thread pools sized to --threads" '[[ "$(cut -f3 "$SAI_LOG" | sort -u)" == "3" ]]'
chk "A: the wide block sits right after the curated spliceai_* columns, before info_*" \
    '[[ "$(hdr "$CALLS" | paste -sd, -)" == "trio_id,mode,chrom,pos,ref,alt,cadd,spliceai_ds,spliceai_event,spliceai_wide_ds,spliceai_wide_event,spliceai_wide_event_pos,spliceai_wide_event2,spliceai_wide_event2_ds,spliceai_wide_event2_pos,spliceai_wide_shift_nt,spliceai_wide_shift_frame,spliceai_wide_effect,spliceai_wide_symbol,spliceai_wide_distal,spliceai_wide_minus_precomputed,info_vep_HGVSc" ]]'
chk "A: even-position variant: donor loss 0.7 @1002 + donor gain @1014 -> 12-nt IN-FRAME cryptic shift, GENEA (the strongest gene) wins" \
    '[[ "$(cell "$CALLS" 1 spliceai_wide_ds)" == "0.7" && "$(cell "$CALLS" 1 spliceai_wide_event)" == "donor_loss" && "$(cell "$CALLS" 1 spliceai_wide_event_pos)" == "1002" && "$(cell "$CALLS" 1 spliceai_wide_event2)" == "donor_gain" && "$(cell "$CALLS" 1 spliceai_wide_shift_nt)" == "12" && "$(cell "$CALLS" 1 spliceai_wide_shift_frame)" == "in_frame" && "$(cell "$CALLS" 1 spliceai_wide_effect)" == "cryptic_shift_in_frame" && "$(cell "$CALLS" 1 spliceai_wide_symbol)" == "GENEA" ]]'
chk "A: within the precomputed window -> distal=0; wide minus precomputed = 0.7-0.55 = 0.15" \
    '[[ "$(cell "$CALLS" 1 spliceai_wide_distal)" == "0" && "$(cell "$CALLS" 1 spliceai_wide_minus_precomputed)" == "0.15" ]]'
chk "A: odd-position variant: acceptor loss 0.6 @-120 -> site_loss at 881, DISTAL (beyond +/-50 bp), no delta without a precomputed score" \
    '[[ "$(cell "$CALLS" 2 spliceai_wide_ds)" == "0.6" && "$(cell "$CALLS" 2 spliceai_wide_event)" == "acceptor_loss" && "$(cell "$CALLS" 2 spliceai_wide_event_pos)" == "881" && "$(cell "$CALLS" 2 spliceai_wide_effect)" == "site_loss" && "$(cell "$CALLS" 2 spliceai_wide_distal)" == "1" && "$(cell "$CALLS" 2 spliceai_wide_minus_precomputed)" == "" ]]'
chk "A: the same variant in a second trio carries identical wide cells" \
    '[[ "$(cell "$CALLS" 3 spliceai_wide_effect)" == "cryptic_shift_in_frame" && "$(cell "$CALLS" 3 spliceai_wide_ds)" == "0.7" ]]'
chk "A: a symbolic ALT is never scored: blank wide cells, not zeros" \
    '[[ "$(cell "$CALLS" 4 spliceai_wide_ds)" == "" && "$(cell "$CALLS" 4 spliceai_wide_effect)" == "" && "$(cell "$CALLS" 4 spliceai_wide_distal)" == "" ]]'
chk "A: the indel row (precomputed 0.1) reads wide 0.7 and a 0.6 delta" \
    '[[ "$(cell "$CALLS" 5 spliceai_wide_ds)" == "0.7" && "$(cell "$CALLS" 5 spliceai_wide_minus_precomputed)" == "0.6" ]]'
chk "A: the curated columns are untouched (row count, spliceai_ds, info_ block)" \
    '[[ "$(wc -l < "$CALLS" | tr -d " ")" -eq 6 && "$(cell "$CALLS" 1 spliceai_ds)" == "0.55" && "$(cell "$CALLS" 1 info_vep_HGVSc)" == "ENST1:c.1A>T" ]]'
chk "A: the per-variant cache holds 3 variants at -D 4999 and the chunk dir is cleaned" \
    'grep -q "^#spliceai_rescore.distance=4999" "$T/w/spliceai_rescore/scores.tsv" && [[ "$(grep -cv "^#" "$T/w/spliceai_rescore/scores.tsv")" -eq 4 ]] && [[ ! -f "$T/w/spliceai_rescore/manifest.txt" ]]'
chk "A: audit records the run" 'grep -q "05b_spliceai_rescore.*variants_pending.3" "$HPRV_AUDIT_DIR/counts.tsv"'
cp "$CALLS" "$T/serial.tsv"

# ---------------------------------------------------------------------------
# B. Step 5 re-ran (the calls table is rewritten WITHOUT the block) -> 5b must re-merge from the
#    cache and score nothing; and a re-merge on the augmented table must not duplicate the block.
# ---------------------------------------------------------------------------
cp "$T/calls.orig.tsv" "$CALLS"
rc=0; bash "$STEP" --calls "$CALLS" --ref "$T/ref.fa" --outdir "$T/w/spliceai_rescore" \
        --distance 4999 --chunk-size 2 --threads 1 > "$T/B.log" 2>&1 || rc=$?
chk "B: re-run after a Step-5 rewrite exits 0 and scores NOTHING new (cache hit)" \
    '[[ "$rc" -eq 0 && "$(wc -l < "$SAI_LOG" | tr -d " ")" -eq 2 ]]'
chk "B: the re-merged table is byte-identical to the first" 'cmp -s "$CALLS" "$T/serial.tsv"'
rc=0; bash "$STEP" --calls "$CALLS" --ref "$T/ref.fa" --outdir "$T/w/spliceai_rescore" \
        --distance 4999 --chunk-size 2 --threads 1 > "$T/B2.log" 2>&1 || rc=$?
chk "B: a run on the ALREADY-augmented table replaces the block, never duplicates it" \
    '[[ "$rc" -eq 0 ]] && cmp -s "$CALLS" "$T/serial.tsv"'

# ---------------------------------------------------------------------------
# C. distributed sub-modes on a FRESH cache: emit-manifest (plan only) -> --chunk 0/1 -> --gather,
#    must reproduce the serial table byte for byte.
# ---------------------------------------------------------------------------
cp "$T/calls.orig.tsv" "$T/w/calls2.tsv"; : > "$SAI_LOG"
rc=0; bash "$STEP" --calls "$T/w/calls2.tsv" --ref "$T/ref.fa" --outdir "$T/w/r2" --distance 4999 \
        --chunk-size 2 --emit-manifest "$T/w/r2/manifest.copy.txt" > "$T/C1.log" 2>&1 || rc=$?
chk "C: --emit-manifest plans 2 chunks and scores nothing" \
    '[[ "$rc" -eq 0 && "$(grep -c . "$T/w/r2/manifest.copy.txt")" -eq 2 && ! -s "$SAI_LOG" ]]'
chk "C: the chunk VCFs exist, carry contig lines from the .fai and one record per pending variant" \
    'c0="$(sed -n 1p "$T/w/r2/manifest.txt")"; [[ -f "$c0" ]] && grep -q "^##contig=<ID=chr1,length=100000>" "$c0" && [[ "$(grep -cv "^#" "$c0")" -eq 2 ]]'
# gather BEFORE scoring: must refuse (a partial merge is the failure nobody notices)
rc=0; bash "$STEP" --calls "$T/w/calls2.tsv" --ref "$T/ref.fa" --outdir "$T/w/r2" --distance 4999 --gather > "$T/C2.log" 2>&1 || rc=$?
chk "C: gather with unscored chunks fails loudly and leaves the calls table untouched" \
    '[[ "$rc" -ne 0 ]] && grep -q "PARTIAL" "$T/C2.log" && cmp -s "$T/w/calls2.tsv" "$T/calls.orig.tsv"'
for i in 0 1; do
    rc=0; bash "$STEP" --calls "$T/w/calls2.tsv" --ref "$T/ref.fa" --outdir "$T/w/r2" --distance 4999 --chunk "$i" > "$T/C3.$i.log" 2>&1 || rc=$?
    chk "C: --chunk $i scores exactly one chunk" '[[ "$rc" -eq 0 ]]'
done
chk "C: two array tasks -> two model calls" '[[ "$(wc -l < "$SAI_LOG" | tr -d " ")" -eq 2 ]]'
rc=0; bash "$STEP" --calls "$T/w/calls2.tsv" --ref "$T/ref.fa" --outdir "$T/w/r2" --distance 4999 --chunk 0 > "$T/C4.log" 2>&1 || rc=$?
chk "C: re-running a chunk is a cached no-op (a requeued array element re-scores nothing)" \
    '[[ "$rc" -eq 0 && "$(wc -l < "$SAI_LOG" | tr -d " ")" -eq 2 ]]'
rc=0; bash "$STEP" --calls "$T/w/calls2.tsv" --ref "$T/ref.fa" --outdir "$T/w/r2" --distance 4999 --gather > "$T/C5.log" 2>&1 || rc=$?
chk "C: gather after scoring exits 0" '[[ "$rc" -eq 0 ]]'
chk "C: distributed output is BYTE-IDENTICAL to the serial run" 'cmp -s "$T/w/calls2.tsv" "$T/serial.tsv"'

# ---------------------------------------------------------------------------
# D. a changed window invalidates the cache (a 50 bp score and a 4999 bp score are different
#    quantities); the cache header records the new window.
# ---------------------------------------------------------------------------
: > "$SAI_LOG"
rc=0; bash "$STEP" --calls "$CALLS" --ref "$T/ref.fa" --outdir "$T/w/spliceai_rescore" \
        --distance 500 --chunk-size 10 --threads 1 > "$T/D.log" 2>&1 || rc=$?
chk "D: a different -D re-scores every variant (1 chunk of 3) and re-stamps the cache" \
    '[[ "$rc" -eq 0 && "$(wc -l < "$SAI_LOG" | tr -d " ")" -eq 1 ]] && grep -q "distance=500" "$T/w/spliceai_rescore/scores.tsv"'

# ---------------------------------------------------------------------------
# E. the isolated env is missing -> exit 3 (UNAVAILABLE), the caller's halt contract.
# ---------------------------------------------------------------------------
rc=0; HPRV_SPLICEAI_ENV="$T/no_such_env" bash "$STEP" --calls "$CALLS" --ref "$T/ref.fa" \
        --outdir "$T/w/r3" --distance 4999 > "$T/E.log" 2>&1 || rc=$?
chk "E: a missing spliceai env exits 3 with an actionable message" \
    '[[ "$rc" -eq 3 ]] && grep -q "spliceai_rescore.enabled: false" "$T/E.log"'

[[ "$fail" -eq 0 ]] && echo "ALL STEP-5B RESCORE ASSERTIONS PASSED" || { echo "STEP-5B RESCORE TEST FAILED"; exit 1; }
