#!/usr/bin/env bash
# =============================================================================
# 02b_spliceai_backfill.sh — OPTIONAL: score cohort-union variants that carry NO precomputed
# SpliceAI delta score (mostly novel indels), live, with the stock Illumina model, and merge the
# scores into cohort.sites.annotated.vcf.gz — into the SAME vep_SpliceAI_pred_DS_* fields Step 2's
# SpliceAI plugin produces, so Step 3 selection reads precomputed and backfilled scores identically.
#
# Runs AFTER Step 2 (the plugin has scored everything the precomputed files cover) and BEFORE Step 3
# (so a backfilled score is a keep-path). The model + GENCODE annotation are bundled in the isolated
# `spliceai` conda env baked into the image — NO data download beyond the reference FASTA. TensorFlow
# inference is heavy, so this is OFF by default and gated on the env being present (graceful skip).
#
# Usage:
#   02b_spliceai_backfill.sh --annotated cohort.sites.annotated.vcf.gz --ref GRCh38.fa \
#       [--distance 500] [--indels-only 1] [--tmpdir DIR]
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"
export PYTHONPATH="${HPRV_HOME:-$(cd "$HERE/.." && pwd)}/src${PYTHONPATH:+:$PYTHONPATH}"

ANNOT="" REF="${HPRV_REF_FASTA:-}" DISTANCE=500 INDELS_ONLY=1
SPLICEAI_ENV="${HPRV_SPLICEAI_ENV:-/opt/conda/envs/spliceai}"   # isolated env baked into the image
while [[ $# -gt 0 ]]; do
    case "$1" in
        --annotated) ANNOT="$2"; shift 2;;
        --ref) REF="$2"; shift 2;;
        --distance) DISTANCE="$2"; shift 2;;
        --indels-only) INDELS_ONLY="$2"; shift 2;;
        --tmpdir) HPRV_TMPDIR="$2"; shift 2;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done
[[ -n "$ANNOT" && -f "$ANNOT" ]] || die "need --annotated <cohort.sites.annotated.vcf.gz>"
[[ -n "$REF" && -f "$REF" ]] || die "need --ref <GRCh38 FASTA> (the same reference the cohort was called against)"
[[ "$DISTANCE" =~ ^[0-9]+$ ]] || die "--distance must be a non-negative integer"

# Gate 1: the isolated spliceai env must be present (image-only). On a bare host / CI it is absent —
# warn and skip so the pipeline still runs on the precomputed scores alone.
if ! hprv_run -- test -x "$SPLICEAI_ENV/bin/spliceai"; then
    warn "Step 2b: the isolated 'spliceai' env ($SPLICEAI_ENV) is not present in the runtime — skipping live backfill (precomputed scores only). It ships only in the container image."
    exit 0
fi

# Idempotency: skip if a completed backfill is NEWER than the annotated union (i.e. the union has not
# been regenerated since). We touch the marker AFTER replacing the union, so the marker is the newest.
DONE="$ANNOT.spliceai_backfill.done"
if [[ -f "$DONE" && "$DONE" -nt "$ANNOT" ]]; then
    log "Step 2b: SpliceAI backfill already applied to $ANNOT (rm $DONE to force) — skipping"
    exit 0
fi

mkdir -p "$HPRV_TMPDIR"
sub="$HPRV_TMPDIR/spliceai_backfill.subset.vcf"
scored="$HPRV_TMPDIR/spliceai_backfill.scored.vcf"
annot="$HPRV_TMPDIR/spliceai_backfill.annot.vcf"
# Bind every dir the wrapped tools must see (union, ref, tmp).
HPRV_BIND="$(printf '%s\n' "$(abspath_dir "$ANNOT")" "$(abspath_dir "$REF")" "$HPRV_TMPDIR" $HPRV_BIND | sort -u | tr '\n' ' ')"; export HPRV_BIND

# 1) select variants with NO precomputed SpliceAI score (default: indels only — SNVs are complete).
io=(); [[ "$INDELS_ONLY" != "0" ]] && io=(--indels-only)
python3 -m hprv.spliceai_backfill select --in "$ANNOT" --out "$sub" "${io[@]}"
n_unscored="$(grep -cve '^#' "$sub" 2>/dev/null || true)"; n_unscored=${n_unscored:-0}
if [[ "$n_unscored" -eq 0 ]]; then
    log "Step 2b: no unscored variants to backfill — nothing to do"
    touch "$DONE"; rm -f "$sub"; exit 0
fi
log "Step 2b: $n_unscored variant(s) lack a precomputed SpliceAI score — scoring live (spliceai -D $DISTANCE)"
# Loud heads-up if the set is large: stock SpliceAI does ~1 variant/sec/CPU (no multiprocessing) and
# ~7/sec on GPU, so a big set is slow. (A big set usually means the precomputed plugin was NOT
# configured — configure resources.vep.spliceai_snv/indel first so this only fills the indel gap.)
[[ "$n_unscored" -gt 50000 ]] && warn "Step 2b: $n_unscored unscored variants is a LOT — live SpliceAI is ~1 var/s/CPU (7/s GPU). Did you configure the precomputed resources.vep.spliceai_snv/indel? Backfill is meant to fill only the small novel-indel gap."

# 2) score them live in the isolated env (bundled Illumina model + GENCODE grch38 annotation).
# A backfill failure is NON-fatal: it is an optional enhancement, so warn loudly and leave the union
# on its precomputed scores rather than aborting a long pipeline run.
if ! hprv_run -- micromamba run -n spliceai spliceai -I "$sub" -O "$scored" -R "$REF" -A grch38 -D "$DISTANCE"; then
    warn "Step 2b: spliceai scoring FAILED — the annotated union keeps its PRECOMPUTED scores only (backfill NOT applied). \
Check the reference matches the cohort build/contigs and the spliceai env, then remove $DONE and re-run Step 2 to retry."
    audit 02b_spliceai_backfill scoring_failed 1
    rm -f "$sub" "$scored"
    exit 0
fi

# 3) SpliceAI= -> vep_SpliceAI_pred_DS_* annotation VCF (per-event max over genes), bgzip + index.
python3 -m hprv.spliceai_backfill annotate --in "$scored" --out "$annot"
n_scored="$(grep -cve '^#' "$annot" 2>/dev/null || true)"; n_scored=${n_scored:-0}
if [[ "$n_scored" -eq 0 ]]; then
    log "Step 2b: SpliceAI returned no scores for the unscored set (all below/outside gene windows) — nothing to merge"
    touch "$DONE"; rm -f "$sub" "$scored" "$annot"; exit 0
fi
bgzip -f "$annot"; index_vcf "$annot.gz"

# 4) merge the backfilled deltas onto the union (only the previously-unscored records match), then
#    atomically replace the union and re-index. `-c INFO/<field>` names the exact 4 fields.
merged="$HPRV_TMPDIR/cohort.sites.annotated.backfilled.vcf.gz"
bcftools annotate -a "$annot.gz" \
    -c INFO/vep_SpliceAI_pred_DS_AG,INFO/vep_SpliceAI_pred_DS_AL,INFO/vep_SpliceAI_pred_DS_DG,INFO/vep_SpliceAI_pred_DS_DL \
    -Oz -o "$merged" "$ANNOT"
require_intact_bgzip "$merged"
mv -f "$merged" "$ANNOT"
rm -f "$ANNOT".tbi "$ANNOT".csi
index_vcf "$ANNOT"
touch "$DONE"   # AFTER replacing $ANNOT, so the marker is newer than the union (idempotency check)
audit 02b_spliceai_backfill unscored_selected "$n_unscored"
audit 02b_spliceai_backfill scores_merged "$n_scored"
rm -f "$sub" "$scored" "$annot.gz" "$annot.gz".{tbi,csi} 2>/dev/null || true
log "Step 2b complete: backfilled $n_scored SpliceAI score(s) into $ANNOT (of $n_unscored unscored)"
