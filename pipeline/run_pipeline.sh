#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh  —  end-to-end orchestrator (resolve + Steps 0-8)
#
# Runs the whole screen from a single config. Designed to run INSIDE the container
# (tools + python env native on PATH):
#
#   apptainer exec --cleanenv \
#       --bind "$REF_DIR" --bind "$VEP_CACHE" --bind "$WORK" --bind "$DATA" \
#       hprv.sif run_pipeline.sh --config config/config.yaml
#
# Inputs (from config): a kid/dad/mom trios file + a VCF source (dir and/or list).
# The resolver maps each trio to the VCF containing all three members and generates
# the internal manifest + PEDs. Each step is idempotent (.done files) and records
# counts to the run audit; a summary is assembled at the end.
# See docs/pipeline_design.md.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"
export PYTHONPATH="${PYTHONPATH:-}:${HPRV_HOME:-$(cd "$HERE/.." && pwd)}/src"

CFG="" FROM=0 TO=8
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CFG="$2"; shift 2;;
        --from)   FROM="$2"; shift 2;;
        --to)     TO="$2"; shift 2;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done
[[ -n "$CFG" && -f "$CFG" ]] || die "need --config <config.yaml>"

is_set() { [[ -n "${1:-}" && "$1" != *'${'* ]]; }
cfg_get() { python3 -m hprv.config get --config "$CFG" --key "$1" --default "${2:-}"; }
run_step() { local n="$1"; [[ "$FROM" -le "$n" && "$n" -le "$TO" ]]; }

# Resolve config -> environment (HPRV_* vars). Warnings for unset placeholders -> stderr.
eval "$(python3 -m hprv.config sh --config "$CFG")"

is_set "${HPRV_OUTPUT_DIR:-}"  || die "project.output_dir is unresolved — set the env var it references"
is_set "${HPRV_REF_FASTA:-}"   || die "reference.fasta is unresolved — set the env var it references"
is_set "${HPRV_TRIOS_FILE:-}"  || die "inputs.trios_file is unresolved — set the env var it references"
[[ -f "$HPRV_TRIOS_FILE" ]]    || die "trios file not found: $HPRV_TRIOS_FILE"
is_set "${HPRV_VCF_DIR:-}" || is_set "${HPRV_VCF_LIST:-}" || die "set inputs.vcf_dir and/or inputs.vcf_list"

W="$HPRV_OUTPUT_DIR"; mkdir -p "$W"
export HPRV_TMPDIR="${HPRV_TMPDIR:-$W/tmp}"; mkdir -p "$HPRV_TMPDIR"
export HPRV_AUDIT_DIR="$W/audit"; mkdir -p "$HPRV_AUDIT_DIR"
RESOLVED="$W/trios.resolved.tsv"

# ---------------------------------------------------------------------------
# Preflight: resolve trios -> VCFs + PEDs. Runs when starting from <=1 or when
# the resolved manifest is missing (needed by steps 0/1/4).
# ---------------------------------------------------------------------------
if [[ "$FROM" -le 1 || ! -f "$RESOLVED" ]]; then
    log "== Resolve: mapping trios to VCFs + generating PEDs =="
    rargs=(--trios "$HPRV_TRIOS_FILE" --outdir "$W")
    is_set "${HPRV_VCF_DIR:-}"  && [[ -d "$HPRV_VCF_DIR" ]]  && rargs+=(--vcf-dir "$HPRV_VCF_DIR")
    is_set "${HPRV_VCF_LIST:-}" && [[ -f "$HPRV_VCF_LIST" ]] && rargs+=(--vcf-list "$HPRV_VCF_LIST")
    python3 "$HERE/resolve_trios.py" "${rargs[@]}"
fi
[[ -f "$RESOLVED" ]] || die "resolved manifest missing: $RESOLVED"
n_trios=$(($(grep -cve '^[[:space:]]*$' "$RESOLVED") - 1))
log "Pipeline: $n_trios resolved trios, steps $FROM..$TO, work=$W"

if run_step 0; then
    log "== Step 0: per-trio QC =="
    python3 "$HERE/00_qc.py" --manifest "$RESOLVED" --config "$CFG" --out "$W/qc_report.tsv"
fi

if run_step 1; then
    log "== Step 1: cohort site union =="
    bash "$HERE/01_make_cohort_sites.sh" --manifest "$RESOLVED" \
        --ref "$HPRV_REF_FASTA" --out "$W/cohort.sites.vcf.gz"
fi

if run_step 2; then
    log "== Step 2: annotate sites (VEP once) =="
    bash "$HERE/02_annotate_sites.sh" --sites "$W/cohort.sites.vcf.gz" \
        --ref "$HPRV_REF_FASTA" --out "$W/cohort.sites.annotated.vcf.gz"
fi

if run_step 3; then
    log "== Step 3: select plausible =="
    python3 "$HERE/03_select_plausible.py" --in "$W/cohort.sites.annotated.vcf.gz" \
        --out "$W/plausible.sites.vcf.gz" --config "$CFG"
fi

if run_step 4; then
    log "== Step 4: per-trio subset + annotate =="
    bash "$HERE/04_subset_and_annotate_trios.sh" --manifest "$RESOLVED" \
        --plausible "$W/plausible.sites.vcf.gz" --ref "$HPRV_REF_FASTA" --outdir "$W"
fi

if run_step 5; then
    log "== Step 5: inheritance screen =="
    python3 "$HERE/05_inheritance_screen.py" --manifest "$W/trios.candidates.tsv" \
        --config "$CFG" --out "$W/candidates.calls.tsv" --qc-report "$W/qc_report.tsv"
fi

if run_step 6; then
    log "== Step 6: cross-pedigree gene burden =="
    mut="$(cfg_get resources.mutation_rate_table)"
    con="$(cfg_get resources.constraint.gnomad_v2_constraint)"
    extra=()
    is_set "$mut" && [[ -e "$mut" ]] && extra+=(--mutrate "$mut")
    is_set "$con" && [[ -e "$con" ]] && extra+=(--constraint "$con")
    python3 "$HERE/06_gene_burden.py" --calls "$W/candidates.calls.tsv" \
        --out "$W/genes.ranked.tsv" --config "$CFG" --n-trios "$n_trios" "${extra[@]}"
fi

if run_step 7 && [[ "$(cfg_get outputs.xlsx true)" != "false" ]]; then
    log "== Step 7: consolidated xlsx summary =="
    python3 "$HERE/07_report_xlsx.py" --work "$W" --config "$CFG" \
        --out "$W/hprv_summary.xlsx" --label "$(basename "$W")"
fi

if run_step 8 && [[ "$(cfg_get outputs.igv.enabled true)" != "false" ]]; then
    log "== Step 8: igv.js variant-review export =="
    pad="$(cfg_get outputs.igv.padding 1000)"
    gen="$(cfg_get outputs.igv.genome hg38)"
    ig=(--work "$W" --ref "$HPRV_REF_FASTA" --padding "$pad" --genome "$gen")
    cm="$(cfg_get resources.cram_map)"
    is_set "$cm" && [[ -f "$cm" ]] && ig+=(--cram-map "$cm")
    bash "$HERE/08_igv_export.sh" "${ig[@]}"
fi

# Assemble the run audit summary (what went where, and why).
python3 -m hprv.audit --dir "$HPRV_AUDIT_DIR" --out "$HPRV_AUDIT_DIR/summary.md" >/dev/null || true

log "Pipeline complete. Key outputs in $W:"
log "  trios.resolved.tsv  trio_resolution.tsv  qc_report.tsv"
log "  candidates.calls.tsv  genes.ranked.tsv  hprv_summary.xlsx"
log "  igv/variants.tsv (+ crams/ vcfs/ trios.tsv curation.json)"
log "  audit/summary.md  audit/counts.tsv"
