#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh  —  end-to-end orchestrator (Steps 0-6)
#
# Runs the whole screen from a single config. Designed to run INSIDE the container
# (tools + python env native on PATH):
#
#   apptainer exec --cleanenv \
#       --bind "$REF_DIR" --bind "$VEP_CACHE" --bind "$WORK" --bind "$DATA" \
#       hprv.sif run_pipeline.sh --config config/config.yaml
#
# Each step is idempotent (.done files), so a re-run resumes where it stopped.
# See docs/pipeline_design.md for the flow and artifacts.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"
export PYTHONPATH="${PYTHONPATH:-}:${HPRV_HOME:-$(cd "$HERE/.." && pwd)}/src"

CFG="" FROM=0 TO=6
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

# Resolve config -> environment (HPRV_* vars). Warnings for unset placeholders go to stderr.
eval "$(python3 -m hprv.config sh --config "$CFG")"

is_set "${HPRV_OUTPUT_DIR:-}"   || die "project.output_dir is unresolved — set the env var it references"
is_set "${HPRV_REF_FASTA:-}"    || die "reference.fasta is unresolved — set the env var it references"
is_set "${HPRV_TRIO_MANIFEST:-}" || die "inputs.trio_manifest is unresolved — set the env var it references"
[[ -f "$HPRV_TRIO_MANIFEST" ]]  || die "trio manifest not found: $HPRV_TRIO_MANIFEST"

W="$HPRV_OUTPUT_DIR"; mkdir -p "$W"
export HPRV_TMPDIR="${HPRV_TMPDIR:-$W/tmp}"; mkdir -p "$HPRV_TMPDIR"
n_trios=$(($(grep -cve '^[[:space:]]*$' "$HPRV_TRIO_MANIFEST") - 1))
log "Pipeline start: $n_trios trios, steps $FROM..$TO, work=$W"

if run_step 0; then
    log "== Step 0: per-trio QC =="
    python3 "$HERE/00_qc.py" --manifest "$HPRV_TRIO_MANIFEST" --config "$CFG" --out "$W/qc_report.tsv"
fi

if run_step 1; then
    log "== Step 1: cohort site union =="
    bash "$HERE/01_make_cohort_sites.sh" --manifest "$HPRV_TRIO_MANIFEST" \
        --ref "$HPRV_REF_FASTA" --out "$W/cohort.sites.vcf.gz"
fi

if run_step 2; then
    log "== Step 2: annotate sites =="
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
    bash "$HERE/04_subset_and_annotate_trios.sh" --manifest "$HPRV_TRIO_MANIFEST" \
        --plausible "$W/plausible.sites.vcf.gz" --ref "$HPRV_REF_FASTA" --outdir "$W"
fi

if run_step 5; then
    log "== Step 5: inheritance screen =="
    python3 "$HERE/05_inheritance_screen.py" --manifest "$W/trios.candidates.tsv" \
        --config "$CFG" --out "$W/candidates.calls.tsv"
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

log "Pipeline complete. Key outputs in $W:"
log "  qc_report.tsv  candidates.calls.tsv  genes.ranked.tsv"
