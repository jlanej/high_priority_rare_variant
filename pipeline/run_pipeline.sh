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
# prepend our src; avoid a leading ':' (which would put CWD on the import path) when unset
export PYTHONPATH="${HPRV_HOME:-$(cd "$HERE/.." && pwd)}/src${PYTHONPATH:+:$PYTHONPATH}"

CFG="" FROM=0 TO=8
# Distributed Step-2 pass-throughs (used by the SLURM orchestration in pipeline/slurm/). These
# only change how Step 2 runs and are forwarded verbatim to 02_annotate_sites.sh; every other step
# is unaffected. Typically paired with `--from 2 --to 2` so one job does one Step-2 sub-task.
S2_PASSTHRU=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CFG="$2"; shift 2;;
        --from)   FROM="$2"; shift 2;;
        --to)     TO="$2"; shift 2;;
        --annotate-emit-manifest) S2_PASSTHRU+=(--emit-shard-manifest "$2"); shift 2;;
        --annotate-shard-contig)  S2_PASSTHRU+=(--shard-contig "$2"); shift 2;;
        --annotate-gather)        S2_PASSTHRU+=(--gather); shift;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done
[[ -n "$CFG" && -f "$CFG" ]] || die "need --config <config.yaml>"

is_set() { [[ -n "${1:-}" && "$1" != *'${'* ]]; }
cfg_get() { python3 -m hprv.config get --config "$CFG" --key "$1" --default "${2:-}"; }
run_step() { local n="$1"; [[ "$FROM" -le "$n" && "$n" -le "$TO" ]]; }

# Resolve config -> environment (HPRV_* vars). Warnings for unset placeholders -> stderr.
# Capture first so a config parse/emit failure aborts loudly (a bare `eval "$(...)"` under
# set -e swallows the non-zero exit of the substitution).
_cfg_sh="$(python3 -m hprv.config sh --config "$CFG")" || die "failed to resolve config: $CFG"
eval "$_cfg_sh"

is_set "${HPRV_OUTPUT_DIR:-}"  || die "project.output_dir is unresolved — set the env var it references"
is_set "${HPRV_REF_FASTA:-}"   || die "reference.fasta is unresolved — set the env var it references"
# existence, not just placeholder-resolution: a wrong bind-mount must fail here, not mid-Step-1
[[ -f "$HPRV_REF_FASTA" ]]     || die "reference FASTA not found: $HPRV_REF_FASTA (check the bind mount)"
[[ -f "$HPRV_REF_FASTA.fai" ]] || warn "reference index $HPRV_REF_FASTA.fai missing — bcftools norm/samtools will samtools-faidx or fail; run 'samtools faidx $HPRV_REF_FASTA'"
is_set "${HPRV_TRIOS_FILE:-}"  || die "inputs.trios_file is unresolved — set the env var it references"
[[ -f "$HPRV_TRIOS_FILE" ]]    || die "trios file not found: $HPRV_TRIOS_FILE"
is_set "${HPRV_VCF_DIR:-}" || is_set "${HPRV_VCF_LIST:-}" || die "set inputs.vcf_dir and/or inputs.vcf_list"

# ---------------------------------------------------------------------------
# Resource preflight. Under the VEP-only contract the entire surface is: a VEP 115 GRCh38
# cache (which carries gnomAD v4.1 frequencies + ClinVar itself) and the CADD plugin files.
# No gnomAD / ClinVar / dbNSFP / SpliceAI / LOFTEE download exists to check.
# Only enforced when Step 2 actually runs — a `--from 3` re-run reads annotations that are
# already in the VCF and needs none of this. If resources.vep.annotated_vcf is set, VEP is
# not invoked at all, so only that file has to exist.
# ---------------------------------------------------------------------------
if run_step 2; then
    r_missing=()
    _need() { local v="$1"; if ! is_set "${!v:-}" || [[ ! -e "${!v:-}" ]]; then r_missing+=("$2 -> \$$v='${!v:-}'"); fi; }
    _opt()  { local v="$1"; if ! is_set "${!v:-}" || [[ ! -e "${!v:-}" ]]; then warn "resource DEGRADED: $2 missing (\$$v) — that evidence will be unavailable"; fi; }
    if is_set "${HPRV_VEP_ANNOTATED_VCF:-}"; then
        _need HPRV_VEP_ANNOTATED_VCF "pre-annotated VEP VCF (resources.vep.annotated_vcf)"
    else
        _need HPRV_VEP_CACHE "VEP cache (resources.vep.cache_dir) — supplies transcripts, gnomAD v4.1 AFs, and ClinVar"
        # CADD is the ONLY functional predictor here, and the only evidence that can keep a
        # variant VEP rates below MODERATE. Without it the screen is impact-only: every
        # intronic / synonymous / UTR / regulatory candidate is unreachable. Loud, not fatal —
        # an impact-only screen is still a coherent (if narrower) run.
        _opt HPRV_CADD_SNV   "CADD SNV (primary non-coding functional evidence)"
        _opt HPRV_CADD_INDEL "CADD indel (indel-capable functional score)"
        # SpliceAI: optional splice keep-path (deep-intronic / exonic-synonymous). Degrades
        # gracefully — its absence just leaves those splice classes to CADD's weak proxy.
        _opt HPRV_SPLICEAI_SNV   "SpliceAI SNV scores (deep-intronic + synonymous splice detection)"
        _opt HPRV_SPLICEAI_INDEL "SpliceAI indel scores (splice detection for indels)"
    fi
    if [[ ${#r_missing[@]} -gt 0 ]]; then
        warn "Required resources are missing:"
        for m in "${r_missing[@]}"; do warn "  - $m"; done
        die "point resources.vep.cache_dir at a VEP ${HPRV_VEP_VERSION:-115} GRCh38 cache (or set resources.vep.annotated_vcf to a VEP VCF you already have), then re-run. See docs/resources.md."
    fi
fi

W="$HPRV_OUTPUT_DIR"
export HPRV_TMPDIR="${HPRV_TMPDIR:-$W/tmp}"
export HPRV_AUDIT_DIR="$W/audit"
# Fail at second 1, not mid-run: every step writes scratch to the tmpdir and outputs to the work
# dir, and a read-only/unbound path otherwise surfaces hours later as a cryptic tool error (e.g.
# Step 4's `bcftools sort -T` -> "mkdtemp(...) failed: Read-only file system"). A plain
# `mkdir -p` cannot catch it — it succeeds on an existing dir even when the FS is read-only.
require_writable_dir "$W" "project.output_dir"
require_writable_dir "$HPRV_TMPDIR" "runtime.tmpdir (scratch for norm/sort/VEP)"
require_writable_dir "$HPRV_AUDIT_DIR" "audit dir"
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
    # Idempotent like every other step. Step 0 rescans every trio's genotypes (Mendelian errors,
    # chrX sex inference, contamination) — an HOUR on a real cohort — and it was the one step
    # with no .done guard, so it re-ran on every invocation. That is pure waste under the resume
    # workflow: a walltime kill and re-submit, or any SLURM `prep` re-run (--from 0 --to 1),
    # repeated the whole scan even though qc_report.tsv was already complete.
    # Re-run it by removing the sentinel: rm "$W/qc_report.tsv.done"  (do that after changing the
    # trio set or any filters.genotype_qc / qc.* threshold, which this cache cannot detect).
    if is_done "$W/qc_report.tsv"; then
        log "== Step 0: per-trio QC — cached, skipping (rm $W/qc_report.tsv.done to force) =="
    else
        log "== Step 0: per-trio QC =="
        python3 "$HERE/00_qc.py" --manifest "$RESOLVED" --config "$CFG" --out "$W/qc_report.tsv"
        mark_done "$W/qc_report.tsv"
    fi
fi

if run_step 1; then
    log "== Step 1: cohort site union =="
    bash "$HERE/01_make_cohort_sites.sh" --manifest "$RESOLVED" \
        --ref "$HPRV_REF_FASTA" --out "$W/cohort.sites.vcf.gz"
fi

if run_step 2; then
    s2_args=(--sites "$W/cohort.sites.vcf.gz" --ref "$HPRV_REF_FASTA"
             --out "$W/cohort.sites.annotated.vcf.gz")
    if [[ ${#S2_PASSTHRU[@]} -gt 0 ]]; then
        # A distributed Step-2 sub-task (emit-manifest / shard-contig / gather) from the SLURM
        # orchestration. Not compatible with the ingest bypass (which produces no shards).
        log "== Step 2: distributed sub-task (${S2_PASSTHRU[*]}) =="
        s2_args+=("${S2_PASSTHRU[@]}")
    elif is_set "${HPRV_VEP_ANNOTATED_VCF:-}"; then
        log "== Step 2: ingest pre-annotated VEP VCF (VEP not run) =="
        s2_args+=(--vep-vcf "$HPRV_VEP_ANNOTATED_VCF")
    else
        log "== Step 2: annotate sites (VEP once) =="
    fi
    bash "$HERE/02_annotate_sites.sh" "${s2_args[@]}"
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
    # CRAMs are reference-compressed: decoding needs the reference they were ENCODED
    # against, which may differ from the variant-calling reference.fasta (e.g. an
    # alignment build with different decoy/HLA/PAR-masking). Use resources.cram_ref
    # when set, else fall back to reference.fasta.
    cref="${HPRV_CRAM_REF:-}"; is_set "$cref" && [[ -f "$cref" ]] || cref="$HPRV_REF_FASTA"
    # samtools threads per slice. Slicing is serial (one CRAM at a time) so a flaky
    # FUSE/SBFS mount stays healthy; this only sets per-slice compression threads.
    jobs="$(cfg_get outputs.igv.extract_jobs "$(cfg_get runtime.threads 4)")"
    ig=(--work "$W" --ref "$cref" --padding "$pad" --genome "$gen" --jobs "$jobs")
    cm="$(cfg_get resources.cram_map)"
    is_set "$cm" && [[ -f "$cm" ]] && ig+=(--cram-map "$cm")
    # Step 8b (non-human-fraction). Default ON, but activates only when a kraken2 DB is provided;
    # otherwise warn here and leave 8b off (its columns stay blank). run_pipeline is the single
    # enable gate — 08 never reads $HPRV_KRAKEN2_DB on its own.
    if [[ "$(cfg_get outputs.igv.nonhuman_screen.enabled true)" != "false" ]]; then
        kdb="${HPRV_KRAKEN2_DB:-}"
        if is_set "$kdb" && [[ -d "$kdb" ]]; then
            ig+=(--kraken2-db "$kdb"
                 --nhf-members    "$(cfg_get outputs.igv.nonhuman_screen.members carriers)"
                 --nhf-confidence "$(cfg_get outputs.igv.nonhuman_screen.confidence 0.05)"
                 --nhf-min-reads  "$(cfg_get outputs.igv.nonhuman_screen.min_reads 5)")
            [[ "$(cfg_get outputs.igv.nonhuman_screen.memory_mapping true)" != "false" ]] \
                && ig+=(--nhf-memory-mapping)
        else
            warn "outputs.igv.nonhuman_screen.enabled but resources.kraken2_db is unset/not a dir ('$kdb') — skipping NHF annotation (variants.tsv NHF columns will be blank)"
        fi
    fi
    bash "$HERE/08_igv_export.sh" "${ig[@]}"
fi

# Assemble the run audit summary (what went where, and why).
python3 -m hprv.audit --dir "$HPRV_AUDIT_DIR" --out "$HPRV_AUDIT_DIR/summary.md" >/dev/null || true

log "Pipeline complete. Key outputs in $W:"
log "  trios.resolved.tsv  trio_resolution.tsv  qc_report.tsv"
log "  candidates.calls.tsv  genes.ranked.tsv  hprv_summary.xlsx"
log "  igv/variants.tsv (+ crams/ vcfs/ trios.tsv curation.json)"
log "  audit/summary.md  audit/counts.tsv"
