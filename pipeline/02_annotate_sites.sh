#!/usr/bin/env bash
# =============================================================================
# 02_annotate_sites.sh  —  Pipeline Step 2: annotate the cohort site list ONCE
#
# Annotate the (deduplicated) cohort site-only union with the full stack, so the
# expensive VEP pass runs once over distinct sites rather than per trio. See
# docs/functional_annotation.md and docs/pipeline_design.md (Step 2).
#
# Pipeline: VEP (cache + plugins) -> bcftools
#           +split-vep (lift key fields to INFO) -> transfer gnomAD faf95/AF and
#           ClinVar CLNSIG/CLNREVSTAT from external sites VCFs.
#
# Resource paths + tag names come from the environment (exported by the config;
# see src/hprv/config.py). Each plugin/annotation is added ONLY if its resource is
# configured; a missing REQUIRED resource (VEP cache) fails loudly, optional ones
# warn and are skipped so the pipeline still runs with whatever you have.
#
# Usage:
#   02_annotate_sites.sh --sites cohort.sites.vcf.gz --ref GRCh38.fa \
#       --out cohort.sites.annotated.vcf.gz [--threads N]
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"

SITES="" REF="${HPRV_REF_FASTA:-}" OUT="" THREADS="${HPRV_THREADS:-4}"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sites)   SITES="$2"; shift 2;;
        --ref)     REF="$2"; shift 2;;
        --out)     OUT="$2"; shift 2;;
        --threads) THREADS="$2"; shift 2;;
        --tmpdir)  HPRV_TMPDIR="$2"; shift 2;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done
[[ -n "$SITES" && -n "$REF" && -n "$OUT" ]] || die "need --sites, --ref, and --out"

# A value is "set" if non-empty and not a leftover ${...} placeholder.
is_set() { [[ -n "${1:-}" && "$1" != *'${'* ]]; }
require_path() { is_set "$1" || die "$2 is not configured"; [[ -e "$1" ]] || die "$2 not found: $1"; }

require_path "${HPRV_VEP_CACHE:-}" "VEP cache (resources.vep.cache_dir)"
VEP_VERSION="${HPRV_VEP_VERSION:-115}"

if is_done "$OUT"; then log "Step 2 already complete: $OUT (skipping)"; exit 0; fi

outdir="$(abspath_dir "$OUT")"; mkdir -p "$outdir"
# Bind every resource dir so tool calls see them inside/outside the container.
binds="$outdir $(abspath_dir "$SITES") $(abspath_dir "$REF") $HPRV_VEP_CACHE"
for r in "${HPRV_VEP_PLUGINS:-}" "${HPRV_CADD_SNV:-}" "${HPRV_CADD_INDEL:-}" \
         "${HPRV_DBNSFP:-}" "${HPRV_SPLICEAI_SNV:-}" "${HPRV_SPLICEAI_INDEL:-}" \
         "${HPRV_LOFTEE_DATA:-}" "${HPRV_GNOMAD_SITES:-}" "${HPRV_CLINVAR_VCF:-}"; do
    is_set "$r" && [[ -e "$r" ]] && binds+=" $(abspath_dir "$r")"
done
HPRV_BIND="$(printf '%s\n' $binds | sort -u | tr '\n' ' ')"; export HPRV_BIND

vep_vcf="$HPRV_TMPDIR/vep.vcf.gz"; mkdir -p "$HPRV_TMPDIR"

# --- assemble the VEP command (targeted; MANE-first single pick) ---
vep_args=(
    --cache --offline --dir_cache "$HPRV_VEP_CACHE" --cache_version "$VEP_VERSION"
    --species homo_sapiens --assembly GRCh38 --fasta "$REF"
    --vcf --compress_output bgzip --force_overwrite --no_stats
    --symbol --biotype --numbers --hgvs --canonical --mane
    --pick --pick_order mane_select,mane_plus_clinical,canonical,rank
    --sift b --polyphen b
    --fork "$THREADS"
)
is_set "${HPRV_VEP_PLUGINS:-}" && vep_args+=(--dir_plugins "$HPRV_VEP_PLUGINS")

# CADD via the dedicated plugin = the COMPLETE CADD source: whole_genome_SNVs scores every
# possible SNV genome-wide (coding AND non-coding) + the precomputed indel set. This is the
# authoritative CADD field (vep_CADD_PHRED); dbNSFP's CADD_phred (coding-only) is intentionally
# NOT requested to avoid a partial duplicate. Matches the group's existing annotate setup.
if is_set "${HPRV_CADD_SNV:-}" && is_set "${HPRV_CADD_INDEL:-}"; then
    vep_args+=(--plugin "CADD,snv=${HPRV_CADD_SNV},indels=${HPRV_CADD_INDEL}")
else warn "CADD plugin not configured — no genome-wide CADD (dbNSFP CADD_phred is coding-only and not requested)"; fi

# dbNSFP (REVEL / AlphaMissense / MPC / MetaRNN — the calibrated missense predictors). CADD comes
# from the dedicated plugin above (more complete), so CADD_phred is deliberately omitted here.
if is_set "${HPRV_DBNSFP:-}"; then
    vep_args+=(--plugin "dbNSFP,${HPRV_DBNSFP},REVEL_score,AlphaMissense_score,AlphaMissense_pred,MPC_score,MetaRNN_score")
else warn "dbNSFP not configured — REVEL/AlphaMissense/MPC will be unavailable"; fi

# SpliceAI (masked scores recommended for interpretation)
if is_set "${HPRV_SPLICEAI_SNV:-}" && is_set "${HPRV_SPLICEAI_INDEL:-}"; then
    vep_args+=(--plugin "SpliceAI,snv=${HPRV_SPLICEAI_SNV},indel=${HPRV_SPLICEAI_INDEL}")
else warn "SpliceAI not configured — splice predictions will be unavailable"; fi

# LOFTEE (GRCh38 fork). Data-file layout varies by build; allow a full override via
# HPRV_LOF_PLUGIN, else construct a best-effort default the user should verify.
if is_set "${HPRV_LOF_PLUGIN:-}"; then
    vep_args+=(--plugin "$HPRV_LOF_PLUGIN")
elif is_set "${HPRV_LOFTEE_DATA:-}" && is_set "${HPRV_VEP_PLUGINS:-}"; then
    vep_args+=(--plugin "LoF,loftee_path:${HPRV_VEP_PLUGINS},human_ancestor_fa:${HPRV_LOFTEE_DATA}/human_ancestor.fa.gz,conservation_file:${HPRV_LOFTEE_DATA}/loftee.sql,gerp_bigwig:${HPRV_LOFTEE_DATA}/gerp_conservation_scores.homo_sapiens.GRCh38.bw")
    warn "Using default LOFTEE data-file names under HPRV_LOFTEE_DATA — verify they exist, or set HPRV_LOF_PLUGIN to override the whole plugin string"
else warn "LOFTEE not configured — pLoF HC/LC confidence will be unavailable"; fi

n_sites="$(count_variants "$SITES")"
# VEP runs EXACTLY ONCE here, on the deduplicated cohort union — never per trio.
# Per-trio steps (Step 4) transfer these annotations with `bcftools annotate`, they
# do not re-run VEP. This is the single most expensive operation in the pipeline.
log "Step 2: VEP-annotating the cohort union ONCE — $n_sites sites (VEP r${VEP_VERSION}). VEP is NOT run per trio."
audit 02_annotate input_sites "$n_sites"
hprv_run -- vep "${vep_args[@]}" -i "$SITES" -o "$vep_vcf"
require_intact_bgzip "$vep_vcf"

# --- split-vep: lift only the CSQ fields that are actually present to INFO ---
# Capture the header ONCE, then parse in bash: piping bcftools into `grep -m1` lets grep
# close the pipe early and SIGPIPE-abort bcftools under `set -o pipefail`.
vep_header="$(hprv_run -- bcftools view -h "$vep_vcf")"
csq_line=""
while IFS= read -r _hl; do
    case "$_hl" in '##INFO=<ID=CSQ'*) csq_line="$_hl"; break ;; esac
done <<< "$vep_header"
[[ -n "$csq_line" ]] || die "VEP output has no CSQ header — annotation failed"
csq_fmt="${csq_line##*Format: }"; csq_fmt="${csq_fmt%%\">*}"
[[ -n "$csq_fmt" ]] || die "VEP CSQ header has no Format description"
want="Consequence IMPACT SYMBOL Gene Feature BIOTYPE HGVSc HGVSp MANE_SELECT CANONICAL \
      REVEL_score AlphaMissense_score AlphaMissense_pred MPC_score MetaRNN_score CADD_PHRED \
      SpliceAI_pred_DS_AG SpliceAI_pred_DS_AL SpliceAI_pred_DS_DG SpliceAI_pred_DS_DL \
      LoF LoF_filter LoF_flags"
have_fields=""
for w in $want; do
    [[ "|$csq_fmt|" == *"|$w|"* ]] && have_fields+="${have_fields:+,}$w"
done
[[ -n "$have_fields" ]] || die "none of the desired CSQ fields are present"
log "Step 2: split-vep lifting fields: $have_fields"

# core VEP fields must always be present; their absence means the annotation itself is broken
_have() { case ",$have_fields," in *",$1,"*) return 0;; *) return 1;; esac; }
for core in Consequence IMPACT SYMBOL; do
    _have "$core" || die "VEP CSQ is missing core field '$core' — annotation is broken (bad cache/plugins?)"
done
# functional predictors: warn LOUDLY (not fatal — a mode can proceed on other evidence) when a
# CONFIGURED plugin's field did not land, so a misconfigured dbNSFP/SpliceAI/CADD/LOFTEE is caught
# rather than silently dropping that evidence on the whole cohort.
if is_set "${HPRV_DBNSFP:-}"      && ! _have REVEL_score;         then warn "dbNSFP configured but no vep_REVEL_score lifted — check the dbNSFP plugin columns"; fi
if is_set "${HPRV_SPLICEAI_SNV:-}" && ! _have SpliceAI_pred_DS_AG; then warn "SpliceAI configured but no vep_SpliceAI_pred_DS lifted — check the SpliceAI plugin"; fi
if is_set "${HPRV_CADD_SNV:-}"    && ! _have CADD_PHRED;          then warn "CADD configured but no vep_CADD_PHRED lifted — check the CADD plugin"; fi
if { is_set "${HPRV_LOFTEE_DATA:-}" || is_set "${HPRV_LOF_PLUGIN:-}"; } && ! _have LoF; then warn "LOFTEE configured but no vep_LoF lifted — check the LoF plugin / data files"; fi

split_vcf="$HPRV_TMPDIR/split.vcf.gz"
hprv_run -- bcftools +split-vep -c "$have_fields" -s worst -p vep_ \
    --threads "$THREADS" -Oz -o "$split_vcf" "$vep_vcf"
index_vcf "$split_vcf"

# --- transfer external gnomAD (faf95/AF/nhomalt) and ClinVar (CLNSIG/CLNREVSTAT) ---
# gnomAD faf95 is the SOLE population-frequency oracle (golden rule #2): a configured-but-missing
# path is a hard error, never a silent skip. Only a genuinely unset resource warns-and-continues.
cur="$split_vcf"
if is_set "${HPRV_GNOMAD_SITES:-}"; then
    [[ -e "$HPRV_GNOMAD_SITES" ]] || die "gnomAD sites configured but not found: $HPRV_GNOMAD_SITES (resources.gnomad.sites_vcf)"
    log "Step 2: transferring gnomAD frequency (faf95 = ${HPRV_GNOMAD_FAF95_TAG:-fafmax_faf95_max_joint})"
    gn="$HPRV_TMPDIR/gnomad.vcf.gz"
    hprv_run -- bcftools annotate -a "$HPRV_GNOMAD_SITES" \
        -c "INFO/hprv_gnomad_af:=INFO/${HPRV_GNOMAD_AF_TAG:-AF_joint},INFO/hprv_gnomad_grpmax_af:=INFO/${HPRV_GNOMAD_GRPMAX_AF_TAG:-AF_grpmax_joint},INFO/hprv_gnomad_faf95:=INFO/${HPRV_GNOMAD_FAF95_TAG:-fafmax_faf95_max_joint},INFO/hprv_gnomad_nhomalt:=INFO/${HPRV_GNOMAD_NHOMALT_TAG:-nhomalt_joint}" \
        --threads "$THREADS" -Oz -o "$gn" "$cur"
    index_vcf "$gn"; cur="$gn"
    # sanity: a contig-naming mismatch (chr1 vs 1) or wrong tag names annotates 0 records but exits 0
    n_gn="$(hprv_run -- bcftools view -H -i 'INFO/hprv_gnomad_af!="." || INFO/hprv_gnomad_faf95!="."' "$gn" | wc -l | tr -d '[:space:]')"
    log "Step 2: gnomAD annotated $n_gn / $n_sites sites"
    [[ "$n_gn" -gt 0 || "$n_sites" -eq 0 ]] || \
        die "gnomAD transfer matched 0/$n_sites sites — check contig naming (chr1 vs 1) and the resources.gnomad.*_tag names"
else warn "gnomAD sites not configured (resources.gnomad.sites_vcf) — rarity filtering will be unavailable downstream"; fi

if is_set "${HPRV_CLINVAR_VCF:-}"; then
    [[ -e "$HPRV_CLINVAR_VCF" ]] || die "ClinVar VCF configured but not found: $HPRV_CLINVAR_VCF (resources.clinvar.vcf)"
    log "Step 2: transferring ClinVar (${HPRV_CLINVAR_SIG_TAG:-CLNSIG})"
    cv="$HPRV_TMPDIR/clinvar.vcf.gz"
    cln_cols="INFO/hprv_clnsig:=INFO/${HPRV_CLINVAR_SIG_TAG:-CLNSIG},INFO/hprv_clnrevstat:=INFO/${HPRV_CLINVAR_REVSTAT_TAG:-CLNREVSTAT}"
    # CLNSIGCONF is absent in some ClinVar builds; requesting a missing source tag aborts the whole
    # transfer (losing CLNSIG too), so only add it when present in the source header (bash match, no pipe).
    sigconf_tag="${HPRV_CLINVAR_SIGCONF_TAG:-CLNSIGCONF}"
    cln_header="$(hprv_run -- bcftools view -h "$HPRV_CLINVAR_VCF")"
    if [[ "$cln_header" == *"##INFO=<ID=${sigconf_tag},"* ]]; then
        cln_cols+=",INFO/hprv_clnsigconf:=INFO/${sigconf_tag}"
    else warn "ClinVar source lacks ${sigconf_tag}; transferring CLNSIG/CLNREVSTAT only"; fi
    hprv_run -- bcftools annotate -a "$HPRV_CLINVAR_VCF" -c "$cln_cols" \
        --threads "$THREADS" -Oz -o "$cv" "$cur"
    index_vcf "$cv"; cur="$cv"
else warn "ClinVar not configured — clinical evidence will be unavailable downstream"; fi

cp "$cur" "$OUT"; index_vcf "$OUT"
require_intact_bgzip "$OUT"; mark_done "$OUT"
audit 02_annotate annotated_sites "$(count_variants "$OUT")"
log "Step 2 complete: $OUT"
