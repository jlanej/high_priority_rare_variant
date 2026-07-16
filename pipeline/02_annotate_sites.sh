#!/usr/bin/env bash
# =============================================================================
# 02_annotate_sites.sh  —  Pipeline Step 2: annotate the cohort site list ONCE
#
# Annotate the (deduplicated) cohort site-only union, so the expensive VEP pass runs
# once over distinct sites rather than per trio. See docs/functional_annotation.md
# and docs/pipeline_design.md (Step 2).
#
# Pipeline: VEP (cache + CADD plugin) -> bcftools +split-vep (lift CSQ -> INFO).
#
# **VEP-only contract.** Everything downstream reads comes out of the VEP cache: the
# gnomAD v4.1 per-population AFs (--af_gnomade/--af_gnomadg), ClinVar CLIN_SIG
# (--check_existing), and CADD from its plugin. NOTHING is transferred in from an
# external sites VCF — no gnomAD, ClinVar, dbNSFP, SpliceAI or LOFTEE download exists.
# The cost is real and deliberate; see docs/allele_frequency.md for the ledger:
#   - no faf95 (the cache has no AC/AN, so the CI correction is not reconstructible)
#   - no nhomalt, no SpliceAI, no LOFTEE, no ClinVar review status/stars
# Re-adding any of them = one bcftools annotate here + its INFO field in annotations.F.
#
# Already have a VEP VCF? Pass --vep-vcf (or set resources.vep.annotated_vcf) and the
# VEP call is skipped entirely; the file is split-vep'd as-is. It must be VEP 115
# GRCh38 and carry the fields above — this script verifies that rather than trusting it.
#
# Usage:
#   02_annotate_sites.sh --sites cohort.sites.vcf.gz --ref GRCh38.fa \
#       --out cohort.sites.annotated.vcf.gz [--vep-vcf pre_annotated.vcf.gz] [--threads N]
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"

SITES="" REF="${HPRV_REF_FASTA:-}" OUT="" THREADS="${HPRV_THREADS:-4}"
PRE_VEP="${HPRV_VEP_ANNOTATED_VCF:-}"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sites)   SITES="$2"; shift 2;;
        --ref)     REF="$2"; shift 2;;
        --out)     OUT="$2"; shift 2;;
        --vep-vcf) PRE_VEP="$2"; shift 2;;
        --threads) THREADS="$2"; shift 2;;
        --tmpdir)  HPRV_TMPDIR="$2"; shift 2;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done

# A value is "set" if non-empty and not a leftover ${...} placeholder.
is_set() { [[ -n "${1:-}" && "$1" != *'${'* ]]; }
require_path() { is_set "$1" || die "$2 is not configured"; [[ -e "$1" ]] || die "$2 not found: $1"; }

# Two modes: annotate $SITES with VEP ourselves, or ingest a VEP VCF someone else made.
# The cache/ref/sites are only needed in the former.
if is_set "$PRE_VEP"; then
    [[ -e "$PRE_VEP" ]] || die "pre-annotated VEP VCF not found: $PRE_VEP (resources.vep.annotated_vcf)"
    [[ -n "$OUT" ]] || die "need --out"
else
    [[ -n "$SITES" && -n "$REF" && -n "$OUT" ]] || die "need --sites, --ref, and --out (or --vep-vcf)"
    require_path "${HPRV_VEP_CACHE:-}" "VEP cache (resources.vep.cache_dir)"
fi
VEP_VERSION="${HPRV_VEP_VERSION:-115}"

if is_done "$OUT"; then log "Step 2 already complete: $OUT (skipping)"; exit 0; fi

outdir="$(abspath_dir "$OUT")"; mkdir -p "$outdir"
# Bind every resource dir so tool calls see them inside/outside the container.
binds="$outdir"
for r in "$SITES" "$REF" "$PRE_VEP" "${HPRV_VEP_CACHE:-}" "${HPRV_VEP_PLUGINS:-}" \
         "${HPRV_CADD_SNV:-}" "${HPRV_CADD_INDEL:-}"; do
    is_set "$r" && [[ -e "$r" ]] && binds+=" $(abspath_dir "$r")"
done
HPRV_BIND="$(printf '%s\n' $binds | sort -u | tr '\n' ' ')"; export HPRV_BIND

vep_vcf="$HPRV_TMPDIR/vep.vcf.gz"; mkdir -p "$HPRV_TMPDIR"

if is_set "$PRE_VEP"; then
    # Ingest mode: someone else already ran VEP. Verify it is the right build rather than
    # trusting the filename — a GRCh37 or older-release VCF would annotate "successfully"
    # and be wrong everywhere downstream.
    n_sites="$(count_variants "$PRE_VEP")"
    log "Step 2: ingesting pre-annotated VEP VCF (VEP is NOT run) — $n_sites sites: $PRE_VEP"
    pre_header="$(hprv_run -- bcftools view -h "$PRE_VEP")"
    case "$pre_header" in
        *'##VEP='*) ;;
        *) die "--vep-vcf has no ##VEP header line — it is not a VEP-annotated VCF: $PRE_VEP";;
    esac
    case "$pre_header" in
        *"assembly=\"GRCh38"*|*'assembly=GRCh38'*) ;;
        *) warn "--vep-vcf does not declare assembly=GRCh38 — verify the build; every coordinate downstream assumes GRCh38";;
    esac
    case "$pre_header" in
        *"##VEP=\"v${VEP_VERSION}"*) ;;
        *) warn "--vep-vcf was not made by VEP v${VEP_VERSION} (resources.vep.version) — transcript models and cached frequencies may differ from what the config documents";;
    esac
    vep_vcf="$PRE_VEP"
    audit 02_annotate input_sites "$n_sites"
else
    # --- assemble the VEP command ---
    # The cache supplies gnomAD v4.1 per-population AFs and ClinVar CLIN_SIG; CADD is the one
    # plugin. --flag_pick (not --pick) keeps EVERY consequence block and merely marks the chosen
    # one PICK=1, which lets split-vep below choose the selection rule — and keeps this output
    # shaped identically to an externally-produced --flag_pick VCF, so both modes take one path.
    vep_args=(
        --cache --offline --dir_cache "$HPRV_VEP_CACHE" --cache_version "$VEP_VERSION"
        --species homo_sapiens --assembly GRCh38 --fasta "$REF"
        --vcf --compress_output bgzip --force_overwrite --no_stats
        --symbol --biotype --numbers --hgvs --canonical --mane
        # The rarity oracle and the clinical evidence, straight from the cache. Without
        # --af_gnomade/--af_gnomadg there is NO population frequency anywhere in this pipeline;
        # without --check_existing there is no CLIN_SIG. Both are load-bearing, not extras.
        --af_gnomade --af_gnomadg --max_af --check_existing
        --flag_pick --pick_order mane_select,mane_plus_clinical,canonical,rank
        --fork "$THREADS"
    )
    is_set "${HPRV_VEP_PLUGINS:-}" && vep_args+=(--dir_plugins "$HPRV_VEP_PLUGINS")

    # CADD via the dedicated plugin: whole_genome_SNVs scores every possible SNV genome-wide
    # (coding AND non-coding) + the precomputed indel set. Under this contract it is the ONLY
    # functional predictor, and therefore the only keep-path for anything below MODERATE impact.
    if is_set "${HPRV_CADD_SNV:-}" && is_set "${HPRV_CADD_INDEL:-}"; then
        vep_args+=(--plugin "CADD,snv=${HPRV_CADD_SNV},indels=${HPRV_CADD_INDEL}")
    else warn "CADD plugin not configured — with no CADD there is NO functional evidence for any variant VEP rates below MODERATE (intronic/synonymous/UTR/regulatory); the screen degrades to an impact-only filter"; fi

    n_sites="$(count_variants "$SITES")"
    # VEP runs EXACTLY ONCE here, on the deduplicated cohort union — never per trio.
    # Per-trio steps (Step 4) transfer these annotations with `bcftools annotate`, they
    # do not re-run VEP. This is the single most expensive operation in the pipeline.
    log "Step 2: VEP-annotating the cohort union ONCE — $n_sites sites (VEP r${VEP_VERSION}). VEP is NOT run per trio."
    audit 02_annotate input_sites "$n_sites"
    hprv_run -- vep "${vep_args[@]}" -i "$SITES" -o "$vep_vcf"
    require_intact_bgzip "$vep_vcf"
fi

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
# The grpmax-ELIGIBLE gnomAD populations (annotations.GRPMAX_POPS) — the rarity oracle.
# Deliberately NOT the full population set and NOT MAX_AF: see the frequency check below.
GRPMAX_AF_FIELDS="gnomADe_AFR_AF gnomADe_AMR_AF gnomADe_EAS_AF gnomADe_NFE_AF gnomADe_SAS_AF \
                  gnomADg_AFR_AF gnomADg_AMR_AF gnomADg_EAS_AF gnomADg_NFE_AF gnomADg_SAS_AF"
want="Consequence IMPACT SYMBOL Gene Feature BIOTYPE HGVSc HGVSp MANE_SELECT \
      CADD_PHRED CLIN_SIG gnomADe_AF gnomADg_AF MAX_AF MAX_AF_POPS $GRPMAX_AF_FIELDS"
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

# The grpmax-eligible AF fields ARE the rarity oracle now (annotations.frequency()). If none of
# them landed, every rarity gate reads None => "not in gnomAD" => rarest, and the screen silently
# keeps everything and calls BA1-common polymorphisms candidates. That is the same class of
# silent-catastrophe the old gnomAD-transfer 0-match guard caught, so it dies the same way.
_n_grpmax=0
for p in $GRPMAX_AF_FIELDS; do _have "$p" && _n_grpmax=$((_n_grpmax + 1)); done
[[ "$_n_grpmax" -gt 0 ]] || die "VEP CSQ carries none of the gnomAD grpmax AF fields ($GRPMAX_AF_FIELDS) — \
population frequency would be silently absent and NOTHING would be filtered as common. Re-run VEP with \
--af_gnomade --af_gnomadg (or --everything)."
log "Step 2: rarity oracle = $_n_grpmax/10 gnomAD grpmax-eligible AF fields"
if ! _have CADD_PHRED; then
    warn "no vep_CADD_PHRED lifted — CADD is the ONLY functional predictor under this contract, so \
without it nothing below MODERATE impact can ever be kept (all intronic/synonymous/UTR evidence is gone)"
fi
if ! _have CLIN_SIG; then warn "no vep_CLIN_SIG lifted — the ClinVar P/LP override is inactive (re-run VEP with --check_existing)"; fi

# CSQ selection. VEP --pick emits ONE block; --flag_pick emits ALL blocks with the chosen one
# marked PICK=1. Under --flag_pick, `-s worst` would silently select the worst consequence ACROSS
# transcripts instead of the picked one — a different annotation than the producer intended. So
# honor PICK when it exists and fall back to worst only when it does not. Override with
# HPRV_CSQ_SELECT (bcftools accepts: all, worst, primary, mane, pick, or an EXPRESSION).
# Caveat worth knowing: a --pick_order beginning with `rank` picks the WORST-consequence
# transcript, so SYMBOL can name a non-MANE/readthrough gene, which Step 6 then aggregates
# carriers under. Set HPRV_CSQ_SELECT=mane if gene attribution matters more than recall.
if is_set "${HPRV_CSQ_SELECT:-}"; then sel="$HPRV_CSQ_SELECT"
elif [[ "|$csq_fmt|" == *"|PICK|"* ]]; then sel="pick"
else sel="worst"; fi

# `pick` / `mane` / the EXPRESSION form were ALL added to split-vep in bcftools 1.20 (2024-04-15,
# commit 944fc93); 1.10-1.19 accept only all/worst/primary and die with "the transcript selection
# key ... is not recognised" — AFTER VEP has already run, which on real data is hours wasted. So
# check up front. Gated on the VERSION, not the usage text: an earlier revision of this guard
# scraped the selector list out of `--help` and was wrong twice over — the wording is
# "TR, transcript:" through 1.21.2 and "TR, filter transcripts:" only from 1.22, so it read empty
# on precisely the old versions it was meant to catch (failing open), while 1.20/1.21 support the
# selectors but print the old wording (failing closed). A release number is the durable fact.
# Relevant if you run steps on a host rather than in the image (CLAUDE.md "Dev/host"): Ubuntu
# jammy ships bcftools 1.13 and noble 1.19, so `apt-get install bcftools` CANNOT run this path.
# There is no faithful pre-1.20 workaround: `-i 'PICK=1'` filters whole sites, not CSQ blocks, and
# `-c PICK -d -i 'PICK=1'` drops sites with no PICK'd transcript and duplicates multiallelics.
# NOT auto-downgraded to `worst`: on a --flag_pick VCF that silently swaps "the transcript VEP
# chose" for "the worst consequence across all transcripts" — a different annotation, and exactly
# the kind of quiet science change this pipeline exists to prevent.
case "$sel" in
    all|worst|primary) ;;   # available in every bcftools that ships split-vep (since 1.10)
    *)
        _bcf_v="$(hprv_run -- bcftools --version 2>/dev/null | head -1 | awk '{print $2}')"
        _bcf_maj="${_bcf_v%%.*}"; _bcf_rest="${_bcf_v#*.}"; _bcf_min="${_bcf_rest%%.*}"
        if [[ ! "$_bcf_maj" =~ ^[0-9]+$ || ! "$_bcf_min" =~ ^[0-9]+$ ]]; then
            # Unparseable version (odd build): warn rather than block. bcftools' own error is
            # survivable; a false die here would be worse than the message we are improving on.
            warn "could not parse a bcftools version from '${_bcf_v:-?}' to pre-check '-s $sel'; continuing (bcftools will error if unsupported)"
        elif (( _bcf_maj < 1 || (_bcf_maj == 1 && _bcf_min < 20) )); then
            die "bcftools $_bcf_v is too old for the '$sel' transcript selector — +split-vep gained pick/mane in \
1.20, and 1.10-1.19 offer only all/worst/primary. The image pins a bcftools that supports it, so run inside the \
image or upgrade (note Ubuntu's apt bcftools is 1.13 on jammy / 1.19 on noble). You may instead set \
resources.vep.csq_select to all|worst|primary — but do NOT reach for 'worst' on a --flag_pick VCF unless you mean \
it: it takes the worst consequence ACROSS transcripts rather than the one VEP picked, changing every downstream \
SYMBOL/IMPACT."
        fi
        ;;
esac
log "Step 2: split-vep transcript selection: -s $sel"

split_vcf="$HPRV_TMPDIR/split.vcf.gz"
# NB: no --threads. `bcftools +split-vep` is a plugin and does NOT accept it ("unrecognized
# option `--threads'"); passing it aborts the step. $THREADS applies to vep --fork above.
hprv_run -- bcftools +split-vep -c "$have_fields" -s "$sel" -p vep_ \
    -Oz -o "$split_vcf" "$vep_vcf"
index_vcf "$split_vcf"

# --- no external transfers: gnomAD/ClinVar came from the cache with the CSQ above ---
# Presence of the FIELD only proves VEP emitted the column, not that the cache actually had
# frequencies to put in it. A cache built without the frequency data, or an input whose alleles
# are all un-accessioned, yields a fully-populated header over entirely empty values — and every
# rarity gate then reads None ("rarest") and keeps everything. So assert on the VALUES.
cur="$split_vcf"
freq_expr='INFO/vep_gnomADe_AF!="." || INFO/vep_gnomADg_AF!="." || INFO/vep_MAX_AF!="."'
n_freq="$(hprv_run -- bcftools view -H -i "$freq_expr" "$cur" | wc -l | tr -d '[:space:]')"
log "Step 2: gnomAD frequency present on $n_freq / $n_sites sites"
[[ "$n_freq" -gt 0 || "$n_sites" -eq 0 ]] || \
    die "0/$n_sites sites carry ANY gnomAD frequency — the cache lookup is not working. A cohort union \
always contains some known alleles, so this is a broken cache/build, not a very rare cohort. Rarity \
filtering would silently pass everything. Check --af_gnomade/--af_gnomadg and the cache version."

cp "$cur" "$OUT"; index_vcf "$OUT"
require_intact_bgzip "$OUT"; mark_done "$OUT"
audit 02_annotate annotated_sites "$(count_variants "$OUT")"
log "Step 2 complete: $OUT"
