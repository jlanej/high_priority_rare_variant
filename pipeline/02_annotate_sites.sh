#!/usr/bin/env bash
# =============================================================================
# 02_annotate_sites.sh  —  Pipeline Step 2: annotate the cohort site list ONCE
#
# Annotate the (deduplicated) cohort site-only union, so the expensive VEP pass runs
# once over distinct sites rather than per trio. See docs/functional_annotation.md
# and docs/pipeline_design.md (Step 2).
#
# Pipeline: VEP (cache + CADD and SpliceAI plugins) -> bcftools +split-vep (lift CSQ -> INFO).
#
# **VEP-centric contract.** Almost everything downstream reads comes out of the VEP cache: the
# gnomAD v4.1 per-population AFs (--af_gnomade/--af_gnomadg), ClinVar CLIN_SIG
# (--check_existing), and CADD, SpliceAI, REVEL + AlphaMissense from their plugins. A plugin
# score file is a PLUGIN input, not a bcftools transfer, so those do not breach the contract.
#
# There are exactly TWO bcftools transfers, both because the cache cannot supply the field at any
# price, both prefixed with their own namespace so the different oracle is visible at a glance,
# and both guarded by a 0-match check below:
#   - clinvar_* : ClinVar review status / GOLD STARS (CLNREVSTAT is absent from the cache).
#   - gnomad_*  : faf95 + nhomalt from the gnomAD v4.1 JOINT slim (OPTIONAL). faf95's CI
#                 correction needs AC/AN, which the cache omits. Without it the rarity oracle
#                 falls back to the grpmax point-estimate proxy — see annotations.frequency().
#
# The remaining cost is real and deliberate; see docs/allele_frequency.md for the ledger:
#   - no LOFTEE (pLoF confidence), no exome/genome discordance flag
# Re-adding either = one bcftools annotate here + its INFO field in annotations.F.
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
# Distributed (SLURM) sub-modes. Default (none set) = the self-contained in-process run.
#   --emit-shard-manifest F : write the ordered contig list one-per-line to F, then exit. A SLURM
#                             job array indexes into it ($SLURM_ARRAY_TASK_ID -> line).
#   --shard-contig CHR      : VEP-annotate ONLY contig CHR -> a shard VCF + .done, then exit. One
#                             array task per contig. Skips split-vep / gather / output.
#   --gather                : skip VEP; verify every expected shard is complete, concat them, then
#                             run split-vep + the guards + write $OUT. The dependent gather job.
# See docs/pipeline_design.md (Step 2, distributed) and pipeline/slurm/.
EMIT_MANIFEST="" SHARD_CONTIG="" GATHER=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sites)   SITES="$2"; shift 2;;
        --ref)     REF="$2"; shift 2;;
        --out)     OUT="$2"; shift 2;;
        --vep-vcf) PRE_VEP="$2"; shift 2;;
        --threads) THREADS="$2"; shift 2;;
        --tmpdir)  HPRV_TMPDIR="$2"; shift 2;;
        --emit-shard-manifest) EMIT_MANIFEST="$2"; shift 2;;
        --shard-contig)        SHARD_CONTIG="$2"; shift 2;;
        --gather)              GATHER=1; shift;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done

# A value is "set" if non-empty and not a leftover ${...} placeholder.
is_set() { [[ -n "${1:-}" && "$1" != *'${'* ]]; }
require_path() { is_set "$1" || die "$2 is not configured"; [[ -e "$1" ]] || die "$2 not found: $1"; }

# The distributed sub-modes are mutually exclusive and never combine with --vep-vcf (ingest).
_nmode=0
is_set "$EMIT_MANIFEST" && _nmode=$((_nmode + 1)); is_set "$SHARD_CONTIG" && _nmode=$((_nmode + 1))
[[ "$GATHER" -eq 1 ]] && _nmode=$((_nmode + 1))
[[ "$_nmode" -le 1 ]] || die "--emit-shard-manifest / --shard-contig / --gather are mutually exclusive"
{ [[ "$_nmode" -eq 1 ]] && is_set "$PRE_VEP"; } && die "the distributed sub-modes run VEP; they cannot combine with --vep-vcf"

# Two run styles: annotate $SITES with VEP ourselves, or ingest a VEP VCF someone else made — the
# cache/ref are only needed in the former. --emit-shard-manifest additionally needs neither.
if is_set "$PRE_VEP"; then
    [[ -e "$PRE_VEP" ]] || die "pre-annotated VEP VCF not found: $PRE_VEP (resources.vep.annotated_vcf)"
    [[ -n "$OUT" ]] || die "need --out"
else
    [[ -n "$SITES" && -n "$OUT" ]] || die "need --sites and --out (or --vep-vcf)"
    if ! is_set "$EMIT_MANIFEST"; then
        [[ -n "$REF" ]] || die "need --ref"
        require_path "${HPRV_VEP_CACHE:-}" "VEP cache (resources.vep.cache_dir)"
    fi
fi
VEP_VERSION="${HPRV_VEP_VERSION:-115}"

# The $OUT-complete short-circuit applies only to run styles that PRODUCE $OUT (default, ingest,
# gather). --shard-contig produces a per-contig shard and --emit-shard-manifest produces a manifest.
if ! is_set "$SHARD_CONTIG" && ! is_set "$EMIT_MANIFEST"; then
    # Keyed to the INPUT that produced $OUT, so a re-unioned cohort (Step 1 rebuilt because the trio
    # set changed) correctly invalidates the annotated union instead of being masked by a bare
    # existence marker. Without this, Step 1's keyed invalidation would not propagate: Step 2 would
    # skip and Step 3 would select from stale annotations.
    # The key source differs by run style: the sites union normally, the ingested VEP VCF under
    # --vep-vcf (where --sites is not required at all, so $SITES may be empty).
    _kin="$SITES"; is_set "$PRE_VEP" && _kin="$PRE_VEP"
    _skey=""; [[ -f "$_kin" ]] && _skey="$(cksum < "$_kin" | awk '{print $1"-"$2}')"
    # The ANNOTATION RESOURCES are part of the key, not just the input sites. A new ClinVar
    # release reclassifies variants and a newly-supplied REVEL/AlphaMissense file adds columns —
    # both change this output while leaving the sites union byte-identical, so a sites-only key
    # reports "already complete" and serves the old annotation forever. (Exactly the staleness
    # Step 9's content key was rewritten to close.) Identity is path+size+mtime, NOT a content
    # hash: these are 0.2-80 GB of static, version-pinned reference data whose bytes never change
    # in place, and cksum'ing the ~80 GB of CADD on every invocation to detect a swap that shows
    # up in the stat anyway would dominate the step's startup.
    if [[ -n "$_skey" ]]; then
        # EVERY annotation resource, not just the new ones. SpliceAI matters MOST here and was
        # the omission that made this a bug rather than a nicety: it is a SELECTION keep-path
        # (selection.py checks it before CADD), so a run made with spliceai_required:false and
        # re-run after downloading the scores would hit `is_done`, match the key, and screen
        # forever against a union that has no SpliceAI — silently, since the "configured but
        # nothing lifted" warning lives inside the skipped path. CADD has the same shape and
        # needs no non-default config at all.
        for _res in "${HPRV_CLINVAR_VCF:-}" "${HPRV_GNOMAD_SITES:-}" \
                    "${HPRV_REVEL:-}" "${HPRV_ALPHAMISSENSE:-}" \
                    "${HPRV_SPLICEAI_SNV:-}" "${HPRV_SPLICEAI_INDEL:-}" \
                    "${HPRV_CADD_SNV:-}" "${HPRV_CADD_INDEL:-}"; do
            if is_set "$_res" && [[ -e "$_res" ]]; then
                _skey+="-$(cksum <<<"$_res$(stat -c '%s-%Y' "$_res" 2>/dev/null \
                          || stat -f '%z-%m' "$_res" 2>/dev/null)" | awk '{print $1}')"
            else
                _skey+="-0"   # absent is itself a state: supplying the file later must invalidate
            fi
        done
    fi
    # An unreadable input leaves the key empty -> never skip; recomputing is the safe direction
    # (the missing input then fails loudly below rather than silently reusing a stale annotation).
    if [[ -n "$_skey" ]] && is_done "$OUT" && [[ "$(cat "$OUT.done" 2>/dev/null)" == "$_skey" ]]; then
        log "Step 2 already complete: $OUT (skipping)"; exit 0
    fi
fi

outdir="$(abspath_dir "$OUT")"; mkdir -p "$outdir"
# Bind every resource dir so tool calls see them inside/outside the container.
binds="$outdir"
for r in "$SITES" "$REF" "$PRE_VEP" "${HPRV_VEP_CACHE:-}" "${HPRV_VEP_PLUGINS:-}" \
         "${HPRV_CADD_SNV:-}" "${HPRV_CADD_INDEL:-}" \
         "${HPRV_SPLICEAI_SNV:-}" "${HPRV_SPLICEAI_INDEL:-}" \
         "${HPRV_REVEL:-}" "${HPRV_ALPHAMISSENSE:-}" "${HPRV_CLINVAR_VCF:-}" \
         "${HPRV_GNOMAD_SITES:-}"; do
    is_set "$r" && [[ -e "$r" ]] && binds+=" $(abspath_dir "$r")"
done
HPRV_BIND="$(printf '%s\n' $binds | sort -u | tr '\n' ' ')"; export HPRV_BIND

vep_vcf="$HPRV_TMPDIR/vep.vcf.gz"; mkdir -p "$HPRV_TMPDIR"

# --- shard helpers: shared by the in-process loop AND the SLURM sub-modes -------------------
# The per-contig VEP OUTPUTS land next to $OUT — i.e. on the SHARED filesystem, visible to every
# node — NOT under HPRV_TMPDIR (often node-local scratch). A SLURM gather job runs on a different
# node than the scatter tasks and must read their shards, and resume across job submissions needs
# them to persist. Only the transient region-subset is written to node-local scratch.
shard_dir() { printf '%s/annotate_shards' "$outdir"; }
# Contigs that carry variants, in FILE (= reference/coordinate) order, from the index — no full
# scan. Scatter, gather and the manifest all call this, so they agree exactly on the shard set.
enum_contigs() { hprv_run -- bcftools index -s "$SITES" | awk -F'\t' '($3+0)>0{print $1}'; }
# Annotate exactly one contig -> shard VCF + .done. Idempotent: skips a complete shard (resume).
annotate_one_contig() {  # $1 = contig ; requires vep_args set
    local c="$1" sd; sd="$(shard_dir)"; mkdir -p "$sd"
    local out="$sd/vep.${c}.vcf.gz"
    if is_done "$out"; then log "  [$c] cached (resume)"; return 0; fi
    local sub="${HPRV_TMPDIR:-/tmp}/sites.${c}.$$.vcf.gz"
    # Brace-quote the contig as a whole region: htslib's region parser splits an unquoted
    # `contig:pos` on ':', so a GRCh38 ALT_HLA contig (e.g. `HLA-A*01:01:01:01`, present in the
    # standard Broad reference and carried untouched through Step 1) fails region-init with
    # "Could not parse the region(s)" and aborts the whole (multi-hour) annotation. `{$c}` forces
    # the entire string to be read as a contig name — verified safe + index-scoped for plain
    # contigs (chr1) too. (`-t "$c"` does NOT help — it colon-splits identically.)
    hprv_run -- bcftools view -r "{$c}" -Oz -o "$sub" "$SITES"
    log "  [$c] VEP ($(count_variants "$sub") sites)"
    hprv_run -- vep "${vep_args[@]}" -i "$sub" -o "$out"
    require_intact_bgzip "$out"; mark_done "$out"; rm -f "$sub"
}
# Verify EVERY expected shard is complete, then concat (file order = sorted) into $vep_vcf. The
# verify IS the coherence guarantee: gather never builds a PARTIAL call set even if a scheduler
# dependency misfires. Dies loudly, naming what is missing.
gather_shards() {  # writes $vep_vcf
    local sd; sd="$(shard_dir)"
    local list="$sd/concat.list"; : > "$list"
    local c so missing=0
    while IFS= read -r c; do
        so="$sd/vep.${c}.vcf.gz"
        if is_done "$so"; then printf '%s\n' "$so" >> "$list"
        else warn "gather: shard for contig '$c' is missing/incomplete: $so"; missing=$((missing + 1)); fi
    done < <(enum_contigs)
    [[ -s "$list" ]] || die "gather: no completed shards under $sd — run the scatter first"
    [[ "$missing" -eq 0 ]] || die "gather: $missing contig shard(s) missing/incomplete — the call set \
would be PARTIAL. Re-run the scatter (only contigs without a .done re-run), then gather again."
    log "Step 2: gather — concatenating $(wc -l < "$list" | tr -d ' ') annotated contig shards"
    hprv_run -- bcftools concat -f "$list" -Oz -o "$vep_vcf"
    require_intact_bgzip "$vep_vcf"
}

# --emit-shard-manifest: write the ordered contig list (the array indexes into it) and exit.
if is_set "$EMIT_MANIFEST"; then
    [[ -f "$SITES.tbi" || -f "$SITES.csi" ]] || die "--emit-shard-manifest needs an indexed --sites"
    enum_contigs > "$EMIT_MANIFEST"
    _n="$(wc -l < "$EMIT_MANIFEST" | tr -d ' ')"
    [[ "$_n" -gt 0 ]] || die "no contigs with variants in $SITES"
    log "Step 2: wrote $_n-contig shard manifest -> $EMIT_MANIFEST"
    exit 0
fi

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
    # Build check must FAIL LOUDLY, not warn: a GRCh37 VEP VCF carries a valid ##VEP header and valid
    # vep_gnomAD*_AF, so nothing else stops it — Step 3 would select on GRCh37 coordinates and Step 4
    # would transfer annotations by POSITION onto GRCh38 trios, silently losing/mis-annotating calls.
    # Scope the assembly token to the ##VEP= line so a stray `assembly=` on a ##contig line isn't
    # misread. die only on an EXPLICIT non-GRCh38 assembly; warn (don't fail) when the token is absent
    # — a legitimate GRCh38 VEP VCF may omit it, and both hard guards below still apply.
    vep_line="$(printf '%s\n' "$pre_header" | grep -m1 '^##VEP=' || true)"
    case "$vep_line" in
        *'assembly="GRCh38'*|*'assembly=GRCh38'*) ;;                              # explicitly GRCh38 — good
        *assembly=*) die "--vep-vcf explicitly declares a non-GRCh38 assembly in its ##VEP header, but \
every coordinate downstream (Step 3 selection, Step 4 position-keyed transfer) assumes GRCh38 — ingesting \
it would silently mis-map. Re-annotate against a GRCh38 cache. ##VEP: ${vep_line#*assembly=}" ;;
        *) warn "--vep-vcf's ##VEP header declares no assembly= token — cannot verify the build; every \
coordinate downstream assumes GRCh38. Proceeding on the assumption it is GRCh38." ;;
    esac
    case "$pre_header" in
        *"##VEP=\"v${VEP_VERSION}"*) ;;
        *) warn "--vep-vcf was not made by VEP v${VEP_VERSION} (resources.vep.version) — transcript models and cached frequencies may differ from what the config documents";;
    esac
    vep_vcf="$PRE_VEP"
    audit 02_annotate input_sites "$n_sites"
elif [[ "$GATHER" -eq 1 ]]; then
    # --- gather-only (the SLURM dependent job): no VEP here ---
    # Verify every contig shard the scatter was supposed to produce, concat them into the
    # whole-union VCF, then fall through to split-vep + the guards + output exactly as a single
    # run would. gather_shards() dies loudly if any shard is missing, so this never builds a
    # partial call set even if the scheduler dependency misfired.
    n_sites="$(count_variants "$SITES")"
    audit 02_annotate input_sites "$n_sites"
    gather_shards
else
    # --- assemble the VEP command (shared by --shard-contig and the in-process scatter) ---
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

    # SpliceAI plugin: precomputed raw genome-wide splice delta scores (Illumina), the one signal
    # that reaches deep-intronic cryptic sites + exonic-synonymous splice disruptions VEP's
    # positional terms and CADD both under-call. Emits vep_SpliceAI_pred_DS_* (+ DP_*) subfields.
    # NB the plugin's own `cutoff=` is deliberately NOT passed — we lift the raw delta scores and
    # threshold in selection.py (filters.functional.spliceai_ds_min) so nothing is masked pre-review.
    # Files must be bgzipped + tabix-indexed. Optional + graceful (same as CADD).
    # Check EXISTENCE, not just set-ness: emit-env exports the SpliceAI paths unconditionally, so a
    # configured-but-not-yet-downloaded file must skip the plugin gracefully rather than make VEP
    # abort on an unopenable score file. Both snv= and indel= are required by the plugin.
    if is_set "${HPRV_SPLICEAI_SNV:-}" && [[ -e "$HPRV_SPLICEAI_SNV" ]] \
       && is_set "${HPRV_SPLICEAI_INDEL:-}" && [[ -e "$HPRV_SPLICEAI_INDEL" ]]; then
        vep_args+=(--plugin "SpliceAI,snv=${HPRV_SPLICEAI_SNV},indel=${HPRV_SPLICEAI_INDEL}")
    else warn "SpliceAI plugin inactive (score files unset or not present) — deep-intronic and exonic-synonymous splice-disrupting variants that VEP's positional terms miss will be invisible (CADD is only a weak, lossy proxy for the splice signal). See docs/resources.md to add spliceai_snv/spliceai_indel."; fi

    # --- calibrated MISSENSE predictors (optional + graceful, same contract as CADD) -----------
    # These change the SCREEN by nothing at all — missense is IMPACT=MODERATE and selection.py
    # keeps it at the impact rung before any predictor runs (docs/limitations.md #7). They exist
    # for Step 9's missense TIER, which without them can only report an off-label CADD rank.
    #
    # REVEL emits the CSQ key `REVEL` -> vep_REVEL. The score file must be sorted on the GRCh38
    # column and tabix'd -s 1 -b 3 -e 3; a file sorted on the GRCh37 column indexes without
    # complaint and then matches NOTHING — see scripts/prepare_resources.sh (prep_revel).
    if is_set "${HPRV_REVEL:-}" && [[ -e "$HPRV_REVEL" ]]; then
        vep_args+=(--plugin "REVEL,${HPRV_REVEL}")
    else warn "REVEL plugin inactive — Step 9's missense tier falls back to an OFF-LABEL CADD rank (missense_evidence_source=cadd_offlabel). No effect on the screen. See docs/resources.md#revel."; fi
    # AlphaMissense takes `file=` (unlike REVEL's bare path) and emits `am_pathogenicity`/`am_class`
    # — NOT `AlphaMissense_score`, which is dbNSFP's column name for the same quantity. Copying the
    # dbNSFP name into the want list below yields a plugin that runs and a column that is never
    # populated. transcript_match=1 makes the plugin require the CSQ transcript to match the row's
    # transcript instead of taking any row at the position: AlphaMissense is per-transcript, and
    # without it a variant can be scored against a transcript VEP did not pick.
    if is_set "${HPRV_ALPHAMISSENSE:-}" && [[ -e "$HPRV_ALPHAMISSENSE" ]]; then
        vep_args+=(--plugin "AlphaMissense,file=${HPRV_ALPHAMISSENSE},transcript_match=1")
    else warn "AlphaMissense plugin inactive — no effect on the screen; Step 9's missense tier loses the SVI-endorsed predictor that reaches Strong on constrained genes. See docs/resources.md#alphamissense."; fi

    if is_set "$SHARD_CONTIG"; then
        # --- scatter (one SLURM array task): annotate a single contig -> shard + .done, then stop.
        # No split-vep, no gather, no $OUT — the dependent gather job assembles the whole. Idempotent
        # via the shard .done, so requeuing the array re-runs only the contigs that did not finish.
        log "Step 2: scatter — VEP for contig '$SHARD_CONTIG' only (VEP r${VEP_VERSION})."
        annotate_one_contig "$SHARD_CONTIG"
        exit 0
    fi

    n_sites="$(count_variants "$SITES")"
    # VEP runs over the deduplicated cohort union — once per contig, never per trio. Per-trio
    # steps (Step 4) transfer these annotations with `bcftools annotate`; they do not re-run VEP.
    # This is the single most expensive operation in the pipeline (WGS: ~57M sites, ~a day).
    audit 02_annotate input_sites "$n_sites"

    # --- in-process shard the VEP call BY CONTIG for walltime-resumability ---------------------
    # One --fork run over a WGS union is single-node AND un-resumable — a job killed at the
    # scheduler walltime re-annotates all 57M sites from scratch. Sharding by contig makes each an
    # independent VEP run with its own `.done`, so a killed job RESUMES (finished contigs skipped);
    # for multi-node parallelism drive the SAME shard functions as a SLURM array (--shard-contig /
    # --gather; see pipeline/slurm/). CORRECTNESS: this stays byte-identical to a single run —
    # only the vep call is sharded, split-vep and the guards run once on the reassembled whole.
    # A sharded==single equivalence test guards this (tests/integration).
    case "${HPRV_VEP_SHARD_BY_CONTIG:-1}" in 1|true|yes|on) vshard=1;; *) vshard=0;; esac
    if [[ "$vshard" -eq 1 && ( -f "$SITES.tbi" || -f "$SITES.csi" ) ]]; then
        shard_tags=()
        while IFS= read -r _c; do shard_tags+=("$_c"); done < <(enum_contigs)
        [[ ${#shard_tags[@]} -gt 0 ]] || die "Step 2: no contigs with variants in $SITES"
        log "Step 2: VEP sharded by contig — ${#shard_tags[@]} shards, $n_sites sites (VEP \
r${VEP_VERSION}). Each shard has its own .done, so a walltime kill RESUMES."
        for _c in "${shard_tags[@]}"; do annotate_one_contig "$_c"; done
        gather_shards
    else
        [[ "$vshard" -eq 1 ]] && warn "Step 2: $SITES is unindexed — cannot shard by contig; annotating in one un-resumable pass"
        log "Step 2: VEP-annotating the union in ONE pass — $n_sites sites (VEP r${VEP_VERSION})."
        hprv_run -- vep "${vep_args[@]}" -i "$SITES" -o "$vep_vcf"
        require_intact_bgzip "$vep_vcf"
    fi
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
SPLICEAI_FIELDS="SpliceAI_pred_DS_AG SpliceAI_pred_DS_AL SpliceAI_pred_DS_DG SpliceAI_pred_DS_DL \
                 SpliceAI_pred_DP_AG SpliceAI_pred_DP_AL SpliceAI_pred_DP_DG SpliceAI_pred_DP_DL \
                 SpliceAI_pred_SYMBOL"
# NB `am_pathogenicity`/`am_class` are the AlphaMissense PLUGIN's key names. dbNSFP calls the same
# quantity `AlphaMissense_score`; using that name here would lift nothing, silently.
want="Consequence IMPACT SYMBOL Gene Feature BIOTYPE HGVSc HGVSp MANE_SELECT \
      CADD_PHRED CLIN_SIG gnomADe_AF gnomADg_AF MAX_AF MAX_AF_POPS $GRPMAX_AF_FIELDS $SPLICEAI_FIELDS \
      REVEL am_pathogenicity am_class"
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
    warn "no vep_CADD_PHRED lifted — CADD is a primary functional predictor here, so without it a \
variant below MODERATE impact can only be kept via SpliceAI (intronic/synonymous/UTR CADD evidence is gone)"
fi
# SpliceAI is required by default (enforced in run_pipeline's preflight, not here); warn only if it
# was CONFIGURED (plugin data given) yet no score lifted — that
# means the plugin/data is misconfigured, not merely unused, and the splice keep-path would be silently dead.
if is_set "${HPRV_SPLICEAI_SNV:-}" && is_set "${HPRV_SPLICEAI_INDEL:-}" && ! _have SpliceAI_pred_DS_AG; then
    warn "SpliceAI score files are configured but no vep_SpliceAI_pred_DS_* was lifted — the plugin \
produced no scores (wrong build/filenames, un-indexed .tbi, or a VEP/plugin mismatch). The splice \
keep-path is inactive; verify the raw snv/indel VCFs are GRCh38, bgzipped and tabix-indexed."
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
        # `|| true`: `head -1` can SIGPIPE bcftools (exit 141) which pipefail+set -e would turn
        # into a silent abort. bcftools' version banner is small enough to usually fit the pipe
        # buffer, so this is latent rather than firing — but it is the same trap as
        # scripts/download_spliceai.sh, and the guard below already handles an empty value.
        _bcf_v="$(hprv_run -- bcftools --version 2>/dev/null | head -1 | awk '{print $2}' || true)"
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
# $HPRV_TMPDIR defaults to the PERSISTENT $W/tmp and is never cleaned, while split.vcf.gz is
# REWRITTEN every run. index_vcf() is a no-op when any index exists, so a stale index left by a run
# that died after this point (e.g. at the frequency guard below) would be reused here and then
# installed as the authoritative index of $OUT by the `mv` at the end — with the record count and
# offsets of the PREVIOUS run. Force a fresh index to match the fresh data.
rm -f "$split_vcf".tbi "$split_vcf".csi
index_vcf "$split_vcf"

# --- ClinVar review status: the ONE external transfer -----------------------------------------
# Everything else on this VCF is a CSQ field lifted by split-vep. This is not, and cannot be: the
# VEP cache exposes CLIN_SIG but carries NO CLNREVSTAT, so gold stars are unreachable from the
# cache at any price. Transferred under a `clinvar_` prefix (not `vep_`) so a reader can tell at a
# glance which oracle a field came from.
#
# Optional and graceful: no configured/present ClinVar VCF -> skip with a warning, and
# clinvar_stars reads UNAVAILABLE downstream exactly as it did before this existed.
if is_set "${HPRV_CLINVAR_VCF:-}" && [[ -e "${HPRV_CLINVAR_VCF:-}" ]]; then
    cv_out="$HPRV_TMPDIR/clinvar_annotated.vcf.gz"
    # -c with `:=` renames on transfer. CLNREVSTAT is Number=. and its VALUES contain commas
    # ("criteria_provided,_multiple_submitters,_no_conflicts"), so it arrives as an array; the
    # canonicalising parser in annotations._canon_revstat is what makes that harmless.
    if hprv_run -- bcftools annotate -a "$HPRV_CLINVAR_VCF" \
            -c "INFO/clinvar_CLNREVSTAT:=INFO/CLNREVSTAT,INFO/clinvar_CLNSIG:=INFO/CLNSIG" \
            --threads "$THREADS" -Oz -o "$cv_out" "$split_vcf"; then
        rm -f "$cv_out".tbi "$cv_out".csi
        index_vcf "$cv_out"
        # THE 0-MATCH GUARD. This is the whole reason the transfer is not a one-liner. A ClinVar
        # VCF on the wrong build, with the wrong contig naming (chr1 vs 1), or simply not indexed
        # transfers ZERO records while bcftools exits 0 — leaving a fully-populated header over
        # entirely empty values. Every star then reads UNAVAILABLE and the run looks like a
        # cohort ClinVar happens to know nothing about. A cohort union always overlaps ClinVar
        # somewhere, so 0 is a broken join, not a rare cohort. Same class of silent catastrophe
        # as the frequency guard below, and it dies the same way.
        n_cv="$(hprv_run -- bcftools query -i 'INFO/clinvar_CLNREVSTAT!="."' -f '\n' "$cv_out" \
                | wc -l | tr -d '[:space:]')"
        log "Step 2: ClinVar review status transferred to $n_cv / $n_sites sites"
        if [[ "$n_cv" -eq 0 && "$n_sites" -gt 0 ]]; then
            die "0/$n_sites sites received a ClinVar review status from $HPRV_CLINVAR_VCF — the \
transfer matched NOTHING. A cohort union always overlaps ClinVar somewhere, so this is a broken join, \
not a rare cohort. Check: (a) the ClinVar VCF is GRCh38, (b) its contig naming matches this cohort \
(ClinVar ships bare '1', a GRCh38 cohort union is usually 'chr1' — see docs/resources.md#clinvar), \
(c) it is tabix-indexed. Set resources.clinvar.vcf to '' to run without ClinVar stars."
        fi
        split_vcf="$cv_out"
    else
        warn "ClinVar transfer failed (bcftools annotate returned non-zero) — continuing WITHOUT review status; clinvar_stars will read UNAVAILABLE and Step 9's clinical term stays ungated"
        rm -f "$cv_out" "$cv_out".tbi "$cv_out".csi
    fi
else
    warn "no ClinVar VCF configured (resources.clinvar.vcf) — review status/GOLD STARS unavailable, so a 1-star single-submitter assertion and a 3-star expert-panel one are indistinguishable downstream. See docs/resources.md#clinvar"
fi

# --- gnomAD joint slim: faf95 + nhomalt (the SECOND external transfer) ---------------------
# THE most consequential optional resource. It replaces the rarity oracle's grpmax point-estimate
# PROXY with gnomAD's real faf95 (the 95% CI lower bound), which is the quantity ACMG/ClinGen
# specify for frequency filtering. The cache cannot supply it at any price: the CI correction
# needs AC/AN, which the cache omits. Prefixed `gnomad_`.
#
# EXPECT MORE CANDIDATES, not fewer. faf95 <= the point estimate, so the same cutoffs stop
# discarding low-count alleles whose confidence interval never justified the call.
if is_set "${HPRV_GNOMAD_SITES:-}" && [[ -e "${HPRV_GNOMAD_SITES:-}" ]]; then
    gn_out="$HPRV_TMPDIR/gnomad_annotated.vcf.gz"
    # Renamed on transfer so the INFO names are stable regardless of which gnomAD release the
    # slim came from (v4.1 spells them *_joint; a future release may not).
    _gn_cols="INFO/gnomad_faf95:=INFO/fafmax_faf95_max_joint"
    _gn_cols+=",INFO/gnomad_faf95_group:=INFO/fafmax_faf95_max_gen_anc_joint"
    _gn_cols+=",INFO/gnomad_nhomalt:=INFO/nhomalt_joint"
    _gn_cols+=",INFO/gnomad_AF_joint:=INFO/AF_joint"
    _gn_cols+=",INFO/gnomad_AF_grpmax:=INFO/AF_grpmax_joint"
    if hprv_run -- bcftools annotate -a "$HPRV_GNOMAD_SITES" -c "$_gn_cols" \
            --threads "$THREADS" -Oz -o "$gn_out" "$split_vcf"; then
        rm -f "$gn_out".tbi "$gn_out".csi
        index_vcf "$gn_out"
        # 0-MATCH GUARD. Same class as ClinVar's, and it matters more here because this feeds the
        # rarity gate: a slim on the wrong build or with the wrong contig naming transfers ZERO
        # records and exits 0, every faf95 reads absent, frequency() silently falls back to the
        # proxy for EVERY variant, and the run looks like a successful faf95 run that never was.
        #
        # NB the denominator is deliberately "sites with ANY gnomAD INFO", not "sites with faf95".
        # A faf95 count of 0 is legitimate — gnomAD emits fafmax as MISSING wherever no group's CI
        # lower bound clears 0, which is 80% of a chr22 sample. Guarding on faf95 alone would abort
        # correct runs; guarding on the JOIN proves the join. That same AF_joint field is also what
        # lets annotations.frequency() tell "gnomAD says faf95 is 0" (rarest) from "gnomAD has
        # never seen this" (use the proxy) — it is a witness, not just a reporting column.
        n_gn="$(hprv_run -- bcftools query -i 'INFO/gnomad_AF_joint!="."' -f '\n' "$gn_out" \
                | wc -l | tr -d '[:space:]')"
        n_faf="$(hprv_run -- bcftools query -i 'INFO/gnomad_faf95!="."' -f '\n' "$gn_out" \
                 | wc -l | tr -d '[:space:]')"
        log "Step 2: gnomAD joint matched $n_gn / $n_sites sites (faf95 present on $n_faf; absent faf95 = no group has a confidently non-zero AF, NOT an error)"
        if [[ "$n_gn" -eq 0 && "$n_sites" -gt 0 ]]; then
            die "0/$n_sites sites matched the gnomAD joint slim ($HPRV_GNOMAD_SITES) — the transfer \
matched NOTHING, so faf95 would be absent everywhere and every rarity gate would silently fall back to \
the grpmax proxy while the run looked like a faf95 run. A cohort union always overlaps gnomAD, so this \
is a broken join: check the slim is GRCh38, its contigs match this cohort (chr-prefixed), and it is \
tabix-indexed. Set resources.gnomad.sites_slim to '' to run on the proxy deliberately."
        fi
        split_vcf="$gn_out"
    else
        warn "gnomAD joint transfer failed (bcftools annotate returned non-zero) — continuing on the grpmax PROXY oracle; faf95/nhomalt will be absent"
        rm -f "$gn_out" "$gn_out".tbi "$gn_out".csi
    fi
else
    warn "no gnomAD joint slim configured (resources.gnomad.sites_slim) — rarity uses the grpmax POINT-ESTIMATE proxy, which sits ~one CI-width stringent on low-count alleles (errs toward dropping). No faf95, no nhomalt. See docs/allele_frequency.md"
fi

# --- no other external transfers ---
# Presence of the FIELD only proves VEP emitted the column, not that the cache actually had
# frequencies to put in it. A cache built without the frequency data, or an input whose alleles
# are all un-accessioned, yields a fully-populated header over entirely empty values — and every
# rarity gate then reads None ("rarest") and keeps everything. So assert on the VALUES.
cur="$split_vcf"
# Build the presence expression over ONLY the frequency tags VEP actually emitted. bcftools
# validates EVERY referenced INFO tag at filter-init, regardless of `||` short-circuit — so naming
# a tag the header lacks aborts the query even when a present tag carries a real value. The self-run
# path always emits all three global/max fields, but a `--vep-vcf` ingest made with, say,
# `--af_gnomadg` but no `--max_af` carries valid frequencies yet lacks vep_MAX_AF; hardcoding all
# three would refuse it with a misleading "broken cache" die. The grpmax-eligible fields are
# guaranteed present by the header guard above (_n_grpmax>0), so freq_expr is never empty.
freq_expr=""
for _f in gnomADe_AF gnomADg_AF MAX_AF $GRPMAX_AF_FIELDS; do
    _have "$_f" && freq_expr+="${freq_expr:+ || }INFO/vep_${_f}!=\".\""
done
[[ -n "$freq_expr" ]] || die "no gnomAD frequency fields present to value-check (unreachable after the grpmax header guard)"
# `query -f '\n'`, not `view -H | wc -l`: we want a COUNT, and view -H reconstructs every matching
# record — all 25 vep_* INFO fields of it — just to throw it away at wc. query emits one byte per
# match. Measured on the integration data: 15,521 bytes vs 22 for the identical count of 22; that
# ratio is what scales, and at WGS it is tens of GB of formatting through a pipe for a number.
n_freq="$(hprv_run -- bcftools query -i "$freq_expr" -f '\n' "$cur" | wc -l | tr -d '[:space:]')"
log "Step 2: gnomAD frequency present on $n_freq / $n_sites sites"
[[ "$n_freq" -gt 0 || "$n_sites" -eq 0 ]] || \
    die "0/$n_sites sites carry ANY gnomAD frequency — the cache lookup is not working. A cohort union \
always contains some known alleles, so this is a broken cache/build, not a very rare cohort. Rarity \
filtering would silently pass everything. Check --af_gnomade/--af_gnomadg and the cache version."

# `mv`, not `cp`+re-index: $cur is the split-vep output under HPRV_TMPDIR — the same filesystem as
# $OUT under run_pipeline.sh — and at WGS scale it is tens of GB that are ALREADY indexed. Copying
# and re-indexing is two more full passes over the data for nothing. mv degrades to a copy across
# filesystems, so this is never worse than what it replaces.
mv "$cur" "$OUT"
if   [[ -f "$cur.tbi" ]]; then mv "$cur.tbi" "$OUT.tbi"
elif [[ -f "$cur.csi" ]]; then mv "$cur.csi" "$OUT.csi"
else index_vcf "$OUT"; fi
# Keyed marker (see the resume guard at the top): records WHICH input this annotation came from, so
# a re-unioned cohort invalidates it. Falls back to a bare marker if the key could not be computed.
require_intact_bgzip "$OUT"
if [[ -n "${_skey:-}" ]]; then printf '%s\n' "$_skey" > "$OUT.done"; else mark_done "$OUT"; fi
audit 02_annotate annotated_sites "$(count_variants "$OUT")"
log "Step 2 complete: $OUT"
