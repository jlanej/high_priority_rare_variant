#!/usr/bin/env bash
# =============================================================================
# prepare_resources.sh — download & prepare the GRCh38 annotation resources the
# hprv pipeline bind-mounts at runtime. The container ships the SOFTWARE; this
# fetches the DATA. Nothing is installed — only downloaded, verified, indexed,
# and wired to the config's ${ENV} placeholders.
#
# Run it INSIDE the image so bcftools/tabix/vep are on PATH, e.g.
#   apptainer exec --bind /data hprv.sif \
#       scripts/prepare_resources.sh --dir /data/hprv_resources fetch
#
# Modes:
#   fetch     (default) download+prepare free resources; validate gated ones you provided
#   verify    check every expected target exists, is indexed, and gnomAD tags match config
#   emit-env  print `export VAR=...` lines for the config ${ENV} placeholders
#
# Flags: --dir DIR (required)  --only id,id  --accept-license  --out FILE (emit-env)
# See docs/resources.md. Pinned URLs/versions/checksums live in resources/manifest.env.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${HPRV_RESOURCE_MANIFEST:-$HERE/../resources/manifest.env}"

DIR="" MODE="fetch" ONLY="" ACCEPT_LICENSE=0 EMIT_OUT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dir)            DIR="$2"; shift 2;;
        --only)           ONLY="$2"; shift 2;;
        --accept-license) ACCEPT_LICENSE=1; shift;;
        --out)            EMIT_OUT="$2"; shift 2;;
        fetch|verify|emit-env) MODE="$1"; shift;;
        -h|--help)        grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) echo "unknown arg: $1" >&2; exit 2;;
    esac
done

log()  { printf '[resources] %s\n' "$*" >&2; }
warn() { printf '[resources] WARNING: %s\n' "$*" >&2; }
die()  { printf '[resources] ERROR: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1; }

[[ -n "$DIR" ]] || die "need --dir DIR (where resources are prepared)"
[[ -f "$MANIFEST" ]] || die "manifest not found: $MANIFEST (set HPRV_RESOURCE_MANIFEST)"
# shellcheck disable=SC1090
source "$MANIFEST"
mkdir -p "$DIR"
DIR="$(cd "$DIR" && pwd)"

# resource id -> the prepared path under $DIR (also what emit-env exports)
REF_FASTA_OUT="$DIR/reference/GRCh38.fa"
VEP_CACHE_OUT="$DIR/vep_cache"
# NB: VEP plugin CODE (CADD/dbNSFP/SpliceAI + LOFTEE grch38) ships IN the image at /plugins —
# it is not fetched here. VEP_PLUGINS defaults to /plugins via the image ENV; not emitted below.
GNOMAD_OUT="$DIR/gnomad/gnomad.joint.v${GNOMAD_VERSION}.sites.slim.vcf.gz"
CLINVAR_OUT="$DIR/clinvar/clinvar_${CLINVAR_DATE}.GRCh38.vcf.gz"
DBNSFP_OUT="$DIR/dbnsfp/${DBNSFP_EXPECT}"
CADD_SNV_OUT="$DIR/cadd/whole_genome_SNVs.tsv.gz"
CADD_INDEL_OUT="$DIR/cadd/gnomad.genomes.r4.0.indel.tsv.gz"
SPLICEAI_SNV_OUT="$DIR/spliceai/${SPLICEAI_SNV_EXPECT}"
SPLICEAI_INDEL_OUT="$DIR/spliceai/${SPLICEAI_INDEL_EXPECT}"
LOFTEE_OUT="$DIR/loftee"
CONSTRAINT_OUT="$DIR/constraint/constraint.by_gene.tsv"
MUTRATE_OUT="$DIR/constraint/mutation_rate.by_gene.tsv"

MISSING=() PREPARED=() SKIPPED=()

selected() {  # honor --only id,id; empty ONLY = all
    [[ -z "$ONLY" ]] && return 0
    case ",$ONLY," in *",$1,"*) return 0;; *) return 1;; esac
}

# -------- helpers ----------------------------------------------------------- #
sha_ok() {  # sha_ok FILE SHA256  (empty/REPLACE_* sha -> skip check, warn once)
    local f="$1" want="$2" got
    [[ -f "$f" ]] || return 1
    case "$want" in ""|REPLACE_*) return 0;; esac
    if need sha256sum; then got="$(sha256sum "$f" | awk '{print $1}')"
    elif need shasum;  then got="$(shasum -a 256 "$f" | awk '{print $1}')"
    else warn "no sha256 tool; cannot verify $f"; return 0; fi
    [[ "$got" == "$want" ]] || { warn "checksum mismatch for $f (got $got want $want)"; return 1; }
}

fetch_to() {  # fetch_to URL OUT  — dispatch by protocol; resumable; retries
    local url="$1" out="$2"
    case "$url" in ""|REPLACE_*) return 3;; esac      # not pinned yet
    mkdir -p "$(dirname "$out")"
    case "$url" in
        gs://*)    need gsutil || return 4; gsutil -q cp "$url" "$out";;
        s3://*)    need aws || return 4; aws s3 cp --no-sign-request "$url" "$out";;
        http://*|https://*|ftp://*)
            if need curl; then curl -fSL --retry 5 --retry-delay 5 -C - -o "$out" "$url"
            elif need wget; then wget -c -O "$out" "$url"
            else return 4; fi;;
        *) return 5;;
    esac
}

index_vcf() { need tabix && tabix -f -p vcf "$1" 2>/dev/null || warn "could not tabix $1"; }

record() { case "$1" in ok) PREPARED+=("$2");; skip) SKIPPED+=("$2");; miss) MISSING+=("$2");; esac; }

# A free download+prepare step, idempotent + checksum-guarded.
# get_free ID URL OUT SHA  -> 0 prepared/cached, 1 not-pinned/failed
get_free() {
    local id="$1" url="$2" out="$3" sha="${4:-}"
    if [[ -f "$out" ]] && sha_ok "$out" "$sha"; then log "[$id] cached"; record skip "$id"; return 0; fi
    log "[$id] downloading -> $out"
    # Capture fetch_to's REAL exit code. Inside `if ! fetch_to ...; then`, `$?` is the status of
    # the negated compound (always 0), so the per-code diagnostics were unreachable and every
    # failure printed the generic "download failed". `if ...; else rc=$?` reads the true code.
    local rc=0
    if fetch_to "$url" "$out"; then :; else rc=$?; fi
    if (( rc )); then
        case $rc in 3) warn "[$id] URL not pinned in manifest.env — skipping";;
                    4) warn "[$id] required downloader (curl/wget/gsutil/aws) missing";;
                    *) warn "[$id] download failed";; esac
        record miss "$id"; return 1
    fi
    sha_ok "$out" "$sha" || { record miss "$id"; return 1; }
    record ok "$id"; return 0
}

# A manual, login-gated resource: never downloaded. Validate if the user placed it, else instruct.
gated() {  # gated ID OUT "acquisition instructions"
    local id="$1" out="$2" how="$3"
    if [[ -f "$out" ]]; then
        [[ "$ACCEPT_LICENSE" == 1 ]] || warn "[$id] present but --accept-license not set; re-run with it to acknowledge terms"
        log "[$id] provided by user: $out"; record ok "$id"; return 0
    fi
    warn "[$id] LICENSE-GATED (login required) — not auto-downloadable. Provide it, then re-run:"
    printf '%s\n' "$how" | sed 's/^/    /' >&2
    warn "[$id] expected at: $out"
    record miss "$id"; return 1
}

# A public-but-non-commercial resource: direct URL exists, but require --accept-license first.
get_licensed() {  # get_licensed ID URL OUT "license note"
    local id="$1" url="$2" out="$3" note="$4"
    [[ -f "$out" ]] && { log "[$id] cached"; record skip "$id"; return 0; }
    if [[ "$ACCEPT_LICENSE" != 1 ]]; then
        warn "[$id] non-commercial license — re-run with --accept-license to download. Terms: $note"
        record miss "$id"; return 1
    fi
    get_free "$id" "$url" "$out" ""
}

# -------- per-resource preparation ------------------------------------------ #
prep_reference() {
    selected reference || return 0
    # DETERMINISTIC output path: fetch/verify/emit-env are separate invocations, so the prepared
    # file must land at exactly $REF_FASTA_OUT (line 51) — never a process-local reassignment.
    if [[ -f "$REF_FASTA_OUT" ]]; then log "[reference] cached"; record skip reference; return 0; fi
    get_free reference "$REF_FASTA_URL" "$REF_FASTA_OUT.dl" "$REF_FASTA_SHA256" || return 0
    # Ensembl serves gzip; decompress to a PLAIN faidx'd FASTA (universally accepted by
    # bcftools norm -f, VEP --fasta, samtools -T; no bgzf-random-access ambiguity).
    if zcat "$REF_FASTA_OUT.dl" > "$REF_FASTA_OUT" 2>/dev/null || gunzip -c "$REF_FASTA_OUT.dl" > "$REF_FASTA_OUT"; then
        rm -f "$REF_FASTA_OUT.dl"
    else
        warn "[reference] decompress failed; raw download left at $REF_FASTA_OUT.dl"; record miss reference; return 0
    fi
    need samtools && samtools faidx "$REF_FASTA_OUT" 2>/dev/null || warn "[reference] samtools faidx failed (need the .fai next to $REF_FASTA_OUT)"
}

prep_vep_cache() {
    selected vep_cache || return 0
    local tar="$DIR/vep_cache/cache.tar.gz"
    if [[ -d "$VEP_CACHE_OUT/homo_sapiens" ]]; then log "[vep_cache] cached"; record skip vep_cache; return 0; fi
    get_free vep_cache "$VEP_CACHE_URL" "$tar" "$VEP_CACHE_SHA256" || return 0
    log "[vep_cache] extracting (indexed cache, r${VEP_CACHE_VERSION})"
    mkdir -p "$VEP_CACHE_OUT"; tar -xzf "$tar" -C "$VEP_CACHE_OUT" && rm -f "$tar"
}

prep_gnomad() {
    selected gnomad_sites || return 0
    if [[ -f "$GNOMAD_OUT" ]]; then log "[gnomad_sites] cached"; record skip gnomad_sites; return 0; fi
    case "${GNOMAD_FILE_TMPL:-}" in ""|REPLACE_*) warn "[gnomad_sites] not pinned in manifest.env — skipping"; record miss gnomad_sites; return 0;; esac
    local base="" proto=""
    if   [[ "${GNOMAD_BASE_GS:-}"    != REPLACE_* && -n "${GNOMAD_BASE_GS:-}"    ]] && need gsutil; then base="$GNOMAD_BASE_GS"; proto=gs
    elif [[ "${GNOMAD_BASE_S3:-}"    != REPLACE_* && -n "${GNOMAD_BASE_S3:-}"    ]] && need aws;    then base="$GNOMAD_BASE_S3"; proto=s3
    elif [[ "${GNOMAD_BASE_HTTPS:-}" != REPLACE_* && -n "${GNOMAD_BASE_HTTPS:-}" ]];                then base="$GNOMAD_BASE_HTTPS"; proto=https
    else warn "[gnomad_sites] no usable bucket base / downloader (need gsutil, aws, or an https mirror)"; record miss gnomad_sites; return 0; fi
    need bcftools || die "bcftools required to slim gnomAD"
    local parts=() list="$DIR/gnomad/parts.txt"; mkdir -p "$DIR/gnomad"
    for c in $GNOMAD_CHROMS; do
        # shellcheck disable=SC2059
        local fn; fn="$(printf "$GNOMAD_FILE_TMPL" "$c")"
        local raw="$DIR/gnomad/raw.${fn}" slim="$DIR/gnomad/slim.chr${c}.vcf.gz" tmp="$DIR/gnomad/.slim.chr${c}.tmp.vcf.gz"
        # reuse a cached slim ONLY if it is a complete, index-valid bgzip — a truncated one from an
        # interrupted run must be rebuilt, never concatenated into the faf95 oracle (silent corruption)
        if [[ -f "$slim" ]] && bcftools index -n "$slim" >/dev/null 2>&1; then parts+=("$slim"); continue; fi
        rm -f "$slim" "$slim".tbi "$slim".csi "$tmp" "$tmp".tbi 2>/dev/null || true
        log "[gnomad_sites] chr${c}: fetch + slim to {$GNOMAD_KEEP_INFO}"
        fetch_to "${base}/${fn}" "$raw" || { warn "[gnomad_sites] chr${c} download failed"; record miss gnomad_sites; return 0; }
        # keep ONLY the INFO fields the pipeline reads; write to a temp and ATOMICALLY move it into
        # place only after tabix indexes it cleanly (which validates the BGZF EOF block) — so a
        # partial/truncated slim can never be trusted or reused.
        if bcftools annotate -x "^INFO/${GNOMAD_KEEP_INFO//,/,INFO/}" -Oz -o "$tmp" "$raw" \
               && tabix -f -p vcf "$tmp" 2>/dev/null; then
            mv -f "$tmp" "$slim"; mv -f "$tmp.tbi" "$slim.tbi"; rm -f "$raw"; parts+=("$slim")
        else
            warn "[gnomad_sites] chr${c} slim failed/truncated — not trusting it"; rm -f "$tmp" "$tmp.tbi" "$raw"; record miss gnomad_sites; return 0
        fi
    done
    [[ ${#parts[@]} -gt 0 ]] || { record miss gnomad_sites; return 0; }
    printf '%s\n' "${parts[@]}" > "$list"
    log "[gnomad_sites] concatenating ${#parts[@]} slim per-chrom files"
    bcftools concat -n -f "$list" -Oz -o "$GNOMAD_OUT" || bcftools concat -f "$list" -Oz -o "$GNOMAD_OUT"
    index_vcf "$GNOMAD_OUT"
    bcftools index -n "$GNOMAD_OUT" >/dev/null 2>&1 || { warn "[gnomad_sites] concatenated oracle failed its integrity check"; record miss gnomad_sites; return 0; }
    record ok gnomad_sites
}

prep_clinvar() {
    selected clinvar || return 0
    get_free clinvar "$CLINVAR_URL" "$CLINVAR_OUT" "" || return 0
    index_vcf "$CLINVAR_OUT"
}

prep_loftee() {
    selected loftee || return 0
    mkdir -p "$LOFTEE_OUT"
    # human_ancestor bgzip-FASTA ships WITH its .fai/.gzi — download all three, do NOT re-bgzip
    get_free loftee_human_ancestor "$LOFTEE_HUMAN_ANCESTOR_URL" "$LOFTEE_OUT/human_ancestor.fa.gz" "" || true
    get_free loftee_human_ancestor_fai "$LOFTEE_HUMAN_ANCESTOR_FAI_URL" "$LOFTEE_OUT/human_ancestor.fa.gz.fai" "" || true
    get_free loftee_human_ancestor_gzi "$LOFTEE_HUMAN_ANCESTOR_GZI_URL" "$LOFTEE_OUT/human_ancestor.fa.gz.gzi" "" || true
    # conservation_file must be the UNCOMPRESSED sqlite db
    if [[ ! -f "$LOFTEE_OUT/loftee.sql" ]]; then
        get_free loftee_conservation "$LOFTEE_CONSERVATION_URL" "$LOFTEE_OUT/loftee.sql.gz" "" \
            && gunzip -f "$LOFTEE_OUT/loftee.sql.gz" || true
    fi
    get_free loftee_gerp "$LOFTEE_GERP_URL" "$LOFTEE_OUT/gerp_conservation_scores.homo_sapiens.GRCh38.bw" "" || true  # 12 GB
    log "[loftee] use the LoF plugin CODE from branch '${LOFTEE_PLUGIN_BRANCH}' (master is GRCh37-only)"
}

prep_constraint() {
    selected constraint || return 0
    mkdir -p "$DIR/constraint"
    local gn="$DIR/constraint/gnomad_lof_metrics.txt.bgz" sh="$DIR/constraint/shet_zeng2024.tsv" ph="$DIR/constraint/phaplo_collins2022.tsv.gz"
    get_free constraint_gnomad "$CONSTRAINT_GNOMAD_URL" "$gn" "" || true
    get_free constraint_shet "$CONSTRAINT_SHET_URL" "$sh" "" || true
    get_free constraint_phaplo "$CONSTRAINT_PHAPLO_URL" "$ph" "" || true
    # left-join into ONE per-gene TSV (oe_lof_upper, pli, s_het, phaplo) — priors/tiers only
    if [[ -f "$gn" ]] && need python3; then
        log "[constraint] joining gnomAD LOEUF/pLI + s_het + pHaplo -> $CONSTRAINT_OUT"
        python3 "$HERE/join_constraint.py" --gnomad "$gn" --shet "$sh" --phaplo "$ph" --out "$CONSTRAINT_OUT" \
            && record ok constraint || warn "[constraint] join failed; join manually (see docs/gene_constraint.md)"
    fi
    get_free mutation_rate "$MUTRATE_URL" "$MUTRATE_OUT" "" || warn "[mutation_rate] optional (Step-6 de novo, secondary) — supply a Samocha per-gene rate TSV if wanted"
}

prep_dbnsfp() {
    selected dbnsfp || return 0
    [[ -f "$DBNSFP_OUT" ]] && { log "[dbnsfp] cached"; record skip dbnsfp; return 0; }
    get_licensed dbnsfp "$DBNSFP_URL" "$DIR/dbnsfp/dbNSFP${DBNSFP_VERSION}.zip" \
        "dbNSFP ${DBNSFP_VERSION} academic build — non-commercial + cite Liu et al. Use the 'a' (academic) build, NOT 'c' (commercial, omits REVEL/CADD/AlphaMissense)." || return 0
    need bcftools || die "need htslib tools to prepare dbNSFP"
    log "[dbnsfp] building GRCh38 table (VEP-plugin recipe; heavy sort under $DIR/dbnsfp/_sort)"
    ( cd "$DIR/dbnsfp" && unzip -o "dbNSFP${DBNSFP_VERSION}.zip" >/dev/null \
        && zcat "dbNSFP${DBNSFP_VERSION}_variant.chr1.gz" | head -n1 > h \
        && mkdir -p _sort \
        && zgrep -h -v '^#chr' dbNSFP${DBNSFP_VERSION}_variant.chr* | LC_ALL=C sort -T _sort -k1,1 -k2,2n \
             | cat h - | bgzip -c > "$DBNSFP_OUT" \
        && tabix -s 1 -b 2 -e 2 "$DBNSFP_OUT" \
        && rm -rf _sort h dbNSFP${DBNSFP_VERSION}_variant.chr* "dbNSFP${DBNSFP_VERSION}.zip" ) \
        && record ok dbnsfp || warn "[dbnsfp] prep failed — run the recipe in docs/resources.md manually"
}

prep_spliceai() {
    selected spliceai || return 0
    # SNV: free no-login Ensembl MANE mirror (MANE-select transcripts only). Full genome-wide set
    # is Illumina BaseSpace (login-gated) — provide it manually to cover non-MANE transcripts.
    if [[ ! -f "$SPLICEAI_SNV_OUT" ]]; then
        get_licensed spliceai_snv "$SPLICEAI_MANE_SNV_URL" "$SPLICEAI_SNV_OUT" \
            "SpliceAI non-commercial (Illumina). Free source is Ensembl MANE-only; full set needs BaseSpace login." \
            && { [[ -f "$SPLICEAI_SNV_OUT.tbi" ]] || index_vcf "$SPLICEAI_SNV_OUT"; }
    fi
    gated spliceai_indel "$SPLICEAI_INDEL_OUT" \
        "SpliceAI indel scores: download spliceai_scores.raw.indel.hg38.vcf.gz from Illumina BaseSpace
(${SPLICEAI_ILLUMINA_URL}; 'genome_scores_v1.3'), then: tabix -p vcf <file>. (No no-login mirror for indels.)" || true
}

prep_cadd() {
    selected cadd || return 0
    # CADD via the dedicated plugin is the COMPLETE CADD source (whole-genome SNVs, coding AND
    # non-coding, + the precomputed indel set) — more complete than dbNSFP's coding-only CADD_phred.
    # ~81 GB SNV file, non-commercial -> requires --accept-license. If your existing VEP setup
    # already holds these files, point CADD_SNV/CADD_INDEL at them and skip (--only without cadd).
    log "[cadd] fetching CADD v${CADD_VERSION} GRCh38 — whole-genome SNVs (~81 GB) + indels (the complete CADD source)"
    get_licensed cadd_snv "$CADD_SNV_URL" "$CADD_SNV_OUT" "CADD ${CADD_VERSION} non-commercial (UW)." \
        && get_free cadd_snv_tbi "$CADD_SNV_URL.tbi" "$CADD_SNV_OUT.tbi" "" || true
    get_licensed cadd_indel "$CADD_INDEL_URL" "$CADD_INDEL_OUT" "CADD ${CADD_VERSION} non-commercial (UW)." \
        && get_free cadd_indel_tbi "$CADD_INDEL_URL.tbi" "$CADD_INDEL_OUT.tbi" "" || true
}

# -------- verify ------------------------------------------------------------ #
verify_one() { [[ -e "$1" ]] && log "  ok   $2" || { warn "  MISS $2 ($1)"; MISSING+=("$2"); }; }
# Existence is NOT enough for the big UNCHECKSUMMED downloads (CADD ~81 GB): the manifest carries
# no usable sha for them, so get_free accepts any existing file as "cached" — an interrupted
# download is silently reused. For CADD that is the quiet science-loss path: VEP's tabix returns
# no score past the truncation, and CADD is the ONLY keep-path below MODERATE impact, so those
# variants vanish from the screen with no error. So integrity-check the compressed targets — a
# complete bgzip/gzip stream AND a tabix index — degrading to a warn if no (b)gzip is on PATH.
_gz_intact() {  # 0 if $1 is a complete gz/bgzip stream (or no tool to check with)
    if need bgzip; then bgzip -t "$1" 2>/dev/null && return 0 || return 1; fi
    if need gzip;  then gzip  -t "$1" 2>/dev/null && return 0 || return 1; fi
    warn "no bgzip/gzip to integrity-check $1 — skipping (install htslib or gzip to enable)"; return 0
}
verify_one_intact() {  # PATH ID — existence + (for .gz) complete stream + tabix index
    local f="$1" id="$2"
    if [[ ! -e "$f" ]]; then warn "  MISS $id ($f)"; MISSING+=("$id"); return; fi
    if ! _gz_intact "$f"; then
        warn "  BAD  $id ($f) — truncated/corrupt compressed stream; delete it and re-fetch"; MISSING+=("$id"); return; fi
    if [[ ! -f "$f.tbi" && ! -f "$f.csi" ]]; then
        warn "  BAD  $id ($f) — no .tbi/.csi index; VEP's per-variant tabix reads would fail. Re-fetch the index."; MISSING+=("$id"); return; fi
    log "  ok   $id (stream + index intact)"
}
# Present-but-not-required: report, never fail. These back the roadmap restorations
# (docs/ROADMAP.md "Restoring the VEP-only contract's gaps"); no step reads them today.
verify_extra() { [[ -e "$1" ]] && log "  ok   $2 (extra; not used by the pipeline)" || log "  --   $2 not prepared (not required)"; }
do_verify() {
    log "verifying prepared resources under $DIR"
    # REQUIRED = exactly what the VEP-only contract consumes. Verifying gnomAD/ClinVar/LOFTEE
    # here would fail every correctly-configured user, since nothing fetches or reads them.
    verify_one "$REF_FASTA_OUT" reference
    # Step 1's `norm -f` and VEP's --fasta both need the .fai; a run must not pass verify without it.
    if [[ -e "$REF_FASTA_OUT" && ! -f "$REF_FASTA_OUT.fai" ]]; then
        warn "  BAD  reference (no .fai next to $REF_FASTA_OUT — run: samtools faidx $REF_FASTA_OUT)"; MISSING+=("reference_fai"); fi
    verify_one "$VEP_CACHE_OUT/homo_sapiens" vep_cache
    verify_one_intact "$CADD_SNV_OUT" cadd_snv
    verify_one_intact "$CADD_INDEL_OUT" cadd_indel
    verify_one "$CONSTRAINT_OUT" constraint
    verify_extra "$GNOMAD_OUT" gnomad_sites
    verify_extra "$CLINVAR_OUT" clinvar
    verify_extra "$LOFTEE_OUT/human_ancestor.fa.gz" loftee
    verify_extra "$DBNSFP_OUT" dbnsfp
    verify_extra "$SPLICEAI_SNV_OUT" spliceai_snv
    [[ ${#MISSING[@]} -eq 0 ]] && { log "verify: all required resources present"; return 0; } || die "verify: ${#MISSING[@]} required resource(s) missing: ${MISSING[*]}"
}

# -------- emit-env ---------------------------------------------------------- #
do_emit() {
    local w=/dev/stdout; [[ -n "$EMIT_OUT" ]] && w="$EMIT_OUT"
    {
        echo "# hprv resource env — source before run_pipeline.sh (generated by prepare_resources.sh)"
        echo "# VEP-only contract: these are every \${ENV} placeholder config/config.example.yaml"
        echo "# still has a key for. gnomAD/ClinVar/dbNSFP/SpliceAI/LOFTEE exports were REMOVED —"
        echo "# nothing reads them, and exporting dead vars invites 'why is this ignored?'."
        echo "export REF_FASTA=$REF_FASTA_OUT"
        echo "export VEP_CACHE=$VEP_CACHE_OUT"
        echo "# VEP_PLUGINS is provided by the image (/plugins) — override only if you relocate plugin code"
        echo "export CADD_SNV=$CADD_SNV_OUT"
        echo "export CADD_INDEL=$CADD_INDEL_OUT"
        echo "export GNOMAD_V2_CONSTRAINT=$CONSTRAINT_OUT"
        echo "export MUTRATE_TABLE=$MUTRATE_OUT"
        echo "# Already have a VEP 115 GRCh38 VCF? Set this and Step 2 skips the VEP call entirely:"
        echo "# export VEP_ANNOTATED_VCF=/path/to/your.vep.vcf.gz"
    } > "$w"
    # `if`, not `[[ ]] && log`: as the LAST statement of the function a false test would make
    # do_emit return 1, and `set -e` would kill the script AFTER a successful emit to stdout.
    if [[ -n "$EMIT_OUT" ]]; then log "wrote env exports -> $EMIT_OUT"; fi
}

# -------- dispatch ---------------------------------------------------------- #
case "$MODE" in
    emit-env) do_emit; exit 0;;
    verify)   do_verify; exit 0;;
    fetch)
        # DEFAULT = only what the VEP-only contract consumes. gnomAD/ClinVar/dbNSFP/SpliceAI/
        # LOFTEE are NOT fetched by default: no step reads them, and defaulting them on would
        # start an ~877 GB gnomAD download for data nothing consumes. They remain reachable by
        # explicit `--only gnomad_sites,clinvar,...` so the roadmap restorations
        # (docs/ROADMAP.md) stay one flag away. Note the pinned dbNSFP URL is DEAD upstream
        # (S3 NoSuchBucket -> registration-gated); prefer the dedicated REVEL/AlphaMissense files.
        prep_reference; prep_vep_cache; prep_cadd; prep_constraint
        if [[ -n "$ONLY" ]]; then
            prep_gnomad; prep_clinvar; prep_loftee; prep_dbnsfp; prep_spliceai
        fi
        log "---------------------------------------------------------------"
        log "prepared: ${#PREPARED[@]}  cached: ${#SKIPPED[@]}  missing/gated: ${#MISSING[@]}"
        [[ ${#MISSING[@]} -gt 0 ]] && log "missing/gated: ${MISSING[*]}"
        log "next: prepare_resources.sh --dir $DIR emit-env --out $DIR/resources.env ; source it"
        ;;
    *) die "unknown mode: $MODE";;
esac
