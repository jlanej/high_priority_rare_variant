#!/usr/bin/env bash
# =============================================================================
# pipeline/lib/common.sh — shared helpers for the high_priority_rare_variant
# pipeline. Source this from every step script:
#
#     source "$(dirname "$0")/lib/common.sh"
#
# Design goals (mirroring the project's engineering ethos):
#   * ONE container holds all analysis software. Every tool (bcftools, bedtools,
#     vep, slivar, python + cyvcf2 ...) is invoked through hprv_run so the same
#     pinned image is used everywhere and nothing depends on host-installed tools.
#   * Runs unchanged on a laptop (docker) and on HPC (apptainer/singularity).
#   * No hardcoded paths. Everything comes from the environment / config, with
#     safe, overridable defaults. This is a PUBLIC repo — keep it that way.
#   * Fail loudly, verify before claiming success, and be idempotent (.done files).
# =============================================================================

# Callers should already have `set -euo pipefail`; enforce it defensively.
set -euo pipefail

# --- Container image -----------------------------------------------------------
# The unified analysis image. Override with HPRV_IMAGE to pin a digest or point
# at a locally-pulled .sif. Default is the GHCR image built by CI. We deliberately
# do NOT bake an owner path into scripts beyond this single overridable default.
: "${HPRV_IMAGE:=ghcr.io/jlanej/high_priority_rare_variant:latest}"

# Container runtime: auto-detect unless set (apptainer|singularity|docker|native).
# Honors HPRV_ENGINE (the config key runtime.engine) as an alias. "native" runs tools
# directly on the host / inside the container (tools already on PATH).
: "${HPRV_RUNTIME:=${HPRV_ENGINE:-auto}}"

# Extra bind mounts (space-separated absolute dirs) always applied, e.g. a shared
# scratch or resources root. Per-call binds are added with `hprv_run --bind DIR`.
: "${HPRV_BIND:=}"

# TMPDIR for apptainer image extraction / container /tmp. Pointing this at real
# disk (not the default tmpfs) avoids the SLURM-cgroup tmpfs-OOM failure mode
# documented in the reference annotate script (apptainer --containall parks /tmp
# in RAM, which is billed to --mem and SIGKILLs heavy VEP jobs). Keep it on disk.
: "${HPRV_TMPDIR:=${TMPDIR:-/tmp}/hprv.$$}"

# --------------------------------------------------------------------------- #
# Logging
# --------------------------------------------------------------------------- #
_hprv_ts() { date +'%Y-%m-%dT%H:%M:%S%z'; }
log()   { printf '[%s] %s\n'      "$(_hprv_ts)" "$*" >&2; }
warn()  { printf '[%s] WARN: %s\n' "$(_hprv_ts)" "$*" >&2; }
die()   { printf '[%s] ERROR: %s\n' "$(_hprv_ts)" "$*" >&2; exit 1; }

# --------------------------------------------------------------------------- #
# Runtime detection
# --------------------------------------------------------------------------- #
hprv_detect_runtime() {
    if [[ "$HPRV_RUNTIME" != "auto" ]]; then
        printf '%s' "$HPRV_RUNTIME"; return 0
    fi
    if   command -v apptainer   >/dev/null 2>&1; then printf 'apptainer'
    elif command -v singularity >/dev/null 2>&1; then printf 'singularity'
    elif command -v docker      >/dev/null 2>&1; then printf 'docker'
    else printf 'native'
    fi
}

# --------------------------------------------------------------------------- #
# hprv_run [--bind DIR]... -- CMD [ARGS...]
#
# Execute a command inside the analysis container. Binds are host dirs that the
# command must read/write. For apptainer/singularity we bind + set a disk-backed
# workdir; for docker we -v mount and run as the calling uid to avoid root-owned
# outputs. `native` just execs the command on the host.
# --------------------------------------------------------------------------- #
hprv_run() {
    local -a binds=()
    while [[ $# -gt 0 && "$1" == "--bind" ]]; do
        binds+=("$2"); shift 2
    done
    [[ "${1:-}" == "--" ]] && shift
    [[ $# -gt 0 ]] || die "hprv_run: no command given"

    # Fold in always-on binds (HPRV_BIND) plus each requested dir.
    local -a allbinds=()
    local b
    for b in $HPRV_BIND "${binds[@]:-}"; do
        [[ -n "$b" ]] && allbinds+=("$b")
    done

    local runtime; runtime="$(hprv_detect_runtime)"
    mkdir -p "$HPRV_TMPDIR"

    case "$runtime" in
        apptainer|singularity)
            local -a bindargs=()
            for b in "${allbinds[@]:-}"; do [[ -n "$b" ]] && bindargs+=(--bind "$b"); done
            APPTAINER_TMPDIR="$HPRV_TMPDIR" SINGULARITY_TMPDIR="$HPRV_TMPDIR" \
            "$runtime" exec \
                --workdir "$HPRV_TMPDIR" \
                "${bindargs[@]:-}" \
                "$HPRV_IMAGE" "$@"
            ;;
        docker)
            local -a vargs=()
            for b in "${allbinds[@]:-}"; do [[ -n "$b" ]] && vargs+=(-v "$b:$b"); done
            docker run --rm \
                -u "$(id -u):$(id -g)" \
                -e HOME=/tmp \
                -v "$PWD:$PWD" -w "$PWD" \
                "${vargs[@]:-}" \
                "$HPRV_IMAGE" "$@"
            ;;
        native)
            # `command` bypasses the bcftools()/bedtools() convenience FUNCTIONS below,
            # which would otherwise re-enter hprv_run and recurse infinitely.
            command "$@"
            ;;
        *) die "hprv_run: unknown runtime '$runtime'";;
    esac
}

# --------------------------------------------------------------------------- #
# Convenience wrappers (all go through the one image)
# --------------------------------------------------------------------------- #
bcftools() { hprv_run -- bcftools "$@"; }
bedtools() { hprv_run -- bedtools "$@"; }
samtools() { hprv_run -- samtools "$@"; }
tabix()    { hprv_run -- tabix "$@"; }
bgzip()    { hprv_run -- bgzip "$@"; }

# --------------------------------------------------------------------------- #
# Integrity / idempotency helpers
# --------------------------------------------------------------------------- #

# Fail unless the bgzip file is complete (detects truncation).
require_intact_bgzip() {
    local f="$1"
    [[ -s "$f" ]] || die "missing or empty: $f"
    hprv_run -- bgzip -t "$f" 2>/dev/null \
        || gzip -t "$f" 2>/dev/null \
        || die "corrupt/truncated bgzip: $f"
}

# Index a VCF/BCF (tabix .tbi for .vcf.gz, .csi otherwise). Idempotent.
index_vcf() {
    local vcf="$1" d
    d="$(cd "$(dirname "$vcf")" && pwd)"
    if [[ "$vcf" == *.vcf.gz ]]; then
        [[ -f "${vcf}.tbi" || -f "${vcf}.csi" ]] && return 0
        hprv_run --bind "$d" -- bcftools index -t "$vcf"
    else
        [[ -f "${vcf}.csi" ]] && return 0
        hprv_run --bind "$d" -- bcftools index "$vcf"
    fi
}

# Count variant records (data lines) in a VCF/BCF.
count_variants() {
    local vcf="$1" d
    d="$(cd "$(dirname "$vcf")" && pwd)"
    hprv_run --bind "$d" -- bcftools view -H "$vcf" 2>/dev/null | wc -l | tr -d ' '
}

# .done idempotency: `is_done OUT` returns 0 if OUT and OUT.done both exist and
# OUT is a non-empty, intact bgzip. `mark_done OUT` stamps completion.
is_done() {
    local out="$1"
    [[ -f "${out}.done" && -s "$out" ]] || return 1
    if [[ "$out" == *.gz ]]; then require_intact_bgzip "$out" 2>/dev/null || return 1; fi
    return 0
}
mark_done() { touch "${1}.done"; }

# Resolve the absolute directory of a path (dir must exist).
abspath_dir() { cd "$(dirname "$1")" && pwd; }

# --------------------------------------------------------------------------- #
# Auditing: append (step, scope, metric, value) to $HPRV_AUDIT_DIR/counts.tsv so
# every count is recoverable. scope = "global" or a trio_id. No-op if unset.
#   audit STEP METRIC VALUE [SCOPE]
# --------------------------------------------------------------------------- #
: "${HPRV_AUDIT_DIR:=}"
audit() {
    [[ -n "$HPRV_AUDIT_DIR" ]] || return 0
    mkdir -p "$HPRV_AUDIT_DIR"
    local f="$HPRV_AUDIT_DIR/counts.tsv"
    [[ -f "$f" ]] || printf 'timestamp\tstep\tscope\tmetric\tvalue\n' > "$f"
    printf '%s\t%s\t%s\t%s\t%s\n' "$(_hprv_ts)" "$1" "${4:-global}" "$2" "$3" >> "$f"
}

# --------------------------------------------------------------------------- #
# Cleanup of the per-run tmpdir on exit (best-effort; only what we created).
# --------------------------------------------------------------------------- #
_hprv_cleanup() { [[ -n "${HPRV_TMPDIR:-}" && -d "$HPRV_TMPDIR" ]] && rm -rf "$HPRV_TMPDIR" || true; }
trap _hprv_cleanup EXIT
