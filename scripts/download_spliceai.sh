#!/usr/bin/env bash
# =============================================================================
# download_spliceai.sh — fetch ONLY the raw hg38 SpliceAI score files from the
# Illumina BaseSpace "genome_scores" dataset, using an authenticated `bs` CLI.
#
# The full raw genome-wide SpliceAI scores are login-gated (BaseSpace); there is no
# no-login mirror for the full set (the free Ensembl mirror is MANE-only SNV). This
# automates the fetch once you have (1) `bs` installed + `bs auth`-ed and (2) accepted
# the shared dataset into your account by opening the share link in a browser:
#     https://basespace.illumina.com/s/otSPW8hnhaZR
#
# It downloads ONLY the 4 files this pipeline uses (raw snv/indel hg38 + their .tbi) —
# NOT the whole dataset (which also carries masked + hg19 builds, ~100 GB). Idempotent:
# a file already present and intact is skipped. Verifies bgzip integrity + a tabix query,
# tabix-indexes anything missing an index, and (with --ref) checks contig naming.
#
# Usage:
#   scripts/download_spliceai.sh --dir /path/to/resources/spliceai [--ref GRCh38.fa]
#       [--bs $HOME/bin/bs] [--dataset-id ds.XXXX]
#
# Then:  source <(scripts/download_spliceai.sh --dir DIR --print-env)   # emits export lines
# =============================================================================
set -euo pipefail

BS="${BS:-$HOME/bin/bs}"
DIR="" REF="" DATASET_ID="" PRINT_ENV=0
NAME_RE="genome_scores"           # dataset/project name pattern to auto-resolve

die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log()  { printf '[download_spliceai] %s\n' "$*" >&2; }

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dir) DIR="$2"; shift 2;;
        --ref) REF="$2"; shift 2;;
        --bs) BS="$2"; shift 2;;
        --dataset-id) DATASET_ID="$2"; shift 2;;
        --print-env) PRINT_ENV=1; shift;;
        -h|--help) usage;;
        *) die "unknown arg: $1 (see --help)";;
    esac
done
[[ -n "$DIR" ]] || die "need --dir (where to put the score files)"

# The exact files this pipeline consumes (config resources.vep.spliceai_snv/indel).
SNV="spliceai_scores.raw.snv.hg38.vcf.gz"
INDEL="spliceai_scores.raw.indel.hg38.vcf.gz"
WANT=("$SNV" "$SNV.tbi" "$INDEL" "$INDEL.tbi")

# --print-env just emits the export lines (for `source <(... --print-env)`); no network.
if [[ "$PRINT_ENV" -eq 1 ]]; then
    printf 'export SPLICEAI_SNV=%s/%s\n'   "$DIR" "$SNV"
    printf 'export SPLICEAI_INDEL=%s/%s\n' "$DIR" "$INDEL"
    exit 0
fi

command -v tabix >/dev/null 2>&1 || die "tabix not on PATH (needed to verify/index; e.g. 'conda install -c bioconda htslib')"
[[ -x "$BS" ]] || command -v "$BS" >/dev/null 2>&1 || die "bs CLI not found at '$BS' — set --bs or install it (see docs/resources.md#spliceai)"
"$BS" whoami >/dev/null 2>&1 || die "bs is not authenticated — run: $BS auth   (then re-run this script)"
mkdir -p "$DIR"

# --- resolve the dataset id (auto by name, unless --dataset-id given) -----------------
if [[ -z "$DATASET_ID" ]]; then
    log "resolving the '$NAME_RE' dataset id via '$BS list datasets'..."
    # -f csv gives a header + rows; find the row whose name matches and pull the ds.* / numeric id.
    DATASET_ID="$("$BS" list datasets -f csv 2>/dev/null \
        | awk -F, -v re="$NAME_RE" 'tolower($0) ~ tolower(re){
              for(i=1;i<=NF;i++){gsub(/^ *"?|"? *$/,"",$i); if($i ~ /^(ds\.[0-9a-f]+|[0-9]{5,})$/){print $i; exit}}}')" || true
    if [[ -z "$DATASET_ID" ]]; then
        log "could not auto-find it under 'list datasets'; trying projects..."
        printf '%s\n' "-- your datasets (for reference) --" >&2
        "$BS" list datasets -f csv 2>/dev/null | sed 's/^/    /' >&2 || true
        die "could not auto-resolve the SpliceAI dataset id. Find it manually with:
    $BS list datasets            # look for genome_scores / SpliceAI
then re-run with:  --dataset-id <ds.XXXX or numeric id>"
    fi
fi
log "dataset id: $DATASET_ID"

# --- list its contents once: build a  name<TAB>fileid  map ----------------------------
CONTENTS="$("$BS" contents dataset -i "$DATASET_ID" -f csv 2>/dev/null)" \
    || die "failed to list dataset contents ('$BS contents dataset -i $DATASET_ID'). Is the id right / accepted into your account?"
file_id() {  # $1 = exact filename -> its numeric file id (or empty)
    # Match the Name as a whole CSV FIELD (not a substring): the ".vcf.gz" name is a prefix of the
    # ".vcf.gz.tbi" name, so a substring test could grab the .tbi row's id depending on row order.
    printf '%s\n' "$CONTENTS" | awk -F, -v n="$1" '
        { row=0
          for(i=1;i<=NF;i++){ f=$i; gsub(/^ *"?|"? *$/,"",f); if(f==n) row=1 }
          if(row){ for(i=1;i<=NF;i++){ f=$i; gsub(/^ *"?|"? *$/,"",f); if(f ~ /^[0-9]{4,}$/){print f; exit} } } }'
}

# --- verify helper: present, non-empty, intact bgzip (index files: just present+non-empty) ---
_valid() {  # $1 = path
    [[ -s "$1" ]] || return 1
    case "$1" in *.tbi) return 0;; esac
    bgzip -t "$1" 2>/dev/null || gzip -t "$1" 2>/dev/null
}

# --- download each wanted file individually (skip already-valid) ----------------------
n_new=0
for name in "${WANT[@]}"; do
    out="$DIR/$name"
    if _valid "$out"; then log "cached, skipping: $name"; continue; fi
    fid="$(file_id "$name")"
    if [[ -z "$fid" ]]; then
        # .tbi may be absent from the dataset — we regenerate it below; a missing DATA file is fatal.
        case "$name" in
            *.tbi) log "no $name in the dataset — will tabix-index locally"; continue;;
            *) die "file '$name' not found in dataset $DATASET_ID. Contents were:
$(printf '%s\n' "$CONTENTS" | sed 's/^/    /')";;
        esac
    fi
    log "downloading $name (file id $fid) — this is large (~27 GB for the SNV file)..."
    "$BS" download file -i "$fid" -o "$DIR" || die "download failed for $name (id $fid)"
    # bs may nest the file in a subdir; normalize to $DIR/$name.
    if [[ ! -f "$out" ]]; then
        found="$(find "$DIR" -maxdepth 3 -name "$name" -type f 2>/dev/null | head -1)"
        [[ -n "$found" ]] && mv -f "$found" "$out"
    fi
    _valid "$out" || die "downloaded $name is missing or a corrupt/truncated bgzip — re-run to resume"
    n_new=$((n_new + 1))
done

# --- ensure both DATA files exist, then index any missing .tbi ------------------------
for f in "$DIR/$SNV" "$DIR/$INDEL"; do
    [[ -s "$f" ]] || die "expected data file missing after download: $f"
    if [[ ! -f "$f.tbi" ]]; then log "tabix-indexing $(basename "$f")..."; tabix -p vcf "$f"; fi
done

# --- verify each is queryable (a tabix fetch returns >=1 record) ----------------------
for f in "$DIR/$SNV" "$DIR/$INDEL"; do
    c0="$(tabix -l "$f" 2>/dev/null | head -1)"
    [[ -n "$c0" ]] || die "no contigs in $(basename "$f") index — the file/index is broken"
    n="$(tabix "$f" "$c0" 2>/dev/null | head -1000 | wc -l | tr -d ' ')"
    [[ "$n" -gt 0 ]] || die "tabix query on $(basename "$f") ($c0) returned 0 records — file/index mismatch"
done
log "verified: both raw hg38 score files present, intact, and queryable."

# --- contig-naming trap: the score VCFs are Ensembl-style ('1','X'); if your reference is
#     'chr'-prefixed, the SpliceAI plugin's tabix queries silently return nothing (exit 0).
#     Warn loudly (never rewrite) so you catch it before a whole annotation run comes back empty. ---
if [[ -n "$REF" && -f "$REF.fai" ]]; then
    ref_chr="$(cut -f1 "$REF.fai" | grep -qx chr1 && echo chr || echo nochr)"
    sai_chr="$(tabix -l "$DIR/$SNV" 2>/dev/null | grep -qx chr1 && echo chr || echo nochr)"
    if [[ "$ref_chr" != "$sai_chr" ]]; then
        log "WARNING: CONTIG-NAMING MISMATCH — your reference is '${ref_chr}'-style but the SpliceAI VCF is '${sai_chr}'-style."
        log "         VEP maps to the cache's naming internally, so this is usually fine — BUT if Step 2 warns 'no vep_SpliceAI_pred_DS_* lifted',"
        log "         this is why. Re-key the score VCF to match (bcftools annotate --rename-chrs) if so."
    else
        log "contig naming matches the reference ('${ref_chr}'-style)."
    fi
fi

log "done. $n_new file(s) downloaded. Point the config at them:"
printf 'export SPLICEAI_SNV=%s/%s\n'   "$DIR" "$SNV"
printf 'export SPLICEAI_INDEL=%s/%s\n' "$DIR" "$INDEL"
