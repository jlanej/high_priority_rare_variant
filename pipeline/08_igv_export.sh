#!/usr/bin/env bash
# =============================================================================
# 08_igv_export.sh  —  Pipeline Step 8: export for the igv.js trio review server
#
# Produces, under WORK/igv/ (the server's --data-dir):
#   variants.tsv     one row per candidate call (igv.js variant-review schema:
#                    chrom/pos/ref/alt + inheritance mode + genotypes + annotations
#                    + per-member track columns); any extra column is filterable.
#   crams/<trio>/<sample>.cram(+.crai)   mini-CRAMs sliced around candidate loci
#                    for child/mother/father (only if a sample->CRAM map is given).
#   vcfs/<trio>.vcf.gz(+.tbi)            per-trio candidate VCF track.
#   sample_qc.tsv    trio_id/role/sample_id + QC metrics (optional --sample-qc input).
#   trios.tsv        #kid mom dad.
#   curation.json    empty review state ({}).
#
# Serve with the jlanej/igv.js variant-review server:
#   node server.js --variants WORK/igv/variants.tsv --data-dir WORK/igv --genome hg38
#   (container: ghcr.io/jlanej/igv-variant-review)
#
# Usage:
#   08_igv_export.sh --work WORKDIR --ref GRCh38.fa [--cram-map map.tsv] [--padding 1000]
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"
export PYTHONPATH="${PYTHONPATH:-}:${HPRV_HOME:-$(cd "$HERE/.." && pwd)}/src"

WORK="" REF="${HPRV_REF_FASTA:-}" CRAM_MAP="${HPRV_CRAM_MAP:-}" PAD=1000
while [[ $# -gt 0 ]]; do
    case "$1" in
        --work) WORK="$2"; shift 2;;
        --ref) REF="$2"; shift 2;;
        --cram-map) CRAM_MAP="$2"; shift 2;;
        --padding) PAD="$2"; shift 2;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done
[[ -n "$WORK" ]] || die "need --work"
calls="$WORK/candidates.calls.tsv"; resolved="$WORK/trios.resolved.tsv"
cand="$WORK/trios.candidates.tsv"
[[ -f "$calls" && -f "$resolved" ]] || die "missing $calls or $resolved (run steps 4-6 first)"

is_set() { [[ -n "${1:-}" && "$1" != *'${'* ]]; }
DATA="$WORK/igv"; mkdir -p "$DATA/crams" "$DATA/vcfs"

# --- source-CRAM map (sample -> path); optional. Looked up on demand (no assoc
#     arrays, for bash 3.2 portability), matching the bespoke get_cram approach. ---
HAVE_MAP=0
binds="$DATA"; is_set "$REF" && [[ -f "$REF" ]] && binds+=" $(abspath_dir "$REF")"
cram_for() { [[ -f "${CRAM_MAP:-/nonexistent}" ]] || return 0; awk -F'\t' -v s="$1" '$1==s{print $2; exit}' "$CRAM_MAP"; }
if is_set "$CRAM_MAP" && [[ -f "$CRAM_MAP" ]]; then
    HAVE_MAP=1
    while IFS=$'\t' read -r s pth _; do
        [[ -z "$s" || "$s" == \#* ]] && continue
        [[ -f "$pth" ]] && binds+=" $(abspath_dir "$pth")"
    done < "$CRAM_MAP"
    log "Step 8: using CRAM map $CRAM_MAP"
else
    warn "no --cram-map; variants.tsv will be emitted without alignment tracks"
fi
HPRV_BIND="$(printf '%s\n' $binds | sort -u | tr '\n' ' ')"; export HPRV_BIND

# --- per-trio: BED of candidate loci -> mini-CRAM extraction + VCF track ---
n_extracted=0
while IFS=$'\t' read -r trio vcf ped samples; do
    [[ "$trio" == "trio_id" || -z "$trio" ]] && continue
    IFS=',' read -r kid dad mom <<< "$samples"

    bed="$HPRV_TMPDIR/${trio}.bed"; mkdir -p "$HPRV_TMPDIR"
    awk -F'\t' -v t="$trio" -v pad="$PAD" '
        NR==1{for(i=1;i<=NF;i++){if($i=="trio_id")ti=i; if($i=="chrom")ci=i; if($i=="pos")pi=i}; next}
        $ti==t{s=$pi-1-pad; if(s<0)s=0; print $ci"\t"s"\t"($pi+pad)}' "$calls" \
        | sort -k1,1 -k2,2n > "$bed"
    [[ -s "$bed" ]] || { continue; }
    merged="$HPRV_TMPDIR/${trio}.merged.bed"
    hprv_run -- bedtools merge -i "$bed" > "$merged" 2>/dev/null || cp "$bed" "$merged"

    if is_set "$REF" && [[ -f "$REF" && "$HAVE_MAP" -eq 1 ]]; then
        for pair in "child:$kid" "mother:$mom" "father:$dad"; do
            role="${pair%%:*}"; sample="${pair#*:}"
            src="$(cram_for "$sample")"
            [[ -n "$src" && -f "$src" ]] || { warn "  [$trio] no CRAM for $role $sample"; continue; }
            odir="$DATA/crams/$trio"; mkdir -p "$odir"
            ocram="$odir/${sample}.cram"
            hprv_run -- samtools view -C -T "$REF" --regions-file "$merged" -o "$ocram" "$src"
            hprv_run -- samtools index "$ocram"
            n_extracted=$((n_extracted + 1))
        done
    fi

    # per-trio VCF track (copy + index the candidate VCF under the data-dir)
    if [[ -f "$cand" ]]; then
        cvcf="$(awk -F'\t' -v t="$trio" '$1==t{print $2}' "$cand" | head -1)"
        if [[ -n "$cvcf" && -f "$cvcf" ]]; then
            cp -f "$cvcf" "$DATA/vcfs/${trio}.vcf.gz"
            hprv_run --bind "$DATA" -- bcftools index -t -f "$DATA/vcfs/${trio}.vcf.gz" 2>/dev/null || true
        fi
    fi
done < <(tail -n +2 "$resolved")

# --- assemble variants.tsv + sample_qc.tsv + trios.tsv + curation.json ---
python3 - "$calls" "$resolved" "$DATA" "$WORK/qc_report.tsv" <<'PY'
import sys
from hprv import igv
calls, manifest, data, qc = sys.argv[1:5]
n = igv.build_variants_tsv(calls, manifest, data, f"{data}/variants.tsv")
m = igv.write_sample_qc(qc, manifest, f"{data}/sample_qc.tsv")
sys.stderr.write(f"  variants.tsv rows: {n}; sample_qc rows: {m}\n")
PY

# trios.tsv (#kid mom dad) + empty curation
{ printf '#kid\tmom\tdad\n'
  awk -F'\t' 'NR>1{split($4,s,","); if(s[1]!="") print s[1]"\t"s[3]"\t"s[2]}' "$resolved"; } > "$DATA/trios.tsv"
[[ -f "$DATA/curation.json" ]] || echo '{}' > "$DATA/curation.json"

audit 08_igv variants "$(($(grep -cve '^[[:space:]]*$' "$DATA/variants.tsv") - 1))"
audit 08_igv minicrams "$n_extracted"
log "Step 8 complete: igv review export -> $DATA (variants.tsv, crams/, vcfs/, trios.tsv, curation.json)"
