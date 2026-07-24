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
# Step 8b (optional): if --kraken2-db is given, classify each candidate's ALT-supporting reads
# in the mini-CRAMs with kraken2 (via nonhuman-screen) and fold the non-human fraction into
# variants.tsv (child_/mother_/father_nhf + nhf_flag). A contamination / mis-mapping down-rank
# signal; rides on the mini-CRAMs already sliced above.
#
# Usage:
#   08_igv_export.sh --work WORKDIR --ref GRCh38.fa [--cram-map map.tsv] [--padding 1000]
#       [--kraken2-db DB [--nhf-members carriers|child_only|all] [--nhf-confidence 0.05]
#        [--nhf-min-reads 5] [--nhf-memory-mapping]]
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$HERE/lib/common.sh"
export PYTHONPATH="${PYTHONPATH:-}:${HPRV_HOME:-$(cd "$HERE/.." && pwd)}/src"

# --ref is the CRAM ENCODING reference (needed to decode reference-compressed CRAMs),
# which may differ from the variant-calling reference used elsewhere; run_pipeline.sh
# passes resources.cram_ref (falling back to reference.fasta). --jobs bounds how many
# CRAM slices run concurrently: mini-slices are small index-based reads, but on a flaky
# FUSE/SBFS mount unbounded concurrency destabilizes it, so this is capped deliberately.
WORK="" REF="${HPRV_CRAM_REF:-${HPRV_REF_FASTA:-}}" CRAM_MAP="${HPRV_CRAM_MAP:-}" PAD=1000 GENOME=hg38 JOBS=1
# Step 8b (non-human-fraction) options; empty KRAKEN2_DB => 8b disabled. Supplied ONLY via
# --kraken2-db (run_pipeline.sh passes it only when outputs.igv.nonhuman_screen.enabled AND a DB is
# set) — deliberately NOT defaulted from $HPRV_KRAKEN2_DB, so the enable gate lives in one place.
KRAKEN2_DB="" NHF_MEMBERS=carriers NHF_CONF=0.05 NHF_MIN_READS=5 NHF_MMAP=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --work) WORK="$2"; shift 2;;
        --ref) REF="$2"; shift 2;;
        --cram-map) CRAM_MAP="$2"; shift 2;;
        --padding) PAD="$2"; shift 2;;
        --genome) GENOME="$2"; shift 2;;
        --jobs) JOBS="$2"; shift 2;;
        --kraken2-db) KRAKEN2_DB="$2"; shift 2;;
        --nhf-members) NHF_MEMBERS="$2"; shift 2;;
        --nhf-confidence) NHF_CONF="$2"; shift 2;;
        --nhf-min-reads) NHF_MIN_READS="$2"; shift 2;;
        --nhf-memory-mapping) NHF_MMAP="--memory-mapping"; shift;;
        -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
        *) die "unknown arg: $1";;
    esac
done
[[ "$JOBS" =~ ^[0-9]+$ && "$JOBS" -ge 1 ]] || die "--jobs must be a positive integer"
[[ "$NHF_MIN_READS" =~ ^[0-9]+$ ]] || die "--nhf-min-reads must be a non-negative integer"
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

# --- pre-group per-trio inputs in ONE pass each (was O(trios x calls) re-scanning) ---
mkdir -p "$HPRV_TMPDIR"
# 1) candidate loci -> one padded BED per trio. Tag each call with its trio, sort so a
#    trio's rows are contiguous, then split into per-trio BEDs keeping a single file
#    handle open at a time (scales to many trios without hitting awk's open-file limit).
beddir="$HPRV_TMPDIR/callbeds"; mkdir -p "$beddir"
awk -F'\t' -v pad="$PAD" '
    NR==1{for(i=1;i<=NF;i++){if($i=="trio_id")ti=i; if($i=="chrom")ci=i; if($i=="pos")pi=i}; next}
    ti&&ci&&pi{s=$pi-1-pad; if(s<0)s=0; print $ti"\t"$ci"\t"s"\t"($pi+pad)}' "$calls" \
    | sort -k1,1 -k2,2 -k3,3n \
    | awk -v dir="$beddir" '$1==""{next} {if($1!=prev){if(prev!="")close(pf); pf=dir"/"$1".bed"; prev=$1} print $2"\t"$3"\t"$4 > pf}'
# 2) trio -> candidate VCF path (trios.candidates.tsv has one row per trio). Write the
#    first path per trio to a tiny lookup file so the loop below is O(1) per trio.
vcfmapdir="$HPRV_TMPDIR/candvcf"; mkdir -p "$vcfmapdir"
[[ -f "$cand" ]] && awk -F'\t' -v dir="$vcfmapdir" 'NR>1 && $1!="" && !seen[$1]++{f=dir"/"$1; print $2 > f; close(f)}' "$cand"

# --- Pass 1 (serial, cheap): per trio, build the merged BED + copy the VCF track,
#     and emit one CRAM-slice TASK per member into a task list. No FUSE reads here. ---
tasks="$HPRV_TMPDIR/extract_tasks.tsv"; : > "$tasks"
while IFS=$'\t' read -r trio vcf ped samples; do
    [[ "$trio" == "trio_id" || -z "$trio" ]] && continue
    IFS=',' read -r kid dad mom <<< "$samples"

    bed="$beddir/${trio}.bed"
    if [[ -s "$bed" ]]; then
        merged="$HPRV_TMPDIR/${trio}.merged.bed"
        hprv_run -- bedtools merge -i "$bed" > "$merged" 2>/dev/null || cp "$bed" "$merged"
        if is_set "$REF" && [[ -f "$REF" && "$HAVE_MAP" -eq 1 ]]; then
            odir="$DATA/crams/$trio"; mkdir -p "$odir"
            for pair in "child:$kid" "mother:$mom" "father:$dad"; do
                role="${pair%%:*}"; sample="${pair#*:}"
                [[ -n "$sample" ]] || continue
                src="$(cram_for "$sample")"
                [[ -n "$src" && -f "$src" ]] || { warn "  [$trio] no CRAM for $role $sample"; continue; }
                # trio<TAB>role<TAB>sample<TAB>src<TAB>merged<TAB>ocram
                printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$trio" "$role" "$sample" "$src" "$merged" "$odir/${sample}.cram" >> "$tasks"
            done
        fi
    fi

    # per-trio VCF track (copy + index the candidate VCF under the data-dir)
    cvcf=""
    if [[ -f "$vcfmapdir/$trio" ]]; then IFS= read -r cvcf < "$vcfmapdir/$trio" || true; fi
    if [[ -n "$cvcf" && -f "$cvcf" ]]; then
        cp -f "$cvcf" "$DATA/vcfs/${trio}.vcf.gz"
        hprv_run --bind "$DATA" -- bcftools index -t -f "$DATA/vcfs/${trio}.vcf.gz" 2>/dev/null || true
    fi
done < <(tail -n +2 "$resolved")

# --- Pass 2 (serial): slice each candidate region set from its source CRAM, one CRAM
#     open at a time. Reads are index-based (small); serial access is the gentlest pattern
#     for a flaky FUSE/SBFS mount (the multi-core part is samtools' own -@ $JOBS per slice).
#     Each slice validates (quickcheck) and retries once on a transient read failure; a
#     persistently-bad source CRAM warns and is skipped (never aborts the whole export).
#
#     IDEMPOTENT / RESUMABLE: Step 8 has no step-level .done (it must re-run to pick up a changed
#     candidate set), so each slice self-guards. A re-run — after a walltime kill, or to recover a
#     few slices that failed on a flaky mount — REUSES a valid prior mini-CRAM for the SAME region
#     set WITHOUT re-reading the multi-GB source CRAM off the mount. That reuse is the whole point
#     on a beta FUSE/SBFS mount; it also composes with Step 8b, which already .done-guards NHF. ---
# Content key of the candidate REGION SET, from the per-trio merged BED. Keyed to BED CONTENT, not
# mtime: Pass 1 rebuilds $merged every run, so its mtime always changes even when the candidate set
# is identical — only the content tells us whether the regions actually changed. cksum is POSIX
# (present on the host and in the image) and reads a small LOCAL file. A plain "file exists +
# quickcheck" guard would instead serve a STALE mini-CRAM after a re-run whose candidate set changed
# (curation re-ran): content-keying skips an unchanged set and re-slices a changed one.
_bed_key() { cksum < "$1" 2>/dev/null | awk '{print $1"-"$2}'; }

# Slice one candidate-region set from one source CRAM. Runs in a SUBSHELL (see the loop) so its
# `set +e` is local: we must handle samtools exit codes ourselves and retry/skip a transiently-bad
# read rather than let `set -e` abort the whole export on the first failure. Slices AND indexes in
# one pass (--write-index) because the igv.js server needs the .crai. Returns: 0 = freshly sliced,
# 3 = reused a valid cached slice (no source read), 1 = failed/skipped.
extract_one() {
    set +e
    local trio="$1" role="$2" sample="$3" src="$4" merged="$5" ocram="$6" attempt
    local donef="${ocram}.done" key; key="$(_bed_key "$merged")"
    # Idempotent skip — LOCAL reads only, NEVER $src (that is the point): a prior slice for the
    # CURRENT region key that is present, indexed, and passes quickcheck. Its .done records the key.
    if [[ -f "$donef" && -s "$ocram" && -f "$ocram.crai" && "$(cat "$donef" 2>/dev/null)" == "$key" ]] \
       && hprv_run -- samtools quickcheck "$ocram" 2>/dev/null; then
        return 3
    fi
    # Absent / stale (region set changed) / corrupt: clear any partial artifacts and (re)slice.
    rm -f "$ocram" "$ocram.crai" "$donef"
    for attempt in 1 2; do
        if hprv_run -- samtools view -C -@ "$JOBS" -T "$REF" --regions-file "$merged" \
                --write-index -o "$ocram" "$src" 2>/dev/null \
           && [[ -s "$ocram" && -f "$ocram.crai" ]] \
           && hprv_run -- samtools quickcheck "$ocram" 2>/dev/null; then
            printf '%s\n' "$key" > "$donef"   # stamp completion with the region key we sliced for
            return 0
        fi
        rm -f "$ocram" "$ocram.crai"
        [[ "$attempt" -eq 1 ]] && { warn "  [$trio] slice transient-failed for $role $sample; retrying"; sleep 2; }
    done
    warn "  [$trio] CRAM slice failed for $role $sample ($src); skipping this alignment track"
    return 1
}

# Extract SERIALLY — one source CRAM open at a time. On a flaky FUSE/SBFS mount this
# is the gentlest access pattern (the user's stated requirement: "keep the mount well
# behaved"); the multi-core part is samtools' own --threads (-@ $JOBS) per slice.
# Run in a single interactive job where the CRAM mount is available.
# NB: `grep -c` prints 0 AND exits 1 on an all-blank/empty file, so `|| echo 0` would
# append a second "0" -> a "0\n0" that breaks the `-gt` test. Keep grep's own count.
n_tasks=$(grep -cve '^[[:space:]]*$' "$tasks" 2>/dev/null || true); n_tasks=${n_tasks:-0}
n_extracted=0 n_sliced=0 n_cached=0
if [[ "$n_tasks" -gt 0 ]]; then
    log "Step 8: up to $n_tasks mini-CRAM(s), serial (FUSE-safe; $JOBS thread(s)/slice); valid cached slices are reused"
    while IFS=$'\t' read -r trio role sample src merged ocram; do
        [[ -z "$trio" ]] && continue
        # `( … ) || rc=$?` keeps extract_one's `set +e` in its subshell AND captures its distinct
        # return code (0 sliced / 3 cached / 1 failed) without tripping the loop's `set -e`.
        rc=0; ( extract_one "$trio" "$role" "$sample" "$src" "$merged" "$ocram" ) || rc=$?
        case "$rc" in
            0) n_extracted=$((n_extracted + 1)); n_sliced=$((n_sliced + 1)) ;;
            3) n_extracted=$((n_extracted + 1)); n_cached=$((n_cached + 1)) ;;
            *) : ;;   # failed/skipped-bad: already warned inside extract_one
        esac
    done < "$tasks"
    [[ "$n_cached" -gt 0 ]] && log "Step 8: reused $n_cached valid cached mini-CRAM(s) (idempotent resume — no source-CRAM reads); sliced $n_sliced"
fi

# --- Step 8b: non-human-fraction (NHF) annotation of ALT-supporting reads --------------------
# For each screened (trio, member): classify the mini-CRAM's ALT reads with kraken2 (via
# nonhuman-screen) against the per-trio candidate VCF -> nhf/<trio>/<sample>.variant_nhf.tsv,
# which build_variants_tsv folds into variants.tsv on the 0-based key (pos-1). Rides entirely on
# the mini-CRAMs + per-trio VCFs produced above — no new source-CRAM I/O.
#
# Gated OFF unless a kraken2 DB is given AND mini-CRAMs exist AND nonhuman-screen resolves in the
# runtime (it lives only in the image, so a bare host / CI warns and skips). All failures degrade
# to "leave the NHF columns blank", never abort the export.
#
# Perf: invocations run SERIALLY so that with --memory-mapping the OS page cache stays warm across
# trios and the (dominant) kraken2 DB-load is paid ~once per node, not once per invocation. Put the
# DB on local disk/NVMe. See docs/resources.md.
nhf_dir="$DATA/nhf"
if is_set "$KRAKEN2_DB"; then
    if [[ ! -d "$KRAKEN2_DB" ]]; then
        warn "Step 8b: kraken2 DB '$KRAKEN2_DB' is not a directory — skipping NHF screening"
    elif [[ "$HAVE_MAP" -ne 1 ]]; then
        warn "Step 8b: no mini-CRAMs (no --cram-map) — skipping NHF screening"
    elif ! ( is_set "$REF" && [[ -f "$REF" ]] ); then
        warn "Step 8b: CRAM reference '$REF' missing — skipping NHF screening"
    elif ! hprv_run --bind "$KRAKEN2_DB" -- nonhuman-screen --version >/dev/null 2>&1; then
        warn "Step 8b: 'nonhuman-screen' not available in the runtime (it ships only in the image) — skipping NHF screening"
    else
        db_ok=1
        for f in hash.k2d opts.k2d taxo.k2d; do
            [[ -f "$KRAKEN2_DB/$f" ]] || { warn "Step 8b: kraken2 DB missing $f — skipping NHF screening"; db_ok=0; break; }
        done
        # Taxonomy dumps are not strictly required by kraken2, but WITHOUT them classification
        # degrades to exact-taxid matching and the NHF signal is unreliable in both directions.
        # Warn loudly (do not silently trust); there is no per-variant flag for this in the output.
        if [[ "$db_ok" -eq 1 ]] \
           && ! { [[ -f "$KRAKEN2_DB/taxonomy/nodes.dmp" || -f "$KRAKEN2_DB/nodes.dmp" ]] \
                  && [[ -f "$KRAKEN2_DB/taxonomy/names.dmp" || -f "$KRAKEN2_DB/names.dmp" ]]; }; then
            warn "Step 8b: kraken2 DB lacks taxonomy nodes.dmp/names.dmp — NHF degrades to exact-taxid matching and is UNRELIABLE. Use a DB that ships the dumps (e.g. PrackenDB)."
        fi
        if [[ "$db_ok" -eq 1 ]]; then
            mkdir -p "$nhf_dir"
            # members=carriers: a parent is screened only where it carries the ALT in >= 1 of the
            # trio's calls (a hom-ref member has ~no ALT reads, so screening it just wastes a DB
            # load). The child is ALWAYS screened — it defines the candidate. NOTE: Step 5 writes
            # child_gt/mother_gt/father_gt from cyvcf2 `gt_bases`, so these columns are BASE-form
            # ('A/G', 'AT/A', './.'), NOT allele-index ('0/1') — a parent "carries the ALT" when a
            # GT token equals the row's `alt` base (rows are biallelic post `norm -m-`, so `alt` is
            # a single allele string). Testing for a digit here would NEVER match and silently
            # degrade `carriers` to `child_only` (regression guarded by tests/test_nhf_carriers.sh).
            carriers_list="$HPRV_TMPDIR/nhf_carriers.tsv"
            awk -F'\t' '
                NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
                h["alt"]{
                  t=$(h["trio_id"]); a=$(h["alt"]);
                  if (h["mother_gt"]) { n=split($(h["mother_gt"]),g,/[\/|]/); for(i=1;i<=n;i++) if(g[i]==a){print t"\tmother"; break} }
                  if (h["father_gt"]) { n=split($(h["father_gt"]),g,/[\/|]/); for(i=1;i<=n;i++) if(g[i]==a){print t"\tfather"; break} } }' \
                "$calls" | sort -u > "$carriers_list"
            is_carrier() { awk -F'\t' -v t="$1" -v r="$2" '$1==t&&$2==r{f=1} END{exit f?0:1}' "$carriers_list"; }

            n_nhf=0
            while IFS=$'\t' read -r trio _ _ samples; do   # cols: trio_id vcf ped samples
                [[ "$trio" == "trio_id" || -z "$trio" ]] && continue
                IFS=',' read -r kid dad mom <<< "$samples"
                tvcf="$DATA/vcfs/${trio}.vcf.gz"
                [[ -f "$tvcf" ]] || continue
                for pair in "child:$kid" "mother:$mom" "father:$dad"; do
                    role="${pair%%:*}"; sample="${pair#*:}"
                    [[ -n "$sample" ]] || continue
                    cram="$DATA/crams/$trio/${sample}.cram"
                    [[ -f "$cram" ]] || continue
                    case "$NHF_MEMBERS" in
                        child_only) [[ "$role" == "child" ]] || continue ;;
                        all)        : ;;
                        *)          [[ "$role" == "child" ]] || is_carrier "$trio" "$role" || continue ;;
                    esac
                    outp="$nhf_dir/$trio/${sample}"; out_tsv="${outp}.variant_nhf.tsv"
                    if [[ -f "$out_tsv" && -f "${out_tsv}.done" ]]; then n_nhf=$((n_nhf+1)); continue; fi
                    mkdir -p "$nhf_dir/$trio"
                    errf="$HPRV_TMPDIR/nhf.${trio}.${sample}.err"
                    # shellcheck disable=SC2086  # $NHF_MMAP is an intentional word (empty or --memory-mapping)
                    if hprv_run --bind "$KRAKEN2_DB" -- nonhuman-screen classify \
                            --bam "$cram" --variants "$tvcf" --ref-fasta "$REF" \
                            --kraken2-db "$KRAKEN2_DB" --confidence "$NHF_CONF" \
                            --threads "$JOBS" $NHF_MMAP --out-prefix "$outp" 2>"$errf" \
                       && [[ -f "$out_tsv" ]]; then
                        touch "${out_tsv}.done"; n_nhf=$((n_nhf+1))
                    else
                        warn "  [$trio] NHF screen failed for $role $sample (kraken2 missing, or a read/DB error): $(tail -n1 "$errf" 2>/dev/null); leaving its NHF blank"
                        rm -f "$out_tsv" "${outp}.summary.json"
                    fi
                    rm -f "$errf"
                done
            done < <(tail -n +2 "$resolved")
            log "Step 8b: NHF-screened $n_nhf (trio,member) mini-CRAM(s) (members=$NHF_MEMBERS, confidence=$NHF_CONF, memmap=${NHF_MMAP:-off})"
            audit 08_igv nhf_screened "$n_nhf"
        fi
    fi
fi

# --- assemble variants.tsv + sample_qc.tsv + trios.tsv + curation.json ---
# NHF columns fold in from nhf/<trio>/<sample>.variant_nhf.tsv when Step 8b produced them (the
# join is on the 0-based key, pos-1); absent files => blank NHF columns (legacy behavior).
python3 - "$calls" "$resolved" "$DATA" "$WORK/qc_report.tsv" "$NHF_MIN_READS" <<'PY'
import sys
from hprv import igv
calls, manifest, data, qc, nhf_min = sys.argv[1:6]
n = igv.build_variants_tsv(calls, manifest, data, f"{data}/variants.tsv",
                           nhf_dir=f"{data}/nhf", nhf_min_reads=int(nhf_min))
m = igv.write_sample_qc(qc, manifest, f"{data}/sample_qc.tsv")
sys.stderr.write(f"  variants.tsv rows: {n}; sample_qc rows: {m}\n")
PY

# trios.tsv (#kid mom dad) + empty curation
{ printf '#kid\tmom\tdad\n'
  awk -F'\t' 'NR>1{split($4,s,","); if(s[1]!="") print s[1]"\t"s[3]"\t"s[2]}' "$resolved"; } > "$DATA/trios.tsv"
[[ -f "$DATA/curation.json" ]] || echo '{}' > "$DATA/curation.json"
printf '{"genome": "%s"}\n' "$GENOME" > "$DATA/config.json"   # server --genome hint

audit 08_igv variants "$(($(grep -cve '^[[:space:]]*$' "$DATA/variants.tsv") - 1))"
audit 08_igv minicrams "$n_extracted"
log "Step 8 complete: igv review export -> $DATA (variants.tsv, crams/, vcfs/, trios.tsv, curation.json)"
