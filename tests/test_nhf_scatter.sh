#!/usr/bin/env bash
# =============================================================================
# tests/test_nhf_scatter.sh — Step-8b per-trio scatter/gather regression test.
#
# Runs the REAL pipeline/08_igv_export.sh (no copy of its logic, so no drift) with a
# stubbed `nonhuman-screen`, over a TWO-trio fixture. Asserts the properties that make
# a per-trio job array safe:
#
#   1. --nhf-emit-manifest lists only trios with outstanding work, and SHRINKS after
#      one is screened (so a resubmit schedules a smaller array).
#   2. --nhf-trio screens ONLY that trio (per-trio isolation).
#   3. Two scatter tasks running CONCURRENTLY do not corrupt shared state:
#        * each gets its own nhf_carriers.<trio>.tsv — the shared path is written with
#          `sort -u > file`, which TRUNCATES, so a reader could otherwise see a partial
#          list and a real carrier parent would read as non-carrier (the exact
#          degradation tests/test_nhf_carriers.sh guards);
#        * the per-trio VCF is NOT rewritten (8b's content key hashes it, so a torn
#          `cp -f` would give an unstable key / a .done against a partial write).
#   4. Scattered output is BYTE-IDENTICAL to the serial path.
#   5. --nhf-gather with an unscreened member WARNS and still writes variants.tsv —
#      never dies, never implies completeness (unlike Step 2's --annotate-gather, an
#      unscreened member here is legitimate: `members: carriers` skips hom-ref parents).
#   6. PCR/optical DUPLICATES are excluded from the mini-CRAMs (samtools -F 3844), and the
#      mask is part of the slice .done key so changing it re-slices. nonhuman-screen does
#      no flag filtering of its own and de-dups only by read NAME (which collapses the two
#      mates of one pair, NOT duplicates of one fragment), so without this N copies of a
#      contaminating fragment count as N independent ALT reads and inflate *_nhf_reads.
#
# Needs bcftools/samtools/bgzip/tabix + a python3 (stdlib only); self-skips if absent.
# Run: bash tests/test_nhf_scatter.sh
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

for t in bcftools samtools bgzip tabix python3; do
    command -v "$t" >/dev/null 2>&1 || { echo "SKIP test_nhf_scatter: $t not on PATH"; exit 0; }
done

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
W="$T/work"; mkdir -p "$W" "$T/bin"
cd "$T"

python3 -c 'open("ref.fa","w").write(">chr1\n"+("ACGTACGTAC"*60)+"\n")'
samtools faidx ref.fa

mk_cram() {  # $1 = sample id
    local s="$1" seq qual
    seq="$(python3 -c 'print("A"*40)')"; qual="$(python3 -c 'print(chr(73)*40)')"
    { printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:600\n@RG\tID:%s\tSM:%s\n' "$s" "$s"
      for p in 90 190; do
          printf '%s_%s\t0\tchr1\t%s\t60\t40M\t*\t0\t0\t%s\t%s\tRG:Z:%s\n' "$s" "$p" "$p" "$seq" "$qual" "$s"
      done
    } > "$s.sam"
    samtools view -C -T ref.fa -o "$s.cram" "$s.sam"; samtools index "$s.cram"
}
# two trios, each kid + carrier mother (so `carriers` screens 2 members per trio)
for s in KID1 MOM1 DAD1 KID2 MOM2 DAD2; do mk_cram "$s"; done

# KID1 additionally gets 3 PCR-DUPLICATE reads (flag 0x400) over chr1:100. nonhuman-screen does no
# flag filtering and de-dups only by read NAME, so without the slice-time -F mask these would count
# as 3 extra independent ALT reads and inflate the *_nhf_reads denominator.
mk_cram_with_dups() {
    local s="$1" seq qual
    seq="$(python3 -c 'print("A"*40)')"; qual="$(python3 -c 'print(chr(73)*40)')"
    { printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:600\n@RG\tID:%s\tSM:%s\n' "$s" "$s"
      for p in 90 190; do
          printf '%s_%s\t0\tchr1\t%s\t60\t40M\t*\t0\t0\t%s\t%s\tRG:Z:%s\n' "$s" "$p" "$p" "$seq" "$qual" "$s"
      done
      for d in 1 2 3; do
          printf '%s_dup%s\t1024\tchr1\t90\t60\t40M\t*\t0\t0\t%s\t%s\tRG:Z:%s\n' "$s" "$d" "$seq" "$qual" "$s"
      done
      # a SUPPLEMENTARY alignment (0x800): GATK HC uses these, so the slice must KEEP it
      printf '%s_supp\t2048\tchr1\t90\t60\t40M\t*\t0\t0\t%s\t%s\tRG:Z:%s\n' "$s" "$seq" "$qual" "$s"
    } > "$s.sam"
    # sort: the dup/supp records above are emitted at pos 90 AFTER the pos-190 read, and an
    # unsorted BAM/CRAM cannot be indexed ("Unsorted positions on sequence #1").
    samtools sort -O cram --reference ref.fa -o "$s.cram" "$s.sam"; samtools index "$s.cram"
}
mk_cram_with_dups KID1
dups_in_source=$(samtools view -c -f 1024 KID1.cram)
supp_in_source=$(samtools view -c -f 2048 KID1.cram)
: > map.tsv
for s in KID1 MOM1 DAD1 KID2 MOM2 DAD2; do printf '%s\t%s/%s.cram\n' "$s" "$T" "$s" >> map.tsv; done

mk_vcf() {  # $1 = trio suffix; a 2-variant biallelic candidate VCF
    { printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=600>\n'
      printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
      printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tKID%s\tDAD%s\tMOM%s\n' "$1" "$1" "$1"
      printf 'chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\t0/1\n'
      printf 'chr1\t200\t.\tC\tG\t.\t.\t.\tGT\t0/1\t0/0\t0/0\n'; } > "cand$1.vcf"
    bgzip -f "cand$1.vcf"; bcftools index -t "cand$1.vcf.gz"
}
mk_vcf 1; mk_vcf 2

{ printf 'trio_id\tvcf\tped\tsamples\n'
  printf 'T1\t%s/cand1.vcf.gz\tp\tKID1,DAD1,MOM1\n' "$T"
  printf 'T2\t%s/cand2.vcf.gz\tp\tKID2,DAD2,MOM2\n' "$T"; } > "$W/trios.resolved.tsv"
{ printf 'trio_id\tcandidates_vcf\tped\n'
  printf 'T1\t%s/cand1.vcf.gz\tp\n' "$T"
  printf 'T2\t%s/cand2.vcf.gz\tp\n' "$T"; } > "$W/trios.candidates.tsv"
# BASE-FORM GTs, as Step 5 emits (cyvcf2 gt_bases). MOMn carries alt=T at chr1:100.
{ printf 'trio_id\tchrom\tpos\tref\talt\tmode\tchild_gt\tmother_gt\tfather_gt\n'
  printf 'T1\tchr1\t100\tA\tT\tdominant\tA/T\tA/T\tA/A\n'
  printf 'T1\tchr1\t200\tC\tG\tdominant\tC/G\tC/C\tC/C\n'
  printf 'T2\tchr1\t100\tA\tT\tdominant\tA/T\tA/T\tA/A\n'
  printf 'T2\tchr1\t200\tC\tG\tdominant\tC/G\tC/C\tC/C\n'; } > "$W/candidates.calls.tsv"
printf 'trio_id\n' > "$W/qc_report.tsv"

mkdir -p fakedb/taxonomy
: > fakedb/hash.k2d; : > fakedb/opts.k2d; : > fakedb/taxo.k2d
: > fakedb/taxonomy/nodes.dmp; : > fakedb/taxonomy/names.dmp

cat > bin/nonhuman-screen <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
[[ "${1:-}" == "--version" ]] && { echo "nonhuman-screen 0.0.0-stub"; exit 0; }
vcf="" outp="" bam=""
while [[ $# -gt 0 ]]; do case "$1" in
  classify) shift;; --bam) bam="$2"; shift 2;; --variants) vcf="$2"; shift 2;;
  --out-prefix) outp="$2"; shift 2;; --memory-mapping) shift;;
  --ref-fasta|--kraken2-db|--confidence|--threads) shift 2;; *) shift;; esac; done
# slow enough that two concurrent tasks genuinely overlap (the point of test 3)
sleep 1
# report the ACTUAL read count of the file we were handed, so the test can prove the
# classify-time filter reached nonhuman-screen (a fixed number could not).
nreads=$(samtools view -c "$bam" 2>/dev/null || echo 0)
printf 'variant_key\tsupporting_reads\tnonhuman_fraction\n' > "$outp.variant_nhf.tsv"
bcftools query -f '%CHROM:%POS0:%REF:%ALT\n' "$vcf" | grep -vE '[*<]' | \
  while IFS= read -r k; do printf '%s\t%s\t0.50\n' "$k" "$nreads" >> "$outp.variant_nhf.tsv"; done
echo '{}' > "$outp.summary.json"
STUB
chmod +x bin/nonhuman-screen

run8() { PATH="$T/bin:$PATH" HPRV_RUNTIME=native HPRV_TMPDIR="${TMPD:-$T/tmp}" \
    bash "$REPO/pipeline/08_igv_export.sh" --work "$W" --ref "$T/ref.fa" --cram-map "$T/map.tsv" \
        --kraken2-db "$T/fakedb" --nhf-members carriers --nhf-min-reads 5 "$@"; }

fail=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fail=1; fi; }

# --- serial reference run (also builds the mini-CRAMs + per-trio VCFs the array needs) ---
run8 >/dev/null 2>&1
cp "$W/igv/variants.tsv" "$T/serial.variants.tsv"
serial_md5="$(python3 -c 'import hashlib,sys;print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest())' "$T/serial.variants.tsv")"

# --- 1. manifest lists outstanding work; empty once everything is screened ---
run8 --nhf-emit-manifest "$T/m_done.txt" >/dev/null 2>&1
chk "manifest is empty when all trios are already screened" \
    '[[ ! -s "$T/m_done.txt" ]]'

rm -rf "$W/igv/nhf"
run8 --nhf-emit-manifest "$T/m_all.txt" >/dev/null 2>&1
chk "manifest lists both trios when none are screened" \
    '[[ "$(sort "$T/m_all.txt" | tr "\n" " ")" == "T1 T2 " ]]'

# --- 2. --nhf-trio screens ONLY that trio ---
run8 --nhf-trio T1 >/dev/null 2>&1
chk "--nhf-trio T1 screened only T1" \
    '[[ -d "$W/igv/nhf/T1" && ! -d "$W/igv/nhf/T2" ]]'

# --- 1b. manifest SHRINKS after partial work ---
run8 --nhf-emit-manifest "$T/m_left.txt" >/dev/null 2>&1
chk "manifest shrinks to the unscreened trio after a partial run" \
    '[[ "$(tr -d "[:space:]" < "$T/m_left.txt")" == "T2" ]]'

# --- 5. gather with an unscreened trio: warns, still writes, exits 0 ---
g_rc=0; run8 --nhf-gather > "$T/gather.log" 2>&1 || g_rc=$?
chk "gather with an unscreened member exits 0 (never dies)" '[[ "$g_rc" -eq 0 ]]'
chk "gather warns about members with no NHF row" 'grep -q "no NHF row" "$T/gather.log"'
chk "gather still wrote variants.tsv" '[[ -s "$W/igv/variants.tsv" ]]'

# --- 3 + 4. concurrent scatter: private carriers files, untouched VCFs, identical output ---
rm -rf "$W/igv/nhf"
vcf_md5_before="$(python3 -c 'import hashlib,sys;print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest())' "$W/igv/vcfs/T1.vcf.gz")"
for t in T1 T2; do
    ( mkdir -p "$T/tmp.$t"; TMPD="$T/tmp.$t" run8 --nhf-trio "$t" >/dev/null 2>&1 ) &
done
wait
vcf_md5_after="$(python3 -c 'import hashlib,sys;print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest())' "$W/igv/vcfs/T1.vcf.gz")"
chk "concurrent tasks did NOT rewrite the per-trio VCF (no torn cp -f / unstable key)" \
    '[[ "$vcf_md5_before" == "$vcf_md5_after" ]]'
chk "each concurrent task used its OWN carriers file (no shared truncation)" \
    '[[ -f "$T/tmp.T1/nhf_carriers.T1.tsv" && -f "$T/tmp.T2/nhf_carriers.T2.tsv" ]]'
chk "both trios were screened concurrently" \
    '[[ -d "$W/igv/nhf/T1" && -d "$W/igv/nhf/T2" ]]'

run8 --nhf-gather >/dev/null 2>&1
scatter_md5="$(python3 -c 'import hashlib,sys;print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest())' "$W/igv/variants.tsv")"
chk "scattered variants.tsv is BYTE-IDENTICAL to the serial run" \
    '[[ "$serial_md5" == "$scatter_md5" ]]'

# --- duplicate exclusion: the mini-CRAM must carry NO duplicate-flagged reads ---
# --- the two masks are INDEPENDENT: archive complete, NHF denominator filtered ---
chk "fixture really contains duplicate-flagged reads (guards the test itself)" \
    '[[ "$dups_in_source" -eq 3 ]]'
# The mini-CRAM is the ARCHIVAL review track: the classify-time mask must not touch it.
chk "mini-CRAM ARCHIVES duplicate-flagged reads (slice default 0 = keep everything)" \
    '[[ "$(samtools view -c -f 1024 "$W/igv/crams/T1/KID1.cram")" -eq "$dups_in_source" ]]'
chk "mini-CRAM archives supplementary alignments (GATK HC calls on them)" \
    '[[ "$(samtools view -c -f 2048 "$W/igv/crams/T1/KID1.cram")" -eq "$supp_in_source" ]]'
mini_total="$(samtools view -c "$W/igv/crams/T1/KID1.cram")"
# ...but the NHF denominator must reflect the CLASSIFY-time filter, i.e. fewer reads than the archive.
nhf_reads="$(awk -F'\t' 'NR==2{print $2}' "$W/igv/nhf/T1/KID1.variant_nhf.tsv")"
chk "NHF read count is LOWER than the archived mini-CRAM (duplicates filtered at classify time)" \
    '[[ "$nhf_reads" -lt "$mini_total" ]]'
chk "NHF read count equals the archive minus the duplicates (supplementary kept)" \
    '[[ "$nhf_reads" -eq "$(( mini_total - dups_in_source ))" ]]'
# Keys: the classify mask belongs to nhf_key, NOT to the mini-CRAM key — so revising it must
# recompute NHF without forcing a re-slice off the (expensive, flaky) source CRAM mount.
chk "mini-CRAM .done key does NOT carry the classify-time mask" \
    '! grep -q -- "-F1796" "$W/igv/crams/T1/KID1.cram.done"'
chk "NHF .done key DOES carry the classify-time mask" \
    'grep -q -- "-F1796" "$W/igv/nhf/T1/KID1.variant_nhf.tsv.done"'
# changing the classify mask recomputes NHF and leaves the mini-CRAM byte-identical
cram_md5_before="$(python3 -c 'import hashlib,sys;print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest())' "$W/igv/crams/T1/KID1.cram")"
run8 --nhf-exclude-flags 0 >/dev/null 2>&1
cram_md5_after="$(python3 -c 'import hashlib,sys;print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest())' "$W/igv/crams/T1/KID1.cram")"
nhf_reads_unfiltered="$(awk -F'\t' 'NR==2{print $2}' "$W/igv/nhf/T1/KID1.variant_nhf.tsv")"
chk "changing the classify mask did NOT re-slice the mini-CRAM" \
    '[[ "$cram_md5_before" == "$cram_md5_after" ]]'
chk "changing the classify mask DID recompute NHF (now counts every archived read)" \
    '[[ "$nhf_reads_unfiltered" -eq "$mini_total" ]]'

# --- MIGRATION: which cached mini-CRAMs must be re-sliced off the source mount, and which must not.
# Re-slicing is the most expensive operation in the pipeline (source CRAM over a flaky FUSE/SBFS
# mount), so a cache that is ALREADY correct must survive the upgrade untouched.
cram="$W/igv/crams/T1/KID1.cram"
md5_now="$(python3 -c 'import hashlib,sys;print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest())' "$cram")"

# (a) LEGACY cache: sliced before this feature existed, so its .done holds a bare bed key and the
#     archive is COMPLETE (nothing was filtered). It must be REUSED, not re-fetched.
legacy_key="$(cksum < "$T/tmp/T1.merged.bed" 2>/dev/null | awk '{print $1"-"$2}')"
if [[ -n "$legacy_key" ]]; then
    printf '%s\n' "$legacy_key" > "$cram.done"
    run8 >/dev/null 2>&1
    chk "a LEGACY-keyed mini-CRAM (no -F suffix) is reused, not re-sliced" \
        '[[ "$(python3 -c "import hashlib,sys;print(hashlib.md5(open(sys.argv[1],\"rb\").read()).hexdigest())" "$cram")" == "$md5_now" \
            && "$(cat "$cram.done")" == "$legacy_key" ]]'
else
    echo "SKIP legacy-key reuse (merged BED not retained in tmp)"
fi

# (b) A cache sliced under a NON-ZERO mask is genuinely INCOMPLETE (reads were dropped), so it must
#     be re-sliced under the new archive-everything default.
printf '%s\n' "deadbeef-1-F1796" > "$cram.done"
run8 >/dev/null 2>&1
chk "a mini-CRAM cached under a non-zero -F mask IS re-sliced (its archive was incomplete)" \
    '[[ "$(cat "$cram.done")" != "deadbeef-1-F1796" ]]'
chk "the re-sliced mini-CRAM again archives the duplicates" \
    '[[ "$(samtools view -c -f 1024 "$cram")" -eq "$dups_in_source" ]]'

[[ "$fail" -eq 0 ]] && echo "All Step-8b scatter/gather tests passed." \
    || { echo "test_nhf_scatter FAILED"; exit 1; }
