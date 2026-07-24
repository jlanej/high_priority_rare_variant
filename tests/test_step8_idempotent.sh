#!/usr/bin/env bash
# =============================================================================
# tests/test_step8_idempotent.sh — Step-8 mini-CRAM slicing idempotency.
#
# Runs the REAL pipeline/08_igv_export.sh (no logic copy) and asserts the resumable
# slicing contract, which matters on a flaky FUSE/SBFS mount where re-reading multi-GB
# source CRAMs is the cost to avoid:
#
#   A) RESUME REUSES a valid prior slice WITHOUT re-reading the source CRAM. Proven by
#      corrupting the source between runs: if the slice is reused the source is never
#      touched and the mini-CRAM is byte-identical; if it re-sliced it would read the
#      garbage source, fail, and drop the mini-CRAM.
#   B) A CHANGED candidate set (different merged BED -> different content key) RE-SLICES,
#      not serve a stale mini-CRAM.
#   C) A CORRUPTED mini-CRAM (quickcheck fails) is REGENERATED.
#
# Needs bcftools/samtools/bedtools/bgzip/tabix/cksum + python3 (stdlib); self-skips if absent.
# Run: bash tests/test_step8_idempotent.sh
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

for t in bcftools samtools bedtools bgzip tabix cksum python3; do
    command -v "$t" >/dev/null 2>&1 || { echo "SKIP test_step8_idempotent: $t not on PATH"; exit 0; }
done

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
W="$T/work"; mkdir -p "$W"
cd "$T"

python3 -c 'open("ref.fa","w").write(">chr1\n"+("ACGTACGTAC"*80)+"\n")'
samtools faidx ref.fa

mk_cram() {  # $1=sample id; reads spanning chr1:100, :300, :500
    local s="$1" seq qual
    seq="$(python3 -c 'print("A"*40)')"; qual="$(python3 -c 'print(chr(73)*40)')"
    { printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:800\n@RG\tID:%s\tSM:%s\n' "$s" "$s"
      for p in 90 290 490; do
          printf '%s_%s\t0\tchr1\t%s\t60\t40M\t*\t0\t0\t%s\t%s\tRG:Z:%s\n' "$s" "$p" "$p" "$seq" "$qual" "$s"
      done
    } > "$s.sam"
    samtools view -C -T ref.fa -o "$s.cram" "$s.sam"; samtools index "$s.cram"
}
mk_cram KID
printf 'KID\t%s/KID.cram\n' "$T" > map.tsv

# per-trio candidate VCF + manifests (single trio, child only screened for CRAM tracks)
printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=800>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tKID\nchr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\n' > cand.vcf
bgzip -f cand.vcf; bcftools index -t cand.vcf.gz
printf 'trio_id\tvcf\tped\tsamples\nT1\t%s/cand.vcf.gz\tp\tKID,KID,KID\n' "$T" > "$W/trios.resolved.tsv"
printf 'trio_id\tcandidates_vcf\tped\nT1\t%s/cand.vcf.gz\tp\n' "$T" > "$W/trios.candidates.tsv"
write_calls() {  # $1 = candidate POS (drives the merged BED)
    { printf 'trio_id\tchrom\tpos\tref\talt\tmode\tchild_gt\n'
      printf 'T1\tchr1\t%s\tA\tT\tdominant\t0/1\n' "$1"; } > "$W/candidates.calls.tsv"
}
write_calls 100
printf 'trio_id\n' > "$W/qc_report.tsv"

run8() { HPRV_RUNTIME=native HPRV_TMPDIR="$T/tmp" bash "$REPO/pipeline/08_igv_export.sh" \
             --work "$W" --ref "$T/ref.fa" --cram-map "$T/map.tsv" >/dev/null 2>&1; }

ocram="$W/igv/crams/T1/KID.cram"
fail=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fail=1; fi; }
ck() { cksum < "$1" 2>/dev/null | awk '{print $1"-"$2}'; }

# --- Run 1: produce the mini-CRAM ---
run8
chk "run1: mini-CRAM sliced"          '[[ -s "$ocram" ]]'
chk "run1: .done sentinel written"    '[[ -f "$ocram.done" ]]'
chk "run1: mini-CRAM validates"       'samtools quickcheck "$ocram"'
c1="$(ck "$ocram")"; key1="$(cat "$ocram.done")"

# --- A) RESUME with a CORRUPT source: must reuse the slice, never re-read the source ---
cp "$T/KID.cram" "$T/KID.cram.bak"
printf 'GARBAGE-not-a-cram' > "$T/KID.cram"      # file still exists -> Pass 1 emits the task
run8
chk "A: export still succeeds with a corrupt source" '[[ -s "$ocram" ]]'
chk "A: mini-CRAM byte-identical (source NOT re-read)" '[[ "$(ck "$ocram")" == "$c1" ]]'
chk "A: mini-CRAM still validates"    'samtools quickcheck "$ocram"'
chk "A: variants.tsv child_file populated" 'python3 -c "import csv,sys; r=next(csv.DictReader(open(sys.argv[1]),delimiter=chr(9))); sys.exit(0 if r[\"child_file\"] else 1)" "$W/igv/variants.tsv"'
mv -f "$T/KID.cram.bak" "$T/KID.cram"             # restore a valid source

# --- B) CHANGED candidate set (different BED -> different key): must RE-SLICE ---
write_calls 300                                    # candidate moved chr1:100 -> chr1:300
run8
key2="$(cat "$ocram.done")"
chk "B: region key changed on a changed candidate set" '[[ "$key2" != "$key1" ]]'
chk "B: re-sliced mini-CRAM validates" 'samtools quickcheck "$ocram"'

# --- C) CORRUPTED mini-CRAM (quickcheck fails) is REGENERATED ---
printf 'CORRUPT' > "$ocram"                        # leave .done in place (key still matches BED)
run8
chk "C: corrupted mini-CRAM regenerated + validates" 'samtools quickcheck "$ocram"'

[[ "$fail" -eq 0 ]] && echo "All Step-8 idempotency tests passed." || { echo "test_step8_idempotent FAILED"; exit 1; }
