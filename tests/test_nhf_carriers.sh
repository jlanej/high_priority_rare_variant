#!/usr/bin/env bash
# =============================================================================
# tests/test_nhf_carriers.sh — Step-8b `members: carriers` regression test.
#
# Runs the REAL pipeline/08_igv_export.sh (no copy of its awk, so no drift) with a
# stubbed `nonhuman-screen` on PATH, over a fixture whose calls.tsv uses BASE-FORM
# genotypes ('A/T', 'A/A', './.') exactly as Step 5 writes them (cyvcf2 gt_bases).
#
# Guards the defect where the carrier awk tested `mother_gt ~ /[1-9]/` — which never
# matches base-form GTs, silently degrading `carriers` to `child_only`. Asserts:
#   * the child (always) AND the carrier mother (mother_gt=A/T at an alt=T call) ARE screened,
#   * the NON-carrier father (father_gt=A/A everywhere) is NOT screened,
#   * variants.tsv reflects that (mother_nhf populated on the shared locus; father_nhf blank).
#
# Needs bcftools/samtools/bgzip/tabix + a python3 (stdlib only); self-skips if absent.
# Run: bash tests/test_nhf_carriers.sh
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

for t in bcftools samtools bgzip tabix python3; do
    command -v "$t" >/dev/null 2>&1 || { echo "SKIP test_nhf_carriers: $t not on PATH"; exit 0; }
done

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
W="$T/work"; mkdir -p "$W" "$T/bin"
cd "$T"

python3 -c 'open("ref.fa","w").write(">chr1\n"+("ACGTACGTAC"*60)+"\n")'
samtools faidx ref.fa

mk_cram() {  # $1 = sample id; a few reads over chr1:100 and chr1:200
    local s="$1" seq qual
    seq="$(python3 -c 'print("A"*40)')"; qual="$(python3 -c 'print(chr(73)*40)')"
    { printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:600\n@RG\tID:%s\tSM:%s\n' "$s" "$s"
      for p in 90 190; do
          printf '%s_%s\t0\tchr1\t%s\t60\t40M\t*\t0\t0\t%s\t%s\tRG:Z:%s\n' "$s" "$p" "$p" "$seq" "$qual" "$s"
      done
    } > "$s.sam"
    samtools view -C -T ref.fa -o "$s.cram" "$s.sam"; samtools index "$s.cram"
}
mk_cram KID; mk_cram DAD; mk_cram MOM
printf 'KID\t%s/KID.cram\nDAD\t%s/DAD.cram\nMOM\t%s/MOM.cram\n' "$T" "$T" "$T" > map.tsv

# per-trio candidate VCF (2 biallelic variants)
cat > cand.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##FORMAT=<ID=GT,Number=1,Type=String,Description="g">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	KID	DAD	MOM
chr1	100	.	A	T	.	.	.	GT	0/1	0/0	0/1
chr1	200	.	C	G	.	.	.	GT	0/1	0/0	0/0
EOF
bgzip -f cand.vcf; bcftools index -t cand.vcf.gz

printf 'trio_id\tvcf\tped\tsamples\nT1\t%s/cand.vcf.gz\tp\tKID,DAD,MOM\n' "$T" > "$W/trios.resolved.tsv"
printf 'trio_id\tcandidates_vcf\tped\nT1\t%s/cand.vcf.gz\tp\n' "$T" > "$W/trios.candidates.tsv"
# calls.tsv with BASE-FORM GTs (as Step 5 emits): MOM carries alt=T at chr1:100; DAD never carries.
{ printf 'trio_id\tchrom\tpos\tref\talt\tmode\tchild_gt\tmother_gt\tfather_gt\n'
  printf 'T1\tchr1\t100\tA\tT\tdominant\tA/T\tA/T\tA/A\n'
  printf 'T1\tchr1\t200\tC\tG\tdominant\tC/G\tC/C\tC/C\n'; } > "$W/candidates.calls.tsv"
printf 'trio_id\n' > "$W/qc_report.tsv"

mkdir -p fakedb/taxonomy
: > fakedb/hash.k2d; : > fakedb/opts.k2d; : > fakedb/taxo.k2d
: > fakedb/taxonomy/nodes.dmp; : > fakedb/taxonomy/names.dmp

# stub nonhuman-screen: emits 0-based keys via bcftools %POS0, uniform per-sample fraction.
cat > bin/nonhuman-screen <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
[[ "${1:-}" == "--version" ]] && { echo "nonhuman-screen 0.0.0-stub"; exit 0; }
bam="" vcf="" outp=""
while [[ $# -gt 0 ]]; do case "$1" in
  classify) shift;; --bam) bam="$2"; shift 2;; --variants) vcf="$2"; shift 2;;
  --out-prefix) outp="$2"; shift 2;; --memory-mapping) shift;;
  --ref-fasta|--kraken2-db|--confidence|--threads) shift 2;; *) shift;; esac; done
printf 'variant_key\tsupporting_reads\tnonhuman_fraction\n' > "$outp.variant_nhf.tsv"
bcftools query -f '%CHROM:%POS0:%REF:%ALT\n' "$vcf" | grep -vE '[*<]' | \
  while IFS= read -r k; do printf '%s\t7\t0.50\n' "$k" >> "$outp.variant_nhf.tsv"; done
echo '{}' > "$outp.summary.json"
STUB
chmod +x bin/nonhuman-screen

PATH="$T/bin:$PATH" HPRV_RUNTIME=native HPRV_TMPDIR="$T/tmp" \
    bash "$REPO/pipeline/08_igv_export.sh" --work "$W" --ref "$T/ref.fa" --cram-map "$T/map.tsv" \
        --kraken2-db "$T/fakedb" --nhf-members carriers --nhf-min-reads 5 --nhf-memory-mapping \
        >/dev/null 2>&1

fail=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fail=1; fi; }

chk "child (KID) screened"                 '[[ -f "$W/igv/nhf/T1/KID.variant_nhf.tsv" ]]'
chk "carrier mother (MOM, A/T) screened"    '[[ -f "$W/igv/nhf/T1/MOM.variant_nhf.tsv" ]]'
chk "non-carrier father (DAD, A/A) NOT screened" '[[ ! -f "$W/igv/nhf/T1/DAD.variant_nhf.tsv" ]]'

python3 - "$W/igv/variants.tsv" <<'PY' || fail=1
import csv, sys
rows = {(r["chrom"], r["pos"]): r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
r = rows[("chr1", "100")]
ok = (r["child_nhf"] == "0.50" and r["mother_nhf"] == "0.50" and r["father_nhf"] == ""
      and r["nhf_flag"] == "1")
print(("PASS" if ok else "FAIL") + " variants.tsv join (child+mother populated, father blank, flag=1)")
sys.exit(0 if ok else 1)
PY

[[ "$fail" -eq 0 ]] && echo "All Step-8b carrier tests passed." || { echo "test_nhf_carriers FAILED"; exit 1; }
