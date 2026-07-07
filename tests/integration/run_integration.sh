#!/usr/bin/env bash
# =============================================================================
# run_integration.sh — end-to-end pipeline test over generated mock data.
#
# Exercises resolve + Steps 0,1,3,4,5,6 for real (bcftools + the python steps).
# Only Step 2's VEP invocation is mocked (VEP cache is too heavy for CI) — the
# gnomAD/ClinVar transfer runs real `bcftools annotate`.
#
# Requires on PATH: bcftools, samtools, bgzip, tabix, and a python with
# cyvcf2/pysam/numpy/scipy/pyyaml. Runs on a laptop (conda env) or in CI.
#
# Usage: run_integration.sh [WORK_ROOT]   (defaults to a temp dir)
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"

W="${1:-$(mktemp -d)}"
mkdir -p "$W"
echo "== mock data root: $W =="

python3 "$HERE/make_mock_data.py" --out "$W"
samtools faidx "$W/reference.fa"
for f in "$W"/vcfs/*.vcf; do bgzip -f "$f"; tabix -f -p vcf "$f.gz"; done
# source SAMs -> sorted+indexed CRAMs (for Step 8 mini-CRAM slicing)
for s in "$W"/crams_src/*.sam; do
    b="${s%.sam}"
    samtools sort -O cram --reference "$W/reference.fa" -o "$b.cram" "$s"
    samtools index "$b.cram"
done
bgzip -f "$W/gnomad.sites.vcf"; tabix -f -p vcf "$W/gnomad.sites.vcf.gz"
bgzip -f "$W/clinvar.vcf";      tabix -f -p vcf "$W/clinvar.vcf.gz"

CFG="$W/config.mock.yaml"
WORK="$W/work"

# --- resolve + Step 0 (QC) + Step 1 (cohort sites) ---
bash "$REPO/pipeline/run_pipeline.sh" --config "$CFG" --from 0 --to 1

# --- mock Step 2: VEP fields + real gnomAD/ClinVar transfer ---
python3 "$HERE/mock_annotate.py" --in "$WORK/cohort.sites.vcf.gz" \
    --lookup "$W/annot.tsv" --out "$WORK/cohort.sites.vep.vcf"
bgzip -f "$WORK/cohort.sites.vep.vcf"; tabix -f -p vcf "$WORK/cohort.sites.vep.vcf.gz"
bcftools annotate -a "$W/gnomad.sites.vcf.gz" \
    -c 'INFO/hprv_gnomad_af:=INFO/AF,INFO/hprv_gnomad_grpmax_af:=INFO/AF_grpmax,INFO/hprv_gnomad_faf95:=INFO/faf95,INFO/hprv_gnomad_nhomalt:=INFO/nhomalt' \
    -Oz -o "$WORK/cohort.sites.gn.vcf.gz" "$WORK/cohort.sites.vep.vcf.gz"
tabix -f -p vcf "$WORK/cohort.sites.gn.vcf.gz"
bcftools annotate -a "$W/clinvar.vcf.gz" \
    -c 'INFO/hprv_clnsig:=INFO/CLNSIG,INFO/hprv_clnrevstat:=INFO/CLNREVSTAT' \
    -Oz -o "$WORK/cohort.sites.annotated.vcf.gz" "$WORK/cohort.sites.gn.vcf.gz"
tabix -f -p vcf "$WORK/cohort.sites.annotated.vcf.gz"
touch "$WORK/cohort.sites.annotated.vcf.gz.done"

# --- Steps 3-8 (includes xlsx summary + igv export) ---
bash "$REPO/pipeline/run_pipeline.sh" --config "$CFG" --from 3 --to 8

# --- assertions ---
python3 "$HERE/assert_integration.py" --work "$WORK"
echo "== integration run complete: $WORK =="
