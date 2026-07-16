#!/usr/bin/env bash
# =============================================================================
# run_integration.sh — end-to-end pipeline test over generated mock data.
#
# Exercises resolve + Steps 0-8 for real (bcftools + the python steps). ONLY the `vep`
# binary itself is mocked (its cache is 24 GB — too heavy for CI): mock_vep.py writes a
# VEP-shaped CSQ, and Step 2 then runs for real via its --vep-vcf ingest path, so the
# build checks, the split-vep lift, the transcript selector and the frequency guard are
# all under test. There is no gnomAD/ClinVar transfer to exercise — under the VEP-only
# contract those annotations come from the CSQ.
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
CFG="$W/config.mock.yaml"
WORK="$W/work"

# --- resolve + Step 0 (QC) + Step 1 (cohort sites) ---
bash "$REPO/pipeline/run_pipeline.sh" --config "$CFG" --from 0 --to 1

# --- stand in for the `vep` binary only: produce a VEP-shaped CSQ over the cohort union ---
python3 "$HERE/mock_vep.py" --in "$WORK/cohort.sites.vcf.gz" \
    --lookup "$W/annot.tsv" --out "$W/cohort.sites.vep.vcf"
bgzip -f "$W/cohort.sites.vep.vcf"; tabix -f -p vcf "$W/cohort.sites.vep.vcf.gz"

# --- Steps 2-8 for real. Step 2 ingests the VEP VCF above (resources.vep.annotated_vcf in the
#     mock config), so its build checks + split-vep + selector + frequency guard all execute. ---
bash "$REPO/pipeline/run_pipeline.sh" --config "$CFG" --from 2 --to 8

# --- assertions ---
python3 "$HERE/assert_integration.py" --work "$WORK"
echo "== integration run complete: $WORK =="
