#!/usr/bin/env bash
# =============================================================================
# assert_shard_equivalence.sh — Step 2's by-contig VEP sharding must produce
# output IDENTICAL to a single un-sharded pass.
#
# This is load-bearing: Step 2 feeds the published call set, so "sharded == single"
# is a correctness guarantee, not a nicety. The test exercises Step 2's REAL
# (non-ingest) code path by putting a `vep` shim on PATH (mock_vep.py in shim mode),
# then runs Step 2 twice on the same multi-contig union — once sharded, once not —
# and asserts the record bodies and the split-vep INFO headers are byte-identical.
#
# Usage: assert_shard_equivalence.sh <cohort.sites.vcf.gz> <annot.tsv> <reference.fa>
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"

SITES="$1"; LOOKUP="$2"; REF="$3"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# A `vep` shim on PATH -> Step 2's real path runs and calls it once per contig shard.
shimdir="$WORK/shim"; mkdir -p "$shimdir"
cat > "$shimdir/vep" <<EOF
#!/usr/bin/env bash
exec python3 "$HERE/mock_vep.py" "\$@"
EOF
chmod +x "$shimdir/vep"

# Step 2 requires an existing VEP cache dir (path check only; the shim ignores it).
mkdir -p "$WORK/cache/homo_sapiens"

run_step2() {  # $1 = shard flag (0|1), $2 = output path
    PATH="$shimdir:$PATH" \
    HPRV_ENGINE=native \
    HPRV_MOCK_VEP_LOOKUP="$LOOKUP" \
    HPRV_VEP_CACHE="$WORK/cache" \
    HPRV_TMPDIR="$WORK/tmp.$1" \
    HPRV_VEP_SHARD_BY_CONTIG="$1" \
        bash "$REPO/pipeline/02_annotate_sites.sh" --sites "$SITES" --ref "$REF" --out "$2" >/dev/null 2>&1
}

echo "== shard-equivalence: annotating the union sharded, then single =="
run_step2 1 "$WORK/sharded.vcf.gz"
run_step2 0 "$WORK/single.vcf.gz"

fail=0

# The union must actually span >=2 contigs, or the test proves nothing.
ncontig="$(bcftools index -s "$SITES" | awk '$3>0' | wc -l | tr -d ' ')"
if [[ "${ncontig:-0}" -ge 2 ]]; then
    echo "PASS union spans $ncontig contigs (sharding was genuinely exercised)"
else
    echo "FAIL union has <2 contigs with variants — equivalence test is vacuous"; fail=1
fi

# Record bodies (CHROM..INFO incl. CSQ + every vep_* field) must be byte-identical:
# same sites, same annotation, same order.
if diff <(bcftools view -H "$WORK/sharded.vcf.gz") <(bcftools view -H "$WORK/single.vcf.gz") >/dev/null; then
    echo "PASS shard==single: record bodies byte-identical"
else
    echo "FAIL shard==single: record bodies DIFFER"
    diff <(bcftools view -H "$WORK/sharded.vcf.gz") <(bcftools view -H "$WORK/single.vcf.gz") | head -10
    fail=1
fi

# The INFO definitions split-vep created (the vep_* fields) must match too.
if diff <(bcftools view -h "$WORK/sharded.vcf.gz" | grep '^##INFO' | sort) \
        <(bcftools view -h "$WORK/single.vcf.gz" | grep '^##INFO' | sort) >/dev/null; then
    echo "PASS shard==single: INFO header definitions identical"
else
    echo "FAIL shard==single: INFO header definitions DIFFER"; fail=1
fi

[[ "$fail" -eq 0 ]] && echo "ALL SHARD-EQUIVALENCE ASSERTIONS PASSED" || { echo "SHARD-EQUIVALENCE FAILED"; exit 1; }
