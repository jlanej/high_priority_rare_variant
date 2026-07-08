#!/usr/bin/env python3
"""Generate a tiny, self-consistent mock dataset exercising the whole pipeline.

Writes under --out:
  reference.fa                mini GRCh38-like genome (chr1, chr2, chrX-with-nonPAR)
  vcfs/fileA.vcf              trio A (CH_A/FA_A/MO_A) — autosomal modes + filter cases
  vcfs/fileB.vcf              FAMILY VCF: [MO_B, SIB_B, CH_B, FA_B] (extra sibling,
                              shuffled order) — de novo (recurrent gene) + X-linked
  vcfs/fileC.vcf              duo CH_C/FA_C (mom absent) — resolver "unresolved" case
  gnomad.sites.vcf            external AF (faf95/AF/AF_grpmax/nhomalt) for some sites
  clinvar.vcf                 CLNSIG/CLNREVSTAT for a P/LP site
  trios.tsv                   #kid dad mom (A, B resolvable; C unresolvable)
  annot.tsv                   per-site VEP-like lookup for mock_annotate.py
  mutrate.tsv, constraint.tsv gene tables for Step 6
  config.mock.yaml            config pointing at the above (concrete paths, ephemeral)

All variants are SNVs (no indel left-align ambiguity). REF matches the reference at
each position. Not committed data — generated into the git-ignored work dir.
"""
from __future__ import annotations

import argparse
import os

CONTIGS = {"chr1": 20000, "chr2": 20000, "chrX": 2782200}  # chrX > PAR1 end (2,781,479)
BASES = "ACGT"


def refbase(pos):
    return BASES[pos % 4]


def altbase(pos):
    return BASES[(pos + 1) % 4]


def ad(gt, dp):
    if gt == "0/0":
        return f"{dp},0"
    if gt == "1/1":
        return f"0,{dp}"
    h = dp // 2
    return f"{dp - h},{h}"


# Sample groupings per VCF file (note fileB order is shuffled and has an extra sib).
FILES = {
    "A": ["CH_A", "FA_A", "MO_A"],
    "B": ["MO_B", "SIB_B", "CH_B", "FA_B"],
    "C": ["CH_C", "FA_C"],
}

# Each variant: file, chrom, pos, gene, csq, impact, revel, am, spliceai, loftee,
# faf95 (None=absent from gnomAD), clnsig, clnrevstat, filter, hidenovo, gts.
V = []


def add(**k):
    k.setdefault("revel", ""); k.setdefault("am", ""); k.setdefault("spliceai", "")
    k.setdefault("loftee", ""); k.setdefault("faf95", None); k.setdefault("clnsig", "")
    k.setdefault("clnrevstat", ""); k.setdefault("filter", "PASS"); k.setdefault("hidenovo", "")
    V.append(k)


# --- Trio A (autosomal), CH_A female ---------------------------------------
# 1) de novo, HIGH LoF, absent -> expect mode=denovo
add(file="A", chrom="chr1", pos=5000, gene="GENE1", csq="stop_gained", impact="HIGH",
    loftee="HC", hidenovo="CH_A",
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/0", 99, 40)})
# 2) homozygous recessive, MODERATE missense, rare -> mode=hom_recessive
add(file="A", chrom="chr1", pos=8000, gene="GENE2", csq="missense_variant", impact="MODERATE",
    revel="0.95", faf95=5e-4,
    gts={"CH_A": ("1/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/1", 99, 40)})
# 3+4) compound het in GENE3 (var3 maternal, var4 paternal) -> mode=compound_het
add(file="A", chrom="chr2", pos=5000, gene="GENE3", csq="missense_variant", impact="MODERATE",
    revel="0.80", faf95=1e-3,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/1", 99, 40)})
add(file="A", chrom="chr2", pos=6000, gene="GENE3", csq="missense_variant", impact="MODERATE",
    revel="0.80", faf95=1e-3,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# 5) common (BA1) -> dropped at Step 3 (never a candidate)
add(file="A", chrom="chr1", pos=12000, gene="GENE4", csq="missense_variant", impact="MODERATE",
    revel="0.70", faf95=0.2,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# 6) low-GQ de-novo-looking (no hiConfDeNovo) -> passes Step 3, FAILS Step 5 QC
add(file="A", chrom="chr1", pos=15000, gene="GENE5", csq="stop_gained", impact="HIGH", loftee="HC",
    gts={"CH_A": ("0/1", 12, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/0", 99, 40)})
# 7) non-PASS -> dropped at Step 1
add(file="A", chrom="chr1", pos=17000, gene="GENE6", csq="stop_gained", impact="HIGH", loftee="HC",
    filter="VQSRTrancheSNP99.00to99.90+", hidenovo="CH_A",
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/0", 99, 40)})
# 8) ClinVar P/LP but LOW impact (synonymous) -> kept at Step 3 via clinvar_plp override
add(file="A", chrom="chr2", pos=8000, gene="GENE7", csq="synonymous_variant", impact="LOW",
    faf95=2e-4, clnsig="Pathogenic",  # > dominant_max: kept via ClinVar but not a dominant call
    clnrevstat="criteria_provided,_multiple_submitters,_no_conflicts",
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- Trio B (CH_B male), family VCF with extra sibling ----------------------
# 9) de novo in GENE1 (recurrent across trios A+B) -> mode=denovo; drives Step-6 burden
add(file="B", chrom="chr1", pos=5100, gene="GENE1", csq="stop_gained", impact="HIGH",
    loftee="HC", hidenovo="CH_B",
    gts={"CH_B": ("0/1", 99, 40), "FA_B": ("0/0", 99, 40), "MO_B": ("0/0", 99, 40),
         "SIB_B": ("0/0", 99, 40)})
# 10) X-linked recessive (male hemizygous), carrier mother -> mode=x_linked_recessive
add(file="B", chrom="chrX", pos=2781600, gene="GENEX", csq="missense_variant", impact="MODERATE",
    revel="0.90", faf95=1e-4,
    gts={"CH_B": ("1/1", 99, 40), "FA_B": ("0/0", 99, 40), "MO_B": ("0/1", 99, 40),
         "SIB_B": ("0/0", 99, 40)})
# 11-13) chrX filler (common) so CH_B is inferred MALE (hemizygous alt -> low het ratio)
for i, p in enumerate((2781700, 2781800, 2781900)):
    add(file="B", chrom="chrX", pos=p, gene=f"XFILL{i}", csq="missense_variant", impact="MODERATE",
        faf95=0.3,
        gts={"CH_B": ("1/1", 99, 40), "FA_B": ("1/1", 99, 40), "MO_B": ("0/1", 99, 40),
             "SIB_B": ("0/1", 99, 40)})

# --- DOMINANT RECURRENCE: same rare functional inherited het in GENED, in BOTH trios
#     (CH_A inherits from dad, CH_B from mom) -> Step 6 nominates GENED (n_dominant=2) ---
add(file="A", chrom="chr2", pos=10000, gene="GENED", csq="missense_variant", impact="MODERATE",
    revel="0.90", faf95=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
add(file="B", chrom="chr2", pos=10000, gene="GENED", csq="missense_variant", impact="MODERATE",
    revel="0.90", faf95=5e-5,
    gts={"CH_B": ("0/1", 99, 40), "MO_B": ("0/1", 99, 40), "FA_B": ("0/0", 99, 40),
         "SIB_B": ("0/0", 99, 40)})

# --- Trio C: duo only (mom MO_C absent everywhere) -> resolver unresolved ---
add(file="C", chrom="chr1", pos=9000, gene="GENE8", csq="missense_variant", impact="MODERATE",
    faf95=1e-3, gts={"CH_C": ("0/1", 99, 40), "FA_C": ("0/1", 99, 40)})


def write_reference(path):
    with open(path, "w") as fh:
        for c, n in CONTIGS.items():
            seq = bytearray(b"A" * n)
            for v in V:
                if v["chrom"] == c:
                    seq[v["pos"] - 1] = ord(refbase(v["pos"]))
            fh.write(f">{c}\n")
            s = seq.decode()
            for i in range(0, n, 60):
                fh.write(s[i:i + 60] + "\n")


VCF_HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=lowGQ,Description="GQ < 20.0">
##FILTER=<ID=VQSRTrancheSNP99.00to99.90+,Description="VQSR tranche">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=hiConfDeNovo,Number=1,Type=String,Description="High-confidence de novo child list">
##INFO=<ID=loConfDeNovo,Number=1,Type=String,Description="Low-confidence de novo child list">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
"""


def write_vcf(path, samples, variants):
    with open(path, "w") as fh:
        fh.write(VCF_HEADER)
        for c, n in CONTIGS.items():
            fh.write(f"##contig=<ID={c},length={n}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for v in sorted(variants, key=lambda x: (list(CONTIGS).index(x["chrom"]), x["pos"])):
            ref, alt = refbase(v["pos"]), altbase(v["pos"])
            ac = sum(g[0].count("1") for g in v["gts"].values())
            an = 2 * len(samples)
            info = f"AC={ac};AN={an};AF={ac / an:.4g}"
            if v["hidenovo"]:
                info += f";hiConfDeNovo={v['hidenovo']}"
            cells = []
            for s in samples:
                gt, gq, dp = v["gts"][s]
                cells.append(f"{gt}:{ad(gt, dp)}:{dp}:{gq}")
            fh.write(f"{v['chrom']}\t{v['pos']}\t.\t{ref}\t{alt}\t100\t{v['filter']}\t"
                     f"{info}\tGT:AD:DP:GQ\t" + "\t".join(cells) + "\n")


def write_sites(path, header_info, want):
    """Write a sites-only VCF with the given INFO records (want: list of (v, infostr))."""
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for line in header_info:
            fh.write(line + "\n")
        for c, n in CONTIGS.items():
            fh.write(f"##contig=<ID={c},length={n}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        seen = set()
        for v, info in sorted(want, key=lambda x: (list(CONTIGS).index(x[0]["chrom"]), x[0]["pos"])):
            key = (v["chrom"], v["pos"])
            if key in seen:
                continue
            seen.add(key)
            ref, alt = refbase(v["pos"]), altbase(v["pos"])
            fh.write(f"{v['chrom']}\t{v['pos']}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")


def write_source_sams(W):
    """Write a tiny per-sample SAM (a few reads over each of the sample's variants,
    coordinate-sorted) and return {sample: cram_path}. run_integration.sh converts each
    SAM to a sorted+indexed CRAM so Step 8 can slice mini-CRAMs from it."""
    sam_dir = os.path.join(W, "crams_src")
    os.makedirs(sam_dir, exist_ok=True)
    L = 60
    per_sample = {}  # sample -> list of (contig_idx, pos, samline)
    counter = 0
    for fk, samples in FILES.items():
        for v in [x for x in V if x["file"] == fk]:
            cidx = list(CONTIGS).index(v["chrom"])
            start = max(1, v["pos"] - 30)
            for s in samples:
                for _ in range(4):
                    counter += 1
                    line = (f"r{counter}\t0\t{v['chrom']}\t{start}\t60\t{L}M\t*\t0\t0\t"
                            f"{'A' * L}\t{'I' * L}")
                    per_sample.setdefault(s, []).append((cidx, start, line))
    cram_map = {}
    for s, reads in per_sample.items():
        reads.sort(key=lambda x: (x[0], x[1]))
        path = os.path.join(sam_dir, f"{s}.sam")
        with open(path, "w") as fh:
            fh.write("@HD\tVN:1.6\tSO:coordinate\n")
            for c, n in CONTIGS.items():
                fh.write(f"@SQ\tSN:{c}\tLN:{n}\n")
            for _, _, line in reads:
                fh.write(line + "\n")
        cram_map[s] = os.path.join(sam_dir, f"{s}.cram")  # produced from the SAM by the runner
    with open(os.path.join(W, "cram_map.tsv"), "w") as fh:
        for s, p in sorted(cram_map.items()):
            fh.write(f"{s}\t{p}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    W = os.path.abspath(args.out)
    os.makedirs(os.path.join(W, "vcfs"), exist_ok=True)

    write_reference(os.path.join(W, "reference.fa"))
    write_source_sams(W)
    for fk, samples in FILES.items():
        write_vcf(os.path.join(W, "vcfs", f"file{fk}.vcf"), samples,
                  [v for v in V if v["file"] == fk])

    # gnomAD sites (only variants with an faf95)
    gnomad_hdr = [
        '##INFO=<ID=AF,Number=1,Type=Float,Description="af">',
        '##INFO=<ID=AF_grpmax,Number=1,Type=Float,Description="grpmax af">',
        '##INFO=<ID=faf95,Number=1,Type=Float,Description="faf95">',
        '##INFO=<ID=nhomalt,Number=1,Type=Integer,Description="nhomalt">',
    ]
    gwant = [(v, f"AF={v['faf95']:.4g};AF_grpmax={v['faf95']:.4g};faf95={v['faf95']:.4g};nhomalt=0")
             for v in V if v["faf95"] is not None]
    write_sites(os.path.join(W, "gnomad.sites.vcf"), gnomad_hdr, gwant)

    # ClinVar sites
    cv_hdr = [
        '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="clinical significance">',
        '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="review status">',
    ]
    cwant = [(v, f"CLNSIG={v['clnsig']};CLNREVSTAT={v['clnrevstat']}")
             for v in V if v["clnsig"]]
    write_sites(os.path.join(W, "clinvar.vcf"), cv_hdr, cwant)

    # annotation lookup for mock_annotate.py
    with open(os.path.join(W, "annot.tsv"), "w") as fh:
        fh.write("chrom\tpos\tref\talt\tgene\tcsq\timpact\trevel\tam\tspliceai\tloftee\n")
        seen = set()
        for v in V:
            key = (v["chrom"], v["pos"])
            if key in seen:
                continue
            seen.add(key)
            fh.write(f"{v['chrom']}\t{v['pos']}\t{refbase(v['pos'])}\t{altbase(v['pos'])}\t"
                     f"{v['gene']}\t{v['csq']}\t{v['impact']}\t{v['revel']}\t{v['am']}\t"
                     f"{v['spliceai']}\t{v['loftee']}\n")

    # Step-6 tables
    with open(os.path.join(W, "mutrate.tsv"), "w") as fh:
        fh.write("gene\tmut_lof\tmut_mis\n")
        for g in ("GENE1", "GENE2", "GENE3", "GENE5", "GENEX"):
            fh.write(f"{g}\t1e-6\t1e-5\n")
    with open(os.path.join(W, "constraint.tsv"), "w") as fh:
        fh.write("gene\toe_lof_upper\tpli\ts_het\n")
        fh.write("GENE1\t0.2\t0.98\t0.15\n")
        fh.write("GENED\t0.25\t0.95\t0.12\n")

    # trios file (#kid dad mom); C is unresolvable (MO_C absent everywhere)
    with open(os.path.join(W, "trios.tsv"), "w") as fh:
        fh.write("#kid\tdad\tmom\nCH_A\tFA_A\tMO_A\nCH_B\tFA_B\tMO_B\nCH_C\tFA_C\tMO_C\n")

    # minimal config (defaults fill in thresholds); concrete ephemeral paths
    with open(os.path.join(W, "config.mock.yaml"), "w") as fh:
        fh.write(f"""project: {{name: mock, genome_build: GRCh38, output_dir: {W}/work}}
runtime: {{image: none, engine: native, tmpdir: {W}/work/tmp, threads: 1}}
reference: {{fasta: {W}/reference.fa}}
resources:
  mutation_rate_table: {W}/mutrate.tsv
  constraint: {{gnomad_v2_constraint: {W}/constraint.tsv}}
  cram_map: {W}/cram_map.tsv
inputs:
  trios_file: {W}/trios.tsv
  vcf_dir: {W}/vcfs
  vcf_list: ""
qc:
  sex_min_sites: 2            # tiny mock has only a handful of chrX sites (prod default is 20)
outputs:
  xlsx: true
  igv: {{enabled: true, padding: 200, genome: hg38}}
""")
    print(f"mock data written to {W}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
