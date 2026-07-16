#!/usr/bin/env python3
"""Generate a tiny, self-consistent mock dataset exercising the whole pipeline.

Writes under --out:
  reference.fa                mini GRCh38-like genome (chr1, chr2, chrX-with-nonPAR)
  vcfs/fileA.vcf              trio A (CH_A/FA_A/MO_A) — autosomal modes + filter cases
  vcfs/fileB.vcf              FAMILY VCF: [MO_B, SIB_B, CH_B, FA_B] (extra sibling,
                              shuffled order) — de novo (recurrent gene) + X-linked
  vcfs/fileC.vcf              duo CH_C/FA_C (mom absent) — resolver "unresolved" case
  trios.tsv                   #kid dad mom (A, B resolvable; C unresolvable)
  annot.tsv                   per-site lookup mock_vep.py turns into a VEP CSQ (frequency +
                              CLIN_SIG + CADD included — there is no external sites VCF)
  mutrate.tsv, constraint.tsv gene tables for Step 6
  config.mock.yaml            config pointing at the above (concrete paths, ephemeral)

All variants are SNVs (no indel left-align ambiguity). REF matches the reference at
each position. Not committed data — generated into the git-ignored work dir.
"""
from __future__ import annotations

import argparse
import os

# chrM is present ON PURPOSE: Step 1 must EXCLUDE it (it is out of scope; see
# 01_make_cohort_sites.sh EXCLUDE_CONTIGS). A mock without chrM cannot prove the filter works.
CONTIGS = {"chr1": 20000, "chr2": 20000, "chrX": 2782200, "chrM": 16569}  # chrX > PAR1 end (2,781,479)
BASES = "ACGT"


def refbase(pos):
    return BASES[pos % 4]


def altbase(pos):
    return BASES[(pos + 1) % 4]


def ad(gt, dp, n_alt=1):
    """Biallelic-style AD for a genotype. n_alt > 1 pads the Number=R array to REF + n ALTs.

    A non-ref/non-ref (1/2) genotype must NOT use this — its AD is the whole point of the test
    and is supplied explicitly via `adov` (see the GENECH2 comp-het case).
    """
    pad = ",0" * (n_alt - 1)
    if gt == "0/0":
        return f"{dp},0{pad}"
    if gt == "1/1":
        return f"0,{dp}{pad}"
    h = dp // 2
    return f"{dp - h},{h}{pad}"


# Sample groupings per VCF file (note fileB order is shuffled and has an extra sib).
FILES = {
    "A": ["CH_A", "FA_A", "MO_A"],
    "B": ["MO_B", "SIB_B", "CH_B", "FA_B"],
    "C": ["CH_C", "FA_C"],
}

# Each variant: file, chrom, pos, gene, csq, impact, cadd, af (None = absent from gnomAD),
# af_pop (which gnomAD CSQ population carries `af`), clnsig, filter, hidenovo, gts.
#
# VEP-only contract: there is no gnomAD/ClinVar/dbNSFP/SpliceAI/LOFTEE file to mock, because the
# pipeline no longer reads one. Frequency and CLIN_SIG are emitted INTO the CSQ by mock_vep.py,
# exactly as a real `vep --af_gnomade --af_gnomadg --check_existing` run would.
#
# af_pop defaults to a grpmax-ELIGIBLE group, so `af` drives annotations.frequency(). Set it to a
# bottlenecked group (ami/asj/fin/mid) to model an allele grpmax deliberately ignores.
V = []


def add(**k):
    k.setdefault("cadd", ""); k.setdefault("af", None); k.setdefault("af_pop", "gnomADe_NFE_AF")
    k.setdefault("clnsig", ""); k.setdefault("filter", "PASS"); k.setdefault("hidenovo", "")
    V.append(k)


# --- Trio A (autosomal), CH_A female ---------------------------------------
# 1) de novo, HIGH LoF, absent -> expect mode=denovo
add(file="A", chrom="chr1", pos=5000, gene="GENE1", csq="stop_gained", impact="HIGH",
    hidenovo="CH_A",
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/0", 99, 40)})
# 2) homozygous recessive, MODERATE missense, rare -> mode=hom_recessive
add(file="A", chrom="chr1", pos=8000, gene="GENE2", csq="missense_variant", impact="MODERATE",
    af=5e-4,
    gts={"CH_A": ("1/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/1", 99, 40)})
# 3+4) compound het in GENE3 (var3 maternal, var4 paternal) -> mode=compound_het
add(file="A", chrom="chr2", pos=5000, gene="GENE3", csq="missense_variant", impact="MODERATE",
    af=1e-3,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/1", 99, 40)})
add(file="A", chrom="chr2", pos=6000, gene="GENE3", csq="missense_variant", impact="MODERATE",
    af=1e-3,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# 5) common (BA1) -> dropped at Step 3 (never a candidate)
add(file="A", chrom="chr1", pos=12000, gene="GENE4", csq="missense_variant", impact="MODERATE",
    af=0.2,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# 6) low-GQ de-novo-looking (no hiConfDeNovo) -> passes Step 3, FAILS Step 5 QC
add(file="A", chrom="chr1", pos=15000, gene="GENE5", csq="stop_gained", impact="HIGH", gts={"CH_A": ("0/1", 12, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/0", 99, 40)})
# 7) non-PASS -> dropped at Step 1
add(file="A", chrom="chr1", pos=17000, gene="GENE6", csq="stop_gained", impact="HIGH", filter="VQSRTrancheSNP99.00to99.90+", hidenovo="CH_A",
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/0", 99, 40)})
# 8) ClinVar P/LP but LOW impact (synonymous) -> kept at Step 3 via clinvar_plp override
add(file="A", chrom="chr2", pos=8000, gene="GENE7", csq="synonymous_variant", impact="LOW",
    af=2e-4, clnsig="pathogenic",  # > dominant_max: kept via ClinVar but not a dominant call
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- Trio B (CH_B male), family VCF with extra sibling ----------------------
# 9) de novo in GENE1 (recurrent across trios A+B) -> mode=denovo; drives Step-6 burden
add(file="B", chrom="chr1", pos=5100, gene="GENE1", csq="stop_gained", impact="HIGH",
    hidenovo="CH_B",
    gts={"CH_B": ("0/1", 99, 40), "FA_B": ("0/0", 99, 40), "MO_B": ("0/0", 99, 40),
         "SIB_B": ("0/0", 99, 40)})
# 10) X-linked recessive (male hemizygous), carrier mother -> mode=x_linked_recessive
add(file="B", chrom="chrX", pos=2781600, gene="GENEX", csq="missense_variant", impact="MODERATE",
    af=1e-4,
    gts={"CH_B": ("1/1", 99, 40), "FA_B": ("0/0", 99, 40), "MO_B": ("0/1", 99, 40),
         "SIB_B": ("0/0", 99, 40)})
# 11-13) chrX filler (common) so CH_B is inferred MALE (hemizygous alt -> low het ratio)
for i, p in enumerate((2781700, 2781800, 2781900)):
    add(file="B", chrom="chrX", pos=p, gene=f"XFILL{i}", csq="missense_variant", impact="MODERATE",
        af=0.3,
        gts={"CH_B": ("1/1", 99, 40), "FA_B": ("1/1", 99, 40), "MO_B": ("0/1", 99, 40),
             "SIB_B": ("0/1", 99, 40)})

# --- DOMINANT RECURRENCE: same rare functional inherited het in GENED, in BOTH trios
#     (CH_A inherits from dad, CH_B from mom) -> Step 6 nominates GENED (n_dominant=2) ---
add(file="A", chrom="chr2", pos=10000, gene="GENED", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
add(file="B", chrom="chr2", pos=10000, gene="GENED", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_B": ("0/1", 99, 40), "MO_B": ("0/1", 99, 40), "FA_B": ("0/0", 99, 40),
         "SIB_B": ("0/0", 99, 40)})

# --- CONTAMINATION: a common hom-alt site where FA_B carries reference reads at a hom-alt
#     genotype (verifyBamID-style contamination). Common (af=0.3) so Step 3 drops it and it
#     never becomes a candidate — it only exercises Step 0's CHARR gate. Mendelian-consistent
#     (all 1/1) so it adds no MIE. Expect: CH_B contam_flag=1 (dad), kid/mom clean. ---
add(file="B", chrom="chr2", pos=14000, gene="XCONTAM", csq="missense_variant", impact="MODERATE",
    af=0.3,
    gts={"CH_B": ("1/1", 99, 40), "FA_B": ("1/1", 99, 40), "MO_B": ("1/1", 99, 40),
         "SIB_B": ("1/1", 99, 40)},
    adov={"FA_B": "6,34"})   # 6 ref reads at a hom-alt site -> CHARR 0.15 > 0.02 threshold

# --- X-linked recessive with an AFFECTED FATHER (hom-alt hemizygous) + carrier mother: the son
#     must STILL be called (father transmits Y, not X, to a son) -> tests the father-genotype
#     relaxation. Flag 'father_carries_x_allele' expected. ---
add(file="B", chrom="chrX", pos=2782000, gene="GENEXAF", csq="missense_variant", impact="MODERATE",
    af=1e-4,
    gts={"CH_B": ("1/1", 99, 40), "FA_B": ("1/1", 99, 40), "MO_B": ("0/1", 99, 40),
         "SIB_B": ("0/0", 99, 40)})

# --- autosomal hom-recessive with a HOM-ALT parent (consanguinity-like): FA_A hom-alt, MO_A het,
#     CH_A hom-alt -> hom_recessive via the {HET,HOM_ALT} carrier rule (tests carrier_ok HOM_ALT). ---
add(file="A", chrom="chr1", pos=8500, gene="GENE2H", csq="missense_variant", impact="MODERATE",
    af=5e-4,
    gts={"CH_A": ("1/1", 99, 40), "FA_A": ("1/1", 99, 40), "MO_A": ("0/1", 99, 40)})

# --- DISTINCT-variant dominant recurrence across trios A+B in GENEDD (two DIFFERENT rare hets) ->
#     the stronger 'gene signal'; recurrence_kind=distinct_variant (ranks above same-variant GENED). ---
add(file="A", chrom="chr2", pos=11000, gene="GENEDD", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
add(file="B", chrom="chr2", pos=11100, gene="GENEDD", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_B": ("0/1", 99, 40), "MO_B": ("0/1", 99, 40), "FA_B": ("0/0", 99, 40),
         "SIB_B": ("0/0", 99, 40)})

# --- comp-het CIS rejection: two rare functional hets in GENEC, BOTH inherited from mom
#     (cis). compound_het requires TRANS (mat x pat), so these must NOT be paired; each is a
#     dominant (maternal-origin) call instead. Guards the trans-pairing logic. ---
add(file="A", chrom="chr2", pos=16000, gene="GENEC", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/1", 99, 40)})
add(file="A", chrom="chr2", pos=16100, gene="GENEC", csq="missense_variant", impact="MODERATE",
    af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/0", 99, 40), "MO_A": ("0/1", 99, 40)})

# --- CADD-only keep: a deep-intronic MODIFIER with no impact-based evidence. CADD is the ONLY
#     functional predictor left, so this is the sole path by which any non-coding variant can
#     survive Step 3. If the CADD branch ever breaks, the screen silently goes coding-only and
#     this is the assertion that notices. Inherited het -> also a dominant call in GENEIN. ---
add(file="A", chrom="chr1", pos=18000, gene="GENEIN", csq="intron_variant", impact="MODIFIER",
    cadd="27.5", af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})
# ...and its control: same intronic MODIFIER, CADD BELOW the cutoff -> must be dropped.
add(file="A", chrom="chr1", pos=18500, gene="GENEINLO", csq="intron_variant", impact="MODIFIER",
    cadd="3.0", af=5e-5,
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- FOUNDER-POPULATION allele: frequent ONLY in a bottlenecked group gnomAD's grpmax excludes
#     (mid, ~AN 700). MAX_AF reports 0.002 — 20x over dominant_max — but grpmax-eligible groups
#     report nothing, so annotations.frequency() must return None and the variant must SURVIVE as
#     a dominant candidate. This is the concrete false-negative that using VEP's MAX_AF as the
#     rarity field would cause, and the reason frequency() reads only GRPMAX_POPS. ---
add(file="A", chrom="chr2", pos=17000, gene="GENEFND", csq="missense_variant", impact="MODERATE",
    af=0.002, af_pop="gnomADe_MID_AF",
    gts={"CH_A": ("0/1", 99, 40), "FA_A": ("0/1", 99, 40), "MO_A": ("0/0", 99, 40)})

# --- MULTIALLELIC trans compound het: child is 1/2 (one allele from each parent), which is the
#     textbook presentation of a recessive diagnosis. `bcftools norm -m-` splits this into two
#     records and, WITHOUT --keep-sum AD, discards the other ALT's reads from each leg — leaving
#     ref_ad~0, so allele_balance() reads ~1.0, the het band rejects BOTH legs, and the whole
#     compound_het vanishes with no warning and no audit counter. This case is the regression
#     test for that (Step 4's --keep-sum AD): expect a compound_het pair in GENECH2.
#     AD is explicit because it IS the thing under test:
#       CH_A 1/2 -> 0 ref, 19 for G, 20 for T  (legs must read AB 0.487 / 0.513, not 1.0)
#       FA_A 0/1 -> transmits G ;  MO_A 0/2 -> transmits T   => trans, not cis. ---
add(file="A", chrom="chr1", pos=19000, gene="GENECH2", csq="missense_variant", impact="MODERATE",
    alt2="T", af=5e-5,
    gts={"CH_A": ("1/2", 99, 39), "FA_A": ("0/1", 99, 38), "MO_A": ("0/2", 99, 38)},
    adov={"CH_A": "0,19,20", "FA_A": "18,20,0", "MO_A": "17,0,21"})

# --- chrM: MUST be excluded at Step 1 (out of scope; a dedicated mtDNA pipeline owns it).
#     Modelled on the real failure: m.8860A>G is a near-fixed rCRS haplogroup variant, so the
#     WHOLE TRIO is hom-alt. Left un-excluded it fires hom_recessive in every trio, and with no
#     gnomAD mito AF the rarity gate passes unconditionally -> Step 6 floors q -> p ~ 1e-12 ->
#     it lands in the recurrent exome-wide-significant tier above real nuclear candidates. ---
add(file="A", chrom="chrM", pos=8860, gene="MT-ATP6", csq="missense_variant", impact="MODERATE",
    gts={"CH_A": ("1/1", 99, 400), "FA_A": ("1/1", 99, 400), "MO_A": ("1/1", 99, 400)})
add(file="B", chrom="chrM", pos=8860, gene="MT-ATP6", csq="missense_variant", impact="MODERATE",
    gts={"CH_B": ("1/1", 99, 400), "FA_B": ("1/1", 99, 400), "MO_B": ("1/1", 99, 400),
         "SIB_B": ("1/1", 99, 400)})

# --- Trio C: duo only (mom MO_C absent everywhere) -> resolver unresolved ---
add(file="C", chrom="chr1", pos=9000, gene="GENE8", csq="missense_variant", impact="MODERATE",
    af=1e-3, gts={"CH_C": ("0/1", 99, 40), "FA_C": ("0/1", 99, 40)})


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
            ref = refbase(v["pos"])
            alt = altbase(v["pos"])
            n_alt = 1
            if v.get("alt2"):            # multiallelic: ALT becomes "G,T"
                alt = f"{alt},{v['alt2']}"
                n_alt = 2
            # AC/AF are Number=A: ONE VALUE PER ALT, counted per allele index. A single summed
            # count is not just imprecise, it is malformed for a multiallelic — bcftools rejects
            # it outright ("wrong number of fields in INFO/AF ... expected 2, found 1").
            an = 2 * len(samples)
            acs = [sum(g[0].replace("|", "/").split("/").count(str(i))
                       for g in v["gts"].values())
                   for i in range(1, n_alt + 1)]
            info = ("AC=" + ",".join(str(a) for a in acs) + f";AN={an};AF="
                    + ",".join(f"{a / an:.4g}" for a in acs))
            if v["hidenovo"]:
                info += f";hiConfDeNovo={v['hidenovo']}"
            cells = []
            adov = v.get("adov", {})     # per-sample AD override (contamination; multiallelic AD)
            for s in samples:
                gt, gq, dp = v["gts"][s]
                a = adov.get(s, ad(gt, dp, n_alt))
                cells.append(f"{gt}:{a}:{dp}:{gq}")
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

    # No gnomad.sites.vcf / clinvar.vcf: under the VEP-only contract nothing is transferred from
    # an external sites VCF, so there is nothing to mock. Frequency + CLIN_SIG go into the CSQ
    # (see mock_vep.py), which is where a real `vep --af_gnomade --check_existing` puts them.

    # annotation lookup for mock_vep.py -> becomes the CSQ.
    # Keyed per (chrom, pos, ALT), not per site: Step 1's `norm -m-` splits a multiallelic into
    # one record per ALT, and mock_vep.py matches on the exact allele — so a site with alt2 needs
    # BOTH alleles here or the second leg reaches Step 3 with no CSQ and is dropped as
    # not_functional, quietly destroying the very comp-het the multiallelic case exists to test.
    with open(os.path.join(W, "annot.tsv"), "w") as fh:
        fh.write("chrom\tpos\tref\talt\tgene\tcsq\timpact\tcadd\taf\taf_pop\tclnsig\n")
        seen = set()
        for v in V:
            alts = [altbase(v["pos"])]
            if v.get("alt2"):
                alts.append(v["alt2"])
            for a in alts:
                key = (v["chrom"], v["pos"], a)
                if key in seen:
                    continue
                seen.add(key)
                af = "" if v["af"] is None else f"{v['af']:.6g}"
                fh.write(f"{v['chrom']}\t{v['pos']}\t{refbase(v['pos'])}\t{a}\t"
                         f"{v['gene']}\t{v['csq']}\t{v['impact']}\t{v['cadd']}\t{af}\t"
                         f"{v['af_pop']}\t{v['clnsig']}\n")

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
  # Step 2 ingests this instead of invoking `vep` (mock_vep.py writes it). Everything else in
  # Step 2 — build checks, split-vep, selector, frequency guard — runs for real against it.
  vep: {{annotated_vcf: {W}/cohort.sites.vep.vcf.gz, version: 115}}
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
