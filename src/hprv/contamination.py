"""Contamination QC: verifyBamID FREEMIX ingestion + a VCF-only raw ref-read-fraction proxy.

Two paths, so a contamination estimate is always available for the Step-0 trio gate:
  * PRIMARY — ingest verifyBamID `.selfSM` FREEMIX if a directory of them is configured
    (mirrors the group's DNM freemix QC: read FREEMIX, gate trios at a threshold). This is the
    recommended production path.
  * FALLBACK — a raw (uncorrected) reference-read fraction from the per-trio VCF itself (no BAM
    needed): Σref_AD / Σ(ref_AD+alt_AD) over high-quality homozygous-ALT sites. An uncontaminated
    sample sits near 0 (sequencing error + reference bias); cross-sample contamination introduces
    reference reads at hom-alt sites and raises it.

IMPORTANT — this fallback is a CHARR-*like* proxy, NOT the calibrated CHARR estimator. Published
CHARR (Lu et al., Am J Hum Genet 2023, DOI 10.1016/j.ajhg.2023.10.011) is the per-GENOTYPE mean of
ref_AD/DP, divided by the mean reference-allele availability ≈ mean(1−AF) at those sites, and
baseline-corrected for reference-mapping/error bias. This proxy does NONE of that (no AF
normalization, no baseline subtraction, and it pools by depth rather than per-genotype). Because
hom-alt sites are HWE-weighted toward high-AF variants, mean(1−AF) ≈ 0.2–0.35, so the raw value is
only ~1/3 of the true contamination fraction: it flags only GROSS contamination (empirically ≳5–8%),
NOT the 1–3% band the calibrated CHARR literature targets. The `qc.charr_threshold` default (0.02) is
therefore a heuristic on THIS uncorrected scale, not the CHARR corrected-scale cutoff it resembles.
Calibrating a real correction (needs per-site gnomAD AF, unavailable at Step 0 which runs
pre-annotation) is a roadmap TODO; until then, prefer the FREEMIX path for anything near 1–3%.
"""

from __future__ import annotations

import glob
import os


def read_selfsm(selfsm_dir):
    """Return {sample_id: freemix} parsed from verifyBamID ``*.selfSM`` files."""
    out = {}
    if not selfsm_dir or not os.path.isdir(selfsm_dir):
        return out
    for fp in glob.glob(os.path.join(selfsm_dir, "*.selfSM")):
        try:
            with open(fp) as fh:
                header = None
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if header is None:
                        header = {c.strip("#"): i for i, c in enumerate(parts)}
                        continue
                    si = header.get("SEQ_ID", 0)
                    fi = header.get("FREEMIX", 6)
                    if len(parts) > max(si, fi):
                        try:
                            out[parts[si].strip()] = float(parts[fi])
                        except ValueError:
                            pass
        except OSError:
            continue
    return out


def charr(ref_sum, dp_sum):
    """Raw (uncorrected) reference-read fraction at hom-alt sites = Σref_AD / Σ(ref+alt), or None.

    A CHARR-*like* proxy, not the calibrated Lu-2023 statistic (no /mean(1−AF), no baseline
    subtraction) — see the module docstring. Reads only ~1/3 of the true contamination fraction, so
    it flags gross (≳5–8%) contamination only; the 0.02 default gate is a heuristic on this scale."""
    return (ref_sum / dp_sum) if dp_sum and dp_sum > 0 else None
