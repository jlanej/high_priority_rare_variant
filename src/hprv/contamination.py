"""Contamination QC: verifyBamID FREEMIX ingestion + a VCF-only CHARR-style estimate.

Two paths, so a contamination estimate is always available for the Step-0 trio gate:
  * PRIMARY — ingest verifyBamID `.selfSM` FREEMIX if a directory of them is configured
    (mirrors the group's DNM freemix QC: read FREEMIX, gate trios at a threshold).
  * FALLBACK — CHARR-style estimate from the per-trio VCF itself (no BAM needed): the
    reference-allele read fraction at high-quality homozygous-ALT sites. An uncontaminated
    sample sits near 0 (sequencing error + reference bias); cross-sample contamination
    introduces reference reads at hom-alt sites and raises it.

CHARR: Lu et al., "CHARR ... contamination from homozygous alternate reference reads,"
Am J Hum Genet 2023 (DOI 10.1016/j.ajhg.2023.10.011).
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
    """VCF-only contamination proxy = reference-read fraction at hom-alt sites, or None."""
    return (ref_sum / dp_sum) if dp_sum and dp_sum > 0 else None
