#!/usr/bin/env python3
"""Resolve trios to VCFs and generate PEDs (pipeline preflight).

Input is a simple, order-independent kid/dad/mom file (sample IDs matching the VCFs)
plus a VCF source (directory and/or list). For each trio this:
  1. locates the VCF(s) containing ALL THREE members (matched EXACTLY by sample ID,
     never by column order; additional members in the VCF are fine),
  2. picks the most trio-specific VCF (fewest samples; lexical tie-break),
  3. validates membership and generates a standard PED,
  4. emits a resolved manifest (trio_id, vcf, ped, samples) for the rest of the
     pipeline, plus a resolution audit (what mapped where, and why anything didn't).

Unresolved/ambiguous trios are reported loudly and skipped (never guessed) — mirroring
the group's pedigree.py. The run fails only if NOTHING resolves. See CLAUDE.md.

Usage:
  resolve_trios.py --trios trios.tsv [--vcf-dir DIR] [--vcf-list list.txt] --outdir OUT
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

from cyvcf2 import VCF

from hprv import audit
from hprv.ped import read_trios_file, write_ped


def enumerate_vcfs(vcf_dir, vcf_list):
    vcfs = []
    if vcf_dir:
        for pat in ("*.vcf.gz", "*.vcf", "*.bcf"):
            vcfs += glob.glob(os.path.join(vcf_dir, "**", pat), recursive=True)
    if vcf_list:
        with open(vcf_list) as fh:
            for line in fh:
                p = line.strip()
                if p and not p.startswith("#"):
                    vcfs.append(p)
    # de-dup, keep existing files, stable order
    seen, out = set(), []
    for p in vcfs:
        ap = os.path.abspath(p)
        if ap in seen:
            continue
        seen.add(ap)
        if os.path.exists(p):
            out.append(p)
        else:
            sys.stderr.write(f"WARN: listed VCF not found, skipping: {p}\n")
    return out


def index_samples(vcfs):
    """Return (sample -> set(vcf), vcf -> set(samples)). Header read only."""
    s2v, v2s = {}, {}
    for v in vcfs:
        try:
            samples = set(VCF(v).samples)
        except Exception as e:  # unreadable/corrupt header — report, don't crash
            sys.stderr.write(f"WARN: could not read samples from {v}: {e}\n")
            continue
        v2s[v] = samples
        for s in samples:
            s2v.setdefault(s, set()).add(v)
    return s2v, v2s


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trios", required=True)
    ap.add_argument("--vcf-dir", default="")
    ap.add_argument("--vcf-list", default="")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)

    if not (args.vcf_dir or args.vcf_list):
        sys.stderr.write("ERROR: provide --vcf-dir and/or --vcf-list\n")
        return 2

    trios = read_trios_file(args.trios)
    if not trios:
        sys.stderr.write(f"ERROR: no trios parsed from {args.trios}\n")
        return 2

    vcfs = enumerate_vcfs(args.vcf_dir, args.vcf_list)
    if not vcfs:
        sys.stderr.write("ERROR: no VCFs found from the given source\n")
        return 2
    s2v, v2s = index_samples(vcfs)

    ped_dir = os.path.join(args.outdir, "peds")
    os.makedirs(ped_dir, exist_ok=True)
    manifest = os.path.join(args.outdir, "trios.resolved.tsv")
    res_audit = os.path.join(args.outdir, "trio_resolution.tsv")

    n_res = n_unres = n_ambig = 0
    with open(manifest, "w") as mf, open(res_audit, "w") as af:
        mf.write("trio_id\tvcf\tped\tsamples\n")
        af.write("kid\tdad\tmom\tstatus\tchosen_vcf\tn_candidate_vcfs\t"
                 "missing_members\tn_samples_in_chosen\n")
        for kid, dad, mom in trios:
            members = {"kid": kid, "dad": dad, "mom": mom}
            missing = [f"{role}:{sid}" for role, sid in members.items() if sid not in s2v]
            candidates = (s2v.get(kid, set()) & s2v.get(dad, set()) & s2v.get(mom, set()))
            if not candidates:
                n_unres += 1
                miss = ",".join(missing) if missing else "no VCF has all three together"
                af.write(f"{kid}\t{dad}\t{mom}\tunresolved\t\t0\t{miss}\t\n")
                sys.stderr.write(f"WARN: {kid}: unresolved — {miss}\n")
                continue
            # most trio-specific: fewest samples, then lexical path
            chosen = sorted(candidates, key=lambda v: (len(v2s[v]), v))[0]
            status = "resolved" if len(candidates) == 1 else "resolved_multi"
            if len(candidates) > 1:
                n_ambig += 1
            n_res += 1
            ped = os.path.join(ped_dir, f"{kid}.ped")
            write_ped(ped, kid, dad, mom, kid_sex="0")  # sex unknown; Step 5 infers X-ploidy
            mf.write(f"{kid}\t{chosen}\t{ped}\t{kid},{dad},{mom}\n")
            af.write(f"{kid}\t{dad}\t{mom}\t{status}\t{chosen}\t{len(candidates)}\t"
                     f"\t{len(v2s[chosen])}\n")

    audit.record("resolve", "trios_input", len(trios))
    audit.record("resolve", "vcfs_scanned", len(v2s))
    audit.record("resolve", "samples_indexed", len(s2v))
    audit.record("resolve", "trios_resolved", n_res)
    audit.record("resolve", "trios_unresolved", n_unres)
    audit.record("resolve", "trios_multi_vcf", n_ambig)

    sys.stderr.write(
        f"Resolve complete: {n_res}/{len(trios)} trios resolved "
        f"({n_unres} unresolved, {n_ambig} matched >1 VCF) over {len(v2s)} VCFs.\n"
        f"  manifest: {manifest}\n  resolution audit: {res_audit}\n"
    )
    if n_res == 0:
        sys.stderr.write("ERROR: no trios resolved to a VCF — check IDs and VCF source\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
