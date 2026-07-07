#!/usr/bin/env python3
"""Pipeline Step 7: consolidated .xlsx supplemental-table summary.

Assembles the run's outputs (gene consolidation, candidate calls, trio resolution,
QC, audit) into one documented workbook. See src/hprv/report.py.

Usage:
  07_report_xlsx.py --work WORKDIR --config config.yaml --out summary.xlsx [--label RUN]
"""
from __future__ import annotations

import argparse
import sys

from hprv.config import load_config
from hprv.report import build


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    out = build(args.work, args.out, cfg, run_label=args.label)
    sys.stderr.write(f"Step 7 complete: xlsx summary -> {out}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
