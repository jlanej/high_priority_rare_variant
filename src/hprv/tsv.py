"""The ONE place hprv builds a TSV/CSV reader — so Python's 128 KiB field cap never applies.

The stdlib `csv` parser refuses any field longer than `csv.field_size_limit()`, which defaults to
131,072 characters (`_csv.Error: field larger than field limit (131072)`). Step 5's dropless
`info_<ID>` block carries the raw multi-transcript VEP `CSQ` verbatim as `info_CSQ`: one
consequence block per overlapping transcript (Step 2 runs VEP with `--flag_pick`, which keeps every
block), each dozens of fields wide once `--numbers --total_length`, NMD, SpliceAI, SpliceVault,
REVEL and AlphaMissense ride in it. A locus with enough overlapping transcripts writes a cell past the cap,
and without this module every step that re-reads the calls table — 5b (where it surfaced), 6, 7,
8 and 9 — raises on the first such row.

So the limit is raised to the platform maximum here, and every reader in `src/hprv/` and
`pipeline/` is built through `reader()` / `DictReader()` below — never `csv.reader` /
`csv.DictReader` directly. `tests/test_pure.py` enforces that by scanning the source, and
round-trips a 1 MiB cell through each entry point. (Writers are unaffected: the cap is the
parser's.)

Two details are load-bearing:

* **The overflow-safe loop.** `csv.field_size_limit` takes a C `long`, which is 32-bit on Windows
  (LLP64) even where `sys.maxsize` is 2**63-1, so `sys.maxsize` itself raises OverflowError there.
  Halving until it fits lands on the largest value the platform accepts. It runs once, at import.
* **Every constructor re-applies it.** The limit is process-global state. Raised only at import,
  an entry point would be correct merely because some OTHER module had imported this one first,
  and a test running every entry point in one process could never catch the one that bypassed it.
  Re-applying the resolved value is one C call; it makes each call site self-sufficient, and lets
  the test restore the stdlib default before each entry point and prove it.

The cap was not an integrity check worth keeping: it fires on a legitimate long cell exactly as
on a runaway quoted field, and misses a runaway field shorter than 128 KiB.
"""

from __future__ import annotations

import csv
import sys


def _platform_max_field_size() -> int:
    """Set and return the largest field-size limit this platform's `csv` accepts."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return limit
        except OverflowError:
            limit //= 2


FIELD_SIZE_LIMIT = _platform_max_field_size()


def raise_field_size_limit() -> int:
    """(Re-)apply FIELD_SIZE_LIMIT to the process-global `csv` parser; returns it."""
    csv.field_size_limit(FIELD_SIZE_LIMIT)
    return FIELD_SIZE_LIMIT


def reader(f, delimiter="\t", **kw):
    """`csv.reader` over `f` with the field-size limit raised. Tab-delimited unless told otherwise."""
    raise_field_size_limit()
    return csv.reader(f, delimiter=delimiter, **kw)


class DictReader(csv.DictReader):
    """`csv.DictReader` over `f` with the field-size limit raised. Tab-delimited unless told otherwise."""

    def __init__(self, f, fieldnames=None, delimiter="\t", **kw):
        raise_field_size_limit()
        super().__init__(f, fieldnames=fieldnames, delimiter=delimiter, **kw)
