"""Minimal PED/FAM parsing for trios.

PED columns (whitespace-separated): family, individual, paternal, maternal, sex, [phenotype].
sex: 1=male, 2=female. We take the first individual that has BOTH parents set as the
proband/child of the trio.
"""

from __future__ import annotations

from typing import Optional


def parse_ped(path: Optional[str]) -> Optional[dict]:
    """Return {child, father, mother, sex} sample IDs for the trio, or None."""
    if not path:
        return None
    try:
        rows = []
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                rows.append(line.split())
    except OSError:
        return None
    for f in rows:
        if len(f) >= 5 and f[2] not in ("0", "") and f[3] not in ("0", ""):
            return {"child": f[1], "father": f[2], "mother": f[3], "sex": f[4]}
    return None
