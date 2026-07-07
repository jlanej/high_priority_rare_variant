"""Minimal PED/FAM parsing and trio-pedigree reading.

PED columns (whitespace-separated): family, individual, paternal, maternal, sex, [phenotype].
sex: 1=male, 2=female. We take the first individual that has BOTH parents set as the
proband/child of the trio.

The user-facing trio input is a simpler tab-separated file with a header naming the kid,
dad, and mom columns (any order), whose values are sample IDs matching the VCFs — e.g.

    #kid    dad     mom
    CH1     FA1     MO1

Column names are matched by alias (kid/child/proband, dad/father/paternal,
mom/mother/maternal), mirroring the group's proven pedigree.py so existing files work.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

# Column-name aliases (lower-cased, leading '#' stripped).
_KID_ALIASES = ("kid", "child", "proband", "sample_id")
_DAD_ALIASES = ("dad", "father", "paternal")
_MOM_ALIASES = ("mom", "mother", "maternal")


def read_trios_file(path: str) -> List[Tuple[str, str, str]]:
    """Read a kid/dad/mom trio file. Returns a list of (kid, dad, mom) tuples.

    Columns are located by header name (alias-matched), falling back to the first
    three columns (kid, dad, mom) if the standard names are absent. Rows missing any
    of the three IDs are skipped. Tab-separated; '#' on the header is tolerated.
    """
    out: List[Tuple[str, str, str]] = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        kid = dad = mom = None
        for i, col in enumerate(header):
            c = col.lower().strip("#").strip()
            if c in _KID_ALIASES:
                kid = i
            elif c in _DAD_ALIASES:
                dad = i
            elif c in _MOM_ALIASES:
                mom = i
        if kid is None:
            kid = 0
        if dad is None:
            dad = 1
        if mom is None:
            mom = 2
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(kid, dad, mom):
                continue
            k, d, m = f[kid].strip(), f[dad].strip(), f[mom].strip()
            if k and d and m:
                out.append((k, d, m))
    return out


def write_ped(path: str, kid: str, dad: str, mom: str, kid_sex: str = "0") -> None:
    """Write a standard PED for one trio (father sex=1, mother sex=2, kid affected).

    kid_sex: '1' male, '2' female, '0' unknown (Step 5 infers X-ploidy when unknown).
    """
    fam = f"FAM_{kid}"
    with open(path, "w") as fh:
        fh.write(f"{fam}\t{kid}\t{dad}\t{mom}\t{kid_sex}\t2\n")
        fh.write(f"{fam}\t{dad}\t0\t0\t1\t1\n")
        fh.write(f"{fam}\t{mom}\t0\t0\t2\t1\n")


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
