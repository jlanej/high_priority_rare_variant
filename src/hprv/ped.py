"""Minimal PED/FAM parsing and trio-pedigree reading.

PED columns (whitespace-separated): family, individual, paternal, maternal, sex, [phenotype].
sex: 1=male, 2=female. We take the first individual that has BOTH parents set as the
proband/child of the trio.

The user-facing trio input is a simpler tab-separated file with a header naming the kid,
dad, and mom columns (any order), whose values are sample IDs matching the VCFs — e.g.

    #kid    dad     mom
    CH1     FA1     MO1

Column names are matched by alias (kid/child/proband, dad/father/paternal,
mom/mother/maternal, an `_id` suffix tolerated), mirroring the group's proven pedigree.py so
existing files work — and STRICTLY: a header that names some roles but not all is an error,
never a positional guess (see read_trios_file).
"""

from __future__ import annotations

import sys
from typing import List, Optional, Tuple

# Column-name aliases (lower-cased, leading '#' stripped, an `_id` / `_sample_id` suffix or a
# `sample_` prefix removed first, so `father_id`, `Mother`, `proband_sample_id` all resolve).
_KID_ALIASES = ("kid", "child", "proband", "sample")
_DAD_ALIASES = ("dad", "father", "paternal")
_MOM_ALIASES = ("mom", "mother", "maternal")
_ROLES = (("kid", _KID_ALIASES), ("dad", _DAD_ALIASES), ("mom", _MOM_ALIASES))


def _canon(col: str) -> str:
    c = col.strip().lstrip("#").strip().lower()
    for suf in ("_sample_id", "_sample", "_id"):       # suffix first: `sample_id` -> `sample`
        if c.endswith(suf) and len(c) > len(suf):
            c = c[: -len(suf)]
            break
    if c.startswith("sample_") and len(c) > len("sample_"):   # then `sample_kid` -> `kid`
        c = c[len("sample_"):]
    return c


def _resolve_header(fields) -> dict:
    """Map role -> column index for the roles named in a header line (may be partial/empty)."""
    roles = {}
    for i, col in enumerate(fields):
        c = _canon(col)
        for role, aliases in _ROLES:
            if c in aliases:
                if role in roles:
                    raise ValueError(f"two columns name the {role}: {fields[roles[role]]!r} and {col!r}")
                roles[role] = i
    return roles


def read_trios_file(path: str) -> List[Tuple[str, str, str]]:
    """Read a kid/dad/mom trio file. Returns a list of (kid, dad, mom) tuples.

    Columns are located by header NAME, never by position, when a header is present. A header
    that names only some of the three roles is an ERROR rather than a positional guess: the old
    fallback silently filled the unresolved roles with columns 1 and 2, so a file whose header
    spelled the parents `mother_id`/`father_id` (unrecognised at the time) in that order had
    every trio's parents transposed — inverting every parent-of-origin call and fabricating
    X-linked calls for the whole trio — and Step 0's Mendelian-error gate is mathematically
    symmetric under a father/mother swap, so nothing downstream could see it.

    A file with NO header at all (its first non-comment line is a trio) is read positionally as
    kid/dad/mom, and that first line is a trio, not a header. Leading `#` comment lines that name
    no role are skipped. Rows missing any of the three IDs are skipped. Tab-separated.
    """
    out: List[Tuple[str, str, str]] = []
    with open(path) as fh:
        lines = fh.read().splitlines()
    kid = dad = mom = None
    start = 0
    for i, line in enumerate(lines):
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        try:
            roles = _resolve_header(fields)
        except ValueError as e:
            raise ValueError(f"trios file {path}: {e}") from None
        if len(roles) == 3:
            kid, dad, mom = roles["kid"], roles["dad"], roles["mom"]
            start = i + 1
            break
        if roles:
            missing = [r for r, _ in _ROLES if r not in roles]
            raise ValueError(
                f"trios file {path}: header {fields!r} names {sorted(roles)} but not "
                f"{missing}. Recognised column names (case-insensitive; an `_id` suffix is fine): "
                f"kid/child/proband, dad/father/paternal, mom/mother/maternal. Column ORDER is never "
                f"assumed when a header is present — a transposed mother/father inverts every "
                f"parent-of-origin call for the trio and no downstream check can detect it.")
        if line.lstrip().startswith("#"):
            continue                      # a comment line naming no role; keep looking
        # headerless: the first data line IS a trio; read positionally
        kid, dad, mom = 0, 1, 2
        start = i
        sys.stderr.write(f"WARN: trios file {path} has no header line; reading columns positionally "
                         "as kid, dad, mom\n")
        break
    if kid is None:
        return out
    for line in lines[start:]:
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
