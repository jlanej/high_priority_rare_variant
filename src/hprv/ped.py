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

OPTIONAL sex columns make the trios file the CANONICAL source of sex: `kid_sex` (aliases `sex`,
`child_sex`, `proband_sex`; `gender` spellings tolerated) and `dad_sex` / `mom_sex` (`father_sex`,
`mother_sex`, ...). Values are 1/2/0 or male/female/unknown (m/f; blank, `.`, NA = unknown).
A stated proband sex is written to the PED and is what Step 5 judges chrX ploidy under; Step 0's
chrX heterozygosity inference then only CHECKS it (see resolve_child_sex, the one precedence
rule). Absent columns leave the proband unknown (`0`, inference fills in) and the parents at
their role (father 1, mother 2). A parent sex that CONTRADICTS the role is an error — it is the
signature of a transposed father/mother, and the pipeline never guesses which column is wrong.
"""

from __future__ import annotations

import sys
from typing import List, NamedTuple, Optional, Tuple

# Column-name aliases (lower-cased, leading '#' stripped, an `_id` / `_sample_id` suffix or a
# `sample_` prefix removed first, so `father_id`, `Mother`, `proband_sample_id` all resolve).
_KID_ALIASES = ("kid", "child", "proband", "sample")
_DAD_ALIASES = ("dad", "father", "paternal")
_MOM_ALIASES = ("mom", "mother", "maternal")
_ROLES = (("kid", _KID_ALIASES), ("dad", _DAD_ALIASES), ("mom", _MOM_ALIASES))

# The OPTIONAL sex columns. A bare `sex`/`gender` is the proband's; the parents' are role-prefixed
# (`dad_sex`, `father_sex`, `paternal_gender`, ...). `sample_sex` canonicalises to `sex` (kid).
_SEX_WORDS = ("sex", "gender")
_KID_SEX_ALIASES = _SEX_WORDS + tuple(f"{a}_{w}" for a in _KID_ALIASES for w in _SEX_WORDS)
_DAD_SEX_ALIASES = tuple(f"{a}_{w}" for a in _DAD_ALIASES for w in _SEX_WORDS)
_MOM_SEX_ALIASES = tuple(f"{a}_{w}" for a in _MOM_ALIASES for w in _SEX_WORDS)
_SEX_ROLES = (("kid_sex", _KID_SEX_ALIASES), ("dad_sex", _DAD_SEX_ALIASES),
              ("mom_sex", _MOM_SEX_ALIASES))
# The sex a role ASSERTS: the PED father is male and the PED mother is female by definition —
# every hemizygous model in Step 5 is keyed on the mother's genotype because of it.
ROLE_SEX = {"dad_sex": "1", "mom_sex": "2"}

# Accepted sex spellings -> PED code. Anything else is an ERROR (a typo must not read as unknown
# and silently hand the trio to the chrX inference).
_SEX_VALUES = {"1": "1", "2": "2", "0": "0", "male": "1", "female": "2", "unknown": "0",
               "m": "1", "f": "2", "u": "0", "": "0", ".": "0", "na": "0", "n/a": "0",
               "nan": "0", "none": "0", "unk": "0"}
SEX_VOCABULARY = "1/2/0 or male/female/unknown (m/f; blank, '.', NA = unknown)"


class TrioRow(NamedTuple):
    """One trios-file row: the three sample IDs and each member's PED sex code ('1'/'2'/'0').

    `dad_sex`/`mom_sex` are always their role's code on a successfully parsed file (an absent or
    unknown value defaults to the role; a contradiction raises) — carried explicitly so write_ped
    receives what the file said rather than re-deriving it.
    """
    kid: str
    dad: str
    mom: str
    kid_sex: str = "0"
    dad_sex: str = "1"
    mom_sex: str = "2"


def parse_sex(value) -> str:
    """Normalise a sex value to the PED code '1' (male) / '2' (female) / '0' (unknown).

    Raises ValueError on anything outside SEX_VOCABULARY, so a mis-typed cell fails the file
    rather than reading as unknown.
    """
    key = str(value if value is not None else "").strip().lower()
    if key not in _SEX_VALUES:
        raise ValueError(f"unrecognised sex value {value!r}; accepted: {SEX_VOCABULARY}")
    return _SEX_VALUES[key]


def resolve_child_sex(ped_sex, inferred_sex) -> Tuple[str, str, bool]:
    """THE precedence rule for a proband's sex: (sex, source, discordant).

    * A PED (trios-file) sex of 1/2 is CANONICAL -> (ped_sex, 'ped', discordant), where
      `discordant` is True when Step 0's chrX inference is positive (1/2) and disagrees. The PED
      value is kept regardless — the disagreement is REPORTED (Step 0 `sex_match=0`, Step 5's
      `sex_discordant_inference` flag), never resolved by silently replacing a pedigree with a
      heuristic.
    * Otherwise a positive inference fills in -> (inferred, 'inferred', False).
    * Otherwise ('0', 'none', False): Step 5 skips the sex chromosomes for the trio.

    Every consumer — Step 0's `sex_match`/`sex_source`, Step 5's chrX ploidy, Step 6's male
    proband count — reads this function, so the three can never disagree on who is male.
    """
    p = str(ped_sex if ped_sex is not None else "").strip()
    i = str(inferred_sex if inferred_sex is not None else "").strip()
    if p in ("1", "2"):
        return p, "ped", (i in ("1", "2") and i != p)
    if i in ("1", "2"):
        return i, "inferred", False
    return "0", "none", False


def _canon(col: str) -> str:
    c = col.strip().lstrip("#").strip().lower()
    for suf in ("_sample_id", "_sample", "_id"):       # suffix first: `sample_id` -> `sample`
        if c.endswith(suf) and len(c) > len(suf):
            c = c[: -len(suf)]
            break
    if c.startswith("sample_") and len(c) > len("sample_"):   # then `sample_kid` -> `kid`
        c = c[len("sample_"):]
    return c


def _resolve_header(fields) -> Tuple[dict, dict]:
    """Map role -> column index for the ID roles named in a header line (may be partial/empty),
    and sex-role -> column index for the optional sex columns."""
    roles, sex_cols = {}, {}
    for i, col in enumerate(fields):
        c = _canon(col)
        for role, aliases in _ROLES:
            if c in aliases:
                if role in roles:
                    raise ValueError(f"two columns name the {role}: {fields[roles[role]]!r} and {col!r}")
                roles[role] = i
        for role, aliases in _SEX_ROLES:
            if c in aliases:
                if role in sex_cols:
                    raise ValueError(f"two columns name the {role}: {fields[sex_cols[role]]!r} and {col!r}")
                sex_cols[role] = i
    return roles, sex_cols


def read_trios_file(path: str) -> List[TrioRow]:
    """Read a kid/dad/mom trio file. Returns a list of TrioRow (kid, dad, mom, kid_sex, dad_sex,
    mom_sex); `row[:3]` is the old (kid, dad, mom) tuple.

    Columns are located by header NAME, never by position, when a header is present. A header
    that names only some of the three roles is an ERROR rather than a positional guess: the old
    fallback silently filled the unresolved roles with columns 1 and 2, so a file whose header
    spelled the parents `mother_id`/`father_id` (unrecognised at the time) in that order had
    every trio's parents transposed — inverting every parent-of-origin call and fabricating
    X-linked calls for the whole trio — and Step 0's Mendelian-error gate is mathematically
    symmetric under a father/mother swap, so nothing downstream could see it.

    The OPTIONAL sex columns (`kid_sex`/`sex`, `dad_sex`, `mom_sex`; see the module docstring)
    are read by name too. An absent column or an unknown value leaves the proband `0` and puts
    each parent at its role's sex; a parent value that contradicts the role (a female father, a
    male mother) is a ValueError naming the line — it is exactly the signature of a transposed
    pair, and which column is wrong cannot be guessed. An unrecognised sex spelling is an error
    too (parse_sex), never a silent unknown.

    A file with NO header at all (its first non-comment line is a trio) is read positionally as
    kid/dad/mom, and that first line is a trio, not a header. Leading `#` comment lines that name
    no role are skipped. Rows missing any of the three IDs are skipped. Tab-separated.
    """
    out: List[TrioRow] = []
    with open(path) as fh:
        lines = fh.read().splitlines()
    kid = dad = mom = None
    sex_cols: dict = {}
    start = 0
    for i, line in enumerate(lines):
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        try:
            roles, sex_cols = _resolve_header(fields)
        except ValueError as e:
            raise ValueError(f"trios file {path}: {e}") from None
        if len(roles) == 3:
            kid, dad, mom = roles["kid"], roles["dad"], roles["mom"]
            start = i + 1
            break
        if roles or sex_cols:
            missing = [r for r, _ in _ROLES if r not in roles]
            raise ValueError(
                f"trios file {path}: header {fields!r} names {sorted(roles) + sorted(sex_cols)} "
                f"but not {missing}. Recognised column names (case-insensitive; an `_id` suffix is "
                f"fine): kid/child/proband, dad/father/paternal, mom/mother/maternal, and the "
                f"optional kid_sex (sex), dad_sex, mom_sex. Column ORDER is never assumed when a "
                f"header is present — a transposed mother/father inverts every parent-of-origin "
                f"call for the trio and no downstream check can detect it.")
        if line.lstrip().startswith("#"):
            continue                      # a comment line naming no role; keep looking
        # headerless: the first data line IS a trio; read positionally (no sex columns)
        kid, dad, mom = 0, 1, 2
        sex_cols = {}
        start = i
        sys.stderr.write(f"WARN: trios file {path} has no header line; reading columns positionally "
                         "as kid, dad, mom (no sex columns; proband sex left unknown)\n")
        break
    if kid is None:
        return out
    for n, line in enumerate(lines[start:], start=start + 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) <= max(kid, dad, mom):
            continue
        k, d, m = f[kid].strip(), f[dad].strip(), f[mom].strip()
        if not (k and d and m):
            continue
        sexes = {}
        for role, _aliases in _SEX_ROLES:
            col = sex_cols.get(role)
            raw = f[col] if (col is not None and col < len(f)) else ""
            try:
                sexes[role] = parse_sex(raw)
            except ValueError as e:
                raise ValueError(f"trios file {path} line {n} ({role} of proband {k!r}): {e}") from None
        for role, want in ROLE_SEX.items():
            stated = sexes[role]
            if stated == "0":
                sexes[role] = want                  # absent/unknown -> the role's sex
            elif stated != want:
                who = "father" if role == "dad_sex" else "mother"
                raise ValueError(
                    f"trios file {path} line {n}: {role}={stated} contradicts the {who} role of "
                    f"{(d if role == 'dad_sex' else m)!r} (a PED {who} is sex {want}). This is the "
                    f"signature of a transposed father/mother — or a wrong sex cell — and the "
                    f"pipeline will not guess which. Fix the roles or the sex column.")
        out.append(TrioRow(k, d, m, sexes["kid_sex"], sexes["dad_sex"], sexes["mom_sex"]))
    return out


def write_ped(path: str, kid: str, dad: str, mom: str, kid_sex: str = "0",
              dad_sex: str = "1", mom_sex: str = "2") -> None:
    """Write a standard PED for one trio (kid affected; parents unaffected).

    kid_sex: '1' male, '2' female, '0' unknown — the trios file's stated sex when it has one
    (CANONICAL for Step 5's chrX ploidy; Step 0 only checks it against the chrX inference), else
    '0' and Step 0's inference fills in downstream (resolve_child_sex). dad_sex/mom_sex default
    to the role (1/2); read_trios_file refuses a value that contradicts it.
    """
    fam = f"FAM_{kid}"
    with open(path, "w") as fh:
        fh.write(f"{fam}\t{kid}\t{dad}\t{mom}\t{kid_sex}\t2\n")
        fh.write(f"{fam}\t{dad}\t0\t0\t{dad_sex}\t1\n")
        fh.write(f"{fam}\t{mom}\t0\t0\t{mom_sex}\t1\n")


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
