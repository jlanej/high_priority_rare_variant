"""Load and resolve the pipeline YAML config.

Responsibilities:
  * parse the YAML,
  * expand ``${ENV_VAR}`` placeholders from the environment (so no real paths are
    ever committed — they are supplied at runtime),
  * expose values by dotted key (``filters.rarity.dominant_max``),
  * emit a curated set of shell ``export`` lines for the bash step scripts.

Usable both as a library (``from hprv.config import load_config``) and as a CLI
(``python -m hprv.config sh --config config.yaml``).
"""

from __future__ import annotations

import argparse
import os
import shlex
import sys


def _expand(value):
    """Recursively expand ${ENV} in strings; leave non-strings alone."""
    if isinstance(value, str):
        return os.path.expandvars(value)
    if isinstance(value, dict):
        return {k: _expand(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_expand(v) for v in value]
    return value


def load_config(path: str) -> dict:
    try:
        import yaml
    except ImportError:  # pragma: no cover
        sys.stderr.write("ERROR: pyyaml is required (it is in the container image).\n")
        raise
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise ValueError(f"config {path!r} did not parse to a mapping")
    return _expand(cfg)


def get(cfg: dict, dotted: str, default=None):
    """Fetch ``a.b.c`` from a nested dict, returning ``default`` if absent."""
    cur = cfg
    for key in dotted.split("."):
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


_TRUE, _FALSE = {"true", "yes", "on", "1"}, {"false", "no", "off", "0"}


def as_bool(value, key="") -> bool:
    """A boolean knob, strictly. `bool(value)` read a quoted or `${ENV}`-templated "false" as
    True — the gates a user had just switched OFF stayed on, silently."""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    text = str(value).strip().lower()
    if text in _TRUE:
        return True
    if text in _FALSE:
        return False
    raise ValueError(f"{key or 'config'}: expected true/false, got {value!r}")


def get_bool(cfg: dict, dotted: str, default: bool) -> bool:
    v = get(cfg, dotted, None)
    return default if v is None else as_bool(v, dotted)


IMPACTS = ("HIGH", "MODERATE", "LOW", "MODIFIER")


def keep_impacts(cfg: dict) -> set:
    """`filters.functional.keep_impacts` as an upper-cased set of VEP IMPACT values.

    `set(get(cfg, key, [...]))` on a YAML SCALAR (`keep_impacts: HIGH`) yielded {'H','I','G'}, on a
    lower-cased list yielded nothing VEP ever emits, and on `[]` an empty set — each of which
    silently disabled the ENTIRE impact rung and filed every stop_gained under `not_functional`.
    A scalar or a comma/space-separated string is accepted and normalised; an unknown value or
    an empty set is an error, because a screen with no impact rung is not a screen.
    """
    raw = get(cfg, "filters.functional.keep_impacts", ["HIGH", "MODERATE"])
    if isinstance(raw, str):
        items = [x for x in raw.replace(",", " ").split() if x]
    elif isinstance(raw, (list, tuple, set)):
        items = [str(x) for x in raw]
    else:
        raise ValueError(f"filters.functional.keep_impacts: expected a list of {list(IMPACTS)}, got {raw!r}")
    out = {x.strip().upper() for x in items if x.strip()}
    bad = sorted(out - set(IMPACTS))
    if bad:
        raise ValueError(f"filters.functional.keep_impacts: unknown IMPACT value(s) {bad}; VEP emits only {list(IMPACTS)}")
    if not out:
        raise ValueError("filters.functional.keep_impacts is empty — that disables the impact rung entirely; list at least HIGH")
    return out


_BOOL_KNOBS = ("filters.genotype_qc.require_pass", "filters.denovo.use_hiconf_tag",
               "filters.denovo.crosscheck_prerefinement_pl", "inheritance.emit_denovo",
               "inheritance.emit_dominant", "burden.rank_by_mutational_target",
               "prioritization.composite.gene_list_prior.enabled")


def validate_filters(cfg: dict):
    """Every problem with the screen's knobs, as a list of messages (empty = sound).

    These are the settings that can turn a rung OFF without an error: a scalar keep_impacts, an
    inverted rarity ladder (dominant_max above recessive_max is silently inert), an inverted or
    percent-scale allele-balance band (`het_ab_min: 25` reduces Step 5 to zero calls with exit
    0), a quoted boolean. Steps 3 and 5 call this at start-up and halt on anything returned.
    """
    problems = []
    try:
        keep_impacts(cfg)
    except ValueError as e:
        problems.append(str(e))

    def num(key, default):
        v = get(cfg, key, default)
        try:
            return float(v)
        except (TypeError, ValueError):
            problems.append(f"{key}: expected a number, got {v!r}")
            return None

    dom = num("filters.rarity.dominant_max", 1e-4)
    rec = num("filters.rarity.recessive_max", 1e-2)
    strict = num("filters.rarity.recessive_strict", 1e-3)
    ba1 = num("filters.rarity.benign_ba1", 0.05)
    if None not in (dom, rec, strict, ba1):
        for k, v in (("dominant_max", dom), ("recessive_max", rec), ("recessive_strict", strict), ("benign_ba1", ba1)):
            if not 0 < v <= 1:
                problems.append(f"filters.rarity.{k}: {v} is not an allele frequency in (0, 1]")
        if dom > rec:
            problems.append(f"filters.rarity.dominant_max ({dom}) exceeds recessive_max ({rec}): hets are pooled at recessive_max, so the dominant gate could never be the binding one")
        if rec > ba1:
            problems.append(f"filters.rarity.recessive_max ({rec}) exceeds benign_ba1 ({ba1}): the BA1 drop would fire before the recessive band is reachable")
        if strict > rec:
            problems.append(f"filters.rarity.recessive_strict ({strict}) exceeds recessive_max ({rec}): high_conf_rarity would tag every recessive call")

    g = "filters.genotype_qc."
    lo, hi = num(g + "het_ab_min", 0.25), num(g + "het_ab_max", 0.75)
    ha, hr = num(g + "homalt_ab_min", 0.90), num(g + "homref_ab_max", 0.10)
    if None not in (lo, hi, ha, hr):
        for k, v in (("het_ab_min", lo), ("het_ab_max", hi), ("homalt_ab_min", ha), ("homref_ab_max", hr)):
            if not 0 <= v <= 1:
                problems.append(f"{g}{k}: {v} is not a fraction in [0, 1] (allele balance is alt/(ref+alt), never a percentage)")
        if lo >= hi:
            problems.append(f"{g}het_ab_min ({lo}) >= het_ab_max ({hi}): no het could pass, and every dominant and compound-het call would vanish with exit 0")
        if ha <= 0:
            problems.append(f"{g}homalt_ab_min ({ha}) must be positive")
    for k, floor in ((g + "min_gq", 0), (g + "min_dp", 1), (g + "denovo_min_dp", 1),
                     ("filters.denovo.parent_min_dp", 1), ("filters.denovo.parent_max_alt_ad", 0)):
        v = get(cfg, k, None)
        if v is not None:
            try:
                if int(v) < floor:
                    problems.append(f"{k}: {v} is below {floor}")
            except (TypeError, ValueError):
                problems.append(f"{k}: expected an integer, got {v!r}")
    try:
        if int(get(cfg, g + "denovo_min_dp", 20)) < int(get(cfg, g + "min_dp", 10)):
            problems.append(f"{g}denovo_min_dp is below {g}min_dp — the de novo floor is meant to be the stricter one")
    except (TypeError, ValueError):
        pass
    for k in _BOOL_KNOBS:
        v = get(cfg, k, None)
        if v is not None:
            try:
                as_bool(v, k)
            except ValueError as e:
                problems.append(str(e) + " — a quoted or ${ENV}-templated string is not a YAML boolean")
    return problems


# Curated map: shell variable -> dotted config key. Only what the bash steps need.
SH_MAP = {
    "HPRV_OUTPUT_DIR": "project.output_dir",
    "HPRV_GENOME_BUILD": "project.genome_build",
    "HPRV_IMAGE": "runtime.image",
    "HPRV_ENGINE": "runtime.engine",
    "HPRV_TMPDIR": "runtime.tmpdir",
    "HPRV_THREADS": "runtime.threads",
    "HPRV_REF_FASTA": "reference.fasta",
    "HPRV_TRIOS_FILE": "inputs.trios_file",
    "HPRV_VCF_DIR": "inputs.vcf_dir",
    "HPRV_VCF_LIST": "inputs.vcf_list",
    # VEP-only contract: the cache + the CADD plugin are the whole resource surface. The
    # gnomAD / ClinVar / dbNSFP / SpliceAI / LOFTEE keys are gone because nothing reads them —
    # population frequency and ClinVar now ride in the CSQ from the cache itself. See
    # src/hprv/annotations.py and docs/allele_frequency.md.
    "HPRV_VEP_CACHE": "resources.vep.cache_dir",
    "HPRV_VEP_VERSION": "resources.vep.version",
    "HPRV_VEP_PLUGINS": "resources.vep.plugins_dir",
    "HPRV_VEP_ANNOTATED_VCF": "resources.vep.annotated_vcf",
    "HPRV_CSQ_SELECT": "resources.vep.csq_select",
    "HPRV_VEP_SHARD_BY_CONTIG": "resources.vep.shard_by_contig",
    "HPRV_CADD_SNV": "resources.vep.cadd_snv",
    "HPRV_CADD_INDEL": "resources.vep.cadd_indel",
    # SpliceAI plugin score files (precomputed raw genome-wide splice deltas). Optional, like CADD.
    "HPRV_SPLICEAI_SNV": "resources.vep.spliceai_snv",
    "HPRV_SPLICEAI_INDEL": "resources.vep.spliceai_indel",
    # Calibrated MISSENSE predictors, both VEP plugins. They buy the SCREEN nothing (every missense
    # is IMPACT=MODERATE and selection.py keeps it at the impact rung before any predictor runs —
    # docs/limitations.md #7); their consumer is Step 9's missense tier, which without them can
    # only report an off-label CADD rank. Optional, like CADD.
    "HPRV_REVEL": "resources.vep.revel",
    "HPRV_ALPHAMISSENSE": "resources.vep.alphamissense",
    # ClinVar sites VCF. One of the two bcftools transfers: the VEP cache exposes CLIN_SIG
    # but no CLNREVSTAT, so review status / gold stars are unavailable from the cache at any price.
    "HPRV_CLINVAR_VCF": "resources.clinvar.vcf",
    # The gnomAD v4.1 JOINT sites slim — the second bcftools transfer, and the one that upgrades
    # the rarity oracle from a grpmax point-estimate PROXY to real faf95 (see annotations.frequency).
    "HPRV_GNOMAD_SITES": "resources.gnomad.sites_slim",
    # The run-level frequency oracle (faf95 | grpmax_proxy). Exported so Step 2 can HALT when a
    # faf95 run loses its slim transfer instead of warning and continuing on no oracle at all.
    "HPRV_GNOMAD_ORACLE": "resources.gnomad.oracle",
    "HPRV_CRAM_MAP": "resources.cram_map",
    "HPRV_CRAM_REF": "resources.cram_ref",
    # kraken2 DB for Step-8b non-human-fraction screening (bind-mounted DATA, never baked).
    "HPRV_KRAKEN2_DB": "resources.kraken2_db",
}


def emit_sh(cfg: dict) -> None:
    """Print ``export VAR='value'`` lines; warn (stderr) on unresolved ${...}."""
    unresolved = []
    for var, key in SH_MAP.items():
        val = get(cfg, key, "")
        # YAML booleans must reach the shell as `true`/`false`: Python's str(True) is "True",
        # which no `case ... in 1|true|yes|on)` guard matches — the shipped
        # `resources.vep.shard_by_contig: true` was read by Step 2 as sharding OFF.
        if isinstance(val, bool):
            val = "true" if val else "false"
        val = "" if val is None else str(val)
        if "${" in val:
            # do NOT export an unresolved placeholder — a non-empty '${FOO}' string would
            # shadow the shell-level `:=`/`:-` defaults in common.sh/run_pipeline.sh. Emit
            # empty so those defaults take effect, and warn.
            unresolved.append((var, val))
            print(f"export {var}=''")
            continue
        print(f"export {var}={shlex.quote(val)}")
    if unresolved:
        sys.stderr.write(
            "WARNING: unresolved ${ENV} placeholders (set these env vars before running):\n"
        )
        for var, val in unresolved:
            sys.stderr.write(f"  {var}={val}\n")


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Resolve the hprv pipeline config.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_sh = sub.add_parser("sh", help="emit shell exports for the bash steps")
    p_sh.add_argument("--config", required=True)

    p_get = sub.add_parser("get", help="print one value by dotted key")
    p_get.add_argument("--config", required=True)
    p_get.add_argument("--key", required=True)
    p_get.add_argument("--default", default="")

    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    if args.cmd == "sh":
        emit_sh(cfg)
    elif args.cmd == "get":
        val = get(cfg, args.key, args.default)
        # normalize YAML booleans to lowercase so shell string compares (`!= "false"`) work;
        # Python's str(True/False) is "True"/"False", which silently defeats those guards.
        if isinstance(val, bool):
            val = "true" if val else "false"
        print("" if val is None else val)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
