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
    "HPRV_VEP_CACHE": "resources.vep.cache_dir",
    "HPRV_VEP_VERSION": "resources.vep.version",
    "HPRV_VEP_PLUGINS": "resources.vep.plugins_dir",
    "HPRV_VEP_FASTA": "resources.vep.fasta",
    "HPRV_CADD_SNV": "resources.vep.cadd_snv",
    "HPRV_CADD_INDEL": "resources.vep.cadd_indel",
    "HPRV_DBNSFP": "resources.vep.dbnsfp",
    "HPRV_SPLICEAI_SNV": "resources.vep.spliceai_snv",
    "HPRV_SPLICEAI_INDEL": "resources.vep.spliceai_indel",
    "HPRV_LOFTEE_DATA": "resources.vep.loftee_data",
    "HPRV_GNOMAD_SITES": "resources.gnomad.sites_vcf",
    "HPRV_GNOMAD_AF_TAG": "resources.gnomad.af_tag",
    "HPRV_GNOMAD_GRPMAX_AF_TAG": "resources.gnomad.grpmax_af_tag",
    "HPRV_GNOMAD_FAF95_TAG": "resources.gnomad.faf95_tag",
    "HPRV_GNOMAD_NHOMALT_TAG": "resources.gnomad.nhomalt_tag",
    "HPRV_CRAM_MAP": "resources.cram_map",
    "HPRV_CLINVAR_VCF": "resources.clinvar.vcf",
    "HPRV_CLINVAR_SIG_TAG": "resources.clinvar.sig_tag",
    "HPRV_CLINVAR_REVSTAT_TAG": "resources.clinvar.revstat_tag",
    "HPRV_CLINVAR_SIGCONF_TAG": "resources.clinvar.sigconf_tag",
}


def emit_sh(cfg: dict) -> None:
    """Print ``export VAR='value'`` lines; warn (stderr) on unresolved ${...}."""
    unresolved = []
    for var, key in SH_MAP.items():
        val = get(cfg, key, "")
        val = "" if val is None else str(val)
        if "${" in val:
            unresolved.append((var, val))
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
        print("" if val is None else val)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
