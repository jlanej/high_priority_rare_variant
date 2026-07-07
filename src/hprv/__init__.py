"""high_priority_rare_variant — screen GMKF Kids First trio VCFs for high-priority rare variants.

Config-driven analysis library used by the ``pipeline/`` step scripts. Methods and the
canonical parameter defaults are documented in ``docs/`` (see ``docs/README.md``).

This is a PUBLIC package: it must contain no hardcoded filesystem paths, no sample or
subject identifiers, and no controlled-access (dbGaP/PHI) data. All inputs, resources,
and thresholds are supplied via configuration.
"""

__version__ = "0.1.0"
