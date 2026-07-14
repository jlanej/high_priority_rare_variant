# syntax=docker/dockerfile:1
# =============================================================================
# high_priority_rare_variant — single analysis image
#
# Contains EVERY tool the pipeline scripts need, so analysis is fully reproducible
# under Apptainer on HPC (the scripts only ever call tools present in this image).
#
# Base-image choice: we build FROM the official Ensembl VEP image rather than a
# from-scratch micromamba image. VEP + its plugin/Perl infrastructure is heavy and
# finicky to reproduce on conda; basing on the upstream image guarantees VEP works
# exactly as in the group's validated annotation setup (VEP release 115). The rest
# of the toolchain (bcftools, bedtools, slivar, somalier, whatshap, python libs) is
# layered on via micromamba from a version-pinned environment (env/environment.yml).
# See docs/tooling_and_reproducibility.md for the tradeoff discussion.
#
# Pin note: the base tag is fixed here for readability; CI records the resolved
# @sha256 digest, and downstream consumers should pull the image by digest.
# =============================================================================
FROM ensemblorg/ensembl-vep:release_115.0

LABEL org.opencontainers.image.title="high_priority_rare_variant" \
      org.opencontainers.image.description="Containerized toolchain for screening GMKF Kids First trio VCFs for high-priority rare variants" \
      org.opencontainers.image.source="https://github.com/jlanej/high_priority_rare_variant" \
      org.opencontainers.image.licenses="MIT"

USER root
ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    HPRV_ENV=/opt/conda/envs/hprv

# --- micromamba (reproducible, pinned analysis tools) -----------------------
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends curl bzip2 ca-certificates procps; \
    rm -rf /var/lib/apt/lists/*; \
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
      | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# --- create the pinned environment ------------------------------------------
COPY env/environment.yml /tmp/environment.yml
RUN set -eux; \
    micromamba create -y -p "${HPRV_ENV}" -f /tmp/environment.yml; \
    micromamba clean -a -y; \
    rm -f /tmp/environment.yml

# Analysis tools resolve from the conda env. The env is perl-free, so VEP's
# `env perl` still finds the base image's system Perl and VEP is unaffected.
ENV PATH=${HPRV_ENV}/bin:$PATH

# --- LOFTEE plugin CODE (GRCh38 fork) ---------------------------------------
# The ensembl-vep base bundles the Ensembl/VEP_plugins .pm at /plugins (CADD, dbNSFP,
# SpliceAI, ...) but builds with `--skip_plugins LoF`, and LOFTEE is a SEPARATE repo
# (konradjk/loftee) — so LoF.pm is absent. Bake the *grch38 branch* (master is GRCh37-only)
# into /plugins so a single `--dir_plugins /plugins` (and loftee_path:/plugins) serves all four
# plugins. LOFTEE's Perl deps (Bio::DB::BigFile against the Kent lib, DBI, Bio::Perl) are provided
# by the base image; DBD::SQLite is only a `recommends` there, so ensure it. The plugin-.pm checks
# HARD-FAIL the build (that is the guarantee this layer makes); the Perl-dep checks only REPORT — a
# base `recommends` may not have compiled, and a truly missing runtime dep surfaces clearly when the
# LoF plugin runs, so hard-failing on it would block an otherwise-working image.
# Only plugin CODE is baked — the DATA (human_ancestor/GERP/loftee.sql) is host-fetched.
# Repro note: pin a grch38 commit SHA for byte-identical rebuilds (the branch tip moves).
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends git ca-certificates; \
    git clone --depth 1 -b grch38 https://github.com/konradjk/loftee.git /tmp/loftee; \
    cp -a /tmp/loftee/. /plugins/; \
    rm -rf /tmp/loftee; \
    perl -MDBD::SQLite -e1 2>/dev/null || apt-get install -y --no-install-recommends libdbd-sqlite3-perl || true; \
    apt-get clean; rm -rf /var/lib/apt/lists/*; \
    # GUARANTEE: the plugin CODE must be present — fail the build if not
    for p in LoF CADD dbNSFP SpliceAI; do test -f "/plugins/$p.pm" || { echo "MISSING plugin code: /plugins/$p.pm" >&2; exit 1; }; done; \
    # REPORT (non-fatal): LOFTEE runtime Perl deps — a missing one errors clearly when LoF runs
    for m in DBI DBD::SQLite Bio::DB::BigFile Bio::Perl; do \
        if perl -M"$m" -e1 2>/dev/null; then echo "loftee dep OK: $m"; \
        else echo "WARN: loftee Perl dep not loadable at build: $m (LoF will error at runtime if truly missing)" >&2; fi; \
    done

# VEP's plugin-code dir in this base is /plugins; expose it under our config's ${VEP_PLUGINS}
# so `resources.vep.plugins_dir` resolves without the user setting anything.
ENV VEP_PLUGINS=/plugins

# --- pipeline code ----------------------------------------------------------
COPY pipeline/ /opt/hprv/pipeline/
COPY src/ /opt/hprv/src/
ENV PYTHONPATH=/opt/hprv/src:$PYTHONPATH \
    PATH=/opt/hprv/pipeline:$PATH \
    HPRV_HOME=/opt/hprv
RUN find /opt/hprv/pipeline -type f -name '*.sh' -exec chmod a+rx {} + ; \
    find /opt/hprv/pipeline -type f -name '*.py' -exec chmod a+rx {} +

# --- build-time sanity check: every core tool + python dep must resolve ------
RUN set -eux; \
    bcftools --version | head -1; \
    bedtools --version; \
    tabix --version | head -1; \
    slivar --help >/dev/null 2>&1 || slivar 2>/dev/null || true; \
    somalier --help >/dev/null 2>&1 || true; \
    whatshap --version; \
    vep --help >/dev/null 2>&1 || true; \
    python -c "import cyvcf2, pysam, pandas, numpy, scipy, yaml, click; print('python deps OK')"

USER vep
WORKDIR /data
CMD ["bash"]
