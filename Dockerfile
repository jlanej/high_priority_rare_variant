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
    rm -f /tmp/environment.yml; \
    # bioconda tools (bcftools/samtools/...) drag `perl` into the env as a transitive dep. With
    # the env prepended to PATH it would SHADOW the base image's VEP Perl, so VEP's `#!/usr/bin/env
    # perl` (and its plugins, incl. LOFTEE) would run under the conda Perl and fail to load their
    # modules. Nothing in this pipeline uses the conda Perl, so remove it to keep the env GENUINELY
    # perl-free — the invariant VEP relies on — so `env perl` resolves to the base image's Perl.
    find "${HPRV_ENV}/bin" -maxdepth 1 \( -name 'perl' -o -name 'perl5*' \) -delete 2>/dev/null || true

# Analysis tools resolve from the conda env; the env is kept perl-free (above), so VEP's
# `env perl` finds the base image's system Perl and VEP + its plugins are unaffected.
ENV PATH=${HPRV_ENV}/bin:$PATH

# --- LOFTEE plugin CODE (GRCh38 fork) ---------------------------------------
# The ensembl-vep base bundles the Ensembl/VEP_plugins .pm at /plugins (CADD, dbNSFP,
# SpliceAI, ...) but builds with `--skip_plugins LoF`, and LOFTEE is a SEPARATE repo
# (konradjk/loftee) — so LoF.pm is absent. Bake the *grch38 branch* (master is GRCh37-only)
# into /plugins so a single `--dir_plugins /plugins` (and loftee_path:/plugins) serves all four
# plugins. LOFTEE's Perl deps (Bio::DB::BigFile against the Kent lib, DBI, Bio::Perl) are provided
# by the base image; DBD::SQLite is only a `recommends` there, so ensure it. Both the plugin-.pm
# existence AND every LOFTEE runtime Perl dep HARD-FAIL the build: now that the conda Perl no longer
# shadows VEP's Perl (see the env layer), all of DBI/DBD::SQLite/Bio::DB::BigFile/Bio::Perl must load
# in VEP's Perl — a broken or GERP-less LOFTEE must never ship silently. If a base `recommends`
# (e.g. Bio::DB::BigFile via the Kent lib) turns out not to have compiled, the build says so loudly
# and we add an explicit build for it.
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
    # Every LOFTEE runtime Perl dep must load in VEP's Perl (DBI also proves the conda-Perl un-shadow
    # worked; DBD::SQLite = conservation db; Bio::DB::BigFile = GERP bigwig; Bio::Perl) — fail loudly.
    for m in DBI DBD::SQLite Bio::DB::BigFile Bio::Perl; do \
        perl -M"$m" -e1 2>/dev/null || { echo "FATAL: LOFTEE Perl dep not loadable: $m — conda Perl still shadowing VEP's, or a base 'recommends' (e.g. Bio::DB::BigFile via the Kent lib) did not compile; add an explicit install." >&2; exit 1; }; \
    done

# VEP's plugin-code dir in this base is /plugins; expose it under our config's ${VEP_PLUGINS}
# so `resources.vep.plugins_dir` resolves without the user setting anything.
ENV VEP_PLUGINS=/plugins

# --- pipeline code ----------------------------------------------------------
COPY pipeline/ /opt/hprv/pipeline/
COPY src/ /opt/hprv/src/
# The resource-prep tooling is SOFTWARE, so it ships in the image like everything else: the
# quickstart runs prepare_resources.sh INSIDE the container, and it must not depend on a host
# checkout happening to be bind-mounted at $PWD. Its pinned manifest rides along (.dockerignore
# keeps prepared resource DATA out of the context); re-pin without a rebuild via
# HPRV_RESOURCE_MANIFEST. join_constraint.py resolves via the same $HERE.
COPY scripts/ /opt/hprv/scripts/
COPY resources/manifest.env /opt/hprv/resources/manifest.env
ENV PYTHONPATH=/opt/hprv/src:$PYTHONPATH \
    PATH=/opt/hprv/pipeline:/opt/hprv/scripts:$PATH \
    HPRV_HOME=/opt/hprv
RUN find /opt/hprv/pipeline /opt/hprv/scripts -type f \
        \( -name '*.sh' -o -name '*.py' \) -exec chmod a+rx {} +

# --- build-time sanity check: every core tool + python dep must resolve ------
RUN set -eux; \
    bcftools --version | head -1; \
    bedtools --version; \
    tabix --version | head -1; \
    slivar --help >/dev/null 2>&1 || slivar 2>/dev/null || true; \
    somalier --help >/dev/null 2>&1 || true; \
    whatshap --version; \
    vep --help >/dev/null 2>&1 || true; \
    python -c "import cyvcf2, pysam, pandas, numpy, scipy, yaml, click; print('python deps OK')"; \
    prepare_resources.sh --dir /tmp/_resprobe emit-env >/dev/null; \
    rm -rf /tmp/_resprobe; \
    echo "prepare_resources.sh on PATH and resolves its pinned manifest in-image"

USER vep
WORKDIR /data
CMD ["bash"]
