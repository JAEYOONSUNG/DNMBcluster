# syntax=docker/dockerfile:1.7
#
# DNMBcluster multi-stage build
#   Stage 1: compile usearch12 from rcedgar/usearch12 (GPL-3)
#   Stage 2: conda env with MMseqs2/DIAMOND/CD-HIT + Python + R
#   Stage 3: slim runtime with dnmbcluster Python pkg + DNMBcluster R pkg
#

########## Stage 1: usearch12 builder ##########
FROM debian:bookworm-slim AS usearch-build
# TARGETARCH is auto-set by buildx (amd64, arm64, …) so the sed-based
# -march patch can branch per target architecture instead of hard-coding
# an x86 baseline that would break the arm64 image.
ARG TARGETARCH

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential ccache python3 git ca-certificates \
 && rm -rf /var/lib/apt/lists/*

ARG USEARCH_REF=master
RUN git clone --depth 1 --branch ${USEARCH_REF} \
        https://github.com/rcedgar/usearch12.git /src/usearch12

WORKDIR /src/usearch12/src
# Replace -march=native with a portable, arch-appropriate baseline so the
# image is redistributable. usearch12's build_linux.py hard-fails if
# `git status --porcelain` is non-empty, so we commit the edit with a
# local CI identity BEFORE calling the build script — otherwise the
# in-place sed leaves the repo dirty and the script aborts with
# "ERROR -- Uncommited changes".
RUN case "${TARGETARCH}" in \
      amd64) MARCH=x86-64-v3 ;; \
      arm64) MARCH=armv8-a   ;; \
      *)     MARCH=native    ;; \
    esac \
 && sed -i "s/-march=native/-march=${MARCH}/g" build_linux.py \
 && git -C /src/usearch12 -c user.email=ci@local -c user.name=ci \
        commit -am "ci: portable -march=${MARCH}" \
 && python3 build_linux.py

RUN strip /src/usearch12/bin/usearch12 || true

########## Stage 2: conda env solve ##########
FROM mambaorg/micromamba:1.5.8 AS conda-build

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml \
 && micromamba clean --all --yes \
 && find /opt/conda -name '*.pyc' -delete \
 && find /opt/conda -name 'tests' -type d -exec rm -rf {} + 2>/dev/null || true

########## Stage 3: runtime ##########
FROM mambaorg/micromamba:1.5.8 AS runtime
USER root

RUN apt-get update \
 && apt-get install -y --no-install-recommends gosu tini ca-certificates \
 && rm -rf /var/lib/apt/lists/*

COPY --from=conda-build /opt/conda /opt/conda
COPY --from=usearch-build /src/usearch12/bin/usearch12 /opt/conda/bin/usearch12
COPY --from=usearch-build /src/usearch12/LICENSE /usr/share/doc/usearch12/LICENSE

# dnmbcluster Python package
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml /app/pyproject.toml
COPY --chown=$MAMBA_USER:$MAMBA_USER src /app/src
RUN /opt/conda/bin/pip install --no-deps /app

# DNMBcluster R package (absorbed BPGAconverter)
COPY --chown=$MAMBA_USER:$MAMBA_USER R /app/R
RUN /opt/conda/bin/R CMD INSTALL /app/R

# r-eulerr is not shipped via conda-forge for linux-aarch64, so install it
# from CRAN source after R is wired up. Compilers are pulled in temporarily
# and left in place — ortho_euler / euler_upset_combined fall back via
# requireNamespace() if the build fails, so we don't want this step to
# abort the whole image.
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential \
 && (/opt/conda/bin/R -e "install.packages(c('polyclip','eulerr'), repos='https://cloud.r-project.org', Ncpus=2)" || true) \
 && rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /data

COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

ENTRYPOINT ["tini", "--", "/usr/local/bin/docker-entrypoint.sh"]
CMD ["dnmbcluster", "--help"]
