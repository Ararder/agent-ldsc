# syntax=docker/dockerfile:1.7
# agent-ldsc image: every scientific and worker dependency, installed at build time only.
# Build context is the repository root (see .dockerignore); platform linux/amd64.

ARG DEBIAN_IMAGE=debian:bookworm-20261005-slim@sha256:7c7b2c966bc9ee8cedfeef67e0e279108992c77681fa595db4a9d65c06ccc587

# --- pinned upstream sources, verified by commit id ------------------------------------
FROM ${DEBIAN_IMAGE} AS sources
ARG LDSC_COMMIT=aa33296abac9569a6422ee6ba7eb4b902422cc74
ARG LDSR_COMMIT=e880755004946770b0d8a048f970d513fafc16c8
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates \
 && rm -rf /var/lib/apt/lists/*
RUN set -eu; \
    fetch() { git init -q "$1"; git -C "$1" fetch -q --depth 1 "$2" "$3"; git -C "$1" checkout -q FETCH_HEAD; \
              test "$(git -C "$1" rev-parse HEAD)" = "$3"; rm -rf "$1/.git"; }; \
    fetch /src/ldsc https://github.com/bulik/ldsc "$LDSC_COMMIT"; \
    fetch /src/ldsR https://github.com/Ararder/ldsR "$LDSR_COMMIT"

# --- locked environments --------------------------------------------------------------
FROM ${DEBIAN_IMAGE} AS envs
ARG MICROMAMBA_VERSION=2.9.0-0
ARG MICROMAMBA_SHA256=366cd9cd8be14df1ab8ed50352a82111082a36686b2d389fdb79a92c3fafb3e3
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl \
 && rm -rf /var/lib/apt/lists/*
RUN curl -fsSL -o /usr/local/bin/micromamba \
      "https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-linux-64" \
 && echo "${MICROMAMBA_SHA256}  /usr/local/bin/micromamba" | sha256sum -c - \
 && chmod +x /usr/local/bin/micromamba
ENV MAMBA_ROOT_PREFIX=/tmp/mamba
COPY container/locks/ldsc-linux-64.lock container/locks/worker-linux-64.lock container/envs/ldsc-requirements.txt /tmp/locks/
# Explicit lockfiles carry URL + md5 per package; no solver runs here.
RUN micromamba create -y -p /opt/envs/ldsc -f /tmp/locks/ldsc-linux-64.lock \
 && micromamba create -y -p /opt/envs/worker -f /tmp/locks/worker-linux-64.lock \
 && /opt/envs/ldsc/bin/pip install --no-cache-dir --no-deps --only-binary :all: --require-hashes \
      -r /tmp/locks/ldsc-requirements.txt \
 && micromamba clean -a -y
COPY --from=sources /src/ldsc /opt/ldsc
COPY --from=sources /src/ldsR /tmp/ldsR
RUN /opt/envs/worker/bin/R CMD INSTALL --no-docs --no-multiarch --no-test-load /tmp/ldsR \
 && rm -rf /tmp/ldsR /opt/envs/*/pkgs \
 && find /opt/envs -name '__pycache__' -prune -exec rm -rf {} + \
 && find /opt/envs -name '*.a' -delete

# --- final image ------------------------------------------------------------------------
FROM ${DEBIAN_IMAGE}
ARG LDSC_COMMIT=aa33296abac9569a6422ee6ba7eb4b902422cc74
ARG LDSR_COMMIT=e880755004946770b0d8a048f970d513fafc16c8
ARG SOURCE_SHA=unknown
COPY --from=envs /opt/envs /opt/envs
COPY --from=envs /opt/ldsc /opt/ldsc
COPY worker/ /opt/agent-ldsc/worker/
COPY R/ /opt/agent-ldsc/R/
COPY schemas/ /opt/agent-ldsc/schemas/
COPY references/references.lock.json /opt/agent-ldsc/references/references.lock.json
COPY profiles/scientific/ /opt/agent-ldsc/profiles/scientific/
COPY container/bin/ /opt/agent-ldsc/bin/
RUN chmod 0755 /opt/agent-ldsc/bin/* \
 && /opt/envs/worker/bin/python -m compileall -q /opt/agent-ldsc/worker
# Single-threaded BLAS by default; the worker sets task concurrency from the allocation.
ENV PATH=/opt/agent-ldsc/bin:/opt/envs/worker/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    PYTHONPATH=/opt/agent-ldsc/worker \
    PYTHONDONTWRITEBYTECODE=1 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    LC_ALL=C.UTF-8 LANG=C.UTF-8 \
    AGENT_LDSC_HOME=/opt/agent-ldsc \
    AGENT_LDSC_LDSC_COMMIT=${LDSC_COMMIT} \
    AGENT_LDSC_LDSR_COMMIT=${LDSR_COMMIT} \
    AGENT_LDSC_SOURCE_SHA=${SOURCE_SHA}
LABEL org.opencontainers.image.source="https://github.com/Ararder/agent-ldsc" \
      org.opencontainers.image.description="agent-ldsc worker: annotations to ldsR LD scores" \
      org.opencontainers.image.licenses="MIT AND GPL-3.0-only AND LicenseRef-UCSC-liftOver" \
      org.opencontainers.image.revision="${SOURCE_SHA}"
ENTRYPOINT ["agent-ldsc-worker"]
