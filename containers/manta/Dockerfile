# syntax=docker/dockerfile:1.7

ARG CLAM_TARGETPLATFORM=linux/amd64
FROM --platform=$CLAM_TARGETPLATFORM quay.io/biocontainers/manta:1.6.0--py27h9948957_6

ARG CLAM_TARGETPLATFORM
ARG BUILDPLATFORM

LABEL org.opencontainers.image.title="CLAM Manta runtime" \
      org.opencontainers.image.description="Pinned Manta 1.6.0 and HTSlib runtime for optional nuclear WES structural-variant calling" \
      org.opencontainers.image.source="https://github.com/aledavini7/clam-mitogenome-pipeline"

RUN echo "Building CLAM Manta tools for ${CLAM_TARGETPLATFORM} from ${BUILDPLATFORM:-unknown}"

# Manta 1.6.0 requires Python 2.7. Use its Bioconda-built runtime instead of
# resolving that legacy interpreter alongside modern HTSlib/OpenSSL packages.
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends tabix && \
    rm -rf /var/lib/apt/lists/*

RUN if ! command -v python2 >/dev/null 2>&1; then \
        python27="$(command -v python2.7)" && \
        test -n "$python27" && \
        ln -sf "$python27" /usr/local/bin/python2; \
    fi

ENV PATH="/opt/conda/bin:/usr/local/bin:${PATH}" \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

RUN python2 --version && \
    configManta.py --help >/dev/null && \
    bgzip --help >/dev/null 2>&1 && \
    tabix --version >/dev/null

WORKDIR /work
