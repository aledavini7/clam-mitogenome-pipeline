# syntax=docker/dockerfile:1.7

ARG CLAM_TARGETPLATFORM=linux/amd64
FROM --platform=$CLAM_TARGETPLATFORM mambaorg/micromamba:1.5.10-jammy

ARG CLAM_TARGETPLATFORM
ARG BUILDPLATFORM

LABEL org.opencontainers.image.title="CLAM Manta runtime" \
      org.opencontainers.image.description="Pinned Manta 1.6.0 and HTSlib runtime for optional nuclear WES structural-variant calling" \
      org.opencontainers.image.source="https://github.com/aledavini7/clam-mitogenome-pipeline"

RUN echo "Building CLAM Manta tools for ${CLAM_TARGETPLATFORM} from ${BUILDPLATFORM:-unknown}"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    ln -sf /opt/conda/bin/python2.7 /opt/conda/bin/python2 && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/usr/local/bin:${PATH}" \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

RUN python2 --version && \
    configManta.py --help >/dev/null && \
    bgzip --help >/dev/null 2>&1 && \
    tabix --version >/dev/null

WORKDIR /work
