# syntax=docker/dockerfile:1.7

ARG CLAM_TARGETPLATFORM=linux/amd64
FROM --platform=$CLAM_TARGETPLATFORM mambaorg/micromamba:1.5.10-jammy

ARG CLAM_TARGETPLATFORM
ARG BUILDPLATFORM

LABEL org.opencontainers.image.title="CLAM nuclear WES runtime" \
      org.opencontainers.image.description="Pinned BWA-MEM2, samtools, CNVkit, and offline VEP runtime for nuclear WES analysis" \
      org.opencontainers.image.source="https://github.com/aledavini7/clam-mitogenome-pipeline"

RUN echo "Building CLAM nuclear tools for ${CLAM_TARGETPLATFORM} from ${BUILDPLATFORM:-unknown}"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/usr/local/bin:${PATH}" \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

RUN python3 --version && \
    bwa-mem2 version && \
    samtools --version | head -n 1 && \
    cnvkit.py version && \
    vep --help >/dev/null

WORKDIR /work
