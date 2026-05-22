# syntax=docker/dockerfile:1.7

ARG CLAM_TARGETPLATFORM=linux/amd64
FROM --platform=$CLAM_TARGETPLATFORM python:3.11-slim

ARG CLAM_TARGETPLATFORM
ARG BUILDPLATFORM

LABEL org.opencontainers.image.title="CLAM annotation runtime" \
      org.opencontainers.image.description="Runtime for CLAM mitochondrial variant annotation and summary tables" \
      org.opencontainers.image.source="https://github.com/aledavini7/clam-mitogenome-pipeline"

RUN echo "Building CLAM annotation for ${CLAM_TARGETPLATFORM} from ${BUILDPLATFORM:-unknown}" && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        procps && \
    rm -rf /var/lib/apt/lists/* && \
    ps --version >/dev/null && \
    python3 --version

WORKDIR /work
