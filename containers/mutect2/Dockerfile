# syntax=docker/dockerfile:1.7

ARG CLAM_TARGETPLATFORM=linux/amd64
FROM --platform=$CLAM_TARGETPLATFORM quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0

ARG CLAM_TARGETPLATFORM
ARG BUILDPLATFORM

LABEL org.opencontainers.image.title="CLAM Mutect2 runtime" \
      org.opencontainers.image.description="GATK4 Mutect2 runtime for CLAM mitochondrial variant calling" \
      org.opencontainers.image.source="https://github.com/aledavini7/clam-mitogenome-pipeline"

RUN echo "Building CLAM Mutect2 for ${CLAM_TARGETPLATFORM} from ${BUILDPLATFORM:-unknown}" && \
    gatk --version

WORKDIR /work
