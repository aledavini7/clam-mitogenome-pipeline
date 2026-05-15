# syntax=docker/dockerfile:1.7

ARG CLAM_TARGETPLATFORM=linux/amd64
FROM --platform=$CLAM_TARGETPLATFORM mambaorg/micromamba:1.5.10-jammy

ARG CLAM_TARGETPLATFORM
ARG BUILDPLATFORM

ARG MITOSCAPE_VERSION=1.0
ARG MITOSCAPE_RESOURCES_REV=118d4d5dfcf66774296017130b8886eee0290f49
ARG HAPLOGREP3_VERSION=3.2.2
ARG MUTSERVE_VERSION=2.0.3
ARG HAPLOGREP3_TREE=phylotree-rcrs@17.2

LABEL org.opencontainers.image.title="CLAM core runtime" \
      org.opencontainers.image.description="Runtime for CLAM mitochondrial genome pipeline core processes" \
      org.opencontainers.image.source="https://github.com/aledavini7/clam-mitogenome-pipeline"

RUN echo "Building CLAM core for ${CLAM_TARGETPLATFORM} from ${BUILDPLATFORM:-unknown}"

USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        coreutils \
        curl \
        findutils \
        gawk \
        grep \
        procps \
        sed \
        tar \
        unzip && \
    rm -rf /var/lib/apt/lists/*

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

USER $MAMBA_USER
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

USER root

ENV PATH="/opt/conda/bin:/usr/local/bin:${PATH}" \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

RUN mkdir -p \
        /opt/clam/haplogrep3 \
        /opt/clam/mutserve \
        /opt/clam/share/mitoscape && \
    curl -fsSL \
        "https://github.com/larryns/MitoScape/releases/download/v${MITOSCAPE_VERSION}/MitoScape-${MITOSCAPE_VERSION}.jar" \
        -o "/opt/clam/share/mitoscape/MitoScape-${MITOSCAPE_VERSION}.jar" && \
    ln -s "MitoScape-${MITOSCAPE_VERSION}.jar" /opt/clam/share/mitoscape/MitoScapeClassify.jar && \
    curl -fsSL \
        "https://raw.githubusercontent.com/larryns/MitoScape/${MITOSCAPE_RESOURCES_REV}/src/universal/mitomap.ld" \
        -o /opt/clam/share/mitoscape/mitomap.ld && \
    curl -fsSL \
        "https://raw.githubusercontent.com/larryns/MitoScape/${MITOSCAPE_RESOURCES_REV}/src/universal/NUMTs_hg38.txt" \
        -o /opt/clam/share/mitoscape/NUMTs_hg38.txt && \
    curl -fsSL \
        "https://raw.githubusercontent.com/larryns/MitoScape/${MITOSCAPE_RESOURCES_REV}/src/universal/MTClassifierModel.RF.tar" \
        -o /tmp/MTClassifierModel.RF.tar && \
    tar -xf /tmp/MTClassifierModel.RF.tar \
        -C /opt/clam/share/mitoscape && \
    curl -fsSL \
        "https://github.com/seppinho/mutserve/releases/download/v${MUTSERVE_VERSION}/mutserve.zip" \
        -o /tmp/mutserve.zip && \
    unzip -q /tmp/mutserve.zip -d /opt/clam/mutserve && \
    chmod +x /opt/clam/mutserve/mutserve && \
    curl -fsSL \
        "https://github.com/genepi/haplogrep3/releases/download/v${HAPLOGREP3_VERSION}/haplogrep3-${HAPLOGREP3_VERSION}-linux.zip" \
        -o /tmp/haplogrep3.zip && \
    unzip -q /tmp/haplogrep3.zip -d /opt/clam/haplogrep3 && \
    chmod +x /opt/clam/haplogrep3/haplogrep3 && \
    rm -rf /tmp/MTClassifierModel.RF.tar /tmp/mutserve.zip /tmp/haplogrep3.zip

COPY bin/haplogrep3 bin/haplogrep bin/mutserve /usr/local/bin/

RUN chmod +x /usr/local/bin/haplogrep3 /usr/local/bin/haplogrep /usr/local/bin/mutserve && \
    cd /opt/clam/haplogrep3 && \
    ./haplogrep3 install-tree "${HAPLOGREP3_TREE}" && \
    gsnap --version && \
    samtools --version | head -n 1 && \
    bgzip --version | head -n 1 && \
    mutserve --version && \
    haplogrep3 trees

WORKDIR /work
