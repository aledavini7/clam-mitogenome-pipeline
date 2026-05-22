# syntax=docker/dockerfile:1.7

ARG CLAM_TARGETPLATFORM=linux/amd64
FROM --platform=$CLAM_TARGETPLATFORM eclipse-temurin:8-jre-jammy

ARG CLAM_TARGETPLATFORM
ARG BUILDPLATFORM

ARG MITOSCAPE_VERSION=1.0
ARG MITOSCAPE_RESOURCES_REV=118d4d5dfcf66774296017130b8886eee0290f49

LABEL org.opencontainers.image.title="CLAM MitoScape runtime" \
      org.opencontainers.image.description="Java 8 runtime for Spark-based MitoScape NUMT classification" \
      org.opencontainers.image.source="https://github.com/aledavini7/clam-mitogenome-pipeline"

RUN echo "Building CLAM MitoScape for ${CLAM_TARGETPLATFORM} from ${BUILDPLATFORM:-unknown}"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        coreutils \
        curl \
        procps \
        tar && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/clam/share/mitoscape && \
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
    rm -rf /tmp/MTClassifierModel.RF.tar && \
    java -version && \
    test -s /opt/clam/share/mitoscape/MitoScapeClassify.jar && \
    test -s /opt/clam/share/mitoscape/mitomap.ld && \
    test -s /opt/clam/share/mitoscape/NUMTs_hg38.txt && \
    test -d /opt/clam/share/mitoscape/MTClassifierModel.RF

WORKDIR /work
