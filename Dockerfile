FROM ubuntu:22.04

LABEL maintainer="illumina_idat_processing"
LABEL description="Illumina IDAT processing pipeline with bcftools and plugins for HPC"

ARG BCFTOOLS_VERSION=1.23
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        bwa \
        bzip2 \
        ca-certificates \
        curl \
        g++ \
        gcc \
        gzip \
        libbz2-dev \
        libcurl4-openssl-dev \
        libdeflate-dev \
        liblzma-dev \
        libssl-dev \
        make \
        python3 \
        python3-numpy \
        samtools \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Build bcftools with gtc2vcf and mocha plugins
WORKDIR /tmp/build
RUN wget -q "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" && \
    tar xjf "bcftools-${BCFTOOLS_VERSION}.tar.bz2" && \
    cd "bcftools-${BCFTOOLS_VERSION}" && \
    for f in idat2gtc.c gtc2vcf.c gtc2vcf.h affy2vcf.c BAFregress.c; do \
        wget -q -P plugins "https://raw.githubusercontent.com/freeseek/gtc2vcf/master/${f}"; \
    done && \
    for f in mocha.h beta_binom.h genome_rules.h mocha.c mochatools.c extendFMT.c; do \
        wget -q -P plugins "https://raw.githubusercontent.com/freeseek/mocha/master/${f}"; \
    done && \
    make -j"$(nproc)" && \
    make install && \
    cp plugins/*.so /usr/local/libexec/bcftools/ && \
    cd / && rm -rf /tmp/build

ENV BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools

# Copy pipeline scripts
COPY scripts/ /opt/scripts/
COPY config/ /opt/config/

WORKDIR /data

ENTRYPOINT ["bash"]
