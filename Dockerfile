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
        pigz \
        python3 \
        python3-matplotlib \
        python3-numpy \
        samtools \
        unzip \
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

# Install plink2 for variant-level QC (missingness, HWE, MAF)
RUN wget -q "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip" -O /tmp/plink2.zip && \
    unzip -o /tmp/plink2.zip plink2 -d /usr/local/bin/ && \
    chmod +x /usr/local/bin/plink2 && \
    rm /tmp/plink2.zip

# Install flashpca2 for fast ancestry PCA computation
RUN apt-get update && apt-get install -y --no-install-recommends \
        libeigen3-dev \
        liblapack-dev \
        libopenblas-dev \
        libboost-program-options-dev \
        git \
    && rm -rf /var/lib/apt/lists/* && \
    cd /tmp && \
    wget -q https://github.com/yixuan/spectra/archive/v0.8.1.tar.gz && \
    tar xzf v0.8.1.tar.gz && \
    git clone https://github.com/gabraham/flashpca.git && \
    cd flashpca && \
    make all EIGEN_INC=/usr/include/eigen3 SPECTRA_INC=/tmp/spectra-0.8.1/include && \
    cp flashpca /usr/local/bin/ && \
    chmod +x /usr/local/bin/flashpca && \
    cd / && rm -rf /tmp/flashpca /tmp/spectra-0.8.1 /tmp/v0.8.1.tar.gz

# Copy pipeline scripts (includes diagnose_qc.py for QC troubleshooting)
COPY scripts/ /opt/scripts/
RUN chmod +x /opt/scripts/*.sh /opt/scripts/*.py
COPY config/ /opt/config/

WORKDIR /data

ENTRYPOINT ["bash"]
