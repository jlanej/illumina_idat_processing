FROM ubuntu:22.04

LABEL maintainer="illumina_idat_processing"
LABEL description="Illumina IDAT processing pipeline with bcftools and plugins for HPC"

ARG BCFTOOLS_VERSION=1.23
ARG DEBIAN_FRONTEND=noninteractive

# --- Pinned dependency versions for reproducibility ---
# plink2: v2.00a6.33 (March 2026)
ARG PLINK2_VERSION=v2.0.0-a.6.33
# flashpca: pinned to latest commit (2020-09-10)
ARG FLASHPCA_COMMIT=b8044f13607a072125828547684fde8b081d6191
# bcftools plugins — pinned commit hashes:
#   gtc2vcf: 2026-02-26
ARG GTC2VCF_COMMIT=cc4898976c11dda6c7bfb3473d13afea11c48a1c
#   mocha: 2025-08-22
ARG MOCHA_COMMIT=95686b7b65f53a490513be76bb120e5fc20a8bcf
#   score (liftover plugin): 2026-02-02
ARG SCORE_COMMIT=909d23019e19aeadf3bf6fe1407fd6afc094592a

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
        python3-pip \
        samtools \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install peddy for pedigree/sex/ancestry QC.
# Pin NumPy below 1.25 to stay compatible with Ubuntu 22.04 scipy/matplotlib
# ABI expectations and avoid runtime scipy numpy-version warnings in peddy.
RUN pip3 install "numpy<1.25" peddy

# Build bcftools with gtc2vcf, mocha, and liftover plugins (pinned commits)
WORKDIR /tmp/build
RUN wget -q "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" && \
    tar xjf "bcftools-${BCFTOOLS_VERSION}.tar.bz2" && \
    cd "bcftools-${BCFTOOLS_VERSION}" && \
    for f in idat2gtc.c gtc2vcf.c gtc2vcf.h affy2vcf.c BAFregress.c; do \
        wget -q -P plugins "https://raw.githubusercontent.com/freeseek/gtc2vcf/${GTC2VCF_COMMIT}/${f}"; \
    done && \
    for f in mocha.h beta_binom.h genome_rules.h mocha.c mochatools.c extendFMT.c; do \
        wget -q -P plugins "https://raw.githubusercontent.com/freeseek/mocha/${MOCHA_COMMIT}/${f}"; \
    done && \
    wget -q -P plugins "https://raw.githubusercontent.com/freeseek/score/${SCORE_COMMIT}/liftover.c" && \
    make -j"$(nproc)" && \
    make install && \
    cp plugins/*.so /usr/local/libexec/bcftools/ && \
    cd / && rm -rf /tmp/build

ENV BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools

# Install plink2 for variant-level QC (missingness, HWE, MAF) — pinned version
RUN wget -q "https://github.com/chrchang/plink-ng/releases/download/${PLINK2_VERSION}/plink2_linux_x86_64.zip" -O /tmp/plink2.zip && \
    unzip -o /tmp/plink2.zip plink2 -d /usr/local/bin/ && \
    chmod +x /usr/local/bin/plink2 && \
    rm /tmp/plink2.zip

# Install flashpca2 for fast ancestry PCA computation — pinned commit
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
    git checkout "${FLASHPCA_COMMIT}" && \
    make all EIGEN_INC=/usr/include/eigen3 SPECTRA_INC=/tmp/spectra-0.8.1/include && \
    cp flashpca /usr/local/bin/ && \
    chmod +x /usr/local/bin/flashpca && \
    cd / && rm -rf /tmp/flashpca /tmp/spectra-0.8.1 /tmp/v0.8.1.tar.gz

# Copy pipeline scripts (includes diagnose_qc.py for QC troubleshooting)
COPY scripts/ /opt/scripts/
RUN chmod +x /opt/scripts/*.sh /opt/scripts/*.py
COPY config/ /opt/config/
COPY tests/data/ /opt/tests/data/

WORKDIR /data

ENTRYPOINT ["bash"]
