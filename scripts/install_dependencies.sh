#!/usr/bin/env bash
#
# install_dependencies.sh
#
# Install bcftools and required plugins (gtc2vcf, idat2gtc, mocha) for
# processing Illumina IDAT data on HPC or command line without Genome Studio.
#
# Usage:
#   ./install_dependencies.sh [--prefix /path/to/install]
#
set -euo pipefail

BCFTOOLS_VERSION="1.23"
INSTALL_PREFIX="${HOME}/bin"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Install bcftools and required plugins for Illumina IDAT processing.

Options:
  --prefix DIR    Installation directory (default: \$HOME/bin)
  --help          Show this help message
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix) INSTALL_PREFIX="$2"; shift 2 ;;
        --help)   usage ;;
        *)        echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "${INSTALL_PREFIX}"

echo "=== Installing bcftools ${BCFTOOLS_VERSION} and plugins ==="
echo "Install prefix: ${INSTALL_PREFIX}"

WORK_DIR=$(mktemp -d)
trap 'rm -rf "${WORK_DIR}"' EXIT
cd "${WORK_DIR}"

# Download and extract bcftools
echo "Downloading bcftools ${BCFTOOLS_VERSION}..."
wget -q "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"
tar xjf "bcftools-${BCFTOOLS_VERSION}.tar.bz2"
cd "bcftools-${BCFTOOLS_VERSION}"

# Download gtc2vcf plugins (idat2gtc, gtc2vcf, affy2vcf, BAFregress)
echo "Downloading gtc2vcf plugins..."
for f in idat2gtc.c gtc2vcf.c gtc2vcf.h affy2vcf.c BAFregress.c; do
    wget -q -P plugins "https://raw.githubusercontent.com/freeseek/gtc2vcf/master/${f}"
done

# Download mocha plugins
echo "Downloading mocha plugins..."
for f in mocha.h beta_binom.h genome_rules.h mocha.c mochatools.c extendFMT.c; do
    wget -q -P plugins "https://raw.githubusercontent.com/freeseek/mocha/master/${f}"
done

# Compile
echo "Compiling bcftools and plugins..."
make -j"$(nproc)" 2>&1 | tail -5

# Install
echo "Installing to ${INSTALL_PREFIX}..."
cp bcftools "${INSTALL_PREFIX}/"
cp plugins/{fill-tags,fixploidy,idat2gtc,gtc2vcf,affy2vcf,BAFregress,mocha,mochatools,extendFMT}.so "${INSTALL_PREFIX}/"

echo ""
echo "=== Installation complete ==="
echo ""
echo "Add the following to your shell profile:"
echo "  export PATH=\"${INSTALL_PREFIX}:\$PATH\""
echo "  export BCFTOOLS_PLUGINS=\"${INSTALL_PREFIX}\""
echo ""

# Install plink2 for variant-level QC
echo "=== Installing plink2 ==="
PLINK2_URL="https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip"
wget -q "${PLINK2_URL}" -O "${WORK_DIR}/plink2.zip"
unzip -o "${WORK_DIR}/plink2.zip" plink2 -d "${INSTALL_PREFIX}/"
chmod +x "${INSTALL_PREFIX}/plink2"
echo "plink2 installed to ${INSTALL_PREFIX}/plink2"
