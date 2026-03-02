#!/usr/bin/env bash
#
# quickstart.sh
#
# WDL quickstart for the Illumina IDAT processing pipeline.
# Downloads Cromwell, creates a Cromwell config for Apptainer,
# and generates a template inputs JSON — ready to run with minimal setup.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CROMWELL_VERSION="87"
OUTPUT_DIR="."
SKIP_CROMWELL="false"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

WDL quickstart for the Illumina IDAT processing pipeline.
Downloads Cromwell and creates everything needed to run the pipeline via WDL.

Options:
  --output-dir DIR        Directory for Cromwell setup files (default: current directory)
  --cromwell-version VER  Cromwell version to download (default: ${CROMWELL_VERSION})
  --skip-cromwell         Skip downloading Cromwell (if already present)
  --help                  Show this help message

This script will:
  1. Check for Java 11+ and Apptainer
  2. Download Cromwell ${CROMWELL_VERSION} JAR
  3. Copy a Cromwell config for Apptainer (cromwell_apptainer.conf)
  4. Generate a template inputs JSON (inputs_template.json)
  5. Print the command to run the pipeline

Requirements:
  - Java 11+     (for running Cromwell)
  - Apptainer    (for running the container; 'singularity' also works)
  - wget or curl (for downloading Cromwell)

Examples:
  # Run quickstart in the current directory
  bash wdl/quickstart.sh

  # Specify an output directory for setup files
  bash wdl/quickstart.sh --output-dir ./cromwell_run

  # After filling in inputs_template.json, run the pipeline:
  java -Dconfig.file=cromwell_apptainer.conf \\
    -jar cromwell-${CROMWELL_VERSION}.jar run \\
    wdl/illumina_idat_processing.wdl \\
    --inputs inputs.json
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)       OUTPUT_DIR="$2"; shift 2 ;;
        --cromwell-version) CROMWELL_VERSION="$2"; shift 2 ;;
        --skip-cromwell)    SKIP_CROMWELL="true"; shift ;;
        --help)             usage ;;
        *)                  echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

CROMWELL_JAR="cromwell-${CROMWELL_VERSION}.jar"
CROMWELL_URL="https://github.com/broadinstitute/cromwell/releases/download/${CROMWELL_VERSION}/${CROMWELL_JAR}"

mkdir -p "${OUTPUT_DIR}"

echo "======================================================"
echo "  Illumina IDAT Processing — WDL Quickstart"
echo "======================================================"
echo ""

# ---------------------------------------------------------------
# Check prerequisites
# ---------------------------------------------------------------
echo "Checking prerequisites..."

if ! command -v java &>/dev/null; then
    echo "ERROR: Java not found. Install Java 11+ before running Cromwell." >&2
    echo "  Ubuntu/Debian: sudo apt-get install openjdk-17-jre-headless" >&2
    echo "  RHEL/Rocky:    sudo yum install java-17-openjdk-headless" >&2
    exit 1
fi
JAVA_VERSION=$(java -version 2>&1 | awk -F '"' '/version/ {print $2}' | cut -d. -f1)
echo "  Java ${JAVA_VERSION}: OK"

if command -v apptainer &>/dev/null; then
    echo "  Apptainer: OK"
elif command -v singularity &>/dev/null; then
    echo "  Singularity: OK (update cromwell_apptainer.conf to use 'singularity' instead of 'apptainer')"
else
    echo "  WARNING: apptainer/singularity not found."
    echo "    Install Apptainer: https://apptainer.org/docs/admin/main/installation.html"
fi

echo ""

# ---------------------------------------------------------------
# Download Cromwell
# ---------------------------------------------------------------
CROMWELL_DEST="${OUTPUT_DIR}/${CROMWELL_JAR}"

if [[ "${SKIP_CROMWELL}" == "true" ]] || [[ -f "${CROMWELL_DEST}" ]]; then
    echo "Cromwell: ${CROMWELL_DEST} (already present, skipping download)"
else
    echo "Downloading Cromwell ${CROMWELL_VERSION}..."
    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "${CROMWELL_DEST}" "${CROMWELL_URL}"
    elif command -v curl &>/dev/null; then
        curl -fL --progress-bar -o "${CROMWELL_DEST}" "${CROMWELL_URL}"
    else
        echo "ERROR: Neither wget nor curl found. Install one to download Cromwell." >&2
        exit 1
    fi
    echo "  Downloaded: ${CROMWELL_DEST}"
fi

echo ""

# ---------------------------------------------------------------
# Copy Cromwell Apptainer config
# ---------------------------------------------------------------
CONF_SRC="${SCRIPT_DIR}/cromwell_apptainer.conf"
CONF_DEST="${OUTPUT_DIR}/cromwell_apptainer.conf"

if [[ -f "${CONF_SRC}" ]]; then
    cp "${CONF_SRC}" "${CONF_DEST}"
    echo "Cromwell config: ${CONF_DEST}"
else
    echo "WARNING: ${CONF_SRC} not found — skipping config copy." >&2
fi

# ---------------------------------------------------------------
# Generate template inputs JSON
# ---------------------------------------------------------------
INPUTS_DEST="${OUTPUT_DIR}/inputs_template.json"
cat > "${INPUTS_DEST}" << 'INPUTS_JSON'
{
  "illumina_idat_processing.idat_files": [
    "/path/to/sample1_Grn.idat",
    "/path/to/sample1_Red.idat",
    "/path/to/sample2_Grn.idat",
    "/path/to/sample2_Red.idat"
  ],
  "illumina_idat_processing.bpm_file":      "/path/to/array.bpm",
  "illumina_idat_processing.egt_file":      "/path/to/array.egt",
  "illumina_idat_processing.csv_file":      "/path/to/array.csv",
  "illumina_idat_processing.ref_fasta":     "/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
  "illumina_idat_processing.ref_fasta_fai": "/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
  "illumina_idat_processing.min_call_rate": 0.97,
  "illumina_idat_processing.max_lrr_sd":    0.35,
  "illumina_idat_processing.threads":       4,
  "illumina_idat_processing.skip_stage2":   false,
  "illumina_idat_processing.docker":        "ghcr.io/jlanej/illumina_idat_processing:main"
}
INPUTS_JSON

echo "Inputs template: ${INPUTS_DEST}"
echo ""
echo "======================================================"
echo "  Setup Complete"
echo "======================================================"
echo ""
echo "Next steps:"
echo ""
echo "1. Download manifests and reference (if not already available):"
echo "   # Using the pre-built container with Apptainer:"
echo "   apptainer exec docker://ghcr.io/jlanej/illumina_idat_processing:main \\"
echo "     bash /opt/scripts/download_manifests.sh --array-name GSA-24v3-0_A1 --output-dir manifests/"
echo "   apptainer exec docker://ghcr.io/jlanej/illumina_idat_processing:main \\"
echo "     bash /opt/scripts/download_reference.sh --genome GRCh38 --output-dir reference/"
echo ""
echo "2. Fill in the inputs JSON with your data paths:"
echo "   cp ${INPUTS_DEST} ${OUTPUT_DIR}/inputs.json"
echo "   # Edit ${OUTPUT_DIR}/inputs.json"
echo ""
echo "3. Run the pipeline:"
WDL_PATH="${SCRIPT_DIR}/illumina_idat_processing.wdl"
echo "   java -Dconfig.file=${CONF_DEST} \\"
echo "     -jar ${CROMWELL_DEST} run \\"
echo "     ${WDL_PATH} \\"
echo "     --inputs ${OUTPUT_DIR}/inputs.json"
echo ""
echo "The pipeline will pull and run the container:"
echo "  ghcr.io/jlanej/illumina_idat_processing:main"
echo ""
echo "For HPC schedulers (SLURM, SGE), see:"
echo "  https://cromwell.readthedocs.io/en/stable/backends/HPC/"
