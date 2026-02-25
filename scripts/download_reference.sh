#!/usr/bin/env bash
#
# download_reference.sh
#
# Download reference genome and supporting resources required for IDAT processing.
#
# Usage:
#   ./download_reference.sh --genome GRCh38 --output-dir /path/to/ref
#
set -euo pipefail

GENOME="GRCh38"
OUTPUT_DIR=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Download reference genome and supporting resources.

Options:
  --genome NAME    Genome build: GRCh37 or GRCh38 (default: GRCh38)
  --output-dir DIR Directory to save reference files (required)
  --help           Show this help message
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome)     GENOME="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --help)       usage ;;
        *)            echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output-dir is required" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

if [[ "${GENOME}" == "GRCh38" ]]; then
    FASTA_NAME="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    MAP_URL="https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"
    MAP_NAME="genetic_map_hg38_withX.txt.gz"
    CYTO_URL="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
elif [[ "${GENOME}" == "GRCh37" ]]; then
    FASTA_NAME="human_g1k_v37.fasta"
    FASTA_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
    MAP_URL="https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz"
    MAP_NAME="genetic_map_hg19_withX.txt.gz"
    CYTO_URL="https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz"
else
    echo "Error: Unsupported genome build: ${GENOME}. Use GRCh37 or GRCh38." >&2
    exit 1
fi

echo "=== Downloading reference resources for ${GENOME} ==="
echo "Output directory: ${OUTPUT_DIR}"

# Reference FASTA
if [[ ! -f "${OUTPUT_DIR}/${FASTA_NAME}" ]]; then
    echo "Downloading reference genome..."
    wget -q -O- "${FASTA_URL}" | gzip -d > "${OUTPUT_DIR}/${FASTA_NAME}"
    echo "Indexing reference genome..."
    samtools faidx "${OUTPUT_DIR}/${FASTA_NAME}"
else
    echo "Reference genome already exists: ${OUTPUT_DIR}/${FASTA_NAME}"
fi

# BWA index (needed for manifest realignment)
if [[ ! -f "${OUTPUT_DIR}/${FASTA_NAME}.bwt" ]]; then
    echo "Creating BWA index (this may take a while)..."
    if command -v bwa &>/dev/null; then
        bwa index "${OUTPUT_DIR}/${FASTA_NAME}"
    else
        echo "Warning: bwa not found. BWA index not created (needed only for manifest realignment)."
    fi
else
    echo "BWA index already exists."
fi

# Genetic map
if [[ ! -f "${OUTPUT_DIR}/${MAP_NAME}" ]]; then
    echo "Downloading genetic map..."
    wget -q -P "${OUTPUT_DIR}" "${MAP_URL}"
else
    echo "Genetic map already exists."
fi

# Cytoband file
if [[ ! -f "${OUTPUT_DIR}/cytoBand.txt.gz" ]]; then
    echo "Downloading cytoband file..."
    wget -q -O "${OUTPUT_DIR}/cytoBand.txt.gz" "${CYTO_URL}"
else
    echo "Cytoband file already exists."
fi

echo ""
echo "=== Reference download complete ==="
echo "Files in ${OUTPUT_DIR}:"
ls -lh "${OUTPUT_DIR}"
