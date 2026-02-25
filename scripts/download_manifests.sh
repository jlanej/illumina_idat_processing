#!/usr/bin/env bash
#
# download_manifests.sh
#
# Download Illumina array manifest files (BPM, EGT, CSV) required for
# IDAT processing. Can auto-detect the array type from IDAT files.
#
# Usage:
#   ./download_manifests.sh --array-name GSA-24v3-0_A1 --output-dir manifests/
#   ./download_manifests.sh --detect-from /path/to/idats --output-dir manifests/
#   ./download_manifests.sh --bpm-url URL --egt-url URL --csv-url URL --output-dir manifests/
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/../config"

ARRAY_NAME=""
DETECT_FROM=""
BPM_URL=""
EGT_URL=""
CSV_URL=""
OUTPUT_DIR=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Download Illumina array manifest files (BPM, EGT, CSV).

Options:
  --array-name NAME       Known array name (see config/manifest_urls.tsv)
  --detect-from DIR       Auto-detect array type from IDAT files in DIR
                          (requires bcftools with gtc2vcf plugin)
  --bpm-url URL           Direct URL for BPM manifest file
  --egt-url URL           Direct URL for EGT cluster file
  --csv-url URL           Direct URL for CSV manifest file
  --output-dir DIR        Directory to save manifest files (required)
  --help                  Show this help message

If --array-name is provided, URLs are looked up from config/manifest_urls.tsv.
If --detect-from is provided, the array type is auto-detected from IDAT files.
Direct URLs (--bpm-url, --egt-url, --csv-url) override any looked-up URLs.
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --array-name)  ARRAY_NAME="$2"; shift 2 ;;
        --detect-from) DETECT_FROM="$2"; shift 2 ;;
        --bpm-url)     BPM_URL="$2"; shift 2 ;;
        --egt-url)     EGT_URL="$2"; shift 2 ;;
        --csv-url)     CSV_URL="$2"; shift 2 ;;
        --output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
        --help)        usage ;;
        *)             echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output-dir is required" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# Auto-detect array type from IDAT files
if [[ -n "${DETECT_FROM}" ]]; then
    if ! command -v bcftools &>/dev/null; then
        echo "Error: bcftools is required for array detection. Run install_dependencies.sh first." >&2
        exit 1
    fi
    echo "Detecting array type from IDAT files in ${DETECT_FROM}..."

    # Use gtc2vcf plugin to identify chip type from IDAT files
    IDAT_FILE=$(find "${DETECT_FROM}" -name "*_Grn.idat" -o -name "*_Grn.idat.gz" | head -1)
    if [[ -z "${IDAT_FILE}" ]]; then
        echo "Error: No green IDAT files found in ${DETECT_FROM}" >&2
        exit 1
    fi
    DETECTED=$(bcftools +gtc2vcf -i -g "$(dirname "${IDAT_FILE}")" 2>/dev/null | \
        awk -F'\t' 'NR>1 {print $NF}' | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')
    if [[ -n "${DETECTED}" && -z "${ARRAY_NAME}" ]]; then
        ARRAY_NAME="${DETECTED}"
        echo "Detected array type: ${ARRAY_NAME}"
    fi
fi

# Look up URLs from config file
MANIFEST_FILE="${CONFIG_DIR}/manifest_urls.tsv"
if [[ -n "${ARRAY_NAME}" && -f "${MANIFEST_FILE}" ]]; then
    echo "Looking up manifest URLs for ${ARRAY_NAME}..."
    MATCH=$(grep -v '^#' "${MANIFEST_FILE}" | awk -F'\t' -v name="${ARRAY_NAME}" '$1==name {print}')
    if [[ -n "${MATCH}" ]]; then
        LOOKED_UP_URL=$(echo "${MATCH}" | cut -f2)
        # Check if URLs are Illumina Support pages (not direct download links)
        if [[ "${LOOKED_UP_URL}" == *"support.illumina.com"* ]]; then
            echo ""
            echo "NOTE: Illumina has moved manifest downloads to their Support portal." >&2
            echo "Direct download URLs are no longer available for ${ARRAY_NAME}." >&2
            echo "Please download BPM, EGT, and CSV files from:" >&2
            echo "  ${LOOKED_UP_URL}" >&2
            echo "" >&2
            echo "Then re-run with --bpm-url, --egt-url, and --csv-url pointing to the" >&2
            echo "downloaded files, or place them in the output directory." >&2
            exit 1
        fi
        [[ -z "${BPM_URL}" ]] && BPM_URL=$(echo "${MATCH}" | cut -f2)
        [[ -z "${EGT_URL}" ]] && EGT_URL=$(echo "${MATCH}" | cut -f3)
        [[ -z "${CSV_URL}" ]] && CSV_URL=$(echo "${MATCH}" | cut -f4)
    else
        echo "Warning: Array '${ARRAY_NAME}' not found in ${MANIFEST_FILE}" >&2
        echo "Available arrays:" >&2
        grep -v '^#' "${MANIFEST_FILE}" | cut -f1 >&2
    fi
fi

if [[ -z "${BPM_URL}" || -z "${EGT_URL}" || -z "${CSV_URL}" ]]; then
    echo "Error: Could not determine all manifest URLs." >&2
    echo "Provide --bpm-url, --egt-url, and --csv-url, or a valid --array-name." >&2
    exit 1
fi

download_file() {
    local url="$1"
    local dest_dir="$2"
    local filename
    filename=$(basename "${url}" | sed 's/?.*//')
    local dest="${dest_dir}/${filename}"

    if [[ -f "${dest}" ]]; then
        echo "  Already exists: ${dest}"
    else
        echo "  Downloading: ${url}"
        wget -q --show-progress -O "${dest}" "${url}"
        # Decompress if gzipped
        if [[ "${dest}" == *.gz ]] && file "${dest}" | grep -q gzip; then
            echo "  Decompressing: ${dest}"
            gunzip -f "${dest}"
        fi
    fi
}

echo ""
echo "Downloading BPM manifest..."
download_file "${BPM_URL}" "${OUTPUT_DIR}"

echo "Downloading EGT cluster file..."
download_file "${EGT_URL}" "${OUTPUT_DIR}"

echo "Downloading CSV manifest..."
download_file "${CSV_URL}" "${OUTPUT_DIR}"

echo ""
echo "Manifest files saved to: ${OUTPUT_DIR}"
ls -lh "${OUTPUT_DIR}"
