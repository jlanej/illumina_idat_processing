#!/usr/bin/env bash
#
# process_1000g.sh
#
# Download and process 1000 Genomes Project Omni2.5 IDAT data.
#
# Downloads a configurable number of samples (default: 200) from the 1000G
# HD genotype chip data, along with the required array manifests (BPM, EGT,
# CSV). Processes the data using the pipeline with GRCh38 coordinates.
#
# The 1000G Omni2.5 data is available from:
#   https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/broad_intensities/
#
# Note: The original manifests are on GRCh37. The pipeline's built-in
# manifest realignment step (realign_manifest.sh) handles GRCh37→GRCh38
# coordinate remapping by aligning probe flank sequences against the
# GRCh38 reference using bwa mem, following the MoChA/gtc2vcf recommended
# approach.
#
# Usage:
#   ./process_1000g.sh --output-dir ./1000g_output --num-samples 200
#   ./process_1000g.sh --output-dir ./1000g_output --num-samples all
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 1000G FTP base URL
FTP_BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/broad_intensities"
IDAT_TGZ="Omni25_idats_gtcs_2141_samples.tgz"

# Defaults
NUM_SAMPLES=200
OUTPUT_DIR=""
THREADS=1
SKIP_DOWNLOAD="false"
SKIP_STAGE2="false"
KEEP_ARCHIVE="false"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Download and process 1000 Genomes Omni2.5 IDAT data on GRCh38.

Required:
  --output-dir DIR       Output directory (needs ~50 GB for all samples)

Options:
  --num-samples N|all    Number of samples to process (default: ${NUM_SAMPLES})
                         Use "all" for all ~2141 samples
  --threads INT          Number of threads (default: ${THREADS})
  --skip-download        Skip downloading data (use existing files)
  --skip-stage2          Skip Stage 2 reclustering
  --keep-archive         Keep the downloaded .tgz archive after extraction
  --help                 Show this help message

Examples:
  # Process 200 samples (quick test)
  $(basename "$0") --output-dir ./1000g_output --num-samples 200

  # Process all ~2141 samples
  $(basename "$0") --output-dir ./1000g_output --num-samples all --threads 8

  # Use with Apptainer/Singularity on HPC
  apptainer exec --bind \$PWD docker://ghcr.io/jlanej/illumina_idat_processing:main \\
      bash /opt/scripts/process_1000g.sh --output-dir \$PWD/1000g_output
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
        --num-samples)   NUM_SAMPLES="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD="true"; shift ;;
        --skip-stage2)   SKIP_STAGE2="true"; shift ;;
        --keep-archive)  KEEP_ARCHIVE="true"; shift ;;
        --help)          usage ;;
        *)               echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output-dir is required" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

DOWNLOAD_DIR="${OUTPUT_DIR}/downloads"
IDAT_DIR="${OUTPUT_DIR}/idats"
MANIFEST_DIR="${OUTPUT_DIR}/manifests"
PIPELINE_DIR="${OUTPUT_DIR}/pipeline_output"

echo "======================================================"
echo "  1000 Genomes Omni2.5 IDAT Processing"
echo "======================================================"
echo ""
echo "Output dir:     ${OUTPUT_DIR}"
echo "Num samples:    ${NUM_SAMPLES}"
echo "Threads:        ${THREADS}"
echo ""

# ---------------------------------------------------------------
# Step 1: Download array manifests (BPM, EGT, CSV)
# ---------------------------------------------------------------
if [[ "${SKIP_DOWNLOAD}" != "true" ]]; then
    mkdir -p "${MANIFEST_DIR}"

    echo "--- Step 1: Downloading array manifests ---"
    echo ""

    # The 1000G Omni2.5 uses InfiniumOmni2-5-8v1-3 manifests
    # NOTE: The 1000G FTP no longer hosts these proprietary Illumina files.
    # Download them from the Illumina Support portal (free registration required):
    #   https://support.illumina.com/downloads/infinium-omni2-5-8-v1-3-product-files.html
    # Then place BPM, EGT, and CSV files in the manifest directory, or re-run
    # with --skip-download and ensure the files are already present.
    BPM_URL="${BPM_URL:-${FTP_BASE}/HumanOmni2.5-8v1-3_A.bpm}"
    EGT_URL="${EGT_URL:-${FTP_BASE}/HumanOmni2.5-8v1-3_A.egt}"
    CSV_URL="${CSV_URL:-${FTP_BASE}/HumanOmni2.5-8v1-3_A.csv}"

    for url in "${BPM_URL}" "${EGT_URL}" "${CSV_URL}"; do
        fname=$(basename "${url}")
        dest="${MANIFEST_DIR}/${fname}"
        if [[ -f "${dest}" ]]; then
            echo "  Already exists: ${dest}"
        else
            echo "  Downloading: ${url}"
            if ! wget -q --show-progress -O "${dest}" "${url}"; then
                echo "  WARNING: Download failed for ${url}" >&2
                echo "  The 1000G FTP may no longer host Illumina manifest files." >&2
                echo "  Please download manually from Illumina Support:" >&2
                echo "    https://support.illumina.com/downloads/infinium-omni2-5-8-v1-3-product-files.html" >&2
                echo "  Then place the file in: ${MANIFEST_DIR}/" >&2
                rm -f "${dest}"
            fi
        fi
    done
    echo ""
fi

BPM=$(find "${MANIFEST_DIR}" -name "*.bpm" -print -quit)
EGT=$(find "${MANIFEST_DIR}" -name "*.egt" -print -quit)
CSV=$(find "${MANIFEST_DIR}" -name "*.csv" -print -quit)

if [[ -z "${BPM}" || -z "${EGT}" || -z "${CSV}" ]]; then
    echo "Error: Manifest files not found in ${MANIFEST_DIR}" >&2
    echo "Run without --skip-download to fetch them." >&2
    exit 1
fi

echo "Manifests:"
echo "  BPM: ${BPM}"
echo "  EGT: ${EGT}"
echo "  CSV: ${CSV}"
echo ""

# ---------------------------------------------------------------
# Step 2: Download and extract IDAT files
# ---------------------------------------------------------------
if [[ "${SKIP_DOWNLOAD}" != "true" ]]; then
    mkdir -p "${DOWNLOAD_DIR}" "${IDAT_DIR}"

    echo "--- Step 2: Downloading 1000G IDAT archive ---"
    echo "  Source: ${FTP_BASE}/${IDAT_TGZ}"
    echo "  This is a large file (~20 GB) and may take a while."
    echo ""

    ARCHIVE="${DOWNLOAD_DIR}/${IDAT_TGZ}"
    if [[ -f "${ARCHIVE}" ]]; then
        echo "  Archive already downloaded: ${ARCHIVE}"
    else
        wget -q --show-progress -O "${ARCHIVE}" "${FTP_BASE}/${IDAT_TGZ}"
    fi

    echo ""
    echo "  Extracting IDAT files..."

    # Extract only IDAT files (skip GTC files to save space)
    # The archive contains files like: broad_intensities/*_Grn.idat, *_Red.idat
    if [[ "${NUM_SAMPLES}" == "all" ]]; then
        tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
            --wildcards '*.idat' 2>/dev/null || \
        tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat' 2>/dev/null || \
        tar xzf "${ARCHIVE}" -C "${IDAT_DIR}"
    else
        # Extract a subset of samples
        # List all IDAT files in the archive, take N sample pairs
        echo "  Selecting ${NUM_SAMPLES} samples from archive..."

        # Get unique sample prefixes (barcode_position) from Grn IDATs
        SAMPLE_LIST=$(tar tzf "${ARCHIVE}" 2>/dev/null | \
            grep -i '_Grn\.idat' | \
            head -n "${NUM_SAMPLES}")

        if [[ -z "${SAMPLE_LIST}" ]]; then
            echo "  Could not list archive contents. Extracting all and subsetting..."
            tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                --wildcards '*.idat' 2>/dev/null || \
            tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat' 2>/dev/null || \
            tar xzf "${ARCHIVE}" -C "${IDAT_DIR}"
        else
            # Build extraction list: both Grn and Red for selected samples
            EXTRACT_LIST=""
            while IFS= read -r grn_file; do
                red_file="${grn_file/_Grn.idat/_Red.idat}"
                red_file="${red_file/_Grn.IDAT/_Red.IDAT}"
                EXTRACT_LIST="${EXTRACT_LIST} ${grn_file} ${red_file}"
            done <<< "${SAMPLE_LIST}"

            # shellcheck disable=SC2086
            tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                ${EXTRACT_LIST} 2>/dev/null || \
            tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" \
                ${EXTRACT_LIST} 2>/dev/null || {
                echo "  Selective extraction failed. Extracting all IDATs..."
                tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                    --wildcards '*.idat' 2>/dev/null || \
                tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat'
            }
        fi
    fi

    # Move any nested IDAT files to the top-level IDAT directory
    find "${IDAT_DIR}" -mindepth 2 -name '*.idat' -exec mv {} "${IDAT_DIR}/" \; 2>/dev/null || true
    find "${IDAT_DIR}" -mindepth 2 -name '*.IDAT' -exec mv {} "${IDAT_DIR}/" \; 2>/dev/null || true
    # Clean up empty subdirectories
    find "${IDAT_DIR}" -mindepth 1 -type d -empty -delete 2>/dev/null || true

    # Remove archive if not keeping it
    if [[ "${KEEP_ARCHIVE}" != "true" && -f "${ARCHIVE}" ]]; then
        echo "  Removing archive to save disk space..."
        rm -f "${ARCHIVE}"
    fi

    echo ""
fi

# Count available samples
N_GRN=$(find "${IDAT_DIR}" -name "*_Grn.idat" -o -name "*_Grn.IDAT" 2>/dev/null | wc -l)
echo "IDAT files found: ${N_GRN} samples (Green+Red pairs)"

if [[ "${N_GRN}" -eq 0 ]]; then
    echo "Error: No IDAT files found in ${IDAT_DIR}" >&2
    exit 1
fi

# If we extracted all but only want a subset, limit the IDAT directory
if [[ "${NUM_SAMPLES}" != "all" && "${N_GRN}" -gt "${NUM_SAMPLES}" ]]; then
    echo "  Subsetting to ${NUM_SAMPLES} samples..."
    SUBSET_DIR="${OUTPUT_DIR}/idats_subset"
    mkdir -p "${SUBSET_DIR}"

    find "${IDAT_DIR}" -name "*_Grn.idat" -o -name "*_Grn.IDAT" | \
        sort | head -n "${NUM_SAMPLES}" | while read -r grn; do
        red="${grn/_Grn.idat/_Red.idat}"
        red="${red/_Grn.IDAT/_Red.IDAT}"
        ln -sf "${grn}" "${SUBSET_DIR}/"
        [[ -f "${red}" ]] && ln -sf "${red}" "${SUBSET_DIR}/"
    done

    IDAT_DIR="${SUBSET_DIR}"
    N_GRN=$(find "${IDAT_DIR}" -name "*_Grn.idat" -o -name "*_Grn.IDAT" 2>/dev/null | wc -l)
    echo "  Using ${N_GRN} samples for processing"
fi

echo ""

# ---------------------------------------------------------------
# Step 3: Run the pipeline
# ---------------------------------------------------------------
echo "--- Step 3: Running processing pipeline on GRCh38 ---"
echo ""

PIPELINE_ARGS=(
    --idat-dir "${IDAT_DIR}"
    --output-dir "${PIPELINE_DIR}"
    --bpm "${BPM}"
    --egt "${EGT}"
    --csv "${CSV}"
    --genome GRCh38
    --threads "${THREADS}"
)

if [[ "${SKIP_STAGE2}" == "true" ]]; then
    PIPELINE_ARGS+=(--skip-stage2)
fi

bash "${SCRIPT_DIR}/run_pipeline.sh" "${PIPELINE_ARGS[@]}"

echo ""
echo "======================================================"
echo "  1000G Processing Complete"
echo "======================================================"
echo ""
echo "Output directory: ${PIPELINE_DIR}"
echo "Samples processed: ${N_GRN}"
echo ""
echo "To process on an HPC cluster with Apptainer:"
echo "  apptainer pull docker://ghcr.io/jlanej/illumina_idat_processing:main"
echo "  apptainer exec --bind \$PWD illumina_idat_processing_main.sif \\"
echo "      bash /opt/scripts/process_1000g.sh --output-dir \$PWD/1000g_output"
