#!/usr/bin/env bash
#
# process_1000g.sh
#
# Download and process 1000 Genomes Project Omni2.5 IDAT data.
#
# Downloads a configurable number of samples (default: 200) from the 1000G
# HD genotype chip data, along with the required array manifests (BPM, EGT,
# CSV: HumanOmni2.5-4v1). Processes the data using the pipeline with GRCh38
# coordinates.
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
USER_BPM=""
USER_EGT=""
USER_CSV=""
USER_MANIFEST_DIR=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Download and process 1000 Genomes Omni2.5 IDAT data on GRCh38.

Required:
  --output-dir DIR       Output directory (needs ~50 GB for all samples)

Manifest options (override auto-downloaded manifests):
  --bpm FILE             Path to BPM manifest file
  --egt FILE             Path to EGT cluster file
  --csv FILE             Path to CSV manifest file
  --manifest-dir DIR     Directory containing BPM, EGT, CSV files

Options:
  --num-samples N|all    Number of samples to process (default: ${NUM_SAMPLES})
                         Use "all" for all ~2141 samples
  --threads INT          Number of threads (default: ${THREADS})
  --skip-download        Skip downloading data (use existing files)
  --skip-stage2          Skip Stage 2 reclustering
  --keep-archive         Keep the downloaded .tgz archive after extraction
  --help                 Show this help message

Examples:
  # Process 200 samples (quick test - manifests downloaded automatically)
  $(basename "$0") --output-dir ./1000g_output --num-samples 200

  # Process all ~2141 samples
  $(basename "$0") --output-dir ./1000g_output --num-samples all --threads 8

  # Use with Apptainer/Singularity on HPC
  apptainer exec --bind \$PWD docker://ghcr.io/jlanej/illumina_idat_processing:main \\
      bash /opt/scripts/process_1000g.sh --output-dir \$PWD/1000g_output

  # Provide manifest files manually (if auto-download fails)
  $(basename "$0") --output-dir ./1000g_output \\
      --bpm /path/to/HumanOmni2.5-4v1-Multi_B.bpm \\
      --egt /path/to/HumanOmni2.5-4v1-Multi_B.egt \\
      --csv /path/to/HumanOmni2.5-4v1_B.csv
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
        --num-samples)   NUM_SAMPLES="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --bpm)           USER_BPM="$2"; shift 2 ;;
        --egt)           USER_EGT="$2"; shift 2 ;;
        --csv)           USER_CSV="$2"; shift 2 ;;
        --manifest-dir)  USER_MANIFEST_DIR="$2"; shift 2 ;;
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
# Step 1: Resolve array manifests (BPM, EGT, CSV)
# ---------------------------------------------------------------
mkdir -p "${MANIFEST_DIR}"

# Handle user-provided manifest directory
if [[ -n "${USER_MANIFEST_DIR}" ]]; then
    [[ -z "${USER_BPM}" ]] && USER_BPM=$(find "${USER_MANIFEST_DIR}" -name "*.bpm" -print -quit)
    [[ -z "${USER_EGT}" ]] && USER_EGT=$(find "${USER_MANIFEST_DIR}" -name "*.egt" -print -quit)
    [[ -z "${USER_CSV}" ]] && USER_CSV=$(find "${USER_MANIFEST_DIR}" -name "*.csv" -print -quit)
fi

# Copy user-provided manifest files to the manifest directory
if [[ -n "${USER_BPM}" || -n "${USER_EGT}" || -n "${USER_CSV}" ]]; then
    echo "--- Step 1: Using provided manifest files ---"
    echo ""
    for src_var in USER_BPM USER_EGT USER_CSV; do
        src="${!src_var}"
        if [[ -n "${src}" ]]; then
            if [[ ! -f "${src}" ]]; then
                echo "Error: Manifest file not found: ${src}" >&2
                exit 1
            fi
            dest="${MANIFEST_DIR}/$(basename "${src}")"
            if [[ ! -f "${dest}" ]]; then
                cp "${src}" "${dest}"
            fi
            echo "  $(basename "${src}"): ${dest}"
        fi
    done
    echo ""
elif [[ "${SKIP_DOWNLOAD}" != "true" ]]; then
    echo "--- Step 1: Downloading array manifests ---"
    echo ""

    # The 1000G Omni2.5 uses HumanOmni2.5-4v1 manifests hosted on the 1000G FTP.
    # Note: BPM and EGT use the "Multi_B" suffix, while the CSV uses "_B".
    BPM_URL="${BPM_URL:-${FTP_BASE}/HumanOmni2.5-4v1-Multi_B.bpm}"
    EGT_URL="${EGT_URL:-${FTP_BASE}/HumanOmni2.5-4v1-Multi_B.egt}"
    CSV_URL="${CSV_URL:-${FTP_BASE}/HumanOmni2.5-4v1_B.csv}"

    for url in "${BPM_URL}" "${EGT_URL}" "${CSV_URL}"; do
        fname=$(basename "${url}")
        dest="${MANIFEST_DIR}/${fname}"
        if [[ -f "${dest}" ]]; then
            echo "  Already exists: ${dest}"
        else
            echo "  Downloading: ${url}"
            if ! wget -q --show-progress -O "${dest}" "${url}"; then
                echo "  WARNING: Download failed for ${url}" >&2
                echo "  You can download manifests manually from the 1000G FTP:" >&2
                echo "    ${FTP_BASE}/" >&2
                echo "  Then provide the files with --bpm, --egt, --csv options." >&2
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
    echo "" >&2
    echo "Download them from the 1000G FTP:" >&2
    echo "  ${FTP_BASE}/HumanOmni2.5-4v1-Multi_B.bpm" >&2
    echo "  ${FTP_BASE}/HumanOmni2.5-4v1-Multi_B.egt" >&2
    echo "  ${FTP_BASE}/HumanOmni2.5-4v1_B.csv" >&2
    echo "" >&2
    echo "Then re-run with the manifest file paths:" >&2
    echo "  $(basename "$0") --output-dir ${OUTPUT_DIR} \\" >&2
    echo "      --bpm /path/to/HumanOmni2.5-4v1-Multi_B.bpm \\" >&2
    echo "      --egt /path/to/HumanOmni2.5-4v1-Multi_B.egt \\" >&2
    echo "      --csv /path/to/HumanOmni2.5-4v1_B.csv" >&2
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

    ARCHIVE="${DOWNLOAD_DIR}/${IDAT_TGZ}"
    ARCHIVE_URL="${FTP_BASE}/${IDAT_TGZ}"

    echo "--- Step 2: Downloading 1000G IDAT files ---"
    echo "  Source: ${ARCHIVE_URL}"
    echo ""

    # Extract only IDAT files (skip GTC files to save space)
    # The archive contains files like: broad_intensities/*_Grn.idat, *_Red.idat
    if [[ "${NUM_SAMPLES}" == "all" ]]; then
        # Download full archive for all samples
        echo "  This is a large file (~20 GB) and may take a while."
        echo ""
        if [[ -f "${ARCHIVE}" ]]; then
            echo "  Archive already downloaded: ${ARCHIVE}"
        else
            wget -q --show-progress -O "${ARCHIVE}" "${ARCHIVE_URL}"
        fi

        echo ""
        echo "  Extracting IDAT files..."
        EXTRACT_START=${SECONDS}
        tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
            --wildcards '*.idat' 2>/dev/null || \
        tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat' 2>/dev/null || \
        tar xzf "${ARCHIVE}" -C "${IDAT_DIR}"
        EXTRACT_ELAPSED=$(( SECONDS - EXTRACT_START ))
        N_EXTRACTED=$(find "${IDAT_DIR}" -name '*.idat' -o -name '*.IDAT' 2>/dev/null | wc -l)
        echo "  Extraction complete: ${N_EXTRACTED} IDAT files in ${EXTRACT_ELAPSED}s"
    else
        # Stream a subset of samples directly, avoiding full archive download
        echo "  Selecting ${NUM_SAMPLES} samples..."
        echo "  Mode: subset streaming (avoids downloading full ~20 GB archive)"
        echo ""

        # Get green IDAT paths: use local archive if present, otherwise stream
        LISTING_START=${SECONDS}
        if [[ -f "${ARCHIVE}" ]]; then
            echo "  [1/3] Listing archive contents (local)..."
            echo "        Archive: ${ARCHIVE}"
            SAMPLE_LIST=$(tar tzf "${ARCHIVE}" 2>/dev/null | \
                grep -i '_Grn\.idat' | \
                head -n "${NUM_SAMPLES}")
        else
            echo "  [1/3] Listing archive contents (streaming from remote)..."
            echo "        URL: ${ARCHIVE_URL}"
            echo "        The full ~20 GB archive must be streamed to list its contents."
            echo "        This typically takes 5-15 minutes depending on connection speed."
            # Background monitor: periodically report elapsed time during listing
            (
                while true; do
                    sleep 30
                    ELAPSED=$(( SECONDS - LISTING_START ))
                    MINS=$(( ELAPSED / 60 ))
                    SECS=$(( ELAPSED % 60 ))
                    printf "        [%dm %02ds elapsed] Still scanning archive contents...\n" \
                        "${MINS}" "${SECS}"
                done
            ) &
            LISTING_MONITOR_PID=$!
            SAMPLE_LIST=$(curl -sL "${ARCHIVE_URL}" | \
                tar -tzf - 2>/dev/null | \
                grep -i '_Grn\.idat' | \
                head -n "${NUM_SAMPLES}") || true
            kill "${LISTING_MONITOR_PID}" 2>/dev/null || true
            wait "${LISTING_MONITOR_PID}" 2>/dev/null || true
        fi
        LISTING_ELAPSED=$(( SECONDS - LISTING_START ))

        N_FOUND=$(echo "${SAMPLE_LIST}" | grep -c '.' || true)
        echo "        Found ${N_FOUND} green IDAT entries in ${LISTING_ELAPSED}s"

        if [[ -z "${SAMPLE_LIST}" ]]; then
            echo ""
            echo "  [!] Could not list archive contents. Falling back to full download..."
            if [[ ! -f "${ARCHIVE}" ]]; then
                wget -q --show-progress -O "${ARCHIVE}" "${ARCHIVE_URL}"
            fi
            tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                --wildcards '*.idat' 2>/dev/null || \
            tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat' 2>/dev/null || \
            tar xzf "${ARCHIVE}" -C "${IDAT_DIR}"
        else
            # Build extraction list: both Grn and Red for selected samples
            echo ""
            echo "  [2/3] Building extraction list..."
            EXTRACT_LIST=""
            N_PAIRS=0
            while IFS= read -r grn_file; do
                red_file="${grn_file/_Grn.idat/_Red.idat}"
                red_file="${red_file/_Grn.IDAT/_Red.IDAT}"
                EXTRACT_LIST="${EXTRACT_LIST} ${grn_file} ${red_file}"
                (( N_PAIRS++ )) || true
            done <<< "${SAMPLE_LIST}"
            echo "        Will extract ${N_PAIRS} sample pairs (${N_PAIRS} Grn + ${N_PAIRS} Red = $(( N_PAIRS * 2 )) files)"

            # Extract from local archive or stream from remote
            echo ""
            EXTRACT_START=${SECONDS}
            if [[ -f "${ARCHIVE}" ]]; then
                echo "  [3/3] Extracting selected samples from local archive..."
                # shellcheck disable=SC2086
                tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                    ${EXTRACT_LIST} 2>/dev/null || \
                tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" \
                    ${EXTRACT_LIST} 2>/dev/null || {
                    echo "        Selective extraction failed. Extracting all IDATs..."
                    tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                        --wildcards '*.idat' 2>/dev/null || \
                    tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat'
                }
            else
                echo "  [3/3] Streaming extraction of ${N_PAIRS} samples from remote..."
                echo "        URL: ${ARCHIVE_URL}"
                echo "        Archive is ~20 GB — streaming to extract $(( N_PAIRS * 2 )) files."
                echo "        tar must scan through the full archive stream to find requested files."
                echo ""
                # Background monitor: periodically report elapsed time and extraction progress
                (
                    TARGET_FILES=$(( N_PAIRS * 2 ))
                    while true; do
                        sleep 15
                        N_DONE=$(find "${IDAT_DIR}" \( -name '*.idat' -o -name '*.IDAT' \) 2>/dev/null | wc -l)
                        ELAPSED=$(( SECONDS - EXTRACT_START ))
                        MINS=$(( ELAPSED / 60 ))
                        SECS=$(( ELAPSED % 60 ))
                        printf "        [%dm %02ds elapsed] %d / %d files extracted...\n" \
                            "${MINS}" "${SECS}" "${N_DONE}" "${TARGET_FILES}"
                    done
                ) &
                MONITOR_PID=$!
                # shellcheck disable=SC2086
                curl -sL "${ARCHIVE_URL}" | \
                    tar xzf - -C "${IDAT_DIR}" --strip-components=1 \
                        ${EXTRACT_LIST} 2>/dev/null || \
                curl -sL "${ARCHIVE_URL}" | \
                    tar xzf - -C "${IDAT_DIR}" \
                        ${EXTRACT_LIST} 2>/dev/null || {
                    echo "        Streaming extraction failed. Falling back to full download..."
                    wget -q --show-progress -O "${ARCHIVE}" "${ARCHIVE_URL}"
                    tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                        --wildcards '*.idat' 2>/dev/null || \
                    tar xzf "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat'
                }
                kill "${MONITOR_PID}" 2>/dev/null || true
                wait "${MONITOR_PID}" 2>/dev/null || true
            fi
            EXTRACT_ELAPSED=$(( SECONDS - EXTRACT_START ))
            N_EXTRACTED=$(find "${IDAT_DIR}" -name '*.idat' -o -name '*.IDAT' 2>/dev/null | wc -l)
            echo "        Extraction complete: ${N_EXTRACTED} files in ${EXTRACT_ELAPSED}s"
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
