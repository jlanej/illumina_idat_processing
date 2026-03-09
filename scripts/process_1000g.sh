#!/usr/bin/env bash
#
# process_1000g.sh
#
# Download and process 1000 Genomes Project Omni2.5 IDAT data.
#
# Downloads a configurable number of samples (default: 200) from the 1000G
# HD genotype chip data, along with the required array manifests (BPM, EGT,
# CSV: HumanOmni2.5-4v1). Processes the data using the pipeline on a selected
# genome build (default: GRCh38).
#
# The 1000G Omni2.5 data is available from:
#   https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/broad_intensities/
#
# Note: The original manifests are on GRCh37. The pipeline's built-in
# manifest realignment step (realign_manifest.sh) handles remapping probe
# coordinates onto the selected reference build using bwa mem, following the
# MoChA/gtc2vcf recommended approach.
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
GENOME="GRCh38"
SKIP_DOWNLOAD="false"
SKIP_STAGE2="false"
SKIP_FAILURES="true"   # default on for 1000G: some samples may have incomplete downloads
KEEP_ARCHIVE="false"
FORCE="false"
USER_ARCHIVE=""
USER_BPM=""
USER_EGT=""
USER_CSV=""
USER_MANIFEST_DIR=""
USER_SAMPLE_NAME_MAP=""
USER_IDAT_DIR=""
FORCE_RENAME="false"

# Default 1000G sample name map (bundled with the repo)
DEFAULT_1000G_NAME_MAP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/tests/data/sample_name_map.txt"

# Prefer pigz for faster multi-threaded gzip decompression when available.
USE_PIGZ="false"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Download and process 1000 Genomes Omni2.5 IDAT data.

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
  --archive FILE         Use a pre-downloaded Omni25_idats_gtcs_2141_samples.tgz
                         archive instead of downloading it
  --idat-dir DIR         Use a pre-extracted IDAT directory instead of
                         downloading/extracting the 1000G archive
  --sample-name-map FILE Two-column tab-delimited file mapping IDAT root names
                         to desired sample names.  Defaults to the bundled
                         1000G map (tests/data/sample_name_map.txt).
                         Use --sample-name-map "" to disable.
  --force-rename         Allow renaming even when fewer than 50% of samples in
                         the name map match the data
  --threads INT          Number of threads (default: ${THREADS})
  --genome NAME          Genome build to use for processing
                         (CHM13, GRCh37, or GRCh38; default: ${GENOME})
  --skip-download        Skip downloading data (use existing files)
  --skip-stage2          Skip Stage 2 reclustering
  --skip-failures        Continue past corrupt/truncated IDAT files (default for
                         1000G data, where a small number of samples may have
                         incomplete downloads).  Disable with --no-skip-failures.
  --no-skip-failures     Halt on IDAT read errors instead of skipping
  --keep-archive         Keep downloaded archive after extraction
  --force                Force re-run of all pipeline steps, ignoring checkpoints
  --help                 Show this help message

Examples:
  # Process 200 samples (quick test - manifests downloaded automatically)
  $(basename "$0") --output-dir ./1000g_output --num-samples 200

  # Process all ~2141 samples
  $(basename "$0") --output-dir ./1000g_output --num-samples all --threads 8

  # Use a pre-downloaded 1000G archive (will never be deleted by this script)
  $(basename "$0") --output-dir ./1000g_output \\
      --archive /data/Omni25_idats_gtcs_2141_samples.tgz --num-samples all

  # Use a pre-extracted IDAT directory
  $(basename "$0") --output-dir ./1000g_output --idat-dir /data/1000g_idats

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
        --archive)       USER_ARCHIVE="$2"; shift 2 ;;
        --idat-dir)      USER_IDAT_DIR="$2"; shift 2 ;;
        --sample-name-map) USER_SAMPLE_NAME_MAP="$2"; shift 2 ;;
        --force-rename)  FORCE_RENAME="true"; shift ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --genome)        GENOME="$2"; shift 2 ;;
        --bpm)           USER_BPM="$2"; shift 2 ;;
        --egt)           USER_EGT="$2"; shift 2 ;;
        --csv)           USER_CSV="$2"; shift 2 ;;
        --manifest-dir)  USER_MANIFEST_DIR="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD="true"; shift ;;
        --skip-stage2)   SKIP_STAGE2="true"; shift ;;
        --skip-failures) SKIP_FAILURES="true"; shift ;;
        --no-skip-failures) SKIP_FAILURES="false"; shift ;;
        --keep-archive)  KEEP_ARCHIVE="true"; shift ;;
        --force)         FORCE="true"; shift ;;
        --help)          usage ;;
        *)               echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output-dir is required" >&2
    exit 1
fi
if [[ -n "${USER_IDAT_DIR}" && ! -d "${USER_IDAT_DIR}" ]]; then
    echo "Error: --idat-dir directory not found: ${USER_IDAT_DIR}" >&2
    exit 1
fi

if command -v pigz >/dev/null 2>&1; then
    USE_PIGZ="true"
fi

# Tar helper wrappers so we can consistently enable pigz acceleration.
tar_list_file() {
    local archive="$1"
    if [[ "${USE_PIGZ}" == "true" ]]; then
        tar --use-compress-program="pigz -dc -p ${THREADS}" -tf "${archive}"
    else
        tar tzf "${archive}"
    fi
}

tar_extract_file() {
    local archive="$1"
    shift
    if [[ "${USE_PIGZ}" == "true" ]]; then
        tar --use-compress-program="pigz -dc -p ${THREADS}" -xf "${archive}" "$@"
    else
        tar xzf "${archive}" "$@"
    fi
}

tar_list_stream() {
    if [[ "${USE_PIGZ}" == "true" ]]; then
        pigz -dc -p "${THREADS}" | tar -tf -
    else
        tar tzf -
    fi
}

tar_extract_stream() {
    if [[ "${USE_PIGZ}" == "true" ]]; then
        pigz -dc -p "${THREADS}" | tar -xf - "$@"
    else
        tar xzf - "$@"
    fi
}

mkdir -p "${OUTPUT_DIR}"

DOWNLOAD_DIR="${OUTPUT_DIR}/downloads"
IDAT_DIR="${OUTPUT_DIR}/idats"
if [[ -n "${USER_IDAT_DIR}" ]]; then
    IDAT_DIR="${USER_IDAT_DIR}"
fi
MANIFEST_DIR="${OUTPUT_DIR}/manifests"
PIPELINE_DIR="${OUTPUT_DIR}/pipeline_output"

echo "======================================================"
echo "  1000 Genomes Omni2.5 IDAT Processing"
echo "======================================================"
echo ""
echo "Output dir:     ${OUTPUT_DIR}"
echo "Num samples:    ${NUM_SAMPLES}"
echo "Threads:        ${THREADS}"
echo "Genome:         ${GENOME}"
if [[ "${USE_PIGZ}" == "true" ]]; then
    echo "Gzip backend:   pigz (parallel)"
else
    echo "Gzip backend:   gzip (tar default)"
fi
if [[ -n "${USER_ARCHIVE}" ]]; then
    echo "Archive file:   ${USER_ARCHIVE}"
fi
if [[ -n "${USER_IDAT_DIR}" ]]; then
    echo "IDAT dir:       ${USER_IDAT_DIR}"
fi
echo ""

# ---------------------------------------------------------------
# Step 1: Resolve array manifests (BPM, EGT, CSV)
# ---------------------------------------------------------------
mkdir -p "${MANIFEST_DIR}"

# Handle user-provided manifest directory
if [[ -n "${USER_MANIFEST_DIR}" ]]; then
    [[ -z "${USER_BPM}" ]] && USER_BPM=$(find "${USER_MANIFEST_DIR}" -name "*.bpm" -print -quit 2>/dev/null || true)
    [[ -z "${USER_EGT}" ]] && USER_EGT=$(find "${USER_MANIFEST_DIR}" -name "*.egt" -print -quit 2>/dev/null || true)
    [[ -z "${USER_CSV}" ]] && USER_CSV=$(find "${USER_MANIFEST_DIR}" -name "*.csv" -print -quit 2>/dev/null || true)
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

BPM=$(find "${MANIFEST_DIR}" -name "*.bpm" -print -quit 2>/dev/null || true)
EGT=$(find "${MANIFEST_DIR}" -name "*.egt" -print -quit 2>/dev/null || true)
CSV=$(find "${MANIFEST_DIR}" -name "*.csv" -print -quit 2>/dev/null || true)

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
if [[ -n "${USER_IDAT_DIR}" ]]; then
    echo "--- Step 2: Using pre-extracted IDAT directory ---"
    echo "  IDAT dir: ${IDAT_DIR}"
    echo "  Skipping 1000G archive download/extraction."
    echo ""
elif [[ "${SKIP_DOWNLOAD}" != "true" || -n "${USER_ARCHIVE}" ]]; then
    mkdir -p "${DOWNLOAD_DIR}" "${IDAT_DIR}"

    ARCHIVE_URL="${FTP_BASE}/${IDAT_TGZ}"
    ARCHIVE_FROM_USER="false"
    REMOTE_ARCHIVE_ALLOWED="true"
    ARCHIVE="${DOWNLOAD_DIR}/${IDAT_TGZ}"

    if [[ -n "${USER_ARCHIVE}" ]]; then
        if [[ ! -f "${USER_ARCHIVE}" ]]; then
            echo "Error: --archive file not found: ${USER_ARCHIVE}" >&2
            exit 1
        fi
        ARCHIVE="${USER_ARCHIVE}"
        ARCHIVE_FROM_USER="true"
        REMOTE_ARCHIVE_ALLOWED="false"
    fi

    # Check if the correct number of IDAT files are already present
    EXISTING_GRN=$(find "${IDAT_DIR}" -name "*_Grn.idat" -o -name "*_Grn.IDAT" 2>/dev/null | wc -l)
    if [[ "${NUM_SAMPLES}" != "all" && "${EXISTING_GRN}" -ge "${NUM_SAMPLES}" ]]; then
        echo "--- Step 2: IDAT files already present ---"
        echo "  Found ${EXISTING_GRN} samples in ${IDAT_DIR} (need ${NUM_SAMPLES})"
        echo "  Skipping download."
        echo ""
    else

    echo "--- Step 2: Downloading 1000G IDAT files ---"
    if [[ "${ARCHIVE_FROM_USER}" == "true" ]]; then
        echo "  Source: local archive ${ARCHIVE}"
    else
        echo "  Source: ${ARCHIVE_URL}"
    fi
    echo ""

    # Extract only IDAT files (skip GTC files to save space)
    # The archive contains files like: broad_intensities/*_Grn.idat, *_Red.idat
    if [[ "${NUM_SAMPLES}" == "all" ]]; then
        # Download full archive for all samples
        echo "  This is a large file (~20 GB) and may take a while."
        echo ""
        if [[ -f "${ARCHIVE}" ]]; then
            echo "  Archive already available: ${ARCHIVE}"
        elif [[ "${REMOTE_ARCHIVE_ALLOWED}" == "true" ]]; then
            wget -q --show-progress -O "${ARCHIVE}" "${ARCHIVE_URL}"
        else
            echo "Error: Local archive not available and remote download is disabled." >&2
            exit 1
        fi

        echo ""
        echo "  Extracting IDAT files..."
        EXTRACT_START=${SECONDS}
        tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
            --wildcards '*.idat' 2>/dev/null || \
        tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat' 2>/dev/null || \
        tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}"
        EXTRACT_ELAPSED=$(( SECONDS - EXTRACT_START ))
        N_EXTRACTED=$(find "${IDAT_DIR}" -name '*.idat' -o -name '*.IDAT' 2>/dev/null | wc -l)
        echo "  Extraction complete: ${N_EXTRACTED} IDAT files in ${EXTRACT_ELAPSED}s"
    else
        # Stream a subset of samples directly, avoiding full archive download
        echo "  Selecting ${NUM_SAMPLES} samples..."
        if [[ "${REMOTE_ARCHIVE_ALLOWED}" == "true" ]]; then
            echo "  Mode: subset streaming (avoids downloading full ~20 GB archive)"
        else
            echo "  Mode: subset extraction from local archive"
        fi
        echo ""

        # Get green IDAT paths: use local archive if present, otherwise stream
        LISTING_START=${SECONDS}
        if [[ -f "${ARCHIVE}" ]]; then
            echo "  [1/3] Listing archive contents (local)..."
            echo "        Archive: ${ARCHIVE}"
            SAMPLE_LIST=$(tar_list_file "${ARCHIVE}" 2>/dev/null | \
                grep -i '_Grn\.idat' | \
                head -n "${NUM_SAMPLES}")
        else
            if [[ "${REMOTE_ARCHIVE_ALLOWED}" != "true" ]]; then
                echo "Error: --archive was provided but not found: ${ARCHIVE}" >&2
                exit 1
            fi
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
                tar_list_stream 2>/dev/null | \
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
            if [[ "${REMOTE_ARCHIVE_ALLOWED}" == "true" ]]; then
                echo "  [!] Could not list archive contents. Falling back to full download..."
                if [[ ! -f "${ARCHIVE}" ]]; then
                    wget -q --show-progress -O "${ARCHIVE}" "${ARCHIVE_URL}"
                fi
                tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                    --wildcards '*.idat' 2>/dev/null || \
                tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat' 2>/dev/null || \
                tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}"
            else
                echo "Error: Could not list IDAT entries from local archive: ${ARCHIVE}" >&2
                exit 1
            fi
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
                tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                    ${EXTRACT_LIST} 2>/dev/null || \
                tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}" \
                    ${EXTRACT_LIST} 2>/dev/null || {
                    echo "        Selective extraction failed. Extracting all IDATs..."
                    tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                        --wildcards '*.idat' 2>/dev/null || \
                    tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat'
                }
            else
                echo "  [3/3] Streaming extraction of ${N_PAIRS} samples from remote..."
                echo "        URL: ${ARCHIVE_URL}"
                echo "        Archive is ~20 GB — streaming to extract $(( N_PAIRS * 2 )) files."
                echo "        This may take several minutes as the full archive is streamed."
                echo ""
                BASELINE_FILES=$(find "${IDAT_DIR}" \( -name '*.idat' -o -name '*.IDAT' \) 2>/dev/null | wc -l)
                TARGET_FILES=$(( N_PAIRS * 2 ))

                # Monitor a background extraction pipeline, reporting progress
                # and terminating the stream early once all files are extracted.
                _monitor_extraction() {
                    local extract_pid="$1"
                    while kill -0 "${extract_pid}" 2>/dev/null; do
                        sleep 15
                        local n_total n_done elapsed mins secs
                        n_total=$(find "${IDAT_DIR}" \( -name '*.idat' -o -name '*.IDAT' \) 2>/dev/null | wc -l)
                        n_done=$(( n_total - BASELINE_FILES ))
                        elapsed=$(( SECONDS - EXTRACT_START ))
                        mins=$(( elapsed / 60 ))
                        secs=$(( elapsed % 60 ))
                        printf "        [%dm %02ds elapsed] %d / %d files extracted...\n" \
                            "${mins}" "${secs}" "${n_done}" "${TARGET_FILES}"
                        if [[ "${n_done}" -ge "${TARGET_FILES}" ]]; then
                            printf "        All %d files extracted. Stopping archive stream.\n" "${TARGET_FILES}"
                            kill "${extract_pid}" 2>/dev/null || true
                            break
                        fi
                    done
                    wait "${extract_pid}" 2>/dev/null || true
                }

                # Attempt 1: stream with --strip-components=1
                # shellcheck disable=SC2086
                curl -sL "${ARCHIVE_URL}" | \
                    tar_extract_stream -C "${IDAT_DIR}" --strip-components=1 \
                        ${EXTRACT_LIST} 2>/dev/null &
                EXTRACT_PID=$!
                _monitor_extraction "${EXTRACT_PID}"

                # Check result; fall back to alternative extraction strategies
                N_DONE=$(( $(find "${IDAT_DIR}" \( -name '*.idat' -o -name '*.IDAT' \) 2>/dev/null | wc -l) - BASELINE_FILES ))

                if [[ "${N_DONE}" -eq 0 ]]; then
                    # Attempt 2: stream without --strip-components
                    echo "        Retrying without --strip-components..."
                    # shellcheck disable=SC2086
                    curl -sL "${ARCHIVE_URL}" | \
                        tar_extract_stream -C "${IDAT_DIR}" \
                            ${EXTRACT_LIST} 2>/dev/null &
                    EXTRACT_PID=$!
                    _monitor_extraction "${EXTRACT_PID}"
                    N_DONE=$(( $(find "${IDAT_DIR}" \( -name '*.idat' -o -name '*.IDAT' \) 2>/dev/null | wc -l) - BASELINE_FILES ))
                fi

                if [[ "${N_DONE}" -eq 0 ]]; then
                    # Attempt 3: full download fallback
                    echo "        Streaming extraction failed. Falling back to full download..."
                    wget -q --show-progress -O "${ARCHIVE}" "${ARCHIVE_URL}"
                    tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}" --strip-components=1 \
                        --wildcards '*.idat' 2>/dev/null || \
                    tar_extract_file "${ARCHIVE}" -C "${IDAT_DIR}" --wildcards '*.idat'
                fi
            fi
            EXTRACT_ELAPSED=$(( SECONDS - EXTRACT_START ))
            N_EXTRACTED=$(find "${IDAT_DIR}" -name '*.idat' -o -name '*.IDAT' 2>/dev/null | wc -l)
            echo "        Extraction complete: ${N_EXTRACTED} files in ${EXTRACT_ELAPSED}s"
        fi
    fi

    fi  # end of "already have enough IDATs" check

    # Move any nested IDAT files to the top-level IDAT directory
    find "${IDAT_DIR}" -mindepth 2 -name '*.idat' -exec mv {} "${IDAT_DIR}/" \; 2>/dev/null || true
    find "${IDAT_DIR}" -mindepth 2 -name '*.IDAT' -exec mv {} "${IDAT_DIR}/" \; 2>/dev/null || true
    # Clean up empty subdirectories
    find "${IDAT_DIR}" -mindepth 1 -type d -empty -delete 2>/dev/null || true

    # Remove archive if not keeping it (never remove a user-provided archive).
    if [[ "${KEEP_ARCHIVE}" != "true" && "${ARCHIVE_FROM_USER}" != "true" && -f "${ARCHIVE}" ]]; then
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

# If more samples are present than requested, process them all
if [[ "${NUM_SAMPLES}" != "all" && "${N_GRN}" -gt "${NUM_SAMPLES}" ]]; then
    echo "  Found ${N_GRN} samples (more than requested ${NUM_SAMPLES}). Processing all."
fi

echo ""

# ---------------------------------------------------------------
# Step 3: Run the pipeline
# ---------------------------------------------------------------
echo "--- Step 3: Running processing pipeline on ${GENOME} ---"
echo ""

# Resolve sample name map: use user-provided, or default 1000G map
SAMPLE_NAME_MAP=""
if [[ -n "${USER_SAMPLE_NAME_MAP}" ]]; then
    SAMPLE_NAME_MAP="${USER_SAMPLE_NAME_MAP}"
elif [[ -f "${DEFAULT_1000G_NAME_MAP}" ]]; then
    SAMPLE_NAME_MAP="${DEFAULT_1000G_NAME_MAP}"
    echo "  Using default 1000G sample name map: ${SAMPLE_NAME_MAP}"
else
    echo "  WARNING: Default 1000G sample name map not found at:"
    echo "    ${DEFAULT_1000G_NAME_MAP}"
    echo "  Samples will retain Sentrix barcode IDs instead of being renamed"
    echo "  to 1000G sample names (e.g., NA12878)."
    echo "  To enable renaming, provide --sample-name-map /path/to/map.txt"
    echo ""
fi

PIPELINE_ARGS=(
    --idat-dir "${IDAT_DIR}"
    --output-dir "${PIPELINE_DIR}"
    --bpm "${BPM}"
    --egt "${EGT}"
    --csv "${CSV}"
    --genome "${GENOME}"
    --threads "${THREADS}"
)

if [[ -n "${SAMPLE_NAME_MAP}" ]]; then
    PIPELINE_ARGS+=(--sample-name-map "${SAMPLE_NAME_MAP}")
fi

if [[ "${FORCE_RENAME}" == "true" ]]; then
    PIPELINE_ARGS+=(--force-rename)
fi

if [[ "${SKIP_STAGE2}" == "true" ]]; then
    PIPELINE_ARGS+=(--skip-stage2)
fi

if [[ "${SKIP_FAILURES}" == "true" ]]; then
    PIPELINE_ARGS+=(--skip-failures)
fi

if [[ "${FORCE}" == "true" ]]; then
    PIPELINE_ARGS+=(--force)
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
