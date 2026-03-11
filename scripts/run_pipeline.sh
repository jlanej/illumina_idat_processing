#!/usr/bin/env bash
#
# run_pipeline.sh
#
# Main entry point for the Illumina IDAT processing pipeline.
# Processes raw Illumina array IDAT data in two stages:
#
#   Stage 1: Initial genotyping - IDAT → GTC → VCF with QC metrics
#   Stage 2: Reclustering - filter high-quality samples, recluster, reprocess
#
# Designed to require minimal arguments and automate manifest/resource
# acquisition where possible.
#
# Inspired by mocha (https://github.com/freeseek/mocha) and mochawdl
# (https://github.com/freeseek/mochawdl).
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults
GENOME="CHM13"
THREADS=1
MIN_CALL_RATE=0.97
MAX_LRR_SD=0.35
MIN_ANCESTRY_SAMPLES=100
SKIP_STAGE2="false"
SKIP_DOWNLOAD="false"
SKIP_FAILURES="false"
SKIP_PEDDY="false"
FORCE="false"
FORCE_RENAME="false"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Process raw Illumina IDAT files to genotyped VCF with two-stage QC.

Required (minimal):
  --idat-dir DIR         Directory containing IDAT files
  --output-dir DIR       Output directory

Manifest options (one of the following):
  --array-name NAME      Known array name (auto-downloads manifests)
  --manifest-dir DIR     Directory containing BPM, EGT, and CSV files
  --bpm FILE             BPM manifest file (with --egt and --csv)
  --egt FILE             EGT cluster file (with --bpm and --csv)
  --csv FILE             CSV manifest file (with --bpm and --egt)

Reference options:
  --ref-dir DIR          Directory with reference genome (auto-downloads if missing)
  --ref-fasta FILE       Reference FASTA file (overrides --ref-dir)
  --genome NAME          Genome build: CHM13, GRCh37, or GRCh38 (default: ${GENOME})

Processing options:
  --threads INT          Number of threads (default: ${THREADS})
  --min-call-rate FLOAT  Min call rate for Stage 2 HQ samples (default: ${MIN_CALL_RATE})
  --max-lrr-sd FLOAT     Max LRR SD for Stage 2 HQ samples (default: ${MAX_LRR_SD})
  --min-ancestry-samples INT
                         Minimum samples per ancestry group for ancestry-
                         stratified QC (variant QC and PCA within each
                         ancestry). Uses peddy predictions. (default: ${MIN_ANCESTRY_SAMPLES})
  --sample-name-map FILE Two-column tab-delimited file mapping IDAT root names
                         (column 1) to desired sample names (column 2).
                         Renames samples in GTC files, VCFs, and QC reports.
  --force-rename         Allow renaming even when fewer than 50% of samples in
                         the name map match the data (default: halt)
  --ped-file FILE        Pedigree file (.ped/.fam) for peddy relationship
                         validation. If not provided, a minimal PED is
                         generated (unrelated singletons) so peddy can still
                         predict sex, ancestry, and discover relationships.
  --skip-peddy           Skip peddy pedigree/sex/ancestry QC step
  --skip-stage2          Skip Stage 2 (reclustering)
  --skip-download        Do not auto-download manifests or reference
  --skip-failures        Continue past corrupt/truncated IDAT files instead of
                         halting.  Failed samples are logged and excluded from
                         downstream outputs.
  --force                Force re-run of all steps, ignoring checkpoints

  --help                 Show this help message

Examples:
  # Minimal: auto-download manifests and reference for a known array
  $(basename "$0") --idat-dir /data/idats --array-name GSA-24v3-0_A1 --output-dir /data/output

  # With pre-downloaded manifests and reference
  $(basename "$0") --idat-dir /data/idats --manifest-dir /data/manifests \\
      --ref-fasta /data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \\
      --output-dir /data/output

  # Provide manifest files directly
  $(basename "$0") --idat-dir /data/idats --bpm array.bpm --egt array.egt \\
      --csv array.csv --ref-dir /data/ref --output-dir /data/output
EOF
    exit 0
}

IDAT_DIR=""
OUTPUT_DIR=""
ARRAY_NAME=""
MANIFEST_DIR=""
BPM=""
EGT=""
CSV=""
REF_DIR=""
REF_FASTA=""
SAMPLE_NAME_MAP=""
PED_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --idat-dir)       IDAT_DIR="$2"; shift 2 ;;
        --output-dir)     OUTPUT_DIR="$2"; shift 2 ;;
        --array-name)     ARRAY_NAME="$2"; shift 2 ;;
        --manifest-dir)   MANIFEST_DIR="$2"; shift 2 ;;
        --bpm)            BPM="$2"; shift 2 ;;
        --egt)            EGT="$2"; shift 2 ;;
        --csv)            CSV="$2"; shift 2 ;;
        --ref-dir)        REF_DIR="$2"; shift 2 ;;
        --ref-fasta)      REF_FASTA="$2"; shift 2 ;;
        --genome)         GENOME="$2"; shift 2 ;;
        --threads)        THREADS="$2"; shift 2 ;;
        --min-call-rate)  MIN_CALL_RATE="$2"; shift 2 ;;
        --max-lrr-sd)     MAX_LRR_SD="$2"; shift 2 ;;
        --min-ancestry-samples) MIN_ANCESTRY_SAMPLES="$2"; shift 2 ;;
        --sample-name-map) SAMPLE_NAME_MAP="$2"; shift 2 ;;
        --ped-file)       PED_FILE="$2"; shift 2 ;;
        --force-rename)   FORCE_RENAME="true"; shift ;;
        --skip-stage2)    SKIP_STAGE2="true"; shift ;;
        --skip-download)  SKIP_DOWNLOAD="true"; shift ;;
        --skip-failures)  SKIP_FAILURES="true"; shift ;;
        --skip-peddy)     SKIP_PEDDY="true"; shift ;;
        --force)          FORCE="true"; shift ;;
        --help)           usage ;;
        *)                echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${IDAT_DIR}" ]]; then
    echo "Error: --idat-dir is required" >&2
    exit 1
fi
if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output-dir is required" >&2
    exit 1
fi
if [[ ! -d "${IDAT_DIR}" ]]; then
    echo "Error: IDAT directory not found: ${IDAT_DIR}" >&2
    exit 1
fi

# Check bcftools is available
if ! command -v bcftools &>/dev/null; then
    echo "Error: bcftools not found. Run scripts/install_dependencies.sh first." >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

any_input_newer_than() {
    local output="$1"
    shift
    local input
    for input in "$@"; do
        if [[ -n "${input}" && -e "${input}" && "${input}" -nt "${output}" ]]; then
            return 0
        fi
    done
    return 1
}

echo "======================================================"
echo "  Illumina IDAT Processing Pipeline"
echo "======================================================"
echo ""

# ---------------------------------------------------------------
# Resolve manifest files
# ---------------------------------------------------------------
if [[ -n "${MANIFEST_DIR}" ]]; then
    # Auto-detect BPM, EGT, CSV from manifest directory
    [[ -z "${BPM}" ]] && BPM=$(find "${MANIFEST_DIR}" -name "*.bpm" -print -quit 2>/dev/null || true)
    [[ -z "${EGT}" ]] && EGT=$(find "${MANIFEST_DIR}" -name "*.egt" -print -quit 2>/dev/null || true)
    [[ -z "${CSV}" ]] && CSV=$(find "${MANIFEST_DIR}" -name "*.csv" -print -quit 2>/dev/null || true)
fi

if [[ -z "${BPM}" || -z "${EGT}" || -z "${CSV}" ]]; then
    if [[ "${SKIP_DOWNLOAD}" == "true" ]]; then
        echo "Error: Manifest files not found and --skip-download is set." >&2
        echo "Provide --bpm, --egt, --csv or --manifest-dir." >&2
        exit 1
    fi

    MANIFEST_DIR="${OUTPUT_DIR}/manifests"
    echo "Downloading manifest files..."
    DOWNLOAD_ARGS=(--output-dir "${MANIFEST_DIR}")

    if [[ -n "${ARRAY_NAME}" ]]; then
        DOWNLOAD_ARGS+=(--array-name "${ARRAY_NAME}")
    else
        echo "  Attempting to auto-detect array type from IDAT files..."
        DOWNLOAD_ARGS+=(--detect-from "${IDAT_DIR}")
    fi

    bash "${SCRIPT_DIR}/download_manifests.sh" "${DOWNLOAD_ARGS[@]}"
    echo ""

    [[ -z "${BPM}" ]] && BPM=$(find "${MANIFEST_DIR}" -name "*.bpm" -print -quit 2>/dev/null || true)
    [[ -z "${EGT}" ]] && EGT=$(find "${MANIFEST_DIR}" -name "*.egt" -print -quit 2>/dev/null || true)
    [[ -z "${CSV}" ]] && CSV=$(find "${MANIFEST_DIR}" -name "*.csv" -print -quit 2>/dev/null || true)

    if [[ -z "${BPM}" || -z "${EGT}" || -z "${CSV}" ]]; then
        echo "Error: Could not find all manifest files after download." >&2
        echo "Please provide --bpm, --egt, and --csv directly." >&2
        exit 1
    fi
fi

# ---------------------------------------------------------------
# Resolve reference genome
# ---------------------------------------------------------------
if [[ -z "${REF_FASTA}" ]]; then
    if [[ -z "${REF_DIR}" ]]; then
        REF_DIR="${OUTPUT_DIR}/reference"
    fi

    # Try to find existing reference
    if [[ "${GENOME}" == "CHM13" ]]; then
        REF_FASTA="${REF_DIR}/chm13v2.0.fa"
    elif [[ "${GENOME}" == "GRCh38" ]]; then
        REF_FASTA="${REF_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    else
        REF_FASTA="${REF_DIR}/human_g1k_v37.fasta"
    fi

    if [[ ! -f "${REF_FASTA}" ]]; then
        if [[ "${SKIP_DOWNLOAD}" == "true" ]]; then
            echo "Error: Reference genome not found and --skip-download is set." >&2
            exit 1
        fi
        echo "Downloading reference genome (${GENOME})..."
        bash "${SCRIPT_DIR}/download_reference.sh" \
            --genome "${GENOME}" \
            --output-dir "${REF_DIR}"
        echo ""
    fi
fi

if [[ ! -f "${REF_FASTA}" ]]; then
    echo "Error: Reference FASTA not found: ${REF_FASTA}" >&2
    exit 1
fi

# ---------------------------------------------------------------
# Print configuration
# ---------------------------------------------------------------
echo "Configuration:"
echo "  IDAT dir:       ${IDAT_DIR}"
echo "  Output dir:     ${OUTPUT_DIR}"
echo "  BPM:            ${BPM}"
echo "  EGT:            ${EGT}"
echo "  CSV:            ${CSV}"
echo "  Reference:      ${REF_FASTA}"
echo "  Genome:         ${GENOME}"
echo "  Threads:        ${THREADS}"
echo "  Min call rate:  ${MIN_CALL_RATE}"
echo "  Max LRR SD:     ${MAX_LRR_SD}"
if [[ -n "${SAMPLE_NAME_MAP}" ]]; then
    echo "  Name map:       ${SAMPLE_NAME_MAP}"
fi
if [[ -n "${PED_FILE}" ]]; then
    echo "  Pedigree file:  ${PED_FILE}"
fi
echo ""

# ---------------------------------------------------------------
# Realign manifest probe sequences against reference genome
# ---------------------------------------------------------------
# Following MoChA/gtc2vcf best practice, always realign the CSV manifest
# flank sequences against the reference. This validates probe positions
# and enables correct coordinate mapping regardless of the original build.
REALIGN_DIR="${OUTPUT_DIR}/realigned_manifests"

if command -v bwa &>/dev/null && [[ -f "${REF_FASTA}.bwt" ]]; then
    # Check if realignment was already completed
    EXISTING_REALIGNED=$(find "${REALIGN_DIR}" -name "*.realigned.csv" -print -quit 2>/dev/null || true)
    if [[ -n "${EXISTING_REALIGNED}" && -f "${EXISTING_REALIGNED}" && "${FORCE}" != "true" ]]; then
        echo "======================================================"
        echo "  Manifest Probe Realignment (cached)"
        echo "======================================================"
        echo ""
        echo "Using existing realigned CSV: ${EXISTING_REALIGNED}"
        CSV="${EXISTING_REALIGNED}"
    else
        echo "======================================================"
        echo "  Realigning Manifest Probes Against Reference"
        echo "======================================================"
        echo ""
        echo "This step validates probe coordinates by aligning flank sequences"
        echo "against the reference genome. For large arrays (>1M probes) this"
        echo "may take 15-45 minutes."
        echo ""

        REALIGN_START=${SECONDS}
        bash "${SCRIPT_DIR}/realign_manifest.sh" \
            --csv "${CSV}" \
            --ref-fasta "${REF_FASTA}" \
            --output-dir "${REALIGN_DIR}" \
            --threads "${THREADS}"
        REALIGN_ELAPSED=$(( SECONDS - REALIGN_START ))
        echo "Manifest realignment completed in $(( REALIGN_ELAPSED / 60 ))m $(( REALIGN_ELAPSED % 60 ))s"

        # Use the realigned CSV for all downstream steps
        REALIGNED_CSV=$(find "${REALIGN_DIR}" -name "*.realigned.csv" -print -quit 2>/dev/null || true)
        if [[ -n "${REALIGNED_CSV}" && -f "${REALIGNED_CSV}" ]]; then
            echo "Using realigned CSV: ${REALIGNED_CSV}"
            CSV="${REALIGNED_CSV}"
        fi
    fi

    # Realignment quality gate: warn if too many probes are unmapped
    REALIGN_SUMMARY=$(find "${REALIGN_DIR}" -name "realign_summary.txt" -print -quit 2>/dev/null || true)
    if [[ -n "${REALIGN_SUMMARY}" && -f "${REALIGN_SUMMARY}" ]]; then
        MAPPED_PCT=$(awk '/Mapped:/ && /%/ { gsub(/[^0-9.]/, "", $NF); print $NF; exit }' "${REALIGN_SUMMARY}")
        if [[ -n "${MAPPED_PCT}" ]]; then
            MAPPED_OK=$(awk -v p="${MAPPED_PCT}" 'BEGIN { print (p+0 >= 90) ? "yes" : "no" }')
            if [[ "${MAPPED_OK}" == "no" ]]; then
                echo ""
                echo "WARNING: Only ${MAPPED_PCT}% of manifest probes mapped to the reference."
                echo "  This may indicate a build mismatch or poor manifest quality."
                echo "  Expect reduced call rates for unmapped probes."
                echo "  Review: ${REALIGN_SUMMARY}"
                echo ""
            else
                echo "Realignment quality: ${MAPPED_PCT}% of probes mapped (OK)"
            fi
        fi
    fi
    echo ""
else
    echo "Note: Skipping manifest realignment (bwa or BWA index not available)."
    echo "  Probe positions will use coordinates from the original CSV manifest."
    echo "  For best results, install bwa and create an index:"
    echo "    bwa index ${REF_FASTA}"
    echo ""
fi

# ---------------------------------------------------------------
# Stage 1: Initial genotyping
# ---------------------------------------------------------------
STAGE1_DIR="${OUTPUT_DIR}/stage1"

echo "======================================================"
echo "  Running Stage 1: Initial Genotyping"
echo "======================================================"
echo ""

STAGE1_ARGS=(
    --idat-dir "${IDAT_DIR}"
    --bpm "${BPM}"
    --egt "${EGT}"
    --csv "${CSV}"
    --ref-fasta "${REF_FASTA}"
    --output-dir "${STAGE1_DIR}"
    --threads "${THREADS}"
)
if [[ -n "${SAMPLE_NAME_MAP}" ]]; then
    STAGE1_ARGS+=(--sample-name-map "${SAMPLE_NAME_MAP}")
fi
if [[ "${FORCE_RENAME}" == "true" ]]; then
    STAGE1_ARGS+=(--force-rename)
fi
if [[ "${FORCE}" == "true" ]]; then
    STAGE1_ARGS+=(--force)
fi
if [[ "${SKIP_FAILURES}" == "true" ]]; then
    STAGE1_ARGS+=(--skip-failures)
fi

bash "${SCRIPT_DIR}/stage1_initial_genotyping.sh" "${STAGE1_ARGS[@]}"

echo ""

# ---------------------------------------------------------------
# Stage 2: Recluster and reprocess
# ---------------------------------------------------------------
if [[ "${SKIP_STAGE2}" == "true" ]]; then
    echo "Skipping Stage 2 (--skip-stage2 specified)"
    FINAL_VCF="${STAGE1_DIR}/vcf/stage1_initial.bcf"
    FINAL_QC="${STAGE1_DIR}/qc/stage1_sample_qc.tsv"
else
    STAGE2_DIR="${OUTPUT_DIR}/stage2"

    echo "======================================================"
    echo "  Running Stage 2: Recluster and Reprocess"
    echo "======================================================"
    echo ""

    STAGE2_ARGS=(
        --idat-dir "${IDAT_DIR}"
        --bpm "${BPM}"
        --egt "${EGT}"
        --csv "${CSV}"
        --ref-fasta "${REF_FASTA}"
        --stage1-dir "${STAGE1_DIR}"
        --output-dir "${STAGE2_DIR}"
        --min-call-rate "${MIN_CALL_RATE}"
        --max-lrr-sd "${MAX_LRR_SD}"
        --threads "${THREADS}"
    )
    if [[ -n "${SAMPLE_NAME_MAP}" ]]; then
        STAGE2_ARGS+=(--sample-name-map "${SAMPLE_NAME_MAP}")
    fi
    if [[ "${FORCE}" == "true" ]]; then
        STAGE2_ARGS+=(--force)
    fi
    if [[ "${SKIP_FAILURES}" == "true" ]]; then
        STAGE2_ARGS+=(--skip-failures)
    fi

    bash "${SCRIPT_DIR}/stage2_recluster.sh" "${STAGE2_ARGS[@]}"

    FINAL_VCF="${STAGE2_DIR}/vcf/stage2_reclustered.bcf"
    FINAL_QC="${STAGE2_DIR}/qc/stage2_sample_qc.tsv"
    COMPARISON_DIR="${STAGE2_DIR}/qc/comparison"
fi

echo ""
echo "======================================================"
echo "  Pipeline Complete (Stages 1-2)"
echo "======================================================"
echo ""
echo "Final outputs:"
echo "  VCF:            ${FINAL_VCF}"
echo "  QC metrics:     ${FINAL_QC}"
if [[ -n "${COMPARISON_DIR:-}" && -d "${COMPARISON_DIR:-}" ]]; then
    echo "  QC comparison:  ${COMPARISON_DIR}/"
fi
echo ""

# ---------------------------------------------------------------
# Sex check plot
# ---------------------------------------------------------------
echo "======================================================"
echo "  Generating Sex Check Plot"
echo "======================================================"
echo ""

SEX_CHECK_DIR="${OUTPUT_DIR}/sex_check"
SEX_CHECK_PNG="${SEX_CHECK_DIR}/sex_check_chrXY_lrr.png"
SEX_CHECK_TSV="${SEX_CHECK_DIR}/sex_check_chrXY_lrr.tsv"
if command -v python3 &>/dev/null && [[ -f "${FINAL_VCF}" && -f "${FINAL_QC}" ]]; then
    if [[ "${FORCE}" != "true" && -s "${SEX_CHECK_PNG}" && -s "${SEX_CHECK_TSV}" ]]; then
        echo "Sex check outputs already exist. Skipping (use --force to regenerate)."
    else
        python3 "${SCRIPT_DIR}/plot_sex_check.py" \
            --vcf "${FINAL_VCF}" \
            --sample-qc "${FINAL_QC}" \
            --output-dir "${SEX_CHECK_DIR}" \
            --threads "${THREADS}" 2>&1 || true
    fi
else
    echo "Note: Skipping sex check plot (python3, VCF, or QC file not available)."
fi
echo ""

# ---------------------------------------------------------------
# Peddy: pedigree/sex/ancestry QC
# ---------------------------------------------------------------
PEDDY_DIR="${OUTPUT_DIR}/peddy"

if [[ "${SKIP_PEDDY}" == "true" ]]; then
    echo "Skipping peddy QC (--skip-peddy specified)"
else
    echo "======================================================"
    echo "  Running Peddy (Pedigree/Sex/Ancestry QC)"
    echo "======================================================"
    echo ""

    if [[ -f "${FINAL_VCF}" ]]; then
        echo "Using pipeline final stage VCF as peddy base input: ${FINAL_VCF}"
        PEDDY_FORCE="${FORCE}"
        PEDDY_HET_CHECK="${PEDDY_DIR}/peddy.het_check.csv"
        if [[ "${PEDDY_FORCE}" != "true" && -s "${PEDDY_HET_CHECK}" ]]; then
            if any_input_newer_than "${PEDDY_HET_CHECK}" "${FINAL_VCF}" "${FINAL_QC}" "${PED_FILE}"; then
                echo "Peddy outputs exist but inputs are newer. Regenerating."
                PEDDY_FORCE="true"
            else
                echo "Peddy outputs already exist and are up to date. Skipping (use --force to regenerate)."
            fi
        fi

        if [[ "${PEDDY_FORCE}" == "true" || ! -s "${PEDDY_HET_CHECK}" ]]; then
            PEDDY_ARGS=(
                --vcf "${FINAL_VCF}"
                --output-dir "${PEDDY_DIR}"
                --genome "${GENOME}"
                --threads "${THREADS}"
            )
            if [[ -n "${REF_FASTA}" && -f "${REF_FASTA}" ]]; then
                PEDDY_ARGS+=(--src-fasta "${REF_FASTA}")
            fi
            if [[ -n "${REF_DIR}" ]]; then
                PEDDY_ARGS+=(--ref-dir "${REF_DIR}")
            fi
            if [[ -n "${PED_FILE}" && -f "${PED_FILE}" ]]; then
                PEDDY_ARGS+=(--ped-file "${PED_FILE}")
            fi
            if [[ -f "${FINAL_QC}" ]]; then
                PEDDY_ARGS+=(--sample-qc "${FINAL_QC}")
            fi
            if [[ "${PEDDY_FORCE}" == "true" ]]; then
                PEDDY_ARGS+=(--force)
            fi

            bash "${SCRIPT_DIR}/run_peddy.sh" "${PEDDY_ARGS[@]}" 2>&1 || true
        fi
    else
        echo "Note: Skipping peddy (final VCF not available)."
    fi
fi
echo ""

# ---------------------------------------------------------------
# Ancestry PCA: stringent QC + PCA computation
# ---------------------------------------------------------------
PCA_DIR="${OUTPUT_DIR}/ancestry_pca"

echo "======================================================"
echo "  Running Ancestry PCA"
echo "======================================================"
echo ""

if [[ -f "${FINAL_VCF}" ]]; then
    PCA_FORCE="${FORCE}"
    PCA_PROJ="${PCA_DIR}/pca_projections.tsv"
    if [[ "${PCA_FORCE}" != "true" && -s "${PCA_PROJ}" ]]; then
        if any_input_newer_than "${PCA_PROJ}" "${FINAL_VCF}"; then
            echo "Ancestry PCA outputs exist but inputs are newer. Regenerating."
            PCA_FORCE="true"
        else
            echo "Ancestry PCA outputs already exist and are up to date. Skipping (use --force to regenerate)."
        fi
    fi

    if [[ "${PCA_FORCE}" == "true" || ! -s "${PCA_PROJ}" ]]; then
        PCA_ARGS=(
            --vcf "${FINAL_VCF}"
            --output-dir "${PCA_DIR}"
            --threads "${THREADS}"
        )
        if [[ "${PCA_FORCE}" == "true" ]]; then
            PCA_ARGS+=(--force)
        fi

        bash "${SCRIPT_DIR}/ancestry_pca.sh" "${PCA_ARGS[@]}" 2>&1 || true
    fi
else
    echo "Note: Skipping ancestry PCA (final VCF not available)."
fi
echo ""

# ---------------------------------------------------------------
# Ancestry-Stratified QC (variant QC and PCA per ancestry)
# ---------------------------------------------------------------
ANCESTRY_QC_DIR="${OUTPUT_DIR}/ancestry_stratified_qc"

if [[ "${SKIP_PEDDY}" == "true" ]]; then
    echo "Skipping ancestry-stratified QC (peddy was skipped — ancestry predictions required)"
else
    echo "======================================================"
    echo "  Running Ancestry-Stratified QC"
    echo "======================================================"
    echo ""

    PEDDY_HET="${PEDDY_DIR}/peddy.het_check.csv"
    if [[ -f "${FINAL_VCF}" && -f "${PEDDY_HET}" ]]; then
        STRAT_FORCE="${FORCE}"
        STRAT_COLLATED="${ANCESTRY_QC_DIR}/collated_variant_qc.tsv"
        if [[ "${STRAT_FORCE}" != "true" && -s "${STRAT_COLLATED}" ]]; then
            if any_input_newer_than "${STRAT_COLLATED}" "${FINAL_VCF}" "${PEDDY_HET}"; then
                echo "Ancestry-stratified QC outputs exist but inputs are newer. Regenerating."
                STRAT_FORCE="true"
            else
                echo "Ancestry-stratified QC outputs already exist and are up to date. Skipping (use --force to regenerate)."
            fi
        fi

        if [[ "${STRAT_FORCE}" == "true" || ! -s "${STRAT_COLLATED}" ]]; then
            STRAT_ARGS=(
                --vcf "${FINAL_VCF}"
                --peddy-het-check "${PEDDY_HET}"
                --output-dir "${ANCESTRY_QC_DIR}"
                --min-samples "${MIN_ANCESTRY_SAMPLES}"
                --threads "${THREADS}"
            )
            if [[ "${STRAT_FORCE}" == "true" ]]; then
                STRAT_ARGS+=(--force)
            fi

            bash "${SCRIPT_DIR}/ancestry_stratified_qc.sh" "${STRAT_ARGS[@]}" 2>&1 || true
        fi
    else
        echo "Note: Skipping ancestry-stratified QC (final VCF or peddy output not available)."
    fi
fi
echo ""

# ---------------------------------------------------------------
# Compile unified sample sheet
# ---------------------------------------------------------------
echo "======================================================"
echo "  Compiling Sample Sheet"
echo "======================================================"
echo ""

COMPILED_SHEET="${OUTPUT_DIR}/compiled_sample_sheet.tsv"
if command -v python3 &>/dev/null && [[ -f "${FINAL_QC}" ]]; then
    if [[ "${FORCE}" != "true" && -s "${COMPILED_SHEET}" ]]; then
        echo "Compiled sample sheet already exists. Skipping (use --force to regenerate)."
    else
        COMPILE_ARGS=(
            --sample-qc "${FINAL_QC}"
            --output "${COMPILED_SHEET}"
        )
        if [[ -f "${PCA_DIR}/pca_projections.tsv" ]]; then
            COMPILE_ARGS+=(--pca-projections "${PCA_DIR}/pca_projections.tsv")
        fi

        # Include inbreeding coefficient from variant QC if available
        for het_file in "${OUTPUT_DIR}/stage2/qc/variant_qc/variant_qc.het" \
                        "${OUTPUT_DIR}/stage1/qc/variant_qc/variant_qc.het"; do
            if [[ -f "${het_file}" ]]; then
                COMPILE_ARGS+=(--het-file "${het_file}")
                break
            fi
        done

        # Include peddy sample-level metrics if available
        if [[ -f "${PEDDY_DIR}/peddy.het_check.csv" ]]; then
            COMPILE_ARGS+=(--peddy-het-check "${PEDDY_DIR}/peddy.het_check.csv")
        fi
        if [[ -f "${PEDDY_DIR}/peddy.sex_check.csv" ]]; then
            COMPILE_ARGS+=(--peddy-sex-check "${PEDDY_DIR}/peddy.sex_check.csv")
        fi

        # Include ancestry-specific PCA projections if available
        if [[ -d "${ANCESTRY_QC_DIR}" ]]; then
            for anc_dir in "${ANCESTRY_QC_DIR}"/*/pca; do
                if [[ -f "${anc_dir}/pca_projections.tsv" ]]; then
                    anc_name=$(basename "$(dirname "${anc_dir}")")
                    COMPILE_ARGS+=(--ancestry-pca "${anc_name}:${anc_dir}/pca_projections.tsv")
                fi
            done
        fi

        python3 "${SCRIPT_DIR}/compile_sample_sheet.py" "${COMPILE_ARGS[@]}" 2>&1 || true
    fi
else
    echo "Note: Skipping sample sheet compilation (python3 or QC file not available)."
fi
echo ""

echo "The VCF contains genotypes with BAF and LRR intensities"
echo "ready for downstream analysis (phasing, MoChA, imputation)."

# ---------------------------------------------------------------
# Run automated QC diagnostics
# ---------------------------------------------------------------
echo ""
echo "======================================================"
echo "  Running QC Diagnostics"
echo "======================================================"
echo ""

DIAG_REPORT="${OUTPUT_DIR}/qc_diagnostic_report.txt"
if command -v python3 &>/dev/null; then
    if [[ "${FORCE}" != "true" && -s "${DIAG_REPORT}" ]]; then
        echo "QC diagnostic report already exists. Skipping (use --force to regenerate)."
    else
        python3 "${SCRIPT_DIR}/diagnose_qc.py" \
            --output-dir "${OUTPUT_DIR}" \
            --report "${DIAG_REPORT}" 2>&1 || true
        echo ""
        echo "Diagnostic report: ${DIAG_REPORT}"
    fi
else
    echo "Note: python3 not found; skipping automated QC diagnostics."
    echo "Run manually: python3 scripts/diagnose_qc.py --output-dir ${OUTPUT_DIR}"
fi

echo ""
echo "======================================================"
echo "  Generating Pipeline Report"
echo "======================================================"
echo ""

REPORT_PATH="${OUTPUT_DIR}/pipeline_report.html"
if command -v python3 &>/dev/null; then
    if [[ "${FORCE}" != "true" && -s "${REPORT_PATH}" && -s "${OUTPUT_DIR}/summary_statistics.tsv" && -s "${OUTPUT_DIR}/methods_text.txt" && -d "${OUTPUT_DIR}/summary" ]]; then
        echo "Pipeline report outputs already exist. Skipping (use --force to regenerate)."
    else
        python3 "${SCRIPT_DIR}/generate_report.py" \
            --output-dir "${OUTPUT_DIR}" \
            --genome "${GENOME}" 2>&1 || true
        echo ""
        if [[ -f "${REPORT_PATH}" ]]; then
            echo "Pipeline report:     ${REPORT_PATH}"
        fi
        if [[ -f "${OUTPUT_DIR}/summary_statistics.tsv" ]]; then
            echo "Summary statistics:  ${OUTPUT_DIR}/summary_statistics.tsv"
        fi
        if [[ -f "${OUTPUT_DIR}/methods_text.txt" ]]; then
            echo "Methods text:        ${OUTPUT_DIR}/methods_text.txt"
        fi
        if [[ -d "${OUTPUT_DIR}/summary" ]]; then
            echo "Consolidated summary: ${OUTPUT_DIR}/summary/"
        fi
    fi
else
    echo "Note: python3 not found; skipping report generation."
fi

mkdir -p "${OUTPUT_DIR}/summary"
for f in \
    pipeline_report.html \
    summary_statistics.tsv \
    methods_text.txt \
    citations_summary.tsv \
    compiled_sample_sheet.tsv \
    qc_diagnostic_report.txt; do
    if [[ -f "${OUTPUT_DIR}/${f}" ]]; then
        cp -f "${OUTPUT_DIR}/${f}" "${OUTPUT_DIR}/summary/${f}"
    fi
done
if [[ -f "${OUTPUT_DIR}/sex_check/sex_check_chrXY_lrr.png" ]]; then
    cp -f "${OUTPUT_DIR}/sex_check/sex_check_chrXY_lrr.png" "${OUTPUT_DIR}/summary/sex_check_chrXY_lrr.png"
fi
if [[ -f "${OUTPUT_DIR}/sex_check/sex_check_chrXY_lrr.tsv" ]]; then
    cp -f "${OUTPUT_DIR}/sex_check/sex_check_chrXY_lrr.tsv" "${OUTPUT_DIR}/summary/sex_check_chrXY_lrr.tsv"
fi
# Copy peddy outputs to summary
for pf in peddy.het_check.csv peddy.sex_check.csv peddy.ped_check.csv peddy.peddy.ped peddy_final.ped; do
    if [[ -f "${OUTPUT_DIR}/peddy/${pf}" ]]; then
        cp -f "${OUTPUT_DIR}/peddy/${pf}" "${OUTPUT_DIR}/summary/${pf}"
    fi
done
# Copy ancestry-stratified QC outputs to summary
if [[ -f "${ANCESTRY_QC_DIR}/collated_variant_qc.tsv" ]]; then
    cp -f "${ANCESTRY_QC_DIR}/collated_variant_qc.tsv" "${OUTPUT_DIR}/summary/collated_variant_qc.tsv"
fi
if [[ -f "${ANCESTRY_QC_DIR}/ancestry_stratified_summary.txt" ]]; then
    cp -f "${ANCESTRY_QC_DIR}/ancestry_stratified_summary.txt" "${OUTPUT_DIR}/summary/ancestry_stratified_summary.txt"
fi

echo ""
echo "======================================================"
echo "  All Steps Complete"
echo "======================================================"
echo ""
echo "Final outputs:"
echo "  VCF:              ${FINAL_VCF}"
echo "  QC metrics:       ${FINAL_QC}"
if [[ -f "${COMPILED_SHEET}" ]]; then
    echo "  Sample sheet:     ${COMPILED_SHEET}"
fi
if [[ -d "${SEX_CHECK_DIR}" ]]; then
    echo "  Sex check:        ${SEX_CHECK_DIR}/"
fi
if [[ -d "${PEDDY_DIR}" ]]; then
    echo "  Peddy QC:         ${PEDDY_DIR}/"
fi
if [[ -d "${PCA_DIR}" ]]; then
    echo "  Ancestry PCA:     ${PCA_DIR}/"
fi
if [[ -d "${ANCESTRY_QC_DIR}" ]]; then
    echo "  Ancestry QC:      ${ANCESTRY_QC_DIR}/"
fi
if [[ -f "${REPORT_PATH}" ]]; then
    echo "  Pipeline report:  ${REPORT_PATH}"
fi
if [[ -f "${OUTPUT_DIR}/summary_statistics.tsv" ]]; then
    echo "  Summary stats:    ${OUTPUT_DIR}/summary_statistics.tsv"
fi
if [[ -f "${OUTPUT_DIR}/methods_text.txt" ]]; then
    echo "  Methods text:     ${OUTPUT_DIR}/methods_text.txt"
fi
if [[ -d "${OUTPUT_DIR}/summary" ]]; then
    echo "  Summary bundle:   ${OUTPUT_DIR}/summary/"
fi
