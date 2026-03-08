#!/usr/bin/env bash
#
# run_peddy.sh
#
# Run peddy for pedigree/sex/ancestry QC on genotyped VCF output.
#
# Peddy validates sex, ancestry, and relatedness from VCF + PED input.
# If no pedigree file is provided, a minimal PED is generated (all
# samples as unrelated singletons) so that peddy can still predict
# sex, ancestry, and discover relationships.
#
# After running, a final pedigree incorporating peddy's discovered
# relationships is written.
#
# Idempotent: skips execution when expected outputs already exist
# unless --force is specified.
#
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run peddy pedigree/sex/ancestry QC on pipeline VCF output.

Required:
  --vcf FILE             Input VCF/BCF file
  --output-dir DIR       Output directory for peddy results

Optional:
  --ped-file FILE        Pedigree file (.ped/.fam). If not provided, a minimal
                         PED with unrelated singletons is generated.
  --sample-qc FILE       Sample QC TSV for sex annotation in generated PED
  --threads INT          Number of threads (default: 4)
  --force                Re-run even if outputs exist
  --help                 Show this help message
EOF
    exit 0
}

VCF=""
OUTPUT_DIR=""
PED_FILE=""
SAMPLE_QC=""
THREADS=4
FORCE="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)         VCF="$2"; shift 2 ;;
        --output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
        --ped-file)    PED_FILE="$2"; shift 2 ;;
        --sample-qc)   SAMPLE_QC="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        --force)       FORCE="true"; shift ;;
        --help)        usage ;;
        *)             echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${VCF}" ]]; then
    echo "Error: --vcf is required" >&2
    exit 1
fi
if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output-dir is required" >&2
    exit 1
fi
if [[ ! -f "${VCF}" ]]; then
    echo "Error: VCF file not found: ${VCF}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

PREFIX="${OUTPUT_DIR}/peddy"

# --- Idempotency check ---
EXPECTED_OUTPUTS=(
    "${PREFIX}.het_check.csv"
    "${PREFIX}.sex_check.csv"
    "${PREFIX}.ped_check.csv"
    "${PREFIX}.peddy.ped"
)

if [[ "${FORCE}" != "true" ]]; then
    all_exist=true
    for f in "${EXPECTED_OUTPUTS[@]}"; do
        if [[ ! -f "$f" ]]; then
            all_exist=false
            break
        fi
    done
    if [[ "${all_exist}" == "true" ]]; then
        echo "Peddy outputs already exist in ${OUTPUT_DIR}. Skipping (use --force to re-run)."
        exit 0
    fi
fi

# --- Prepare input VCF ---
# peddy requires a VCF (not BCF) with a .vcf.gz extension and tabix index.
# Convert BCF to VCF.gz if necessary.
INPUT_VCF="${VCF}"
TMP_DIR="${OUTPUT_DIR}/tmp"
mkdir -p "${TMP_DIR}"

if [[ "${VCF}" == *.bcf ]]; then
    echo "Converting BCF to VCF.gz for peddy..."
    INPUT_VCF="${TMP_DIR}/peddy_input.vcf.gz"
    if [[ ! -f "${INPUT_VCF}" || "${FORCE}" == "true" ]]; then
        bcftools view "${VCF}" -Oz -o "${INPUT_VCF}" --threads "${THREADS}"
        bcftools index -t "${INPUT_VCF}"
    fi
elif [[ "${VCF}" == *.vcf && ! "${VCF}" == *.vcf.gz ]]; then
    echo "Compressing VCF to VCF.gz for peddy..."
    INPUT_VCF="${TMP_DIR}/peddy_input.vcf.gz"
    if [[ ! -f "${INPUT_VCF}" || "${FORCE}" == "true" ]]; then
        bcftools view "${VCF}" -Oz -o "${INPUT_VCF}" --threads "${THREADS}"
        bcftools index -t "${INPUT_VCF}"
    fi
elif [[ "${VCF}" == *.vcf.gz ]]; then
    # Ensure tabix index exists
    if [[ ! -f "${VCF}.tbi" ]]; then
        bcftools index -t "${VCF}"
    fi
fi

# --- Generate PED file if not provided ---
USED_PED="${PED_FILE}"
GENERATED_PED="${OUTPUT_DIR}/generated_input.ped"

if [[ -z "${PED_FILE}" || ! -f "${PED_FILE}" ]]; then
    echo "No pedigree file provided. Generating minimal PED from VCF samples..."

    # Extract sample names from VCF
    bcftools query -l "${INPUT_VCF}" > "${TMP_DIR}/sample_list.txt"

    # Build sex map from sample QC if available
    declare -A sex_map
    if [[ -n "${SAMPLE_QC}" && -f "${SAMPLE_QC}" ]]; then
        while IFS=$'\t' read -r sid _ _ _ _ _ _ gender _; do
            case "${gender}" in
                M|male|1)   sex_map["${sid}"]="1" ;;
                F|female|2) sex_map["${sid}"]="2" ;;
                *)          sex_map["${sid}"]="0" ;;
            esac
        done < <(tail -n +2 "${SAMPLE_QC}")
    fi

    # Write minimal PED: family_id individual_id father mother sex phenotype
    # All samples are unrelated singletons with unknown phenotype
    : > "${GENERATED_PED}"
    while read -r sample; do
        sex="${sex_map[${sample}]:-0}"
        echo -e "${sample}\t${sample}\t0\t0\t${sex}\t-9" >> "${GENERATED_PED}"
    done < "${TMP_DIR}/sample_list.txt"

    USED_PED="${GENERATED_PED}"
    echo "  Generated PED: ${GENERATED_PED} ($(wc -l < "${GENERATED_PED}") samples)"
else
    echo "Using provided pedigree file: ${PED_FILE}"
    # Copy user PED so outputs are self-contained
    cp "${PED_FILE}" "${OUTPUT_DIR}/input.ped"
fi

# --- Run peddy ---
echo ""
echo "Running peddy..."
echo "  VCF:     ${INPUT_VCF}"
echo "  PED:     ${USED_PED}"
echo "  Prefix:  ${PREFIX}"
echo "  Threads: ${THREADS}"
echo ""

PEDDY_START=${SECONDS}

python3 -m peddy \
    -p "${THREADS}" \
    --plot \
    --prefix "${PREFIX}" \
    "${INPUT_VCF}" \
    "${USED_PED}" 2>&1 || {
        echo "Warning: peddy exited with non-zero status. Some outputs may be incomplete." >&2
    }

PEDDY_ELAPSED=$(( SECONDS - PEDDY_START ))
echo ""
echo "Peddy completed in $(( PEDDY_ELAPSED / 60 ))m $(( PEDDY_ELAPSED % 60 ))s"

# --- Create final pedigree with discovered relationships ---
FINAL_PED="${OUTPUT_DIR}/peddy_final.ped"

if [[ -f "${PREFIX}.peddy.ped" ]]; then
    echo "Creating final pedigree from peddy results..."
    cp "${PREFIX}.peddy.ped" "${FINAL_PED}"
    echo "  Final pedigree: ${FINAL_PED}"
fi

# --- Report output summary ---
echo ""
echo "Peddy output files:"
for f in "${EXPECTED_OUTPUTS[@]}"; do
    if [[ -f "$f" ]]; then
        row_count=$(wc -l < "$f" | tr -d ' ')
        echo "  $(basename "$f"): ${row_count} lines"
    else
        echo "  $(basename "$f"): MISSING"
    fi
done

if [[ -f "${PREFIX}.html" ]]; then
    echo "  peddy.html: interactive report"
fi
if [[ -f "${PREFIX}.background_pca.json" ]]; then
    echo "  background_pca.json: ancestry PCA data"
fi

# Cleanup temporary files
rm -rf "${TMP_DIR}"

echo ""
echo "Done."
