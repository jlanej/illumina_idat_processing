#!/usr/bin/env bash
#
# compute_sex_chr_variant_qc.sh
#
# Compute sex-aware variant QC for chrX and chrY using PLINK2.
#
# Sex chromosomes require special handling in variant QC because males
# are hemizygous for chrX (one copy outside PAR) and carry chrY, while
# females are diploid for chrX and carry no chrY.  Best-practice
# recommendations (Anderson et al. 2010; Marees et al. 2018;
# Laurie et al. 2012 Genet Epidemiol 36:384-91) include:
#
#   - chrX HWE computed on females only (diploid genotypes)
#   - chrX missingness and allele frequency computed on all samples
#   - Per-sex missingness for chrX (separate male and female call rates)
#   - chrY variant QC computed on males only
#   - PAR/XTR regions excluded (they behave like autosomes)
#
# Outputs are written to a structured directory layout that
# collate_variant_qc.py can consume via --sex-chr-qc-dir.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/utils.sh"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Compute sex-aware variant QC for chrX and chrY.

Required:
  --vcf FILE          Input VCF/BCF file
  --sample-qc FILE    Sample QC TSV with computed_gender column (M/F/1/2)
  --output-dir DIR    Output directory for sex-chr variant QC results

Options:
  --genome STRING     Genome build for PAR/XTR exclusion (default: CHM13)
  --threads INT       Number of threads (default: all available)
  --help              Show this help message

Output files:
  chrX/chrX_all.vmiss          chrX missingness (all samples, non-PAR/XTR)
  chrX/chrX_all.afreq          chrX allele frequencies (all samples)
  chrX/chrX_female.hardy       chrX HWE (females only — diploid chrX)
  chrX/chrX_female.vmiss       chrX missingness (females only)
  chrX/chrX_male.vmiss         chrX missingness (males only)
  chrY/chrY_male.vmiss         chrY missingness (males only, non-PAR)
  chrY/chrY_male.afreq         chrY allele frequencies (males only)
EOF
    exit 0
}

VCF=""
SAMPLE_QC=""
OUTPUT_DIR=""
GENOME="CHM13"
THREADS=$(nproc 2>/dev/null || echo 1)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)        VCF="$2"; shift 2 ;;
        --sample-qc)  SAMPLE_QC="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --genome)     GENOME="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --help)       usage ;;
        *)            echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${VCF}" ]]; then
    echo "Error: --vcf is required" >&2
    exit 1
fi
if [[ -z "${SAMPLE_QC}" ]]; then
    echo "Error: --sample-qc is required" >&2
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
if [[ ! -f "${SAMPLE_QC}" ]]; then
    echo "Error: Sample QC file not found: ${SAMPLE_QC}" >&2
    exit 1
fi

# Check required tools
if ! command -v plink2 &>/dev/null; then
    echo "Warning: plink2 not found. Skipping sex-chromosome variant QC." >&2
    exit 0
fi
if ! command -v bcftools &>/dev/null; then
    echo "Warning: bcftools not found. Skipping sex-chromosome variant QC." >&2
    exit 0
fi

mkdir -p "${OUTPUT_DIR}"
TMP_DIR="${OUTPUT_DIR}/tmp"
mkdir -p "${TMP_DIR}"

echo "Computing sex-chromosome variant QC..."
echo "  VCF:     ${VCF}"
echo "  Genome:  ${GENOME}"
echo "  Threads: ${THREADS}"

# -----------------------------------------------------------------
# Extract male and female sample lists from the sample QC TSV.
# The computed_gender column uses M/F or 1/2 coding.
# -----------------------------------------------------------------
MALE_SAMPLES="${TMP_DIR}/male_samples.txt"
FEMALE_SAMPLES="${TMP_DIR}/female_samples.txt"

python3 -c "
import sys
males, females = [], []
with open('${SAMPLE_QC}') as f:
    header = f.readline().strip().split('\t')
    sex_idx = None
    sid_idx = 0
    for i, col in enumerate(header):
        if col == 'computed_gender':
            sex_idx = i
        if col == 'sample_id':
            sid_idx = i
    if sex_idx is None:
        print('Warning: computed_gender column not found in sample QC', file=sys.stderr)
        sys.exit(0)
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) <= sex_idx:
            continue
        sid = fields[sid_idx]
        sex = fields[sex_idx].upper().strip()
        if sex in ('M', '1', 'MALE'):
            males.append(sid)
        elif sex in ('F', '2', 'FEMALE'):
            females.append(sid)

with open('${MALE_SAMPLES}', 'w') as out:
    for s in males:
        out.write(s + '\n')
with open('${FEMALE_SAMPLES}', 'w') as out:
    for s in females:
        out.write(s + '\n')

print(f'  Males: {len(males)}, Females: {len(females)}')
"

N_MALES=$(wc -l < "${MALE_SAMPLES}" | tr -d ' ')
N_FEMALES=$(wc -l < "${FEMALE_SAMPLES}" | tr -d ' ')

if [[ "${N_MALES}" -eq 0 && "${N_FEMALES}" -eq 0 ]]; then
    echo "Warning: No sex assignments found. Skipping sex-chromosome QC."
    exit 0
fi

# -----------------------------------------------------------------
# Write PAR/XTR exclusion BED
# -----------------------------------------------------------------
PAR_BED="${TMP_DIR}/par_xtr.bed"
get_par_xtr_bed "${GENOME}" "${PAR_BED}"

# -----------------------------------------------------------------
# Filter VCF to remove intensity-only probes
# -----------------------------------------------------------------
FILTERED_VCF=$(mktemp -p "${TMP_DIR}" --suffix=.bcf filtered_sexchr.XXXXXX)
bcftools view -e 'INFO/INTENSITY_ONLY=1' --threads "${THREADS}" \
    -Ob -o "${FILTERED_VCF}" "${VCF}" 2>/dev/null
bcftools index --threads "${THREADS}" "${FILTERED_VCF}" 2>/dev/null

# -----------------------------------------------------------------
# Determine chromosome naming convention from VCF (chrX vs X)
# -----------------------------------------------------------------
CHR_PREFIX=""
if bcftools view -h "${FILTERED_VCF}" 2>/dev/null | grep -q '##contig=<ID=chrX'; then
    CHR_PREFIX="chr"
fi

CHRX="${CHR_PREFIX}X"
CHRY="${CHR_PREFIX}Y"

# Build bcftools target regions for non-PAR chrX and non-PAR chrY.
# plink2 v2.00a6+ replaced --exclude-range with --exclude bed0.
# We use --exclude bed0 for PAR/XTR exclusion in variant QC.

# =================================================================
# chrX variant QC (non-PAR/XTR)
# =================================================================
CHRX_DIR="${OUTPUT_DIR}/chrX"
mkdir -p "${CHRX_DIR}"

echo ""
echo "--- chrX Variant QC ---"

# Check if VCF contains chrX variants (use -H to skip headers efficiently)
N_CHRX=$(bcftools view -H -r "${CHRX}" "${FILTERED_VCF}" 2>/dev/null | wc -l | tr -d ' ')
echo "  chrX variants in VCF: ${N_CHRX}"

if [[ "${N_CHRX}" -gt 0 ]]; then
    # Subset to chrX
    CHRX_VCF="${TMP_DIR}/chrX_subset.bcf"
    bcftools view -r "${CHRX}" --threads "${THREADS}" \
        -Ob -o "${CHRX_VCF}" "${FILTERED_VCF}" 2>/dev/null
    bcftools index --threads "${THREADS}" "${CHRX_VCF}" 2>/dev/null

    # chrX all-samples: missingness and allele frequency (PAR/XTR excluded)
    echo "  Computing chrX all-sample missingness and allele frequency..."
    plink2 \
        --bcf "${CHRX_VCF}" \
        --allow-extra-chr \
        --exclude bed0 "${PAR_BED}" \
        --threads "${THREADS}" \
        --missing variant-only \
        --freq \
        --out "${CHRX_DIR}/chrX_all" 2>&1 | tail -3 || true

    # chrX females-only: HWE and missingness (females are diploid on chrX)
    # HWE is only valid for diploid genotypes, so males are excluded.
    # Reference: Laurie et al. 2012, Genet Epidemiol 36:384-91
    if [[ "${N_FEMALES}" -gt 0 ]]; then
        echo "  Computing chrX female-only HWE and missingness (${N_FEMALES} females)..."
        # Create plink2 keep file (FID IID format; use sample_id for both)
        FEMALE_KEEP="${TMP_DIR}/female_keep.txt"
        awk '{print $1, $1}' "${FEMALE_SAMPLES}" > "${FEMALE_KEEP}"

        plink2 \
            --bcf "${CHRX_VCF}" \
            --allow-extra-chr \
            --exclude bed0 "${PAR_BED}" \
            --keep "${FEMALE_KEEP}" \
            --threads "${THREADS}" \
            --missing variant-only \
            --hardy midp \
            --out "${CHRX_DIR}/chrX_female" 2>&1 | tail -3 || true
    fi

    # chrX males-only: missingness
    if [[ "${N_MALES}" -gt 0 ]]; then
        echo "  Computing chrX male-only missingness (${N_MALES} males)..."
        MALE_KEEP="${TMP_DIR}/male_keep.txt"
        awk '{print $1, $1}' "${MALE_SAMPLES}" > "${MALE_KEEP}"

        plink2 \
            --bcf "${CHRX_VCF}" \
            --allow-extra-chr \
            --exclude bed0 "${PAR_BED}" \
            --keep "${MALE_KEEP}" \
            --threads "${THREADS}" \
            --missing variant-only \
            --out "${CHRX_DIR}/chrX_male" 2>&1 | tail -3 || true
    fi
else
    echo "  No chrX variants found — skipping chrX QC."
fi

# =================================================================
# chrY variant QC (non-PAR, males only)
# =================================================================
CHRY_DIR="${OUTPUT_DIR}/chrY"
mkdir -p "${CHRY_DIR}"

echo ""
echo "--- chrY Variant QC ---"

N_CHRY=$(bcftools view -H -r "${CHRY}" "${FILTERED_VCF}" 2>/dev/null | wc -l | tr -d ' ')
echo "  chrY variants in VCF: ${N_CHRY}"

if [[ "${N_CHRY}" -gt 0 && "${N_MALES}" -gt 0 ]]; then
    CHRY_VCF="${TMP_DIR}/chrY_subset.bcf"
    bcftools view -r "${CHRY}" --threads "${THREADS}" \
        -Ob -o "${CHRY_VCF}" "${FILTERED_VCF}" 2>/dev/null
    bcftools index --threads "${THREADS}" "${CHRY_VCF}" 2>/dev/null

    # chrY males-only: missingness and allele frequency (PAR excluded)
    echo "  Computing chrY male-only missingness and allele frequency (${N_MALES} males)..."
    MALE_KEEP="${TMP_DIR}/male_keep.txt"
    awk '{print $1, $1}' "${MALE_SAMPLES}" > "${MALE_KEEP}"

    plink2 \
        --bcf "${CHRY_VCF}" \
        --allow-extra-chr \
        --exclude bed0 "${PAR_BED}" \
        --keep "${MALE_KEEP}" \
        --threads "${THREADS}" \
        --missing variant-only \
        --freq \
        --out "${CHRY_DIR}/chrY_male" 2>&1 | tail -3 || true
elif [[ "${N_MALES}" -eq 0 ]]; then
    echo "  No male samples — skipping chrY QC."
else
    echo "  No chrY variants found — skipping chrY QC."
fi

# Clean up temp files
rm -rf "${TMP_DIR}"

echo ""
echo "Sex-chromosome variant QC complete."
echo "  Output: ${OUTPUT_DIR}"
for d in chrX chrY; do
    if [[ -d "${OUTPUT_DIR}/${d}" ]]; then
        for f in "${OUTPUT_DIR}/${d}"/*; do
            [[ -f "${f}" ]] && echo "    ${f}"
        done
    fi
done
