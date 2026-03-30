#!/usr/bin/env bash
#
# realign_manifest.sh
#
# Realign CSV manifest probe flank sequences against the reference genome.
#
# Following the MoChA/gtc2vcf recommended approach, this script always
# realigns probe flanking sequences against the target reference genome,
# regardless of whether the manifest claims to be on the same build.
# This ensures probe positions are validated and correctly placed, catching
# any discrepancies between the manifest coordinates and the actual
# reference sequence.
#
# Workflow:
#   1. Extract flank sequences from CSV manifest (bcftools +gtc2vcf --fasta-flank)
#   2. Align flanks to reference genome (bwa mem -M)
#   3. Generate realigned CSV with updated coordinates (bcftools +gtc2vcf --sam-flank)
#   4. Output detailed mapping summary
#
# Usage:
#   ./realign_manifest.sh --csv manifest.csv --ref-fasta ref.fna --output-dir manifests/
#
set -euo pipefail

CSV=""
REF_FASTA=""
OUTPUT_DIR=""
THREADS=1

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Realign CSV manifest probe flank sequences against the reference genome.

Required:
  --csv FILE             Illumina CSV manifest file
  --ref-fasta FILE       Reference genome FASTA file (must have BWA index)
  --output-dir DIR       Directory for realigned manifest and reports

Options:
  --threads INT          Number of threads for BWA alignment (default: 1)
  --help                 Show this help message

The BWA index (.bwt, .sa, .ann, .amb, .pac) must exist alongside the
reference FASTA. If not present, run: bwa index <ref.fasta>

Output:
  <output-dir>/<csv_basename>.realigned.csv   Realigned CSV manifest
  <output-dir>/realign_summary.txt            Mapping summary report
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --csv)        CSV="$2"; shift 2 ;;
        --ref-fasta)  REF_FASTA="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --help)       usage ;;
        *)            echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${CSV}" ]]; then
    echo "Error: --csv is required" >&2
    exit 1
fi
if [[ -z "${REF_FASTA}" ]]; then
    echo "Error: --ref-fasta is required" >&2
    exit 1
fi
if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output-dir is required" >&2
    exit 1
fi

if [[ ! -f "${CSV}" ]]; then
    echo "Error: CSV file not found: ${CSV}" >&2
    exit 1
fi
if [[ ! -f "${REF_FASTA}" ]]; then
    echo "Error: Reference FASTA not found: ${REF_FASTA}" >&2
    exit 1
fi

# Check for BWA index
if [[ ! -f "${REF_FASTA}.bwt" ]]; then
    echo "Error: BWA index not found for ${REF_FASTA}" >&2
    echo "Run: bwa index ${REF_FASTA}" >&2
    exit 1
fi

# Check required tools
for cmd in bcftools bwa; do
    if ! command -v "${cmd}" &>/dev/null; then
        echo "Error: ${cmd} not found in PATH" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}"

CSV_BASENAME=$(basename "${CSV}" .csv)
CSV_BASENAME="${CSV_BASENAME%.gz}"
REALIGNED_CSV="${OUTPUT_DIR}/${CSV_BASENAME}.realigned.csv"
SUMMARY_FILE="${OUTPUT_DIR}/realign_summary.txt"
FLANK_FASTA="${OUTPUT_DIR}/flanks.fasta"
FLANK_SAM="${OUTPUT_DIR}/flanks.sam"
FLANK_LOG="${OUTPUT_DIR}/fasta_flank.log"
BWA_LOG="${OUTPUT_DIR}/bwa_mem.log"

# Idempotency: skip if realigned CSV already exists and is non-empty
if [[ -f "${REALIGNED_CSV}" && -s "${REALIGNED_CSV}" && -f "${SUMMARY_FILE}" ]]; then
    echo "============================================"
    echo "  Manifest Probe Realignment (cached)"
    echo "============================================"
    echo "  Realigned CSV already exists: ${REALIGNED_CSV}"
    echo "  Summary: ${SUMMARY_FILE}"
    echo "  Skipping realignment. Delete these files to force re-run."
    echo ""
    exit 0
fi

echo "============================================"
echo "  Manifest Probe Realignment"
echo "============================================"
echo "  CSV:       ${CSV}"
echo "  Reference: ${REF_FASTA}"
echo "  Output:    ${REALIGNED_CSV}"
echo "============================================"
echo ""

# ---------------------------------------------------------------
# Step 1: Extract flank sequences from CSV manifest
# ---------------------------------------------------------------
echo "--- Step 1/4: Extracting flank sequences from CSV manifest ---"

FLANK_RC=0
bcftools +gtc2vcf \
    --csv "${CSV}" \
    --fasta-flank \
    -o "${FLANK_FASTA}" 2>"${FLANK_LOG}" || FLANK_RC=$?

# Show first few lines of plugin output for diagnostics
if [[ -f "${FLANK_LOG}" ]]; then
    head -20 "${FLANK_LOG}" | sed 's/^/  /'
fi

if [[ "${FLANK_RC}" -ne 0 && ! -s "${FLANK_FASTA}" ]]; then
    echo "ERROR: bcftools +gtc2vcf --fasta-flank failed (exit code ${FLANK_RC}) and produced no output." >&2
    echo "  Check that the CSV file contains an [Assay] section with a SourceSeq column." >&2
    echo "  Verify with: head -20 '${CSV}'" >&2
    exit 1
fi

# Count sequence headers: FASTA uses '>' prefix, gtc2vcf --fasta-flank uses '@' prefix
N_FLANKS=$(grep -c '^[>@]' "${FLANK_FASTA}" 2>/dev/null || true)
N_FLANKS="${N_FLANKS:-0}"
echo "  Extracted ${N_FLANKS} flank sequences"
if [[ -f "${FLANK_FASTA}" ]]; then
    FASTA_SIZE=$(wc -c < "${FLANK_FASTA}")
    echo "  FASTA file size: ${FASTA_SIZE} bytes"
    echo "  FASTA first 3 lines:"
    head -3 "${FLANK_FASTA}" | sed 's/^/    /'
fi

if [[ "${N_FLANKS}" -eq 0 ]]; then
    echo "ERROR: No flank sequences extracted from CSV manifest. Cannot realign." >&2
    echo "  Check that the CSV file contains an [Assay] section with a SourceSeq column." >&2
    echo "  Verify with: head -20 '${CSV}'" >&2
    exit 1
fi
echo ""

# ---------------------------------------------------------------
# Step 2: Align flank sequences to reference genome
# ---------------------------------------------------------------
echo "--- Step 2/4: Aligning flank sequences to reference ---"
echo "  This step aligns ${N_FLANKS} probe flanks against the reference genome."
echo "  With ${THREADS} threads, this may take 10-30+ minutes for large arrays."
echo "  Running: bwa mem -M -t ${THREADS} ..."

BWA_START=${SECONDS}
bwa mem -M -t "${THREADS}" "${REF_FASTA}" "${FLANK_FASTA}" > "${FLANK_SAM}" 2>"${BWA_LOG}"
BWA_ELAPSED=$(( SECONDS - BWA_START ))
BWA_MINS=$(( BWA_ELAPSED / 60 ))
BWA_SECS=$(( BWA_ELAPSED % 60 ))

echo "  Alignment complete (${BWA_MINS}m ${BWA_SECS}s)"
SAM_LINES=$(grep -cv '^@' "${FLANK_SAM}" 2>/dev/null || true)
SAM_LINES="${SAM_LINES:-0}"
echo "  SAM alignment records: ${SAM_LINES}"
echo "  BWA log tail:"
tail -3 "${BWA_LOG}" | sed 's/^/    /'
echo ""

# ---------------------------------------------------------------
# Step 3: Parse SAM alignment and generate mapping summary
# ---------------------------------------------------------------
echo "--- Step 3/4: Analyzing alignment results ---"

# Count alignment categories from SAM (skip header lines)
awk -v summary="${SUMMARY_FILE}" '
BEGIN {
    total = 0; mapped = 0; unmapped = 0; secondary = 0
    mapq_0 = 0; mapq_low = 0; mapq_high = 0
}
/^@/ { next }
{
    total++
    flag = $2
    mapq = $5 + 0

    # Check if secondary/supplementary alignment (POSIX-portable bitwise check)
    if (int(flag / 256) % 2 == 1 || int(flag / 2048) % 2 == 1) {
        secondary++
        next
    }

    # Check if unmapped
    if (int(flag / 4) % 2 == 1) {
        unmapped++
        next
    }

    mapped++

    if (mapq == 0) mapq_0++
    else if (mapq < 10) mapq_low++
    else mapq_high++
}
END {
    primary = total - secondary
    pct_mapped = (primary > 0) ? sprintf("%.1f", mapped * 100.0 / primary) : "0.0"
    pct_unmapped = (primary > 0) ? sprintf("%.1f", unmapped * 100.0 / primary) : "0.0"
    pct_unique = (primary > 0) ? sprintf("%.1f", mapq_high * 100.0 / primary) : "0.0"

    report = "============================================\n"
    report = report "  Manifest Realignment Summary\n"
    report = report "============================================\n"
    report = report sprintf("  Total flank sequences:     %d\n", primary)
    report = report sprintf("  Mapped:                    %d (%s%%)\n", mapped, pct_mapped)
    report = report sprintf("    Unique (MAPQ >= 10):     %d (%s%%)\n", mapq_high, pct_unique)
    report = report sprintf("    Low MAPQ (1-9):          %d\n", mapq_low)
    report = report sprintf("    Multi-mapped (MAPQ 0):   %d\n", mapq_0)
    report = report sprintf("  Unmapped:                  %d (%s%%)\n", unmapped, pct_unmapped)
    report = report sprintf("  Secondary alignments:      %d\n", secondary)
    report = report "============================================\n"

    printf "%s", report > summary
    printf "%s", report
}
' "${FLANK_SAM}"

echo ""

# ---------------------------------------------------------------
# Step 4: Generate realigned CSV with updated coordinates
# ---------------------------------------------------------------
echo "--- Step 4/4: Generating realigned CSV manifest ---"

bcftools +gtc2vcf \
    --csv "${CSV}" \
    --sam-flank "${FLANK_SAM}" \
    --verbose \
    -o "${REALIGNED_CSV}"

echo "  Realigned CSV: ${REALIGNED_CSV}"

# Strict validation: abort if realigned CSV is missing, empty, or malformed
if [[ ! -f "${REALIGNED_CSV}" || ! -s "${REALIGNED_CSV}" ]]; then
    echo "ERROR: Realigned CSV file was not created or is empty." >&2
    echo "  The bcftools +gtc2vcf --sam-flank step may have failed." >&2
    echo "  Check the SAM alignment file: ${FLANK_SAM}" >&2
    exit 1
fi

RCSV_SIZE=$(wc -c < "${REALIGNED_CSV}")
RCSV_LINES=$(wc -l < "${REALIGNED_CSV}")
echo "  Realigned CSV size: ${RCSV_SIZE} bytes, ${RCSV_LINES} lines"

if ! grep -q '^\[Assay\]' "${REALIGNED_CSV}"; then
    echo "ERROR: [Assay] section NOT found in realigned CSV." >&2
    echo "  The realigned manifest is malformed and cannot be used." >&2
    echo "  First 10 lines of realigned CSV:" >&2
    head -10 "${REALIGNED_CSV}" | sed 's/^/    /' >&2
    exit 1
fi

echo "  [Assay] section: found"
ASSAY_HEADER=$(grep -A1 '^\[Assay\]' "${REALIGNED_CSV}" | tail -1)
echo "  Column header: ${ASSAY_HEADER:0:120}..."
DATA_LINES=$(awk '/^\[Assay\]/{found=1; next} found && !/^\[/{n++} found && /^\[/{exit} END{print n+0}' "${REALIGNED_CSV}")
echo "  Data lines after [Assay] header: ${DATA_LINES}"

if [[ "${DATA_LINES}" -le 2 ]]; then
    echo "ERROR: Realigned CSV contains ${DATA_LINES} data lines (expected thousands)." >&2
    echo "  The realignment likely failed to produce valid probe coordinates." >&2
    awk '/^\[Assay\]/{found=1; next} found{print "    " $0; if(++n>=5) exit}' "${REALIGNED_CSV}" >&2
    exit 1
fi

# ---------------------------------------------------------------
# Compare original and realigned CSV to report coordinate changes
# ---------------------------------------------------------------
echo ""
echo "--- Coordinate comparison (original vs realigned) ---"

# Extract IlmnID, Chr, MapInfo columns from both CSVs and compare
# CSV format has a header section before [Assay]; data starts after [Assay] header row
compare_csvs() {
    local original="$1"
    local realigned="$2"

    # Parse data section of Illumina CSV (after [Assay] header, columns vary)
    # Find the header row and extract Name, Chr, MapInfo columns
    python3 -c "
import csv
import sys

def parse_illumina_csv(filepath):
    \"\"\"Parse Illumina manifest CSV, returning dict of Name -> (Chr, MapInfo).\"\"\"
    probes = {}
    in_data = False
    chr_col = None
    pos_col = None
    name_col = None

    with open(filepath, 'r', errors='replace') as f:
        for line in f:
            line = line.strip()
            if line.startswith('[Assay]'):
                in_data = True
                continue
            if in_data and name_col is None:
                # This is the header row
                cols = line.split(',')
                for i, c in enumerate(cols):
                    c_lower = c.strip().lower()
                    if c_lower in ('ilmnid', 'name'):
                        if name_col is None:
                            name_col = i
                    if c_lower in ('chr', 'chromosome'):
                        chr_col = i
                    if c_lower in ('mapinfo', 'position'):
                        pos_col = i
                continue
            if in_data and name_col is not None:
                if line.startswith('['):
                    break
                cols = line.split(',')
                if len(cols) > max(name_col, chr_col or 0, pos_col or 0):
                    name = cols[name_col].strip()
                    chrom = cols[chr_col].strip() if chr_col is not None else ''
                    pos = cols[pos_col].strip() if pos_col is not None else ''
                    probes[name] = (chrom, pos)
    return probes

orig = parse_illumina_csv('$original')
real = parse_illumina_csv('$realigned')

def normalize_chrom(c):
    \"\"\"Normalize chromosome name for comparison (strip 'chr' prefix, uppercase).\"\"\"
    c = c.strip()
    if c.lower().startswith('chr'):
        c = c[3:]
    return c.upper()

total = len(orig)
matched = 0
chr_changed = 0
pos_changed = 0
newly_placed = 0
lost = 0
unchanged = 0

# Diagnostic: show first few chromosome values for debugging
import itertools
orig_chroms = set(c for c, p in itertools.islice(orig.values(), 100))
real_chroms = set(c for c, p in itertools.islice(real.values(), 100))
print(f'  [diag] Original CSV sample chroms: {sorted(orig_chroms)[:10]}')
print(f'  [diag] Realigned CSV sample chroms: {sorted(real_chroms)[:10]}')

for name, (o_chr, o_pos) in orig.items():
    if name in real:
        r_chr, r_pos = real[name]
        o_chr_norm = normalize_chrom(o_chr)
        r_chr_norm = normalize_chrom(r_chr)
        if o_chr_norm == r_chr_norm and o_pos == r_pos:
            unchanged += 1
        elif (o_chr == '' or o_chr == '0') and r_chr not in ('', '0'):
            newly_placed += 1
        elif o_chr_norm != r_chr_norm:
            chr_changed += 1
        elif o_pos != r_pos:
            pos_changed += 1
    else:
        lost += 1

gained = len(real) - len(orig) + lost

pct = lambda n: f'{n/total*100:.1f}' if total > 0 else '0.0'

summary = []
summary.append(f'  Total probes in original:  {total}')
summary.append(f'  Total probes in realigned: {len(real)}')
summary.append(f'  Unchanged coordinates:     {unchanged} ({pct(unchanged)}%)')
summary.append(f'  Position changed:          {pos_changed} ({pct(pos_changed)}%)')
summary.append(f'  Chromosome changed:        {chr_changed} ({pct(chr_changed)}%)')
summary.append(f'  Newly placed:              {newly_placed} ({pct(newly_placed)}%)')
summary.append(f'  Lost (not in realigned):   {lost} ({pct(lost)}%)')

for line in summary:
    print(line)

# Append to summary file
with open('$SUMMARY_FILE', 'a') as sf:
    sf.write('\n--- Coordinate Comparison ---\n')
    for line in summary:
        sf.write(line + '\n')
" 2>&1 || echo "  (CSV comparison requires Python 3)"
}

compare_csvs "${CSV}" "${REALIGNED_CSV}"

# Append summary file location
{
    echo ""
    echo "Summary saved to: ${SUMMARY_FILE}"
} >> "${SUMMARY_FILE}"

echo ""
echo "  Full summary: ${SUMMARY_FILE}"
echo ""

# Clean up intermediate files
rm -f "${FLANK_FASTA}" "${FLANK_SAM}" "${BWA_LOG}" "${FLANK_LOG}"
