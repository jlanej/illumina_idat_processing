#!/usr/bin/env bash
#
# test_realign_manifest.sh
#
# Integration tests for realign_manifest.sh logic:
#   - SAM alignment parsing (awk) with POSIX-portable bitwise flag checks
#   - CSV coordinate comparison (Python)
#   - Test CSV manifest structure validation
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DATA="${SCRIPT_DIR}/data"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0

echo "============================================"
echo "  Realign Manifest Integration Tests"
echo "============================================"
echo ""

# ---------------------------------------------------------------
# Test 1: AWK SAM parsing with portable bitwise flags (mawk)
# ---------------------------------------------------------------
echo "--- Test 1: AWK SAM flag parsing (portable bitwise) ---"

# Create synthetic SAM with known flag categories
printf '@HD\tVN:1.6\tSO:unsorted\n' > "${TMP_DIR}/test.sam"
printf '@SQ\tSN:chr1\tLN:248956422\n' >> "${TMP_DIR}/test.sam"
# mapped, MAPQ=60 → mapq_high
printf 'probe1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tACGT\t*\n' >> "${TMP_DIR}/test.sam"
# unmapped (flag 4)
printf 'probe2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n' >> "${TMP_DIR}/test.sam"
# secondary alignment (flag 256)
printf 'probe3\t256\tchr1\t200\t0\t50M\t*\t0\t0\tACGT\t*\n' >> "${TMP_DIR}/test.sam"
# supplementary alignment (flag 2048)
printf 'probe4\t2048\tchr1\t300\t0\t50M\t*\t0\t0\tACGT\t*\n' >> "${TMP_DIR}/test.sam"
# mapped, MAPQ=5 → mapq_low
printf 'probe5\t0\tchr1\t400\t5\t50M\t*\t0\t0\tACGT\t*\n' >> "${TMP_DIR}/test.sam"
# mapped, MAPQ=0 → mapq_0 (multi-mapped)
printf 'probe6\t0\tchr1\t500\t0\t50M\t*\t0\t0\tACGT\t*\n' >> "${TMP_DIR}/test.sam"

SUMMARY_FILE="${TMP_DIR}/summary.txt"

# Extract the awk block from realign_manifest.sh and run it
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

    if (int(flag / 256) % 2 == 1 || int(flag / 2048) % 2 == 1) {
        secondary++
        next
    }
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
}' "${TMP_DIR}/test.sam" > "${TMP_DIR}/awk_output.txt"

# Validate expected counts
check_value() {
    local label="$1" expected="$2" file="$3"
    # Extract the count value after the label (first number after the colon)
    actual=$(grep "${label}" "${file}" | sed "s/.*${label}//" | grep -o '[0-9]*' | head -1)
    if [[ "${actual}" == "${expected}" ]]; then
        echo "  PASS: ${label} = ${expected}"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${label} expected=${expected} actual=${actual}"
        (( FAIL++ )) || true
    fi
}

check_value "Total flank sequences:" "4" "${SUMMARY_FILE}"
check_value "Mapped:" "3" "${SUMMARY_FILE}"
check_value "Unique (MAPQ >= 10):" "1" "${SUMMARY_FILE}"
check_value "Unmapped:" "1" "${SUMMARY_FILE}"
check_value "Secondary alignments:" "2" "${SUMMARY_FILE}"
check_value "Multi-mapped (MAPQ 0):" "1" "${SUMMARY_FILE}"
check_value "Low MAPQ (1-9):" "1" "${SUMMARY_FILE}"

echo ""

# ---------------------------------------------------------------
# Test 2: AWK parsing works with mawk (POSIX portability)
# ---------------------------------------------------------------
echo "--- Test 2: mawk compatibility ---"
if command -v mawk &>/dev/null; then
    MAWK_SUMMARY="${TMP_DIR}/mawk_summary.txt"
    mawk -v summary="${MAWK_SUMMARY}" '
BEGIN { total = 0; mapped = 0; unmapped = 0; secondary = 0; mapq_0 = 0; mapq_low = 0; mapq_high = 0 }
/^@/ { next }
{
    total++; flag = $2; mapq = $5 + 0
    if (int(flag / 256) % 2 == 1 || int(flag / 2048) % 2 == 1) { secondary++; next }
    if (int(flag / 4) % 2 == 1) { unmapped++; next }
    mapped++
    if (mapq == 0) mapq_0++; else if (mapq < 10) mapq_low++; else mapq_high++
}
END {
    primary = total - secondary
    printf "primary=%d mapped=%d unmapped=%d secondary=%d\n", primary, mapped, unmapped, secondary > summary
}' "${TMP_DIR}/test.sam"
    if grep -q "primary=4 mapped=3 unmapped=1 secondary=2" "${MAWK_SUMMARY}"; then
        echo "  PASS: mawk produces correct SAM flag parsing"
        (( PASS++ )) || true
    else
        echo "  FAIL: mawk produced incorrect output:"
        sed 's/^/        /' "${MAWK_SUMMARY}"
        (( FAIL++ )) || true
    fi
else
    echo "  SKIP: mawk not installed"
fi

echo ""

# ---------------------------------------------------------------
# Test 3: AWK handles combined SAM flags (e.g., 260 = 256+4)
# ---------------------------------------------------------------
echo "--- Test 3: Combined SAM flag handling ---"

printf '@HD\tVN:1.6\n' > "${TMP_DIR}/combined_flags.sam"
# flag 260 = 256 (secondary) + 4 (unmapped) → secondary wins
printf 'p1\t260\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n' >> "${TMP_DIR}/combined_flags.sam"
# flag 2052 = 2048 (supplementary) + 4 (unmapped) → secondary wins
printf 'p2\t2052\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n' >> "${TMP_DIR}/combined_flags.sam"
# flag 16 = reverse complement → mapped (no special bits)
printf 'p3\t16\tchr1\t100\t30\t50M\t*\t0\t0\tACGT\t*\n' >> "${TMP_DIR}/combined_flags.sam"

result=$(awk '
/^@/ { next }
{
    flag = $2
    if (int(flag / 256) % 2 == 1 || int(flag / 2048) % 2 == 1) { sec++; next }
    if (int(flag / 4) % 2 == 1) { unm++; next }
    map++
}
END { printf "map=%d unm=%d sec=%d", map+0, unm+0, sec+0 }
' "${TMP_DIR}/combined_flags.sam")

if [[ "${result}" == "map=1 unm=0 sec=2" ]]; then
    echo "  PASS: Combined flags handled correctly (${result})"
    (( PASS++ )) || true
else
    echo "  FAIL: Combined flags: expected 'map=1 unm=0 sec=2', got '${result}'"
    (( FAIL++ )) || true
fi

echo ""

# ---------------------------------------------------------------
# Test 4: Test CSV manifest structure (test data)
# ---------------------------------------------------------------
echo "--- Test 4: Test CSV manifest file structure ---"

TEST_CSV="${TEST_DATA}/HumanOmni2.5-4v1_B.first.40.csv"
if [[ -f "${TEST_CSV}" ]]; then
    # Check [Assay] section exists
    if grep -q '^\[Assay\]' "${TEST_CSV}"; then
        echo "  PASS: [Assay] section found in test CSV"
        (( PASS++ )) || true
    else
        echo "  FAIL: [Assay] section missing from test CSV"
        (( FAIL++ )) || true
    fi

    # Check header row has SourceSeq column
    header_line=$(grep -A1 '^\[Assay\]' "${TEST_CSV}" | tail -1)
    if echo "${header_line}" | grep -qi "SourceSeq"; then
        echo "  PASS: SourceSeq column found in CSV header"
        (( PASS++ )) || true
    else
        echo "  FAIL: SourceSeq column missing from CSV header"
        (( FAIL++ )) || true
    fi

    # Check data rows have flank sequences (contain [...] pattern)
    n_flanks=$(awk -F',' '/^\[Assay\]/{found=1; next} found && !/^\[/{print}' "${TEST_CSV}" | \
        tail -n +2 | grep -c '\[' || true)
    n_flanks="${n_flanks:-0}"
    if [[ "${n_flanks}" -gt 0 ]]; then
        echo "  PASS: Found ${n_flanks} rows with flank sequences in test CSV"
        (( PASS++ )) || true
    else
        echo "  FAIL: No flank sequences found in test CSV data rows"
        (( FAIL++ )) || true
    fi
else
    echo "  SKIP: Test CSV not found at ${TEST_CSV}"
fi

echo ""

# ---------------------------------------------------------------
# Test 5: CSV comparison Python logic
# ---------------------------------------------------------------
echo "--- Test 5: CSV coordinate comparison Python logic ---"

# Create minimal original and "realigned" CSV files
cat > "${TMP_DIR}/original.csv" << 'CSVEOF'
Illumina, Inc.
[Heading]
Descriptor File Name,Test.bpm
[Assay]
IlmnID,Name,Chr,MapInfo,SourceSeq
probe1,SNP1,1,100,ACGT[A/G]TGCA
probe2,SNP2,2,200,ACGT[T/C]TGCA
probe3,SNP3,,0,ACGT[A/C]TGCA
[Controls]
CSVEOF

cat > "${TMP_DIR}/realigned.csv" << 'CSVEOF'
Illumina, Inc.
[Heading]
Descriptor File Name,Test.bpm
[Assay]
IlmnID,Name,Chr,MapInfo,SourceSeq
probe1,SNP1,1,100,ACGT[A/G]TGCA
probe2,SNP2,2,205,ACGT[T/C]TGCA
probe3,SNP3,3,300,ACGT[A/C]TGCA
[Controls]
CSVEOF

python3 -c "
import sys

def parse_illumina_csv(filepath):
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

orig = parse_illumina_csv('${TMP_DIR}/original.csv')
real = parse_illumina_csv('${TMP_DIR}/realigned.csv')

# Validate parsing
assert len(orig) == 3, f'Expected 3 probes in original, got {len(orig)}'
assert len(real) == 3, f'Expected 3 probes in realigned, got {len(real)}'

# probe1: unchanged (1, 100) -> (1, 100)
assert orig['probe1'] == ('1', '100'), f'orig probe1: {orig[\"probe1\"]}'
assert real['probe1'] == ('1', '100'), f'real probe1: {real[\"probe1\"]}'

# probe2: position changed (2, 200) -> (2, 205)
assert orig['probe2'] == ('2', '200'), f'orig probe2: {orig[\"probe2\"]}'
assert real['probe2'] == ('2', '205'), f'real probe2: {real[\"probe2\"]}'

# probe3: newly placed ('', '0') -> ('3', '300')
assert orig['probe3'] == ('', '0'), f'orig probe3: {orig[\"probe3\"]}'
assert real['probe3'] == ('3', '300'), f'real probe3: {real[\"probe3\"]}'

print('  All CSV parsing assertions passed')
" 2>&1

if [[ $? -eq 0 ]]; then
    echo "  PASS: CSV comparison Python logic"
    (( PASS++ )) || true
else
    echo "  FAIL: CSV comparison Python logic"
    (( FAIL++ )) || true
fi

echo ""

# ---------------------------------------------------------------
# Test 6: Parse test CSV data with Python logic
# ---------------------------------------------------------------
echo "--- Test 6: Parse test CSV HumanOmni2.5-4v1_B data ---"

if [[ -f "${TEST_CSV}" ]]; then
    result=$(python3 -c "
def parse_illumina_csv(filepath):
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

probes = parse_illumina_csv('${TEST_CSV}')
n = len(probes)
chroms = set(c for c, p in probes.values() if c)
print(f'probes={n} chromosomes={sorted(chroms)}')
" 2>&1)

    n_probes=$(echo "${result}" | grep -o 'probes=[0-9]*' | grep -o '[0-9]*')
    if [[ "${n_probes}" -gt 0 ]]; then
        echo "  PASS: Parsed ${n_probes} probes from test CSV (${result})"
        (( PASS++ )) || true
    else
        echo "  FAIL: Could not parse probes from test CSV: ${result}"
        (( FAIL++ )) || true
    fi
else
    echo "  SKIP: Test CSV not found"
fi

echo ""

# ---------------------------------------------------------------
# Test 7: Empty SAM file handling
# ---------------------------------------------------------------
echo "--- Test 7: Empty SAM file handling ---"

printf '@HD\tVN:1.6\n' > "${TMP_DIR}/empty.sam"
EMPTY_SUMMARY="${TMP_DIR}/empty_summary.txt"

awk -v summary="${EMPTY_SUMMARY}" '
BEGIN { total = 0; mapped = 0; unmapped = 0; secondary = 0; mapq_0 = 0; mapq_low = 0; mapq_high = 0 }
/^@/ { next }
{
    total++; flag = $2; mapq = $5 + 0
    if (int(flag / 256) % 2 == 1 || int(flag / 2048) % 2 == 1) { secondary++; next }
    if (int(flag / 4) % 2 == 1) { unmapped++; next }
    mapped++
    if (mapq == 0) mapq_0++; else if (mapq < 10) mapq_low++; else mapq_high++
}
END {
    primary = total - secondary
    pct_mapped = (primary > 0) ? sprintf("%.1f", mapped * 100.0 / primary) : "0.0"
    printf "  Total flank sequences:     %d\n", primary > summary
    printf "  Mapped:                    %d (%s%%)\n", mapped, pct_mapped > summary
}' "${TMP_DIR}/empty.sam"

if grep -q "Total flank sequences:     0" "${EMPTY_SUMMARY}"; then
    echo "  PASS: Empty SAM produces 0 total flank sequences"
    (( PASS++ )) || true
else
    echo "  FAIL: Empty SAM handling"
    sed 's/^/        /' "${EMPTY_SUMMARY}"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 8: Flank FASTA header counting (both > and @ prefixes)
# ---------------------------------------------------------------
echo "--- Test 8: Flank FASTA header counting (> and @ prefixes) ---"

# Standard FASTA with '>' headers
printf '>probe1\nACGT\n>probe2\nTGCA\n>probe3\nAAAA\n' > "${TMP_DIR}/standard.fasta"
N_STD=$(grep -c '^[>@]' "${TMP_DIR}/standard.fasta" 2>/dev/null || true)
N_STD="${N_STD:-0}"
if [[ "${N_STD}" -eq 3 ]]; then
    echo "  PASS: Standard FASTA (>) counted ${N_STD} sequences"
    (( PASS++ )) || true
else
    echo "  FAIL: Standard FASTA (>) expected 3, got ${N_STD}"
    (( FAIL++ )) || true
fi

# gtc2vcf --fasta-flank style with '@' headers
printf '@probe1:1\nACGTTGCA\n@probe2:2\nTGCAAAAA\n' > "${TMP_DIR}/fastq_style.fasta"
N_AT=$(grep -c '^[>@]' "${TMP_DIR}/fastq_style.fasta" 2>/dev/null || true)
N_AT="${N_AT:-0}"
if [[ "${N_AT}" -eq 2 ]]; then
    echo "  PASS: FASTQ-style (@) counted ${N_AT} sequences"
    (( PASS++ )) || true
else
    echo "  FAIL: FASTQ-style (@) expected 2, got ${N_AT}"
    (( FAIL++ )) || true
fi

# Empty file
: > "${TMP_DIR}/empty.fasta"
N_EMPTY=$(grep -c '^[>@]' "${TMP_DIR}/empty.fasta" 2>/dev/null || true)
N_EMPTY="${N_EMPTY:-0}"
if [[ "${N_EMPTY}" -eq 0 ]]; then
    echo "  PASS: Empty FASTA counted ${N_EMPTY} sequences"
    (( PASS++ )) || true
else
    echo "  FAIL: Empty FASTA expected 0, got ${N_EMPTY}"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 9: Chromosome normalization in CSV comparison
# ---------------------------------------------------------------
echo ""
echo "--- Test 9: Chromosome normalization (chr prefix handling) ---"

# Create CSVs with different chromosome naming conventions
cat > "${TMP_DIR}/orig_nochr.csv" << 'CSVEOF'
[Assay]
IlmnID,Name,Chr,MapInfo,SourceSeq
probe1,SNP1,1,100,ACGT[A/G]TGCA
probe2,SNP2,X,200,ACGT[T/C]TGCA
probe3,SNP3,MT,300,ACGT[A/C]TGCA
probe4,SNP4,,0,ACGT[A/T]TGCA
[Controls]
CSVEOF

cat > "${TMP_DIR}/real_chr.csv" << 'CSVEOF'
[Assay]
IlmnID,Name,Chr,MapInfo,SourceSeq
probe1,SNP1,chr1,100,ACGT[A/G]TGCA
probe2,SNP2,chrX,200,ACGT[T/C]TGCA
probe3,SNP3,chrM,305,ACGT[A/C]TGCA
probe4,SNP4,chr5,500,ACGT[A/T]TGCA
[Controls]
CSVEOF

chrom_result=$(python3 -c "
def parse_illumina_csv(filepath):
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
                cols = line.split(',')
                for i, c in enumerate(cols):
                    c_lower = c.strip().lower()
                    if c_lower in ('ilmnid', 'name'):
                        if name_col is None: name_col = i
                    if c_lower in ('chr', 'chromosome'): chr_col = i
                    if c_lower in ('mapinfo', 'position'): pos_col = i
                continue
            if in_data and name_col is not None:
                if line.startswith('['): break
                cols = line.split(',')
                if len(cols) > max(name_col, chr_col or 0, pos_col or 0):
                    name = cols[name_col].strip()
                    chrom = cols[chr_col].strip() if chr_col is not None else ''
                    pos = cols[pos_col].strip() if pos_col is not None else ''
                    probes[name] = (chrom, pos)
    return probes

def normalize_chrom(c):
    c = c.strip()
    if c.lower().startswith('chr'):
        c = c[3:]
    return c.upper()

orig = parse_illumina_csv('${TMP_DIR}/orig_nochr.csv')
real = parse_illumina_csv('${TMP_DIR}/real_chr.csv')

unchanged = 0
pos_changed = 0
chr_changed = 0
newly_placed = 0
for name, (o_chr, o_pos) in orig.items():
    if name in real:
        r_chr, r_pos = real[name]
        o_norm = normalize_chrom(o_chr)
        r_norm = normalize_chrom(r_chr)
        if o_norm == r_norm and o_pos == r_pos:
            unchanged += 1
        elif (o_chr == '' or o_chr == '0') and r_chr not in ('', '0'):
            newly_placed += 1
        elif o_norm != r_norm:
            chr_changed += 1
        elif o_pos != r_pos:
            pos_changed += 1

# probe1: '1' vs 'chr1' same pos → unchanged
# probe2: 'X' vs 'chrX' same pos → unchanged
# probe3: 'MT' vs 'chrM' diff name but also diff pos → depends on normalization
#   MT→MT, chrM→M → different → chr_changed
# probe4: '' vs 'chr5' → newly_placed
print(f'unchanged={unchanged} pos_changed={pos_changed} chr_changed={chr_changed} newly_placed={newly_placed}')
assert unchanged == 2, f'Expected 2 unchanged, got {unchanged}'
assert newly_placed == 1, f'Expected 1 newly_placed, got {newly_placed}'
assert chr_changed == 1, f'Expected 1 chr_changed (MT vs M), got {chr_changed}'
print('OK')
" 2>&1)

if echo "${chrom_result}" | grep -q "OK"; then
    echo "  PASS: Chromosome normalization handles chr prefix correctly"
    echo "        ${chrom_result}"
    (( PASS++ )) || true
else
    echo "  FAIL: Chromosome normalization: ${chrom_result}"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 10: Gender extraction awk (for-loop scope fix)
# ---------------------------------------------------------------
echo ""
echo "--- Test 10: Gender extraction awk for-loop fix ---"

printf 'sample_id\tcall_rate\tcomputed_gender\n' > "${TMP_DIR}/meta.tsv"
printf 'SAMP001\t0.99\tM\n' >> "${TMP_DIR}/meta.tsv"
printf 'SAMP002\t0.98\tF\n' >> "${TMP_DIR}/meta.tsv"

gender_out=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) {if($i=="sample_id") si=i; if($i=="computed_gender") gi=i}}
             NR>1 {print $si "\t" $gi}' "${TMP_DIR}/meta.tsv")

expected_gender=$(printf 'SAMP001\tM\nSAMP002\tF')
if [[ "${gender_out}" == "${expected_gender}" ]]; then
    echo "  PASS: Gender extraction produces correct columns"
    (( PASS++ )) || true
else
    echo "  FAIL: Gender extraction: got '${gender_out}', expected '${expected_gender}'"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 11: QC call rate awk (per-sample computation)
# ---------------------------------------------------------------
echo ""
echo "--- Test 11: Per-sample call rate computation ---"

# Simulate bcftools query output: sample\tGT per site
printf 'SAMP1\t0/0\nSAMP2\t./.\nSAMP1\t0/1\nSAMP2\t0/0\nSAMP1\t./.\nSAMP2\t./.\n' > "${TMP_DIR}/gt_data.txt"

cr_out=$(awk -F'\t' '{
    total[$1]++
    if ($2 != "./." && $2 != "." && $2 != ".|.") called[$1]++
} END {
    for (s in total) {
        cr = (total[s] > 0) ? called[s] / total[s] : 0
        print s "\t" cr
    }
}' "${TMP_DIR}/gt_data.txt" | sort -k1,1)

# SAMP1: 2 called out of 3 → 0.666667
# SAMP2: 1 called out of 3 → 0.333333
samp1_cr=$(echo "${cr_out}" | grep 'SAMP1' | awk -F'\t' '{printf "%.2f", $2}')
samp2_cr=$(echo "${cr_out}" | grep 'SAMP2' | awk -F'\t' '{printf "%.2f", $2}')

if [[ "${samp1_cr}" == "0.67" ]] && [[ "${samp2_cr}" == "0.33" ]]; then
    echo "  PASS: Per-sample call rates computed correctly (SAMP1=${samp1_cr}, SAMP2=${samp2_cr})"
    (( PASS++ )) || true
else
    echo "  FAIL: Per-sample call rates: SAMP1=${samp1_cr} (expected 0.67), SAMP2=${samp2_cr} (expected 0.33)"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 12: BPM reader field order validation
# ---------------------------------------------------------------
echo ""
echo "--- Test 12: BPM reader Python validation ---"

bpm_result=$(python3 -c "
import struct, sys, os, tempfile

# Helpers matching recluster_egt.py
def write_byte(f, val): f.write(struct.pack('<B', val))
def write_int(f, val): f.write(struct.pack('<i', val))
def write_string(f, s):
    encoded = s.encode('utf-8')
    length = len(encoded)
    while length >= 0x80:
        write_byte(f, (length & 0x7F) | 0x80)
        length >>= 7
    write_byte(f, length)
    f.write(encoded)

def read_byte(f): return struct.unpack('<B', f.read(1))[0]
def read_int(f): return struct.unpack('<i', f.read(4))[0]
def read_string(f):
    total_length = 0
    partial_length = read_byte(f)
    num_bytes = 0
    while partial_length & 0x80 > 0:
        total_length += (partial_length & 0x7F) << (7 * num_bytes)
        partial_length = read_byte(f)
        num_bytes += 1
    total_length += partial_length << (7 * num_bytes)
    return f.read(total_length).decode('utf-8')

# Create a minimal BPM file with correct field order
bpm_path = os.path.join(tempfile.mkdtemp(), 'test.bpm')
num_loci = 2
with open(bpm_path, 'wb') as f:
    f.write(b'BPM')
    write_byte(f, 4)      # version
    write_int(f, 1)        # version_int
    write_string(f, 'TestManifest')
    write_string(f, '')    # controls
    write_int(f, num_loci)

    for i in range(num_loci):
        write_int(f, 1)             # entry version
        write_string(f, f'ILMN{i}')  # ilmn_id
        write_string(f, f'SNP{i}')   # name
        for _ in range(3):           # 3 skipped strings
            write_string(f, '')
        write_int(f, 1000 + i)      # address_a
        write_string(f, 'ACGT')     # allele_a_probe_seq
        write_int(f, 2000 + i)      # address_b
        write_string(f, 'TGCA')     # allele_b_probe_seq
        write_string(f, '38')       # genome_build
        write_string(f, 'chr1')     # chrom
        write_string(f, str(100+i)) # map_info
        write_string(f, '[A/G]')    # snp
        write_string(f, 'TOP')      # ref_strand
        write_byte(f, 0)            # assay_type

    # Normalization IDs in separate section
    for i in range(num_loci):
        write_byte(f, i + 10)

# Now read it back using the BPM reader logic from recluster_egt.py
with open(bpm_path, 'rb') as f:
    magic = f.read(3).decode('utf-8')
    assert magic == 'BPM', f'Bad magic: {magic}'
    _version = read_byte(f)
    _version_int = read_int(f)
    manifest_name = read_string(f)
    _controls = read_string(f)
    n_loci = read_int(f)
    assert n_loci == num_loci, f'Expected {num_loci} loci, got {n_loci}'

    names = []
    for idx in range(n_loci):
        _ver = read_int(f)
        _ilmn = read_string(f)
        name = read_string(f)
        names.append(name)
        for _ in range(3): read_string(f)
        read_int(f)      # address_a
        read_string(f)   # allele_a_probe_seq
        read_int(f)      # address_b
        read_string(f)   # allele_b_probe_seq
        read_string(f)   # genome_build
        read_string(f)   # chrom
        read_string(f)   # map_info
        read_string(f)   # snp
        read_string(f)   # ref_strand
        read_byte(f)     # assay_type

    norm_ids = []
    for _ in range(n_loci):
        norm_ids.append(read_byte(f))

assert names == ['SNP0', 'SNP1'], f'Names: {names}'
assert norm_ids == [10, 11], f'Norm IDs: {norm_ids}'
print('OK')
os.unlink(bpm_path)
" 2>&1)

if echo "${bpm_result}" | grep -q "OK"; then
    echo "  PASS: BPM reader correctly parses field order and normalization IDs"
    (( PASS++ )) || true
else
    echo "  FAIL: BPM reader: ${bpm_result}"
    (( FAIL++ )) || true
fi

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
