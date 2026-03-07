#!/usr/bin/env bash
#
# test_ancestry_pca_projection.sh
#
# Deterministic integration test for ancestry_pca.sh projection behavior.
# Ensures training samples project back to their PCA values and projected-only
# samples have finite, non-degenerate, reasonable PC coordinates.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PASS=0
FAIL=0

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

FAKE_BIN="${TMP_DIR}/bin"
OUT_DIR="${TMP_DIR}/out"
mkdir -p "${FAKE_BIN}" "${OUT_DIR}"

echo "============================================"
echo "  ancestry_pca projection integration test"
echo "============================================"
echo ""

cat > "${FAKE_BIN}/bcftools" << 'EOF'
#!/usr/bin/env bash
set -euo pipefail

cmd="${1:-}"
shift || true

case "${cmd}" in
    view)
        if [[ "${1:-}" == "-H" ]]; then
            exit 0
        fi
        out=""
        while [[ "$#" -gt 0 ]]; do
            case "$1" in
                -o) out="$2"; shift 2 ;;
                *) shift ;;
            esac
        done
        [[ -n "${out}" ]] && : > "${out}"
        ;;
    index)
        if [[ "${1:-}" == "-n" ]]; then
            file="${2:-}"
            if [[ "${file}" == *"filtered_ancestry_pca"* ]]; then
                echo "5"
            else
                echo "6"
            fi
            exit 0
        fi
        file="${1:-}"
        [[ -n "${file}" ]] && : > "${file}.csi"
        ;;
    query)
        if [[ "${1:-}" == "-l" ]]; then
            printf "HG00096\nHG00097\nHG00098\nHG00099\nHG00100\n"
        fi
        ;;
    *)
        echo "Unsupported fake bcftools command: ${cmd}" >&2
        exit 1
        ;;
esac
EOF

cat > "${FAKE_BIN}/plink2" << 'EOF'
#!/usr/bin/env bash
set -euo pipefail

out=""
is_make_bed="false"
is_ld="false"

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --out) out="$2"; shift 2 ;;
        --make-bed) is_make_bed="true"; shift ;;
        --indep-pairwise) is_ld="true"; shift 4 ;;
        *) shift ;;
    esac
done

if [[ "${is_ld}" == "true" ]]; then
    cat > "${out}.prune.in" << 'PRUNE'
rs1
rs2
rs3
rs4
PRUNE
fi

if [[ "${is_make_bed}" == "true" ]]; then
    : > "${out}.bed"
    cat > "${out}.bim" << 'BIM'
1	rs1	0	101	A	G
1	rs2	0	102	C	T
1	rs3	0	103	G	A
1	rs4	0	104	T	C
BIM

    if [[ "${out}" == *"_all_ldpruned" ]]; then
        cat > "${out}.fam" << 'FAM'
0 HG00096 0 0 0 -9
0 HG00097 0 0 0 -9
0 HG00098 0 0 0 -9
0 HG00099 0 0 0 -9
0 HG00100 0 0 0 -9
FAM
    else
        cat > "${out}.fam" << 'FAM'
0 HG00096 0 0 0 -9
0 HG00097 0 0 0 -9
0 HG00099 0 0 0 -9
FAM
    fi
fi
EOF

cat > "${FAKE_BIN}/flashpca" << 'EOF'
#!/usr/bin/env bash
set -euo pipefail

project="false"
outpc=""
outval=""
outload=""
outmeansd=""
inload=""
outproj=""

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --project) project="true"; shift ;;
        --outpc) outpc="$2"; shift 2 ;;
        --outval) outval="$2"; shift 2 ;;
        --outload) outload="$2"; shift 2 ;;
        --outmeansd) outmeansd="$2"; shift 2 ;;
        --inload) inload="$2"; shift 2 ;;
        --outproj) outproj="$2"; shift 2 ;;
        *) shift ;;
    esac
done

if [[ "${project}" == "false" ]]; then
    if [[ -z "${outload}" ]]; then
        echo "Missing required --outload argument" >&2
        exit 2
    fi

    cat > "${outpc}" << 'PCS'
0 HG00096 -0.0600 0.2500 -0.0600 0.0200
0 HG00097 -0.0620 0.2480 -0.0610 0.0210
0 HG00099 -0.0660 0.2490 -0.0580 0.0180
PCS
    cat > "${outval}" << 'EIG'
1.0
0.8
0.6
0.4
EIG
    cat > "${outload}" << 'LOAD'
0.1 0.2 0.3 0.4
0.2 0.1 0.4 0.3
0.3 0.4 0.1 0.2
0.4 0.3 0.2 0.1
LOAD
    cat > "${outmeansd}" << 'MS'
0.0 1.0
0.0 1.0
0.0 1.0
0.0 1.0
MS
else
    if [[ ! -f "${inload}" ]]; then
        echo "Loadings file not found: ${inload}" >&2
        exit 3
    fi

    cat > "${outproj}" << 'PROJ'
0 HG00096 -0.0600 0.2500 -0.0600 0.0200
0 HG00097 -0.0620 0.2480 -0.0610 0.0210
0 HG00098 -0.0610 0.2490 -0.0590 0.0190
0 HG00099 -0.0660 0.2490 -0.0580 0.0180
0 HG00100 -0.0640 0.2510 -0.0600 0.0220
PROJ
fi
EOF

chmod +x "${FAKE_BIN}/bcftools" "${FAKE_BIN}/plink2" "${FAKE_BIN}/flashpca"

echo "--- Test 1: ancestry_pca.sh deterministic projection ---"
: > "${TMP_DIR}/input.bcf"

if PATH="${FAKE_BIN}:${PATH}" bash "${REPO_DIR}/scripts/ancestry_pca.sh" \
    --vcf "${TMP_DIR}/input.bcf" \
    --output-dir "${OUT_DIR}" \
    --n-pcs 4 \
    --threads 1 >/dev/null 2>&1; then
    echo "  PASS: ancestry_pca.sh completed with mocked tools"
    (( PASS++ )) || true
else
    echo "  FAIL: ancestry_pca.sh failed with mocked tools"
    (( FAIL++ )) || true
fi

echo "--- Test 2: overlap samples match PCA values within tolerance ---"
PCS_FILE="${OUT_DIR}/ancestry_pca_pcs.txt"
PROJ_FILE="${OUT_DIR}/pca_projections.tsv"

if [[ -f "${PCS_FILE}" && -f "${PROJ_FILE}" ]] && \
   python3 - "${PCS_FILE}" "${PROJ_FILE}" << 'PY'
import csv
import math
import sys

pcs_path, proj_path = sys.argv[1], sys.argv[2]
n_pcs = 4
tolerance = 1e-9

pcs = {}
with open(pcs_path, "r", encoding="utf-8") as f:
    for line in f:
        fields = line.strip().split()
        if len(fields) < 2 + n_pcs:
            continue
        sid = fields[1]
        pcs[sid] = [float(x) for x in fields[2:2 + n_pcs]]

proj = {}
with open(proj_path, "r", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    for row in reader:
        sid = row[1]
        proj[sid] = [float(x) for x in row[2:2 + n_pcs]]

overlap = sorted(set(pcs).intersection(proj))
if len(overlap) < 3:
    raise SystemExit("Expected >=3 overlapping samples between PCs and projections")

for sid in overlap:
    for i, (a, b) in enumerate(zip(pcs[sid], proj[sid]), start=1):
        if abs(a - b) > tolerance:
            raise SystemExit(
                f"Sample {sid} PC{i} mismatch: base={a} projected={b} tolerance={tolerance}"
            )

for sid, vals in proj.items():
    if not all(math.isfinite(v) for v in vals):
        raise SystemExit(f"Sample {sid} has non-finite projected values: {vals}")
    if any(abs(v) > 10 for v in vals):
        raise SystemExit(f"Sample {sid} has unreasonable projected magnitude: {vals}")

projected_only = sorted(set(proj).difference(pcs))
if len(projected_only) < 2:
    raise SystemExit("Expected projected-only samples in projection output")

for sid in projected_only:
    vals = proj[sid]
    if max(vals) - min(vals) <= 1e-8:
        raise SystemExit(
            f"Sample {sid} projected PCs appear degenerate (all nearly identical): {vals}"
        )
PY
then
    echo "  PASS: overlapping and projected-only sample checks passed"
    (( PASS++ )) || true
else
    echo "  FAIL: projected values failed overlap/reasonableness checks"
    (( FAIL++ )) || true
fi

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
