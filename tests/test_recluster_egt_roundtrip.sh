#!/usr/bin/env bash
#
# test_recluster_egt_roundtrip.sh
#
# Regression test for the Stage 2 reclustering EGT writer.
#
# The reference gtc2vcf/idat2gtc reader (gtc2vcf.c) requires that
# EGT files with data_block_version == 9 contain a trailing
# per-record float block of size num_records * sizeof(float) after
# the cluster counts.  When this block is missing, idat2gtc fails
# with errors of the form:
#
#     Failed to reposition stream forward <num_records*4> bytes
#
# for every sample, exactly as reported in the upstream issue
# ("Stage 2 egt failure"), where the chip in question
# (InfiniumCoreExome-24v1-3, ~551,004 probes) yields the canonical
# 2,204,016-byte shortfall.
#
# This test:
#   1. Builds a synthetic EGT with data_block_version 9 mimicking
#      the layout consumed by gtc2vcf, including the trailing float
#      block.
#   2. Round-trips it through EGTFile.read / EGTFile.write.
#   3. Asserts the written file's byte size and trailing block
#      contents match the original (reproducing the bug fix).
#   4. Repeats the round-trip for data_block_version 5, 7, 8 to make
#      sure we don't accidentally write a trailing block for the
#      versions that don't have one.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PASS=0
FAIL=0

check() {
    local label="$1"
    local result="$2"
    if [[ "${result}" == "0" ]]; then
        echo "  PASS: ${label}"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${label}"
        (( FAIL++ )) || true
    fi
}

echo "============================================"
echo "  Reclustering EGT Round-Trip Tests"
echo "============================================"
echo ""

for VERSION in 9 8 7 5; do
    echo "--- Round-trip test: data_block_version=${VERSION} ---"
    python3 - "${REPO_DIR}" "${VERSION}" <<'PYEOF'
import os
import struct
import sys
import tempfile

sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
data_block_version = int(sys.argv[2])

import numpy as np
from recluster_egt import EGTFile, write_byte, write_int, write_float, write_string


NUM_RECORDS = 37  # small but >1 so all loops exercise


def build_synthetic_egt(path, version):
    """Write a synthetic EGT file matching the layout consumed by
    gtc2vcf/idat2gtc (gtc2vcf.c egt_init).  This is the inverse of
    EGTFile.read and includes the trailing per-record float block
    only when version == 9.
    """
    rng = np.random.default_rng(seed=42 + version)
    locus_names = [f"rs{1000 + i:06d}" for i in range(NUM_RECORDS)]

    aa_n = rng.integers(2, 100, size=NUM_RECORDS).astype(np.int32)
    ab_n = rng.integers(2, 100, size=NUM_RECORDS).astype(np.int32)
    bb_n = rng.integers(2, 100, size=NUM_RECORDS).astype(np.int32)
    addresses = rng.integers(1_000_000, 9_999_999, size=NUM_RECORDS).astype(np.int32)

    # Trailing per-record float block (only meaningful for v9).
    trailing = rng.standard_normal(NUM_RECORDS).astype(np.float32)

    extra_per_rec = []  # 14 floats per record for v >= 7
    for _ in range(NUM_RECORDS):
        extra_per_rec.append(rng.standard_normal(14).astype(np.float32))

    with open(path, "wb") as f:
        write_int(f, 3)                 # version
        write_string(f, "GenCall-7.0")  # gencall_version
        write_string(f, "ClusterAlgo")  # cluster_version
        write_string(f, "CallAlgo")     # call_version
        write_string(f, "NormAlgo")     # normalization_version
        write_string(f, "1/1/2024")     # date_created
        write_byte(f, 1)                # is_wgt
        write_string(f, "TestManifest") # manifest_name
        write_int(f, version)           # data_block_version
        write_string(f, "")             # opa
        write_int(f, NUM_RECORDS)

        # Cluster records
        for i in range(NUM_RECORDS):
            # N counts
            write_int(f, int(aa_n[i]))
            write_int(f, int(ab_n[i]))
            write_int(f, int(bb_n[i]))
            # 12 floats: r_dev x3, r_mean x3, theta_dev x3, theta_mean x3
            for v in rng.standard_normal(12).astype(np.float32):
                write_float(f, float(v))
            if version >= 7:
                # intensity_threshold + 14 extra floats
                write_float(f, float(rng.standard_normal()))
                for v in extra_per_rec[i]:
                    write_float(f, float(v))

        # Cluster scores: 3 floats + 1 byte each
        for _ in range(NUM_RECORDS):
            write_float(f, float(rng.standard_normal()))
            write_float(f, float(rng.standard_normal()))
            write_float(f, float(rng.standard_normal()))
            write_byte(f, 0)

        # Genotype strings
        for _ in range(NUM_RECORDS):
            write_string(f, "AA/AB/BB")

        # Locus names
        for name in locus_names:
            write_string(f, name)

        # Addresses
        for a in addresses:
            write_int(f, int(a))

        # Cluster counts (mirrors N values)
        for i in range(NUM_RECORDS):
            write_int(f, int(aa_n[i]))
            write_int(f, int(ab_n[i]))
            write_int(f, int(bb_n[i]))

        # v9 only: trailing per-record float block
        if version == 9:
            for v in trailing:
                write_float(f, float(v))

    return trailing


with tempfile.TemporaryDirectory() as tmp:
    src = os.path.join(tmp, "src.egt")
    dst = os.path.join(tmp, "dst.egt")

    expected_trailing = build_synthetic_egt(src, data_block_version)
    src_size = os.path.getsize(src)

    egt = EGTFile.read(src)
    assert egt.data_block_version == data_block_version
    assert len(egt.records) == NUM_RECORDS

    if data_block_version == 9:
        assert len(egt.trailing_floats) == NUM_RECORDS, (
            f"Expected {NUM_RECORDS} trailing floats, got {len(egt.trailing_floats)}"
        )
        np.testing.assert_array_equal(
            np.asarray(egt.trailing_floats, dtype=np.float32),
            expected_trailing,
        )
    else:
        assert len(egt.trailing_floats) == 0, (
            f"v{data_block_version} should have no trailing floats, "
            f"got {len(egt.trailing_floats)}"
        )

    egt.write(dst)
    dst_size = os.path.getsize(dst)

    # The two files must be identical in length.  Specifically, when the
    # source is v9, the destination must NOT be num_records*4 bytes
    # short of the source — that shortfall is exactly the bug reported
    # upstream as "Failed to reposition stream forward <N> bytes".
    assert dst_size == src_size, (
        f"v{data_block_version}: rewritten EGT size {dst_size} != original {src_size} "
        f"(diff = {src_size - dst_size} bytes; expected 0)"
    )

    # Stronger: read the file back and confirm the trailing block
    # round-trips byte-exact.
    egt2 = EGTFile.read(dst)
    if data_block_version == 9:
        np.testing.assert_array_equal(
            np.asarray(egt2.trailing_floats, dtype=np.float32),
            expected_trailing,
        )

    print(f"  v{data_block_version}: src={src_size}B dst={dst_size}B ✓")
    sys.exit(0)
PYEOF
    check "data_block_version=${VERSION} round-trip preserves byte size" "$?"
    echo ""
done

# ---------------------------------------------------------------
# Targeted regression: simulate the issue's exact failure mode by
# stripping the trailing block from a v9 file and confirming the
# size shortfall is num_records*4, matching the reported error.
# ---------------------------------------------------------------
echo "--- Failure-mode reproduction: stripped v9 trailing block ---"
python3 - "${REPO_DIR}" <<'PYEOF'
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
from recluster_egt import EGTFile

# Reuse a tiny v9 EGT and verify that an unfixed writer would have
# produced a file num_records*4 bytes shorter.  We do this by writing
# with trailing_floats forcibly cleared, then comparing sizes.

with tempfile.TemporaryDirectory() as tmp:
    src = os.path.join(tmp, "src.egt")
    fixed = os.path.join(tmp, "fixed.egt")
    broken = os.path.join(tmp, "broken.egt")

    # Build src by re-using the helper above via a small inline write.
    # Easier: just call the round-trip path (we already trust it).
    # Build a v9 EGT inline (smaller copy of helper)
    import numpy as np
    from recluster_egt import write_byte, write_int, write_float, write_string

    NUM = 11
    rng = np.random.default_rng(0)
    with open(src, "wb") as f:
        write_int(f, 3)
        for s in ("g", "c", "ca", "n", "d"):
            write_string(f, s)
        write_byte(f, 1)
        write_string(f, "m")
        write_int(f, 9)
        write_string(f, "")
        write_int(f, NUM)
        for _ in range(NUM):
            write_int(f, 5); write_int(f, 5); write_int(f, 5)
            for v in rng.standard_normal(12).astype(np.float32):
                write_float(f, float(v))
            write_float(f, 0.0)
            for v in rng.standard_normal(14).astype(np.float32):
                write_float(f, float(v))
        for _ in range(NUM):
            write_float(f, 0.0); write_float(f, 0.0); write_float(f, 0.0); write_byte(f, 0)
        for _ in range(NUM): write_string(f, "")
        for i in range(NUM): write_string(f, f"r{i}")
        for i in range(NUM): write_int(f, i + 1)
        for _ in range(NUM):
            write_int(f, 5); write_int(f, 5); write_int(f, 5)
        for v in rng.standard_normal(NUM).astype(np.float32):
            write_float(f, float(v))

    src_size = os.path.getsize(src)
    egt = EGTFile.read(src)
    egt.write(fixed)
    fixed_size = os.path.getsize(fixed)
    assert fixed_size == src_size, f"Fixed writer size mismatch: {fixed_size} vs {src_size}"

    # Simulate the bug: blank out trailing_floats and write
    egt.trailing_floats = []
    # Temporarily monkey-patch the version check to skip writing the block
    orig_dbv = egt.data_block_version
    egt.data_block_version = 8  # fool the writer into omitting the block
    egt.write(broken)
    egt.data_block_version = orig_dbv
    broken_size = os.path.getsize(broken)

    shortfall = src_size - broken_size
    expected = NUM * 4
    assert shortfall == expected, (
        f"Expected shortfall of {expected} bytes (NUM*4), got {shortfall}"
    )
    print(f"  Confirmed: missing v9 trailing block produces "
          f"{shortfall}-byte shortfall (NUM={NUM}, NUM*4={expected})")
sys.exit(0)
PYEOF
check "v9 trailing-block omission produces NUM*4-byte shortfall" "$?"

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

exit "${FAIL}"
