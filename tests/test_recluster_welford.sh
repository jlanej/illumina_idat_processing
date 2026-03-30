#!/usr/bin/env bash
#
# test_recluster_welford.sh
#
# Unit tests for the streaming Welford's algorithm in recluster_egt.py.
# Verifies that _welford_update_pair produces mean and sample-std (ddof=1)
# identical to np.nanmean / np.nanstd for various edge cases.
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
echo "  Reclustering Welford Algorithm Tests"
echo "============================================"
echo ""

# ---------------------------------------------------------------
# Test 1: Welford vs batch NumPy — random data, multiple genotypes
# ---------------------------------------------------------------
echo "--- Test 1: Welford matches batch NumPy (random, multi-genotype) ---"

python3 - "${REPO_DIR}" <<'PYEOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
import numpy as np

from recluster_egt import _welford_update_pair

np.random.seed(42)
num_probes = 200
num_samples = 50

# Generate random data
all_theta = np.random.rand(num_probes, num_samples).astype(np.float32)
all_r = np.random.rand(num_probes, num_samples).astype(np.float32) * 2.0
all_geno = np.random.choice([0, 1, 2, 3], size=(num_probes, num_samples),
                             p=[0.05, 0.30, 0.35, 0.30]).astype(np.uint8)

# --- Batch reference (old approach) ---
ref_counts = {}
ref_theta_mean = {}
ref_theta_std = {}
ref_r_mean = {}
ref_r_std = {}

for gc in (1, 2, 3):
    mask = all_geno == gc
    ref_counts[gc] = mask.sum(axis=1)
    theta_masked = np.where(mask, all_theta, np.nan)
    r_masked = np.where(mask, all_r, np.nan)
    with np.errstate(all="ignore"):
        ref_theta_mean[gc] = np.nan_to_num(np.nanmean(theta_masked, axis=1), nan=0.0)
        ref_theta_std[gc] = np.nan_to_num(np.nanstd(theta_masked, axis=1, ddof=1), nan=0.0)
        ref_r_mean[gc] = np.nan_to_num(np.nanmean(r_masked, axis=1), nan=0.0)
        ref_r_std[gc] = np.nan_to_num(np.nanstd(r_masked, axis=1, ddof=1), nan=0.0)

# --- Streaming Welford ---
w_counts = {gc: np.zeros(num_probes, dtype=np.int64) for gc in (1, 2, 3)}
w_theta_mean = {gc: np.zeros(num_probes, dtype=np.float64) for gc in (1, 2, 3)}
w_theta_m2 = {gc: np.zeros(num_probes, dtype=np.float64) for gc in (1, 2, 3)}
w_r_mean = {gc: np.zeros(num_probes, dtype=np.float64) for gc in (1, 2, 3)}
w_r_m2 = {gc: np.zeros(num_probes, dtype=np.float64) for gc in (1, 2, 3)}

for sidx in range(num_samples):
    theta_col = all_theta[:, sidx]
    r_col = all_r[:, sidx]
    geno_col = all_geno[:, sidx]
    for gc in (1, 2, 3):
        mask = geno_col == gc
        _welford_update_pair(w_counts[gc], w_theta_mean[gc], w_theta_m2[gc],
                             w_r_mean[gc], w_r_m2[gc],
                             mask, theta_col, r_col)

# Finalize std
w_theta_std = {}
w_r_std = {}
for gc in (1, 2, 3):
    n = w_counts[gc].astype(np.float64)
    w_theta_std[gc] = np.where(n > 1, np.sqrt(w_theta_m2[gc] / (n - 1)), 0.0)
    w_r_std[gc] = np.where(n > 1, np.sqrt(w_r_m2[gc] / (n - 1)), 0.0)

# --- Compare ---
max_err = 0.0
for gc in (1, 2, 3):
    assert np.array_equal(ref_counts[gc], w_counts[gc]), f"Count mismatch for gc={gc}"
    for label, ref, wel in [
        ("theta_mean", ref_theta_mean[gc], w_theta_mean[gc]),
        ("theta_std", ref_theta_std[gc], w_theta_std[gc]),
        ("r_mean", ref_r_mean[gc], w_r_mean[gc]),
        ("r_std", ref_r_std[gc], w_r_std[gc]),
    ]:
        err = np.max(np.abs(ref - wel))
        max_err = max(max_err, err)
        assert err < 1e-5, f"gc={gc} {label}: max error {err:.2e} exceeds 1e-5"

print(f"  Max absolute error across all stats: {max_err:.2e}")
sys.exit(0)
PYEOF
check "Welford matches batch NumPy (200 probes, 50 samples)" "$?"

# ---------------------------------------------------------------
# Test 2: Edge case — single sample per genotype (n=1 → std=0)
# ---------------------------------------------------------------
echo ""
echo "--- Test 2: n=1 yields std=0 (Bessel undefined → 0.0) ---"

python3 - "${REPO_DIR}" <<'PYEOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
import numpy as np
from recluster_egt import _welford_update_pair

num_probes = 5
counts = np.zeros(num_probes, dtype=np.int64)
t_mean = np.zeros(num_probes, dtype=np.float64)
t_m2 = np.zeros(num_probes, dtype=np.float64)
r_mean = np.zeros(num_probes, dtype=np.float64)
r_m2 = np.zeros(num_probes, dtype=np.float64)

mask = np.ones(num_probes, dtype=bool)
theta = np.array([0.1, 0.5, 0.9, 0.3, 0.7], dtype=np.float32)
r = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
_welford_update_pair(counts, t_mean, t_m2, r_mean, r_m2, mask, theta, r)

assert np.all(counts == 1), f"Expected all counts=1, got {counts}"
assert np.allclose(t_mean, theta, atol=1e-7), f"Mean should equal value for n=1"

std = np.where(counts > 1, np.sqrt(t_m2 / (counts - 1)), 0.0)
assert np.all(std == 0.0), f"Expected std=0 for n=1, got {std}"

print("  n=1: counts, mean, std all correct")
sys.exit(0)
PYEOF
check "n=1 yields std=0" "$?"

# ---------------------------------------------------------------
# Test 3: Edge case — no samples (n=0 → mean=0, std=0)
# ---------------------------------------------------------------
echo ""
echo "--- Test 3: n=0 yields mean=0, std=0 ---"

python3 - "${REPO_DIR}" <<'PYEOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
import numpy as np
from recluster_egt import _welford_update_pair

num_probes = 3
counts = np.zeros(num_probes, dtype=np.int64)
t_mean = np.zeros(num_probes, dtype=np.float64)
t_m2 = np.zeros(num_probes, dtype=np.float64)
r_mean = np.zeros(num_probes, dtype=np.float64)
r_m2 = np.zeros(num_probes, dtype=np.float64)

mask = np.zeros(num_probes, dtype=bool)
theta = np.zeros(num_probes, dtype=np.float32)
r = np.zeros(num_probes, dtype=np.float32)
_welford_update_pair(counts, t_mean, t_m2, r_mean, r_m2, mask, theta, r)

assert np.all(counts == 0)
assert np.all(t_mean == 0.0)
assert np.all(t_m2 == 0.0)
std = np.where(counts > 1, np.sqrt(t_m2 / (counts - 1)), 0.0)
assert np.all(std == 0.0)

print("  n=0: counts, mean, std all zero")
sys.exit(0)
PYEOF
check "n=0 yields mean=0, std=0" "$?"

# ---------------------------------------------------------------
# Test 4: Two samples — verifies ddof=1 Bessel correction
# ---------------------------------------------------------------
echo ""
echo "--- Test 4: n=2 gives correct sample std (ddof=1) ---"

python3 - "${REPO_DIR}" <<'PYEOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
import numpy as np
from recluster_egt import _welford_update_pair

num_probes = 1
counts = np.zeros(num_probes, dtype=np.int64)
t_mean = np.zeros(num_probes, dtype=np.float64)
t_m2 = np.zeros(num_probes, dtype=np.float64)
r_mean = np.zeros(num_probes, dtype=np.float64)
r_m2 = np.zeros(num_probes, dtype=np.float64)

mask = np.ones(num_probes, dtype=bool)
_welford_update_pair(counts, t_mean, t_m2, r_mean, r_m2, mask,
                     np.array([3.0], dtype=np.float32),
                     np.array([10.0], dtype=np.float32))
_welford_update_pair(counts, t_mean, t_m2, r_mean, r_m2, mask,
                     np.array([5.0], dtype=np.float32),
                     np.array([20.0], dtype=np.float32))

assert counts[0] == 2
assert abs(t_mean[0] - 4.0) < 1e-10, f"Expected theta mean=4.0, got {t_mean[0]}"
assert abs(r_mean[0] - 15.0) < 1e-10, f"Expected r mean=15.0, got {r_mean[0]}"

t_std = np.sqrt(t_m2[0] / (counts[0] - 1))
r_std = np.sqrt(r_m2[0] / (counts[0] - 1))
expected_t_std = np.sqrt(2.0)  # std([3,5], ddof=1)
expected_r_std = np.sqrt(50.0)  # std([10,20], ddof=1)
assert abs(t_std - expected_t_std) < 1e-10
assert abs(r_std - expected_r_std) < 1e-10

print(f"  n=2: theta_mean={t_mean[0]}, theta_std={t_std:.6f}, r_mean={r_mean[0]}, r_std={r_std:.6f}")
sys.exit(0)
PYEOF
check "n=2 Bessel-corrected std" "$?"

# ---------------------------------------------------------------
# Test 5: Large-scale accuracy — 1M probes, 100 samples
# ---------------------------------------------------------------
echo ""
echo "--- Test 5: Large-scale accuracy (1M probes × 100 samples) ---"

python3 - "${REPO_DIR}" <<'PYEOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "scripts"))
import numpy as np
from recluster_egt import _welford_update_pair

np.random.seed(123)
num_probes = 1_000_000
num_samples = 100
geno_code = 2

geno_all = np.random.choice([0, 1, 2, 3], size=(num_probes, num_samples),
                             p=[0.1, 0.3, 0.3, 0.3]).astype(np.uint8)
theta_all = np.random.rand(num_probes, num_samples).astype(np.float32)
r_all = np.random.rand(num_probes, num_samples).astype(np.float32) * 2.0

# Batch reference
mask_all = geno_all == geno_code
ref_counts = mask_all.sum(axis=1)
theta_masked = np.where(mask_all, theta_all, np.nan)
r_masked = np.where(mask_all, r_all, np.nan)
with np.errstate(all="ignore"):
    ref_t_mean = np.nan_to_num(np.nanmean(theta_masked, axis=1), nan=0.0)
    ref_t_std = np.nan_to_num(np.nanstd(theta_masked, axis=1, ddof=1), nan=0.0)
    ref_r_mean = np.nan_to_num(np.nanmean(r_masked, axis=1), nan=0.0)
    ref_r_std = np.nan_to_num(np.nanstd(r_masked, axis=1, ddof=1), nan=0.0)

# Streaming Welford
w_counts = np.zeros(num_probes, dtype=np.int64)
w_t_mean = np.zeros(num_probes, dtype=np.float64)
w_t_m2 = np.zeros(num_probes, dtype=np.float64)
w_r_mean = np.zeros(num_probes, dtype=np.float64)
w_r_m2 = np.zeros(num_probes, dtype=np.float64)

for sidx in range(num_samples):
    mask = geno_all[:, sidx] == geno_code
    _welford_update_pair(w_counts, w_t_mean, w_t_m2, w_r_mean, w_r_m2,
                         mask, theta_all[:, sidx], r_all[:, sidx])

w_t_std = np.where(w_counts > 1, np.sqrt(w_t_m2 / (w_counts - 1)), 0.0)
w_r_std = np.where(w_counts > 1, np.sqrt(w_r_m2 / (w_counts - 1)), 0.0)

assert np.array_equal(ref_counts, w_counts), "Count mismatch"
t_mean_err = np.max(np.abs(ref_t_mean - w_t_mean))
t_std_err = np.max(np.abs(ref_t_std - w_t_std))
r_mean_err = np.max(np.abs(ref_r_mean - w_r_mean))
r_std_err = np.max(np.abs(ref_r_std - w_r_std))

print(f"  1M probes × 100 samples: theta_mean_err={t_mean_err:.2e}, "
      f"theta_std_err={t_std_err:.2e}, r_mean_err={r_mean_err:.2e}, "
      f"r_std_err={r_std_err:.2e}")
assert t_mean_err < 1e-5, f"Theta mean error {t_mean_err:.2e} too large"
assert t_std_err < 1e-5, f"Theta std error {t_std_err:.2e} too large"
assert r_mean_err < 1e-5, f"R mean error {r_mean_err:.2e} too large"
assert r_std_err < 1e-5, f"R std error {r_std_err:.2e} too large"
sys.exit(0)
PYEOF
check "Large-scale accuracy (1M probes × 100 samples)" "$?"

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

exit "${FAIL}"
