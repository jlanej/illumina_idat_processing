#!/usr/bin/env python3
"""
generate_mock_data.py

Generate deterministic mock pipeline output for integration testing and
GitHub Pages demo report.  All random values use a fixed seed so output
is reproducible across runs.

Usage:
    python3 tests/generate_mock_data.py --output-dir /tmp/mock_output
"""

import argparse
import math
import os
import random
import sys

# Fixed seed for reproducibility
SEED = 42


def _rng():
    """Return a seeded RNG."""
    return random.Random(SEED)


def generate_sample_qc(output_dir, stage, n_samples=30):
    """Write a stage sample_qc.tsv with deterministic mock data.

    The mock cohort is designed to exercise every QC feature:
      - Most samples pass all thresholds
      - 2-3 samples fail call rate (< 0.97)
      - 2-3 samples fail LRR SD (> 0.35)
      - 1-2 samples flagged for BAF SD (> 0.15)
      - 1-2 samples are het rate outliers (very high/low het)
      - Mix of M/F predicted sex
    """
    rng = _rng()

    qc_dir = os.path.join(output_dir, stage, 'qc')
    os.makedirs(qc_dir, exist_ok=True)

    rows = []
    for i in range(n_samples):
        sid = f"SAMPLE_{i + 1:03d}"

        # --- Call Rate ---
        if i == 3:
            cr = 0.9521  # fail: below 0.97
        elif i == 12:
            cr = 0.9388  # fail: below 0.97
        elif i == 22:
            cr = 0.9655  # fail: below 0.97
        else:
            cr = round(rng.uniform(0.9700, 0.9999), 4)

        # --- LRR SD ---
        if i == 5:
            lsd = 0.4123  # fail: above 0.35
        elif i == 14:
            lsd = 0.3891  # fail: above 0.35
        elif i == 24:
            lsd = 0.3601  # fail: above 0.35
        else:
            lsd = round(rng.uniform(0.08, 0.30), 4)

        # --- LRR Mean / Median ---
        lrr_mean = round(rng.uniform(-0.02, 0.02), 4)
        lrr_median = round(lrr_mean + rng.uniform(-0.005, 0.005), 4)

        # --- BAF SD (het sites) ---
        if i == 7:
            bsd = 0.1823  # flag: above 0.15
        elif i == 18:
            bsd = 0.1645  # flag: above 0.15
        else:
            bsd = round(rng.uniform(0.02, 0.09), 4)

        # --- Het Rate ---
        if i == 9:
            hr = 0.3412  # very high → outlier
        elif i == 20:
            hr = 0.1534  # very low → outlier
        else:
            hr = round(rng.uniform(0.22, 0.28), 4)

        # --- Gender ---
        gender = 'M' if i % 2 == 0 else 'F'

        rows.append({
            'sample_id': sid,
            'call_rate': f"{cr:.4f}",
            'lrr_sd': f"{lsd:.4f}",
            'lrr_mean': f"{lrr_mean:.4f}",
            'lrr_median': f"{lrr_median:.4f}",
            'baf_sd': f"{bsd:.4f}",
            'het_rate': f"{hr:.4f}",
            'computed_gender': gender,
        })

    filepath = os.path.join(qc_dir, f'{stage}_sample_qc.tsv')
    cols = ['sample_id', 'call_rate', 'lrr_sd', 'lrr_mean',
            'lrr_median', 'baf_sd', 'het_rate', 'computed_gender']
    with open(filepath, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(r[c] for c in cols) + '\n')

    return filepath


def generate_variant_qc_summary(output_dir, stage):
    """Write a mock variant QC summary."""
    vqc_dir = os.path.join(output_dir, stage, 'qc', 'variant_qc')
    os.makedirs(vqc_dir, exist_ok=True)

    text = """==================================================
  Variant QC Summary (autosomes)
==================================================

  Input variants:         650,427
  After call-rate filter: 642,185 (≥ 0.98)
  After HWE filter:       639,412 (p ≥ 1e-6)
  After MAF filter:       531,208 (≥ 0.01)

  Ti/Tv ratio:            2.07
  Mean inbreeding F:      0.0012
  Inbreeding F range:     -0.032 to 0.041

  Missingness histogram:
    0.00-0.01:  ████████████████████████  598,127
    0.01-0.02:  ████                       44,058
    0.02-0.05:  ██                          6,923
    0.05-0.10:  ▏                           1,062
    0.10-1.00:  ▏                             257

  MAF distribution:
    0.00-0.01:  ████████████████           108,977
    0.01-0.05:  ████████████████████       132,245
    0.05-0.10:  ██████████████             98,456
    0.10-0.25:  ████████████████████████  164,123
    0.25-0.50:  ████████████████████████  146,626
"""
    filepath = os.path.join(vqc_dir, 'variant_qc_summary.txt')
    with open(filepath, 'w') as f:
        f.write(text)
    return filepath


def generate_comparison_report(output_dir):
    """Write a mock QC comparison report for stage1→stage2."""
    comp_dir = os.path.join(output_dir, 'stage2', 'qc', 'comparison')
    os.makedirs(comp_dir, exist_ok=True)

    text = """==================================================
  Stage 1 → Stage 2 QC Comparison
==================================================

  Metric              Stage 1     Stage 2     Change
  ─────────────────── ─────────── ─────────── ──────────
  Mean Call Rate       0.9841      0.9893      +0.0052
  Mean LRR SD          0.2312      0.2045      -0.0267
  Mean BAF SD          0.0612      0.0534      -0.0078
  Mean Het Rate        0.2489      0.2478      -0.0011

  Variant Missingness  0.0134      0.0098      -0.0036
  HWE failures (1e-6)  1,245       823         -422
  Mean MAF             0.1723      0.1734      +0.0011
  Ti/Tv ratio          2.05        2.07        +0.02

  Reclustering improved call rate for 28/30 samples
  and reduced LRR SD for 26/30 samples.
"""
    filepath = os.path.join(comp_dir, 'qc_comparison_report.txt')
    with open(filepath, 'w') as f:
        f.write(text)
    return filepath


def generate_realign_summary(output_dir):
    """Write a mock manifest realignment summary."""
    realign_dir = os.path.join(output_dir, 'realigned_manifests')
    os.makedirs(realign_dir, exist_ok=True)

    text = """==================================================
  Manifest Realignment Summary
==================================================

  Manifest:         GSA-24v3-0_A1
  Reference:        chm13v2.0 (T2T-CHM13v2.0)
  Aligner:          bwa 0.7.18

  Total probes:     654,027
  Mapped:           650,427 (99.45%)
  Unmapped:         2,134 (0.33%)
  Multi-mapped:     1,466 (0.22%)

  Position changes: 3,812 (0.58%)
    Same chromosome, different position:  3,701
    Different chromosome:                    111

  Strand changes:   245
  New placements:   89 (previously unmapped)
"""
    filepath = os.path.join(realign_dir, 'realign_summary.txt')
    with open(filepath, 'w') as f:
        f.write(text)
    return filepath


def generate_diagnostic_report(output_dir):
    """Write a mock QC diagnostic report."""
    text = """==================================================
  QC Diagnostic Report
==================================================

  ✓ Call rate distribution appears normal
  ✓ LRR SD distribution is within expected range
  ⚠ 2 samples flagged for BAF SD > 0.15 (possible contamination)
  ⚠ 2 samples are heterozygosity rate outliers (± 3 SD)
  ✓ Gender predictions consistent (15 M, 15 F)
  ✓ Ti/Tv ratio 2.07 is within expected range (1.8-2.2)
  ✓ No build mismatch detected between manifest and reference

  Recommendations:
  - Review BAF SD-flagged samples (SAMPLE_008, SAMPLE_019) for
    possible contamination or mosaicism
  - Review het rate outliers (SAMPLE_010, SAMPLE_021) for
    population stratification or sample swaps
"""
    filepath = os.path.join(output_dir, 'qc_diagnostic_report.txt')
    with open(filepath, 'w') as f:
        f.write(text)
    return filepath


def generate_pca_data(output_dir, n_samples=30):
    """Write mock PCA projections and QC summary."""
    pca_dir = os.path.join(output_dir, 'ancestry_pca')
    os.makedirs(pca_dir, exist_ok=True)

    rng = _rng()

    # PCA projections (20 PCs)
    header = ['FID', 'IID'] + [f'PC{i}' for i in range(1, 21)]
    filepath = os.path.join(pca_dir, 'pca_projections.tsv')
    with open(filepath, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for i in range(n_samples):
            sid = f"SAMPLE_{i + 1:03d}"
            # Create 3 rough clusters
            cluster = i % 3
            base_pc1 = [-0.02, 0.01, 0.03][cluster]
            base_pc2 = [0.01, -0.02, 0.02][cluster]
            pcs = [round(base_pc1 + rng.gauss(0, 0.005), 6),
                   round(base_pc2 + rng.gauss(0, 0.005), 6)]
            for _ in range(18):
                pcs.append(round(rng.gauss(0, 0.003), 6))
            row = [sid, sid] + [str(p) for p in pcs]
            f.write('\t'.join(row) + '\n')

    # QC summary
    qc_text = """==================================================
  PCA QC Summary
==================================================

  Input samples:          30
  After call rate ≥ 0.98: 28
  After mind ≤ 0.02:      28

  Input variants:         650,427
  After geno ≤ 0.02:      642,185
  After HWE (1e-6):       639,412
  After MAF ≥ 0.05:       485,234
  After LD pruning:       98,412

  PCs computed: 20
  Eigenvalues: 12.34, 8.56, 5.23, 3.89, 2.67, ...
"""
    with open(os.path.join(pca_dir, 'qc_summary.txt'), 'w') as f:
        f.write(qc_text)

    return filepath


def main():
    parser = argparse.ArgumentParser(
        description="Generate deterministic mock pipeline output for testing"
    )
    parser.add_argument("--output-dir", required=True,
                        help="Directory to write mock output")
    parser.add_argument("--n-samples", type=int, default=30,
                        help="Number of mock samples (default: 30)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    n = args.n_samples

    print(f"Generating mock pipeline output in: {args.output_dir}")
    print(f"  Samples: {n}")

    # Stage 1 QC (same structure, deterministic from fixed seed)
    s1 = generate_sample_qc(args.output_dir, 'stage1', n)
    print(f"  Stage 1 sample QC: {s1}")

    # Stage 2 QC (same samples, slightly better metrics)
    s2 = generate_sample_qc(args.output_dir, 'stage2', n)
    print(f"  Stage 2 sample QC: {s2}")

    # Variant QC
    vqc = generate_variant_qc_summary(args.output_dir, 'stage2')
    print(f"  Variant QC:        {vqc}")

    # Comparison report
    comp = generate_comparison_report(args.output_dir)
    print(f"  Comparison report: {comp}")

    # Realignment summary
    realign = generate_realign_summary(args.output_dir)
    print(f"  Realignment:       {realign}")

    # Diagnostic report
    diag = generate_diagnostic_report(args.output_dir)
    print(f"  Diagnostics:       {diag}")

    # PCA data
    pca = generate_pca_data(args.output_dir, n)
    print(f"  PCA projections:   {pca}")

    print("\nDone. Run generate_report.py --output-dir to build the HTML report.")


if __name__ == "__main__":
    main()
