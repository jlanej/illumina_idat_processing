#!/usr/bin/env python3
"""
compute_sample_qc.py

Compute per-sample QC metrics from a VCF/BCF stream in a SINGLE pass.

Reads bcftools query output (GT:LRR:BAF per sample per variant) from stdin
and computes:
  - call_rate: fraction of non-missing genotypes (autosomes only)
  - lrr_sd: standard deviation of LRR (autosomes only)
  - lrr_mean: mean of LRR (autosomes only)
  - lrr_median: median of LRR (autosomes only) via numpy.partition (O(n))
  - baf_sd: standard deviation of BAF for heterozygous SNPs (contamination metric)
  - het_rate: fraction of heterozygous genotype calls

This replaces the previous two-pass approach (awk for call_rate/lrr_sd + Python
for lrr_median), halving the BCF decompression I/O.

Usage (called by collect_qc_metrics.sh):
    bcftools view -e 'INFO/INTENSITY_ONLY=1' -t ^chrX,chrY,chrM,X,Y,MT file.bcf |
    bcftools query -f '[\\t%GT:%LRR:%BAF]\\n' |
    python3 compute_sample_qc.py --num-samples N --output metrics.tsv

Performance:
    - Single pass over the BCF stream
    - O(n) median via numpy.partition (no full sort)
    - Streaming mean/variance via Welford's algorithm (constant memory per sample
      for call rate, LRR SD, BAF SD, het rate)
    - LRR values stored for exact median (bounded by #autosomal variants)
"""

import argparse
import sys

import numpy as np


_INVALID_LRR = frozenset({'.', '', 'nan', '-nan', 'inf', '-inf'})
_HET_GT = frozenset({'0/1', '1/0', '0|1', '1|0'})
_MISSING_GT = frozenset({'./.', '.', '.|.'})


def compute_metrics(num_samples, output_path):
    """Read GT:LRR:BAF triples from stdin, compute per-sample QC metrics."""

    # Per-sample accumulators
    total = np.zeros(num_samples, dtype=np.int64)
    called = np.zeros(num_samples, dtype=np.int64)
    het_count = np.zeros(num_samples, dtype=np.int64)

    # LRR: Welford's online algorithm for mean/variance
    lrr_n = np.zeros(num_samples, dtype=np.int64)
    lrr_mean = np.zeros(num_samples, dtype=np.float64)
    lrr_m2 = np.zeros(num_samples, dtype=np.float64)

    # BAF SD for het sites: Welford's online algorithm
    baf_n = np.zeros(num_samples, dtype=np.int64)
    baf_mean = np.zeros(num_samples, dtype=np.float64)
    baf_m2 = np.zeros(num_samples, dtype=np.float64)

    # Collect LRR values for median (use list of lists, convert to numpy at end)
    lrr_values = [[] for _ in range(num_samples)]

    line_count = 0
    for line in sys.stdin:
        line_count += 1
        fields = line.rstrip('\n').split('\t')
        # fields[0] is empty (leading tab from bcftools query format)
        samples = fields[1:]

        if len(samples) != num_samples:
            if line_count == 1:
                # Auto-detect num_samples if it doesn't match
                num_samples = len(samples)
                total = np.zeros(num_samples, dtype=np.int64)
                called = np.zeros(num_samples, dtype=np.int64)
                het_count = np.zeros(num_samples, dtype=np.int64)
                lrr_n = np.zeros(num_samples, dtype=np.int64)
                lrr_mean = np.zeros(num_samples, dtype=np.float64)
                lrr_m2 = np.zeros(num_samples, dtype=np.float64)
                baf_n = np.zeros(num_samples, dtype=np.int64)
                baf_mean = np.zeros(num_samples, dtype=np.float64)
                baf_m2 = np.zeros(num_samples, dtype=np.float64)
                lrr_values = [[] for _ in range(num_samples)]
            else:
                continue

        for i, sample_field in enumerate(samples):
            parts = sample_field.split(':')
            gt = parts[0] if len(parts) > 0 else '.'
            lrr = parts[1] if len(parts) > 1 else '.'
            baf = parts[2] if len(parts) > 2 else '.'

            # Call rate
            total[i] += 1
            is_called = gt not in _MISSING_GT
            if is_called:
                called[i] += 1

            # Heterozygosity
            is_het = gt in _HET_GT
            if is_het:
                het_count[i] += 1

            # LRR statistics (Welford's online algorithm)
            if lrr.lower() not in _INVALID_LRR:
                try:
                    lrr_val = float(lrr)
                    lrr_n[i] += 1
                    delta = lrr_val - lrr_mean[i]
                    lrr_mean[i] += delta / lrr_n[i]
                    delta2 = lrr_val - lrr_mean[i]
                    lrr_m2[i] += delta * delta2
                    lrr_values[i].append(lrr_val)
                except ValueError:
                    pass

            # BAF SD for het sites only
            if is_het and baf.lower() not in _INVALID_LRR:
                try:
                    baf_val = float(baf)
                    baf_n[i] += 1
                    delta = baf_val - baf_mean[i]
                    baf_mean[i] += delta / baf_n[i]
                    delta2 = baf_val - baf_mean[i]
                    baf_m2[i] += delta * delta2
                except ValueError:
                    pass

        if line_count % 500000 == 0:
            print(f"  Processed {line_count} variants...", file=sys.stderr)

    print(f"  Processed {line_count} total variants across {num_samples} samples",
          file=sys.stderr)

    # Compute final statistics
    with open(output_path, 'w') as out:
        out.write("idx\tcall_rate\tlrr_sd\tlrr_mean\tlrr_median\tbaf_sd\thet_rate\n")

        for i in range(num_samples):
            # Call rate
            cr = called[i] / total[i] if total[i] > 0 else 0.0

            # LRR SD (population variance, matching original implementation)
            if lrr_n[i] > 1:
                variance = lrr_m2[i] / lrr_n[i]
                if variance < 0:
                    variance = 0
                lrr_sd = variance ** 0.5
                mean_val = lrr_mean[i]
            else:
                lrr_sd = float('nan')
                mean_val = float('nan')

            # LRR median via numpy.partition (O(n) selection, no full sort)
            if lrr_values[i]:
                arr = np.array(lrr_values[i], dtype=np.float32)
                n = len(arr)
                if n % 2 == 1:
                    k = n // 2
                    np.partition(arr, k)
                    median_val = float(arr[k])
                else:
                    k1, k2 = n // 2 - 1, n // 2
                    np.partition(arr, k2)
                    # After partition around k2, arr[:k2] are all <= arr[k2]
                    # The k1-th element is the max of arr[:k2]
                    median_val = float((np.max(arr[:k2]) + arr[k2]) / 2)
            else:
                median_val = float('nan')

            # BAF SD (population variance)
            if baf_n[i] > 1:
                baf_var = baf_m2[i] / baf_n[i]
                if baf_var < 0:
                    baf_var = 0
                baf_sd_val = baf_var ** 0.5
            else:
                baf_sd_val = float('nan')

            # Heterozygosity rate
            het_rate = het_count[i] / called[i] if called[i] > 0 else 0.0

            # Format values
            def fmt(v, decimals=6):
                if v != v:  # NaN check
                    return "NA"
                return f"{v:.{decimals}f}"

            out.write(f"{i}\t{fmt(cr, 6)}\t{fmt(lrr_sd, 6)}\t{fmt(mean_val, 6)}"
                      f"\t{fmt(median_val, 6)}\t{fmt(baf_sd_val, 6)}\t{fmt(het_rate, 6)}\n")

    print(f"  Metrics written to {output_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Compute per-sample QC metrics in a single BCF pass"
    )
    parser.add_argument("--num-samples", type=int, default=0,
                        help="Number of samples (auto-detected if 0)")
    parser.add_argument("--output", required=True,
                        help="Output TSV file for metrics")
    args = parser.parse_args()

    compute_metrics(args.num_samples, args.output)


if __name__ == "__main__":
    main()

