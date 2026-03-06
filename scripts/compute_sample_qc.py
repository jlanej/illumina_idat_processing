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
    - Chunked numpy-vectorized processing (~5-10× faster than per-sample loop)
    - O(n) median via numpy.partition (no full sort)
    - Streaming sum/sum-of-squares for mean/variance (constant memory per sample
      for call rate, LRR SD, BAF SD, het rate)
    - LRR values stored for exact median (bounded by #autosomal variants)
"""

import argparse
import sys

import numpy as np


_INVALID_LRR = frozenset({'.', '', 'nan', '-nan', 'inf', '-inf'})
_HET_GT = frozenset({'0/1', '1/0', '0|1', '1|0'})
_MISSING_GT = frozenset({'./.', '.', '.|.'})


def compute_metrics(num_samples, output_path, chunk_size=10000):
    """Read GT:LRR:BAF triples from stdin, compute per-sample QC metrics.

    Uses chunked numpy-vectorized processing for ~5-10× speedup over the
    previous pure-Python per-sample inner loop.  Lines are accumulated into
    chunks of *chunk_size* variants, then parsed and processed with
    vectorized numpy operations.
    """

    # Per-sample accumulators (initialised lazily after first chunk)
    total = None
    called = None
    het_count = None

    lrr_n = None
    lrr_sum = None
    lrr_sum_sq = None

    baf_n = None
    baf_sum = None
    baf_sum_sq = None

    # Collect LRR values for median (per-sample list of float32 chunks)
    lrr_values = None

    def _init_accumulators(ns):
        """(Re-)initialise all accumulators for *ns* samples."""
        nonlocal total, called, het_count
        nonlocal lrr_n, lrr_sum, lrr_sum_sq
        nonlocal baf_n, baf_sum, baf_sum_sq
        nonlocal lrr_values, num_samples

        num_samples = ns
        total = np.zeros(ns, dtype=np.int64)
        called = np.zeros(ns, dtype=np.int64)
        het_count = np.zeros(ns, dtype=np.int64)
        lrr_n = np.zeros(ns, dtype=np.int64)
        lrr_sum = np.zeros(ns, dtype=np.float64)
        lrr_sum_sq = np.zeros(ns, dtype=np.float64)
        baf_n = np.zeros(ns, dtype=np.int64)
        baf_sum = np.zeros(ns, dtype=np.float64)
        baf_sum_sq = np.zeros(ns, dtype=np.float64)
        lrr_values = [[] for _ in range(ns)]

    def _process_chunk(lines):
        """Parse and accumulate a chunk of variant lines (vectorised)."""
        nonlocal num_samples

        n_lines = len(lines)
        if n_lines == 0:
            return

        # --- Parse all GT:LRR:BAF fields into parallel arrays ----------
        # Each line: \tGT:LRR:BAF\tGT:LRR:BAF\t...
        gt_arr = np.empty((n_lines, num_samples), dtype=object)
        lrr_arr = np.full((n_lines, num_samples), np.nan, dtype=np.float64)
        baf_arr = np.full((n_lines, num_samples), np.nan, dtype=np.float64)

        for row_idx, line in enumerate(lines):
            fields = line.rstrip('\n').split('\t')
            samples = fields[1:]  # skip leading empty field
            if len(samples) != num_samples:
                gt_arr[row_idx, :] = '.'
                continue
            for col_idx, sf in enumerate(samples):
                parts = sf.split(':')
                gt_arr[row_idx, col_idx] = parts[0] if parts else '.'
                if len(parts) > 1:
                    try:
                        lrr_arr[row_idx, col_idx] = float(parts[1])
                    except (ValueError, IndexError):
                        pass
                if len(parts) > 2:
                    try:
                        baf_arr[row_idx, col_idx] = float(parts[2])
                    except (ValueError, IndexError):
                        pass

        # --- Vectorised call-rate / het counting ----------------------
        is_missing = np.isin(gt_arr, list(_MISSING_GT))
        is_het = np.isin(gt_arr, list(_HET_GT))
        total[:] += n_lines
        called[:] += np.sum(~is_missing, axis=0)
        het_count[:] += np.sum(is_het, axis=0)

        # --- Vectorised LRR accumulation ------------------------------
        lrr_valid = np.isfinite(lrr_arr)
        lrr_n[:] += np.sum(lrr_valid, axis=0)
        lrr_clean = np.where(lrr_valid, lrr_arr, 0.0)
        lrr_sum[:] += np.sum(lrr_clean, axis=0)
        lrr_sum_sq[:] += np.sum(lrr_clean ** 2, axis=0)

        # Collect LRR values for median (keep compact float32 chunk copies)
        for col_idx in range(num_samples):
            col_valid = lrr_valid[:, col_idx]
            if np.any(col_valid):
                lrr_values[col_idx].append(
                    lrr_arr[col_valid, col_idx].astype(np.float32, copy=True)
                )

        # --- Vectorised BAF accumulation (het-only) -------------------
        baf_valid = is_het & np.isfinite(baf_arr)
        baf_n[:] += np.sum(baf_valid, axis=0)
        baf_clean = np.where(baf_valid, baf_arr, 0.0)
        baf_sum[:] += np.sum(baf_clean, axis=0)
        baf_sum_sq[:] += np.sum(baf_clean ** 2, axis=0)

    # --- Main streaming loop ------------------------------------------
    line_count = 0
    chunk_lines = []

    for line in sys.stdin:
        line_count += 1
        chunk_lines.append(line)

        # Auto-detect / validate num_samples on first line
        if line_count == 1:
            first_fields = line.rstrip('\n').split('\t')
            detected = len(first_fields) - 1  # skip leading empty field
            if num_samples <= 0 or detected != num_samples:
                num_samples = detected
            _init_accumulators(num_samples)

        if len(chunk_lines) >= chunk_size:
            _process_chunk(chunk_lines)
            chunk_lines = []
            if line_count % 500000 < chunk_size:
                print(f"  Processed {line_count} variants...", file=sys.stderr)

    # Flush remaining lines
    if chunk_lines:
        _process_chunk(chunk_lines)

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
                mean_val = lrr_sum[i] / lrr_n[i]
                variance = lrr_sum_sq[i] / lrr_n[i] - mean_val ** 2
                if variance < 0:
                    variance = 0
                lrr_sd = variance ** 0.5
            else:
                lrr_sd = float('nan')
                mean_val = lrr_sum[i] / lrr_n[i] if lrr_n[i] > 0 else float('nan')

            # LRR median via in-place numpy partition (O(n) selection)
            if lrr_values[i]:
                if len(lrr_values[i]) == 1:
                    arr = lrr_values[i][0].copy()
                else:
                    arr = np.concatenate(lrr_values[i]).astype(np.float32, copy=False)
                n = arr.size
                if n % 2 == 1:
                    k = n // 2
                    arr.partition(k)
                    median_val = float(arr[k])
                else:
                    k2 = n // 2
                    k1 = k2 - 1
                    arr.partition((k1, k2))
                    median_val = float((arr[k1] + arr[k2]) / 2.0)
            else:
                median_val = float('nan')

            # BAF SD (population variance)
            if baf_n[i] > 1:
                baf_mean_val = baf_sum[i] / baf_n[i]
                baf_var = baf_sum_sq[i] / baf_n[i] - baf_mean_val ** 2
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
