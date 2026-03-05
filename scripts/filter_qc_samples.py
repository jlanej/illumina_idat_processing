#!/usr/bin/env python3
"""
filter_qc_samples.py

Filter samples from a QC TSV file based on call rate and LRR SD thresholds.
Produces lists of high-quality and excluded samples for Stage 2 reclustering.

Uses header-based column lookup so the script works regardless of whether
optional columns (e.g. original_name) are present.

Usage:
    python3 filter_qc_samples.py \
        --qc-file stage1/qc/stage1_sample_qc.tsv \
        --hq-output high_quality_samples.txt \
        --excluded-output excluded_samples.txt \
        --min-call-rate 0.97 \
        --max-lrr-sd 0.35
"""

import argparse
import sys


def filter_samples(qc_file, hq_output, excluded_output, min_call_rate, max_lrr_sd):
    """Read QC TSV, split samples into high-quality and excluded lists."""

    hq_samples = []
    excluded_samples = []
    excluded_details = []

    with open(qc_file) as f:
        header_line = f.readline().rstrip('\n')
        columns = header_line.split('\t')

        # Build column index map from header names
        col_idx = {name: i for i, name in enumerate(columns)}

        required = ['sample_id', 'call_rate', 'lrr_sd']
        missing = [c for c in required if c not in col_idx]
        if missing:
            print(f"Error: QC file missing required columns: {missing}",
                  file=sys.stderr)
            print(f"  Header: {header_line}", file=sys.stderr)
            print(f"  Available columns: {columns}", file=sys.stderr)
            sys.exit(1)

        sid_idx = col_idx['sample_id']
        cr_idx = col_idx['call_rate']
        sd_idx = col_idx['lrr_sd']

        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= max(sid_idx, cr_idx, sd_idx):
                continue

            sample_id = fields[sid_idx]
            call_rate_str = fields[cr_idx]
            lrr_sd_str = fields[sd_idx]

            # Parse call rate
            try:
                call_rate = float(call_rate_str)
            except (ValueError, TypeError):
                # Non-numeric call_rate → exclude
                excluded_samples.append(sample_id)
                excluded_details.append(
                    (sample_id, call_rate_str, lrr_sd_str,
                     'FAIL(parse)', 'FAIL(parse)'))
                continue

            # Parse LRR SD
            if lrr_sd_str in ('NA', '', 'nan', 'NaN'):
                lrr_sd = None
            else:
                try:
                    lrr_sd = float(lrr_sd_str)
                except (ValueError, TypeError):
                    excluded_samples.append(sample_id)
                    excluded_details.append(
                        (sample_id, call_rate_str, lrr_sd_str,
                         'ok', 'FAIL(parse)'))
                    continue

            # Apply thresholds
            cr_pass = call_rate >= min_call_rate
            sd_pass = lrr_sd is None or lrr_sd <= max_lrr_sd

            if cr_pass and sd_pass:
                hq_samples.append(sample_id)
            else:
                excluded_samples.append(sample_id)
                cr_status = 'ok' if cr_pass else 'FAIL'
                sd_status = 'ok' if sd_pass else 'FAIL'
                excluded_details.append(
                    (sample_id, call_rate_str, lrr_sd_str,
                     cr_status, sd_status))

    # Write output files
    with open(hq_output, 'w') as f:
        for s in hq_samples:
            f.write(s + '\n')

    with open(excluded_output, 'w') as f:
        for s in excluded_samples:
            f.write(s + '\n')

    # Print summary
    n_total = len(hq_samples) + len(excluded_samples)
    print(f"  Total samples:            {n_total}")
    print(f"  High-quality samples:     {len(hq_samples)}")
    print(f"  Excluded samples:         {len(excluded_samples)}")
    print(f"  [diag] QC thresholds: call_rate >= {min_call_rate}, "
          f"lrr_sd <= {max_lrr_sd}")

    # Show first 5 lines of QC file
    with open(qc_file) as f:
        lines = [f.readline().rstrip('\n') for _ in range(5)]
    print("  [diag] Stage 1 QC file first 5 lines:")
    for line in lines:
        print(f"    [diag]   {line}")

    # Show first 5 excluded samples with reasons
    if excluded_details:
        print("  [diag] First 5 excluded samples with values:")
        for sid, cr, sd, cr_s, sd_s in excluded_details[:5]:
            print(f"    [diag]   {sid}: call_rate={cr}({cr_s}) "
                  f"lrr_sd={sd}({sd_s})")


def main():
    parser = argparse.ArgumentParser(
        description="Filter samples by QC thresholds for Stage 2 reclustering"
    )
    parser.add_argument("--qc-file", required=True,
                        help="Stage 1 sample QC TSV file")
    parser.add_argument("--hq-output", required=True,
                        help="Output file for high-quality sample IDs")
    parser.add_argument("--excluded-output", required=True,
                        help="Output file for excluded sample IDs")
    parser.add_argument("--min-call-rate", type=float, default=0.97,
                        help="Minimum call rate threshold (default: 0.97)")
    parser.add_argument("--max-lrr-sd", type=float, default=0.35,
                        help="Maximum LRR SD threshold (default: 0.35)")
    args = parser.parse_args()

    filter_samples(args.qc_file, args.hq_output, args.excluded_output,
                   args.min_call_rate, args.max_lrr_sd)


if __name__ == "__main__":
    main()

