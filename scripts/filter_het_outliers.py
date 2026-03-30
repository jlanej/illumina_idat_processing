#!/usr/bin/env python3
"""
filter_het_outliers.py

Identify and exclude samples whose heterozygosity rate is more than a
configurable number of standard deviations (default 3) from the mean,
computed **within** the purported ancestry group.

Ancestry assignments are read from peddy's het_check.csv
(``ancestry-prediction`` column).  Sample heterozygosity rates are read
from the pipeline's sample QC TSV (``het_rate`` column).

Samples that fall outside the per-ancestry mean ± N SD are written to
an exclusion list.  These samples are removed prior to PCA, HWE, and
ancestry-stratified QC to avoid biased population structure estimation
and inflated HWE deviation.

Usage:
    python3 filter_het_outliers.py \
        --sample-qc stage2_sample_qc.tsv \
        --peddy-het-check peddy.het_check.csv \
        --output-excluded het_outlier_excluded.txt \
        --output-log het_outlier_details.tsv \
        --sd-threshold 3

References:
    - Anderson et al. 2010, Nat Protoc 5:1564-73
    - Marees et al. 2018, Int J Methods Psychiatr Res 27:e1608
"""

import argparse
import csv
import statistics
import sys


def load_het_rates(qc_file):
    """Read het_rate from sample QC TSV. Returns dict of sample_id -> het_rate."""
    het_rates = {}
    with open(qc_file) as f:
        header = f.readline().rstrip('\n').split('\t')
        col_idx = {name: i for i, name in enumerate(header)}

        sid_idx = col_idx.get('sample_id')
        hr_idx = col_idx.get('het_rate')

        if sid_idx is None or hr_idx is None:
            print("Error: QC file missing required columns "
                  "(sample_id, het_rate)", file=sys.stderr)
            sys.exit(1)

        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= max(sid_idx, hr_idx):
                continue
            sid = fields[sid_idx]
            try:
                het_rates[sid] = float(fields[hr_idx])
            except (ValueError, TypeError):
                pass

    return het_rates


def load_ancestry_assignments(peddy_het_check):
    """Read ancestry-prediction from peddy het_check.csv.

    Returns dict of sample_id -> ancestry_prediction.
    """
    assignments = {}
    with open(peddy_het_check) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sid = row.get('sample_id', '')
            anc = row.get('ancestry-prediction', 'UNKNOWN')
            if sid:
                assignments[sid] = anc
    return assignments


def filter_het_outliers(qc_file, peddy_het_check, output_excluded, output_log,
                        sd_threshold):
    """Flag samples with het rate > N SD from their ancestry-group mean."""
    het_rates = load_het_rates(qc_file)
    ancestry = load_ancestry_assignments(peddy_het_check)

    # Group het rates by ancestry
    ancestry_groups = {}
    for sid, hr in het_rates.items():
        anc = ancestry.get(sid, 'UNKNOWN')
        ancestry_groups.setdefault(anc, []).append((sid, hr))

    excluded = set()
    details = []  # (sample_id, ancestry, het_rate, ancestry_mean, ancestry_sd, direction)

    for anc, samples in sorted(ancestry_groups.items()):
        rates = [hr for _, hr in samples]
        if len(rates) < 2:
            continue

        mean_hr = statistics.mean(rates)
        sd_hr = statistics.stdev(rates)

        if sd_hr == 0:
            continue

        for sid, hr in samples:
            deviation = abs(hr - mean_hr)
            if deviation > sd_threshold * sd_hr:
                excluded.add(sid)
                direction = 'high' if hr > mean_hr else 'low'
                details.append((sid, anc, hr, mean_hr, sd_hr, direction))

    # Write excluded sample IDs
    with open(output_excluded, 'w') as f:
        for sid in sorted(excluded):
            f.write(sid + '\n')

    # Write detailed log
    if output_log:
        with open(output_log, 'w') as f:
            f.write('sample_id\tancestry\thet_rate\tancestry_mean\t'
                    'ancestry_sd\tdirection\n')
            for sid, anc, hr, mean_hr, sd_hr, direction in details:
                f.write(f'{sid}\t{anc}\t{hr:.6f}\t{mean_hr:.6f}\t'
                        f'{sd_hr:.6f}\t{direction}\n')

    # Summary
    print(f"  Heterozygosity outlier filtering (± {sd_threshold} SD):")
    print(f"    Total samples with het_rate:  {len(het_rates)}")
    print(f"    Ancestry groups evaluated:    {len(ancestry_groups)}")
    print(f"    Outlier samples excluded:     {len(excluded)}")
    for anc, samples in sorted(ancestry_groups.items()):
        n_exc = sum(1 for sid, _ in samples if sid in excluded)
        if n_exc > 0:
            print(f"      {anc}: {n_exc} outlier(s) of {len(samples)} samples")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Exclude samples whose heterozygosity rate is > N SD from the "
            "per-ancestry-group mean (GWAS QC best practice)."
        )
    )
    parser.add_argument("--sample-qc", required=True,
                        help="Sample QC TSV file (must have sample_id, het_rate)")
    parser.add_argument("--peddy-het-check", required=True,
                        help="Peddy het_check.csv with ancestry-prediction")
    parser.add_argument("--output-excluded", required=True,
                        help="Output file: one excluded sample ID per line")
    parser.add_argument("--output-log",
                        help="Output TSV log of outlier details")
    parser.add_argument("--sd-threshold", type=float, default=3.0,
                        help="Number of SDs from mean to flag (default: 3)")
    args = parser.parse_args()

    filter_het_outliers(args.sample_qc, args.peddy_het_check,
                        args.output_excluded, args.output_log,
                        args.sd_threshold)


if __name__ == "__main__":
    main()
