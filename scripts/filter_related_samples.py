#!/usr/bin/env python3
"""
filter_related_samples.py

Identify related sample pairs from peddy ped_check.csv output and select
samples to exclude prior to PCA, HWE testing, and ancestry-stratified QC.

For each pair with estimated relatedness >= threshold (default 0.125,
second-degree or closer), the sample with lower call rate (or, when tied,
higher LRR SD) is selected for exclusion.  This implements the GWAS QC
best practice of restricting PCA and HWE analyses to an unrelated sample
set (Anderson et al. 2010; Marees et al. 2018).

Excluded samples are still projected onto PCs downstream (flashpca2
--project) so they retain ancestry estimates.

Usage:
    python3 filter_related_samples.py \
        --ped-check peddy.ped_check.csv \
        --sample-qc stage2_sample_qc.tsv \
        --output-excluded relatedness_excluded.txt \
        --output-log relatedness_exclusions.tsv \
        --min-relatedness 0.125

References:
    - Anderson et al. 2010, Nat Protoc 5:1564-73
    - Marees et al. 2018, Int J Methods Psychiatr Res 27:e1608
"""

import argparse
import csv
import sys


def load_sample_qc(qc_file):
    """Read sample QC TSV and return dict of sample_id -> {call_rate, lrr_sd}."""
    qc = {}
    with open(qc_file) as f:
        header = f.readline().rstrip('\n').split('\t')
        col_idx = {name: i for i, name in enumerate(header)}

        sid_idx = col_idx.get('sample_id')
        cr_idx = col_idx.get('call_rate')
        sd_idx = col_idx.get('lrr_sd')

        if sid_idx is None:
            print("Error: QC file missing 'sample_id' column", file=sys.stderr)
            sys.exit(1)

        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= sid_idx:
                continue
            sid = fields[sid_idx]

            call_rate = None
            if cr_idx is not None and cr_idx < len(fields):
                try:
                    call_rate = float(fields[cr_idx])
                except (ValueError, TypeError):
                    pass

            lrr_sd = None
            if sd_idx is not None and sd_idx < len(fields):
                val = fields[sd_idx]
                if val not in ('NA', '', 'nan', 'NaN'):
                    try:
                        lrr_sd = float(val)
                    except (ValueError, TypeError):
                        pass

            qc[sid] = {'call_rate': call_rate, 'lrr_sd': lrr_sd}

    return qc


def choose_exclusion(sample_a, sample_b, qc_data):
    """Choose which sample to exclude from a related pair.

    Prefer to exclude the sample with:
      1. Lower call rate (worse genotyping quality)
      2. Higher LRR SD (more intensity noise) as tiebreaker
    If QC data is missing for one sample, exclude that one.
    """
    qa = qc_data.get(sample_a, {})
    qb = qc_data.get(sample_b, {})

    cr_a = qa.get('call_rate')
    cr_b = qb.get('call_rate')
    sd_a = qa.get('lrr_sd')
    sd_b = qb.get('lrr_sd')

    # If one sample has no QC data at all, exclude it
    if cr_a is None and cr_b is not None:
        return sample_a
    if cr_b is None and cr_a is not None:
        return sample_b

    # Compare call rates (lower is worse → exclude)
    if cr_a is not None and cr_b is not None:
        if cr_a < cr_b:
            return sample_a
        if cr_b < cr_a:
            return sample_b

    # Tie on call rate → compare LRR SD (higher is worse → exclude)
    if sd_a is not None and sd_b is not None:
        if sd_a > sd_b:
            return sample_a
        if sd_b > sd_a:
            return sample_b

    # Complete tie or both missing → exclude first alphabetically for
    # deterministic behavior
    return min(sample_a, sample_b)


def filter_related(ped_check_file, qc_file, output_excluded, output_log,
                   min_relatedness):
    """Identify related pairs and write exclusion list."""
    qc_data = {}
    if qc_file:
        qc_data = load_sample_qc(qc_file)

    # Read related pairs from peddy ped_check.csv
    related_pairs = []
    with open(ped_check_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                rel = float(row.get('rel', 0))
            except (ValueError, TypeError):
                continue

            if rel >= min_relatedness:
                sample_a = row.get('sample_a', '')
                sample_b = row.get('sample_b', '')
                if sample_a and sample_b:
                    related_pairs.append((sample_a, sample_b, rel))

    # Greedy exclusion: for each pair, exclude the worse-quality sample.
    # Track exclusions to avoid double-excluding.
    excluded = set()
    exclusion_log = []

    for sample_a, sample_b, rel in related_pairs:
        # If one is already excluded, the pair is already broken
        if sample_a in excluded or sample_b in excluded:
            continue

        to_exclude = choose_exclusion(sample_a, sample_b, qc_data)
        excluded.add(to_exclude)
        kept = sample_b if to_exclude == sample_a else sample_a
        exclusion_log.append((to_exclude, kept, rel))

    # Write excluded sample IDs
    with open(output_excluded, 'w') as f:
        for sid in sorted(excluded):
            f.write(sid + '\n')

    # Write detailed log
    if output_log:
        with open(output_log, 'w') as f:
            f.write('excluded_sample\tkept_sample\trelatedness\n')
            for exc, kept, rel in exclusion_log:
                f.write(f'{exc}\t{kept}\t{rel:.4f}\n')

    # Summary
    print(f"  Relatedness filtering (threshold: {min_relatedness}):")
    print(f"    Related pairs found:    {len(related_pairs)}")
    print(f"    Samples excluded:       {len(excluded)}")
    if exclusion_log:
        print(f"    First excluded samples:")
        for exc, kept, rel in exclusion_log[:5]:
            print(f"      {exc} (rel={rel:.4f} with {kept})")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Exclude one sample per related pair (relatedness >= threshold) "
            "prior to PCA and HWE testing. Selects the sample with lower "
            "call rate / higher LRR SD for exclusion."
        )
    )
    parser.add_argument("--ped-check", required=True,
                        help="Peddy ped_check.csv file")
    parser.add_argument("--sample-qc",
                        help="Sample QC TSV with call_rate and lrr_sd columns")
    parser.add_argument("--output-excluded", required=True,
                        help="Output file: one excluded sample ID per line")
    parser.add_argument("--output-log",
                        help="Output TSV log of exclusion decisions")
    parser.add_argument("--min-relatedness", type=float, default=0.125,
                        help="Minimum relatedness to flag a pair (default: 0.125)")
    args = parser.parse_args()

    filter_related(args.ped_check, args.sample_qc, args.output_excluded,
                   args.output_log, args.min_relatedness)


if __name__ == "__main__":
    main()
