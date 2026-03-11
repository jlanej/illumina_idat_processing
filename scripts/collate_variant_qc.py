#!/usr/bin/env python3
"""
collate_variant_qc.py

Collate per-variant QC metrics from the full cohort and ancestry-specific
subsets into a single unified TSV file. Each variant has columns for
call rate, HWE p-value, and MAF from each ancestry group alongside the
full-cohort values. The output also includes cross-ancestry pass flags
indicating whether each variant passes call-rate/HWE/MAF thresholds in
all ancestry subsets that were run, plus an overall pass flag.

This enables analysts to compare variant quality across ancestry strata,
which is essential because HWE tests assume panmixia and can produce
misleading results in mixed-ancestry samples (Wahlund effect).

Usage:
    python3 collate_variant_qc.py \\
        --all-variant-qc-dir stage2/qc/variant_qc \\
        --ancestry-variant-qc EUR:ancestry_qc/EUR/variant_qc \\
        --ancestry-variant-qc AFR:ancestry_qc/AFR/variant_qc \\
        --output collated_variant_qc.tsv
"""

import argparse
import os
import sys

CALL_RATE_THRESHOLD = 0.98
HWE_P_THRESHOLD = 1e-6
MAF_THRESHOLD = 0.01


def _safe_float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def read_plink2_vmiss(filepath):
    """Read plink2 .vmiss file. Returns dict of variant_id -> missingness."""
    data = {}
    if not os.path.exists(filepath):
        return data
    with open(filepath) as f:
        header = f.readline().strip().split()
        # plink2 .vmiss columns: #CHROM  ID  MISSING_CT  OBS_CT  F_MISS
        id_idx = header.index('ID') if 'ID' in header else 1
        fmiss_idx = header.index('F_MISS') if 'F_MISS' in header else len(header) - 1
        for line in f:
            fields = line.strip().split()
            if len(fields) > max(id_idx, fmiss_idx):
                data[fields[id_idx]] = fields[fmiss_idx]
    return data


def read_plink2_hardy(filepath):
    """Read plink2 .hardy file. Returns dict of variant_id -> HWE p-value."""
    data = {}
    if not os.path.exists(filepath):
        return data
    with open(filepath) as f:
        header = f.readline().strip().split()
        id_idx = header.index('ID') if 'ID' in header else 1
        # P column is typically the last one
        p_idx = header.index('P') if 'P' in header else len(header) - 1
        for line in f:
            fields = line.strip().split()
            if len(fields) > max(id_idx, p_idx):
                data[fields[id_idx]] = fields[p_idx]
    return data


def read_plink2_afreq(filepath):
    """Read plink2 .afreq file. Returns dict of variant_id -> MAF."""
    data = {}
    if not os.path.exists(filepath):
        return data
    with open(filepath) as f:
        header = f.readline().strip().split()
        id_idx = header.index('ID') if 'ID' in header else 1
        af_idx = header.index('ALT_FREQS') if 'ALT_FREQS' in header else 4
        for line in f:
            fields = line.strip().split()
            if len(fields) > max(id_idx, af_idx):
                try:
                    af = float(fields[af_idx])
                    maf = min(af, 1 - af)
                    data[fields[id_idx]] = f'{maf:.6f}'
                except (ValueError, IndexError):
                    data[fields[id_idx]] = 'NA'
    return data


def load_variant_qc(qc_dir, prefix='variant_qc'):
    """Load all variant QC metrics from a plink2 output directory.

    Returns dict of variant_id -> {call_rate, hwe_p, maf}.
    """
    vmiss = read_plink2_vmiss(os.path.join(qc_dir, f'{prefix}.vmiss'))
    hardy = read_plink2_hardy(os.path.join(qc_dir, f'{prefix}.hardy'))
    afreq = read_plink2_afreq(os.path.join(qc_dir, f'{prefix}.afreq'))

    all_variants = sorted(set(list(vmiss.keys()) + list(hardy.keys()) + list(afreq.keys())))

    data = {}
    for vid in all_variants:
        fmiss = vmiss.get(vid, 'NA')
        try:
            cr = f'{1 - float(fmiss):.6f}' if fmiss != 'NA' else 'NA'
        except ValueError:
            cr = 'NA'
        data[vid] = {
            'call_rate': cr,
            'hwe_p': hardy.get(vid, 'NA'),
            'maf': afreq.get(vid, 'NA'),
        }
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Collate variant QC metrics across ancestry groups"
    )
    parser.add_argument("--all-variant-qc-dir", default=None,
                        help="Directory with full-cohort variant QC files")
    parser.add_argument("--ancestry-variant-qc", action='append', default=[],
                        help="Ancestry-specific QC in format ANCESTRY:DIR (repeatable)")
    parser.add_argument("--output", required=True,
                        help="Output collated TSV file")
    args = parser.parse_args()

    # Parse ancestry:dir pairs
    ancestry_dirs = []
    for spec in args.ancestry_variant_qc:
        if ':' not in spec:
            print(f"Warning: Invalid ancestry spec (expected ANCESTRY:DIR): {spec}",
                  file=sys.stderr)
            continue
        anc, qc_dir = spec.split(':', 1)
        ancestry_dirs.append((anc, qc_dir))

    # Load full-cohort QC
    all_data = {}
    if args.all_variant_qc_dir and os.path.isdir(args.all_variant_qc_dir):
        all_data = load_variant_qc(args.all_variant_qc_dir)
        print(f"Full-cohort variant QC: {len(all_data)} variants")

    # Load per-ancestry QC
    ancestry_data = {}
    for anc, qc_dir in ancestry_dirs:
        if os.path.isdir(qc_dir):
            ancestry_data[anc] = load_variant_qc(qc_dir)
            print(f"  {anc} variant QC: {len(ancestry_data[anc])} variants")

    # Collect all variant IDs
    all_variants = set(all_data.keys())
    for anc_data in ancestry_data.values():
        all_variants.update(anc_data.keys())
    all_variants = sorted(all_variants)

    if not all_variants:
        print("Warning: No variant QC data found.", file=sys.stderr)
        # Write empty header
        with open(args.output, 'w') as f:
            f.write("variant_id\n")
        return

    # Build output columns
    columns = ['variant_id']
    if all_data:
        columns.extend(['all_call_rate', 'all_hwe_p', 'all_maf'])
    ancestries_sorted = sorted(ancestry_data.keys())
    for anc in ancestries_sorted:
        columns.extend([
            f'{anc}_call_rate',
            f'{anc}_hwe_p',
            f'{anc}_maf',
        ])
    columns.extend([
        'all_ancestries_call_rate_pass',
        'all_ancestries_hwe_pass',
        'all_ancestries_maf_pass',
        'all_ancestries_qc_pass',
    ])

    # Write collated output
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, 'w') as f:
        f.write('\t'.join(columns) + '\n')
        for vid in all_variants:
            row = [vid]
            if all_data:
                vd = all_data.get(vid, {})
                row.extend([
                    vd.get('call_rate', 'NA'),
                    vd.get('hwe_p', 'NA'),
                    vd.get('maf', 'NA'),
                ])
            for anc in ancestries_sorted:
                vd = ancestry_data.get(anc, {}).get(vid, {})
                row.extend([
                    vd.get('call_rate', 'NA'),
                    vd.get('hwe_p', 'NA'),
                    vd.get('maf', 'NA'),
                ])

            if ancestries_sorted:
                cr_pass = True
                hwe_pass = True
                maf_pass = True
                for anc in ancestries_sorted:
                    anc_variant = ancestry_data.get(anc, {}).get(vid, {})
                    cr_val = _safe_float(anc_variant.get('call_rate'))
                    hwe_val = _safe_float(anc_variant.get('hwe_p'))
                    maf_val = _safe_float(anc_variant.get('maf'))
                    cr_pass = cr_pass and cr_val is not None and cr_val >= CALL_RATE_THRESHOLD
                    hwe_pass = hwe_pass and hwe_val is not None and hwe_val >= HWE_P_THRESHOLD
                    maf_pass = maf_pass and maf_val is not None and maf_val >= MAF_THRESHOLD
                all_qc_pass = cr_pass and hwe_pass and maf_pass
                row.extend([
                    '1' if cr_pass else '0',
                    '1' if hwe_pass else '0',
                    '1' if maf_pass else '0',
                    '1' if all_qc_pass else '0',
                ])
            else:
                row.extend(['NA', 'NA', 'NA', 'NA'])
            f.write('\t'.join(row) + '\n')

    print(f"Collated variant QC: {args.output}")
    print(f"  Total variants: {len(all_variants)}")
    print(f"  Columns: {', '.join(columns)}")


if __name__ == "__main__":
    main()
