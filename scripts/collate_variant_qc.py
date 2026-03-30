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
    """Read plink2 .vmiss file. Returns dict of variant_id -> row dict.

    Returned fields per variant: F_MISS, OBS_CT, MISSING_CT.
    """
    data = {}
    if not os.path.exists(filepath):
        return data
    with open(filepath) as f:
        header = f.readline().strip().split()
        id_idx = header.index('ID') if 'ID' in header else 1
        fmiss_idx = header.index('F_MISS') if 'F_MISS' in header else len(header) - 1
        obs_idx = header.index('OBS_CT') if 'OBS_CT' in header else None
        miss_idx = header.index('MISSING_CT') if 'MISSING_CT' in header else None
        for line in f:
            fields = line.strip().split()
            if len(fields) > max(id_idx, fmiss_idx):
                row = {'F_MISS': fields[fmiss_idx]}
                if obs_idx is not None and obs_idx < len(fields):
                    row['OBS_CT'] = fields[obs_idx]
                if miss_idx is not None and miss_idx < len(fields):
                    row['MISSING_CT'] = fields[miss_idx]
                data[fields[id_idx]] = row
    return data


def read_plink2_hardy(filepath):
    """Read plink2 .hardy file. Returns dict of variant_id -> row dict.

    Returned fields per variant: P, HOM_A1_CT, HET_A1_AX_CT, TWO_AX_CT.
    """
    data = {}
    if not os.path.exists(filepath):
        return data
    with open(filepath) as f:
        header = f.readline().strip().split()
        id_idx = header.index('ID') if 'ID' in header else 1
        p_idx = header.index('P') if 'P' in header else len(header) - 1
        hom_a1_idx = header.index('HOM_A1_CT') if 'HOM_A1_CT' in header else None
        het_idx = header.index('HET_A1_AX_CT') if 'HET_A1_AX_CT' in header else None
        two_ax_idx = header.index('TWO_AX_CT') if 'TWO_AX_CT' in header else None
        for line in f:
            fields = line.strip().split()
            if len(fields) > max(id_idx, p_idx):
                row = {'P': fields[p_idx]}
                if hom_a1_idx is not None and hom_a1_idx < len(fields):
                    row['HOM_A1_CT'] = fields[hom_a1_idx]
                if het_idx is not None and het_idx < len(fields):
                    row['HET_A1_AX_CT'] = fields[het_idx]
                if two_ax_idx is not None and two_ax_idx < len(fields):
                    row['TWO_AX_CT'] = fields[two_ax_idx]
                data[fields[id_idx]] = row
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

    Returns dict of variant_id -> {call_rate, hwe_p, maf, obs_ct, missing_ct,
                                    hom_a1_ct, het_ct, hom_a2_ct}.
    """
    vmiss = read_plink2_vmiss(os.path.join(qc_dir, f'{prefix}.vmiss'))
    hardy = read_plink2_hardy(os.path.join(qc_dir, f'{prefix}.hardy'))
    afreq = read_plink2_afreq(os.path.join(qc_dir, f'{prefix}.afreq'))

    all_variants = sorted(set(list(vmiss.keys()) + list(hardy.keys()) + list(afreq.keys())))

    data = {}
    for vid in all_variants:
        vmiss_row = vmiss.get(vid, {})
        fmiss = vmiss_row.get('F_MISS', 'NA') if isinstance(vmiss_row, dict) else 'NA'
        try:
            cr = f'{1 - float(fmiss):.6f}' if fmiss != 'NA' else 'NA'
        except ValueError:
            cr = 'NA'

        hardy_row = hardy.get(vid, {})
        hwe_p = hardy_row.get('P', 'NA') if isinstance(hardy_row, dict) else 'NA'

        data[vid] = {
            'call_rate': cr,
            'hwe_p': hwe_p,
            'maf': afreq.get(vid, 'NA'),
            'obs_ct': vmiss_row.get('OBS_CT', 'NA') if isinstance(vmiss_row, dict) else 'NA',
            'missing_ct': vmiss_row.get('MISSING_CT', 'NA') if isinstance(vmiss_row, dict) else 'NA',
            'hom_a1_ct': hardy_row.get('HOM_A1_CT', 'NA') if isinstance(hardy_row, dict) else 'NA',
            'het_ct': hardy_row.get('HET_A1_AX_CT', 'NA') if isinstance(hardy_row, dict) else 'NA',
            'hom_a2_ct': hardy_row.get('TWO_AX_CT', 'NA') if isinstance(hardy_row, dict) else 'NA',
        }
    return data


def parse_tstv_file(filepath):
    """Parse bcftools stats tstv output for project-level Ti/Tv ratio.

    Returns the Ti/Tv ratio as a string, or 'NA' if not available.
    """
    if not filepath or not os.path.exists(filepath):
        return 'NA'
    with open(filepath) as f:
        for line in f:
            # TSTV line format: TSTV\t0\tts\ttv\tts/tv\t...
            if line.startswith('TSTV\t'):
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    return fields[4]
    return 'NA'


def main():
    parser = argparse.ArgumentParser(
        description="Collate variant QC metrics across ancestry groups"
    )
    parser.add_argument("--all-variant-qc-dir", default=None,
                        help="Directory with full-cohort variant QC files")
    parser.add_argument("--ancestry-variant-qc", action='append', default=[],
                        help="Ancestry-specific QC in format ANCESTRY:DIR (repeatable)")
    parser.add_argument("--tstv-file", default=None,
                        help="bcftools stats tstv output file for project-level "
                             "Ti/Tv ratio (tstv_stats.txt)")
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

    # Parse project-level Ti/Tv ratio
    tstv_ratio = parse_tstv_file(args.tstv_file)

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
        columns.extend([
            'all_call_rate', 'all_hwe_p', 'all_maf',
            'all_obs_ct', 'all_missing_ct',
            'all_hom_a1_ct', 'all_het_ct', 'all_hom_a2_ct',
        ])
    ancestries_sorted = sorted(ancestry_data.keys())
    for anc in ancestries_sorted:
        columns.extend([
            f'{anc}_call_rate',
            f'{anc}_hwe_p',
            f'{anc}_maf',
            f'{anc}_obs_ct',
            f'{anc}_missing_ct',
            f'{anc}_hom_a1_ct',
            f'{anc}_het_ct',
            f'{anc}_hom_a2_ct',
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
        # Write project-level Ti/Tv as a header comment if available
        if tstv_ratio != 'NA':
            f.write(f'#tstv_ratio={tstv_ratio}\n')
        f.write('\t'.join(columns) + '\n')
        for vid in all_variants:
            row = [vid]
            if all_data:
                vd = all_data.get(vid, {})
                row.extend([
                    vd.get('call_rate', 'NA'),
                    vd.get('hwe_p', 'NA'),
                    vd.get('maf', 'NA'),
                    vd.get('obs_ct', 'NA'),
                    vd.get('missing_ct', 'NA'),
                    vd.get('hom_a1_ct', 'NA'),
                    vd.get('het_ct', 'NA'),
                    vd.get('hom_a2_ct', 'NA'),
                ])
            for anc in ancestries_sorted:
                vd = ancestry_data.get(anc, {}).get(vid, {})
                row.extend([
                    vd.get('call_rate', 'NA'),
                    vd.get('hwe_p', 'NA'),
                    vd.get('maf', 'NA'),
                    vd.get('obs_ct', 'NA'),
                    vd.get('missing_ct', 'NA'),
                    vd.get('hom_a1_ct', 'NA'),
                    vd.get('het_ct', 'NA'),
                    vd.get('hom_a2_ct', 'NA'),
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
    if tstv_ratio != 'NA':
        print(f"  Project Ti/Tv ratio: {tstv_ratio}")


if __name__ == "__main__":
    main()
