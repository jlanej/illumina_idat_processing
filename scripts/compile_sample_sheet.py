#!/usr/bin/env python3
"""
compile_sample_sheet.py

Compile a unified sample sheet combining all QC metrics and ancestry PCs
from the Illumina IDAT processing pipeline.

Merges:
  - Sample QC metrics (call rate, LRR SD/mean/median, BAF mean/SD, het rate)
  - Sex check cross-tabulation (chrX/Y LRR medians, F-stat, peddy sex,
    concordance status)
  - Pre-PCA exclusion flags (relatedness, het outlier, combined)
  - Inbreeding coefficient (F) from plink2 .het
  - Ancestry PCA projections (default 20 PCs)
  - Peddy het_check and sex_check metrics

Usage:
    python3 compile_sample_sheet.py \\
        --sample-qc stage2/qc/stage2_sample_qc.tsv \\
        --pca-projections pca/pca_projections.tsv \\
        --sex-check sex_check/sex_check_chrXY_lrr.tsv \\
        --relatedness-excluded relatedness_excluded_samples.txt \\
        --het-outlier-excluded het_outlier_excluded_samples.txt \\
        --pre-pca-excluded pre_pca_excluded_samples.txt \\
        --output compiled_sample_sheet.tsv
"""

import argparse
import os
import sys


def read_tsv(filepath):
    """Read a TSV file. Returns (header_list, dict of id -> row_dict)."""
    rows = {}
    header = []
    if not os.path.exists(filepath):
        return header, rows
    with open(filepath) as f:
        header_line = f.readline().strip()
        if not header_line:
            return header, rows
        header = header_line.split('\t')
        for line in f:
            fields = line.strip().split('\t')
            if not fields:
                continue
            row = {}
            for i, col in enumerate(header):
                row[col] = fields[i] if i < len(fields) else 'NA'
            # Use first column as ID, or IID if present
            row_id = fields[0]
            rows[row_id] = row
    return header, rows


def read_pca_projections(filepath):
    """Read PCA projections TSV (from flashpca2 or plink2).

    Handles both tab-separated and space-separated formats.

    Returns (pc_columns, dict of sample_id -> {PC1: val, PC2: val, ...}).
    """
    pc_data = {}
    pc_cols = []
    if not os.path.exists(filepath):
        return pc_cols, pc_data
    with open(filepath) as f:
        header_line = f.readline().strip()
        if not header_line:
            return pc_cols, pc_data
        # Detect delimiter: use tab if present, otherwise whitespace
        if '\t' in header_line:
            header = header_line.split('\t')
            delim = '\t'
        else:
            header = header_line.split()
            delim = None  # split on any whitespace
        # Find ID column: prefer IID, then second column, then first
        id_idx = 1 if len(header) > 1 else 0
        for i, col in enumerate(header):
            if col.upper() in ('IID', '#IID'):
                id_idx = i
                break
        # PC columns: everything after FID/IID
        fid_iid_indices = [i for i, h in enumerate(header)
                           if h.upper() in ('FID', 'IID', '#FID', '#IID')]
        pc_start = max(fid_iid_indices) + 1 if fid_iid_indices else 2
        pc_cols = header[pc_start:]

        for line in f:
            fields = line.strip().split(delim) if delim else line.strip().split()
            if len(fields) <= id_idx:
                continue
            sample_id = fields[id_idx]
            row = {}
            for j, col in enumerate(pc_cols):
                idx = pc_start + j
                row[col] = fields[idx] if idx < len(fields) else 'NA'
            pc_data[sample_id] = row
    return pc_cols, pc_data


def read_peddy_csv(filepath):
    """Read a peddy output CSV file.

    Returns (columns, dict of sample_id -> row_dict).
    """
    import csv
    data = {}
    columns = []
    if not filepath or not os.path.exists(filepath):
        return columns, data
    with open(filepath) as f:
        reader = csv.DictReader(f)
        columns = reader.fieldnames or []
        for row in reader:
            sid = row.get('sample_id', '')
            if sid:
                data[sid] = row
    return columns, data


def read_exclusion_list(filepath):
    """Read a sample exclusion list (one sample ID per line).

    Returns a set of excluded sample IDs.
    """
    excluded = set()
    if not filepath or not os.path.exists(filepath):
        return excluded
    with open(filepath) as f:
        for line in f:
            sid = line.strip()
            if sid:
                excluded.add(sid)
    return excluded


def main():
    parser = argparse.ArgumentParser(
        description="Compile unified sample sheet with QC metrics and PCs"
    )
    parser.add_argument("--sample-qc", required=True,
                        help="Sample QC TSV file (from collect_qc_metrics.sh)")
    parser.add_argument("--pca-projections", default=None,
                        help="PCA projections TSV (from ancestry_pca.sh)")
    parser.add_argument("--het-file", default=None,
                        help="Plink2 .het file with per-sample inbreeding coefficient")
    parser.add_argument("--peddy-het-check", default=None,
                        help="Peddy het_check.csv with ancestry/het QC metrics")
    parser.add_argument("--peddy-sex-check", default=None,
                        help="Peddy sex_check.csv with sex prediction metrics")
    parser.add_argument("--sex-check", default=None,
                        help="Sex check cross-tabulation TSV "
                             "(sex_check_chrXY_lrr.tsv from plot_sex_check.py)")
    parser.add_argument("--relatedness-excluded", default=None,
                        help="Relatedness-excluded samples list "
                             "(relatedness_excluded_samples.txt)")
    parser.add_argument("--het-outlier-excluded", default=None,
                        help="Het-outlier-excluded samples list "
                             "(het_outlier_excluded_samples.txt)")
    parser.add_argument("--pre-pca-excluded", default=None,
                        help="Combined pre-PCA excluded samples list "
                             "(pre_pca_excluded_samples.txt)")
    parser.add_argument("--ancestry-pca", action='append', default=[],
                        help="Ancestry-specific PCA projections in format "
                             "ANCESTRY:FILE (repeatable). Samples not of a "
                             "given ancestry will have NaN for those PCs.")
    parser.add_argument("--output", required=True,
                        help="Output compiled sample sheet TSV")
    parser.add_argument("--n-pcs", type=int, default=20,
                        help="Number of PCs to include (default: 20)")
    args = parser.parse_args()

    if not os.path.exists(args.sample_qc):
        print(f"Error: Sample QC file not found: {args.sample_qc}",
              file=sys.stderr)
        sys.exit(1)

    # Read QC metrics
    qc_header, qc_data = read_tsv(args.sample_qc)
    if not qc_data:
        print(f"Error: No data in QC file: {args.sample_qc}", file=sys.stderr)
        sys.exit(1)

    # Read PCA projections if available
    pc_cols, pc_data = [], {}
    if args.pca_projections and os.path.exists(args.pca_projections):
        pc_cols, pc_data = read_pca_projections(args.pca_projections)
        # Limit to requested number of PCs
        pc_cols = pc_cols[:args.n_pcs]

    # Read inbreeding coefficient (F) from plink2 .het file if available
    het_data = {}
    if args.het_file and os.path.exists(args.het_file):
        with open(args.het_file) as f:
            het_header = f.readline().strip().split()
            # plink2 .het columns: #FID  IID  OBS_CT  E(HOM_CT)  HOM_CT  F
            iid_idx = 1 if len(het_header) > 1 else 0
            f_idx = len(het_header) - 1  # F is the last column
            for line in f:
                fields = line.strip().split()
                if len(fields) > f_idx:
                    het_data[fields[iid_idx]] = fields[f_idx]

    # Read sex check cross-tabulation if available
    SEX_CHECK_FIELDS = [
        'chrx_lrr_median', 'chry_lrr_median',
        'peddy_sex', 'sex_status',
    ]
    sex_check_cols = []
    sex_check_data = {}
    if args.sex_check and os.path.exists(args.sex_check):
        _sc_header, _sc_data = read_tsv(args.sex_check)
        for field in SEX_CHECK_FIELDS:
            if field in _sc_header:
                sex_check_cols.append(field)
        for sid, row in _sc_data.items():
            sex_check_data[sid] = {col: row.get(col, 'NA')
                                   for col in sex_check_cols}

    # Read pre-PCA exclusion lists if available
    rel_excluded = read_exclusion_list(args.relatedness_excluded)
    het_excluded = read_exclusion_list(args.het_outlier_excluded)
    pca_excluded = read_exclusion_list(args.pre_pca_excluded)
    has_exclusion_flags = bool(
        args.relatedness_excluded or args.het_outlier_excluded
        or args.pre_pca_excluded
    )

    # Read peddy sample-level metrics if available
    # Columns are prefixed with 'peddy_' to denote their source.
    peddy_het_cols = []
    peddy_het_data = {}
    PEDDY_HET_FIELDS = [
        'het_ratio', 'mean_depth', 'ancestry-prediction',
        'ancestry-prob', 'PC1', 'PC2', 'PC3', 'PC4',
    ]
    if args.peddy_het_check and os.path.exists(args.peddy_het_check):
        _cols, peddy_het_raw = read_peddy_csv(args.peddy_het_check)
        for field in PEDDY_HET_FIELDS:
            if field in _cols:
                peddy_het_cols.append(f'peddy_{field}')
        for sid, row in peddy_het_raw.items():
            peddy_het_data[sid] = {}
            for field in PEDDY_HET_FIELDS:
                if field in row:
                    peddy_het_data[sid][f'peddy_{field}'] = row[field]

    peddy_sex_cols = []
    peddy_sex_data = {}
    PEDDY_SEX_FIELDS = [
        'predicted_sex', 'het_ratio', 'error',
    ]
    if args.peddy_sex_check and os.path.exists(args.peddy_sex_check):
        _cols, peddy_sex_raw = read_peddy_csv(args.peddy_sex_check)
        for field in PEDDY_SEX_FIELDS:
            if field in _cols:
                peddy_sex_cols.append(f'peddy_sex_{field}')
        for sid, row in peddy_sex_raw.items():
            peddy_sex_data[sid] = {}
            for field in PEDDY_SEX_FIELDS:
                if field in row:
                    peddy_sex_data[sid][f'peddy_sex_{field}'] = row[field]

    # Read ancestry-specific PCA projections if provided
    # Each entry is ANCESTRY:FILE. Samples not belonging to that ancestry
    # will have NaN for those PCs, following best practice for ancestry-
    # stratified analyses (Peterson et al. 2019, AJHG 105:921-935).
    ancestry_pca_data = {}  # ancestry -> (pc_cols, {sample_id -> {col: val}})
    ancestry_pca_order = []
    for spec in args.ancestry_pca:
        if ':' not in spec:
            print(f"Warning: Invalid ancestry-pca spec: {spec}",
                  file=sys.stderr)
            continue
        anc, pca_file = spec.split(':', 1)
        if os.path.exists(pca_file):
            anc_pc_cols, anc_pc_data = read_pca_projections(pca_file)
            anc_pc_cols = anc_pc_cols[:args.n_pcs]
            # Prefix PC columns with ancestry name
            prefixed_cols = [f'{anc}_{col}' for col in anc_pc_cols]
            prefixed_data = {}
            for sid, row in anc_pc_data.items():
                prefixed_data[sid] = {}
                for orig, pref in zip(anc_pc_cols, prefixed_cols):
                    prefixed_data[sid][pref] = row.get(orig, 'NaN')
            ancestry_pca_data[anc] = (prefixed_cols, prefixed_data)
            ancestry_pca_order.append(anc)

    # Build unified header
    output_cols = list(qc_header)
    if het_data:
        output_cols.append('inbreeding_F')
    output_cols.extend(sex_check_cols)
    if has_exclusion_flags:
        output_cols.extend([
            'excluded_relatedness', 'excluded_het_outlier', 'pre_pca_excluded',
        ])
    output_cols.extend(pc_cols)
    # Add ancestry-specific PC columns
    for anc in ancestry_pca_order:
        anc_cols, _ = ancestry_pca_data[anc]
        output_cols.extend(anc_cols)
    output_cols.extend(peddy_het_cols)
    output_cols.extend(peddy_sex_cols)

    # Collect all sample IDs
    all_samples = sorted(set(
        list(qc_data.keys()) + list(pc_data.keys())
        + list(peddy_het_data.keys()) + list(peddy_sex_data.keys())
    ))

    # Write output
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, 'w') as f:
        f.write('\t'.join(output_cols) + '\n')
        for sample_id in all_samples:
            row = []
            # QC columns
            qc_row = qc_data.get(sample_id, {})
            for col in qc_header:
                row.append(qc_row.get(col, 'NA'))
            # Override sample_id if it's the first column
            if qc_header and qc_header[0] == 'sample_id':
                row[0] = sample_id
            # Inbreeding coefficient
            if het_data:
                row.append(het_data.get(sample_id, 'NA'))
            # Sex check cross-tabulation columns
            sc_row = sex_check_data.get(sample_id, {})
            for col in sex_check_cols:
                row.append(sc_row.get(col, 'NA'))
            # Pre-PCA exclusion flags
            if has_exclusion_flags:
                row.append('1' if sample_id in rel_excluded else '0')
                row.append('1' if sample_id in het_excluded else '0')
                row.append('1' if sample_id in pca_excluded else '0')
            # PC columns (full-cohort)
            pc_row = pc_data.get(sample_id, {})
            for col in pc_cols:
                row.append(pc_row.get(col, 'NA'))
            # Ancestry-specific PC columns (NaN for non-members)
            for anc in ancestry_pca_order:
                anc_cols, anc_data_dict = ancestry_pca_data[anc]
                anc_row = anc_data_dict.get(sample_id, {})
                for col in anc_cols:
                    row.append(anc_row.get(col, 'NaN'))
            # Peddy het_check columns
            peddy_het_row = peddy_het_data.get(sample_id, {})
            for col in peddy_het_cols:
                row.append(peddy_het_row.get(col, 'NA'))
            # Peddy sex_check columns
            peddy_sex_row = peddy_sex_data.get(sample_id, {})
            for col in peddy_sex_cols:
                row.append(peddy_sex_row.get(col, 'NA'))
            f.write('\t'.join(row) + '\n')

    n_samples = len(all_samples)
    n_with_pcs = sum(1 for s in all_samples if s in pc_data)
    n_with_peddy = sum(1 for s in all_samples if s in peddy_het_data)
    n_with_sex_check = sum(1 for s in all_samples if s in sex_check_data)
    print(f"Compiled sample sheet: {args.output}")
    print(f"  Total samples:     {n_samples}")
    print(f"  QC columns:        {len(qc_header)}")
    print(f"  PC columns:        {len(pc_cols)}")
    print(f"  Samples with PCs:  {n_with_pcs}")
    if sex_check_cols:
        print(f"  Sex check columns: {len(sex_check_cols)}"
              f" ({n_with_sex_check} samples)")
    if has_exclusion_flags:
        print(f"  Exclusion flags:   excluded_relatedness,"
              f" excluded_het_outlier, pre_pca_excluded")
        print(f"    Relatedness excluded: {len(rel_excluded)}")
        print(f"    Het outlier excluded: {len(het_excluded)}")
        print(f"    Pre-PCA excluded:     {len(pca_excluded)}")
    if ancestry_pca_order:
        for anc in ancestry_pca_order:
            anc_cols, anc_data_dict = ancestry_pca_data[anc]
            n_with_anc_pcs = sum(1 for s in all_samples if s in anc_data_dict)
            print(f"  {anc} PC columns:  {len(anc_cols)} ({n_with_anc_pcs} samples)")
    if peddy_het_cols or peddy_sex_cols:
        print(f"  Peddy columns:     {len(peddy_het_cols) + len(peddy_sex_cols)}")
        print(f"  Samples w/ peddy:  {n_with_peddy}")


if __name__ == "__main__":
    main()
