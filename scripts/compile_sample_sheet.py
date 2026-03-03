#!/usr/bin/env python3
"""
compile_sample_sheet.py

Compile a unified sample sheet combining all QC metrics and ancestry PCs
from the Illumina IDAT processing pipeline.

Merges:
  - Sample QC metrics (call rate, LRR SD, LRR mean, LRR median, predicted sex)
  - Ancestry PCA projections (default 20 PCs)

Usage:
    python3 compile_sample_sheet.py \\
        --sample-qc stage2/qc/stage2_sample_qc.tsv \\
        --pca-projections pca/pca_projections.tsv \\
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


def main():
    parser = argparse.ArgumentParser(
        description="Compile unified sample sheet with QC metrics and PCs"
    )
    parser.add_argument("--sample-qc", required=True,
                        help="Sample QC TSV file (from collect_qc_metrics.sh)")
    parser.add_argument("--pca-projections", default=None,
                        help="PCA projections TSV (from ancestry_pca.sh)")
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

    # Build unified header
    output_cols = list(qc_header)
    output_cols.extend(pc_cols)

    # Collect all sample IDs
    all_samples = sorted(set(list(qc_data.keys()) + list(pc_data.keys())))

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
            # PC columns
            pc_row = pc_data.get(sample_id, {})
            for col in pc_cols:
                row.append(pc_row.get(col, 'NA'))
            f.write('\t'.join(row) + '\n')

    n_samples = len(all_samples)
    n_with_pcs = sum(1 for s in all_samples if s in pc_data)
    print(f"Compiled sample sheet: {args.output}")
    print(f"  Total samples:     {n_samples}")
    print(f"  QC columns:        {len(qc_header)}")
    print(f"  PC columns:        {len(pc_cols)}")
    print(f"  Samples with PCs:  {n_with_pcs}")


if __name__ == "__main__":
    main()
