#!/usr/bin/env python3
"""
enrich_qc_with_names.py

Add an original_name column to a QC TSV file using a sample name mapping.

The name mapping file has the format (with header):
    sample_id\toriginal_name
    NA19152\t5426529111_R04C01

The QC file has the format (with header):
    sample_id\tcall_rate\tlrr_sd\t...

This script inserts the original_name column after sample_id, producing:
    sample_id\toriginal_name\tcall_rate\tlrr_sd\t...

If a sample_id is not found in the mapping, original_name defaults to sample_id.

Usage:
    python3 enrich_qc_with_names.py \
        --qc-file stage1_sample_qc.tsv \
        --name-map sample_name_mapping.tsv \
        --output stage1_sample_qc_enriched.tsv
"""

import argparse
import sys


def enrich_qc(qc_file, name_map_file, output_file):
    """Add original_name column to QC file using the name mapping."""

    # Read name mapping (sample_id -> original_name)
    name_map = {}
    with open(name_map_file) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 2:
                # In the mapping file: col1 = sample_id (new name),
                # col2 = original_name (Sentrix ID)
                name_map[fields[0]] = fields[1]

    # Read QC file and insert original_name column
    with open(qc_file) as fin, open(output_file, 'w') as fout:
        qc_header = fin.readline().rstrip('\n').split('\t')

        # Insert 'original_name' after 'sample_id'
        new_header = [qc_header[0], 'original_name'] + qc_header[1:]
        fout.write('\t'.join(new_header) + '\n')

        for line in fin:
            fields = line.rstrip('\n').split('\t')
            if not fields:
                continue
            sample_id = fields[0]
            original_name = name_map.get(sample_id, sample_id)
            new_fields = [sample_id, original_name] + fields[1:]
            fout.write('\t'.join(new_fields) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Add original_name column to QC TSV from name mapping"
    )
    parser.add_argument("--qc-file", required=True,
                        help="Input QC TSV file")
    parser.add_argument("--name-map", required=True,
                        help="Name mapping TSV file (sample_id → original_name)")
    parser.add_argument("--output", required=True,
                        help="Output enriched QC TSV file")
    args = parser.parse_args()

    enrich_qc(args.qc_file, args.name_map, args.output)


if __name__ == "__main__":
    main()

