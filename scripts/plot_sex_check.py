#!/usr/bin/env python3
"""
plot_sex_check.py

Generate a scatter plot of median per-sample LRR on chrX vs median per-sample
LRR on chrY, colored by predicted sex.  This is a standard sex-check
visualization for genotyping array QC.

Males are expected to have lower chrX LRR (hemizygous) and higher chrY LRR,
while females should have higher chrX LRR and lower chrY LRR.

Usage:
    python3 plot_sex_check.py \\
        --vcf stage2/vcf/stage2_reclustered.bcf \\
        --sample-qc stage2/qc/stage2_sample_qc.tsv \\
        --output-dir stage2/qc
"""

import argparse
import os
import subprocess
import sys


def extract_chr_lrr_medians(vcf_file, chrom, threads=1):
    """Extract median LRR per sample for a specific chromosome.

    Returns dict of sample_index -> median_lrr.
    """
    # Determine correct chromosome name by checking what's in the VCF
    try:
        idx_output = subprocess.run(
            ['bcftools', 'index', '-s', vcf_file],
            capture_output=True, text=True
        )
        contigs = [line.split('\t')[0] for line in idx_output.stdout.strip().split('\n') if line]
    except Exception:
        contigs = []

    # Try chr-prefixed first, then non-prefixed
    target_chrom = None
    for candidate in [f'chr{chrom}', chrom]:
        if candidate in contigs:
            target_chrom = candidate
            break

    if target_chrom is None:
        return {}

    # Extract LRR values per sample for this chromosome
    cmd = (
        f'bcftools view -e "INFO/INTENSITY_ONLY=1" -r {target_chrom} '
        f'--threads {threads} "{vcf_file}" 2>/dev/null | '
        f'bcftools query -f "[\\t%LRR]\\n" 2>/dev/null'
    )
    try:
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = proc.stdout.strip().split('\n')
    except Exception:
        return {}

    if not lines or lines == ['']:
        return {}

    # Collect LRR values per sample column
    data = {}
    for line in lines:
        fields = line.split('\t')[1:]  # skip leading empty
        for i, val in enumerate(fields):
            if val in ('.', '', 'nan', '-nan', 'inf', '-inf'):
                continue
            try:
                v = float(val)
                data.setdefault(i, []).append(v)
            except ValueError:
                pass

    # Compute median per sample
    medians = {}
    for i, vals in data.items():
        vals.sort()
        n = len(vals)
        if n == 0:
            continue
        if n % 2 == 1:
            medians[i] = vals[n // 2]
        else:
            medians[i] = (vals[n // 2 - 1] + vals[n // 2]) / 2.0
    return medians


def read_sample_names(vcf_file):
    """Read sample names from VCF. Returns list of sample names."""
    try:
        proc = subprocess.run(
            ['bcftools', 'query', '-l', vcf_file],
            capture_output=True, text=True
        )
        return proc.stdout.strip().split('\n')
    except Exception:
        return []


def read_predicted_sex(sample_qc_file):
    """Read predicted sex from sample QC TSV. Returns dict of sample_id -> sex."""
    sex_map = {}
    with open(sample_qc_file) as f:
        header = f.readline().strip().split('\t')
        sex_idx = None
        for i, col in enumerate(header):
            if col == 'computed_gender':
                sex_idx = i
                break
        if sex_idx is None:
            return sex_map
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) > sex_idx:
                sex_map[fields[0]] = fields[sex_idx]
    return sex_map


def create_sex_check_plot(chrx_medians, chry_medians, sample_names,
                          sex_map, output_dir):
    """Create scatter plot of median chrX LRR vs median chrY LRR."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not available. Skipping sex check plot.",
              file=sys.stderr)
        return None

    # Gather data points with their sex labels
    points = []  # list of (x, y, sex_code)
    for i, name in enumerate(sample_names):
        if i in chrx_medians and i in chry_medians:
            sex = sex_map.get(name, 'NA')
            points.append((chrx_medians[i], chry_medians[i], sex))

    if not points:
        print("Warning: No samples with both chrX and chrY LRR data.",
              file=sys.stderr)
        return None

    fig, ax = plt.subplots(figsize=(8, 7))

    # Plot by sex category for legend (deduplicate Male codes 1/M, Female 2/F)
    sex_groups = {
        'Male': {'codes': {'1', 'M'}, 'color': '#4393C3'},
        'Female': {'codes': {'2', 'F'}, 'color': '#D6604D'},
        'Unknown': {'codes': {'NA', 'U', ''}, 'color': '#999999'},
    }
    for label, group in sex_groups.items():
        xi = [p[0] for p in points if p[2] in group['codes']]
        yi = [p[1] for p in points if p[2] in group['codes']]
        if xi:
            ax.scatter(xi, yi, c=group['color'], label=label, alpha=0.6,
                       s=20, edgecolors='none')

    ax.set_xlabel('Median LRR (chrX)', fontsize=12)
    ax.set_ylabel('Median LRR (chrY)', fontsize=12)
    ax.set_title('Sex Check: Median chrX vs chrY LRR', fontsize=14,
                 fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, 'sex_check_chrXY_lrr.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    return plot_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate sex check plot: median chrX LRR vs median chrY LRR"
    )
    parser.add_argument("--vcf", required=True,
                        help="Input VCF/BCF file with LRR FORMAT field")
    parser.add_argument("--sample-qc", required=True,
                        help="Sample QC TSV with computed_gender column")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for the plot")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads for bcftools (default: 1)")
    args = parser.parse_args()

    if not os.path.exists(args.vcf):
        print(f"Error: VCF file not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.sample_qc):
        print(f"Error: Sample QC file not found: {args.sample_qc}",
              file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    print("Generating sex check plot...")

    sample_names = read_sample_names(args.vcf)
    if not sample_names:
        print("Error: No samples found in VCF.", file=sys.stderr)
        sys.exit(1)

    print(f"  Extracting median chrX LRR ({len(sample_names)} samples)...")
    chrx_medians = extract_chr_lrr_medians(args.vcf, 'X', args.threads)

    print(f"  Extracting median chrY LRR ({len(sample_names)} samples)...")
    chry_medians = extract_chr_lrr_medians(args.vcf, 'Y', args.threads)

    sex_map = read_predicted_sex(args.sample_qc)

    plot_path = create_sex_check_plot(
        chrx_medians, chry_medians, sample_names, sex_map, args.output_dir
    )

    if plot_path:
        print(f"  Sex check plot: {plot_path}")
    else:
        print("  Warning: Sex check plot could not be generated.")


if __name__ == "__main__":
    main()
