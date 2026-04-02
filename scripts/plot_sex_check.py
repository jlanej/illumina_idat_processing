#!/usr/bin/env python3
"""
plot_sex_check.py

Generate a scatter plot of median per-sample LRR on chrX vs median per-sample
LRR on chrY, colored by predicted sex.  This is a standard sex-check
visualization for genotyping array QC.

LRR-based sex predictions are cross-tabulated with peddy's sex check
(when available) to flag discordances, potential aneuploidy, sample swaps,
or sex mismatches.

Males are expected to have lower chrX LRR (hemizygous) and higher chrY LRR,
while females should have higher chrX LRR and lower chrY LRR.

Usage:
    python3 plot_sex_check.py \\
        --vcf stage2/vcf/stage2_reclustered.bcf \\
        --sample-qc stage2/qc/stage2_sample_qc.tsv \\
        --output-dir stage2/qc \\
        [--peddy-sex-check peddy/peddy.sex_check.csv]
"""

import argparse
import csv
import os
import subprocess
import sys

# Allow importing par_xtr_regions from the same directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from par_xtr_regions import get_par_xtr_bed  # noqa: E402


def _get_par_xtr_bed_path(genome, output_dir=None):
    """Return path to a PAR/XTR exclusion BED for the given genome build.

    When *output_dir* is provided the BED is cached there so it is written
    once per run.  When *output_dir* is ``None`` a temporary file is created
    (caller should clean up or rely on OS temp-file reaping).
    """
    if output_dir:
        bed_path = os.path.join(output_dir, 'par_xtr_exclusion.bed')
        if not os.path.exists(bed_path):
            get_par_xtr_bed(genome, bed_path)
        return bed_path
    # No output_dir — write to a temp file so PAR/XTR exclusion always works.
    return get_par_xtr_bed(genome)


def extract_chrXY_lrr_medians(vcf_file, threads=1, genome='CHM13',
                               output_dir=None):
    """Extract median LRR per sample for chrX and chrY in a SINGLE VCF pass.

    Instead of two separate passes (one for chrX, one for chrY), this
    extracts both chromosomes in one bcftools view -r chrX,chrY call,
    then splits by chromosome in Python. Halves I/O.

    PAR1, PAR2, and XTR regions are excluded from both chrX and chrY via
    ``-T ^bed`` so that the medians reflect only non-PAR/XTR loci.

    If a cached result exists at ``output_dir/chrXY_lrr_cache.tsv`` it is
    returned directly, avoiding the expensive BCF scan on the second
    invocation (post-peddy).

    Returns (dict of sample_index -> median_lrr_X, dict of sample_index -> median_lrr_Y).
    """
    # Check for cached results
    cache_path = None
    if output_dir:
        cache_path = os.path.join(output_dir, 'chrXY_lrr_cache.tsv')
        if os.path.exists(cache_path):
            cached = _read_lrr_cache(cache_path, expected_genome=genome)
            if cached is not None:
                return cached

    # Determine correct chromosome names by checking what's in the VCF
    try:
        idx_output = subprocess.run(
            ['bcftools', 'index', '-s', vcf_file],
            capture_output=True, text=True
        )
        contigs = [line.split('\t')[0] for line in idx_output.stdout.strip().split('\n') if line]
    except Exception:
        contigs = []

    # Find chrX and chrY names
    chrom_x = None
    chrom_y = None
    for candidate in ['chrX', 'X']:
        if candidate in contigs:
            chrom_x = candidate
            break
    for candidate in ['chrY', 'Y']:
        if candidate in contigs:
            chrom_y = candidate
            break

    if chrom_x is None and chrom_y is None:
        return {}, {}

    # Build region string for both chromosomes
    regions = ','.join(c for c in [chrom_x, chrom_y] if c is not None)

    # Write PAR/XTR exclusion BED — always generated so that chrX and chrY
    # LRR medians reflect only non-PAR/XTR loci.  The BED contains both
    # chrX and chrY entries, so ``-T ^bed`` with ``-r chrX,chrY`` correctly
    # excludes PAR/XTR from both chromosomes in a single pass.
    bed_path = _get_par_xtr_bed_path(genome, output_dir)
    par_xtr_arg = f' -T ^"{bed_path}"'

    # Single pass: extract CHROM and LRR per sample
    cmd = (
        f'bcftools view -e "INFO/INTENSITY_ONLY=1" -r {regions}'
        f'{par_xtr_arg} '
        f'--threads {threads} "{vcf_file}" 2>/dev/null | '
        f'bcftools query -f "%CHROM[\\t%LRR]\\n" 2>/dev/null'
    )
    try:
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = proc.stdout.strip().split('\n')
    except Exception:
        return {}, {}

    if not lines or lines == ['']:
        return {}, {}

    # Collect LRR values per sample column, split by chromosome
    data_x = {}  # sample_idx -> list of float
    data_y = {}
    for line in lines:
        fields = line.split('\t')
        if len(fields) < 2:
            continue
        chrom = fields[0]
        for i, val in enumerate(fields[1:]):
            if val in ('.', '', 'nan', '-nan', 'inf', '-inf'):
                continue
            try:
                v = float(val)
                if chrom == chrom_x:
                    data_x.setdefault(i, []).append(v)
                elif chrom == chrom_y:
                    data_y.setdefault(i, []).append(v)
            except ValueError:
                pass

    # Compute median per sample
    def _compute_medians(data):
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

    medians_x = _compute_medians(data_x)
    medians_y = _compute_medians(data_y)

    # Cache for reuse within the same run (second invocation post-peddy)
    if cache_path:
        _write_lrr_cache(cache_path, medians_x, medians_y, genome)

    return medians_x, medians_y


def _write_lrr_cache(cache_path, medians_x, medians_y, genome='CHM13'):
    """Write LRR medians to a TSV cache file."""
    with open(cache_path, 'w') as f:
        f.write(f'# genome={genome}\n')
        f.write('sample_idx\tchrx_lrr_median\tchry_lrr_median\n')
        for idx in sorted(set(medians_x.keys()) | set(medians_y.keys())):
            x_val = f'{medians_x[idx]:.6f}' if idx in medians_x else 'NA'
            y_val = f'{medians_y[idx]:.6f}' if idx in medians_y else 'NA'
            f.write(f'{idx}\t{x_val}\t{y_val}\n')


def _read_lrr_cache(cache_path, expected_genome=None):
    """Read LRR medians from a TSV cache file.

    Returns (medians_x, medians_y) or None on genome mismatch.
    """
    medians_x, medians_y = {}, {}
    with open(cache_path) as f:
        first_line = f.readline().strip()
        if first_line.startswith('# genome='):
            cached_genome = first_line.split('=', 1)[1]
            if expected_genome and cached_genome != expected_genome:
                print(f"  LRR cache genome mismatch ({cached_genome} vs "
                      f"{expected_genome}); recomputing.", file=sys.stderr)
                return None
            f.readline()  # skip header
        else:
            # No genome comment; first_line is the header — skip it.
            pass
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            try:
                idx = int(fields[0])
            except ValueError:
                continue
            if fields[1] != 'NA':
                medians_x[idx] = float(fields[1])
            if fields[2] != 'NA':
                medians_y[idx] = float(fields[2])
    return medians_x, medians_y


def read_peddy_sex_check(peddy_sex_check_file):
    """Read peddy sex_check.csv and return dict of sample_id -> peddy sex info.

    Returns dict of sample_id -> {'peddy_sex': str, 'error': bool}.
    """
    results = {}
    if not peddy_sex_check_file or not os.path.exists(peddy_sex_check_file):
        return results
    try:
        with open(peddy_sex_check_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                sid = row.get('sample_id', '')
                if not sid:
                    continue
                predicted = row.get('predicted_sex', 'NA')
                error = row.get('error', 'False').strip().lower() == 'true'
                results[sid] = {
                    'peddy_sex': predicted,
                    'error': error,
                }
    except Exception:
        pass
    return results


def cross_tabulate_sex(sample_names, sex_map, chrx_medians, chry_medians,
                       peddy_sex_results):
    """Cross-tabulate sex determinations from multiple methods.

    For each sample, compare:
      1. LRR-based sex (computed_gender from GTC metadata)
      2. Peddy sex (genotype-based prediction, if available)

    Status values:
      CONCORDANT   – ≥2 methods agree on M or F
      DISCORDANT   – ≥2 methods disagree
      SINGLE_METHOD – only 1 method produced a definitive M/F call
      UNDETERMINED – no method produced a definitive call

    Returns list of dicts, one per sample.
    """
    _normalize = {'1': 'M', 'M': 'M', '2': 'F', 'F': 'F',
                  'male': 'M', 'female': 'F'}

    rows = []
    for i, name in enumerate(sample_names):
        if i not in chrx_medians and i not in chry_medians:
            continue

        lrr_sex = _normalize.get(sex_map.get(name, ''), 'NA')

        peddy_info = peddy_sex_results.get(name, {})
        peddy_sex = _normalize.get(peddy_info.get('peddy_sex', ''), 'NA')
        peddy_error = peddy_info.get('error', False)

        # Determine concordance using LRR-based and peddy sex only.
        methods = []
        if lrr_sex in ('M', 'F'):
            methods.append(lrr_sex)
        if peddy_sex in ('M', 'F'):
            methods.append(peddy_sex)

        if len(methods) >= 2 and len(set(methods)) > 1:
            status = 'DISCORDANT'
        elif len(methods) >= 2 and len(set(methods)) == 1:
            status = 'CONCORDANT'
        elif len(methods) == 1:
            status = 'SINGLE_METHOD'
        else:
            status = 'UNDETERMINED'

        rows.append({
            'sample_id': name,
            'lrr_sex': lrr_sex,
            'peddy_sex': peddy_sex,
            'peddy_error': peddy_error,
            'status': status,
            'chrx_lrr_median': chrx_medians.get(i),
            'chry_lrr_median': chry_medians.get(i),
        })

    return rows


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
    except (ImportError, AttributeError):
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


def write_sex_check_table(chrx_medians, chry_medians, sample_names, sex_map,
                          output_dir, peddy_sex_results=None):
    """Write per-sample sex check table with all available methods.

    Columns include LRR-based and peddy sex determinations along with a
    concordance status flag.
    """
    if peddy_sex_results is None:
        peddy_sex_results = {}

    cross_tab = cross_tabulate_sex(
        sample_names, sex_map, chrx_medians, chry_medians,
        peddy_sex_results,
    )

    table_path = os.path.join(output_dir, 'sex_check_chrXY_lrr.tsv')
    with open(table_path, 'w') as out:
        out.write("sample_id\tchrx_lrr_median\tchry_lrr_median\t"
                  "computed_gender\t"
                  "peddy_sex\tsex_status\n")
        for row in cross_tab:
            chrx = f'{row["chrx_lrr_median"]:.6f}' if row['chrx_lrr_median'] is not None else 'NA'
            chry = f'{row["chry_lrr_median"]:.6f}' if row['chry_lrr_median'] is not None else 'NA'
            out.write(f"{row['sample_id']}\t{chrx}\t{chry}\t"
                      f"{row['lrr_sex']}\t"
                      f"{row['peddy_sex']}\t{row['status']}\n")

    # Write discordance summary
    n_discordant = sum(1 for r in cross_tab if r['status'] == 'DISCORDANT')
    n_concordant = sum(1 for r in cross_tab if r['status'] == 'CONCORDANT')
    n_single = sum(1 for r in cross_tab if r['status'] == 'SINGLE_METHOD')
    n_undetermined = sum(1 for r in cross_tab if r['status'] == 'UNDETERMINED')

    summary_path = os.path.join(output_dir, 'sex_check_summary.txt')
    with open(summary_path, 'w') as out:
        out.write("Sex Check Cross-Tabulation Summary\n")
        out.write("=" * 50 + "\n\n")
        out.write(f"Total samples evaluated:       {len(cross_tab)}\n")
        out.write(f"Concordant (≥2 methods agree): {n_concordant}\n")
        out.write(f"Discordant (methods disagree): {n_discordant}\n")
        out.write(f"Single method only:            {n_single}\n")
        out.write(f"Undetermined (no calls):       {n_undetermined}\n\n")
        out.write("Status definitions:\n")
        out.write("  CONCORDANT   – ≥2 methods agree on M or F\n")
        out.write("  DISCORDANT   – ≥2 methods disagree (possible "
                  "sample swap / aneuploidy)\n")
        out.write("  SINGLE_METHOD – only 1 method produced a "
                  "definitive M/F call\n")
        out.write("  UNDETERMINED – no method produced a definitive "
                  "call\n\n")
        out.write("Methods used:\n")
        out.write("  1. LRR-based (computed_gender): sex from GTC "
                  "metadata / intensity separation\n")
        if peddy_sex_results:
            out.write("  2. Peddy: Genotype-based sex prediction via "
                      "peddy\n")
        out.write("\nDiscordant samples may indicate:\n")
        out.write("  - Sex chromosome aneuploidy (e.g. XXY, X0)\n")
        out.write("  - Sample swaps or contamination\n")
        out.write("  - Data quality issues\n")

    return table_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate sex check plot and cross-tabulate LRR "
                    "and peddy sex determinations"
    )
    parser.add_argument("--vcf", required=True,
                        help="Input VCF/BCF file with LRR FORMAT field")
    parser.add_argument("--sample-qc", required=True,
                        help="Sample QC TSV with computed_gender column")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for the plot")
    parser.add_argument("--peddy-sex-check", default=None,
                        help="Peddy sex_check.csv for cross-tabulation "
                             "(optional)")
    parser.add_argument("--genome", default='CHM13',
                        help="Genome build for PAR/XTR exclusion: "
                             "CHM13, GRCh38, or GRCh37 (default: CHM13)")
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
    print(f"  Genome build: {args.genome} (PAR/XTR regions excluded from chrX)")

    sample_names = read_sample_names(args.vcf)
    if not sample_names:
        print("Error: No samples found in VCF.", file=sys.stderr)
        sys.exit(1)

    print(f"  Extracting median chrX + chrY LRR ({len(sample_names)} samples, single pass)...")
    chrx_medians, chry_medians = extract_chrXY_lrr_medians(
        args.vcf, args.threads, genome=args.genome,
        output_dir=args.output_dir)

    sex_map = read_predicted_sex(args.sample_qc)

    peddy_sex_results = read_peddy_sex_check(args.peddy_sex_check)
    if peddy_sex_results:
        print(f"  Peddy sex check: {len(peddy_sex_results)} samples loaded")

    plot_path = create_sex_check_plot(
        chrx_medians, chry_medians, sample_names, sex_map, args.output_dir
    )
    table_path = write_sex_check_table(
        chrx_medians, chry_medians, sample_names, sex_map, args.output_dir,
        peddy_sex_results=peddy_sex_results,
    )

    if plot_path:
        print(f"  Sex check plot: {plot_path}")
    else:
        print("  Warning: Sex check plot could not be generated.")
    if table_path:
        print(f"  Sex check table: {table_path}")
    summary_path = os.path.join(args.output_dir, 'sex_check_summary.txt')
    if os.path.exists(summary_path):
        print(f"  Sex check summary: {summary_path}")


if __name__ == "__main__":
    main()
