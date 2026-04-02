#!/usr/bin/env python3
"""
plot_sex_check.py

Generate a scatter plot of median per-sample LRR on chrX vs median per-sample
LRR on chrY, colored by predicted sex.  This is a standard sex-check
visualization for genotyping array QC.

Additionally computes X-chromosome F-statistic (genotype-based inbreeding
coefficient on chrX) for a complementary genotype-based sex determination.
Males (hemizygous X) have F close to 1.0; females (diploid X) have F close
to 0.0.  Results from both methods (LRR-based and F-statistic) are
cross-tabulated with peddy's sex check (when available) to flag discordances,
potential aneuploidy, sample swaps, or sex mismatches.

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
import shutil
import subprocess
import sys
import tempfile

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


def compute_chrx_f_statistic(vcf_file, sample_names, threads=1,
                             genome='CHM13', output_dir=None):
    """Compute X-chromosome inbreeding F-statistic per sample.

    F = 1 - (observed_het / expected_het).  For chrX:
      - Males (hemizygous): F ≈ 1.0  (very few heterozygous calls)
      - Females (diploid):  F ≈ 0.0  (normal heterozygosity)

    Thresholds (following plink --check-sex convention):
      - F > 0.8  → likely male
      - F < 0.2  → likely female
      - 0.2 ≤ F ≤ 0.8 → ambiguous (potential aneuploidy or data issue)

    PAR1, PAR2, and XTR regions are excluded so that the F-statistic
    reflects only non-pseudoautosomal chrX genotypes.

    When plink2 is available it is used for an allele-frequency-aware
    F-statistic; otherwise a simplified bcftools-based approximation is
    used as a fallback.

    If a cached result file exists at ``output_dir/chrx_f_stat_cache.tsv``,
    results are read from cache instead of recomputing.

    Returns dict of sample_name -> {'f_stat': float, 'n_het': int,
    'n_called': int, 'f_sex': str}.
    """
    # Check for cached results (cache includes genome build for invalidation)
    cache_path = None
    if output_dir:
        cache_path = os.path.join(output_dir, 'chrx_f_stat_cache.tsv')
        if os.path.exists(cache_path):
            cached = _read_f_stat_cache(cache_path, expected_genome=genome)
            if cached is not None:
                return cached

    if shutil.which('plink2'):
        results = _compute_f_stat_plink2(
            vcf_file, sample_names, threads, genome, output_dir)
    else:
        results = _compute_f_stat_bcftools(
            vcf_file, sample_names, threads, genome, output_dir)

    # Write cache for reuse within the same run
    if cache_path and results:
        _write_f_stat_cache(cache_path, results, genome)

    return results


def _compute_f_stat_plink2(vcf_file, sample_names, threads, genome,
                           output_dir):
    """Compute chrX F-statistic via plink2 --check-sex.

    PAR/XTR regions are excluded via ``--exclude-range``.  An alternative
    approach would be plink2's ``--split-par`` flag, which reclassifies PAR
    variants as autosomal — but ``--exclude-range`` is more explicit and
    also handles XTR, which ``--split-par`` does not cover.

    Note: some newer plink2 alpha builds renamed ``--check-sex`` to
    ``--sex-check``.  We try ``--check-sex`` first (compatible with the
    Dockerfile-pinned v2.00a6.33 and most stable builds) and fall back to
    ``--sex-check`` if the first attempt fails with an unrecognized-flag
    error.
    """
    work_dir = tempfile.mkdtemp(prefix='plink2_sexcheck_')
    try:
        prefix = os.path.join(work_dir, 'sexcheck')
        bed_path = _get_par_xtr_bed_path(genome, output_dir)
        base_cmd = [
            'plink2', '--bcf', vcf_file,
            '--allow-extra-chr', '--threads', str(threads),
            '--exclude-range', bed_path,
        ]

        # Try --check-sex first; fall back to --sex-check for newer alphas
        for sex_flag in ['--check-sex', '--sex-check']:
            cmd = base_cmd + [sex_flag, '0.2', '0.8', '--out', prefix]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            sexcheck_file = prefix + '.sexcheck'
            if proc.returncode == 0 and os.path.exists(sexcheck_file):
                break
            # If the error is specifically about an unrecognized flag name,
            # try the next flag variant.  plink2 uses phrases like
            # "Unrecognized flag" or "Invalid flag" in its error output.
            stderr_lower = proc.stderr.lower()
            if 'unrecognized flag' in stderr_lower or \
               'invalid flag' in stderr_lower:
                continue
            # Other error — report and return empty
            if proc.stderr:
                print(f"  plink2 failed (rc={proc.returncode}): "
                      f"{proc.stderr.strip()[:200]}", file=sys.stderr)
            return {}
        else:
            # Neither flag worked
            if proc.stderr:
                print(f"  plink2 failed (rc={proc.returncode}): "
                      f"{proc.stderr.strip()[:200]}", file=sys.stderr)
            return {}

        # Copy raw .sexcheck to output_dir for diagnostics
        if output_dir:
            shutil.copy2(sexcheck_file,
                         os.path.join(output_dir, 'plink2_sexcheck.tsv'))

        # Parse .sexcheck: columns are #FID IID PEDSEX SNPSEX STATUS F
        # plink2 may also provide OBS_CT (non-missing genotype count).
        results = {}
        sample_set = set(sample_names)
        with open(sexcheck_file) as f:
            header = f.readline().strip().split()
            col_map = {col.lstrip('#'): i for i, col in enumerate(header)}
            iid_idx = col_map.get('IID', 1)
            f_idx = col_map.get('F', len(header) - 1)
            obs_ct_idx = col_map.get('OBS_CT')
            for line in f:
                fields = line.strip().split()
                if len(fields) <= max(iid_idx, f_idx):
                    continue
                name = fields[iid_idx]
                if name not in sample_set:
                    continue
                try:
                    f_stat = float(fields[f_idx])
                except ValueError:
                    continue
                if f_stat > 0.8:
                    f_sex = 'M'
                elif f_stat < 0.2:
                    f_sex = 'F'
                else:
                    f_sex = 'ambiguous'

                # Extract OBS_CT if present in the header
                n_called = None
                if obs_ct_idx is not None and obs_ct_idx < len(fields):
                    try:
                        n_called = int(fields[obs_ct_idx])
                    except ValueError:
                        pass

                results[name] = {
                    'f_stat': round(f_stat, 6),
                    'n_het': None,
                    'n_called': n_called,
                    'f_sex': f_sex,
                }
        return results
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def _compute_f_stat_bcftools(vcf_file, sample_names, threads, genome,
                             output_dir):
    """Compute chrX F-statistic via bcftools genotype counting (fallback).

    This is a simplified approximation used when plink2 is not available.
    It estimates F ≈ 1 - 2 * observed_het_rate rather than computing a
    proper allele-frequency-aware inbreeding coefficient.  The result is
    less accurate than plink2 --check-sex but sufficient for flagging
    obvious sex mismatches on genotyping array data.
    """
    # Determine chrX contig name
    try:
        idx_output = subprocess.run(
            ['bcftools', 'index', '-s', vcf_file],
            capture_output=True, text=True
        )
        contigs = [line.split('\t')[0]
                   for line in idx_output.stdout.strip().split('\n') if line]
    except Exception:
        contigs = []

    chrom_x = None
    for candidate in ['chrX', 'X']:
        if candidate in contigs:
            chrom_x = candidate
            break
    if chrom_x is None:
        return {}

    # PAR/XTR exclusion BED — always generated so chrX F-stat reflects
    # only non-PAR/XTR loci.
    bed_path = _get_par_xtr_bed_path(genome, output_dir)
    par_xtr_arg = f' -T ^"{bed_path}"'

    # Extract genotypes on chrX (biallelic SNPs only, non-intensity-only)
    cmd = (
        f'bcftools view -e "INFO/INTENSITY_ONLY=1" -r {chrom_x}'
        f'{par_xtr_arg} '
        f'-v snps -m2 -M2 --threads {threads} "{vcf_file}" 2>/dev/null | '
        f'bcftools query -f "[%GT\\t]\\n" 2>/dev/null'
    )
    try:
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = proc.stdout.strip().split('\n')
    except Exception:
        return {}

    if not lines or lines == ['']:
        return {}

    n_samples = len(sample_names)
    het_counts = [0] * n_samples
    called_counts = [0] * n_samples

    for line in lines:
        gts = line.rstrip('\t').split('\t')
        for i, gt in enumerate(gts):
            if i >= n_samples:
                break
            if gt in ('.', './.', '.|.', ''):
                continue
            called_counts[i] += 1
            # Heterozygous: 0/1, 0|1, 1/0, 1|0
            alleles = gt.replace('|', '/').split('/')
            if len(alleles) == 2 and alleles[0] != alleles[1]:
                het_counts[i] += 1

    # Simplified F-statistic approximation for sex chromosome ploidy.
    # The standard plink --check-sex computes F = 1 - (obs_het / E[het])
    # where E[het] is the per-locus expected heterozygosity from allele
    # frequencies summed across all loci.  This fallback uses a proxy:
    #   F ≈ 1 - 2 * observed_het_rate
    # which maps het_rate ∈ [0, 0.5] → F ∈ [1, 0].
    #   Males (hemizygous, het_rate ≈ 0)    → F ≈ 1.0  (classified M)
    #   Females (diploid, het_rate ≈ 0.4+)  → F ≈ 0.2  (classified F)
    # Note: this approximation is less accurate than plink2's allele-
    # frequency-aware F, particularly for datasets with many low-MAF
    # variants where expected heterozygosity differs from the uniform
    # assumption.  Prefer plink2 --check-sex when available.
    results = {}
    for i, name in enumerate(sample_names):
        n_called = called_counts[i]
        n_het = het_counts[i]
        if n_called == 0:
            continue
        het_rate = n_het / n_called
        f_stat = 1.0 - (2.0 * het_rate) if het_rate <= 0.5 else 0.0

        if f_stat > 0.8:
            f_sex = 'M'
        elif f_stat < 0.2:
            f_sex = 'F'
        else:
            f_sex = 'ambiguous'

        results[name] = {
            'f_stat': round(f_stat, 6),
            'n_het': n_het,
            'n_called': n_called,
            'f_sex': f_sex,
        }

    return results


def _write_f_stat_cache(cache_path, results, genome='CHM13'):
    """Write F-statistic results to a TSV cache file."""
    with open(cache_path, 'w') as f:
        f.write(f'# genome={genome}\n')
        f.write('sample_id\tf_stat\tn_het\tn_called\tf_sex\n')
        for name, info in results.items():
            n_het = info['n_het'] if info['n_het'] is not None else 'NA'
            n_called = info['n_called'] if info['n_called'] is not None else 'NA'
            f.write(f"{name}\t{info['f_stat']}\t{n_het}\t"
                    f"{n_called}\t{info['f_sex']}\n")


def _read_f_stat_cache(cache_path, expected_genome=None):
    """Read F-statistic results from a TSV cache file.

    Returns None if the cache genome build does not match expected_genome.
    """
    results = {}
    with open(cache_path) as f:
        first_line = f.readline().strip()
        if first_line.startswith('# genome='):
            cached_genome = first_line.split('=', 1)[1]
            if expected_genome and cached_genome != expected_genome:
                print(f"  Cache genome mismatch ({cached_genome} vs "
                      f"{expected_genome}); recomputing.", file=sys.stderr)
                return None
            header = f.readline()  # skip header
        else:
            header = first_line  # no genome comment; treat as header
        for lineno, line in enumerate(f, start=3):
            fields = line.strip().split('\t')
            if len(fields) < 5:
                print(f"  Warning: skipping malformed cache line {lineno}: "
                      f"{line.strip()}", file=sys.stderr)
                continue
            results[fields[0]] = {
                'f_stat': float(fields[1]),
                'n_het': int(fields[2]) if fields[2] != 'NA' else None,
                'n_called': int(fields[3]) if fields[3] != 'NA' else None,
                'f_sex': fields[4],
            }
    return results


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
                       f_stat_results, peddy_sex_results):
    """Cross-tabulate sex determinations from multiple methods.

    For each sample, compare:
      1. LRR-based sex (from computed_gender in sample QC)
      2. F-statistic sex (from chrX genotype heterozygosity)
      3. Peddy sex (if available)

    Flag discordances and potential issues.

    Status values:
      CONCORDANT   – ≥2 methods agree on M or F
      DISCORDANT   – ≥2 methods disagree
      AMBIGUOUS    – chrX F-statistic in the 0.2–0.8 indeterminate zone
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

        f_info = f_stat_results.get(name, {})
        raw_f_sex = f_info.get('f_sex', '')
        f_sex = _normalize.get(raw_f_sex, 'NA')
        f_stat = f_info.get('f_stat', None)
        # Display value preserves 'ambiguous' so users can distinguish
        # an indeterminate F-stat from a missing one (plain 'NA').
        f_sex_display = raw_f_sex if raw_f_sex else 'NA'

        peddy_info = peddy_sex_results.get(name, {})
        peddy_sex = _normalize.get(peddy_info.get('peddy_sex', ''), 'NA')
        peddy_error = peddy_info.get('error', False)

        # Determine concordance
        methods = []
        if lrr_sex in ('M', 'F'):
            methods.append(lrr_sex)
        if f_sex in ('M', 'F'):
            methods.append(f_sex)
        if peddy_sex in ('M', 'F'):
            methods.append(peddy_sex)

        if len(methods) >= 2 and len(set(methods)) > 1:
            status = 'DISCORDANT'
        elif raw_f_sex == 'ambiguous':
            status = 'AMBIGUOUS'
        elif len(methods) >= 2 and len(set(methods)) == 1:
            status = 'CONCORDANT'
        elif len(methods) == 1:
            status = 'SINGLE_METHOD'
        else:
            status = 'UNDETERMINED'

        rows.append({
            'sample_id': name,
            'lrr_sex': lrr_sex,
            'f_stat': f_stat,
            'f_sex': f_sex_display,
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
                          output_dir, f_stat_results=None,
                          peddy_sex_results=None):
    """Write per-sample sex check table with all available methods.

    Columns include LRR-based, F-statistic, and peddy sex determinations
    along with a concordance status flag.
    """
    if f_stat_results is None:
        f_stat_results = {}
    if peddy_sex_results is None:
        peddy_sex_results = {}

    cross_tab = cross_tabulate_sex(
        sample_names, sex_map, chrx_medians, chry_medians,
        f_stat_results, peddy_sex_results,
    )

    table_path = os.path.join(output_dir, 'sex_check_chrXY_lrr.tsv')
    with open(table_path, 'w') as out:
        out.write("sample_id\tchrx_lrr_median\tchry_lrr_median\t"
                  "computed_gender\tchrx_f_stat\tf_sex\t"
                  "peddy_sex\tsex_status\n")
        for row in cross_tab:
            f_val = f'{row["f_stat"]:.6f}' if row['f_stat'] is not None else 'NA'
            chrx = f'{row["chrx_lrr_median"]:.6f}' if row['chrx_lrr_median'] is not None else 'NA'
            chry = f'{row["chry_lrr_median"]:.6f}' if row['chry_lrr_median'] is not None else 'NA'
            out.write(f"{row['sample_id']}\t{chrx}\t{chry}\t"
                      f"{row['lrr_sex']}\t{f_val}\t{row['f_sex']}\t"
                      f"{row['peddy_sex']}\t{row['status']}\n")

    # Write discordance summary
    n_discordant = sum(1 for r in cross_tab if r['status'] == 'DISCORDANT')
    n_ambiguous = sum(1 for r in cross_tab if r['status'] == 'AMBIGUOUS')
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
        out.write(f"Ambiguous (F-stat 0.2-0.8):    {n_ambiguous}\n")
        out.write(f"Single method only:            {n_single}\n")
        out.write(f"Undetermined (no calls):       {n_undetermined}\n\n")
        out.write("Status definitions:\n")
        out.write("  CONCORDANT   – ≥2 methods agree on M or F\n")
        out.write("  DISCORDANT   – ≥2 methods disagree (possible "
                  "sample swap / aneuploidy)\n")
        out.write("  AMBIGUOUS    – chrX F-statistic in the 0.2–0.8 "
                  "indeterminate zone\n")
        out.write("  SINGLE_METHOD – only 1 method produced a "
                  "definitive M/F call\n")
        out.write("  UNDETERMINED – no method produced a definitive "
                  "call\n\n")
        out.write("Methods used:\n")
        out.write("  1. LRR-based (computed_gender): sex from GTC "
                  "metadata / intensity separation\n")
        out.write("  2. F-statistic: chrX genotype heterozygosity "
                  "(F>0.8=M, F<0.2=F, PAR/XTR excluded)\n")
        if peddy_sex_results:
            out.write("  3. Peddy: Genotype-based sex prediction via "
                      "peddy\n")
        out.write("\nColumn notes:\n")
        out.write("  f_sex values: M, F, ambiguous (F in 0.2–0.8), "
                  "or NA (not computed)\n")
        out.write("\nDiscordant samples may indicate:\n")
        out.write("  - Sex chromosome aneuploidy (e.g. XXY, X0)\n")
        out.write("  - Sample swaps or contamination\n")
        out.write("  - Data quality issues\n")

    return table_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate sex check plot and cross-tabulate LRR, "
                    "F-statistic, and peddy sex determinations"
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
    if shutil.which('plink2'):
        print("  Using plink2 --check-sex for allele-frequency-aware F-statistic")
    else:
        print("  plink2 not found; using bcftools-based F-statistic approximation")

    sample_names = read_sample_names(args.vcf)
    if not sample_names:
        print("Error: No samples found in VCF.", file=sys.stderr)
        sys.exit(1)

    print(f"  Extracting median chrX + chrY LRR ({len(sample_names)} samples, single pass)...")
    chrx_medians, chry_medians = extract_chrXY_lrr_medians(
        args.vcf, args.threads, genome=args.genome,
        output_dir=args.output_dir)

    sex_map = read_predicted_sex(args.sample_qc)

    print("  Computing chrX F-statistic for genotype-based sex check...")
    f_stat_results = compute_chrx_f_statistic(
        args.vcf, sample_names, args.threads,
        genome=args.genome, output_dir=args.output_dir)
    if f_stat_results:
        n_m = sum(1 for v in f_stat_results.values() if v['f_sex'] == 'M')
        n_f = sum(1 for v in f_stat_results.values() if v['f_sex'] == 'F')
        n_a = sum(1 for v in f_stat_results.values() if v['f_sex'] == 'ambiguous')
        print(f"    F-statistic sex: {n_m} male, {n_f} female, {n_a} ambiguous")
    else:
        print("    Warning: Could not compute chrX F-statistic.")

    peddy_sex_results = read_peddy_sex_check(args.peddy_sex_check)
    if peddy_sex_results:
        print(f"  Peddy sex check: {len(peddy_sex_results)} samples loaded")

    plot_path = create_sex_check_plot(
        chrx_medians, chry_medians, sample_names, sex_map, args.output_dir
    )
    table_path = write_sex_check_table(
        chrx_medians, chry_medians, sample_names, sex_map, args.output_dir,
        f_stat_results=f_stat_results,
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
