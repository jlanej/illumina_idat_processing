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

When a sex-chromosome QC directory is provided (--sex-chr-qc-dir), chrX
and chrY variant metrics are also included in the output.  Sex-chromosome
QC follows best practices (Laurie et al. 2012, Genet Epidemiol 36:384-91):
  - chrX HWE is computed on females only (diploid genotypes)
  - chrX call rate / MAF are reported for all samples, plus per-sex
  - chrY metrics are computed on males only
  - PAR/XTR regions are excluded

Usage:
    python3 collate_variant_qc.py \\
        --all-variant-qc-dir stage2/qc/variant_qc \\
        --ancestry-variant-qc EUR:ancestry_qc/EUR/variant_qc \\
        --ancestry-variant-qc AFR:ancestry_qc/AFR/variant_qc \\
        --sex-chr-qc-dir ancestry_qc/sex_chr_qc \\
        --ancestry-sex-chr-qc EUR:ancestry_qc/EUR/sex_chr_qc \\
        --ancestry-sex-chr-qc AFR:ancestry_qc/AFR/sex_chr_qc \\
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
        fmiss = vmiss_row.get('F_MISS', 'NA')
        try:
            cr = f'{1 - float(fmiss):.6f}' if fmiss != 'NA' else 'NA'
        except ValueError:
            cr = 'NA'

        hardy_row = hardy.get(vid, {})

        data[vid] = {
            'call_rate': cr,
            'hwe_p': hardy_row.get('P', 'NA'),
            'maf': afreq.get(vid, 'NA'),
            'obs_ct': vmiss_row.get('OBS_CT', 'NA'),
            'missing_ct': vmiss_row.get('MISSING_CT', 'NA'),
            'hom_a1_ct': hardy_row.get('HOM_A1_CT', 'NA'),
            'het_ct': hardy_row.get('HET_A1_AX_CT', 'NA'),
            'hom_a2_ct': hardy_row.get('TWO_AX_CT', 'NA'),
        }
    return data


def _fmiss_to_call_rate(vmiss_row):
    """Convert a plink2 vmiss row to a formatted call-rate string.

    Returns 'NA' if F_MISS is missing or unparseable.
    """
    fmiss = vmiss_row.get('F_MISS', 'NA')
    try:
        return f'{1 - float(fmiss):.6f}' if fmiss != 'NA' else 'NA'
    except (TypeError, ValueError):
        return 'NA'


def load_sex_chr_qc(sex_chr_dir):
    """Load sex-chromosome variant QC from the structured output directory.

    Expected layout (produced by compute_sex_chr_variant_qc.sh):
      sex_chr_dir/chrX/chrX_all.vmiss     - all-sample missingness
      sex_chr_dir/chrX/chrX_all.afreq     - all-sample allele frequency
      sex_chr_dir/chrX/chrX_female.hardy   - female-only HWE
      sex_chr_dir/chrX/chrX_female.vmiss   - female-only missingness
      sex_chr_dir/chrX/chrX_male.vmiss     - male-only missingness
      sex_chr_dir/chrY/chrY_male.vmiss     - male-only missingness
      sex_chr_dir/chrY/chrY_male.afreq     - male-only allele frequency

    Returns a dict with keys 'chrX' and 'chrY', each mapping variant_id
    to a dict of QC metrics with column-name-ready keys.
    """
    result = {'chrX': {}, 'chrY': {}}

    if not sex_chr_dir or not os.path.isdir(sex_chr_dir):
        return result

    # ---- chrX ----
    chrx_dir = os.path.join(sex_chr_dir, 'chrX')
    if os.path.isdir(chrx_dir):
        # All-sample call rate & MAF
        chrx_all_vmiss = read_plink2_vmiss(
            os.path.join(chrx_dir, 'chrX_all.vmiss'))
        chrx_all_afreq = read_plink2_afreq(
            os.path.join(chrx_dir, 'chrX_all.afreq'))
        # Female-only HWE & missingness
        chrx_female_hardy = read_plink2_hardy(
            os.path.join(chrx_dir, 'chrX_female.hardy'))
        chrx_female_vmiss = read_plink2_vmiss(
            os.path.join(chrx_dir, 'chrX_female.vmiss'))
        # Male-only missingness
        chrx_male_vmiss = read_plink2_vmiss(
            os.path.join(chrx_dir, 'chrX_male.vmiss'))

        all_chrx_ids = sorted(set(
            list(chrx_all_vmiss.keys()) +
            list(chrx_all_afreq.keys()) +
            list(chrx_female_hardy.keys()) +
            list(chrx_female_vmiss.keys()) +
            list(chrx_male_vmiss.keys())
        ))

        for vid in all_chrx_ids:
            fem_hardy = chrx_female_hardy.get(vid, {})
            result['chrX'][vid] = {
                'chrX_call_rate': _fmiss_to_call_rate(
                    chrx_all_vmiss.get(vid, {})),
                'chrX_maf': chrx_all_afreq.get(vid, 'NA'),
                'chrX_female_call_rate': _fmiss_to_call_rate(
                    chrx_female_vmiss.get(vid, {})),
                'chrX_female_hwe_p': fem_hardy.get('P', 'NA'),
                'chrX_male_call_rate': _fmiss_to_call_rate(
                    chrx_male_vmiss.get(vid, {})),
            }

    # ---- chrY ----
    chry_dir = os.path.join(sex_chr_dir, 'chrY')
    if os.path.isdir(chry_dir):
        chry_male_vmiss = read_plink2_vmiss(
            os.path.join(chry_dir, 'chrY_male.vmiss'))
        chry_male_afreq = read_plink2_afreq(
            os.path.join(chry_dir, 'chrY_male.afreq'))

        all_chry_ids = sorted(set(
            list(chry_male_vmiss.keys()) +
            list(chry_male_afreq.keys())
        ))

        for vid in all_chry_ids:
            result['chrY'][vid] = {
                'chrY_male_call_rate': _fmiss_to_call_rate(
                    chry_male_vmiss.get(vid, {})),
                'chrY_male_maf': chry_male_afreq.get(vid, 'NA'),
            }

    return result


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
    parser.add_argument("--sex-chr-qc-dir", default=None,
                        help="Directory with sex-chromosome variant QC from "
                             "compute_sex_chr_variant_qc.sh (chrX/ and chrY/ "
                             "subdirectories)")
    parser.add_argument("--ancestry-sex-chr-qc", action='append', default=[],
                        help="Ancestry-specific sex-chr QC in format "
                             "ANCESTRY:DIR (repeatable)")
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

    # Load full-cohort sex-chromosome QC if available
    sex_chr_data = load_sex_chr_qc(args.sex_chr_qc_dir)
    has_chrx = bool(sex_chr_data.get('chrX'))
    has_chry = bool(sex_chr_data.get('chrY'))
    if has_chrx:
        print(f"  chrX variant QC: {len(sex_chr_data['chrX'])} variants")
    if has_chry:
        print(f"  chrY variant QC: {len(sex_chr_data['chrY'])} variants")

    # Load per-ancestry sex-chromosome QC
    ancestry_sex_chr_data = {}
    for spec in args.ancestry_sex_chr_qc:
        if ':' not in spec:
            print(f"Warning: Invalid ancestry-sex-chr-qc spec "
                  f"(expected ANCESTRY:DIR): {spec}", file=sys.stderr)
            continue
        anc, scqc_dir = spec.split(':', 1)
        if os.path.isdir(scqc_dir):
            anc_scqc = load_sex_chr_qc(scqc_dir)
            if anc_scqc.get('chrX') or anc_scqc.get('chrY'):
                ancestry_sex_chr_data[anc] = anc_scqc
                n_x = len(anc_scqc.get('chrX', {}))
                n_y = len(anc_scqc.get('chrY', {}))
                print(f"  {anc} chrX: {n_x}, chrY: {n_y} sex-chr variants")

    # Collect all variant IDs (autosomes + sex chromosomes)
    all_variants = set(all_data.keys())
    for anc_data in ancestry_data.values():
        all_variants.update(anc_data.keys())
    if has_chrx:
        all_variants.update(sex_chr_data['chrX'].keys())
    if has_chry:
        all_variants.update(sex_chr_data['chrY'].keys())
    for anc_scqc in ancestry_sex_chr_data.values():
        for chrom_data in anc_scqc.values():
            all_variants.update(chrom_data.keys())
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

    # Sex-chromosome columns (appended when data is available)
    chrx_columns = []
    if has_chrx:
        chrx_columns = [
            'chrX_call_rate', 'chrX_maf',
            'chrX_female_call_rate', 'chrX_female_hwe_p',
            'chrX_male_call_rate',
        ]
        columns.extend(chrx_columns)
    chry_columns = []
    if has_chry:
        chry_columns = [
            'chrY_male_call_rate', 'chrY_male_maf',
        ]
        columns.extend(chry_columns)

    # Per-ancestry sex-chromosome columns
    # These mirror the full-cohort sex-chr columns but are prefixed with
    # the ancestry label, e.g. EUR_chrX_female_hwe_p, AFR_chrY_male_maf.
    anc_sex_chr_sorted = sorted(ancestry_sex_chr_data.keys())
    anc_chrx_col_defs = {}  # ancestry -> list of column names
    anc_chry_col_defs = {}
    for anc in anc_sex_chr_sorted:
        anc_scqc = ancestry_sex_chr_data[anc]
        if anc_scqc.get('chrX'):
            cols = [
                f'{anc}_chrX_call_rate', f'{anc}_chrX_maf',
                f'{anc}_chrX_female_call_rate', f'{anc}_chrX_female_hwe_p',
                f'{anc}_chrX_male_call_rate',
            ]
            anc_chrx_col_defs[anc] = cols
            columns.extend(cols)
        if anc_scqc.get('chrY'):
            cols = [
                f'{anc}_chrY_male_call_rate', f'{anc}_chrY_male_maf',
            ]
            anc_chry_col_defs[anc] = cols
            columns.extend(cols)

    # Write collated output
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, 'w') as f:
        # Write project-level Ti/Tv as a header comment if available
        if tstv_ratio != 'NA':
            f.write(f'#tstv_ratio={tstv_ratio}\n')
        # Document sex-chromosome QC design decisions in header comments
        if has_chrx or has_chry or ancestry_sex_chr_data:
            f.write('#sex_chr_qc_note=chrX HWE computed on females only '
                    '(males are hemizygous; diploid HWE is undefined). '
                    'chrX_male_hwe_p is intentionally omitted. '
                    'chrY metrics are computed on males only. '
                    'PAR/XTR regions excluded. '
                    'Per-ancestry sex-chr columns (e.g. EUR_chrX_female_hwe_p) '
                    'follow the same conventions within each ancestry subset. '
                    'Ref: Laurie et al. 2012 Genet Epidemiol 36:384-91\n')
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

            # Append sex-chromosome columns (full cohort)
            if has_chrx:
                chrx_vd = sex_chr_data['chrX'].get(vid, {})
                for col in chrx_columns:
                    row.append(chrx_vd.get(col, 'NA'))
            if has_chry:
                chry_vd = sex_chr_data['chrY'].get(vid, {})
                for col in chry_columns:
                    row.append(chry_vd.get(col, 'NA'))

            # Append per-ancestry sex-chromosome columns
            for anc in anc_sex_chr_sorted:
                anc_scqc = ancestry_sex_chr_data[anc]
                if anc in anc_chrx_col_defs:
                    chrx_vd = anc_scqc.get('chrX', {}).get(vid, {})
                    # Map ANC_chrX_foo -> chrX_foo for internal lookup
                    for col in anc_chrx_col_defs[anc]:
                        internal_key = col[len(anc) + 1:]  # strip "ANC_"
                        row.append(chrx_vd.get(internal_key, 'NA'))
                if anc in anc_chry_col_defs:
                    chry_vd = anc_scqc.get('chrY', {}).get(vid, {})
                    for col in anc_chry_col_defs[anc]:
                        internal_key = col[len(anc) + 1:]
                        row.append(chry_vd.get(internal_key, 'NA'))

            f.write('\t'.join(row) + '\n')

    print(f"Collated variant QC: {args.output}")
    print(f"  Total variants: {len(all_variants)}")
    print(f"  Columns: {', '.join(columns)}")
    if tstv_ratio != 'NA':
        print(f"  Project Ti/Tv ratio: {tstv_ratio}")


if __name__ == "__main__":
    main()
