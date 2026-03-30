#!/usr/bin/env python3
"""
par_xtr_regions.py

Centralized PAR (pseudoautosomal region) and XTR (X-transposed region)
definitions for chrX and chrY per reference genome build.

These regions must be excluded from sex determination analyses (both
LRR-based and genotype-based F-statistic) because they behave as
autosomal loci:

  - PAR1/PAR2: homologous between chrX and chrY; males are diploid here.
  - XTR: X-transposed region with high X-Y sequence identity; reads may
    mismap between X and Y, inflating apparent heterozygosity in males.

Including PAR/XTR sites inflates male heterozygosity on chrX, biasing
the F-statistic toward 0 and causing false female or ambiguous calls.
On chrY, PAR/XTR sites bias the LRR median by including diploid-behaving
loci.

Supported genome builds:
  - CHM13 (T2T-CHM13v2.0)
  - GRCh38 (hg38)
  - GRCh37 (hg19)

To add a new genome build, add an entry to _PAR_XTR_REGIONS with the
build name as key and a list of (chrom, start, end, label) tuples.

References:
  - T2T CHM13v2.0: https://github.com/marbl/CHM13-annotations
  - T2T chrY PAR/XTR: chm13v2.0_PAR.bed, GIAB genome-stratifications v3.1
  - GRCh38: https://www.ncbi.nlm.nih.gov/grc/human
  - GRCh37: Ensembl GRCh37 assembly report
"""

import os
import sys
import tempfile

# Region definitions: {build: [(chrom, start, end, label), ...]}
# Coordinates are 0-based half-open (BED format).
# Both chrX and chrY PAR/XTR regions are included so that a single BED
# file can exclude diploid-behaving loci from both chromosomes.
_PAR_XTR_REGIONS = {
    'CHM13': [
        # chrX (T2T-CHM13v2.0)
        ('chrX', 0, 2781479, 'PAR1'),
        ('chrX', 2781479, 6400875, 'XTR'),
        ('chrX', 155701382, 156040895, 'PAR2'),
        # chrY (chm13v2.0_PAR.bed; GIAB genome-stratifications v3.1)
        ('chrY', 0, 2458320, 'PAR1'),
        ('chrY', 2458320, 6400875, 'XTR'),
        ('chrY', 62122809, 62460029, 'PAR2'),
    ],
    'GRCh38': [
        # chrX
        ('chrX', 10001, 2781479, 'PAR1'),
        ('chrX', 2781479, 6400000, 'XTR'),
        ('chrX', 155701383, 156030895, 'PAR2'),
        # chrY (NCBI GRCh38 assembly report)
        ('chrY', 10001, 2781479, 'PAR1'),
        # XTR on chrY is not well-defined for GRCh38 and is omitted.
        ('chrY', 56887903, 57217415, 'PAR2'),
    ],
    'GRCh37': [
        # chrX (note: GRCh37 uses non-chr-prefixed contig names)
        ('X', 60001, 2699520, 'PAR1'),
        ('X', 2699520, 6100000, 'XTR'),
        ('X', 154931044, 155260560, 'PAR2'),
        # chrY (Ensembl GRCh37; XTR on Y not clearly defined for this build)
        ('Y', 10001, 2649520, 'PAR1'),
        ('Y', 59034050, 59363566, 'PAR2'),
    ],
}

# Alias for common alternate names
_BUILD_ALIASES = {
    'hg38': 'GRCh38',
    'hg19': 'GRCh37',
    'T2T': 'CHM13',
    'chm13': 'CHM13',
    'chm13v2': 'CHM13',
    'chm13v2.0': 'CHM13',
}


def get_par_xtr_regions(genome='CHM13'):
    """Return list of (chrom, start, end, label) tuples for PAR/XTR regions.

    Args:
        genome: Genome build name (CHM13, GRCh38, GRCh37, or aliases).

    Returns:
        List of (chrom, start, end, label) tuples in BED coordinate space.

    Raises:
        ValueError: If genome build is not recognized.
    """
    build = _BUILD_ALIASES.get(genome, genome)
    if build not in _PAR_XTR_REGIONS:
        raise ValueError(
            f"Unknown genome build '{genome}'. "
            f"Supported: {', '.join(sorted(_PAR_XTR_REGIONS.keys()))}"
        )
    return list(_PAR_XTR_REGIONS[build])


def get_par_xtr_bed(genome='CHM13', output_path=None):
    """Write PAR/XTR regions to a BED file and return the path.

    If output_path is None, writes to a temporary file.

    Args:
        genome: Genome build name.
        output_path: Optional path to write the BED file.

    Returns:
        Path to the BED file.
    """
    regions = get_par_xtr_regions(genome)
    if output_path is None:
        fd, output_path = tempfile.mkstemp(suffix='.bed', prefix='par_xtr_')
        os.close(fd)
    with open(output_path, 'w') as f:
        for chrom, start, end, label in regions:
            f.write(f'{chrom}\t{start}\t{end}\t{label}\n')
    return output_path


def get_nonpar_chrx_regions(genome='CHM13'):
    """Return bcftools-compatible region strings for non-PAR/XTR chrX.

    Computes the complement of PAR/XTR on chrX only (chrY entries in the
    region table are ignored), returning a list of 'chrom:start-end'
    region strings (1-based inclusive, for bcftools -r).

    Args:
        genome: Genome build name.

    Returns:
        List of region strings suitable for bcftools -r.
    """
    all_regions = get_par_xtr_regions(genome)
    # Filter to chrX entries only
    regions = [r for r in all_regions if r[0] in ('chrX', 'X')]
    if not regions:
        return []

    chrom = regions[0][0]

    # Sort by start position
    sorted_regions = sorted(regions, key=lambda r: r[1])

    # Build non-PAR/XTR intervals
    nonpar = []
    prev_end = 0
    for _, start, end, _ in sorted_regions:
        if start > prev_end:
            # bcftools uses 1-based inclusive coordinates
            nonpar.append(f'{chrom}:{prev_end + 1}-{start}')
        prev_end = max(prev_end, end)

    # Add region after last PAR/XTR to end of chromosome
    # Use a generous sentinel: human chrX is ~156 Mb in CHM13 and ~155 Mb
    # in GRCh37/38; 200 Mb safely exceeds all current assemblies.
    chrx_end = 200000000
    if prev_end < chrx_end:
        nonpar.append(f'{chrom}:{prev_end + 1}-{chrx_end}')

    return nonpar


def supported_builds():
    """Return sorted list of supported genome build names."""
    return sorted(_PAR_XTR_REGIONS.keys())


def main():
    """CLI: print PAR/XTR regions in BED format for a given genome build."""
    import argparse
    parser = argparse.ArgumentParser(
        description='Print PAR/XTR region BED for a genome build'
    )
    parser.add_argument(
        '--genome', default='CHM13',
        help=f'Genome build ({", ".join(supported_builds())}). Default: CHM13'
    )
    parser.add_argument(
        '--output', default=None,
        help='Output BED file path (default: stdout)'
    )
    parser.add_argument(
        '--nonpar-regions', action='store_true',
        help='Print non-PAR/XTR chrX region strings (bcftools -r format)'
    )
    args = parser.parse_args()

    if args.nonpar_regions:
        for region in get_nonpar_chrx_regions(args.genome):
            print(region)
    elif args.output:
        path = get_par_xtr_bed(args.genome, args.output)
        print(f'Wrote {path}', file=sys.stderr)
    else:
        for chrom, start, end, label in get_par_xtr_regions(args.genome):
            print(f'{chrom}\t{start}\t{end}\t{label}')


if __name__ == '__main__':
    main()
