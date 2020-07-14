#!/usr/bin/env python3
import argparse
import pathlib
import sys


import utility


class CheckInput(argparse.Action):

    def __call__(self, parser, namespace, filepath, option_string=None):
        if not filepath.exists():
            parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, filepath)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_fp', required=True, type=pathlib.Path,
            help='Filtered BAM filepath', action=CheckInput)
    parser.add_argument('--vcf_q30_fp', required=True, type=pathlib.Path,
            help='q30 VCF filepath', action=CheckInput)
    parser.add_argument('--vcf_hets_fp', required=True, type=pathlib.Path,
            help='Hets VCF filepath', action=CheckInput)
    parser.add_argument('--coverage_depth_fp', required=True, type=pathlib.Path,
            help='Coverage depth metrics filepath', action=CheckInput)
    parser.add_argument('--bam_unmapped_fp', required=True, type=pathlib.Path,
            help='Unmapped reads BAM filepath', action=CheckInput)

    parser.add_argument('--min_depth', type=float, default=10,
            help='Minimum depth to pass')
    parser.add_argument('--min_coverage', type=float, default=50,
            help='Minimum coverage to pass')
    parser.add_argument('--min_mapped_reads', type=float, default=50,
            help='Minimum reads mapped (applies only to largest replicon)')
    return parser.parse_args()


def main():
    # Get command line arguments
    args = get_arguments()

    # Get unmapped total reads
    total_unmapped = get_unmapped_read_count(args.bam_unmapped_fp)

    # Get coverage and depth, and mapped read counts
    print('Calculating coverage and average depth', file=sys.stderr)
    depth_average, coverage, replicon_sizes = get_depth_coverage_and_sizes(args.coverage_depth_fp)
    total_mapped, replicon_mapped = get_read_counts(args.bam_fp)

    # Calculate total mapped
    total_reads = total_mapped + total_unmapped

    # Get count of SNPs, INDELS, and heterozygous SNPs
    print('Collecting SNPs, INDELs, and heterozygous counts', file=sys.stderr)
    replicon_snps, replicon_indels = get_snp_indel_counts(args.vcf_q30_fp)
    replicon_hets = get_het_counts(args.vcf_hets_fp)

    # Get replicon names
    # This is an exhaustive collection of replicon names
    replicons_all = list(coverage) + list(replicon_snps) + list(replicon_hets)
    replicons = set(replicons_all)

    # Get largest replicon (used below for QC)
    replicon_largest = max(replicon_sizes, key=lambda k: replicon_sizes[k])

    # Calculate stats and output
    print('Calculating final statistics', file=sys.stderr)
    header = ('replicon', 'isolate', 'replicon_coverage', 'replicon_average_depth',
              'replicon_reads_mapped', 'total_mapped_reads', 'total_reads',
              'snps', 'snps_heterozygous', 'indels', 'pass_fail')
    print(*header, sep='\t')
    for replicon in sorted(replicons):
        # Check for missing data
        data_check = (
                coverage,
                depth_average,
                replicon_mapped,
                replicon_mapped,
                replicon_snps,
                replicon_hets,
                replicon_indels
                )
        for data in data_check:
            if replicon not in data:
                data[replicon] = 0
        # Aggregate
        isolate = args.bam_fp.stem
        stats = {
                'isolate': isolate,
                'replicon_coverage': round(coverage[replicon], 4),
                'replicon_average_depth': round(depth_average[replicon], 4),
                'replicon_reads_mapped': round(replicon_mapped[replicon] / total_reads * 100, 4),
                'total_mapped_reads': round(total_mapped / total_reads * 100, 4),
                'total_reads': total_reads,
                'snps': replicon_snps[replicon],
                'snps_heterogzyous': replicon_hets[replicon],
                'indels': replicon_indels[replicon],
                'pass_fail': None,
                }

        # Check mapping and SNP call quality for all
        # Percentage of reads mapped is only considered for the largest replicon
        if stats['replicon_average_depth'] < args.min_depth:
            stats['pass_fail'] = 'f'
        elif stats['replicon_coverage'] < args.min_coverage:
            stats['pass_fail'] = 'f'
        elif replicon == replicon_largest and stats['total_mapped_reads'] < args.min_mapped_reads:
            stats['pass_fail'] = 'f'
        else:
            stats['pass_fail'] = 'p'

        # Print
        print(replicon, *stats.values(), sep='\t')


def get_unmapped_read_count(bam_fp):
    # Execute command
    command = f'samtools view {bam_fp} | wc -l'
    result = utility.execute_command(command)
    return int(result.stdout.rstrip())


def get_depth_coverage_and_sizes(coverage_depth_fp):
    depth_averages = dict()
    coverages = dict()
    sizes = dict()
    token_types = (str, float, float, int)
    with coverage_depth_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        for line_tokens in line_token_gen:
            # Cast strings to appropriate types
            data = list()
            for token, token_type in zip(line_tokens, token_types):
                token_cast = token_type(token)
                data.append(token_cast)
            # Sort into variables
            replicon, depth_average, coverage, size = data
            depth_averages[replicon] = depth_average
            coverages[replicon] = coverage
            sizes[replicon] = size
    return depth_averages, coverages, sizes


def get_read_counts(bam_fp):
    # Execute command
    command = f'samtools view {bam_fp} | get_reads_mapped.awk'
    result = utility.execute_command(command)
    # Parse output
    line_token_gen = (line.rstrip().split('\t') for line in result.stdout.split('\n') if line)
    mapped_replicons = dict()
    token_types = (str, int)
    for line_tokens in line_token_gen:
        # Cast strings to appropriate types
        data = list()
        for token, token_type in zip(line_tokens, token_types):
            token_cast = token_type(token)
            data.append(token_cast)
        # Sort
        replicon, mapped = data
        mapped_replicons[replicon] = mapped
    mapped_total = mapped_replicons.pop('total_mapped')
    return mapped_total, mapped_replicons


def get_snp_indel_counts(vcf_q30_fp):
    # Execute command
    command = f'get_snp_indels.awk {vcf_q30_fp}'
    result = utility.execute_command(command)
    # Parse results
    line_token_gen = (line.rstrip().split('\t') for line in result.stdout.split('\n') if line)
    replicon_snps = dict()
    replicon_indels = dict()
    token_types = (str, int, int)
    for line_tokens in line_token_gen:
        # Cast strings to appropriate types
        data = list()
        for token, token_type in zip(line_tokens, token_types):
            token_cast = token_type(token)
            data.append(token_cast)
        # Sort
        replicon, snps, indels = data
        replicon_snps[replicon] = snps
        replicon_indels[replicon] = indels
    return replicon_snps, replicon_indels


def get_het_counts(vcf_hets_fp):
    # Execute command
    command = f'get_het_counts.awk {vcf_hets_fp}'
    result = utility.execute_command(command)
    line_token_gen = (line.rstrip().split('\t') for line in result.stdout.split('\n') if line)
    replicon_hets = {rep: int(hets) for rep, hets in line_token_gen}
    return replicon_hets


if __name__ == '__main__':
    main()
