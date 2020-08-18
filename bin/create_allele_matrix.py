#!/usr/bin/env python3
import argparse
import pathlib
import re
import subprocess
import sys


import utility


dp4_regex = re.compile(r'DP4=([0-9,]+)')


vcf_fields = {'chrom': str,
              'pos': int,
              'id': str,
              'ref': str,
              'alt': str,
              'qual': float,
              'filter': str,
              'info': str,
              'format': str,
              'genotype': str}


class VcfRecord:

    def __init__(self, data):
        for value, (attr, attr_type) in zip(data, vcf_fields.items()):
            setattr(self, attr, attr_type(value))
        self._data = data

    def __str__(self):
        return '\t'.join(self._data)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_fp', required=True, type=pathlib.Path,
            help='Sample BAM filepath')
    parser.add_argument('--sites_fp', required=True, type=pathlib.Path,
            help='Replicon sites filepaths')
    parser.add_argument('--reference_fp', required=True, type=pathlib.Path,
            help='Output directory')
    parser.add_argument('--output_dir', required=True, type=pathlib.Path,
            help='Output directory')

    parser.add_argument('--min_quality', default=20, type=int,
            help='Minimum phread quality to call allele')
    parser.add_argument('--min_support', default=90, type=float,
            help='Minimum proportion of reads support reference allele call')

    args = parser.parse_args()
    if not args.bam_fp.exists():
        parser.error(f'Input file {args.bam_fp} does not exist')
    if not args.reference_fp.exists():
        parser.error(f'Input file {args.reference_fp} does not exist')
    if not args.output_dir.exists():
        parser.error(f'Output directory {args.output_dir} does not exist')
    if args.min_quality <= 0:
        parser.error('--min_quality must be greater than 0, got {args.min_quality}')
    if args.min_support > 100:
        parser.error('--min_support must be less than 100, got {args.min_support}')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Get isolate id
    isolate_id = args.bam_fp.stem.replace('_filtered', '')

    # Get site-specific consensus
    command_mpileup_base = 'bcftools mpileup -q20 -Q20 -B'
    command_mpileup_inputs = f'-f {args.reference_fp} -R {args.sites_fp} {args.bam_fp}'
    command_call = 'bcftools call -V indels -c --ploidy 1 -Ov'
    command = f'{command_mpileup_base} {command_mpileup_inputs} | {command_call}'
    result = utility.execute_command(command)

    # Check we have data to process
    lines = result.stdout.rstrip().split('\n')
    for i, line in enumerate(lines):
        if not line.startswith('#'):
            break

    # Iterate VCF records and print alleles to appropriate replicon file
    line_token_gen = (line.rstrip().split('\t') for line in lines[i:])
    record_gen = (VcfRecord(line_tokens) for line_tokens in line_token_gen)
    replicon_current = None
    fh = None
    for record in record_gen:
        # Skip records with more than one allele
        if ',' in record.alt:
            continue
        # Check if we need to change output file
        if replicon_current != record.chrom:
            if replicon_current:
                fh.close()
            replicon_current = record.chrom
            output_fp = args.output_dir / f'{replicon_current}_{isolate_id}_alleles.csv'
            fh = output_fp.open('w')
            print('Pos', 'Reference', isolate_id, sep=',', file=fh)
        # Check quality and support
        support_ref, support_alt = get_position_support(record.info)
        if record.qual < args.min_quality:
            print(record.pos, record.ref, '-', sep=',', file=fh)
        elif record.alt != '.' and support_alt >= args.min_support:
            print(record.pos, record.ref, record.alt, sep=',', file=fh)
        elif record.alt == '.' and support_ref >= args.min_support:
            print(record.pos, record.ref, record.ref, sep=',', file=fh)
        else:
            print(record.pos, record.ref, '-', sep=',', file=fh)


def get_position_support(record_info):
    # Collect tokens, convert, and sum
    dp4_string = dp4_regex.search(record_info).group(1)
    dp4_counts = [int(token) for token in dp4_string.split(',')]
    reads_ref = sum(dp4_counts[:2])
    reads_alt = sum(dp4_counts[2:])
    # Get support proportions
    reads_sum = reads_ref + reads_alt
    support_ref = reads_ref / reads_sum * 100 if reads_sum else 0
    support_alt = reads_alt / reads_sum * 100 if reads_sum else 0
    return (support_ref, support_alt)


if __name__ == '__main__':
    main()
