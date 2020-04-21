#!/usr/bin/env python3
# TODO: revert base quality (-Q) to 20, reduced for testing purposes


import argparse
import pathlib
import re
import subprocess
import sys


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


class CheckInput(argparse.Action):

    def __call__(self, parser, namespace, filepath, option_string=None):
        if not filepath.exists():
            parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, filepath)


class CheckOutput(argparse.Action):

    def __call__(self, parser, namespace, filepath, option_string=None):
        if not filepath.exists():
            parser.error(f'Output directory {filepath.parent} does not exist')
        setattr(namespace, self.dest, filepath)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_fp', required=True, type=pathlib.Path,
            help='Sample BAM filepath', action=CheckInput)
    parser.add_argument('--sites_fp', required=True, type=pathlib.Path,
            help='Replicon sites filepaths', action=CheckInput)
    parser.add_argument('--reference_fp', required=True, type=pathlib.Path,
            help='Output directory', action=CheckOutput)
    parser.add_argument('--output_dir', required=True, type=pathlib.Path,
            help='Output directory', action=CheckOutput)

    parser.add_argument('--min_quality', default=20, type=int,
            help='Minimum phread quality to call allele')
    parser.add_argument('--min_support', default=0.9, type=float,
            help='Minimum proportion of reads support reference allele call')
    return parser.parse_args()


def main():
    # Get command line arguments
    args = get_arguments()

    # Get isolate id
    isolate_id = args.bam_fp.stem.replace('_filtered', '')

    # Get site-specific consensus
    command_mpileup_base = 'bcftools mpileup -q20 -Q5 -B'
    command_mpileup_inputs = f'-f {args.reference_fp} -R {args.sites_fp} {args.bam_fp}'
    command_call = 'bcftools call -V indels -c --ploidy 1 -Ov'
    command = f'{command_mpileup_base} {command_mpileup_inputs} | {command_call}'
    result = execute_command(command)

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
        # Check if we need to change output file
        if replicon_current != record.chrom:
            if replicon_current:
                fh.close()
            replicon_current = record.chrom
            output_fp = args.output_dir / f'{replicon_current}_{isolate_id}_alleles.tsv'
            fh = output_fp.open('w')
            print('Position', 'Reference', isolate_id, sep='\t', file=fh)
        # Check quality and support
        support_ref, support_alt = get_position_support(record.info)
        if record.qual < args.min_quality:
            print(record.pos, record.ref, '-', sep='\t', file=fh)
        elif record.alt != '.' and support_alt >= args.min_support:
            print(record.pos, record.ref, record.alt, sep='\t', file=fh)
        elif record.alt == '.' and support_ref >= args.min_support:
            print(record.pos, record.ref, record.ref, sep='\t', file=fh)
        else:
            print(record.pos, record.ref, '-', sep='\t', file=fh)


def get_position_support(record_info):
    # Collect tokens, convert, and sum
    dp4_string = dp4_regex.search(record_info).group(1)
    dp4_counts = [int(token) for token in dp4_string.split(',')]
    reads_ref = sum(dp4_counts[:2])
    reads_alt = sum(dp4_counts[2:])
    # Get support proportions
    reads_sum = reads_ref + reads_alt
    support_ref = reads_ref / reads_sum if reads_sum else 0
    support_alt = reads_alt / reads_sum if reads_sum else 0
    return (support_ref, support_alt)


def execute_command(command):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8')
    if result.returncode != 0:
        print('Failed to run command:', result.args, file=sys.stderr)
        print('stdout:', result.stdout, file=sys.stderr)
        print('stderr:', result.stderr, file=sys.stderr)
        sys.exit(1)
    return result


if __name__ == '__main__':
    main()
