#!/usr/bin/env python3
import argparse
import pathlib
import re
import subprocess
import sys


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

    def __call__(self, parser, namespace, value, option_string=None):
        value_check = [value] if isinstance(value, pathlib.Path) else value
        for filepath in value_check:
            if not filepath.exists():
                parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, value)


class CheckOutput(argparse.Action):

    def __call__(self, parser, namespace, filepath, option_string=None):
        if not filepath.exists():
            parser.error(f'Output directory {filepath.parent} does not exist')
        setattr(namespace, self.dest, filepath)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_fps', required=True, type=pathlib.Path,
            help='Input VCF filepaths', nargs='+', action=CheckInput)
    parser.add_argument('--replicon', required=True, type=str,
            help='Replicon name')
    return parser.parse_args()


def main():
    # Get command line arguments
    args = get_arguments()

    # For each VCF get SNPs that pass
    replicons_seen = set()
    ref_alleles = dict()
    set_alleles = list()
    for vcf_fp in args.vcf_fps:
        vcf_alleles = dict()
        with vcf_fp.open('r') as fh:
            # Skip header
            while fh.readline().startswith('#'):
                position = fh.tell()

            # Check if there is any data to process
            if fh.readline():
                fh.seek(position, 0)
            else:
                break

            # Collect SNPs
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            record_gen = (VcfRecord(line_tokens) for line_tokens in line_token_gen)
            for record in record_gen:
                # Skip INDELs
                if record.info.startswith('I'):
                    continue
                # Skip non-SNPs
                if ',' in record.alt:
                    continue
                # Record alt allele
                if record.qual >= 20 and record.alt in {'A', 'T', 'G', 'C'}:
                    vcf_alleles[record.pos] = record.alt
                else:
                    vcf_alleles[record.pos] = '-'
                # Record ref allele
                if record.pos not in ref_alleles:
                    ref_alleles[record.pos] = record.ref
            set_alleles.append(vcf_alleles)

            # Ensure that we're getting the same replicon
            replicons_seen.add(record.chrom)
            if len(replicons_seen) != 1:
                print('error: got duplicate replicon names', file=sys.stderr)
                sys.exit(1)

    # Get SNP positions
    snp_positions = {pos for vcf_alleles in set_alleles for pos in vcf_alleles}
    snp_positions = sorted(snp_positions)

    # Compare seen replicon to provided
    if replicons_seen:
        [replicon] = replicons_seen
        if args.replicon != replicon:
            msg = f'error: got replicon id mismatch, expected {args.replicon} but got {replicon}'
            print(msg, file=sys.stderr)
            sys.exit(1)

    # Header
    # Create matrix
    sample_ids = [fp.stem.replace(f'_{args.replicon}_q30', '') for fp in args.vcf_fps]
    print('Position', 'Reference', *sample_ids, sep='\t')
    # Data
    for snp_position in snp_positions:
        reference_allele = ref_alleles[snp_position]
        position_alleles = [alleles.get(snp_position, reference_allele) for alleles in set_alleles]
        print(snp_position, reference_allele, *position_alleles, sep='\t')


if __name__ == '__main__':
    main()
