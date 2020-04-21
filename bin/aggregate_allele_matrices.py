#!/usr/bin/env python3
import argparse
import pathlib
import sys


class CheckInput(argparse.Action):

    def __call__(self, parser, namespace, value, option_string=None):
        value_check = value if isinstance(value, list) else [value]
        for filepath in value_check:
            if not filepath.exists():
                parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, value)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--allele_fps', required=True, type=pathlib.Path,
            help='Allele matrix filepaths', nargs='+', action=CheckInput)
    parser.add_argument('--sites_fp', required=True, type=pathlib.Path,
            help='SNPs sites filepath', action=CheckInput)
    return parser.parse_args()


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in SNP sites
    alleles = dict()
    with args.sites_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        next(line_token_gen)  # header
        for replicon, position_str in line_token_gen:
            position = int(position_str)
            alleles[position] = list()

    # Set to compare against and discover missing sites
    sites_all = set(alleles)

    # Record reference alleles separately
    ref_alleles = dict()

    # Serially read in individual allele matrices
    isolates = list()
    for allele_fp in args.allele_fps:
        sites_seen = set()
        # Record site alleles
        with allele_fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            isolate = next(line_token_gen)[2]
            isolates.append(isolate)
            for position_str, ref_snp, isolate_snp in line_token_gen:
                position = int(position_str)
                sites_seen.add(position)
                alleles[position].append(isolate_snp)
                if position not in ref_alleles:
                    ref_alleles[position] = ref_snp
        # Complete missing sites
        for position in sites_all ^ sites_seen:
            alleles[position].append('-')

    # Print aggregated matrix
    print('Position', 'Reference', *isolates, sep='\t')
    for position in sorted(sites_all):
        print(position, ref_alleles[position], *alleles[position], sep='\t')


if __name__ == '__main__':
    main()
