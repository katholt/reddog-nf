#!/usr/bin/env python3
import argparse
import pathlib


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--allele_fp', required=True, type=pathlib.Path,
            help='Input allele matrix filepath')
    parser.add_argument('--conservation', default=95, type=int,
            help='Minimum isolates required to have a known allele at a given site (0-100)')
    args = parser.parse_args()
    if not args.allele_fp.exists():
        parser.error(f'Input file {args.allele_fp} does not exist')
    if args.conservation <= 0:
        parser.error('--conservation must be greater than 0, got {args.conversation}')
    if args.conservation > 100:
        parser.error('--conservation must be less than 100, got {args.conversation}')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Filter at n% conservation
    with args.allele_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split(',') for line in fh)
        # Get and print header
        header_tokens = next(line_token_gen)
        print(*header_tokens, sep=',')
        # Process alleles
        isolate_count = len(header_tokens) - 2
        for position, allele_ref, *allele_isolates in line_token_gen:
            # Remove invariant SNPs
            alleles_unique = {allele for allele in allele_isolates if allele != '-'}
            alleles_unique.add(allele_ref)
            if len(alleles_unique) == 1:
                continue
            # Remove SNPs with fewer than n% known alleles
            alleles_known_pct = (1 - allele_isolates.count('-') / isolate_count) * 100
            if alleles_known_pct < args.conservation:
                continue
            # Print passing alleles
            print(position, allele_ref, *allele_isolates, sep=',')


if __name__ == '__main__':
    main()
