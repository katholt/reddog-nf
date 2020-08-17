#!/usr/bin/env python3
import argparse
import pathlib
import sys


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--allele_fp', required=True, type=pathlib.Path,
            help='Input allele matrix filepath')
    args = parser.parse_args()
    if not args.allele_fp.exists():
        parser.error(f'Input file {args.allele_fp} does not exist')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read allele matrix into memory
    with args.allele_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split(',') for line in fh)
        header_tokens = next(line_token_gen)
        isolates = ['Reference', *header_tokens[2:]]
        column_gen = zip(*line_token_gen)
        next(column_gen)  # positions
        for isolate, snps in zip(isolates, column_gen):
            desc = f'>{isolate}'
            snp_alignment = ''.join(snps)
            seq_lines = [snp_alignment[i:i+80] for i in range(0, len(snp_alignment), 80)]
            print(desc)
            print(*seq_lines, sep='\n')


if __name__ == '__main__':
    main()
