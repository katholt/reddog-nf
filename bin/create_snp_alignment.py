#!/usr/bin/env python3
import argparse
import pathlib
import sys


class CheckInput(argparse.Action):

    def __call__(self, parser, namespace, filepath, option_string=None):
        if not filepath.exists():
            parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, filepath)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--allele_fp', required=True, type=pathlib.Path,
            help='Input allele matrix filepath', action=CheckInput)
    return parser.parse_args()


def main():
    # Get command line arguments
    args = get_arguments()

    # Read allele matrix into memory
    with args.allele_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
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
