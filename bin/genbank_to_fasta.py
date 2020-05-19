#!/usr/bin/env python3
import argparse
import contextlib
import pathlib


import Bio.SeqIO


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fp', required=True, type=pathlib.Path,
            help='Input genbank file')
    parser.add_argument('--output_fp', required=True, type=pathlib.Path,
            help='Output file fasta file')
    args = parser.parse_args()
    if not args.input_fp.exists():
        parser.error('Input file {args.input_fp} does not exist')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Convert
    with contextlib.ExitStack() as stack:
        input_fh = stack.enter_context(args.input_fp.open('r'))
        output_fh = stack.enter_context(args.output_fp.open('w'))
        genbank = Bio.SeqIO.parse(input_fh, "genbank")
        Bio.SeqIO.write(genbank, output_fh, "fasta")


if __name__ == '__main__':
    main()
