#!/usr/bin/env python3
import argparse
import contextlib
import pathlib


import Bio.SeqIO


class CheckInput(argparse.Action):

    def __call__(self, parser, namespace, filepath, option_string=None):
        if not filepath.exists():
            parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, filepath)


class CheckOutput(argparse.Action):

    def __call__(self, parser, namespace, filepath, option_string=None):
        if not filepath.parent.exists():
            parser.error(f'Directory {filepath.parent} for {option_string} does not exist')
        setattr(namespace, self.dest, filepath)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fp', required=True, type=pathlib.Path,
            help='Input genbank file', action=CheckInput)
    parser.add_argument('--output_fp', required=True, type=pathlib.Path,
            help='Ouput file fasta file', action=CheckOutput)
    return parser.parse_args()


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
