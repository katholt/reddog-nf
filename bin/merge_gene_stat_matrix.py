#!/usr/bin/env python3
import argparse
import pathlib


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fp_1', required=True, type=pathlib.Path,
            help='First filepath')
    parser.add_argument('--fp_2', required=True, type=pathlib.Path,
            help='Second filepath')
    args = parser.parse_args()
    if not args.fp_1.exists():
        parser.error(f'Input file {args.fp_1} does not exist')
    if not args.fp_1.exists():
        parser.error(f'Input file {args.fp_2} does not exist')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read everything into memory and process
    isolates_1, genes_1 = read_data(args.fp_1)
    isolates_2, genes_2 = read_data(args.fp_2)

    # Iterate all genes and complete data if missing
    print('replicon_id', 'locus_tag', *isolates_1, *isolates_2, sep=',')
    for replicon_id, gene in set(genes_1) | set(genes_2):
        print(replicon_id, gene, sep=',', end=',')
        if (replicon_id, gene) not in genes_1:
            print(*['0.0'] * len(isolates_1), sep=',', end=',')
        else:
            print(genes_1[(replicon_id, gene)], sep=',', end=',')
        if (replicon_id, gene) not in genes_2:
            print(*['0.0'] * len(isolates_2), sep=',')
        else:
            print(genes_2[(replicon_id, gene)], sep=',')


def read_data(filepath):
    genes = dict()
    with filepath.open('r') as fh:
        header_tokens = fh.readline().rstrip().split(',')
        for line in fh:
            replicon_id, gene, line_str = line.rstrip().split(',', maxsplit=2)
            assert (replicon_id, gene) not in genes
            genes[(replicon_id, gene)] = line_str
    return header_tokens[2:], genes


if __name__ == '__main__':
    main()
