#!/usr/bin/env python3
import argparse
import hashlib
import pathlib
import sys


import Bio.SeqIO


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference_fp', required=True, type=pathlib.Path,
            help='Reference filepath')
    args = parser.parse_args()
    if not args.reference_fp.exists():
        parser.error(f'Input file {args.reference_fp} does not exist')
    return args


def main():
    # Get commandline arguments
    args = get_arguments()

    # Check we have genbank by file extension
    reference_extension = args.reference_fp.suffix[1:]
    if reference_extension == 'fasta':
        print('error: reference provided as fasta format but must be in genbank format', file=sys.stderr)
        sys.exit(1)
    elif reference_extension not in {'gbk', 'gb', 'gbff', 'genbank'}:
        print('error: reference does not appear to be in genbank format', file=sys.stderr)
        sys.exit(1)

    # Parse reference
    with args.reference_fp.open('r') as fh:
        try:
            records = [record for record in Bio.SeqIO.parse(fh, 'genbank')]
        except:
            print('error: reference is not in genbank format', file=sys.stderr)
            sys.exit(1)

    # Check we have records
    if len(records) < 1:
        print('error: could not find any records in provided reference', file=sys.stderr)
        sys.exit(1)

    # Get replicon names, size, and m5dsum hash of sequence
    names = list()
    sizes = list()
    hashes = list()
    for record in records:
        seq_bin = record.seq.encode()
        names.append(record.name)
        sizes.append(f'{record.name}:{len(record.seq)}')
        hashes.append(f'{record.name}:{hashlib.md5(seq_bin).hexdigest()}')

    # Print data
    print('reference_name', args.reference_fp.stem, sep='\t')
    print('reference_replicon_names', ','.join(names), sep='\t')
    print('reference_replicon_sizes', ','.join(sizes), sep='\t')
    print('reference_replicon_md5hashes', ','.join(hashes), sep='\t')


if __name__ == '__main__':
    main()
