#!/usr/bin/env python3
import argparse
import pathlib


import Bio.SeqIO


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fp', required=True, type=pathlib.Path,
            help='Input filepath')
    args = parser.parse_args()
    if not args.input_fp.exists():
        parser.error(f'Input file {args.input_fp} does not exist')
    return args


def main():
    # Get commandline arguments
    args = get_arguments()

    # Get list of coding genes
    features = list()
    with args.input_fp.open('r') as fh:
        records = {record.name: record for record in Bio.SeqIO.parse(fh, 'genbank')}
        for record in records.values():
            for feature in record.features:
                if feature.type != 'CDS':
                    continue
                feature.record_name = record.name
                features.append(feature)

    # Print top 10 genes as fasta
    i = 0
    features = sorted(features, key=lambda k: k.location.end - k.location.start, reverse=True)
    for feature in features:
        if 'gene' not in feature.qualifiers:
            continue
        if feature.location.end - feature.location.start + 1 > 2500:
            continue
        i += 1
        if i > 10:
            break
        seqfeature = feature.extract(records[feature.record_name])
        sequence = seqfeature.seq

        [gene] = feature.qualifiers['gene']
        seqlines = [sequence[i:i+80] for i in range(0, len(sequence), 80)]
        print(f'>{gene}')
        print(*seqlines, sep='\n')


if __name__ == '__main__':
    main()
