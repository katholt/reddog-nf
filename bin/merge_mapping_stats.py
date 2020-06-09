#!/usr/bin/env python3
import argparse
import pathlib
import statistics


class Record:

    fields = None

    def __init__(self, data):
        assert self.fields
        for attr, value in zip(self.fields, data):
            setattr(self, attr, value)
        self.ratio = float()

    def __str__(self):
        return '\t'.join(getattr(self, attr) for attr in self.fields)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fp_1', required=True, type=pathlib.Path,
            help='First filepath')
    parser.add_argument('--fp_2', required=True, type=pathlib.Path,
            help='Second filepath')
    parser.add_argument('--stddev_mod', default=2, type=float,
            help='Modifier for outgroup designation')
    args = parser.parse_args()
    if not args.fp_1.exists():
        parser.error(f'Input file {args.fp_1} does not exist')
    if not args.fp_1.exists():
        parser.error(f'Input file {args.fp_2} does not exist')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read data into memory and recalculate ingroup/outgroup
    records = read_data(args.fp_1, args.fp_2)
    # Ratios for all passing isolates
    ratios = list()
    for record in records:
        if record.pass_fail != 'p':
            continue
        record.ratio = int(record.total_reads) / float(record.replicon_coverage) / 100
        ratios.append(record.ratio)
    # If we have sufficient entries, calculate ingroup/outgroup threshold
    if len(ratios) > 1:
        ratio_stddev = statistics.stdev(ratios)
        ratio_mean = statistics.mean(ratios)
        ratio_max = ratio_mean + ratio_stddev * args.stddev_mod
    # Apply groups
    for record in records:
        if not record.ratio:
            continue
        if len(ratios) <= 1:
            record.phylogeny_group = 'undetermined'
        elif record.ratio <= ratio_max:
            record.phylogeny_group = 'i'
        elif record.ratio > ratio_max:
            record.phylogeny_group = 'o'

    # Sort and print data
    print(*Record.fields, sep='\t')
    # Explicit sorting to avoiding missing entries
    records_failed = list()
    records_undetermined = list()
    records_outgroup = list()
    records_ingroup = list()
    for record in records:
        if record.pass_fail == 'f':
            records_failed.append(record)
        elif record.phylogeny_group == 'undetermined':
            records_undetermined.append(record)
        elif record.phylogeny_group == 'o':
            records_outgroup.append(record)
        elif record.phylogeny_group == 'i':
            records_ingroup.append(record)
        else:
            raise ValueError(f'could not sort entry into an output group:\n{stats}')
    print(*records_failed, *records_undetermined, *records_outgroup, *records_ingroup, sep='\n')


def read_data(fp_1, fp_2):
    # First file
    records = list()
    with fp_1.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        Record.fields = next(line_token_gen)
        for line_tokens in line_token_gen:
            record = Record(line_tokens)
            records.append(record)
    # Second file
    with fp_2.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        assert Record.fields == next(line_token_gen)
        for line_tokens in line_token_gen:
            record = Record(line_tokens)
            records.append(record)
    return records


if __name__ == '__main__':
    main()
