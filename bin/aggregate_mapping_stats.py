#!/usr/bin/env python3
import argparse
import pathlib
import statistics
import sys


class Stats:

    fields = {'replicon': str,
              'isolate': str,
              'replicon_coverage': float,
              'replicon_average_depth': float,
              'replicon_reads_mapped': float,
              'total_mapped_reads': float,
              'total_reads': int,
              'snps': int,
              'snps_heterozygous': int,
              'indels': int,
              'pass_fail': str
             }

    fields_out = [field for field in fields if field != 'replicon'] + ['phylogeny_group']


    def __init__(self, values):
        for value, (attr, attr_type) in zip(values, self.fields.items()):
            setattr(self, attr, attr_type(value))
        self._ratio = float()
        self.phylogeny_group = 'n/a'


    @property
    def ratio(self):
        if not self._ratio:
            self._ratio = self.total_reads / self.replicon_coverage / 100
        return self._ratio


    def __str__(self):
        token_gen = (getattr(self, attr) for attr in self.fields_out)
        return '\t'.join(str(value) for value in token_gen)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rep_stats_fps', required=True, type=pathlib.Path,
            help='Replicon stats filepaths', nargs='+')
    parser.add_argument('--output_dir', required=True, type=pathlib.Path,
            help='Output directory')
    parser.add_argument('--stddev_mod', default=2, type=float,
            help='Modifier for outgroup designation')
    args = parser.parse_args()
    for rep_stats_fp in args.rep_stats_fps:
        if not rep_stats_fp.exists():
            parser.error(f'Input file {rep_stats_fp} does not exist')
    if not args.output_dir.exists():
        parser.error(f'Output directory {args.output_dir} does not exist')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read all data in and aggregate by replicon name
    replicon_stats = dict()
    for rep_stats_fp in args.rep_stats_fps:
        with rep_stats_fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            next(line_token_gen)  # header
            for line_tokens in line_token_gen:
                stats = Stats(line_tokens)
                if stats.replicon not in replicon_stats:
                    replicon_stats[stats.replicon] = list()
                replicon_stats[stats.replicon].append(stats)

    # Calculate ingroup/outgroup for each
    for replicon, stats_all in replicon_stats.items():
        # Get all passing entries and set to undetermined
        stats_pass = [stats for stats in stats_all if stats.pass_fail == 'p']
        for stats in stats_pass:
            stats.phylogeny_group = 'undetermined'
        # If we have sufficient entries, calculate ingroup/outgroups
        if len(stats_pass) > 1:
            ratio_stddev = statistics.stdev(stats.ratio for stats in stats_pass)
            ratio_mean = statistics.mean(stats.ratio for stats in stats_pass)
            ratio_max = ratio_mean + ratio_stddev * args.stddev_mod
            for stats in stats_pass:
                stats.phylogeny_group = 'i' if stats.ratio <= ratio_max else 'o'
        # Write to file, order by fail, outgroup, and then ingroup
        output_fp = args.output_dir / f'{replicon}_mapping_stats.tsv'
        with output_fp.open('w') as fh:
            # Explicit sorting to avoiding missing entries
            stats_failed = list()
            stats_undetermined = list()
            stats_outgroup = list()
            stats_ingroup = list()
            for stats in stats_all:
                if stats.pass_fail == 'f':
                    stats_failed.append(stats)
                elif stats.phylogeny_group == 'undetermined':
                    stats_undetermined.append(stats)
                elif stats.phylogeny_group == 'o':
                    stats_outgroup.append(stats)
                elif stats.phylogeny_group == 'i':
                    stats_ingroup.append(stats)
                else:
                    raise ValueError(f'could not sort entry into an output group:\n{stats}')
            print(*Stats.fields_out, sep='\t', file=fh)
            print(*stats_failed, *stats_undetermined, *stats_outgroup, *stats_ingroup, sep='\n', file=fh)


if __name__ == '__main__':
    main()
