#!/usr/bin/env python3
import argparse
import pathlib
import statistics
import sys


class CheckInput(argparse.Action):

    def __call__(self, parser, namespace, filepaths, option_string=None):
        for filepath in filepaths:
            if not filepath.exists():
                parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, filepaths)


class CheckOutput(argparse.Action):

    def __call__(self, parser, namespace, filepath, option_string=None):
        if not filepath.exists():
            parser.error(f'Output {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, filepath)


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
            help='Replicon stats filepaths', nargs='+', action=CheckInput)
    parser.add_argument('--output_dir', required=True, type=pathlib.Path,
            help='Output directory', action=CheckOutput)
    return parser.parse_args()


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
        stats_pass = [stats for stats in stats_all if stats.pass_fail == 'p']
        ratio_stddev = statistics.stdev(stats.ratio for stats in stats_pass)
        ratio_mean = statistics.mean(stats.ratio for stats in stats_pass)
        ratio_max = ratio_mean + ratio_stddev * 2
        for stats in stats_pass:
            stats.phylogeny_group = 'i' if stats.ratio <= ratio_max else 'o'
        # Write to file, order by fail, outgroup, and then ingroup
        output_fp = args.output_dir / f'{replicon}_mapping_stats.tsv'
        with output_fp.open('w') as fh:
            stats_failed = [stats for stats in stats_all if stats.pass_fail == 'f']
            stats_outgroup = [stats for stats in stats_all if stats.phylogeny_group == 'o']
            stats_ingroup = [stats for stats in stats_all if stats.phylogeny_group == 'i']
            print(*Stats.fields_out, sep='\t', file=fh)
            print(*stats_failed, *stats_outgroup, *stats_ingroup, sep='\n', file=fh)


if __name__ == '__main__':
    main()
