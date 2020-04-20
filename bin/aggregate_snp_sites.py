#!/usr/bin/env python3
import argparse
import pathlib
import sys


class CheckInput(argparse.Action):

    def __call__(self, parser, namespace, filepaths, option_string=None):
        for filepath in filepaths:
            if not filepath.exists():
                parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, filepaths)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sites_fps', required=True, type=pathlib.Path,
            help='SNP sites filepaths', nargs='+', action=CheckInput)
    parser.add_argument('--replicon_stats_fps', required=True, type=pathlib.Path,
            help='Replicon statistics filepaths', nargs='+', action=CheckInput)
    return parser.parse_args()


def main():
    # Get command line arguments
    args = get_arguments()

    # Get a set of replicon + isolates that pass
    isolates = set()
    for replicon_stats_fp in args.replicon_stats_fps:
        replicon = replicon_stats_fp.stem.replace('_RepStats', '')
        with replicon_stats_fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            header_tokens = next(line_token_gen)
            for line_tokens in line_token_gen:
                record = {field: value for field, value in zip(header_tokens, line_tokens)}
                if record['pass_fail'] == 'f':
                    continue
                isolates.add((replicon, record['isolate']))

    # Read in SNP sites and combine
    snp_sites = set()
    for sites_fp in args.sites_fps:
        isolate = sites_fp.stem.replace('_sites', '')
        with sites_fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            next(line_token_gen)  # header
            for replicon, site in line_token_gen:
                if (replicon, isolate) not in isolates:
                    continue
                snp_sites.add((replicon, int(site)))

    # Print ordered sites, sort by position then replicon
    print('#CHROM', '#POS', sep='\t')
    for replicon, site in sorted(snp_sites):
        print(replicon, site, sep='\t')


if __name__ == '__main__':
    main()
