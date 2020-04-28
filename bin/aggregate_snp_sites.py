#!/usr/bin/env python3
import argparse
import pathlib
import sys


class CheckInput(argparse.Action):

    def __call__(self, parser, namespace, filepaths, option_string=None):
        if not isinstance(filepaths, list):
            check_value = [filepaths]
        else:
            check_value = filepaths
        for filepath in check_value:
            if not filepath.exists():
                parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, filepaths)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sites_fps', required=True, type=pathlib.Path,
            help='SNP sites filepaths', nargs='+', action=CheckInput)
    parser.add_argument('--isolate_replicons_passing_fp', required=True, type=pathlib.Path,
            help='Isolate replicons passing filepath', action=CheckInput)
    return parser.parse_args()


def main():
    # Get command line arguments
    args = get_arguments()

    # Get a set of replicon + isolates that pass
    isolates = set()
    with args.isolate_replicons_passing_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        for isolate, *replicons in line_token_gen:
            for replicon in replicons:
                isolates.add((replicon, isolate))

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
