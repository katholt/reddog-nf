#!/usr/bin/env python3
import argparse
import pathlib
import sys


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sites_fps', required=True, type=pathlib.Path,
            help='SNP sites filepaths', nargs='+')
    parser.add_argument('--replicons_passing_fp', required=True, type=pathlib.Path,
            help='Isolate replicons passing filepath')
    args = parser.parse_args()
    for sites_fp in args.sites_fps:
        if not sites_fp.exists():
            parser.error(f'Input file {sites_fp} does not exist')
    if not args.replicons_passing_fp.exists():
        parser.error(f'Input file {args.replicons_passing_fp} does not exist')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Get a set of replicon + isolates that pass
    isolates = set()
    replicon_ids = set()
    with args.replicons_passing_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        for isolate_id, replicon_id in line_token_gen:
            replicon_ids.add(replicon_id)
            isolates.add((replicon_id, isolate_id))

    # During merge runs, extra files will be passed and need to be processed differently
    # These files have defined names, which we use to identify them
    merge_filenames = {f'merge_{rid}_sites.tsv' for rid in replicon_ids}

    # Read in SNP sites and combine
    snp_sites = set()
    for sites_fp in args.sites_fps:
        isolate = sites_fp.stem.replace('_sites', '')
        with sites_fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            next(line_token_gen)  # header
            for replicon_id, site in line_token_gen:
                if (replicon_id, isolate) in isolates:
                    snp_sites.add((replicon_id, int(site)))
                elif sites_fp.name in merge_filenames:
                    snp_sites.add((replicon_id, int(site)))

    # Print ordered sites, sort by position then replicon
    print('#CHROM', '#POS', sep='\t')
    for replicon_id, site in sorted(snp_sites):
        print(replicon_id, site, sep='\t')


if __name__ == '__main__':
    main()
