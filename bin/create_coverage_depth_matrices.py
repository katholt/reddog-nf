#!/usr/bin/env python3
import argparse
import contextlib
import pathlib
import sys


import Bio.SeqIO


import utility


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mpileup_fps', required=True, type=pathlib.Path,
            help='samtools mpileup filepaths', nargs='+')
    parser.add_argument('--reference_fp', required=True, type=pathlib.Path,
            help='Reference filepath (genbank format)')
    parser.add_argument('--output_dir', required=True, type=pathlib.Path,
            help='Output directory')
    args = parser.parse_args()
    for mpileup_fp in args.mpileup_fps:
        if not mpileup_fp.exists():
            parser.error(f'Input file {mpileup_fp} does not exist')
    if not args.reference_fp.exists():
        parser.error(f'Input file {args.reference_fp} does not exist')
    if not args.output_dir.exists():
        parser.error(f'Output directory {args.output_dir} does not exist')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in reference and get list of coding features for each replicon
    ref_gbk = {record.id: record for record in Bio.SeqIO.parse(args.reference_fp, 'genbank')}
    replicon_genes = dict()
    for record in Bio.SeqIO.parse(args.reference_fp, 'genbank'):
        replicon_genes[record.id] = list()
        for feature in record.features:
            if feature.type != 'CDS':
                continue
            replicon_genes[record.id].append(feature)

    # Process each isolate mpileup
    gene_stats = dict()
    for mpileup_fp in args.mpileup_fps:
        # Read entire mpileup into memory
        with mpileup_fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            replicon_depths = dict()
            for rep, pos, base, depth_str, reads, qual in line_token_gen:
                if rep not in replicon_depths:
                    replicon_depths[rep] = list()
                depth = int(depth_str)
                replicon_depths[rep].append(depth)
        # For each feature, calculate the mean depth and coverage where depth is greater than n
        replicon_gene_stats = {replicon: dict() for replicon in replicon_depths}
        for replicon in replicon_depths:
            for feature in replicon_genes[replicon]:
                start = feature.location.nofuzzy_start
                end = feature.location.nofuzzy_end
                feature_depths = replicon_depths[replicon][start:end]
                depth_mean = sum(feature_depths) / len(feature_depths)
                coverage = len([d for d in feature_depths if d > 0]) / len(feature_depths) * 100
                locus_tag = utility.get_locus_tag(feature)
                replicon_gene_stats[replicon][locus_tag] = (coverage, depth_mean)
        # Record results
        isolate_id = mpileup_fp.stem.replace('_mpileup', '')
        gene_stats[isolate_id] = replicon_gene_stats

    # Output data
    coverage_fp = args.output_dir / 'gene_coverage.tsv'
    depth_fp = args.output_dir / 'gene_depth.tsv'
    with contextlib.ExitStack() as stack:
        # Get filehandles and allow exit stack to handle closing
        coverage_fh = stack.enter_context(coverage_fp.open('w'))
        depth_fh = stack.enter_context(depth_fp.open('w'))
        # Write results
        print('replicon', 'locus_tag', *gene_stats, sep=',', file=coverage_fh)
        print('replicon', 'locus_tag', *gene_stats, sep=',', file=depth_fh)
        for replicon, features in replicon_genes.items():
            for feature in features:
                locus_tag = utility.get_locus_tag(feature)
                print(replicon, locus_tag, sep=',', end='', file=coverage_fh)
                print(replicon, locus_tag, sep=',', end='', file=depth_fh)
                for isolate_id in gene_stats:
                    if replicon in gene_stats[isolate_id]:
                        coverage, depth_mean = gene_stats[isolate_id][replicon][locus_tag]
                        print(f',{round(coverage, 2)}', end='', file=coverage_fh)
                        print(f',{round(depth_mean, 2)}', end='', file=depth_fh)
                    else:
                        print(f',0.00', end='', file=coverage_fh)
                        print(f',0.00', end='', file=depth_fh)
                print(file=coverage_fh)
                print(file=depth_fh)


if __name__ == '__main__':
    main()
