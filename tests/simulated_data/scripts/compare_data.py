#!/usr/bin/env python3
import argparse
import pathlib
import re
import sys


import Bio.SeqIO.FastaIO


import utility


file_patterns = {
    'allele_table': 'alleles.tsv',
    'allele_core_table': 'alleles_core.tsv',
    'consequences': 'consequences_core.tsv',
    'mapping_stats': 'mapping_stats.tsv',
    'snp_alignment': 'core.mfasta',
    'gene_coverage': 'coverage.tsv',
    'gene_depth': 'depth.tsv'
}


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_dir', required=True, type=pathlib.Path,
            help='Run directory path')
    parser.add_argument('--spec_fp', required=True, type=pathlib.Path,
            help='Dataset specification filepath')
    parser.add_argument('--test_data_dir', type=pathlib.Path,
            help='Test data directory path')
    args = parser.parse_args()
    if not args.run_dir.exists():
        parser.error('Input directory {args.run_dir} does not exist')
    if not args.spec_fp.exists():
        parser.error('Input file {args.spec_fp} does not exist')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Get data
    spec_info, spec_isolates = utility.read_spec_file(args.spec_fp)
    data_input, data_test = get_data(args.run_dir, args.test_data_dir)

    # Run comparisons on input run directory
    run_spec_data_comparison(spec_isolates, spec_info, data_input)

    # If we have a previous run directory to compare with, do now
    if not args.test_data_dir:
        return
    # NOTE: we assume that all filetypes and replicons are present in the test run
    for filetype in file_patterns.keys():
        if filetype in {'gene_coverage', 'gene_depth'}:
            assert data_input[filetype] == data_test[filetype]
            continue
        for replicon_id in spec_info['replicons']:
            assert data_input[filetype][replicon_id] == data_test[filetype][replicon_id]


def get_data(run_dir, test_data_dir):
    data_input = {filetype: dict() for filetype in file_patterns}
    data_test = {filetype: dict() for filetype in file_patterns}
    for filepath_input in run_dir.iterdir():
        for filetype, file_suffix in file_patterns.items():
            if filepath_input.name.endswith(file_suffix):
                if test_data_dir:
                    filepath_test = test_data_dir / filepath_input.name
                    if not filepath_test.exists():
                        print(f'error: missing file {filepath_test}', file=sys.stderr)
                        sys.exit(1)
                    read_file_data(filepath_test, filetype, file_suffix, data_test)
                read_file_data(filepath_input, filetype, file_suffix, data_input)
                break
    return data_input, data_test


def read_file_data(filepath, filetype, suffix, run_data):
    # NOTE: we modify run_data in place to simplify code
    replicon_id_re = re.compile(f'reference_(.+?)_{suffix}')
    re_result = replicon_id_re.match(filepath.name)
    assert re_result
    replicon_id = re_result.group(1)
    if filetype not in {'gene_coverage', 'gene_depth'}:
        assert replicon_id not in run_data[filetype]
    if filetype in {'allele_table', 'allele_core_table'}:
        run_data[filetype][replicon_id] = read_allele_table(filepath)
    elif filetype == 'consequences':
        run_data[filetype][replicon_id] = read_consequences(filepath)
    elif filetype == 'mapping_stats':
        run_data[filetype][replicon_id] = read_mapping_stats(filepath)
    elif filetype == 'snp_alignment':
        run_data[filetype][replicon_id] = read_snp_alignment(filepath)
    elif filetype in {'gene_coverage', 'gene_depth'}:
        run_data[filetype] = read_gene_coverage_depth(filepath)
    else:
        raise ValueError(f'got bad filetype {filetype}')


def read_allele_table(filepath):
    with filepath.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        header_tokens = next(line_token_gen)
        isolates = header_tokens[2:]
        allele_table = {isolate: dict() for isolate in isolates}
        for position_str, reference_allele, *alleles in line_token_gen:
            position = int(position_str)
            for isolate, allele in zip(isolates, alleles):
                allele_table[isolate][position] = allele
    return allele_table


def read_mapping_stats(filepath):
    stats = dict()
    with filepath.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        header_tokens = next(line_token_gen)
        metrics = header_tokens[1:]
        for isolate, *data_values in line_token_gen:
            assert isolate not in stats
            data = {var: val for var, val in zip(metrics, data_values)}
            stats[isolate] = data
    return stats


def read_consequences(filepath):
    consequences = dict()
    with filepath.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        header_tokens = next(line_token_gen)
        for line_tokens in line_token_gen:
            # We currently have to skip intergenic
            if line_tokens[3] == 'intergenic':
                continue
            data = {var: val for var, val in zip(header_tokens, line_tokens)}
            data['position'] = int(data['position'])
            # Must order affected isolates
            data['isolates'] = sorted(data['isolates'].split(','))
            for isolate in data['isolates']:
                if isolate not in consequences:
                    consequences[isolate] = dict()
                assert data['position'] not in consequences[isolate]
                consequences[isolate][data['position']] = data
    return consequences


def read_snp_alignment(filepath):
    with filepath.open('r') as fh:
        return {desc: seq for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh)}


def read_gene_coverage_depth(filepath):
    data = dict()
    with filepath.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        header_tokens = next(line_token_gen)
        isolates = header_tokens[2:]
        for replicon_id, gene, *values in line_token_gen:
            if replicon_id not in data:
                data[replicon_id] = dict()
            if gene not in data[replicon_id]:
                data[replicon_id][gene] = dict()
            for isolate, value in zip(isolates, values):
                assert isolate not in data[replicon_id][gene]
                data[replicon_id][gene][isolate] = value
    return data


def run_spec_data_comparison(spec_isolates, spec_info, data_input):
    # Unpack some data for brevity
    allele_tables = data_input['allele_table']
    consequence_data = data_input['consequences']
    mapping_stats = data_input['mapping_stats']

    mutation_re = re.compile(r'^(?P<gene>[a-zA-Z0-9]+)_(?P<ref>[A-Z*])(?P<pos>[0-9]+)(?P<alt>[A-Z*]).*$')
    for isolate_data in spec_isolates.values():
        # Mapping stats output
        for replicon_id in spec_info['replicons']:
            stats_unmapped = 100 - float(mapping_stats[replicon_id][isolate_data.name]['total_mapped_reads'])
            stats_mean_depth = float(mapping_stats[replicon_id][isolate_data.name]['replicon_average_depth'])
            assert abs(stats_mean_depth - isolate_data.mean_depth) < 2
            assert abs(stats_unmapped - isolate_data.unmapped * 100) < 1
        # Consequence data and allele table output
        for replicon_id, variant_data in isolate_data.variants.items():
            # Homozygous SNPs
            for snp_data in variant_data.homs:
                # Do not compare intergenic or SNPs that are filtered
                if snp_data['note'] == 'intergenic':
                    continue
                if snp_data['note'].endswith('filtered'):
                    continue
                # Get consequence data
                position = int(snp_data['position'])
                cons_data = consequence_data[replicon_id][isolate_data.name][position]
                # Get mutation data from spec file
                re_result = mutation_re.match(snp_data['note'])
                assert re_result
                mutation = re_result.groupdict()
                mutation_type = 'non-synonymous' if mutation['ref'] != mutation['alt'] else 'synonymous'
                # Compare mutation spec data and consequence output data
                assert cons_data['change_type'] == mutation_type
                assert snp_data['alt'].upper() == cons_data['alt']
                assert mutation['ref'] == cons_data['ref_aa']
                assert mutation['alt'] == cons_data['alt_aa']
                assert mutation['pos'] == cons_data['gene_codon_position']
                # Compare spec data to allele table output data
                assert snp_data['alt'].upper() == allele_tables[replicon_id][isolate_data.name][position]
            # Unknown alleles/ missing data ('-')
            for snp_data in variant_data.lquals:
                position = int(snp_data['position'])
                assert '-' == allele_tables[replicon_id][isolate_data.name][position]


if __name__ == '__main__':
    main()
