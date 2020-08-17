#!/usr/bin/env python3
import argparse
import collections
import pathlib
import re
import sys


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--src_dir', required=True, type=pathlib.Path,
            help='Merge source directory')
    parser.add_argument('--dst_dir', required=True, type=pathlib.Path,
            help='New output directory')
    parser.add_argument('--reference_fp', required=True, type=pathlib.Path,
            help='Reference filepath')
    # We take the globs here as they're available immediately (avoids issues with nf async)
    parser.add_argument('--read_globs', required=True, nargs='+', type=pathlib.Path,
            help='Read globs provided to nextflow')
    args = parser.parse_args()
    if not args.src_dir.exists():
        parser.error(f'Input merge directory {args.src_dir} does not exist')
    if not args.dst_dir.exists():
        parser.error(f'New output directory {args.dst_dir} does not exist')
    if not args.reference_fp.exists():
        parser.error(f'Input file {args.reference_fp} does not exist')
    return args


def main():
    # Get commandline arguments
    args = get_arguments()
    status = {'return_code': 0, 'messages': list()}

    # Validate configuration
    compare_run_configs(args.src_dir, args.dst_dir, status)

    # Get reference data from run config and compare
    merge_ref_name, merge_rep_names = compare_references(args.src_dir, args.dst_dir, status)

    # Using BAMs as ground truth - this could be changed to something else
    # Require at least one BAM
    isolate_names_bam = {fp.stem for fp in args.src_dir.glob('bams/*.bam')}
    if len(isolate_names_bam) < 1:
        status['messages'].append(f'error: found no BAM files in {args.src_dir}/bams/')
        status['return_code'] = 1

    # Compare isolate names between BAMs and VCFs
    isolate_names_vcf = set()
    for vcf_fp in args.src_dir.glob('vcfs/*_q30.vcf'):
        isolate_names_vcf.add(vcf_fp.stem.replace('_q30', ''))

    if isolate_names_bam != isolate_names_vcf:
        status['messages'].append('error: isolate names for BAMs and VCFs differ')
        status['return_code'] = 1

    # Compare to fastqc if we have any
    isolate_names_fqc = set()
    fastqc_re = re.compile(r'^(.+)(?:_001_)?(?:_R?[12])(?:_subsampled)?_fastqc.zip$')
    for fastqc_fp in args.src_dir.glob('fastqc/*/*.zip'):
        re_result = fastqc_re.match(fastqc_fp.name)
        assert re_result
        isolate_names_fqc.add(re_result.group(1))
    if isolate_names_fqc and isolate_names_fqc != isolate_names_bam:
        status['messages'].append('error: isolate names in FastQC reports and BAMs differ')
        status['return_code'] = 1

    # Check isolate names from gene coverage and depth
    gene_depth_fp = args.src_dir / f'{merge_ref_name}_gene_depth.tsv'
    gene_coverage_fp = args.src_dir / f'{merge_ref_name}_gene_coverage.tsv'
    if not gene_depth_fp.exists():
        status['messages'].append(f'error: gene depth file {gene_depth_fp} does not exist')
        status['return_code'] = 1
        isolate_names_depth = set()
    else:
        with gene_depth_fp.open('r') as fh:
            isolate_names_depth = set(fh.readline().rstrip().split(',')[2:])
    if not gene_coverage_fp.exists():
        status['messages'].append(f'error: gene coverage file {gene_coverage_fp} does not exist')
        status['return_code'] = 1
        isolate_names_coverage = set()
    else:
        with gene_coverage_fp.open('r') as fh:
            isolate_names_coverage = set(fh.readline().rstrip().split(',')[2:])
    if isolate_names_depth and isolate_names_depth != isolate_names_bam:
        status['messages'].append('error: isolate names in gene depth and BAMs differ')
        status['return_code'] = 1
    if isolate_names_coverage and isolate_names_coverage != isolate_names_bam:
        status['messages'].append('error: isolate names in coverage depth and BAMs differ')
        status['return_code'] = 1

    # Check isolate names from mapping stats
    isolate_names_stats = get_isolate_names_from_mapping_stats(args.src_dir, merge_rep_names, merge_ref_name, status)
    for replicon_name, isolate_names in isolate_names_stats['all'].items():
        if isolate_names != isolate_names_bam:
            status['messages'].append(f'error: isolate names in mapping stats ({replicon_name}) and BAMs differ')
            status['return_code'] = 1

    # Get isolates from allele matrices and compare
    isolate_names_allele = get_isolate_names_from_allele_matrices(args.src_dir, merge_rep_names, merge_ref_name, status)
    allele_names = set(isolate_names_allele)
    stat_names = set(isolate_names_stats['pass'])
    missing_stats = allele_names.difference(stat_names)
    missing_allele = stat_names.difference(allele_names)
    if missing_stats:
        missing_stats_str = ', '.join(missing_stats)
        status['messages'].append(f'error: missing mapping stats files for {missing_stats_str}')
        status['return_code'] = 1
    if missing_allele:
        missing_allele_str = ','.join(missing_allele)
        status['messages'].append(f'error: missing allele matrix files for {missing_allele_str}')
        status['return_code'] = 1
    # Only compare shared replicons and only check that allele isolate names are a subset of
    # mapping stats isolate names - isolates that pass but have no SNPs will not appear in the
    # allele matrix
    rep_names = allele_names & stat_names
    for name in rep_names:
        missing_isolates = isolate_names_allele[name].difference(isolate_names_stats['pass'][name])
        if missing_isolates:
            missing_isolates_str = ', '.join(missing_isolates)
            status['messages'].append(f'error: missing isolates {missing_isolates_str} in mapping stats for {name}')
            status['return_code'] = 1

    # Check for namespace collisions
    isolate_names_new = get_isolate_names_from_reads(args.read_globs)
    isolate_names_collisions = isolate_names_new & isolate_names_bam
    if isolate_names_collisions:
        names_str = ', '.join(isolate_names_collisions)
        status['messages'].append(f'error: detected isolate name collision for {names_str}')
        status['return_code'] = 1

    # Print all messages and exit
    print_messages_and_exit(status)


def compare_run_configs(src_dir, dst_dir, status):
    merge_data = collect_run_config_data(src_dir, status)
    new_data = collect_run_config_data(dst_dir, status)
    if not new_data or not merge_data:
        status['return_code'] = 1
        return
    shared_vars = set(merge_data) & set(new_data)
    for var in shared_vars:
        if merge_data[var] != new_data[var]:
            message = f'error: merge run and current run have different settings for {var}:'
            message += f' merge: {merge_data[var]}; new: {new_data[var]}'
            status['messages'].append(message)
            status['return_code'] = 1


def collect_run_config_data(src_dir, status):
    val_config = 0
    run_config_fp = src_dir / 'run_info/run_config.tsv'
    if not run_config_fp.exists():
        message = f'error: aborting comparison: run config could not be found at {run_config_fp}'
        status['messages'].append(message)
        status['return_code'] = 1
    entries_skip = {
            'reference_name',
            'reference_replicon_names',
            'reference_replicon_sizes',
            'reference_replicon_md5hashes'
        }
    entries_retain = {
            'bt2_max_frag_len',
            'bt2_mode',
            'var_depth_min',
            'mapping_cover_min',
            'mapping_depth_min',
            'mapping_mapped_min',
            'outgroup_mod',
            'allele_matrix_cons'
        }
    run_config = dict()
    with run_config_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        for var, val in line_token_gen:
            if var in entries_skip:
                continue
            elif var in entries_retain:
                assert var not in run_config
                run_config[var] = val
            else:
                status['messages'].append(f'error: unexpected item in run_config {var}')
                status['return_code'] = 1
    entries_missing = set(run_config) ^ entries_retain
    if entries_missing:
        entries_missing_str = ', '.join(entries_missing)
        status['messages'].append(f'error: {run_config_fp} has missing entries: {entries_missing_str}')
        status['return_code'] = 1
    return run_config


def compare_references(src_dir, dst_dir, status):
    merge_data = collect_reference_info(src_dir, status)
    new_data = collect_reference_info(dst_dir, status)
    if not new_data or not merge_data:
        status['return_code'] = 1
    # Reference name
    if new_data['name'] != merge_data['name']:
        message = f'error: reference name mismatch - merge: {merge_data["name"]}; new: {new_data["name"]}'
        status['messages'].append(message)
        status['return_code'] = 1
    # Replicon names
    new_names = set(new_data['names'])
    merge_names = set(merge_data['names'])
    if new_names != merge_names:
        new_missing = ', '.join(merge_names.difference(new_names))
        merge_missing = ', '.join(new_names.difference(merge_names))
        if new_missing:
            status['messages'].append(f'error: missing contigs in new reference {new_missing}')
        if merge_missing:
            status['messages'].append(f'error: missing contigs in merge reference {merge_missing}')
        status['return_code'] = 1
    # Contig sizes and hashes (only those that are shared)
    rep_names = new_names & merge_names
    for name in rep_names:
        new_size = new_data['sizes'][name]
        merge_size = merge_data['sizes'][name]
        if new_size != merge_size:
            status['messages'].append(f'error: replicon size differs for {name} - merge {merge_size}; new {new_size}')
            status['return_code'] = 1
            continue
        new_hash = new_data['hashes'][name]
        merge_hash = merge_data['hashes'][name]
        if new_hash != merge_hash:
            message = f'error: replicon hash (i.e. seq) differs for {name} - merge {merge_hash}; new {new_hash}'
            status['messages'].append(message)
            status['return_code'] = 1
    return merge_data['name'], merge_data['names']


def collect_reference_info(data_dir, status):
    run_config_fp = data_dir / 'run_info/run_config.tsv'
    if not run_config_fp.exists():
        status['messages'].append(f'error: fatal: run config could not be found at {run_config_fp}')
        status['return_code'] = 1
        print_messages_and_exit(status)
    # Using different top level vars for robustness
    name = str()
    names = set()
    sizes = dict()
    hashes = dict()
    with run_config_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        for var, val in line_token_gen:
            if not var.startswith('reference'):
                continue
            if var == 'reference_name':
                assert not name
                name = val
            elif var == 'reference_replicon_names':
                for replicon_name in val.split(','):
                    assert replicon_name not in names
                    names.add(replicon_name)
            elif var == 'reference_replicon_sizes':
                for size_token in val.split(','):
                    replicon_name, replicon_size = size_token.split(':')
                    assert replicon_name not in sizes
                    sizes[replicon_name] = replicon_size
            elif var == 'reference_replicon_md5hashes':
                for hash_token in val.split(','):
                    replicon_name, replicon_hash = hash_token.split(':')
                    assert replicon_name not in hashes
                    hashes[replicon_name] = replicon_hash
    return {'name': name, 'names': names, 'sizes': sizes, 'hashes': hashes}


def get_isolate_names_from_mapping_stats(src_dir, merge_rep_names, merge_ref_name, status):
    fps_count = 0
    isolate_names_all = dict()
    isolate_names_pass = dict()
    for name in merge_rep_names:
        fp = src_dir / f'{merge_ref_name}_{name}_mapping_stats.tsv'
        if not fp.exists():
            continue
        fps_count += 1
        assert name not in isolate_names_all
        assert name not in isolate_names_pass
        isolate_names_all[name] = set()
        isolate_names_pass[name] = set()
        # Get isolate names
        with fp.open('r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            header_tokens = next(line_token_gen)
            Record = collections.namedtuple('Record', [t.lower() for t in header_tokens])
            for line_tokens in line_token_gen:
                record = Record(*line_tokens)
                assert record.isolate not in isolate_names_all[name]
                assert record.isolate not in isolate_names_pass[name]
                isolate_names_all[name].add(record.isolate)
                if record.pass_fail == 'p':
                    isolate_names_pass[name].add(record.isolate)
    if fps_count < 1:
        status['messages'].append('error: did not find any mapping stats file in merge directory')
        status['return_code'] = 1
    else:
        rep_names = list(isolate_names_all.keys())
        for i in range(len(rep_names)):
            for j in range(i+1, len(rep_names)):
                name_i = rep_names[i]
                name_j = rep_names[j]
                isolates_i = isolate_names_all[name_i]
                isolates_j = isolate_names_all[name_j]
                if isolates_i != isolates_j:
                    status['messages'].append(f'error: isolates in mapping stats for {name_i} and {name_j} differ')
                    status['return_code'] = 1
    return {'all': isolate_names_all, 'pass': isolate_names_pass}


def get_isolate_names_from_allele_matrices(src_dir, merge_rep_names, merge_ref_name, status):
    fps_count = 0
    isolate_names = dict()
    for name in merge_rep_names:
        fp = src_dir / f'{merge_ref_name}_{name}_alleles.tsv'
        if not fp.exists():
            continue
        fps_count += 1
        assert name not in isolate_names
        # Get isolate names
        with fp.open('r') as fh:
            isolates = set(fh.readline().rstrip().split(',')[2:])
            isolate_names[name] = isolates
    if fps_count < 1:
        status['messages'].append('error: did not find any allele matrices file in merge directory')
        status['return_code'] = 1
    return isolate_names


def get_isolate_names_from_reads(read_globs):
    # Discover reads
    read_fps = list()
    for read_glob in read_globs:
        read_fps.extend(read_glob.parent.glob(read_glob.name))
    # Get read prefix for each
    regexes = (re.compile(r'^(.+?)_(?:_001_)?R?[12].fastq(?:.gz)?$'), re.compile(r'^(.+?).fastq(?:.gz)?$'))
    isolate_names_new = set()
    for read_fp in read_fps:
        for regex in regexes:
            re_result = regex.match(read_fp.name)
            if re_result:
                break
        else:
            assert False
        isolate_names_new.add(re_result.group(1))
    return isolate_names_new


def print_messages_and_exit(status):
    for message in status['messages']:
        print(message, file=sys.stderr)
    sys.exit(status['return_code'])


if __name__ == '__main__':
    main()
