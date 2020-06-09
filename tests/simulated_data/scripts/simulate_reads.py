#!/usr/bin/env python3
import argparse
import copy
import contextlib
import math
import pathlib
import sys
import tempfile


import Bio.SeqIO
import Bio.Seq


import utility


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference_fp', required=True, type=pathlib.Path,
            help='Genbank format reference filepath')
    parser.add_argument('--spec_fp', required=True, type=pathlib.Path,
            help='Dataset specification filepath')
    parser.add_argument('--output_dir', required=True, type=pathlib.Path,
            help='Output directory')
    args = parser.parse_args()
    if not args.reference_fp.exists():
        parser.error('Input file {args.reference_fp} does not exist')
    if not args.spec_fp.exists():
        parser.error('Input file {args.spec_fp} does not exist')
    if not args.output_dir.exists():
        parser.error('Input file {args.output_dir} does not exist')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Parse readset specification table
    spec_info, readsets = utility.read_spec_file(args.spec_fp)
    # Unpacking for brevity later
    read_len = spec_info['read_length']
    outer_len = spec_info['outer_length']

    # TODO: fix issue where ends of contigs have low mean depth
    # One solution is to add n nucleotides from either end to the other
    # Where n is the size of the read to be simulated

    # Read in reference sequence
    replicon_seqs = dict()
    with args.reference_fp.open('r') as fh:
        for record in Bio.SeqIO.parse(fh, 'genbank'):
            mutable_seq = Bio.Seq.MutableSeq(str(record.seq))
            replicon_seqs[record.name] = mutable_seq

    # Simulate data
    for readset in readsets.values():
        isolate_seqs = copy.deepcopy(replicon_seqs)
        # Add homozygous SNPs and INDELs to each replicon sequence
        # Also arrange heterozygous SNPs and low quality sites by position
        # These variants will be introduced during read simulation
        variant_sites = {replicon: dict() for replicon in replicon_seqs}
        for replicon, variants in readset.variants.items():
            # Modify reference sequence
            for var_hom in variants.homs:
                position = var_hom['position']
                isolate_seqs[replicon][position-1] = var_hom['alt']
            for var_indel in variants.indels:
                position = var_indel['position']
                if var_indel['op'] == 'ins':
                    isolate_seqs[replicon].insert(position-1, var_indel['alt'])
                elif var_indel['op'] == 'del':
                    isolate_seqs[replicon].pop(position-1)
            # Collect heterozygous SNPs and low quality sites by position
            variant_sites[replicon] = dict()
            for var_het in variants.hets:
                position = var_het['position']
                if position in variant_sites[replicon]:
                    raise ValueError(f'got more than one modification at {position} in {replicon}')
                variant_sites[replicon][position] = ('het', var_het)
            for var_lqual in variants.lquals:
                position = var_lqual['position']
                if position in variant_sites[replicon]:
                    raise ValueError(f'got more than one modification at {position} in {replicon}')
                variant_sites[replicon][position] = ('lqual', var_lqual)

        # Output files
        forward_fp = args.output_dir / f'{readset.name}_1.fastq'
        reverse_fp = args.output_dir / f'{readset.name}_2.fastq'

        # Simulate reads
        read_number = 0
        with contextlib.ExitStack() as stack:
            forward_fh = stack.enter_context(forward_fp.open('w'))
            reverse_fh = stack.enter_context(reverse_fp.open('w'))
            for replicon, seq in isolate_seqs.items():
                # Record hets generated at each site
                het_counts = dict()
                # Generate steps for read positions without going past the end
                # Also ensure the last read does include the final nucleotide of contig
                step = math.floor(read_len / readset.mean_depth * 2)
                def start_gen():
                    yield from range(0, len(seq)-outer_len, step)
                    yield from [len(seq)-outer_len]
                # Iterate read start positions
                for start in start_gen():
                    read_number += 1
                    end = start + outer_len
                    # Forward
                    desc = f'@{replicon}_{start+1}_{end}_0:0:0_0:0:0_{read_number}/1'
                    read_seq = seq[start:start+read_len]
                    read_seq, read_qual = add_mods(read_seq, start, start+read_len, het_counts, variant_sites[replicon])
                    print(desc, read_seq, '+', read_qual, sep='\n', file=forward_fh)
                    # Reverse
                    desc = f'@{replicon}_{start+1}_{end}_0:0:0_0:0:0_{read_number}/2'
                    read_seq = seq[end-read_len:end]
                    read_seq, read_qual = add_mods(read_seq, end-read_len, end, het_counts, variant_sites[replicon])
                    read_seq = Bio.Seq.Seq(str(read_seq)).reverse_complement()
                    print(desc, read_seq, '+', read_qual, sep='\n', file=reverse_fh)
            # Create unmappable reads if requested
            if not readset.unmapped:
                continue
            reads_unmapped = round(read_number / (1 - readset.unmapped) * readset.unmapped)
            if reads_unmapped >= 200000:
                print('error: refusing to generate more than 200,000 unmappable reads', file=sys.stderr)
                sys.exit(1)
            for i in range(reads_unmapped):
                read_number += 1
                desc = f'@none_none_none_0:0:0_0:0:0_{read_number}/1'
                print(desc, 'N'*read_len, '+', 'I'*read_len, sep='\n', file=forward_fh)
                desc = f'@none_none_none_0:0:0_0:0:0_{read_number}/2'
                print(desc, 'N'*read_len, '+', 'I'*read_len, sep='\n', file=reverse_fh)


def add_mods(read_seq, start, end, het_counts, variant_sites):
    # Create base read quality
    quality = ['I'] * len(read_seq)
    # Apply other modifications if needed
    overlaps = set(range(start+1, end+1)) & set(variant_sites)
    for position in overlaps:
        var_type, variant = variant_sites[position]
        position_seq = position - start - 1
        if var_type == 'het':
            if position not in het_counts:
                het_counts[position] = {'alt_1': 0, 'alt_2': 0}
            nom = het_counts[position]['alt_1']
            # Catch when denominating is zero
            # Otherwise add appropriate het SNP to balance ratio
            ratio = float(variant['ratio'])
            if nom == 0:
                read_seq[position_seq] = variant['alt_1']
                het_counts[position]['alt_1'] += 1
            elif nom / sum(het_counts[position].values()) >= ratio:
                read_seq[position_seq] = variant['alt_2']
                het_counts[position]['alt_2'] += 1
            else:
                read_seq[position_seq] = variant['alt_1']
                het_counts[position]['alt_1'] += 1
        elif var_type == 'lqual':
            read_seq[position_seq] = 'n'
            quality[position_seq] = '!'
        else:
            raise ValueError(f'got bad modification type {var_type}')
    # TODO: check if read_seq is mutable and modifed in scopes above
    return read_seq, ''.join(quality)


if __name__ == '__main__':
    main()
