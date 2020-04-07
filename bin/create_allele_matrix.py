#!/usr/bin/env python3
import argparse
import pathlib
import sys


from Bio.SeqIO.QualityIO import FastqPhredIterator


vcf_fields = {'chrom': str,
              'pos': int,
              'id': str,
              'ref': str,
              'alt': str,
              'qual': float,
              'filter': str,
              'info': str,
              'format': str,
              'genotype': str}


class VcfRecord:

    def __init__(self, data):
        for value, (attr, attr_type) in zip(data, vcf_fields.items()):
            setattr(self, attr, attr_type(value))
        self._data = data

    def __str__(self):
        return '\t'.join(self._data)


class CheckInput(argparse.Action):

    def __call__(self, parser, namespace, value, option_string=None):
        value_check = [value] if isinstance(value, pathlib.Path) else value
        for filepath in value_check:
            if not filepath.exists():
                parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, value)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_fps', required=True, type=pathlib.Path,
            help='Input VCF filepaths', nargs='+', action=CheckInput)
    parser.add_argument('--consensus_fps', required=True, type=pathlib.Path,
            help='Input consensus filepaths', nargs='+', action=CheckInput)
    parser.add_argument('--replicon', required=True, type=str,
            help='Replicon name')
    return parser.parse_args()


def main():
    # Get command line arguments
    args = get_arguments()

    # Get high quality SNP positions from filtered VCFs
    snp_positions = dict()
    replicons_seen = set()
    for vcf_fp in args.vcf_fps:
        vcf_alleles = dict()
        with vcf_fp.open('r') as fh:
            # Skip header
            while fh.readline().startswith('#'):
                position = fh.tell()

            # Check if there is any data to process
            if fh.readline():
                fh.seek(position, 0)
            else:
                break

            # Collect SNPs
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            record_gen = (VcfRecord(line_tokens) for line_tokens in line_token_gen)
            for record in record_gen:
                # Skip INDELs
                if record.info.startswith('I'):
                    continue
                # Skip non-SNPs
                if ',' in record.alt:
                    continue
                # Skip when ref and alt are the same
                if record.ref == record.alt:
                    continue
                # Note SNPs
                if record.pos in snp_positions and snp_positions[record.pos] != record.ref:
                    print(f'got mismatch reference allele for position {record.pos}', file=sys.stderr)
                    sys.exit(1)
                else:
                    snp_positions[record.pos] = record.ref

            # Ensure that we're getting the same replicon
            replicons_seen.add(record.chrom)
            if len(replicons_seen) != 1:
                print('error: got duplicate replicon names', file=sys.stderr)
                sys.exit(1)


    # Read in consensus sequences and get appropriate replicon
    consensus = dict()
    for consensus_fp in args.consensus_fps:
        sample_id = consensus_fp.name.replace('_consensus.fastq', '')
        with consensus_fp.open('r') as fh:
            for record in FastqPhredIterator(fh):
                if record.id == args.replicon:
                    consensus[sample_id] = record
                    break
            else:
                msg = f'error: could not find replicon {args.replicon} in {consensus_fp}'
                print(msg, file=sys.stderr)
                sys.exit(1)


    # Print SNP calls from consensus sequence (called from raw VCF)
    # The consensus sequence gives us access to low quality calls (allowing us to make missing data
    # designations) and also access to the high quality SNP calls
    print('Position', 'Reference', *consensus, sep='\t')
    for snp_position in sorted(snp_positions):
        ref_allele = snp_positions[snp_position]
        print(snp_position, ref_allele, sep='\t', end='')
        for sample_id in consensus:
            record = consensus[sample_id]
            snp = record.seq[snp_position-1]
            quality = record.letter_annotations['phred_quality'][snp_position-1]
            # Exclude low quality SNPs
            if quality < 20:
                print('\t', '-', sep='', end='')
            # vcfutils vcf2fq uses lowercase nucleotides for low quailty calls
            # Exclude those here as well
            elif snp not in {'A', 'T', 'G', 'C'}:
                print('\t', '-', sep='', end='')
            else:
                print('\t', snp, sep='', end='')
        print()


if __name__ == '__main__':
    main()
