#!/usr/bin/env python3
import argparse
import pathlib
import sys


import Bio.Alphabet.IUPAC
import Bio.SeqFeature
import Bio.SeqIO


import utility


nucleotide_complement = str.maketrans('atgcATGC', 'tacgTACG')
ambiguous_amino_acids = {'B', 'J', 'X', 'Z'}


class Interval:

    def __init__(self, feature, start, end):
        self.feature = feature
        self.start = start
        self.end = end

    @classmethod
    def from_feature(cls, feature):
        start = feature.location.nofuzzy_start + 1
        end = feature.location.nofuzzy_end
        return Interval(feature, start, end)


class Node:

    def __init__(self):
        self.center_start = int()
        self.center_set = set()
        self.node_left = list()
        self.node_right = list()

    def search_position(self, position, result=None):
        if not result:
            result = list()
        # Check if interval is in current node
        for interval in self.center_set:
            if interval.start <= position <= interval.end:
                result.append(interval)
        # Recursively check child nodes
        if position < self.center_start and self.node_left:
            return self.node_left.search_position(position, result)
        elif position > self.center_start and self.node_right:
            return self.node_right.search_position(position, result)
        return result


class AlleleRecord:

    def __init__(self, record_data, isolates):
        position_str, reference, *alleles = record_data
        self.position = int(position_str)
        self.reference = reference
        self.alleles = alleles
        self.isolates = isolates

        self.feature_intervals = list()


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference_fp', required=True, type=pathlib.Path,
            help='Genbank-format reference filepath')
    parser.add_argument('--allele_fp', required=True, type=pathlib.Path,
            help='Allele matrix filepath')
    parser.add_argument('--replicon', required=True, type=str,
            help='Name of replicon')
    args = parser.parse_args()
    if not args.reference_fp.exists():
        parser.error(f'Input file {args.reference_fp} does not exist')
    if not args.allele_fp.exists():
        parser.error(f'Input file {args.allele_fp} does not exist')
    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in replicon from reference
    ref_gbk = {record.id: record for record in Bio.SeqIO.parse(args.reference_fp, 'genbank')}
    if args.replicon not in ref_gbk:
        print(f'error: could not find {args.replicon} in {args.reference_fp}', file=sys.stderr)
        sys.exit(1)
    rep_gbk = ref_gbk[args.replicon]

    # Create interval tree for features
    # Tree is assumed to be balanced, which it should be for this type of data
    intervals = list()
    for feature in rep_gbk.features:
        if feature.type != 'CDS':
            continue
        intervals.append(Interval.from_feature(feature))
    interval_tree = create_interval_tree(intervals)

    # Determine consequences
    header_tokens = (
        'position', 'ref', 'alt', 'change_type', 'gene', 'ref_codon', 'alt_codon', 'ref_aa', 'alt_aa',
        'gene_product', 'gene_nucleotide_position', 'gene_codon_position', 'codon_nucleotide_position',
        'notes'
    )
    print(*header_tokens, sep='\t')
    with args.allele_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split(',') for line in fh)
        header_tokens = next(line_token_gen)
        isolates = header_tokens[2:]
        # Iterate through each entry via AlleleRecords
        allele_record_gen = (AlleleRecord(line_tokens, isolates) for line_tokens in line_token_gen)
        for record in allele_record_gen:
            # Find feature intervals that this record lands in
            record.feature_intervals = interval_tree.search_position(record.position)
            # Process intergenic
            if not record.feature_intervals:
                for allele in get_alleles(record):
                    print(
                        record.position, record.reference, allele, 'intergenic',
                        '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', sep='\t'
                    )
                continue
            # Process genic
            for interval in record.feature_intervals:
                # Check strand is acceptable
                strand = interval.feature.strand
                if not (strand == 1 or strand == -1):
                    print(f'error: got bad strand from {feature}', file=sys.stderr)
                    sys.exit(1)
                # Get SNP location within the gene and codon
                if strand == 1:
                    position_gene = record.position - interval.start + 1
                elif strand == -1:
                    position_gene = interval.end - record.position + 1
                position_codon = (position_gene - 1) % 3 + 1
                # Get the codon which contains the SNP
                codon_number = (position_gene - 1) // 3 + 1
                codon_start = (codon_number - 1) * 3
                # Collect gene and codon sequence
                gene_sequence = interval.feature.extract(rep_gbk.seq)
                codon_sequence = gene_sequence[codon_start:codon_start+3]
                codon = codon_sequence.translate()
                # Get consequences for each allele
                for allele in get_alleles(record):
                    # Skip alleles that fall within compound locations
                    if isinstance(interval.feature.location, Bio.SeqFeature.CompoundLocation):
                        locus_tag = utility.get_locus_tag(interval.feature)
                        msg = (f'warning: skipping position {record.position} in gene {locus_tag} '
                                'as the gene has a compound location')
                        print(msg, file=sys.stderr)
                        print(
                            record.position, record.reference, allele, 'not assessed',
                            '-', '-', '-', '-', '-', '-', '-', '-', '-',
                            f'allele falls within compound location of {locus_tag}',
                            sep='\t'
                        )
                    else:
                        # Must complement allele if on reverse strand
                        # The ref and alt alleles relative to coding strand
                        # Both codon and amino acid details are reported on the appropriate strand
                        # Get consequence
                        allele_cns = allele if strand == 1 else allele.translate(nucleotide_complement)
                        sequence_allele, codon_allele = get_consequence(codon_sequence, allele_cns, position_codon)
                        if codon == codon_allele:
                            change_type = 'synonymous'
                        elif codon in ambiguous_amino_acids or codon_allele in ambiguous_amino_acids:
                            change_type = 'ambiguous'
                        else:
                            change_type = 'non-synonymous'
                        # Get feature info
                        if 'product' in interval.feature.qualifiers:
                            [gene_product] = interval.feature.qualifiers['product']
                        else:
                            gene_product = '-'
                        locus_tag = utility.get_locus_tag(interval.feature)
                        print(
                            record.position, record.reference, allele, change_type, locus_tag, codon_sequence,
                            sequence_allele, codon, codon_allele, gene_product, position_gene, codon_number,
                            position_codon, '-', sep='\t'
                        )


def get_alleles(record):
    alleles = set(record.alleles)
    # Remove alleles that are unknown
    if '-' in alleles:
        alleles.remove('-')
    # Remove alleles that match reference
    if record.reference in alleles:
        alleles.remove(record.reference)
    return alleles


def get_consequence(codon_sequence, allele, position_codon):
    sequence_list = list(codon_sequence)
    sequence_list[position_codon-1] = allele
    sequence = Bio.Seq.Seq(''.join(sequence_list), Bio.Alphabet.IUPAC.ambiguous_dna)
    return sequence, sequence.translate()


def create_interval_tree(intervals):
    if not intervals:
        return None
    node = Node()
    return add_intervals(node, intervals)


def add_intervals(node, intervals):
    # Sort intervals into left or right of center
    center_interval = intervals[len(intervals) // 2]
    node.center_start = center_interval.start
    intervals_left = list()
    intervals_right = list()
    for interval in intervals:
        if interval.end < node.center_start:
            intervals_left.append(interval)
        elif interval.start > node.center_start:
            intervals_right.append(interval)
        else:
            node.center_set.add(interval)
    # Recursively create subtrees
    node.node_left = create_interval_tree(intervals_left)
    node.node_right = create_interval_tree(intervals_right)
    return node


if __name__ == '__main__':
    main()
