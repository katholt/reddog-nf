#!/usr/bin/env python3
import argparse
import pathlib
import sys


import Bio.Alphabet.IUPAC
import Bio.SeqFeature
import Bio.SeqIO


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


class Consequence:

    # Mapping ordered header fields to class attributes
    # header name -> attribute
    fields = {
            'position': 'position',
            'ref': 'reference',
            'alt': 'allele',
            'change_type': 'change_type',
            'gene': 'gene_id',
            'ref_codon': 'sequence',
            'alt_codon': 'sequence_allele',
            'ref_aa': 'codon',
            'alt_aa': 'codon_allele',
            'gene_product': 'gene_product',
            'gene_nucleotide_position': 'gene_nucleotide_position',
            'gene_codon_position': 'gene_codon_position',
            'codon_nucleotide_position': 'codon_nucleotide_position',
            'isolates': 'isolate_str',
            'notes': 'notes',
            }

    def __init__(self):
        for field in self.fields.values():
            setattr(self, field, str())

    def __str__(self):
        data_gen = (getattr(self, field) for field in self.fields.values())
        return '\t'.join(str(d) for d in data_gen)

    @classmethod
    def from_allele_record(cls, record):
        consequence = Consequence()
        consequence.position = record.position
        consequence.reference = record.reference
        return consequence


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

    # Read in reference
    ref_gbk = {record.id: record for record in Bio.SeqIO.parse(args.reference_fp, 'genbank')}

    # Process for replicon
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
    print(*Consequence.fields, sep='\t')
    with args.allele_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        header_tokens = next(line_token_gen)
        isolates = header_tokens[2:]
        # Iterate through each entry via AlleleRecords
        allele_record_gen = (AlleleRecord(line_tokens, isolates) for line_tokens in line_token_gen)
        for record in allele_record_gen:
            # Find feature intervals that this record lands in
            record.feature_intervals = interval_tree.search_position(record.position)
            # Process intergenic
            if not record.feature_intervals:
                for allele, allele_isolates in isolates_by_allele(record).items():
                    # Create new consequence object to hold output data for nicer printing imo
                    consequence = Consequence.from_allele_record(record)
                    consequence.allele = allele
                    consequence.change_type = 'intergenic'
                    print(consequence)
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
                for allele, allele_isolates in isolates_by_allele(record).items():
                    # Create new consequence object to hold output data for nicer printing imo
                    # Doubles to makes intergenic block much cleaner
                    consequence = Consequence.from_allele_record(record)
                    consequence.allele = allele

                    # Must complement allele if on reverse strand
                    # The ref and alt alleles relative to coding strand
                    # Both codon and amino acid details are reported on the appropriate strand
                    if strand == -1:
                        allele = allele.translate(nucleotide_complement)

                    # Skip alleles that fall within compound locations
                    if isinstance(interval.feature.location, Bio.SeqFeature.CompoundLocation):
                        [locus_tag] = interval.feature.qualifiers['locus_tag']
                        msg = (f'warning: skipping position {record.position} in gene {locus_tag} '
                                'as the gene has a compound location')
                        print(msg, file=sys.stderr)
                        consequence.change_type = 'not assessed'
                        consequence.notes = f'allele falls within compound location of {locus_tag}'
                        print(consequence)
                    else:
                        # Get consequence
                        sequence_allele, codon_allele = get_consequence(codon_sequence, allele, position_codon)
                        consequence.codon = codon
                        consequence.codon_allele = codon_allele
                        consequence.sequence = codon_sequence
                        consequence.sequence_allele = sequence_allele
                        consequence.gene_nucleotide_position = position_gene
                        consequence.gene_codon_position = codon_number
                        consequence.codon_nucleotide_position = position_codon
                        if codon == codon_allele:
                            consequence.change_type = 'synonymous'
                        elif codon in ambiguous_amino_acids or codon_allele in ambiguous_amino_acids:
                            consequence.change_type = 'ambiguous'
                        else:
                            consequence.change_type = 'non-synonymous'
                        consequence.isolate_str = ','.join(allele_isolates)
                        [consequence.gene_id] = interval.feature.qualifiers['locus_tag']
                        [consequence.gene_product] = interval.feature.qualifiers['product']
                        # Print data
                        print(consequence)


def isolates_by_allele(record):
    ret = dict()
    for isolate, allele in zip(record.isolates, record.alleles):
        # Skip alleles that are unknown
        if allele == '-':
            continue
        # Skip alleles that match reference
        if allele == record.reference:
            continue
        if allele not in ret:
            ret[allele] = list()
        ret[allele].append(isolate)
    return ret


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
        if interval.end <= node.center_start:
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
