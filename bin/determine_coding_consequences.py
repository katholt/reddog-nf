#!/usr/bin/env python3
# TODO: deal with '-' and 'N'


import argparse
import pathlib
import sys


import Bio.SeqIO
import Bio.Alphabet.IUPAC


nucleotide_complement = str.maketrans('atgcATGC', 'tacgTACG')


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

    def __repr__(self):
        self_type = type(self)
        return f'<{self_type.__qualname__} at {hex(id(self))} with {self.start}:{self.end}>'


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

    fields = {
            'SNP': 'position',
            'ref': 'reference',
            'alt': 'allele',
            'change': 'change_type',
            'gene': 'gene_id',
            'ancestralCodon': 'sequence',
            'derivedCodon': 'sequence_allele',
            'ancestralAA': 'codon',
            'derivedAA': 'codon_allele',
            'product': 'gene_product',
            'isolate_str': 'isolate_str'
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


class CheckInput(argparse.Action):

    def __call__(self, parser, namespace, filepath, option_string=None):
        if not filepath.exists():
            parser.error(f'Filepath {filepath} for {option_string} does not exist')
        setattr(namespace, self.dest, filepath)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference_fp', required=True, type=pathlib.Path,
            help='Genbank-format reference filepath', action=CheckInput)
    parser.add_argument('--allele_matrix', required=True, type=pathlib.Path,
            help='Allele matrix filepath', action=CheckInput)
    return parser.parse_args()


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in reference
    ref_gbk = {record.name: record for record in Bio.SeqIO.parse(args.reference_fp, 'genbank')}

    # Process for replicon
    replicon = 'NZ_UHGC01000001'
    rep_gbk = ref_gbk[replicon]

    # Create interval tree for features
    # NOTE: we assume tree is balanced, which it should be for this type of data
    feature_cds_gen = (feature for feature in rep_gbk.features if feature.type == 'CDS')
    intervals = [Interval.from_feature(feature) for feature in feature_cds_gen]
    interval_tree = create_interval_tree(intervals)

    # Determine consequences
    print(*Consequence.fields, sep='\t')
    with args.allele_matrix.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        header_tokens = next(line_token_gen)
        isolates = header_tokens[2:]
        # Iterate through each entry via AlleleRecords
        allele_record_gen = (AlleleRecord(line_tokens, isolates) for line_tokens in line_token_gen)
        for record in allele_record_gen:
            # Find features intervals that this record lands in
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
                    # Must complement codon if on reverse strand
                    if strand == -1:
                        allele = allele.translate(nucleotide_complement)
                    # Get consequence
                    sequence_allele, codon_allele = get_consequence(codon_sequence, allele, position_codon)
                    # Create new consequence object to hold output data for nicer printing imo
                    # Doubles to makes intergenic block much cleaner
                    consequence = Consequence.from_allele_record(record)
                    consequence.allele = allele
                    consequence.codon = codon
                    consequence.codon_allele = codon_allele
                    consequence.sequence = codon_sequence
                    consequence.sequence_allele = sequence_allele
                    consequence.change_type = 'ns' if codon != codon_allele else 's'
                    consequence.isolate_str = ','.join(allele_isolates)
                    [consequence.gene_id] = interval.feature.qualifiers['locus_tag']
                    [consequence.gene_product] = interval.feature.qualifiers['product']
                    # Print data
                    print(consequence)


def isolates_by_allele(record):
    ret = dict()
    for isolate, allele in zip(record.isolates, record.alleles):
        if allele == record.reference:
            continue
        if allele not in ret:
            ret[allele] = list()
        ret[allele].append(isolate)
    return ret


def get_consequence(codon_sequence, allele, position_codon):
    sequence_list = list(codon_sequence)
    sequence_list[position_codon-1] = allele
    sequence = Bio.Seq.Seq(''.join(sequence_list), Bio.Alphabet.IUPAC.unambiguous_dna)
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
