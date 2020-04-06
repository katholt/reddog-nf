#!/usr/bin/env python3
# TODO: deal with '-' and 'N'
# TODO: clean up printing approach
# TODO: clean up consequence section in general please


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

    # Get consequences for each allele
    with args.allele_matrix.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        header_tokens = next(line_token_gen)
        isolates = header_tokens[2:]
        for position_str, reference_allele, *isolate_alleles in line_token_gen:
            # Get features
            position_snp = int(position_str)
            features = interval_tree.search_position(position_snp)
            # Process intergenic
            if not features:
                continue
                isolates_by_allele = sort_isolates_by_allele(isolates, isolate_alleles, reference_allele)
                empty_fields = [''] * 6  # TODO: please fix this
                for allele, allele_isolates in isolates_by_allele.items():
                    isolate_str = ','.join(allele_isolates)
                    print(position_snp, reference_allele, allele, 'intergenic', *empty_fields, isolate_str, sep='\t')
                continue
            # Process genic
            for interval in interval_tree.search_position(position_snp):
                # Get strand
                feature = interval.feature
                if feature.strand == 1:
                    continue
                    strand = 'forward'
                elif feature.strand == -1:
                    strand = 'reverse'
                else:
                    print(f'error: got bad strand from {feature}', file=sys.stderr)
                    sys.exit(1)
                process_allele(position_snp, feature, interval, strand, isolates, isolate_alleles,
                        reference_allele, rep_gbk)


def process_allele(position_snp, feature, interval, strand, isolates, isolate_alleles, reference_allele, rep_gbk):
    # Get SNP location within the gene and codon
    if strand == 'forward':
        position_gene = position_snp - interval.start + 1
    elif strand == 'reverse':
        position_gene = interval.end - position_snp + 1
    position_codon = (position_gene - 1) % 3 + 1
    # Get the nth codon which contains the SNP
    codon_number = (position_gene - 1) // 3 + 1
    codon_start = (codon_number - 1) * 3
    # Collect gene and codon sequence
    gene_sequence = feature.extract(rep_gbk.seq)
    codon_sequence = gene_sequence[codon_start:codon_start+3]
    codon = codon_sequence.translate()
    # Grab some other information for output
    [gene_id] = feature.qualifiers['locus_tag']
    [gene_product] = feature.qualifiers['product']
    # Sort isolates by allele
    # TODO: check if there's a builtin or stdlib function for this now
    isolates_by_allele = sort_isolates_by_allele(isolates, isolate_alleles, reference_allele)
    # Get consequences for each allele
    for allele, allele_isolates in isolates_by_allele.items():
        if strand == 'forward':
            allele_cons = allele
        elif strand == 'reverse':
            allele_cons = allele.translate(nucleotide_complement)
        sequence, codon_allele = get_consequence(codon_sequence, allele_cons, position_codon)
        change_type = 'ns' if codon != codon_allele else 's'
        isolate_str = ','.join(allele_isolates)
        print(position_snp, reference_allele, allele, change_type, gene_id, codon_sequence,
                sequence, codon, codon_allele, gene_product, isolate_str, sep='\t')


def sort_isolates_by_allele(isolates, isolate_alleles, reference_allele):
    ret = dict()
    for isolate, allele in zip(isolates, isolate_alleles):
        if allele == reference_allele:
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
