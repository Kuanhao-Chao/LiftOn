from Bio.Seq import Seq
from collections import defaultdict
import numpy as np
from liftofftools.annotation import merge_children_intervals


class FeatureDict(defaultdict):
    def __init__(self, annotation, feature_type):
        super().__init__(list)
        self.populate_feature_dict(annotation, feature_type)


    def populate_feature_dict(self, annotation, feature_types):
        db_connection = annotation.db_connection
        for feature in db_connection.all_features():
            if feature.featuretype in feature_types and annotation.is_lowest_child(feature):
                parents = list(db_connection.parents(feature, level=1))
                if len(parents) > 0:
                    self[parents[0].id].append(feature)
                else:
                    self[feature.id].append(feature)


class SequenceDict(dict):
    def __init__(self,annotation, fasta, feature_type, is_protein):
        super().__init__()
        self.annotation = annotation
        self.feature_type = feature_type
        self.feature_dict = FeatureDict(annotation, feature_type)
        self.populate_sequence_dict(fasta, feature_type)
        self.is_protein = is_protein


    def populate_sequence_dict(self, fasta, feature_types):
        for parent, feature_group in self.feature_dict.items():
            feature_types = [feature.featuretype for feature in feature_group]
            if 'CDS' in feature_types and 'exon' not in feature_types:
                sequence = get_protein_sequence(fasta, feature_group)
            elif 'CDS' in feature_types and 'exon' in feature_types:
                CDS_only_group = [feature for feature in feature_group if feature.featuretype == 'CDS']
                sequence = get_dna_sequence(fasta, CDS_only_group)
            else:
                sequence = get_dna_sequence(fasta, feature_group)
            self[parent] = sequence


    def get_longest_isoform_dict(self, gene_list):
        longest_isoform_dict = {}
        for gene in gene_list:
            parent_to_child = self.annotation.make_parent_to_child_dict(self.is_protein, gene)
            if len(parent_to_child) > 0:
                longest_isoform = get_longest_isoform(parent_to_child)
                longest_isoform_dict[gene] = self[longest_isoform]
        return longest_isoform_dict


def get_longest_isoform(child_dict):
    max_length = 0
    longest_tran = ''
    for parent, child_group in child_dict.items():
        length = np.sum([child.end - child.start for child in child_group])
        if length > max_length:
            max_length = length
            longest_tran = parent
    return longest_tran


def get_dna_sequence(fasta, features):
    chrom = features[0].seqid
    strand = features[0].strand
    merged_features = merge_children_intervals(features)
    sequence =''
    for start, end in merged_features:
        sequence += str(fasta[chrom][start -1: end])
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    sequence += 'N' * get_padding_length(len(sequence))
    return sequence

def get_padding_length(sequence_length):
    return (3-sequence_length%3)%3

def get_protein_sequence(fasta, features):
    dna = get_dna_sequence(fasta, features)
    return str(Seq(dna).translate())