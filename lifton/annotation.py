import gffutils
import sys
from collections import defaultdict


class Annotation():

    def __init__(self, file_name, infer_genes, infer_transcripts, merge_strategy, id_spec, force, verbose):
        self.file_name = file_name
        self.infer_genes = infer_genes
        self.infer_transcripts = infer_transcripts
        self.merge_strategy = merge_strategy
        self.id_spec = id_spec
        self.force = force
        self.verbose = verbose
        self.get_db_connnection()
        
    def get_db_connnection(self):
        try:
            feature_db = gffutils.FeatureDB(self.file_name)
        except:
            feature_db = self.build_database()
        self.db_connection = feature_db


    def build_database(self):
        if self.infer_genes:
            disable_genes = False
        else:
            disable_genes = True
        if self.infer_transcripts:
            disable_transcripts = False
        else:
            disable_transcripts = True
        try:
            transform_func = self.get_transform_func()
            feature_db = gffutils.create_db(self.file_name, self.file_name + "_db", 
                                        merge_strategy=self.merge_strategy, 
                                        id_spec=self.id_spec,
                                        force=self.force,
                                        verbose=self.verbose, disable_infer_transcripts=disable_transcripts,
                                            disable_infer_genes=disable_genes, transform=transform_func)
        except Exception as e:
            print("gffutils database build failed with", e)
            feature_db = self.build_database_again()
        return feature_db


    def build_database_again(self):
        if self.infer_genes:
            disable_genes = False
        else:
            disable_genes = True
        if self.infer_transcripts:
            disable_transcripts = False
        else:
            disable_transcripts = True
        try:
            transform_func = self.get_transform_func()
            feature_db = gffutils.create_db(self.file_name, self.file_name + "_db", 
                                        merge_strategy="warning",
                                        id_spec=self.id_spec,
                                        force=True,
                                        verbose=True, disable_infer_transcripts=disable_transcripts,
                                            disable_infer_genes=disable_genes, transform=transform_func)
        except Exception as e:
            print("gffutils database build failed with", e)
            sys.exit()
        return feature_db


    def get_transform_func(self):
        if self.infer_genes is False:
            return None
        else:
            return transform_func



    def get_protein_coding_features(self, feature_types):
        protein_coding_genes = []
        for f_type in feature_types:
            for feature in self.db_connection.features_of_type(featuretype=f_type):
                if self.is_highest_parent(feature):
                    CDS_list = [child for child in self.db_connection.children(feature) if child.featuretype == 'CDS']
                    if len(CDS_list) > 0:
                        protein_coding_genes.append(feature.id)
        return protein_coding_genes


    def get_noncoding_features(self, feature_types):
        noncoding_genes = []
        for f_type in feature_types:
            for feature in self.db_connection.features_of_type(featuretype=f_type):
                if self.is_highest_parent(feature):
                    CDS_list = [child for child in self.db_connection.children(feature) if child.featuretype == 'CDS']
                    if len(CDS_list)==  0:
                        noncoding_genes.append(feature.id)
        return noncoding_genes


    def get_novel_protein_coding_features(self, ref_genes, feature_types):
        novel_protein_coding = []
        all_protein_coding = self.get_protein_coding_features(feature_types)
        for feature in all_protein_coding:
            if feature not in ref_genes:
                novel_protein_coding.append(feature)
        return novel_protein_coding


    def get_novel_noncoding_features(self,ref_genes, feature_types):
        novel_noncoding = []
        all_noncoding = self.get_noncoding_features(feature_types)
        for feature in all_noncoding:
            if feature not in ref_genes:
                novel_noncoding.append(feature)
        return novel_noncoding


    def get_all_parent_feature_ids(self, feature_types):
        feature_ids = []
        for f_type in feature_types:
            feature_ids += [feature.id for feature in self.db_connection.features_of_type(
                featuretype=f_type) if self.is_highest_parent(feature)]
        return feature_ids


    def make_parent_to_child_dict(self, protein_coding, gene):
        child_dict=defaultdict(list)
        child_list = self.db_connection.children(gene)
        child_count = 0
        for child in child_list:
            child_count += 1
            if (protein_coding and child.featuretype == "CDS") or (protein_coding is False and self.is_lowest_child(
                    child)):
                parent = list(self.db_connection.parents(child, level=1))[0]
                child_dict[parent.id].append(child)
        if child_count == 0:
            child_dict[gene] = [self.db_connection[gene]]
        return child_dict


    def get_features_of_type(self,feature_types):
        features_of_type = []
        for f_type in feature_types:
            features_of_type += list(self.db_connection.features_of_type(f_type))
        return features_of_type


    def get_feature_dict(self, feature_types):
        id_to_feature = {}
        features = self.get_features_of_type(feature_types)
        for feature in features:
            id_to_feature[feature.id] = feature
        return id_to_feature


    def get_source_name(self, feature_name):
        feature = self.db_connection[feature_name]
        if 'extra_copy_number' in feature.attributes:
            if feature.attributes['extra_copy_number'][0] == feature.id.split("_")[-1]:
                return "_".join(feature.id.split("_")[:-1])
        return feature.id


    def is_lowest_child(self, feature_name):
        return len(list(self.db_connection.children(feature_name))) == 0


    def is_highest_parent(self, feature_name):
        return len(list(self.db_connection.parents(feature_name))) == 0


    def get_paralog_name(self, feature_name):
        gene_attributes = self.db_connection[feature_name].attributes
        if "extra_copy_number" in gene_attributes:
            paralog_name = "_".join(gene_attributes["ID"][0].split("_")[:-1])
            return paralog_name
        return ""


    def get_num_levels(self, feature_name):
        level = 1
        new_children = list(self.db_connection.children(feature_name, level=level))
        while len(new_children) !=0:
            level += 1
            new_children = list(self.db_connection.children(feature_name, level=level))
        return level



def transform_func(x):
    if 'transcript_id' in x.attributes:
        x.attributes['transcript_id'][0] += '_transcript'
    return x


def get_perc_id(feature):
    if 'sequence_ID' in feature.attributes:
        return float(feature.attributes['sequence_ID'][0])
    return 0.0


def merge_children_intervals(children):
    if len(children) == 0:
        return []
    intervals = [[child.start, child.end] for child in children]
    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged