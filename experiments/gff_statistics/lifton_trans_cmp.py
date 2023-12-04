####################
# Goal of this script is to count the number of 
# 1. genes, 
# 2. transcripts, 
# 3. protein-coding genes,
# 4. protein-coding transcripts
####################

import argparse
import gffutils
import sys

class Lifton_feature:
    def __init__(self, id):
        self.id = id
        self.copy_num = 0
        self.is_protein_coding = False
        self.children = set()

def get_ref_liffover_features(features, ref_db):
    ref_features_dict = {}
    ref_features_reverse_dict = {}
    new_gene_feature = Lifton_feature("Lifton-gene")
    ref_features_dict["LiftOn-gene"] = new_gene_feature

    for f_itr in features:
        for locus in ref_db.db_connection.features_of_type(f_itr):

            CDS_children = list(ref_db.db_connection.children(locus, featuretype='CDS'))

            feature = Lifton_feature(locus.id)
            if len(CDS_children) > 0:
                # This is the protien-coding gene
                feature.is_protein_coding = True
                
            exon_children = list(ref_db.db_connection.children(locus, featuretype='exon', level=1, order_by='start'))

            if len(exon_children) > 0:
                process_ref_liffover_features(locus, ref_db, None)
            else:
                transcripts = ref_db.db_connection.children(locus, level=1)
                for transcript in list(transcripts):
                    # print(transcript)
                    process_ref_liffover_features(transcript, ref_db, feature)
                    ref_features_reverse_dict[transcript.id] = locus.id
            ref_features_dict[locus.id] = feature
    return ref_features_dict, ref_features_reverse_dict


def process_ref_liffover_features(locus, ref_db, feature):
    if feature != None:
        feature.children.add(locus.id)

    
def count_features(ref_db):
    ref_features_dict, ref_features_reverse_dict = get_ref_liffover_features(["gene"], ref_db)


def build_database(infer_genes = True):
    disable_genes = False
    try:
        # transform_func = self.get_transform_func()
        feature_db = gffutils.create_db(args.gff_file, args.gff_file + "_db", 
                                    merge_strategy="create_unique", 
                                        # merge_strategy="create_unique", 
                                    # id_spec='ID',
                                    force=True,
                                    verbose=True, disable_infer_transcripts=True, disable_infer_genes=disable_genes)
        # , transform=transform_func)
        
                                    # id_spec={"gene": ['ID', 'Name'], "mRNA": ['ID', 'Name'], "transcript": ['ID', 'Name'], "lnc_RNA": ['ID', 'Name'], "nc_RNA": ['ID', 'Name']},

    except Exception as e:
        print("gffutils database build failed with", e)
        sys.exit()
    return feature_db



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Count features in a GFF file')
    parser.add_argument('gff_file', help='Input GFF file')
    args = parser.parse_args()

    # db = gffutils.FeatureDB(args.gff_file, keep_order=True)    



    try:
        feature_db = gffutils.FeatureDB(args.gff_file, keep_order=True)
    except:
        feature_db = build_database(args.gff_file)
    
    # # print("feature_db: ", feature_db)
    # # feature_db.execute('ANALYZE features')
    # self.db_connection = feature_db


    count_features(feature_db)

    print(f'Total number of genes: {gene_count}')
    print(f'Total number of transcripts: {trans_count}')
    print(f'Total number of protein-coding genes: {protein_coding_gene_count}')
    print(f'Total number of protein-coding transcripts: {protein_coding_trans_count}')