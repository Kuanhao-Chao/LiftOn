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

    for f_itr in features:
        for locus in ref_db.features_of_type(f_itr):                
            exon_children = list(ref_db.children(locus, featuretype='exon', level=1, order_by='start'))

            if len(exon_children) > 0:
                process_ref_liffover_features(locus, ref_db, None)
            else:
                transcripts = ref_db.children(locus, level=1)
                for transcript in list(transcripts):

                    CDS_children = list(ref_db.children(transcript, featuretype='CDS'))
                    feature = Lifton_feature(transcript.id)
                    if len(CDS_children) > 0:
                        # This is the protien-coding gene
                        feature.is_protein_coding = True

                    process_ref_liffover_features(transcript, ref_db, feature)
                    # ref_features_reverse_dict[transcript.id] = locus.id
                    ref_features_dict[transcript.id] = feature
    return ref_features_dict, ref_features_reverse_dict


def process_ref_liffover_features(locus, ref_db, feature):
    if feature != None:
        feature.children.add(locus.id)


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


def get_ID_base(id):
    # # Regular expression pattern to match the desired substrings
    splits = id.split("_")
    try:
        int(splits[-1])
        id_base = "_".join(splits[:-1])
    except:
        id_base = id
    return id_base

def get_ID(tgt_key, ref_features_dict):
    if tgt_key in ref_features_dict.keys():
        return tgt_key
    else:
        tgt_key_base = get_ID_base(tgt_key)
        return tgt_key_base
    
        # if tgt_key_base in ref_features_dict.keys():
        #     return tgt_key_base
        # else:
        #     return tgt_key



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Count features in a GFF file')
    parser.add_argument('gff_file', help='Input GFF file')
    parser.add_argument('ref_gff_file', help='Input GFF file')


    args = parser.parse_args()

    # db = gffutils.FeatureDB(args.gff_file, keep_order=True)    



    try:
        feature_db = gffutils.FeatureDB(args.gff_file, keep_order=True)
    except:
        feature_db = build_database(args.gff_file)
    
    # # print("feature_db: ", feature_db)
    # # feature_db.execute('ANALYZE features')
    # self.db_connection = feature_db

    # ref_annotation = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/rRNA_removed/NCBI_RefSeq_no_rRNA.gff_db"

    # ref_annotation = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff_db"

    ref_feature_db = gffutils.FeatureDB(args.ref_gff_file, keep_order=True)

    ref_features_dict, ref_features_reverse_dict = get_ref_liffover_features(["gene"], ref_feature_db)

    tgt_features_dict, tgt_features_reverse_dict = get_ref_liffover_features(["gene"], feature_db)

    print("ref_features_dict         : ", len(ref_features_dict))
    print("tgt_features_dict         : ", len(tgt_features_dict))



    ref_trans_count = len(ref_features_dict)
    tgt_trans_count = len(tgt_features_dict)

    lifted_trans_coding_count = 0
    lifted_trans_noncoding_count = 0

    missed_trans_coding_count = 0
    missed_trans_noncoding_count = 0


    for ref_key, ref_val in ref_features_dict.items():
        if ref_key in tgt_features_dict:
            tgt_val = tgt_features_dict[ref_key]
            if tgt_val.is_protein_coding:
                lifted_trans_coding_count += 1
            else:
                lifted_trans_noncoding_count += 1
            #     print("ref_key: ", ref_key)
            #     print("ref_val: ", ref_val.is_protein_coding)
            #     print("tgt_val: ", tgt_val.is_protein_coding)
            #     print("")
        else:
            ref_val = ref_features_dict[ref_key]
            if ref_val.is_protein_coding:
                missed_trans_coding_count += 1
            else:
                missed_trans_noncoding_count += 1



    tgt_trans_coding_single_cp = 0
    tgt_trans_coding_multi_cp = 0
    tgt_trans_coding_multi_cp_count = 0
    tgt_coding_lifted = {}
    

    tgt_trans_noncoding_single_cp = 0
    tgt_trans_noncoding_multi_cp = 0
    tgt_trans_noncoding_multi_cp_count = 0
    tgt_noncoding_lifted = {}

    for tgt_key, tgt_val in tgt_features_dict.items():
        tgt_key = get_ID(tgt_key, ref_features_dict)
        if tgt_key in ref_features_dict.keys():
            if tgt_val.is_protein_coding:
                if tgt_key in tgt_coding_lifted.keys():
                    tgt_coding_lifted[tgt_key] += 1
                else:
                    tgt_coding_lifted[tgt_key] = 1
            else:
                if tgt_key in tgt_noncoding_lifted.keys():
                    tgt_noncoding_lifted[tgt_key] += 1
                else:
                    tgt_noncoding_lifted[tgt_key] = 1
        else:
            pass

        
        # if tgt_key in ref_features_dict.keys():
        #     return tgt_key
        # else:
        #     tgt_key_base = get_ID_base(tgt_key)
        #     if tgt_key_base in ref_features_dict.keys():
        #         return tgt_key_base
        #     else:
        #         return tgt_key



        #     tgt_coding_lifted[tgt_key] = 
        #     if tgt_val.is_protein_coding:

        # else:
        #     tgt_key = get_ID_base(tgt_key)
        #     if tgt_key in ref_features_dict.keys():
        #         tgt_coding_lifted_set.add(tgt_key)


        # if tgt_key_base in ref_features_dict:
        #     # ref_val = ref_features_dict[tgt_key]
            
        #     if tgt_val.is_protein_coding:
        #         lifted_trans_coding_count += 1
        #     else:
        #         lifted_trans_noncoding_count += 1

        # else:
        #     if tgt_val.is_protein_coding:
        #         missed_trans_coding_count += 1
        #     else:
        #         missed_trans_noncoding_count += 1


    print(f'Reference coding count    : {lifted_trans_coding_count + missed_trans_coding_count}')
    # print(f'missed_trans_coding_count    : {missed_trans_coding_count}')
    
    print(f'tgt_coding_lifted            : {len(tgt_coding_lifted)}')

    tgt_coding_lifted_single= 0
    tgt_coding_lifted_extra= 0
    tgt_coding_lifted_extra_sum = 0
    for key, val in tgt_coding_lifted.items():
        # tgt_coding_lifted_sum += val
        if val == 1:
            tgt_coding_lifted_single += 1
        elif val > 1:
            tgt_coding_lifted_extra += 1
            tgt_coding_lifted_extra_sum += val

    print(f'\ttgt_coding_lifted single     : {tgt_coding_lifted_single}')
    print(f'\ttgt_coding_lifted extra      : {tgt_coding_lifted_extra}')
    print(f'\ttgt_coding_lifted extra sum  : {tgt_coding_lifted_extra_sum}')
    print(f'\ttgt_coding_lifted total      : {tgt_coding_lifted_single + tgt_coding_lifted_extra_sum}')

    # print(f'\ttgt_trans_coding_lost         : {tgt_trans_coding_lost}')


    # print(f'lifted_trans_noncoding_count : {lifted_trans_noncoding_count}')
    # print(f'missed_trans_noncoding_count : {missed_trans_noncoding_count}')
    print("\n\n")

    print(f'Reference noncoding count    : {lifted_trans_noncoding_count + missed_trans_noncoding_count}')
    print(f'tgt_noncoding_lifted         : {len(tgt_noncoding_lifted)}')

    tgt_noncoding_lifted_single= 0
    tgt_noncoding_lifted_extra= 0
    tgt_noncoding_lifted_extra_sum = 0
    for key, val in tgt_noncoding_lifted.items():
        # tgt_coding_lifted_sum += val
        if val == 1:
            tgt_noncoding_lifted_single += 1
        elif val > 1:
            tgt_noncoding_lifted_extra += 1
            tgt_noncoding_lifted_extra_sum += val

    print(f'\ttgt_noncoding_lifted single     : {tgt_noncoding_lifted_single}')
    print(f'\ttgt_noncoding_lifted extra      : {tgt_noncoding_lifted_extra}')
    print(f'\ttgt_noncoding_lifted extra sum  : {tgt_noncoding_lifted_extra_sum}')
    print(f'\ttgt_noncoding_lifted total      : {tgt_noncoding_lifted_single + tgt_noncoding_lifted_extra_sum}')

    print("\n\n")

    print(f'Reference total count    : {lifted_trans_coding_count + missed_trans_coding_count + lifted_trans_noncoding_count + missed_trans_noncoding_count}')
    print(f'Target total count    : {tgt_coding_lifted_single + tgt_coding_lifted_extra_sum + tgt_noncoding_lifted_single + tgt_noncoding_lifted_extra_sum}')
