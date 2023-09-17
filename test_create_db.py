import gffutils
from itertools import chain


gff_file_l = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/CHM13_MANE/CHM13_MANE.sort.sorted.gff3"
# feature_db_l = gffutils.create_db(gff_file_l, gff_file_l + "_db", merge_strategy="create_unique", force=True,
#                                 disable_infer_transcripts=False,
#                                 disable_infer_genes=True, verbose=True, keep_order=True, 
#                                 sort_attribute_values=True)

feature_db_l = gffutils.FeatureDB(gff_file_l + "_db")

genes_l = list(feature_db_l.features_of_type('gene', limit=("chr1",0,500000000), strand="+"))
genes_l_len = len(genes_l)

print("genes_l: ", genes_l_len)
# for l_idx in range(len(genes_l)):
#     print(genes_l[l_idx])

gff_file_m = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/CHM13_MANE/CHM13_MANE_miniprot.fix.sorted.gff"
# feature_db_m = gffutils.create_db(gff_file_m, gff_file_m + "_db", merge_strategy="create_unique", force=True,
#                                 disable_infer_transcripts=False,
#                                 disable_infer_genes=True, verbose=True, keep_order=True, 
#                                 sort_attribute_values=True)

feature_db_m = gffutils.FeatureDB(gff_file_m + "_db")

trans_m = list(feature_db_m.features_of_type('mRNA', limit=("chr1",0,500000000), strand="+"))
trans_m_len = len(trans_m)

print("trans_m: ", trans_m_len)

l_idx = 0
m_idx = 0

while l_idx < genes_l_len or m_idx < trans_m_len:
    if l_idx == genes_l_len:
        m_idx += 1
    elif m_idx == trans_m_len:
        l_idx += 1

    liftoff_gene = genes_l[l_idx]
    miniprot_trans = trans_m[m_idx]

    print("Liftoff : ", liftoff_gene.start, liftoff_gene.end)
    print("miniprot: ", miniprot_trans.start, miniprot_trans.end)

    if liftoff_gene.end >= miniprot_trans.end:
        m_idx += 1
    elif liftoff_gene.end < miniprot_trans.end:
        l_idx += 1


# def overlapping(s_1, e_1, s_2, e_2):

#     if e_1 <= s_2:
#         pass
#     elif 


# for l_i, l_feature in enumerate(feature_db_l.all_features()):
#     print(l_i)
#     print(l_feature)


# c = chain(feature_db_m.all_features(), feature_db_l.all_features())  # Chain the generator and a string

# feature_db_m.all_features()

gff_file_merged = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/CHM13_MANE/merged.gff"
# # feature_db_merged = gffutils.create_db(c, gff_file_merged + "_db", merge_strategy="create_unique", force=True,
# #                                 disable_infer_transcripts=False,
# #                                 disable_infer_genes=True, verbose=True, keep_order=True, 
# #                                 sort_attribute_values=True)

# feature_db_merged = gffutils.FeatureDB(gff_file_merged + "_db")

# Liftoff_count = 0
# miniprot_count = 0

# for i, feature in enumerate(feature_db_merged.all_features()):
#     print(i)
#     print(feature)
    
# # for feature in feature_db_merged.all_features():#strand="+", order_by="seqid"):
#     if feature.source == "Liftoff":
#         Liftoff_count += 1
#     elif feature.source == "miniprot":
#         miniprot_count += 1

# print("Liftoff_count: ", Liftoff_count)
        
# print("miniprot_count: ", miniprot_count)

