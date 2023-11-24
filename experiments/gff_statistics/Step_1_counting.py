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

# TARGET = sys.argv[1]

# protein_fa = ""

# if TARGET == "human_to_chimp" or TARGET == "mouse_to_rat" or TARGET == "drosophila" or TARGET == "yeast" or TARGET == "arabadop" or TARGET == "bee" or TARGET == "mouse" or TARGET == "rice" or TARGET == "human_mane" or TARGET == "human_chess"  or TARGET == "human_refseq" or TARGET == "drosophila_erecta" or TARGET == "human_mane_to_mouse" or TARGET == "human_refseq_to_mouse" or TARGET == "human_to_chimp_test" or TARGET == "mouse_to_rat_test" or TARGET == "drosophila_test" or TARGET == "yeast_test" or TARGET == "arabadop_test" or TARGET == "bee_test" or TARGET == "mouse_test" or TARGET == "rice_test" or TARGET == "human_mane_test" or TARGET == "human_chess_test"  or TARGET == "human_refseq_test" or TARGET == "drosophila_erecta_test" or TARGET == "human_mane_to_mouse_test" or TARGET == "human_refseq_to_mouse_test" or TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1":
#     print("Running with ", TARGET)
#     genome = ""

#     if TARGET == "human_to_chimp": 
#         genome = "/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3-v1.1.fna"
    
#     elif TARGET == "mouse_to_rat": 
#         genome = "/ccb/salz2/jheinz3/shared/lifton/cross_species/mouse_to_rat/mRatBN7.2_genomic.fna"

#     elif TARGET == "yeast": 
#         genome = "/ccb/salz2/kh.chao/Lifton/data/yeast/S288C_reference_genome_R64-4-1_20230830/S288C_reference_sequence_R64-4-1_20230830.fsa"

#     elif TARGET == "arabadop": 
#         genome = "/ccb/salz2/jheinz3/shared/lifton/arabadop/Tanz-1.fna"

#     elif TARGET == "bee": 
#         genome = "/ccb/salz2/jheinz3/shared/lifton/bee/ASM1932182v1_genomic.fna"

#     elif TARGET == "mouse": 
#         genome = "/ccb/salz2/jheinz3/shared/lifton/mouse/mMusMuc1.1_genomic.fna"

#     elif TARGET == "rice": 
#         genome = "/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_genomic.fna"


#     elif TARGET == "GRCh38_RefSeq":
#         genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"

#     elif TARGET == "CHM13_MANE" or TARGET == "human_mane": 
#         genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"

#     elif TARGET == "CHM13_RefSeq":
#         genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"

#     elif TARGET == "Han1":
#         genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/Han1/v1.0/Assembly/Han1_v1.2.fasta"

#     elif TARGET == "Ash1":
#         genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/Ash1/v2.2/Assembly/Ash1_v2.2.fa"

#     elif TARGET == "PR1":
#         genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/PR1/v3.0/Assembly/PR1.fa"

#     elif TARGET == "Mus_musculus_MANE":
#         genome = "/ccb/salz3/kh.chao/ref_genome/mus_musculus/NCBI_Refseq/GCF_000001635.27_GRCm39_genomic.fna"
# else:
#     sys.exit(-1)



# def count_features(gff_file):
#     gene_count = 0
#     transcript_count = 0
#     protein_coding_gene_count = 0
#     protein_coding_transcript_count = 0

#     with open(gff_file, 'r') as f:
#         for line in f:
#             if line.startswith('#'):
#                 continue
#             fields = line.strip().split('\t')
#             if len(fields) < 3:
#                 continue

#             feature_type = fields[2]
#             attributes = fields[-1]

#             if feature_type == 'gene':
#                 gene_count += 1
#                 if 'protein_coding' in attributes:
#                     protein_coding_gene_count += 1

#             if feature_type == 'mRNA' or feature_type == 'transcript':
#                 transcript_count += 1
#                 if feature_type == 'mRNA' or 'protein_coding' in attributes:
#                     protein_coding_transcript_count += 1

#     return gene_count, transcript_count, protein_coding_gene_count, protein_coding_transcript_count


gene_count = 0
protein_coding_gene_count = 0

trans_count = 0
protein_coding_trans_count = 0




def count_features(ref_db):
    global gene_count
    global trans_count
    global protein_coding_gene_count
    global protein_coding_trans_count
    
    features = ["gene"]
    for feature in features:
        for locus in ref_db.features_of_type(feature):
            gene_count += 1

            if locus.attributes["gene_biotype"][0] == "protein_coding":
                protein_coding_gene_count += 1
            __inner_count_feature(ref_db, locus)


def __inner_count_feature(ref_db, feature):

    global gene_count
    global trans_count
    global protein_coding_gene_count
    global protein_coding_trans_count

    # If exon is the first level children
    children_exons = list(ref_db.children(feature, featuretype='exon', level=1))
    children_CDSs = list(ref_db.children(feature, featuretype=('start_codon', 'CDS', 'stop_codon'), level=1))

    if len(children_exons) > 0 or len(children_CDSs) > 0:
        # print(f"Parent: {feature.id};  exon: {len(children_exons)}; CDS: {len(children_CDSs)}")
        if len(children_exons) > 0:
            trans_count += 1

        if feature.featuretype == "mRNA":
            protein_coding_trans_count += 1


    else:
        for child in ref_db.children(feature, level=1, order_by='start'):
            __inner_count_feature(ref_db, child)





if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Count features in a GFF file')
    parser.add_argument('gff_file', help='Input GFF file')
    args = parser.parse_args()

    db = gffutils.FeatureDB(args.gff_file, keep_order=True)    
    count_features(db)

    print(f'Total number of genes: {gene_count}')
    print(f'Total number of transcripts: {trans_count}')
    print(f'Total number of protein-coding genes: {protein_coding_gene_count}')
    print(f'Total number of protein-coding transcripts: {protein_coding_trans_count}')