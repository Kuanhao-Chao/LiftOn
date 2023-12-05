from gffutils import FeatureDB
import pybedtools
import os, re
import gffutils

import sys
from gffutils import FeatureDB
from Bio import SeqIO
import parasail

def read_fasta(file_path):
    sequences = {}
    current_sequence = ""
    current_header = ""

    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if current_header != "":
                    sequences[current_header] = current_sequence
                current_header = line[1:]
                current_sequence = ""
            else:
                current_sequence += line

        if current_header != "":
            sequences[current_header] = current_sequence

    return sequences

TARGET = "human_refseq_test"

genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"
protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/protein.fasta"


output_dir = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/ref_chm13_cmp/"

os.makedirs(output_dir, exist_ok=True)

ref_chm13_protein_fa = f"/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0_RefSeq_Liftoff_v5.1_protein.fa"
ref_chm13_gff = f"/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0_RefSeq_Liftoff_v5.1.gff3"
ref_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.gff"

lifton_protein_fa = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/lifton_protein.fa"
lifton_protein_gff = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/lifton.gff3" 


print("protein_fa: ", protein_fa)
print("ref_chm13_protein_fa: ", ref_chm13_protein_fa)
print("lifton_protein_fa: ", lifton_protein_fa)
print("lifton_protein_gff: ", lifton_protein_gff)


# ##############################
# # Creating database
# ##############################
# dbfn_ref_gff = ref_gff + "_db"
# dbfn_ref_chm13_gff = ref_chm13_gff + "_db"
# dbfn_lifton_protein_gff = lifton_protein_gff + "_db"

# # Create the GFF database
# db_ref_chm13_gff = gffutils.create_db(ref_chm13_gff, dbfn=dbfn_ref_chm13_gff,
#                             merge_strategy="create_unique", 
#                                 # merge_strategy="create_unique", 
#                             # id_spec='ID',
#                             force=True,
#                             verbose=True, disable_infer_transcripts=True,
#                                 disable_infer_genes=True)

# db_lifton_protein_gff = gffutils.create_db(lifton_protein_gff, dbfn=dbfn_lifton_protein_gff,
#                             merge_strategy="create_unique", 
#                                 # merge_strategy="create_unique", 
#                             # id_spec='ID',
#                             force=True,
#                             verbose=True, disable_infer_transcripts=True,
#                                 disable_infer_genes=True)

# # # Load the GFF database
# db_ref_gff = gffutils.FeatureDB(dbfn_ref_gff)
# # db_ref_chm13_gff = gffutils.FeatureDB(dbfn_ref_chm13_gff)
# # db_lifton_protein_gff = gffutils.FeatureDB(dbfn_lifton_protein_gff)



# frname_both = output_dir + "identities.txt"
# frname_ref_chm13_only = output_dir + "ref_chm13_only.txt"
# frname_lifton_only = output_dir + "lifton_only.txt"

# fwname_both = output_dir + "identities_coords.txt"
# fwname_ref_chm13_only = output_dir + "ref_chm13_only_coords.txt"
# fwname_lifton_only = output_dir + "lifton_only_coords.txt"


# fr_both = open(frname_both, "r")
# fr_ref_chm13_only = open(frname_ref_chm13_only, "r")
# fr_lifton_only = open(frname_lifton_only, "r")

# frs = [fr_both, fr_ref_chm13_only, fr_lifton_only]
# fwfnames = [fwname_both, fwname_ref_chm13_only, fwname_lifton_only]

# for idx in range(3):
#     fr = frs[idx]
#     fwfname = fwfnames[idx]

#     with open(fwfname, "w") as fw:
#         lines = fr.read().splitlines()
#         for line in lines:
#             # # entry = str(line)
#             # print(line.start)
#             trans_id = line.split("\t")[0]
#             # print(trans_id)
            
#             if idx == 0:
#                 # Liftoff database and ID
#                 eles_ref_db = str(db_ref_gff[trans_id]).split("\t")
#                 # print("eles_MANE_db: ", eles_MANE_db)
#                 coord_MANE = f"\t{eles_ref_db[0]}:{eles_ref_db[3]}-{eles_ref_db[4]}"

#                 eles = str(db_ref_chm13_gff[trans_id]).split("\t")
#                 coord = f"\t{eles[0]}:{eles[3]}-{eles[4]}"
#                 # fw.write(coord)
#                 # print(coord)

#                 # lifton database and ID
#                 coords = ""                
#                 eles = str(db_lifton_protein_gff[trans_id]).split("\t")
#                 coords = coords + f"{eles[0]}:{eles[3]}-{eles[4]};"
#                 fw.write(f"{line}\t{coord}\t{coords}\n")
            
#             elif idx == 1:
#                 # Liftoff database and ID
#                 eles = str(db_ref_chm13_gff[trans_id]).split("\t")
#                 coord = f"{trans_id}\t{eles[0]}:{eles[3]}-{eles[4]}\n"
#                 fw.write(coord)
#                 # print(coord)

#             elif idx == 2:
#                 # lifton database and ID
#                 # print(db_ref_chm13_gff[trans_id])
#                 coords = ""
#                 eles = str(db_lifton_protein_gff[trans_id]).split("\t")
#                 coords = coords + f"{eles[0]}:{eles[3]}-{eles[4]};"
#                 fw.write(f"{trans_id}\t{coords}\n")


# frname_lifton_only = output_dir + "lifton_only.txt"
# fr_lifton_only = open(frname_lifton_only, "r")
