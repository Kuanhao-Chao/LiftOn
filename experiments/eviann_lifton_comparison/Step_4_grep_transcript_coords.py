from gffutils import FeatureDB
import pybedtools
import os, re
import gffutils
import json

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

def get_id_fraction(reference, target):
    matches = 0
    for i, letter in enumerate(reference):
        if letter == target[i]:
            matches += 1
    return matches, max(len(reference), len(target))

TARGET = "CHM13_MANE"
genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"
ref_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff"

output_dir = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/eviann_lifton_cmp/"

os.makedirs(output_dir, exist_ok=True)

eviann_protein_fa = f"/ccb/salz2/kh.chao/Lifton/data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna_AA.fa"
eviann_gff = f"/ccb/salz2/kh.chao/Lifton/data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.chr_fix.gff"

lifton_protein_fa = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_lifton_protein.fa"
lifton_protein_gff = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_lifton.gff3" 


print("protein_fa: ", protein_fa)
print("eviann_protein_fa: ", eviann_protein_fa)
print("lifton_protein_fa: ", lifton_protein_fa)
print("lifton_protein_gff: ", lifton_protein_gff)


##############################
# Creating database
##############################
dbfn_ref_gff = ref_gff + "_db"
dbfn_eviann_gff = eviann_gff + "_db"
dbfn_lifton_protein_gff = lifton_protein_gff + "_db"

# Create the GFF database
# db_MANE_ref_gff = gffutils.create_db(
#     MANE_ref_gff,
#     dbfn=dbfn_eviann_gff,
#     merge_strategy="merge",
#     force=True,  # Set this to True to overwrite the database if it already exists    
# )

# db_eviann_gff = gffutils.create_db(
#     eviann_gff,
#     dbfn=dbfn_eviann_gff,
#     merge_strategy="merge",
#     force=True,  # Set this to True to overwrite the database if it already exists    
# )

db_lifton_protein_gff = gffutils.create_db(
    lifton_protein_gff,
    dbfn=dbfn_lifton_protein_gff,
    merge_strategy="merge",
    force=True,  # Set this to True to overwrite the database if it already exists    
)

# # Load the GFF database
db_ref_gff = gffutils.FeatureDB(dbfn_ref_gff)
db_eviann_gff = gffutils.FeatureDB(dbfn_eviann_gff)
# db_lifton_protein_gff = gffutils.FeatureDB(dbfn_lifton_protein_gff)


frname_both = output_dir + "identities.txt"
frname_eviann_only = output_dir + "eviann_only.txt"
frname_lifton_only = output_dir + "lifton_only.txt"

fwname_both = output_dir + "identities_coords.txt"
fwname_eviann_only = output_dir + "eviann_only_coords.txt"
fwname_lifton_only = output_dir + "lifton_only_coords.txt"


# frname_miss_both = output_dir + "miss_both.txt"

fr_both = open(frname_both, "r")
fr_eviann_only = open(frname_eviann_only, "r")
fr_lifton_only = open(frname_lifton_only, "r")

frs = [fr_both, fr_eviann_only, fr_lifton_only]
fwfnames = [fwname_both, fwname_eviann_only, fwname_lifton_only]

gene_loci_dict_fn = "/ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/eviann_lifton_cmp/gene_loci.json"

with open(gene_loci_dict_fn, 'r') as openfile:
    gene_loci_dict = json.load(openfile)
print("gene_loci_dict: ", gene_loci_dict)


for idx in range(3):
    fr = frs[idx]
    fwfname = fwfnames[idx]

    with open(fwfname, "w") as fw:
        lines = fr.read().splitlines()
        for line in lines:
            # # entry = str(line)
            print(line)
            trans_id = line.split("\t")[0]
            # print(trans_id)
            
            if idx == 0:
                # eviann database and ID
                eles_ref_db = str(db_ref_gff[trans_id]).split("\t")
                # print("eles_MANE_db: ", eles_MANE_db)
                coord_MANE = f"\t{eles_ref_db[0]}:{eles_ref_db[3]}-{eles_ref_db[4]}"

                coords = ""                
                # eviann database and ID
                for eviann_id in gene_loci_dict[trans_id]:
                    eles = str(db_eviann_gff[eviann_id]).split("\t")
                    coords = coords + f"\t{eles[0]}:{eles[3]}-{eles[4]}"
                    # fw.write(coord)
                    # print(coord)

                # lifton database and ID
                eles = str(db_lifton_protein_gff[trans_id]).split("\t")
                coord = f"{eles[0]}:{eles[3]}-{eles[4]};"
                fw.write(f"{line}\t{coord}\t{coords}\n")

            elif idx == 1:
                coords = ""                
                # eviann database and ID
                for eviann_id in gene_loci_dict[trans_id]:
                    eles = str(db_eviann_gff[eviann_id]).split("\t")
                    coords = coords + f"\t{eles[0]}:{eles[3]}-{eles[4]}"
                    # fw.write(coord)
                    # print(coord)
                fw.write(coords)
                # print(coord)

            elif idx == 2:
                # lifton database and ID
                # print(db_eviann_gff[trans_id])
                coords = ""
                eles = str(db_lifton_protein_gff[trans_id]).split("\t")
                coords = coords + f"{eles[0]}:{eles[3]}-{eles[4]};"
                fw.write(f"{trans_id}\t{coords}\n")


frname_lifton_only = output_dir + "lifton_only.txt"
fr_lifton_only = open(frname_lifton_only, "r")
