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

TARGET = sys.argv[1]
protein_fa = ""
if TARGET == "CHM13_MANE" or TARGET == "CHM13_RefSeq" or TARGET == "GRCh38_RefSeq" or TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1" or TARGET == "Mus_musculus_MANE":
    print("Running with ", TARGET)
    genome = ""
    if TARGET == "GRCh38_RefSeq":
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/protein.fasta"
        ref_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.gff"

    elif TARGET == "CHM13_MANE": 
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"
        ref_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff"

    elif TARGET == "CHM13_RefSeq":
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/protein.fasta"
        ref_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.gff"

    elif TARGET == "Han1":
        genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/Han1/v1.0/Assembly/Han1_v1.2.fasta"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"
        ref_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff"

    elif TARGET == "Ash1":
        genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/Ash1/v2.2/Assembly/Ash1_v2.2.fa"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"
        ref_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff"

    elif TARGET == "PR1":
        genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/PR1/v3.0/Assembly/PR1.fa"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"
        ref_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff"

    elif TARGET == "Mus_musculus_MANE":
        genome = "/ccb/salz3/kh.chao/ref_genome/mus_musculus/NCBI_Refseq/GCF_000001635.27_GRCm39_genomic.fna"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"
        ref_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff"

else:
    sys.exit(-1)


output_dir = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/liftoff_lifton_cmp/"

os.makedirs(output_dir, exist_ok=True)

liftoff_protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/{TARGET}/{TARGET}_AA.fa"
liftoff_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/{TARGET}/{TARGET}_liftoff.sort.gff3"

lifton_protein_fa = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_lifton_protein.fa"
lifton_protein_gff = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_lifton.gff3" 


print("protein_fa: ", protein_fa)
print("liftoff_protein_fa: ", liftoff_protein_fa)
print("lifton_protein_fa: ", lifton_protein_fa)
print("lifton_protein_gff: ", lifton_protein_gff)


##############################
# Creating database
##############################
dbfn_ref_gff = ref_gff + "_db"
dbfn_liftoff_gff = liftoff_gff + "_db"
dbfn_lifton_protein_gff = lifton_protein_gff + "_db"

# Create the GFF database
# db_MANE_ref_gff = gffutils.create_db(
#     MANE_ref_gff,
#     dbfn=dbfn_liftoff_gff,
#     merge_strategy="merge",
#     force=True,  # Set this to True to overwrite the database if it already exists    
# )

db_liftoff_gff = gffutils.create_db(
    liftoff_gff,
    dbfn=dbfn_liftoff_gff,
    merge_strategy="merge",
    force=True,  # Set this to True to overwrite the database if it already exists    
)

db_lifton_protein_gff = gffutils.create_db(
    lifton_protein_gff,
    dbfn=dbfn_lifton_protein_gff,
    merge_strategy="merge",
    force=True,  # Set this to True to overwrite the database if it already exists    
)

# # Load the GFF database
db_ref_gff = gffutils.FeatureDB(dbfn_ref_gff)
# db_liftoff_gff = gffutils.FeatureDB(dbfn_liftoff_gff)
# db_lifton_protein_gff = gffutils.FeatureDB(dbfn_lifton_protein_gff)



frname_both = output_dir + "identities.txt"
frname_liftoff_only = output_dir + "liftoff_only.txt"
frname_lifton_only = output_dir + "lifton_only.txt"

fwname_both = output_dir + "identities_coords.txt"
fwname_liftoff_only = output_dir + "liftoff_only_coords.txt"
fwname_lifton_only = output_dir + "lifton_only_coords.txt"



# frname_miss_both = output_dir + "miss_both.txt"

fr_both = open(frname_both, "r")
fr_liftoff_only = open(frname_liftoff_only, "r")
fr_lifton_only = open(frname_lifton_only, "r")

frs = [fr_both, fr_liftoff_only, fr_lifton_only]
fwfnames = [fwname_both, fwname_liftoff_only, fwname_lifton_only]

for idx in range(3):
    fr = frs[idx]
    fwfname = fwfnames[idx]

    with open(fwfname, "w") as fw:
        lines = fr.read().splitlines()
        for line in lines:
            # # entry = str(line)
            # print(line.start)
            trans_id = line.split("\t")[0]
            # print(trans_id)
            
            if idx == 0:
                # Liftoff database and ID
                eles_ref_db = str(db_ref_gff[trans_id]).split("\t")
                # print("eles_MANE_db: ", eles_MANE_db)
                coord_MANE = f"\t{eles_ref_db[0]}:{eles_ref_db[3]}-{eles_ref_db[4]}"

                eles = str(db_liftoff_gff[trans_id]).split("\t")
                coord = f"\t{eles[0]}:{eles[3]}-{eles[4]}"
                # fw.write(coord)
                # print(coord)

                # lifton database and ID
                coords = ""                
                eles = str(db_lifton_protein_gff[trans_id]).split("\t")
                coords = coords + f"{eles[0]}:{eles[3]}-{eles[4]};"
                fw.write(f"{line}\t{coord}\t{coords}\n")
            
            elif idx == 1:
                # Liftoff database and ID
                eles = str(db_liftoff_gff[trans_id]).split("\t")
                coord = f"{trans_id}\t{eles[0]}:{eles[3]}-{eles[4]}\n"
                fw.write(coord)
                # print(coord)

            elif idx == 2:
                # lifton database and ID
                # print(db_liftoff_gff[trans_id])
                coords = ""
                eles = str(db_lifton_protein_gff[trans_id]).split("\t")
                coords = coords + f"{eles[0]}:{eles[3]}-{eles[4]};"
                fw.write(f"{trans_id}\t{coords}\n")


frname_lifton_only = output_dir + "lifton_only.txt"
fr_lifton_only = open(frname_lifton_only, "r")
