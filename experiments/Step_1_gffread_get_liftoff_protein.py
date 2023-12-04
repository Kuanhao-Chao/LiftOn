from gffutils import FeatureDB
import pybedtools
import os

import sys
from gffutils import FeatureDB
from Bio import SeqIO
import subprocess

TARGET = sys.argv[1]

if TARGET == "CHM13_MANE" or TARGET == "CHM13_RefSeq" or TARGET == "GRCh38_RefSeq" or TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1" or TARGET == "Mus_musculus_MANE" or TARGET == "human_refseq_test":
    print("Running with ", TARGET)
    genome = ""

    if TARGET == "GRCh38_RefSeq":
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"

    elif TARGET == "human_refseq_test":
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"

    elif TARGET == "CHM13_MANE": 
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"

    elif TARGET == "CHM13_RefSeq":
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"

    elif TARGET == "Han1":
        genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/Han1/v1.0/Assembly/Han1_v1.2.fasta"

    elif TARGET == "Ash1":
        genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/Ash1/v2.2/Assembly/Ash1_v2.2.fa"

    elif TARGET == "PR1":
        genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/PR1/v3.0/Assembly/PR1.fa"

    elif TARGET == "Mus_musculus_MANE":
        genome = "/ccb/salz3/kh.chao/ref_genome/mus_musculus/NCBI_Refseq/GCF_000001635.27_GRCm39_genomic.fna"
else:
    sys.exit(-1)

refernece_gff3 = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/lifton.gff3" 
output_protein_fa = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/lifton_protein.fa"


print("refernece_gff3: ", refernece_gff3)
print("output_protein_fa: ", output_protein_fa)

command = ['gffread', '-S', '-y', output_protein_fa, '-g', genome, refernece_gff3]
print(command)

gffread_extract = subprocess.call(command)