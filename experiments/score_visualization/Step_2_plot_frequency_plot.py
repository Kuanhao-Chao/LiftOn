import parasail
import os
import pybedtools

import sys
from gffutils import FeatureDB
from Bio import SeqIO
import subprocess
import matplotlib.pyplot as plt
import pandas as pd


TARGET = sys.argv[1]

protein_fa = ""
if TARGET == "human_to_chimp" or TARGET == "mouse_to_rat" or TARGET == "yeast" or TARGET == "arabadop" or TARGET == "bee" or TARGET == "mouse" or TARGET == "rice" or TARGET == "CHM13_MANE" or TARGET == "CHM13_RefSeq" or TARGET == "GRCh38_RefSeq" or TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1" or TARGET == "Mus_musculus_MANE":
    print("Running with ", TARGET)
    genome = ""

    if TARGET == "human_to_chimp": 
        genome = "/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3-v1.1.fna"
    
    elif TARGET == "mouse_to_rat": 
        genome = "/ccb/salz2/jheinz3/shared/lifton/cross_species/mouse_to_rat/mRatBN7.2_genomic.fna"

    elif TARGET == "yeast": 
        genome = "/ccb/salz2/kh.chao/Lifton/data/yeast/S288C_reference_genome_R64-4-1_20230830/S288C_reference_sequence_R64-4-1_20230830.fsa"

    elif TARGET == "arabadop": 
        genome = "/ccb/salz2/jheinz3/shared/lifton/arabadop/Tanz-1.fna"

    elif TARGET == "bee": 
        genome = "/ccb/salz2/jheinz3/shared/lifton/bee/ASM1932182v1_genomic.fna"

    elif TARGET == "mouse": 
        genome = "/ccb/salz2/jheinz3/shared/lifton/mouse/mMusMuc1.1_genomic.fna"

    elif TARGET == "rice": 
        genome = "/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_genomic.fna"


    elif TARGET == "GRCh38_RefSeq":
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"

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



outdir_root = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/"
fname = f"{outdir_root}score.txt"
outdir_root = f"{outdir_root}visualization/"

table = pd.read_csv(fname, sep="\t", header=None)


for target in ["lifton", "miniprot", "liftoff"]:
    # target = "lifton"

    figure_path = f"{outdir_root}{target}_frequency.png"

    if target == "liftoff":
        select_scores = table[1][table[1] <0.99]
    elif target == "miniprot":
        select_scores = table[2][table[2] <0.99]
    elif target == "lifton":
        select_scores = table[3][table[3] <0.99]

    plt.hist(select_scores, bins=100)
    plt.gca().set(title='Score frequency histogram', ylabel='Frequency')

    # plt.xlabel('lifton score')
    # plt.ylabel('miniprot score')
    # plt.title('Comparing lifton vs miniprot protein searching scores')
    plt.tight_layout()
    plt.savefig(figure_path, dpi=300)
    plt.close()
    plt.clf()


    for i in range(20):
        threshold = i/20
        selected_num = len(select_scores[select_scores < threshold]) 
        print(f"* < {threshold}: {selected_num}")

    print(f"* < 1: {len(select_scores)}")