import parasail
import os
import pybedtools

import sys
from gffutils import FeatureDB
from Bio import SeqIO
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

TARGET = sys.argv[1]

protein_fa = ""
if TARGET == "human_to_chimp" or TARGET == "mouse_to_rat" or TARGET == "drosophila" or TARGET == "yeast" or TARGET == "arabadop" or TARGET == "bee" or TARGET == "mouse" or TARGET == "rice" or TARGET == "human_mane" or TARGET == "human_chess"  or TARGET == "human_refseq" or TARGET == "drosophila_erecta" or TARGET == "human_mane_to_mouse" or TARGET == "human_refseq_to_mouse" or TARGET == "human_to_chimp_test" or TARGET == "mouse_to_rat_test" or TARGET == "drosophila_test" or TARGET == "yeast_test" or TARGET == "arabadop_test" or TARGET == "bee_test" or TARGET == "mouse_test" or TARGET == "rice_test" or TARGET == "human_mane_test" or TARGET == "human_chess_test"  or TARGET == "human_refseq_test" or TARGET == "drosophila_erecta_test" or TARGET == "human_mane_to_mouse_test" or TARGET == "human_refseq_to_mouse_test" or TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1":
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

    elif TARGET == "CHM13_MANE" or TARGET == "human_mane": 
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



outdir_root = f"/ccb/salz2/kh.chao/LiftOn/results/{TARGET}/lifton_output/"
fname = f"{outdir_root}score.txt"
outdir_root = f"{outdir_root}visualization/"

table = pd.read_csv(fname, sep="\t", header=None)

upper_threshold = 1.00
for target in ["lifton", "miniprot", "liftoff"]:
    # target = "lifton"

    figure_path = f"{outdir_root}{target}_frequency.png"

    if target == "liftoff":
        select_scores = table[1][(table[1] <= upper_threshold) & (table[1] > 0.0)]
    elif target == "miniprot":
        select_scores = table[2][(table[2] <= upper_threshold) & (table[2] > 0.0)]
    elif target == "lifton":
        select_scores = table[4][(table[4] <= upper_threshold) & (table[4] > 0.0)]

    plt.hist(select_scores, bins=100, alpha=0.9)
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

upper_threshold = 1.00
# upper_threshold = 0.6



count_threshold = 0.4


################################################
# Create a subplot with 1 row and 3 columns
#   taking LOG 
################################################
fig, axes = plt.subplots(1, 3, figsize=(21, 5), sharex='col', sharey='row')

for idx, target in enumerate(["liftoff", "lifton", "miniprot"]):
    select_scores = None

    target_idx = None
    if target == "liftoff":
        target_idx = 1
    elif target == "miniprot":
        target_idx = 2
        # select_scores = table[2][(table[2] <= upper_threshold) & (table[2] > 0.0)]
        # selected_count= len(table[2][(table[2] <= count_threshold) & (table[2] > 0.0)])
    elif target == "lifton":
        target_idx = 4
        # select_scores = table[4][(table[4] <= upper_threshold) & (table[4] > 0.0)]
        # selected_count= len(table[1][(table[1] <= count_threshold) & (table[1] > 0.0)])

    select_scores = table[target_idx][(table[target_idx] <= upper_threshold) & (table[target_idx] > 0.0)]

    selected_count= len(table[target_idx][(table[target_idx] <= count_threshold) & (table[target_idx] > 0.0)])

    print("selected_count: ", selected_count)

    # Plot the histogram on the corresponding subplot
    # axes[idx].hist(select_scores, bins=100, log=True, alpha=0.9)
    axes[idx].hist(select_scores, bins=100, alpha=0.9)

    # Set logarithmic scale for the y-axis
    axes[idx].set_yscale('log')

    # axes[idx].hist(select_scores, bins=np.logspace(np.log10(0.1),np.log10(1.0), 50))


    axes[idx].axvline(x = count_threshold, color = 'r', linestyle='--', linewidth=2)

    axes[idx].annotate(str(selected_count), xy=(0.4, 1000),
                        xytext =(0.2, 850), fontsize=18,   
                       arrowprops = dict(facecolor ='red',  arrowstyle="<-", linewidth=2),) 
    
    axes[idx].set(title=f'{target.capitalize()} Score Frequency Histogram (Log)', xlabel='Log(Score)', ylabel='Frequency')


# Adjust layout
plt.tight_layout()

# Save the figure
figure_path = f"{outdir_root}combined_frequency_log.png"
plt.savefig(figure_path, dpi=300)

# Show the plot
plt.show()

print(f"Length for liftoff: {len(table[1])}")
print(f"Length for miniprot: {len(table[2])}")
print(f"Length for lifton: {len(table[4])}")


print(f"Length for liftoff: {len(table[1][(table[1] <= count_threshold) & (table[1] > 0.0)])}")
print(f"Length for miniprot: {len(table[2][(table[2] <= count_threshold) & (table[2] > 0.0)])}")
print(f"Length for lifton: {len(table[4][(table[4] <= count_threshold) & (table[4] > 0.0)])}")



################################################
# Create a subplot with 1 row and 3 columns
#   not taking LOG 
################################################
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex='col', sharey='row')

for idx, target in enumerate(["liftoff", "lifton", "miniprot"]):
    select_scores = None

    if target == "liftoff":
        select_scores = table[1][(table[1] <= upper_threshold) & (table[1] > 0.0)]
    elif target == "miniprot":
        select_scores = table[2][(table[2] <= upper_threshold) & (table[2] > 0.0)]
    elif target == "lifton":
        select_scores = table[4][(table[4] <= upper_threshold) & (table[4] > 0.0)]

    # Plot the histogram on the corresponding subplot
    axes[idx].hist(select_scores, bins=100, alpha=0.9)
    axes[idx].set(title=f'{target.capitalize()} Score Frequency Histogram', xlabel='Score', ylabel='Frequency')

# Adjust layout
plt.tight_layout()

# Save the figure
figure_path = f"{outdir_root}combined_frequency.png"
plt.savefig(figure_path, dpi=300)

# Show the plot
plt.show()