import os
import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px

import matplotlib.pyplot as plt

def randrange(n, vmin, vmax):
    """
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    """
    return (vmax - vmin)*np.random.rand(n) + vmin


TARGET = sys.argv[1]

protein_fa = ""
if TARGET == "human_to_chimp" or TARGET == "mouse_to_rat" or TARGET == "drosophila" or TARGET == "yeast" or TARGET == "arabadop" or TARGET == "bee" or TARGET == "mouse" or TARGET == "rice" or TARGET == "human_mane" or TARGET == "human_chess"  or TARGET == "human_refseq" or TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1":
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


outdir_root = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/"
fname = f"{outdir_root}score.txt"
outdir_root = f"{outdir_root}visualization/"

table = pd.read_csv(fname, sep="\t", header=None)

# 2D scatter plot
for target in ["Liftoff", "miniprot"]:
    os.makedirs(f"{outdir_root}{target}", exist_ok=True)

    figure_out = f"{outdir_root}{target}/parasail_identities.png"

    if target == "Liftoff":

        select_table = table[(table[1] > 0.0) | (table[3] > 0.0)]
        plt.scatter(select_table[1], select_table[3], s=2)
    elif target == "miniprot":
        select_table = table[(table[2] > 0.0) | (table[3] > 0.0)]
        plt.scatter(table[2], table[3], s=2)
    print(table)
    # # Add labels to the points
    # for i, row in table.iterrows():
    #     if (abs(row[1] - row[2]) > 0.2):
    #         plt.text(row[1], row[2], row[0], fontsize=4, ha='center', va='bottom')

    plt.axline((0, 0), (1, 1), linewidth=1, color='r')
    plt.axis('square')
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.1, 1.1)
    plt.xlabel(f'{target} score')
    plt.ylabel('LiftOn score')
    plt.title(f'Comparing LiftOn vs {target} protein searching scores')
    plt.tight_layout()
    plt.savefig(figure_out, dpi=300)
    plt.close()
    plt.clf()



# Combined dot plots
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 2D scatter plot for Liftoff
select_table_liftoff = table[(table[1] > 0.0) | (table[3] > 0.0)]
axes[0].scatter(select_table_liftoff[1], select_table_liftoff[3], s=2)
axes[0].axline((0, 0), (1, 1), linewidth=1, color='r')
axes[0].axis('square')
axes[0].set_xlim(-0.1, 1.1)
axes[0].set_ylim(-0.1, 1.1)
axes[0].set_xlabel('Liftoff score')
axes[0].set_ylabel('LiftOn score')
axes[0].set_title('Comparing LiftOn vs Liftoff protein searching scores')

# 2D scatter plot for miniprot
select_table_miniprot = table[(table[2] > 0.0) | (table[3] > 0.0)]
axes[1].scatter(select_table_miniprot[2], select_table_miniprot[3], s=2)
axes[1].axline((0, 0), (1, 1), linewidth=1, color='r')
axes[1].axis('square')
axes[1].set_xlim(-0.1, 1.1)
axes[1].set_ylim(-0.1, 1.1)
axes[1].set_xlabel('miniprot score')
axes[1].set_ylabel('LiftOn score')
axes[1].set_title('Comparing LiftOn vs miniprot protein searching scores')

# Adjust layout
plt.tight_layout()

# Save or show the figure
plt.savefig(f"{outdir_root}combined_scatter_plots.png", dpi=300)
plt.show()





# 2D scatter plot
target = "Liftoff_miniprot"
os.makedirs(f"{outdir_root}{target}", exist_ok=True)

figure_out = f"{outdir_root}{target}/parasail_identities.png"

plt.scatter(table[1], table[2], s=2)
print(table)
# # Add labels to the points
# for i, row in table.iterrows():
#     if (abs(row[1] - row[2]) > 0.2):
#         plt.text(row[1], row[2], row[0], fontsize=4, ha='center', va='bottom')

plt.axline((0, 0), (1, 1), linewidth=1, color='r')
plt.axis('square')
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 1.1)
plt.xlabel('Liftoff score')
plt.ylabel('miniprot score')
plt.title(f'Comparing Liftoff vs miniprot protein searching scores')
plt.tight_layout()
plt.savefig(figure_out, dpi=300)
plt.close()
plt.clf()



# 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
n = 100

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
# for m, zlow, zhigh in [('o', -50, -25), ('^', -30, -5)]:

labels = ['Liftoff score', 'miniprot score', 'LiftOn score']
xs = table[1]
ys = table[2]
zs = table[3]

# ax.plot_surface(1, 0, 0, alpha=0.2)

ax.scatter(xs, ys, zs, marker="o", s=2)

ax.set_xlabel(labels[0])
ax.set_ylabel(labels[1])
ax.set_zlabel(labels[2])

plt.tight_layout()
plt.savefig(f"{outdir_root}3d_scatter.png", dpi=300)
plt.close()
plt.clf()
