import os
import pandas as pd
import sys
import matplotlib.pyplot as plt

TARGET = sys.argv[1]

if TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1" or TARGET == "CHM13_MANE" or TARGET == "CHM13_GRCh38" or TARGET == "GRCh38_RefSeq" or TARGET == "GRCh38_MANE" or TARGET == "Mus_musculus":
    print("Running with ", TARGET)
    genome = ""
    if TARGET == "Han1":
        genome = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/Han1/v1.0/Assembly/Han1_v1.2.fasta"
    elif TARGET == "Ash1":
        genome = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/Ash1/v2.2/Assembly/Ash1_v2.2.fa"
    elif TARGET == "PR1":
        genome = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/PR1/v3.0/Assembly/PR1.fa"
    elif TARGET == "GRCh38_RefSeq" or TARGET == "GRCh38_MANE":
        genome = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"
    elif TARGET == "CHM13_MANE" or TARGET == "CHM13_GRCh38":
        genome = "/ccb/salz2/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
    elif TARGET == "Mus_musculus":
        genome = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/mus_musculus/NCBI_Refseq/GCF_000001635.27_GRCm39_genomic.fna"
else:
    sys.exit(-1)

outdir_root = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/{TARGET}/MANE_liftoff_miniprot_cmp/"


os.makedirs(outdir_root + "images/", exist_ok=True)

fname = outdir_root + "identities.txt"
figure_out = outdir_root + "images/parasail_identities.png"
table = pd.read_csv(fname, sep="\t", header=None)
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
plt.xlabel('liftoff score')
plt.ylabel('miniprot score')
plt.title('Comparing liftoff vs miniprot protein searching scores')
plt.tight_layout()
plt.savefig(figure_out, dpi=300)
plt.close()
plt.clf()

