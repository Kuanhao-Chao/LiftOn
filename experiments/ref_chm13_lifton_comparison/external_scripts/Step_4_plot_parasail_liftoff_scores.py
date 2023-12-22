import os
import pandas as pd
import sys
import matplotlib.pyplot as plt

TARGET = "human_refseq_test"


# outdir_root = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/liftoff_lifton_cmp/"
outdir_root = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/ref_chm13_cmp/"


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
plt.ylabel('Lifton score')
plt.title('Comparing liftoff vs lifton protein searching scores')
plt.tight_layout()
plt.savefig(figure_out, dpi=300)
plt.close()
plt.clf()

