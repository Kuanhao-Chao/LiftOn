import parasail
import os
import pybedtools

import sys
from gffutils import FeatureDB
from Bio import SeqIO
import subprocess
import matplotlib.pyplot as plt

def get_id_fraction(reference, target):
    matches = 0
    for i, letter in enumerate(reference):
        if letter == target[i]:
            matches += 1
    return matches, max(len(reference), len(target))

TARGET = "human_refseq_test"


# outdir_root = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/liftoff_lifton_cmp/"
outdir_root = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/ref_chm13_cmp/"


os.makedirs(outdir_root + "images/", exist_ok=True)

fname = outdir_root + "identities.txt"


targets = ["liftoff", "lifton"]
for idx in range(len(targets)):

    nums = []
    with open(fname, "r") as fr:
        lines = fr.read().splitlines()
        for line in lines:
            eles = line.split("\t")
            nums.append(float(eles[idx+1]))
            # print(eles[1])

    figure_path = outdir_root + "images/"+targets[idx]+"_frequency.png"

    plt.hist(nums, bins=100)
    plt.gca().set(title=f'{targets[idx]} score frequency histogram', ylabel='Frequency', xlabel='Protein pairwise alignment identity score')

    # plt.xlabel('liftoff score')
    # plt.ylabel('lifton score')
    # plt.title('Comparing liftoff vs lifton protein searching scores')
    plt.tight_layout()
    plt.savefig(figure_path, dpi=300)
    plt.close()
    plt.clf()