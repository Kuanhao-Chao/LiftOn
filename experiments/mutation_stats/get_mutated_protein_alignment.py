import shutil
import sys, os
import pandas as pd

TARGET = sys.argv[1]

protein_fa = ""
if TARGET == "CHM13_MANE" or TARGET == "CHM13_RefSeq" or TARGET == "GRCh38_RefSeq" or TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1" or TARGET == "Mus_musculus_MANE":
    print("Running with ", TARGET)
    genome = ""
    if TARGET == "GRCh38_RefSeq":
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
mutation_dir = outdir_root + "mutations/" 

fname = f"{outdir_root}score.txt"
# outdir_root = f"{outdir_root}visualization/"

table = pd.read_csv(fname, sep="\t", header=None)

mutation_ls = [
    "identical", 
    "synonymous", 
    "frameshift",
    "5'_truncated",
    "3'_truncated",
    "start_lost",
    "inframe_insertion",
    # "nonsynonymous",
    "stop_missing",
    "stop_codon_gain"]


aln_dir = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/visualization/aln/"

new_aln_dir = aln_dir + "mutations/"

for mutation in mutation_ls:
    mutation_txt = mutation_dir + mutation + ".txt"
    print(">>> mutation_txt: ", mutation_txt)


    mutation_aln_dir = new_aln_dir + mutation + "/"

    os.makedirs(mutation_aln_dir, exist_ok=True)

    with open(mutation_txt, "r") as fr:
        lines = fr.read().splitlines()
        for line in lines:
            src = f"{aln_dir}{line}.png"
            if not os.path.exists(src):
                continue
            dst = f"{mutation_aln_dir}{line}.png"
            shutil.copyfile(src, dst)

            # print(f"\t{line}")
            

