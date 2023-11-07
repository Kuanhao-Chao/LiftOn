from pymsaviz import MsaViz, get_msa_testdata
import os, sys

TARGET = sys.argv[1]

mutation_ls = [
    # "identical", 
    # "synonymous", 
    "frameshift",
    # "5'_truncated",
    "3'_truncated",
    "start_lost",
    "inframe_insertion",
    # "nonsynonymous",
    "stop_missing",
    "stop_codon_gain"]

for target in ["DNA", "AA"]:
    path = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/lifton_{target}/"
    dir_list = os.listdir(path)

    for mutation in mutation_ls:
        mutation_aln_dir = path + mutation + "/"

        for target_fa in os.listdir(mutation_aln_dir):
            target = os.path.splitext(os.path.basename(target_fa))[0]
            output_file = f"{mutation_aln_dir}{target}.png"

            # print("output_file: ", output_file)
            # if os.path.exists(output_file):
            #     continue

            try:
                msa_file = f"{mutation_aln_dir}{target_fa}"
                print("msa_file: ", msa_file)
                mv = MsaViz(msa_file, wrap_length=60, show_count=True, show_consensus= True)
                os.makedirs(mutation_aln_dir, exist_ok=True)
                mv.savefig(output_file)
            except:
                print("Error: ", target_fa)
                continue
        