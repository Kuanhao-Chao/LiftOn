import sys, os
import pandas as pd
import matplotlib.pyplot as plt

TARGET = sys.argv[1]

protein_fa = ""

if TARGET == "human_to_chimp" or TARGET == "mouse_to_rat" or TARGET == "drosophila" or TARGET == "yeast" or TARGET == "arabadop" or TARGET == "bee" or TARGET == "mouse" or TARGET == "rice" or TARGET == "human_mane" or TARGET == "human_chess"  or TARGET == "human_refseq" or TARGET == "human_to_chimp_test" or TARGET == "mouse_to_rat_test" or TARGET == "drosophila_test" or TARGET == "yeast_test" or TARGET == "arabadop_test" or TARGET == "bee_test" or TARGET == "mouse_test" or TARGET == "rice_test" or TARGET == "human_mane_test" or TARGET == "human_chess_test"  or TARGET == "human_refseq_test" or TARGET == "drosophila_erecta_test" or TARGET == "human_mane_to_mouse_test" or TARGET == "human_refseq_to_mouse_test" or TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1":
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
# outdir_root = f"{outdir_root}visualization/"

table = pd.read_csv(fname, sep="\t", header=None)

print("table: ", table)

mutation_ls = [
    "nc_transcript",
    "no_ref_protein",
    "identical",
    "synonymous",
    "nonsynonymous",
    "inframe_insertion",
    "inframe_deletion",
    "frameshift",
    "5'_truncated",
    "3'_truncated",
    "start_lost",
    "stop_missing",
    "stop_codon_gain"
]

# 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 
threasholds = [0.95]

for threashold in threasholds:

    table_lcl = table[table[4] <= threashold]

    for type in ["repeat", "non_repeat"]:

        dict_mutation_count = {
            "nc_transcript": 0,
            "no_ref_protein": 0, 
            "identical": 0,
            "synonymous": 0,
            "frameshift": 0,
            "5'_truncated": 0,
            "3'_truncated": 0,
            "start_lost": 0,
            "inframe_insertion": 0,
            "inframe_deletion": 0,
            "nonsynonymous": 0,
            "stop_missing": 0,
            "stop_codon_gain": 0
        }


        dict_mutation_scores = {
            "nc_transcript": [],
            "no_ref_protein": [], 
            "identical": [],
            "synonymous": [],
            "frameshift": [],
            "5'_truncated": [],
            "3'_truncated": [],
            "start_lost": [],
            "inframe_insertion": [],
            "inframe_deletion": [],
            "nonsynonymous": [],
            "stop_missing": [],
            "stop_codon_gain": []
        }

        dict_mutation_fw = {
            "nc_transcript": None,
            "no_ref_protein": None, 
            "identical": None,
            "synonymous": None,
            "frameshift": None,
            "5'_truncated": None,
            "3'_truncated": None,
            "start_lost": None,
            "inframe_insertion": None,
            "inframe_deletion": None,
            "nonsynonymous": None,
            "stop_missing": None,
            "stop_codon_gain": None
        }

        mutation_dir = f"{outdir_root}mutations/{str(threashold)}/{type}/" 

        os.makedirs(mutation_dir, exist_ok=True)
        for ele in dict_mutation_count:
            outfile = mutation_dir + ele + ".txt"
            dict_mutation_fw[ele] = open(outfile, "w")



        for index, row in table_lcl.iterrows():
            mutations = row[6]
            eles = mutations.split(";")
            # if eles[0] == "nc_transcript":

            if type == "non_repeat":
                max_mutation_idx = 0
                for ele in eles:
                    max_mutation_idx = max(max_mutation_idx, mutation_ls.index(ele))

                    # dict_mutation_count[ele] += 1
                    # dict_mutation_fw[ele].write(row[0] + "\n")

                dict_mutation_count[mutation_ls[max_mutation_idx]] += 1
                dict_mutation_scores[mutation_ls[max_mutation_idx]].append(float(row[4]))

                dict_mutation_fw[mutation_ls[max_mutation_idx]].write(row[0] + "\n")
            
            elif type == "repeat":
                for ele in eles:

                    dict_mutation_count[ele] += 1
                    dict_mutation_scores[ele].append(float(row[4]))
                    dict_mutation_fw[ele].write(row[0] + "\n")




        fw = open(mutation_dir + "summary.txt", "w")

        all_scores = []
        for ele in dict_mutation_count:
            dict_mutation_fw[ele].close()

            fw.write(f"{ele}\t{dict_mutation_count[ele]}\n")

            if ele not in ["nc_transcript", "no_ref_protein", "identical"]:
                all_scores += dict_mutation_scores[ele]

            if threashold == 1:
                figure_path = f"{mutation_dir}/{ele}.png"

                plt.hist(dict_mutation_scores[ele], bins=100)
                plt.gca().set(title='Score frequency histogram', ylabel='Frequency')

                # plt.xlabel('lifton score')
                # plt.ylabel('miniprot score')
                # plt.title('Comparing lifton vs miniprot protein searching scores')
                plt.tight_layout()
                plt.savefig(figure_path, dpi=300)
                plt.close()
                plt.clf()
        
        
        figure_path = f"{outdir_root}mutations/{str(threashold)}/frequency.png"





        # 30 points between [0, 0.2) originally made using np.random.rand(30)*.2
        f, (ax, ax2) = plt.subplots(2, 1, sharex=True)

        # plot the same data on both axes
        ax.hist(all_scores, bins=100, edgecolor='black', alpha=0.7)
        ax2.hist(all_scores, bins=100, edgecolor='black', alpha=0.7)

        # ax.plot(select_scores)
        # ax2.plot(select_scores)

        # zoom-in / limit the view to different portions of the data
        ax.set_ylim(29300, 29500)  # outliers only
        ax2.set_ylim(0, 80)  # most of the data

        # hide the spines between ax and ax2
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop=False)  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()

        # This looks pretty good, and was fairly painless, but you can get that
        # cut-out diagonal lines look with just a bit more work. The important
        # thing to know here is that in axes coordinates, which are always
        # between 0-1, spine endpoints are at these locations (0,0), (0,1),
        # (1,0), and (1,1).  Thus, we just need to put the diagonals in the
        # appropriate corners of each of our axes, and so long as we use the
        # right transform and disable clipping.

        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal






        # plt.hist(all_scores, bins=100)
        # plt.gca().set(title='Score frequency histogram', ylabel='Frequency')

        # # plt.xlabel('lifton score')
        # # plt.ylabel('miniprot score')
        # # plt.title('Comparing lifton vs miniprot protein searching scores')
        plt.tight_layout()
        plt.savefig(figure_path, dpi=300)
        plt.close()
        plt.clf()


        fw.close()

        print("Total number of mutated proteins: ", len(all_scores))