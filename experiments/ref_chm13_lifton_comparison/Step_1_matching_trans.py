import pandas as pd
import matplotlib.pyplot as plt

def main():
    lifton_table = pd.read_csv("/ccb/salz2/kh.chao/Lifton/results/human_refseq_test/eval.txt", sep="\t", header=None)

    ref_table = pd.read_csv("/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/eval.txt", sep="\t", header=None)
    
    print("lifton_table: ", len(lifton_table))
    print("ref_table   : ", len(ref_table))

    print("lifton_table: ", lifton_table)
    print("ref_table: ", ref_table)

    ref_table[0] = 'rna-' + ref_table[0].astype(str)

    # Extract values from the second column of both tables
    lifton_values = lifton_table[0].values
    ref_values = ref_table[0].values

    # Find the intersection of values
    intersection_values = set(lifton_values) & set(ref_values)

    print("len(lifton_values): ", len(set(lifton_values)))
    print("len(ref_values): ", len(set(ref_values)))

    print("len(intersection_values): ", len(intersection_values))
    # print("intersection_values     : ", intersection_values)

    # Filter rows based on the intersection values
    lifton_filtered = lifton_table[lifton_table[0].isin(intersection_values)]
    ref_filtered = ref_table[ref_table[0].isin(intersection_values)]




    print("lifton_table (filtered): ", len(lifton_filtered))
    print("ref_table (filtered)   : ", len(ref_filtered))


    # Merge based on the values in the first column
    merged_df = pd.merge(lifton_filtered, ref_filtered, on=0, how='inner')

    # Print the merged result
    print("Merged DataFrame:")
    print(merged_df)
    print("Merged DataFrame Length:", len(merged_df))

    # for trans in intersection_values:
    #     lifton_selected = lifton_filtered[lifton_filtered[0] == trans]
    #     if len(lifton_selected) > 1:
    #         print("Lifton selected")
    #         print(lifton_selected)


    #     ref_selected = ref_filtered[ref_filtered[0] == trans]
    #     if len(ref_selected) > 1:
    #         print("Ref selected")
    #         print(ref_selected)



    # print("lifton_table (filtered): ", lifton_filtered)
    # print("ref_filtered (filtered): ", ref_filtered)




    targets = ["LiftOn", "Reference"]
    for target in targets:
        for type in ["dna", "protein"]:
            # nums = []
            # with open(fname, "r") as fr:
            #     lines = fr.read().splitlines()
            #     for line in lines:
            #         eles = line.split("\t")
            #         nums.append(float(eles[idx+1]))
            #         # print(eles[1])

            # figure_path = outdir_root + "images/"+targets[idx]+"_frequency.png"
            figure_out = "/ccb/salz2/kh.chao/Lifton/results/human_refseq_test/ref_chm13_cmp/"

            if target == "LiftOn" and type == "dna":
                select_score = merged_df["1_x"]
            elif target == "LiftOn" and type == "protein":
                select_score = merged_df["1_y"]
            elif target == "Reference" and type == "dna":
                select_score = merged_df["2_x"]
            elif target == "Reference" and type == "protein":
                select_score = merged_df["2_y"]

            select_score = select_score[(select_score > 0) & (select_score < 0.98) ]
            plt.hist(select_score, bins=100)
            plt.gca().set(title=f'{target} {type} score frequency histogram', ylabel='Frequency', xlabel='Protein pairwise alignment identity score')

            # plt.xlabel('liftoff score')
            # plt.ylabel('lifton score')
            # plt.title('Comparing liftoff vs lifton protein searching scores')
            plt.tight_layout()
            plt.savefig(f"{figure_out}/{target}_{type}_frequency.png", dpi=300)
            plt.close()
            plt.clf()


    # # Plot scatter plot:
    # figure_out = "/ccb/salz2/kh.chao/Lifton/results/human_refseq_test/ref_chm13_cmp/"
    # for type in ["dna", "protein"]:
    #     print(f"type: {type}")
    #     if type == "dna":
    #         plt.scatter(merged_df["1_x"], merged_df["2_x"], s=2)
    #         print("\tLiftOn > ref: ", len(merged_df[merged_df["1_x"] > merged_df["2_x"]]))
    #         print("\tLiftOn = ref: ", len(merged_df[merged_df["1_x"] == merged_df["2_x"]]))
    #         print("\tLiftOn < ref: ", len(merged_df[merged_df["1_x"] < merged_df["2_x"]]))

    #     elif type == "protein":
    #         plt.scatter(merged_df["1_y"], merged_df["2_y"], s=2)
    #         print("\tLiftOn > ref: ", len(merged_df[merged_df["1_y"] > merged_df["2_y"]]))
    #         print("\tLiftOn = ref: ", len(merged_df[merged_df["1_y"] == merged_df["2_y"]]))
    #         print("\tLiftOn < ref: ", len(merged_df[merged_df["1_y"] < merged_df["2_y"]]))
    #     # # Add labels to the points
    #     # for i, row in table.iterrows():
    #     #     if (abs(row[1] - row[2]) > 0.2):
    #     #         plt.text(row[1], row[2], row[0], fontsize=4, ha='center', va='bottom')

    #     plt.axline((0, 0), (1, 1), linewidth=1, color='r')
    #     plt.axis('square')
    #     plt.xlim(-0.1, 1.1)
    #     plt.ylim(-0.1, 1.1)
    #     plt.xlabel('liftoff score')
    #     plt.ylabel('Lifton score')
    #     plt.title('Comparing liftoff vs lifton protein searching scores')
    #     plt.tight_layout()
    #     plt.savefig(f"{figure_out}/scatter_{type}.png", dpi=300)
    #     plt.close()
    #     plt.clf()




if __name__ == "__main__":
    main()