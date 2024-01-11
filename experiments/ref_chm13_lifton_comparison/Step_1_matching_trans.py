import pandas as pd
import matplotlib.pyplot as plt
import os

def main():
    lifton_table = pd.read_csv("/ccb/salz2/kh.chao/LiftOn/results/human_refseq/lifton_output/eval.txt", sep="\t", header=None)

    ref_table = pd.read_csv("/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/lifton_output/eval.txt", sep="\t", header=None)
    
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

    # Use the mask to filter out rows with chrY
    mask = (merged_df['5_x'].str.contains('chrY') | merged_df['5_y'].str.contains('chrY'))
    merged_df = merged_df[~mask]

    # Print the merged result
    print("Merged DataFrame:")
    print(merged_df)
    print("Merged DataFrame Length:", len(merged_df))

    targets = ["JHU RefSeqv110 + Liftoff v5.1", "LiftOn v1.0.0"]

    THRESHOLD = 0.98   

    # Filtering by DNA scores
    # select_score = merged_df["1_x"]

    # print(f"Filter 0 & {THRESHOLD} => Merged DataFrame Length:", len(merged_df))

    for type in ["DNA", "Protein"]:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7.1, 6), sharey=True)  # Adjust figsize as needed
        # for i, target in enumerate(targets):
        #     figure_out = "/ccb/salz2/kh.chao/Lifton/results/human_refseq_test/ref_chm13_cmp/"

        #     if target == "LiftOn" and type == "dna":
        #         select_score = merged_df["1_x"]
        #     elif target == "LiftOn" and type == "protein":
        #         select_score = merged_df["1_y"]
        #     elif target == "Reference" and type == "dna":
        #         select_score = merged_df["2_x"]
        #     elif target == "Reference" and type == "protein":
        #         select_score = merged_df["2_y"]

        #     select_score = select_score[(select_score > 0) & (select_score < 0.98)]
        #     # select_score = select_score[(select_score > 0)]
        #     axes[i].hist(select_score, bins=100)
        #     axes[i].set(title=f'{target} {type} score frequency histogram',
        #                 ylabel='Frequency', xlabel=f'{type} pairwise alignment identity score')

        # plt.tight_layout()
        # plt.savefig(f"{figure_out}/combined_{type}_frequency.png", dpi=300)
        # plt.clf()

        if type == "DNA":
            select_merged_df = merged_df[((merged_df["1_x"] > 0) & (merged_df["1_y"] > 0))]# & ((merged_df["1_x"] < THRESHOLD) & (merged_df["1_y"] < THRESHOLD))]
        elif type == "Protein":
            select_merged_df = merged_df[((merged_df["2_x"] > 0) & (merged_df["2_y"] > 0))]# & ((merged_df["1_x"] < THRESHOLD) & (merged_df["1_y"] < THRESHOLD))]

        for i, target in enumerate(targets):
            figure_out = "/ccb/salz2/kh.chao/LiftOn/results/human_refseq/ref_chm13_cmp/"
            os.makedirs(figure_out, exist_ok=True)

            if target == "LiftOn v1.0.0" and type == "DNA":
                select_score = select_merged_df["1_x"]
            elif target == "LiftOn v1.0.0" and type == "Protein":
                select_score = select_merged_df["2_x"]
            elif target == "JHU RefSeqv110 + Liftoff v5.1" and type == "DNA":
                select_score = select_merged_df["1_y"]
            elif target == "JHU RefSeqv110 + Liftoff v5.1" and type == "Protein":
                select_score = select_merged_df["2_y"]

            # select_score = select_score[(select_score > 0) & (select_score < THRESHOLD)]

            # Use alpha to make the histograms semi-transparent
            plt.hist(select_score, bins=100, alpha=0.65, label=target, log=True)#, color=colors[i])


        # plt.axvline(THRESHOLD, linestyle='--', color='r', label=f'Threshold {THRESHOLD}')  # Add mean line
        plt.title(f'{type} sequence identity score frequency histogram', fontsize=17)
        plt.ylabel('Frequency', fontsize=16)
        plt.xlabel(f'{type} sequence identity score', fontsize=16)
        plt.legend(fontsize=16)  # Show legend with target labels
        plt.tight_layout()
        plt.savefig(f"{figure_out}/log_combined_{type}_frequency_ovp.png", dpi=300)
        plt.clf()


    # for target in targets:
    #     for type in ["dna", "protein"]:
    #         figure_out = "/ccb/salz2/kh.chao/Lifton/results/human_refseq_test/ref_chm13_cmp/"
    #         if target == "LiftOn" and type == "dna":
    #             select_score = merged_df["1_x"]
    #         elif target == "LiftOn" and type == "protein":
    #             select_score = merged_df["1_y"]
    #         elif target == "Reference" and type == "dna":
    #             select_score = merged_df["2_x"]
    #         elif target == "Reference" and type == "protein":
    #             select_score = merged_df["2_y"]

    #         select_score = select_score[(select_score > 0) & (select_score < 0.98) ]
    #         plt.hist(select_score, bins=100)
    #         plt.gca().set(title=f'{target} {type} score frequency histogram', ylabel='Frequency', xlabel='Protein pairwise alignment identity score')

    #         # plt.xlabel('liftoff score')
    #         # plt.ylabel('lifton score')
    #         # plt.title('Comparing liftoff vs lifton protein searching scores')
    #         plt.tight_layout()
    #         plt.savefig(f"{figure_out}/{target}_{type}_frequency.png", dpi=300)
    #         plt.close()
    #         plt.clf()


    # Plot scatter plot:
    figure_out = "/ccb/salz2/kh.chao/LiftOn/results/human_refseq/ref_chm13_cmp/"
    for type in ["dna", "protein"]:
        print(f"type: {type}")
        if type == "dna":
            plt.scatter(merged_df["1_x"], merged_df["1_y"], s=2)
            print("\tLiftOn > ref: ", len(merged_df[merged_df["1_x"] > merged_df["1_y"]]))
            print("\tLiftOn = ref: ", len(merged_df[merged_df["1_x"] == merged_df["1_y"]]))
            print("\tLiftOn < ref: ", len(merged_df[merged_df["1_x"] < merged_df["1_y"] -0.4]))
            print("\tLiftOn < ref: ", merged_df[merged_df["1_x"] < merged_df["1_y"] -0.4])


        elif type == "protein":
            plt.scatter(merged_df["2_x"], merged_df["2_y"], s=2)
            print("\tLiftOn > ref: ", len(merged_df[merged_df["2_x"] > merged_df["2_y"]]))
            print("\tLiftOn = ref: ", len(merged_df[merged_df["2_x"] == merged_df["2_y"]]))
            print("\tLiftOn < ref: ", len(merged_df[merged_df["2_x"] < merged_df["2_y"]]))
        # # Add labels to the points
        # for i, row in table.iterrows():
        #     if (abs(row[1] - row[2]) > 0.2):
        #         plt.text(row[1], row[2], row[0], fontsize=4, ha='center', va='bottom')

        plt.axline((0, 0), (1, 1), linewidth=1, color='r')
        plt.axis('square')
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.xlabel('LiftOn score')
        plt.ylabel('Reference T2T-CHM13 score')
        plt.title('Comparing liftoff vs lifton protein searching scores')
        plt.tight_layout()
        plt.savefig(f"{figure_out}/scatter_{type}.png", dpi=300)
        plt.close()
        plt.clf()


    merged_df.to_csv(f'{figure_out}/merged_eval.txt', header=False,  index=False, sep="\t") 




if __name__ == "__main__":
    main()