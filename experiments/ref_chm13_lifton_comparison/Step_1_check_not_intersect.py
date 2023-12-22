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



    # Find values only in lifton_values
    only_in_lifton_values = set(lifton_values) - set(ref_values)

    # Find values only in ref_values
    only_in_ref_values = set(ref_values) - set(lifton_values)

    
    print("len(lifton_values): ", len(set(lifton_values)))
    print("len(ref_values): ", len(set(ref_values)))

    print("len(intersection_values): ", len(intersection_values))

    print("Values only in lifton_values:", len(only_in_lifton_values))
    print("Values only in ref_values:", len(only_in_ref_values))

    print("Values only in lifton_values:", only_in_lifton_values)
    print("Values only in ref_values:", only_in_ref_values)


    # Filter rows based on the intersection values
    lifton_filtered = lifton_table[lifton_table[0].isin(intersection_values)]
    ref_filtered = ref_table[ref_table[0].isin(intersection_values)]




    print("lifton_table (filtered): ", len(lifton_filtered))
    print("ref_table (filtered)   : ", len(ref_filtered))



    # # Merge based on the values in the first column
    # merged_df = pd.merge(lifton_filtered, ref_filtered, on=0, how='inner')

    # # Print the merged result
    # print("Merged DataFrame:")
    # print(merged_df)
    # print("Merged DataFrame Length:", len(merged_df))


    # # Left merge to find rows only in lifton_filtered
    # only_in_lifton = pd.merge(lifton_filtered, ref_filtered, on=0, how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)

    # # Right merge to find rows only in ref_filtered
    # only_in_ref = pd.merge(lifton_filtered, ref_filtered, on=0, how='right', indicator=True).query('_merge == "right_only"').drop('_merge', axis=1)

    # # Print the results
    # print("Rows only in lifton_filtered:")
    # print(only_in_lifton)

    # print("\nRows only in ref_filtered:")
    # print(only_in_ref)


if __name__ == "__main__":
    main()