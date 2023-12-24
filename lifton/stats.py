import sys
from lifton import logger

def print_report(ref_features_dict, fw_unmapped, fw_extra_copy, debug=False):
    LIFTED_FEATURES = 0
    MISSED_FEATURES = 0

    LIFTED_SINGLE_CODING_FEATURES = 0
    LIFTED_SINGLE_NONCODING_FEATURES = 0
    
    LIFTED_EXTRA_CODING_FEATURES = 0
    LIFTED_EXTRA_NONCODING_FEATURES = 0

    LIFTED_EXTRA_CODING_SUM_FEATURES = 0
    LIFTED_EXTRA_NONCODING_SUM_FEATURES = 0

    
    
    # SINGLE_COPY_FEATURES = 0
    # EXTRA_COPY_FEATURES = 0
    # EXTRA_COPY_FEATURES_NUM = 0

    # TOTOAL_CHILDREN = 0
    # LIFTED_CHILDREN = 0
    # MISSED_CHILDREN = 0

    # SINGLE_COPY_CHILDREN = 0
    # EXTRA_COPY_CHILDREN = 0
    # EXTRA_COPY_CHILDREN_NUM = 0

    fw_gene = open("gene.txt", "w")
    fw_trans = open("trans.txt", "w")

    for feature in ref_features_dict.keys():
        if feature == "LiftOn-gene":
            continue

        fw_gene.write(f"{feature}\n")
        copy_num = ref_features_dict[feature].copy_num 

        if copy_num >= 1:
            LIFTED_FEATURES += 1

            if ref_features_dict[feature].is_protein_coding and copy_num == 1:
                LIFTED_SINGLE_CODING_FEATURES += 1
            elif ref_features_dict[feature].is_protein_coding and copy_num > 1:
                fw_extra_copy.write(f"{feature}\t{copy_num}\tcoding\n")
                # print(f"{copy_num}\t{feature}")
                LIFTED_EXTRA_CODING_FEATURES += 1
                LIFTED_EXTRA_CODING_SUM_FEATURES += (copy_num)

            elif not ref_features_dict[feature].is_protein_coding and copy_num == 1:
                LIFTED_SINGLE_NONCODING_FEATURES += 1
            elif not ref_features_dict[feature].is_protein_coding and copy_num > 1:
                fw_extra_copy.write(f"{feature}\t{copy_num}\tnon-coding\n")
                # print(f"{copy_num}\t{feature}")
                LIFTED_EXTRA_NONCODING_FEATURES += 1
                LIFTED_EXTRA_NONCODING_SUM_FEATURES += (copy_num)

        elif copy_num == 0:
            MISSED_FEATURES += 1
            fw_unmapped.write(f"{feature}\n")

        # if copy_num > 1:
        #     print(f"{copy_num}\t{feature}")
        #     EXTRA_COPY_FEATURES += 1
        #     EXTRA_COPY_FEATURES_NUM += (copy_num - 1)
        #     LIFTED_FEATURES += 1
        #     fw_extra_copy.write(f"{feature}\t{copy_num}\n")
        # elif copy_num == 1:
        #     SINGLE_COPY_FEATURES += 1
        #     LIFTED_FEATURES += 1
        # elif copy_num == 0:
        #     MISSED_FEATURES += 1
        #     fw_unmapped.write(f"{feature}\n")

        # children_num = len(ref_features_dict[feature].children)
        # for child in ref_features_dict[feature].children:
        #     fw_trans.write(f"{child}\n")

        # if children_num > 0:
        #     TOTOAL_CHILDREN += children_num
        #     if copy_num > 1:
        #         EXTRA_COPY_CHILDREN += children_num
        #         EXTRA_COPY_CHILDREN_NUM += children_num*(copy_num-1)
        #         LIFTED_CHILDREN += children_num
        #     elif copy_num == 1:
        #         SINGLE_COPY_CHILDREN += children_num
        #         LIFTED_CHILDREN += children_num
        #     elif copy_num == 0:
        #         MISSED_CHILDREN += children_num

        # # for trans in ref_features_dict[gene].children:
        # #     print(f"gene: {gene}; trans: {trans}")

    # if TOTOAL_CHILDREN == 0:
    #     print("*********************************************", file=sys.stderr)
    #     # print("LiftOn report:", file=sys.stderr)
    #     print("* Transcript level", file=sys.stderr)
    #     print(f"\t* Total transcript         : {len(ref_features_dict.keys())-1}", file=sys.stderr)
    #     print(f"\t* Lifted transcript        : {LIFTED_FEATURES}", file=sys.stderr)
    #     print(f"\t* Missed transcript        : {MISSED_FEATURES}", file=sys.stderr)
    #     print(f"\t* Single copy transcript   : {SINGLE_COPY_FEATURES}", file=sys.stderr)
    #     print(f"\t* Extra copy transcript    : {EXTRA_COPY_FEATURES}", file=sys.stderr)
    #     print(f"\t   Total extra transcript  : {EXTRA_COPY_FEATURES_NUM}", file=sys.stderr)
    #     print(f"\t* Novel LiftOn transcript  : {ref_features_dict['LiftOn-gene'].copy_num}", file=sys.stderr)
    #     print(f"*********************************************")

    # elif TOTOAL_CHILDREN > 0:
    print("*********************************************", file=sys.stderr)
    # print("LiftOn report:", file=sys.stderr)
    # print("* Gene level", file=sys.stderr)
    print(f"* Total features in reference\t\t: {len(ref_features_dict.keys())-1}", 
    file=sys.stderr)
    print(f"* Lifted feature\t\t\t: {LIFTED_FEATURES} ({LIFTED_SINGLE_CODING_FEATURES + LIFTED_EXTRA_CODING_FEATURES} + {LIFTED_SINGLE_NONCODING_FEATURES + LIFTED_EXTRA_NONCODING_FEATURES})", file=sys.stderr)
    print(f"\t* Protein-coding feature\t: {LIFTED_SINGLE_CODING_FEATURES + LIFTED_EXTRA_CODING_FEATURES}", file=sys.stderr)
    print(f"\t* Non-coding feature\t\t: {LIFTED_SINGLE_NONCODING_FEATURES + LIFTED_EXTRA_NONCODING_FEATURES}", file=sys.stderr)
    print(f"* Missed feature\t\t\t: {MISSED_FEATURES}\n", file=sys.stderr)


    print(f"* Total features in target\t\t: {LIFTED_SINGLE_CODING_FEATURES + LIFTED_EXTRA_CODING_SUM_FEATURES + LIFTED_SINGLE_NONCODING_FEATURES + LIFTED_EXTRA_NONCODING_SUM_FEATURES + ref_features_dict['LiftOn-gene'].copy_num} ({LIFTED_SINGLE_CODING_FEATURES + LIFTED_EXTRA_CODING_SUM_FEATURES} + {LIFTED_SINGLE_NONCODING_FEATURES + LIFTED_EXTRA_NONCODING_SUM_FEATURES} + {ref_features_dict['LiftOn-gene'].copy_num})", file=sys.stderr)
    print(f"\t* Protein-coding feature\t: {LIFTED_SINGLE_CODING_FEATURES + LIFTED_EXTRA_CODING_SUM_FEATURES} ({LIFTED_SINGLE_CODING_FEATURES} + {LIFTED_EXTRA_CODING_SUM_FEATURES})", file=sys.stderr)    
    print(f"\t\t* single copy\t\t: {LIFTED_SINGLE_CODING_FEATURES}", file=sys.stderr)
    print(f"\t\t* > 1 copy\t\t: {LIFTED_EXTRA_CODING_FEATURES}, {LIFTED_EXTRA_CODING_SUM_FEATURES} in total", file=sys.stderr)

    print(f"\t* Non-coding feature\t\t: {LIFTED_SINGLE_NONCODING_FEATURES + LIFTED_EXTRA_NONCODING_SUM_FEATURES} ({LIFTED_SINGLE_NONCODING_FEATURES} + {LIFTED_EXTRA_NONCODING_SUM_FEATURES})", file=sys.stderr)
    print(f"\t\t* single copy\t\t: {LIFTED_SINGLE_NONCODING_FEATURES}", file=sys.stderr)
    print(f"\t\t* > 1 copy\t\t: {LIFTED_EXTRA_NONCODING_FEATURES}, {LIFTED_EXTRA_NONCODING_SUM_FEATURES} in total", file=sys.stderr)

    print(f"\t* Novel LiftOn feature\t\t: {ref_features_dict['LiftOn-gene'].copy_num}", file=sys.stderr)
    
    # # for gene in gene_copy_num_dict.keys():
    # #     if gene_copy_num_dict[gene] > 0:
    # #         print(f"\t{gene}: {gene_copy_num_dict[gene]}")
    
    # print(f"* Transcript level:")
    # print(f"\t* Total transcript          : {TOTOAL_CHILDREN}", file=sys.stderr)
    # print(f"\t* Lifted transcript         : {LIFTED_CHILDREN}", file=sys.stderr)
    # print(f"\t* Missed transcript         : {MISSED_CHILDREN}", file=sys.stderr)
    # print(f"\t* Single copy transcript    : {SINGLE_COPY_CHILDREN}", file=sys.stderr)
    # print(f"\t* Extra copy transcript     : {EXTRA_COPY_CHILDREN}", file=sys.stderr)
    # print(f"\t   Total extra transcript   : {EXTRA_COPY_CHILDREN_NUM}", file=sys.stderr)
    # print(f"\t* Novel LiftOn transcript   : {ref_features_dict['LiftOn-gene'].copy_num}", file=sys.stderr)
    print(f"*********************************************")
