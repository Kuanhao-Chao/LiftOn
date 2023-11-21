import sys
from lifton import logger

def print_report(ref_features_dict, fw_unmapped, fw_extra_copy, debug=False):
    LIFTED_FEATURES = 0
    MISSED_FEATURES = 0

    SINGLE_COPY_FEATURES = 0
    EXTRA_COPY_FEATURES = 0
    EXTRA_COPY_FEATURES_NUM = 0

    TOTOAL_CHILDREN = 0
    LIFTED_CHILDREN = 0
    MISSED_CHILDREN = 0

    SINGLE_COPY_CHILDREN = 0
    EXTRA_COPY_CHILDREN = 0
    EXTRA_COPY_CHILDREN_NUM = 0

    for feature in ref_features_dict.keys():
        if feature == "LiftOn-gene":
            continue
        
        copy_num = ref_features_dict[feature].copy_num 
        if copy_num > 1:
            print(f"{copy_num}\t{feature}")

            EXTRA_COPY_FEATURES += 1
            EXTRA_COPY_FEATURES_NUM += (copy_num - 1)
            LIFTED_FEATURES += 1
            fw_extra_copy.write(f"{feature}\t{copy_num}\n")
        elif copy_num == 1:
            SINGLE_COPY_FEATURES += 1
            LIFTED_FEATURES += 1
        elif copy_num == 0:
            MISSED_FEATURES += 1
            fw_unmapped.write(f"{feature}\n")

        children_num = len(ref_features_dict[feature].children)
        if children_num > 0:
            TOTOAL_CHILDREN += children_num
            if copy_num > 1:
                EXTRA_COPY_CHILDREN += children_num
                EXTRA_COPY_CHILDREN_NUM += children_num*(copy_num-1)
                LIFTED_CHILDREN += children_num
            elif copy_num == 1:
                SINGLE_COPY_CHILDREN += children_num
                LIFTED_CHILDREN += children_num
            elif copy_num == 0:
                MISSED_CHILDREN += children_num

        # for trans in ref_features_dict[gene].children:
        #     print(f"gene: {gene}; trans: {trans}")

    if TOTOAL_CHILDREN == 0:
        print("*********************************************", file=sys.stderr)
        # print("LiftOn report:", file=sys.stderr)
        print("* Transcript level", file=sys.stderr)
        print(f"\t* Total transcript         : {len(ref_features_dict.keys())}", file=sys.stderr)
        print(f"\t* Lifted transcript        : {LIFTED_FEATURES}", file=sys.stderr)
        print(f"\t* Missed transcript        : {MISSED_FEATURES}", file=sys.stderr)
        print(f"\t* Single copy transcript   : {SINGLE_COPY_FEATURES}", file=sys.stderr)
        print(f"\t* Extra copy transcript    : {EXTRA_COPY_FEATURES}", file=sys.stderr)
        print(f"\t   Total extra transcript  : {EXTRA_COPY_FEATURES_NUM}", file=sys.stderr)
        print(f"\t* Novel LiftOn transcript  : {ref_features_dict['LiftOn-gene'].copy_num}", file=sys.stderr)
        print(f"*********************************************")

    elif TOTOAL_CHILDREN > 0:
        print("*********************************************", file=sys.stderr)
        # print("LiftOn report:", file=sys.stderr)
        print("* Gene level", file=sys.stderr)
        print(f"\t* Total gene                : {len(ref_features_dict.keys())}", file=sys.stderr)
        print(f"\t* Lifted gene               : {LIFTED_FEATURES}", file=sys.stderr)
        print(f"\t* Missed gene               : {MISSED_FEATURES}", file=sys.stderr)
        print(f"\t* Single copy gene          : {SINGLE_COPY_FEATURES}", file=sys.stderr)
        print(f"\t* Extra copy gene           : {EXTRA_COPY_FEATURES}", file=sys.stderr)
        print(f"\t   Total extra gene loci    : {EXTRA_COPY_FEATURES_NUM}", file=sys.stderr)
        print(f"\t* Novel LiftOn gene         : {ref_features_dict['LiftOn-gene'].copy_num}", file=sys.stderr)
        # for gene in gene_copy_num_dict.keys():
        #     if gene_copy_num_dict[gene] > 0:
        #         print(f"\t{gene}: {gene_copy_num_dict[gene]}")
        
        print(f"* Transcript level:")
        print(f"\t* Total transcript          : {TOTOAL_CHILDREN}", file=sys.stderr)
        print(f"\t* Lifted transcript         : {LIFTED_CHILDREN}", file=sys.stderr)
        print(f"\t* Missed transcript         : {MISSED_CHILDREN}", file=sys.stderr)
        print(f"\t* Single copy transcript    : {SINGLE_COPY_CHILDREN}", file=sys.stderr)
        print(f"\t* Extra copy transcript     : {EXTRA_COPY_CHILDREN}", file=sys.stderr)
        print(f"\t   Total extra transcript   : {EXTRA_COPY_CHILDREN_NUM}", file=sys.stderr)
        print(f"\t* Novel LiftOn transcript   : {ref_features_dict['LiftOn-gene'].copy_num}", file=sys.stderr)
        print(f"*********************************************")
    