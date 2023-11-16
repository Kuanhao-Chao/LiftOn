import sys
from lifton import logger

def print_report(ref_features_dict, fw_unmapped, fw_extra_copy, debug=False):

    # TOTOAL_GENES = 0
    # TOTOAL_TRANS = 0

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



    # SINGLE_COPY_TRANS = 0
    # EXTRA_COPY_TRANS = 0
    # EXTRA_COPY_TRANS_NUM = 0

    for feature in ref_features_dict.keys():
        if feature == "LiftOn-gene":
            continue
        # TOTOAL_GENES += 1        
        gene_lifted = False
        copy_num = ref_features_dict[feature].copy_num 

        if copy_num > 1:
            EXTRA_COPY_FEATURES += 1
            EXTRA_COPY_FEATURES_NUM += copy_num
            LIFTED_FEATURES += 1
        elif copy_num == 1:
            SINGLE_COPY_FEATURES += 1
            LIFTED_FEATURES += 1
        elif copy_num == 0:
            MISSED_FEATURES += 1


        if len(ref_features_dict[feature].children) > 0:
            children_num = len(ref_features_dict[feature].children)
            TOTOAL_CHILDREN += children_num
            if copy_num > 1:
                EXTRA_COPY_CHILDREN += children_num
                EXTRA_COPY_CHILDREN_NUM += children_num*copy_num
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
        print(f"*********************************************")
    
    


    #         # for trans in ref_features_dict[gene].keys():
    #         TOTOAL_TRANS += 1
    #         if ref_features_dict[gene][trans]:
    #             LIFTED_TRANS += 1
    #             gene_lifted = True
    #         else:
    #             MISSED_TRANS += 1
        
    #     if gene_lifted:
    #         LIFTED_GENES += 1
    #     else:
    #         MISSED_GENES += 1
    #         fw_unmapped.write(f"{gene}\n")

    # SINGLE_COPY_GENES = 0
    # EXTRA_COPY_GENES = 0
    # EXTRA_COPY_GENES_NUM = 0

    # SINGLE_COPY_TRANS = 0
    # EXTRA_COPY_TRANS = 0
    # EXTRA_COPY_TRANS_NUM = 0


    # # for gene in gene_copy_num_dict.keys():
    # #     if gene == "gene-LiftOn":
    # #         continue
    # #     if gene_copy_num_dict[gene] > 0:
    # #         EXTRA_COPY_GENES += 1
    # #         EXTRA_COPY_GENES_NUM += gene_copy_num_dict[gene]
    # #         fw_extra_copy.write(f"{gene}\t{gene_copy_num_dict[gene]}\n")
    # #     else:
    # #         SINGLE_COPY_GENES += 1

    # # for trans in trans_copy_num_dict.keys():
    # #     if trans_copy_num_dict[trans] > 0:
    # #         EXTRA_COPY_TRANS += 1
    # #         EXTRA_COPY_TRANS_NUM += trans_copy_num_dict[trans]
    # #     else:
    # #         SINGLE_COPY_TRANS += 1

    # if (TOTOAL_GENES - MISSED_GENES) != (EXTRA_COPY_GENES + SINGLE_COPY_GENES):
    #     logger.log("Error: ", TOTOAL_GENES-MISSED_GENES, " != ", (EXTRA_COPY_GENES + SINGLE_COPY_GENES), debug=debug)
    # if (TOTOAL_TRANS - MISSED_TRANS) != (EXTRA_COPY_TRANS + SINGLE_COPY_TRANS):
    #     logger.log("Error: ", TOTOAL_TRANS-MISSED_TRANS, " != ", (EXTRA_COPY_TRANS + SINGLE_COPY_TRANS), debug=debug)


    # print("*********************************************", file=sys.stderr)
    # print("LiftOn report:", file=sys.stderr)
    # print("* Gene level", file=sys.stderr)
    # print(f"\t* Total genes            : {TOTOAL_GENES}", file=sys.stderr)
    # print(f"\t* Lifted genes           : {LIFTED_GENES}", file=sys.stderr)
    # print(f"\t* Missed genes           : {MISSED_GENES}", file=sys.stderr)
    # print(f"\t* Single copy genes      : {SINGLE_COPY_GENES}", file=sys.stderr)
    # print(f"\t* Extra copy genes       : {EXTRA_COPY_GENES}", file=sys.stderr)
    # print(f"\t   Total extra gene loci : {EXTRA_COPY_GENES_NUM}", file=sys.stderr)
    # # for gene in gene_copy_num_dict.keys():
    # #     if gene_copy_num_dict[gene] > 0:
    # #         print(f"\t{gene}: {gene_copy_num_dict[gene]}")
    
    # print(f"* Transcript level:")
    # print(f"\t* Total trans            : {TOTOAL_TRANS}", file=sys.stderr)
    # print(f"\t* Lifted trans           : {LIFTED_TRANS}", file=sys.stderr)
    # print(f"\t* Missed trans           : {MISSED_TRANS}", file=sys.stderr)
    # print(f"\t* Single copy trans      : {SINGLE_COPY_TRANS}", file=sys.stderr)
    # print(f"\t* Extra copy trans       : {EXTRA_COPY_TRANS}", file=sys.stderr)
    # print(f"\t   Total extra trans     : {EXTRA_COPY_TRANS_NUM}", file=sys.stderr)
    # print(f"*********************************************")
    
    
    
    # # for trans in trans_copy_num_dict.keys():
    # #     if trans_copy_num_dict[trans] > 0:
    # #         print(f"\t{trans}: {trans_copy_num_dict[trans]}")
