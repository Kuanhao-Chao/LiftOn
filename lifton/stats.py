import sys
from lifton import logger

def print_report(ref_features_dict, gene_copy_num_dict, trans_copy_num_dict, fw_unmapped, fw_extra_copy, debug=False):

    TOTOAL_GENES = 0
    TOTOAL_TRANS = 0

    LIFTED_GENES = 0
    MISSED_GENES = 0

    LIFTED_TRANS = 0
    MISSED_TRANS = 0

    for gene in ref_features_dict.keys():
        TOTOAL_GENES += 1
        
        gene_lifted = False

        for trans in ref_features_dict[gene].keys():

            TOTOAL_TRANS += 1
            if ref_features_dict[gene][trans]:
                LIFTED_TRANS += 1
                gene_lifted = True
            else:
                MISSED_TRANS += 1
        
        if gene_lifted:
            LIFTED_GENES += 1
        else:
            MISSED_GENES += 1
            fw_unmapped.write(f"{gene}\n")

    SINGLE_COPY_GENES = 0
    EXTRA_COPY_GENES = 0
    EXTRA_COPY_GENES_NUM = 0

    SINGLE_COPY_TRANS = 0
    EXTRA_COPY_TRANS = 0
    EXTRA_COPY_TRANS_NUM = 0

    for gene in gene_copy_num_dict.keys():
        if gene == "gene-LiftOn":
            continue
        if gene_copy_num_dict[gene] > 0:
            EXTRA_COPY_GENES += 1
            EXTRA_COPY_GENES_NUM += gene_copy_num_dict[gene]
            fw_extra_copy.write(f"{gene}\t{gene_copy_num_dict[gene]}\n")
        else:
            SINGLE_COPY_GENES += 1

    for trans in trans_copy_num_dict.keys():
        if trans_copy_num_dict[trans] > 0:
            EXTRA_COPY_TRANS += 1
            EXTRA_COPY_TRANS_NUM += trans_copy_num_dict[trans]
        else:
            SINGLE_COPY_TRANS += 1

    if (TOTOAL_GENES - MISSED_GENES) != (EXTRA_COPY_GENES + SINGLE_COPY_GENES):
        logger.log("Error: ", TOTOAL_GENES-MISSED_GENES, " != ", (EXTRA_COPY_GENES + SINGLE_COPY_GENES), debug=debug)
    if (TOTOAL_TRANS - MISSED_TRANS) != (EXTRA_COPY_TRANS + SINGLE_COPY_TRANS):
        logger.log("Error: ", TOTOAL_TRANS-MISSED_TRANS, " != ", (EXTRA_COPY_TRANS + SINGLE_COPY_TRANS), debug=debug)


    print("*********************************************", file=sys.stderr)
    print("LiftOn report:", file=sys.stderr)
    print("* Gene level", file=sys.stderr)
    print(f"\t* Total genes            : {TOTOAL_GENES}", file=sys.stderr)
    print(f"\t* Lifted genes           : {LIFTED_GENES}", file=sys.stderr)
    print(f"\t* Missed genes           : {MISSED_GENES}", file=sys.stderr)
    print(f"\t* Single copy genes      : {SINGLE_COPY_GENES}", file=sys.stderr)
    print(f"\t* Extra copy genes       : {EXTRA_COPY_GENES}", file=sys.stderr)
    print(f"\t   Total extra gene loci : {EXTRA_COPY_GENES_NUM}", file=sys.stderr)
    # for gene in gene_copy_num_dict.keys():
    #     if gene_copy_num_dict[gene] > 0:
    #         print(f"\t{gene}: {gene_copy_num_dict[gene]}")
    
    print(f"* Transcript level:")
    print(f"\t* Total trans            : {TOTOAL_TRANS}", file=sys.stderr)
    print(f"\t* Lifted trans           : {LIFTED_TRANS}", file=sys.stderr)
    print(f"\t* Missed trans           : {MISSED_TRANS}", file=sys.stderr)
    print(f"\t* Single copy trans      : {SINGLE_COPY_TRANS}", file=sys.stderr)
    print(f"\t* Extra copy trans       : {EXTRA_COPY_TRANS}", file=sys.stderr)
    print(f"\t   Total extra trans     : {EXTRA_COPY_TRANS_NUM}", file=sys.stderr)
    print(f"*********************************************")
    # for trans in trans_copy_num_dict.keys():
    #     if trans_copy_num_dict[trans] > 0:
    #         print(f"\t{trans}: {trans_copy_num_dict[trans]}")
