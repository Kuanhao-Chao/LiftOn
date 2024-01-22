import sys
from lifton import logger

def print_report(ref_features_dict, fw_unmapped, fw_extra_copy, debug=False):
    """
        This function prints the report of the lifton results.

        Parameters:
        - ref_features_dict: reference features dictionary
        - fw_unmapped: file writer for unmapped features
        - fw_extra_copy: file writer for extra copy features
        - debug: debug mode
    """
    LIFTED_FEATURES = 0
    MISSED_FEATURES = 0
    LIFTED_SINGLE_CODING_FEATURES = 0
    LIFTED_SINGLE_NONCODING_FEATURES = 0
    LIFTED_EXTRA_CODING_FEATURES = 0
    LIFTED_EXTRA_NONCODING_FEATURES = 0
    LIFTED_EXTRA_CODING_SUM_FEATURES = 0
    LIFTED_EXTRA_NONCODING_SUM_FEATURES = 0
    for feature in ref_features_dict.keys():
        if feature == "LiftOn-gene":
            continue
        copy_num = ref_features_dict[feature].copy_num 
        if copy_num >= 1:
            LIFTED_FEATURES += 1
            if ref_features_dict[feature].is_protein_coding and copy_num == 1:
                LIFTED_SINGLE_CODING_FEATURES += 1
            elif ref_features_dict[feature].is_protein_coding and copy_num > 1:
                fw_extra_copy.write(f"{feature}\t{copy_num}\tcoding\n")
                LIFTED_EXTRA_CODING_FEATURES += 1
                LIFTED_EXTRA_CODING_SUM_FEATURES += (copy_num)
            elif not ref_features_dict[feature].is_protein_coding and copy_num == 1:
                LIFTED_SINGLE_NONCODING_FEATURES += 1
            elif not ref_features_dict[feature].is_protein_coding and copy_num > 1:
                fw_extra_copy.write(f"{feature}\t{copy_num}\tnon-coding\n")
                LIFTED_EXTRA_NONCODING_FEATURES += 1
                LIFTED_EXTRA_NONCODING_SUM_FEATURES += (copy_num)
        elif copy_num == 0:
            MISSED_FEATURES += 1
            fw_unmapped.write(f"{feature}\n")
    print("*********************************************", file=sys.stderr)
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
    print(f"*********************************************")
