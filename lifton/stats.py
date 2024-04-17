import sys
from lifton import logger

def print_report(ref_features_dict, transcripts_stats_dict, fw_unmapped, fw_extra_copy, fw_mapped_feature, fw_mapped_trans, debug=False):
    """
    This function prints the report of the lifton results.

    Parameters:
    - ref_features_dict: reference features dictionary
    - fw_unmapped: file writer for unmapped features
    - fw_extra_copy: file writer for extra copy features
    - fw_mapped_feature: file writer for mapped features
    - fw_mapped_trans: file writer for mapped transcripts
    - debug: debug mode
    """
    LIFTED_FEATURES = 0
    MISSED_FEATURES = 0
    LIFTED_SINGLE_CODING_FEATURES = 0
    LIFTED_SINGLE_NONCODING_FEATURES = 0
    LIFTED_SINGLE_OTHER_FEATURES = 0
    LIFTED_EXTRA_CODING_FEATURES = 0
    LIFTED_EXTRA_NONCODING_FEATURES = 0
    LIFTED_EXTRA_OTHER_FEATURES = 0
    LIFTED_EXTRA_CODING_SUM_FEATURES = 0
    LIFTED_EXTRA_NONCODING_SUM_FEATURES = 0
    LIFTED_EXTRA_OTHER_SUM_FEATURES = 0
    
    # Feature stats (gene)
    for feature in ref_features_dict.keys():
        if feature == "LiftOn-gene":
            continue
        copy_num = ref_features_dict[feature].copy_num
        TYPE = ""
        if copy_num >= 1:
            LIFTED_FEATURES += 1
            if ref_features_dict[feature].is_protein_coding:
                TYPE = "coding"
                if copy_num == 1:
                    LIFTED_SINGLE_CODING_FEATURES += 1
                else:
                    LIFTED_EXTRA_CODING_FEATURES += 1
                    LIFTED_EXTRA_CODING_SUM_FEATURES += copy_num
            elif ref_features_dict[feature].is_non_coding:
                TYPE = "non-coding"
                if copy_num == 1:
                    LIFTED_SINGLE_NONCODING_FEATURES += 1
                else:
                    LIFTED_EXTRA_NONCODING_FEATURES += 1
                    LIFTED_EXTRA_NONCODING_SUM_FEATURES += copy_num
            else:
                TYPE = "other"
                if copy_num == 1:
                    LIFTED_SINGLE_OTHER_FEATURES += 1
                else:
                    LIFTED_EXTRA_OTHER_FEATURES += 1
                    LIFTED_EXTRA_OTHER_SUM_FEATURES += copy_num

            fw_mapped_feature.write(f"{feature}\t{copy_num}\t{TYPE}\n")
            if copy_num > 1:
                fw_extra_copy.write(f"{feature}\t{copy_num}\t{TYPE}\n")
        elif copy_num == 0:
            MISSED_FEATURES += 1
            if ref_features_dict[feature].is_protein_coding:
                TYPE = "coding"
            elif ref_features_dict[feature].is_non_coding:
                TYPE = "non-coding"
            else:
                TYPE = "other"
            fw_unmapped.write(f"{feature}\t{TYPE}\n")

    # Transcript stats        
    for TYPE, transs in transcripts_stats_dict.items():
        for trans, trans_copy_num in transs.items():
            fw_mapped_trans.write(f"{trans}\t{trans_copy_num}\t{TYPE}\n")

    print("\n\n*********************************************", file=sys.stderr)
    print(f"* Total features in reference\t\t: {len(ref_features_dict.keys())-1}", file=sys.stderr)
    print(f"* Lifted feature\t\t\t: {LIFTED_FEATURES} ({LIFTED_SINGLE_CODING_FEATURES + LIFTED_EXTRA_CODING_FEATURES} + {LIFTED_SINGLE_NONCODING_FEATURES + LIFTED_EXTRA_NONCODING_FEATURES} + {LIFTED_SINGLE_OTHER_FEATURES + LIFTED_EXTRA_OTHER_FEATURES})", file=sys.stderr)
    print(f"\t* Protein-coding feature\t: {LIFTED_SINGLE_CODING_FEATURES + LIFTED_EXTRA_CODING_FEATURES}", file=sys.stderr)
    print(f"\t* Non-coding feature\t\t: {LIFTED_SINGLE_NONCODING_FEATURES + LIFTED_EXTRA_NONCODING_FEATURES}", file=sys.stderr)
    print(f"\t* Other feature\t\t\t: {LIFTED_SINGLE_OTHER_FEATURES + LIFTED_EXTRA_OTHER_FEATURES}", file=sys.stderr)
    print(f"* Missed feature\t\t\t: {MISSED_FEATURES}\n", file=sys.stderr)

    print(f"* Total features in target\t\t: {LIFTED_SINGLE_CODING_FEATURES + LIFTED_EXTRA_CODING_SUM_FEATURES + LIFTED_SINGLE_NONCODING_FEATURES + LIFTED_EXTRA_NONCODING_SUM_FEATURES + LIFTED_SINGLE_OTHER_FEATURES + LIFTED_EXTRA_OTHER_SUM_FEATURES} ({LIFTED_SINGLE_CODING_FEATURES + LIFTED_EXTRA_CODING_SUM_FEATURES} + {LIFTED_SINGLE_NONCODING_FEATURES + LIFTED_EXTRA_NONCODING_SUM_FEATURES} + {LIFTED_SINGLE_OTHER_FEATURES + LIFTED_EXTRA_OTHER_SUM_FEATURES})", file=sys.stderr)
    print(f"\t* Protein-coding feature\t: {LIFTED_SINGLE_CODING_FEATURES + LIFTED_EXTRA_CODING_SUM_FEATURES} ({LIFTED_SINGLE_CODING_FEATURES} + {LIFTED_EXTRA_CODING_SUM_FEATURES})", file=sys.stderr)
    print(f"\t\t* single copy\t\t: {LIFTED_SINGLE_CODING_FEATURES}", file=sys.stderr)
    print(f"\t\t* > 1 copy\t\t: {LIFTED_EXTRA_CODING_FEATURES}, {LIFTED_EXTRA_CODING_SUM_FEATURES} in total", file=sys.stderr)
    print(f"\t* Non-coding feature\t\t: {LIFTED_SINGLE_NONCODING_FEATURES + LIFTED_EXTRA_NONCODING_SUM_FEATURES} ({LIFTED_SINGLE_NONCODING_FEATURES} + {LIFTED_EXTRA_NONCODING_SUM_FEATURES})", file=sys.stderr)
    print(f"\t\t* single copy\t\t: {LIFTED_SINGLE_NONCODING_FEATURES}", file=sys.stderr)
    print(f"\t\t* > 1 copy\t\t: {LIFTED_EXTRA_NONCODING_FEATURES}, {LIFTED_EXTRA_NONCODING_SUM_FEATURES} in total", file=sys.stderr)
    print(f"\t* Other feature\t\t\t: {LIFTED_SINGLE_OTHER_FEATURES + LIFTED_EXTRA_OTHER_SUM_FEATURES} ({LIFTED_SINGLE_OTHER_FEATURES} + {LIFTED_EXTRA_OTHER_SUM_FEATURES})", file=sys.stderr)
    print(f"\t\t* single copy\t\t: {LIFTED_SINGLE_OTHER_FEATURES}", file=sys.stderr)
    print(f"\t\t* > 1 copy\t\t: {LIFTED_EXTRA_OTHER_FEATURES}, {LIFTED_EXTRA_OTHER_SUM_FEATURES} in total", file=sys.stderr)
    print(f"*********************************************")
