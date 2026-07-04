import os
import sys
from collections import defaultdict
from lifton import logger


class _Tee:
    """Write the same text to several streams at once (GitHub issue #50: the
    run summary goes to BOTH stderr and stats/summary.txt)."""
    def __init__(self, *streams):
        self._streams = streams

    def write(self, data):
        for s in self._streams:
            s.write(data)

    def flush(self):
        for s in self._streams:
            try:
                s.flush()
            except Exception:
                pass

def print_report(ref_features_dict, transcripts_stats_dict, fw_unmapped, fw_extra_copy, fw_mapped_feature, fw_mapped_trans, debug=False, fw_feature_type=None):
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
    # per raw-GFF-feature-type completeness tally (all-feature tracking)
    by_type = defaultdict(lambda: {"reference": 0, "lifted": 0, "missed": 0,
                                   "extra": 0, "target": 0})

    # Feature stats (gene)
    for feature in ref_features_dict.keys():
        if feature == "LiftOn-gene":
            continue
        copy_num = ref_features_dict[feature].copy_num
        TYPE = ""
        ftype = (getattr(ref_features_dict[feature], "feature_type", None)
                 or getattr(ref_features_dict[feature], "biotype", None) or "unknown")
        bucket = by_type[ftype]
        bucket["reference"] += 1
        if copy_num >= 1:
            bucket["lifted"] += 1
            bucket["target"] += copy_num
            if copy_num > 1:
                bucket["extra"] += 1
        else:
            bucket["missed"] += 1
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

    # All-feature-type completeness side file (additive; never affects the GFF3)
    if fw_feature_type is not None:
        fw_feature_type.write("feature_type\tn_reference\tn_lifted\tn_missed\t"
                              "n_extra_copies\tn_target\tpct_recovered\n")
        for ftype in sorted(by_type):
            b = by_type[ftype]
            pct = (b["lifted"] / b["reference"]) if b["reference"] else 0.0
            fw_feature_type.write(f"{ftype}\t{b['reference']}\t{b['lifted']}\t"
                                  f"{b['missed']}\t{b['extra']}\t{b['target']}\t{pct:.5f}\n")

    # GH #50: persist the summary block (below) to stats/summary.txt as well as
    # stderr. Derive the stats dir from a side-car writer's path so no new argument
    # or caller change is needed; tee stderr for the duration of the block.
    _summary_fh = None
    _stats_dir = os.path.dirname(getattr(fw_mapped_feature, "name", "") or "")
    if _stats_dir:
        try:
            _summary_fh = open(os.path.join(_stats_dir, "summary.txt"), "w")
        except OSError:
            _summary_fh = None
    _prev_stderr = sys.stderr
    if _summary_fh is not None:
        sys.stderr = _Tee(_prev_stderr, _summary_fh)

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
    print("* Completeness by feature type:", file=sys.stderr)
    for ftype in sorted(by_type):
        b = by_type[ftype]
        pct = (100.0 * b["lifted"] / b["reference"]) if b["reference"] else 0.0
        print(f"\t* {ftype}\t: {b['lifted']}/{b['reference']} ({pct:.1f}%)", file=sys.stderr)
    print(f"*********************************************")

    # GH #50: restore stderr and finalise summary.txt (add the trailing rule, which
    # the line above prints to stdout, so the file is self-contained).
    sys.stderr = _prev_stderr
    if _summary_fh is not None:
        _summary_fh.write("*********************************************\n")
        _summary_fh.close()
