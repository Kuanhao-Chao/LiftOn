import re
from Bio.Seq import Seq
from lifton import align, lifton_class

def segments_overlap(segment1, segment2):
    # Check if the segments have valid endpoints
    # print("Checking two segments overlapping.!")

    # print(segment1, segment2)

    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    
    # Sort the segments by their left endpoints
    segment1, segment2 = sorted([segment1, segment2], key=lambda x: x[0])


    # Check if the right endpoint of the first segment is greater than or equal to the left endpoint of the second segment
    # print(segment1[1] >= segment2[0])

    return segment1[1] >= segment2[0]


def custom_bisect_insert(sorted_list, element_to_insert):
    low = 0
    high = len(sorted_list)

    while low < high:
        mid = (low + high) // 2
        if sorted_list[mid].entry.end < element_to_insert.entry.end:
            low = mid + 1
        else:
            high = mid

    sorted_list.insert(low, element_to_insert)

def get_ID_base(id):

    id_base = id.split("_")[0]
    return id_base

def get_trans_ID_base(id):

    # Regular expression pattern to match the desired substrings
    pattern = r'[A-Za-z0-9_]+-([A-Za-z0-9_]+_\d+\.\d+)'

    match = re.search(pattern, id)
    id_base = ""
    if match:
        id_base = match.group(0)  # Full match

    return id_base

def get_parent_features_to_lift(feature_types_file):
    feature_types = ["gene"]
    if feature_types_file is not None:
        f = open(feature_types_file)
        for line in f.readlines():
            feature_types.append(line.rstrip())
    return feature_types

def update_copy(id_base, copy_num_dict):
    if id_base in copy_num_dict.keys():
        copy_num_dict[id_base] += 1
    else:
        copy_num_dict[id_base] = 0

def LiftOn_no_miniprot(lifton_gene, transcript_id, fai, fai_protein, lifon_status, outdir, LIFTON_BAD_PROT_TRANS_COUNT):
    on_lifton_aln, good_trans = lifton_gene.fix_truncated_protein(transcript_id, fai, fai_protein)
    # SETTING LiftOn score
    lifon_status.lifton = max(lifon_status.lifton, on_lifton_aln.identity)
    lifon_status.status = "LiftOff_truncated_protein"

    if on_lifton_aln.identity < 1:
        # Writing out truncated miniprot annotation
        on_lifton_aln.write_alignment(outdir, "lifton", transcript_id)
    if not good_trans:
        LIFTON_BAD_PROT_TRANS_COUNT += 1

    return LIFTON_BAD_PROT_TRANS_COUNT
    