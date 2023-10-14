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


def translate_CDS_2_protein(cds_children, fai, fai_protein, transcript_id):
    ################################
    # Step 2: Iterate through the children and chain the DNA sequence
    ################################
    trans_seq = ""
    cdss_lens = []
    for cds_idx, cds in enumerate(cds_children):
        # print(">> cds: ", cds)

        p_seq = cds.entry.sequence(fai)
        p_seq = Seq(p_seq)

        # Chaining the CDS features
        if cds.entry.strand == '-':
            trans_seq = p_seq + trans_seq
            # cdss_lens.append(cds.entry.end - cds.entry.start + 1)
            cdss_lens.insert(0, cds.entry.end - cds.entry.start + 1)
        elif cds.entry.strand == '+':
            trans_seq = trans_seq + p_seq
            cdss_lens.append(cds.entry.end - cds.entry.start + 1)

    ################################
    # Step 3: Translate the DNA sequence & get the reference protein sequence.
    ################################
    # ref_protein_seq = str(fai_protein[aa_trans_id])

    protein_seq = trans_seq.translate()

    peps = protein_seq.split("*")

    # for pep in peps:
    #     print(str(pep))

    if len(peps) == 2 and str(peps[1]) == "":
        # This is a valid protein
        return True
    else:
        print("Invalid protein: ", peps)
        find_orfs(trans_seq, fai_protein, transcript_id)
        return False

def find_orfs(sequence, fai_protein, transcript_id):
    # sequence = Seq.Seq("ATGTTGCTCTGATGAGGGGTGAAGCAAGGTTACAGTGAGAAGGGCCTGGAGGGAGGAGGTCCTGGAGGAGGGGGG")

    # Find ORFs manually
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]

    orf_list = []

    for frame in range(3):
        orf_idx_s = 0
        for i in range(frame, len(sequence), 3):
            codon = str(sequence[i:i+3])
            
            if codon == start_codon:
                orf_idx_s = i
                orf_idx_e = i
                orf_seq = ""
                for j in range(i, len(sequence), 3):
                    codon = str(sequence[j:j+3])
                    orf_seq += codon
                    if codon in stop_codons:
                        orf_idx_e = j+3
                        break
                if orf_seq:
                    orf = lifton_class.Lifton_ORF(orf_idx_s, orf_idx_e)
                    orf_list.append(orf)

    # Print the ORFs
    final_orf = ""
    max_align_score = 0
    for i, orf in enumerate(orf_list):
        orf_DNA_seq = sequence[orf.start:orf.end]
        orf_protein_seq = orf_DNA_seq.translate()
        
        print(f"\tORF {i+1}: {orf_DNA_seq}")
        print(f"\torf_protein_seq: {orf_protein_seq}")
        extracted_parasail_res, extracted_seq, reference_seq = align.parasail_align_base(orf_protein_seq, str(fai_protein[transcript_id]))

        alignment_score = extracted_parasail_res.score
        if alignment_score > max_align_score:
            max_align_score = alignment_score
            final_orf = orf_protein_seq
            print(f"\talignment_score: {alignment_score}")

    print(f">> max_align_score: {max_align_score}")
    print(f">> final_orf: {final_orf}")
    


