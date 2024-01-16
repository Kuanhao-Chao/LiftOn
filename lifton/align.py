import parasail
from Bio.Seq import Seq
from cigar import Cigar
from lifton import get_id_fraction, lifton_class

# Change the CDS protein boundaries based on CIGAR string
def adjust_cdss_protein_boundary(cdss_protein_aln_boundary, cigar_accum_len, length):
    cds_boundary_shift = 0
    for i in range(len(cdss_protein_aln_boundary)):
        cdss_start = cdss_protein_aln_boundary[i][0] + cds_boundary_shift
        cdss_end = cdss_protein_aln_boundary[i][1] + cds_boundary_shift
        if (cdss_start <= cigar_accum_len) and (cdss_end >= cigar_accum_len):
            # print("### cdss_start: ", cdss_start)
            # print("### cigar_accum_len: ", cigar_accum_len)
            # print("### cdss_end: ", cdss_end)
            # print("### cigar_accum_len: ", cigar_accum_len)
            cds_boundary_shift += length
            cdss_end += length
        cdss_protein_aln_boundary[i] = (cdss_start, cdss_end)
    return cdss_protein_aln_boundary


# Map the CDSs boundaries on to the protein.
def get_cdss_protein_boundary(cdss_lens):
    cdss_cumulative = [sum(cdss_lens[:i+1]) for i in range(len(cdss_lens))]
    cdss_cumulative_div = [x / 3 for x in cdss_cumulative]
    cdss_protein_boundary = {}
    for idx in range(len(cdss_cumulative_div)):
        start = cdss_cumulative_div[idx-1] if idx > 0 else 0
        end = cdss_cumulative_div[idx]
        cdss_protein_boundary[idx] = (start, end)
    # print("cdss_cumulative : ", cdss_cumulative)
    # print("cdss_cumulative_div : ", cdss_cumulative_div)
    # print("cdss_protein_boundary: ", cdss_protein_boundary)
    # print("protein_seq len : ", len(protein_seq))
    return cdss_protein_boundary


def parasail_align_protein_base(protein_seq, ref_protein_seq):
    matrix = parasail.Matrix("blosum62")
    gap_open = 11
    gap_extend = 1
    extracted_seq = str(protein_seq)
    reference_seq = str(ref_protein_seq)
    extracted_seq = "*" if extracted_seq == "" else extracted_seq
    # (Query, Reference)
    extracted_parasail_res = parasail.nw_trace_scan_sat(extracted_seq, reference_seq, gap_open, gap_extend, matrix)
    return extracted_parasail_res, extracted_seq, reference_seq


def protein_align(protein_seq, ref_protein_seq):
    protein_seq = protein_seq.upper()
    extracted_parasail_res, extracted_seq, reference_seq = parasail_align_protein_base(protein_seq, str(ref_protein_seq))
    alignment_query = extracted_parasail_res.traceback.query
    alignment_comp = extracted_parasail_res.traceback.comp
    alignment_ref = extracted_parasail_res.traceback.ref
    # print("alignment_query: ", alignment_query)
    # print("alignment_comp : ", alignment_comp)
    # print("alignment_ref  : ", alignment_ref)
    extracted_matches, extracted_length = get_id_fraction.get_AA_id_fraction(extracted_parasail_res.traceback.ref, extracted_parasail_res.traceback.query)
    extracted_identity = extracted_matches/extracted_length
    lifton_aln = lifton_class.Lifton_Alignment(extracted_identity, None, alignment_query, alignment_comp, alignment_ref, None, None, extracted_seq, reference_seq, None)
    return lifton_aln


def parasail_align_DNA_base(trans_seq, ref_trans_seq):
    matrix = parasail.matrix_create("ACGT*", 1, -3)
    gap_open = 5
    gap_extend = 2
    extracted_seq = str(trans_seq)
    reference_seq = str(ref_trans_seq)
    # (Query, Reference)
    extracted_parasail_res = parasail.nw_trace_scan_sat(extracted_seq, reference_seq, gap_open, gap_extend, matrix)
    return extracted_parasail_res, extracted_seq, reference_seq


def trans_align(trans_seq, ref_trans_seq):
    trans_seq = trans_seq.upper()
    extracted_parasail_res, extracted_seq, reference_seq = parasail_align_DNA_base(trans_seq, str(ref_trans_seq))
    alignment_query = extracted_parasail_res.traceback.query
    alignment_comp = extracted_parasail_res.traceback.comp
    alignment_ref = extracted_parasail_res.traceback.ref
    # print("alignment_query: ", alignment_query) 
    # print("alignment_comp: ", alignment_comp) 
    # print("alignment_ref: ", alignment_ref) 
    extracted_matches, extracted_length = get_id_fraction.get_DNA_id_fraction(extracted_parasail_res.traceback.ref, extracted_parasail_res.traceback.query)
    extracted_identity = extracted_matches/extracted_length
    lifton_aln = lifton_class.Lifton_Alignment(extracted_identity, None, alignment_query, alignment_comp, alignment_ref, None, None, extracted_seq, reference_seq, None)
    return lifton_aln


def LiftOn_translate(tool, lifton_trans, db_entry, fai, ref_proteins, ref_trans_id):
    ref_protein_seq = str(ref_proteins[ref_trans_id])
    coding_seq, cds_children, cdss_lens = lifton_trans.get_coding_seq(fai)
    protein_seq = lifton_trans.translate_coding_seq(coding_seq)
    # print("protein_seq: ", protein_seq)
    return ref_protein_seq, protein_seq, cdss_lens, cds_children


def parasail_align(tool, lifton_trans, db_entry, fai, ref_proteins, ref_trans_id, lifton_status):
    # Step 1: Check if the ref_trans_id is in the reference protein dictionary
    aln = None
    if ref_trans_id not in ref_proteins.keys():
        return aln
    # Step 2: Translate the given annotation
    ref_protein_seq, protein_seq, cdss_lens, cds_children = LiftOn_translate(tool, lifton_trans, db_entry, fai, ref_proteins, ref_trans_id)
    # Step 3: Get the corresponding protein boundaries for each CDS
    cdss_protein_boundary = get_cdss_protein_boundary(cdss_lens)
    # Step 4: Protein alignment
    if protein_seq == None:
        return aln
    protein_seq = protein_seq.upper()
    extracted_parasail_res, extracted_seq, reference_seq = parasail_align_protein_base(protein_seq, ref_protein_seq)
    # Step 5: Extract the alignment information
    alignment_query = extracted_parasail_res.traceback.query
    alignment_comp = extracted_parasail_res.traceback.comp
    alignment_ref = extracted_parasail_res.traceback.ref
    cigar = extracted_parasail_res.cigar
    decoded_cigar = cigar.decode.decode()
    cigar_ls = list(Cigar(decoded_cigar).items())
    # Step 6: Change the CDS protein boundaries based on CIGAR string.
    cigar_accum_len = 0
    cdss_protein_aln_boundary = cdss_protein_boundary.copy()
    for length, symbol in cigar_ls:
        if symbol == "D":
            # print(length, symbol)
            cdss_protein_aln_boundary = adjust_cdss_protein_boundary(cdss_protein_aln_boundary, cigar_accum_len, length)
        cigar_accum_len += length
        # print("cigar_accum_len: ", cigar_accum_len)
    # print("cdss_protein_boundary    : ", cdss_protein_boundary)
    # print("cdss_protein_aln_boundary: ", cdss_protein_aln_boundary)
    # print("\t>> alignment_query: ", len(alignment_query))
    # print("\t>> alignment_comp: ", len(alignment_comp))
    # print("\t>> alignment_ref : ", len(alignment_ref))
    extracted_matches, extracted_length = get_id_fraction.get_AA_id_fraction(extracted_parasail_res.traceback.ref, extracted_parasail_res.traceback.query)
    extracted_identity = extracted_matches/extracted_length
    aln = lifton_class.Lifton_Alignment(extracted_identity, cds_children, alignment_query, alignment_comp, alignment_ref, cdss_protein_boundary, cdss_protein_aln_boundary, extracted_seq, reference_seq, db_entry)
    if tool == "liftoff":
        # SETTING Liftoff identity score
        lifton_status.liftoff = aln.identity
    elif tool == "miniprot":
        # SETTING miniprot identity score
        lifton_status.miniprot = aln.identity
    elif tool == "lifton":
        # SETTING miniprot identity score
        lifton_status.lifton_aa = aln.identity
    return aln
