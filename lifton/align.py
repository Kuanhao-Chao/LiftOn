import parasail
from Bio.Seq import Seq
from cigar import Cigar
from lifton import adjust_cds_boundaries

def get_id_fraction(reference, target):
    matches = 0
    for i, letter in enumerate(reference):
        if letter == target[i]:
            matches += 1
    return matches, max(len(reference), len(target))

def get_cdss_protein_boundary(cdss_len):
    cdss_cumulative = [sum(cdss_len[:i+1]) for i in range(len(cdss_len))]
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


def parasail_align(tool, db, db_entry, fai, fai_protein, aa_trans_id):
    # Get the children of the entry
    # cds_children = [child for child in db.children(db_entry, featuretype='CDS')]
    # print("Before cds_children: ", cds_children)
    cds_children = []
    for child in db.children(db_entry, featuretype='CDS'):
        if len(cds_children) == 0:
            cds_children.append(child)
            continue
        idx_insert = 0
        for idx_c in range(len(cds_children)):
            itr_c = cds_children[idx_c]
            if child.start > itr_c.end:
                idx_insert += 1
        
        cds_children.insert(idx_insert, child)

    # print("After cds_children: ", cds_children)
    # Iterate through the children and print their attributes
    trans_seq = ""
    cdss_len = []
    for cds_idx, cds in enumerate(cds_children):

        # Include the stop coding for the last CDS(+) / first CDS(-) for miniprot 
        if tool == "miniprot" and cds_idx == 0 and cds.strand == '-':
            cds.start = cds.start -3
        if tool == "miniprot" and cds_idx == len(cds_children)-1 and cds.strand == '+':
            cds.end = cds.end + 3

        p_seq = cds.sequence(fai)
        p_seq = Seq(p_seq)

        # Chaining the CDS features
        if cds.strand == '-':
            trans_seq = p_seq + trans_seq
            # cdss_len.append(cds.end - cds.start + 1)
            cdss_len.insert(0, cds.end - cds.start + 1)
        elif cds.strand == '+':
            trans_seq = trans_seq + p_seq
            cdss_len.append(cds.end - cds.start + 1)
        # print('>' + cds.id + '\n' + p_seq)
    # print(cds.strand+' trans_seq: ' + trans_seq)


    protein_seq = str(trans_seq.translate())
    ref_protein_seq = str(fai_protein[aa_trans_id])

    cdss_protein_boundary = get_cdss_protein_boundary(cdss_len)
    matrix = parasail.Matrix("blosum62")
    gap_open = 11
    gap_extend = 1

    extracted_seq = str(protein_seq)
    reference_seq = str(ref_protein_seq) + "*"
    # print("db_entry: ", db_entry)
    # print("\treference_seq: ", reference_seq)
    # print("\textracted_seq: ", extracted_seq)
    
    # (Query, Reference)
    extracted_parasail_res = parasail.nw_trace_scan_sat(extracted_seq, reference_seq, gap_open, gap_extend, matrix)

    # Extract the alignment information
    alignment_score = extracted_parasail_res.score
    alignment_query = extracted_parasail_res.traceback.query
    alignment_comp = extracted_parasail_res.traceback.comp
    alignment_ref = extracted_parasail_res.traceback.ref

    cigar = extracted_parasail_res.cigar
    # alignment_start_query = extracted_parasail_res.traceback.query_begin
    # alignment_end_query = extracted_parasail_res.traceback.query_end
    # alignment_start_comp = extracted_parasail_res.traceback.comp_begin
    # alignment_end_comp = extracted_parasail_res.traceback.comp_end

    # print("extracted_seq: ", extracted_seq)
    # # Print the additional alignment results
    # print("\t>> alignment_score: ", alignment_score)
    # print("\t>> cigar.seq   : ", cigar.seq)
    # for i in cigar.seq:
    #     print("\t\t>> cigar.i   : ", i)
    # use decode attribute to return a decoded cigar string
    print("\t>> cigar.decode: ", cigar.decode.decode())
    decoded_cigar = cigar.decode.decode()
    cigar_ls = list(Cigar(decoded_cigar).items())
    print("\t>> cigar_ls: ", cigar_ls)

    
    cigar_accum_len = 0
    cdss_protein_aln_boundary = cdss_protein_boundary.copy()
    for length, symbol in cigar_ls:
        if symbol == "D":
            # print(length, symbol)
            cdss_protein_aln_boundary = adjust_cds_boundaries.adjust_cdss_protein_boundary(cdss_protein_aln_boundary, cigar_accum_len, length)
        # else:
        #     print("cigar_accum_len: ", cigar_accum_len)
        cigar_accum_len += length
        # print("cigar_accum_len: ", cigar_accum_len)

    print("cdss_protein_boundary    : ", cdss_protein_boundary)
    print("cdss_protein_aln_boundary: ", cdss_protein_aln_boundary)
    print("\t>> alignment_query: ", len(alignment_query))
    print("\t>> alignment_comp: ", len(alignment_comp))
    print("\t>> alignment_ref : ", len(alignment_ref))
    
    extracted_matches, extracted_length = get_id_fraction(extracted_parasail_res.traceback.ref, extracted_parasail_res.traceback.query)
    extracted_identity = extracted_matches/extracted_length

    return extracted_identity, cds_children, alignment_query, alignment_comp, alignment_ref, cdss_protein_aln_boundary, extracted_seq, reference_seq
