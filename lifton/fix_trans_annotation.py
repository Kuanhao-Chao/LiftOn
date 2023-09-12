from lifton import get_id_fraction, lifton_class
import math

def get_protein_boundary(cdss_aln_boundary, c_idx_last, c_idx):
    # for mi in range(m_c_idx_last, m_c_idx):
    #     # print("\t m_i: ", mi)
    #     # m_aa_start, m_aa_end = get_protein_boundary(mi, m_protein)
        
    #     m_c = m_children[mi]
    #     print(f"m_cdss_aln_boundary[{mi}]:  {m_cdss_aln_boundary[mi]}")
    aa_start = cdss_aln_boundary[c_idx_last][0]
    aa_end = cdss_aln_boundary[c_idx-1][1]
    print(f"### protein chunk boundaries: {aa_start} - {aa_end}")
    return aa_start, aa_end


def process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai):
    print("l_c_idx: ", l_c_idx)
    print("m_c_idx: ", m_c_idx)
    print("l_c_idx_last: ", l_c_idx_last)
    print("m_c_idx_last: ", m_c_idx_last)

    # print("## m_cdss_aln_boundary: ", m_cdss_aln_boundary)
    # print("## l_cdss_aln_boundary: ", l_cdss_aln_boundary)

    m_aa_start, m_aa_end = get_protein_boundary(m_lifton_aln.cdss_protein_aln_boundaries, m_c_idx_last, m_c_idx)
    l_aa_start, l_aa_end = get_protein_boundary(l_lifton_aln.cdss_protein_aln_boundaries, l_c_idx_last, l_c_idx)


        # self.identity = extracted_identity
        # self.cds_children = cds_children
        # self.query = alignment_query
        # self.comp = alignment_comp
        # self.ref = alignment_ref
        # self.cdss_protein_boundaries = cdss_protein_boundary
        # self.cdss_protein_aln_boundaries = cdss_protein_aln_boundary
        # self.query_seq = extracted_seq
        # self.ref_seq = reference_seq
    
    m_matches, m_length = get_id_fraction.get_id_fraction(m_lifton_aln.ref_aln, m_lifton_aln.query_aln, math.floor(m_aa_start), math.ceil(m_aa_end))
    l_matches, l_length = get_id_fraction.get_id_fraction(l_lifton_aln.ref_aln, l_lifton_aln.query_aln, math.floor(l_aa_start), math.ceil(l_aa_end))

    print(f"m_matches: {m_matches};  m_length: {m_length}")
    print(f"l_matches: {l_matches};  l_length: {l_length}")

    if (m_matches/m_length <= l_matches/l_length):
        print("Liftoff is doing better")
    elif (m_matches/m_length > l_matches/l_length):
        print("miniprot is doing better")


def fix_transcript_annotation(m_lifton_aln, l_lifton_aln, fai):

    m_children = m_lifton_aln.cds_children
    l_children = l_lifton_aln.cds_children
    print("number of children, m: ", len(m_children))
    print("number of children, l: ", len(l_children))

    m_len_sum = 0
    
    m_c_idx = 0
    l_c_idx = 0
    m_c_idx_last = m_c_idx
    l_c_idx_last = l_c_idx
    while m_c_idx != (len(m_children)-1) or l_c_idx != (len(l_children)-1):
        if (m_c_idx == len(m_children)-1) and (l_c_idx < (len(l_children)-1)):
            l_c_idx += 1
            continue
        if (m_c_idx < len(m_children)-1) and (l_c_idx == (len(l_children)-1)):
            m_c_idx += 1
            continue

        # print(">>> m_c_idx: ", m_c_idx)
        # print(">>> l_c_idx: ", l_c_idx)
        m_c = m_children[m_c_idx]
        l_c = l_children[l_c_idx]

        if m_c.end > l_c.end:
            l_c_idx += 1
        elif m_c.end < l_c.end:
            m_c_idx += 1
        elif m_c.end == l_c.end:
            l_c_idx += 1
            m_c_idx += 1

            process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai)

            m_c_idx_last = m_c_idx
            l_c_idx_last = l_c_idx
            print("\tGroup!")
        
    l_c_idx += 1
    m_c_idx += 1
    process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai)
    print("\tGroup!")


    # for child in m_children:
    #     child_len = child.end - child.start + 1
    #     m_len_sum += child_len
    #     print("start: ", child.start, ";  start: ", child.end, ";  len: ", child_len, ";  len_sum: ", m_len_sum, ";  aa_len: ", m_len_sum/3)
    # print("m_aln_query: ", len(m_aln_query))
    
    # print("number of children: ", len(l_children))
    # l_len_sum = 0
    # for child in l_children:
    #     child_len = child.end - child.start + 1
    #     l_len_sum += child_len
    #     print("start: ", child.start, ";  start: ", child.end, ";  len: ", child_len, ";  len_sum: ", l_len_sum, ";  aa_len: ", l_len_sum/3)
    # print("l_aln_query: ", len(l_aln_query))


