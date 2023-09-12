
def get_protein_boundary():
    pass

def process_m_l_children(m_c_idx, m_c_idx_last, m_children, m_aln_query, m_aln_comp, m_aln_ref, m_protein, l_c_idx, l_c_idx_last, l_children, l_aln_query, l_aln_comp, l_aln_ref, l_protein, fai):
    print("l_c_idx: ", l_c_idx)
    print("m_c_idx: ", m_c_idx)
    print("l_c_idx_last: ", l_c_idx_last)
    print("m_c_idx_last: ", m_c_idx_last)

    m_trans_seq = ""


    m_seq_start = 0
    m_seq_end = 0

    m_group_len = 0
    for mi in range(m_c_idx_last, m_c_idx):
        print("\t m_i: ", mi)
        # m_aa_start, m_aa_end = get_protein_boundary(mi, m_protein)
        
        m_c = m_children[mi]
        m_group_len += (m_c.end - m_c.start + 1)

    m_seq_end = m_seq_start + m_group_len

    #     m_seq = m_c.sequence(fai)
    #     m_seq = Seq(m_seq)

    #     # Chaining the CDS features
    #     if m_c.strand == '-':
    #         m_trans_seq = m_seq + m_trans_seq
    #     elif m_c.strand == '+':
    #         m_trans_seq = m_trans_seq + m_seq
    # m_protein_seq = str(m_trans_seq.translate())
    # print("m_protein_seq: ", m_protein_seq)


    for li in range(l_c_idx_last, l_c_idx):
        print("\t l_i: ", li)
        l_c = l_children[li]

    pass

def fix_transcript_annotation(m_children, m_aln_query, m_aln_comp, m_aln_ref, m_protein, l_children, l_aln_query, l_aln_comp, l_aln_ref, l_protein, fai):
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
            process_m_l_children(m_c_idx, m_c_idx_last, m_children, m_aln_query, m_aln_comp, m_aln_ref, m_protein, l_c_idx, l_c_idx_last, l_children, l_aln_query, l_aln_comp, l_aln_ref, l_protein, fai)

            m_c_idx_last = m_c_idx
            l_c_idx_last = l_c_idx
            print("\tGroup!")
        
    l_c_idx += 1
    m_c_idx += 1
    process_m_l_children(m_c_idx, m_c_idx_last, m_children, m_aln_query, m_aln_comp, m_aln_ref, m_protein, l_c_idx, l_c_idx_last, l_children, l_aln_query, l_aln_comp, l_aln_ref, l_protein, fai)
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


