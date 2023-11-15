from lifton import get_id_fraction, lifton_entry
import math

def get_protein_boundary(cdss_aln_boundary, c_idx_last, c_idx):
    # for mi in range(m_c_idx_last, m_c_idx):
    #     # print("\t m_i: ", mi)
    #     # m_aa_start, m_aa_end = get_protein_boundary(mi, m_protein)
        
    #     m_c = m_children[mi]
    #     print(f"m_cdss_aln_boundary[{mi}]:  {m_cdss_aln_boundary[mi]}")

    aa_start = cdss_aln_boundary[c_idx_last][0]
    aa_end = cdss_aln_boundary[c_idx-1][1]
    # print(f"\t### protein chunk boundaries: {aa_start} - {aa_end}")
    
    return aa_start, aa_end

def get_protein_reference_length_single(lifton_aln, c_idx):
    # for mi in range(m_c_idx_last, m_c_idx):
    #     # print("\t m_i: ", mi)
    #     # m_aa_start, m_aa_end = get_protein_boundary(mi, m_protein)
        
    #     m_c = m_children[mi]
    #     print(f"m_cdss_aln_boundary[{mi}]:  {m_cdss_aln_boundary[mi]}")

    # aa_start = lifton_aln.cdss_protein_aln_boundaries[c_idx][0]
    aa_start = 0
    aa_end = lifton_aln.cdss_protein_aln_boundaries[c_idx][1]
    # print(f"\t### protein chunk boundaries: {aa_start} - {aa_end}")

    # aa_start = math.floor(aa_start)
    aa_end = math.ceil(aa_end)


    ref_count = 0
    # print("\tlifton_aln.ref_seq[aa_start:aa_start]: ", lifton_aln.ref_seq[aa_start:aa_end])
    for i, letter in enumerate(lifton_aln.ref_seq[aa_start:aa_end]):
        if letter != "-":
            ref_count += 1
    return ref_count

def push_cds_idx(c_idx, lifton_aln, ref_aa_count):
    # print("\t>> c_idx: ", c_idx)
    ref_aa_count = get_protein_reference_length_single(lifton_aln, c_idx)
    c_idx += 1
    # ref_aa_count += ref_aa_lcl_count
    return c_idx, ref_aa_count


def process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, fw):
    # print(">> process_m_l_children!! ")
    # print("\tl_c_idx: ", l_c_idx)
    # print("\tm_c_idx: ", m_c_idx)
    # print("\tl_c_idx_last: ", l_c_idx_last)
    # print("\tm_c_idx_last: ", m_c_idx_last)

    # print("## m_cdss_aln_boundary: ", m_lifton_aln.cdss_protein_aln_boundaries)
    # print("## l_cdss_aln_boundary: ", l_lifton_aln.cdss_protein_aln_boundaries)

    m_aa_start, m_aa_end = get_protein_boundary(m_lifton_aln.cdss_protein_aln_boundaries, m_c_idx_last, m_c_idx)
    l_aa_start, l_aa_end = get_protein_boundary(l_lifton_aln.cdss_protein_aln_boundaries, l_c_idx_last, l_c_idx)

    # print(f"\tminiprot: {m_aa_start} - {m_aa_end}")
    # print(f"\tliftoff : {l_aa_start} - {l_aa_end}")

    m_matches, m_length = get_id_fraction.get_id_fraction(m_lifton_aln.ref_aln, m_lifton_aln.query_aln, math.floor(m_aa_start), math.ceil(m_aa_end))
    l_matches, l_length = get_id_fraction.get_id_fraction(l_lifton_aln.ref_aln, l_lifton_aln.query_aln, math.floor(l_aa_start), math.ceil(l_aa_end))

    # print(f"\tm_matches: {m_matches};  m_length: {m_length}")
    # print(f"\tl_matches: {l_matches};  l_length: {l_length}")

    if (m_matches/m_length <= l_matches/l_length):
        # print("\t# => Liftoff is doing better")
        cds_ls = lifton_entry.create_lifton_entries(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, fw, False)
    elif (m_matches/m_length > l_matches/l_length):
        # print("\t# => miniprot is doing better")
        cds_ls = lifton_entry.create_lifton_entries(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, fw, True)

    # print(">> New iteration!\n\n ")

    return cds_ls


def chaining_algorithm(l_lifton_aln, m_lifton_aln, fai, fw):
    l_children = l_lifton_aln.cds_children
    m_children = m_lifton_aln.cds_children

    # print("number of children, m: ", len(m_children))
    # print("number of children, l: ", len(l_children))
    
    m_c_idx = 0
    l_c_idx = 0
    m_c_idx_last = m_c_idx
    l_c_idx_last = l_c_idx

    cds_list = []
    ref_aa_liftoff_count = 0
    ref_aa_miniprot_count = 0

    while m_c_idx != (len(m_children)-1) or l_c_idx != (len(l_children)-1):
        
        if (m_c_idx == len(m_children)-1) and (l_c_idx < (len(l_children)-1)):
            l_c_idx, ref_aa_liftoff_count = push_cds_idx(l_c_idx, l_lifton_aln, ref_aa_liftoff_count)
            continue
        if (m_c_idx < len(m_children)-1) and (l_c_idx == (len(l_children)-1)):
            m_c_idx, ref_aa_miniprot_count = push_cds_idx(m_c_idx, m_lifton_aln, ref_aa_miniprot_count)
            continue

        # print("m_lifton_aln.db_entry.strand: ", m_lifton_aln.db_entry.strand)
        # ref_aa_miniprot_lcl_count = get_protein_reference_length_single(m_lifton_aln, m_c_idx)
        # ref_aa_miniprot_count += ref_aa_miniprot_lcl_count
        # ref_aa_liftoff_lcl_count = get_protein_reference_length_single(l_lifton_aln, l_c_idx)
        # ref_aa_liftoff_count += ref_aa_liftoff_lcl_count

        # print(">> ref_aa_miniprot_count: ", ref_aa_miniprot_count)
        # print(">> ref_aa_liftoff_count: ", ref_aa_liftoff_count)

        if m_lifton_aln.db_entry.strand == "+":
            # print(">>> m_c_idx: ", m_c_idx)
            # print(">>> l_c_idx: ", l_c_idx)
            m_c = m_children[m_c_idx]
            l_c = l_children[l_c_idx]

            # print(f"miniprot : {m_c.start} - {m_c.end}")
            # print(f"liftoff  : {l_c.start} - {l_c.end}")

            if ref_aa_liftoff_count < ref_aa_miniprot_count:
                l_c_idx, ref_aa_liftoff_count = push_cds_idx(l_c_idx, l_lifton_aln, ref_aa_liftoff_count)
            elif ref_aa_liftoff_count > ref_aa_miniprot_count:
                m_c_idx, ref_aa_miniprot_count = push_cds_idx(m_c_idx, m_lifton_aln, ref_aa_miniprot_count)
            elif ref_aa_liftoff_count == ref_aa_miniprot_count:
                # 1. protein references are the same
                # 2. CDS ends at the same position.
                if l_c_idx > 0 and m_c_idx > 0 and m_c.end == l_c.end:
                    cdss = process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, fw)
                    cds_list += cdss
                    m_c_idx_last = m_c_idx
                    l_c_idx_last = l_c_idx

                l_c_idx, ref_aa_liftoff_count = push_cds_idx(l_c_idx, l_lifton_aln, ref_aa_liftoff_count)
                m_c_idx, ref_aa_miniprot_count = push_cds_idx(m_c_idx, m_lifton_aln, ref_aa_miniprot_count)


        elif m_lifton_aln.db_entry.strand == "-":
            m_c = m_children[len(m_children) - m_c_idx - 1]
            l_c = l_children[len(l_children) - l_c_idx - 1]

            # print(f"miniprot : {m_c.start} - {m_c.end}")
            # print(f"liftoff  : {l_c.start} - {l_c.end}")
            
            if ref_aa_liftoff_count < ref_aa_miniprot_count:
                # print(">> Liftoff")
                l_c_idx, ref_aa_liftoff_count = push_cds_idx(l_c_idx, l_lifton_aln, ref_aa_liftoff_count)
            elif ref_aa_liftoff_count > ref_aa_miniprot_count:
                # print(">> miniprot")
                m_c_idx, ref_aa_miniprot_count = push_cds_idx(m_c_idx, m_lifton_aln, ref_aa_miniprot_count)
            elif ref_aa_liftoff_count == ref_aa_miniprot_count:
                # 1. protein references are the same
                # 2. CDS ends at the same position.
                # print(">> Finally comparison now!")
                # print(">> miniprot")
                if l_c_idx > 0 and m_c_idx > 0 and m_c.end == l_c.end:
                    cdss = process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, fw)
                    cds_list += cdss
                    m_c_idx_last = m_c_idx
                    l_c_idx_last = l_c_idx

                l_c_idx, ref_aa_liftoff_count = push_cds_idx(l_c_idx, l_lifton_aln, ref_aa_liftoff_count)
                m_c_idx, ref_aa_miniprot_count = push_cds_idx(m_c_idx, m_lifton_aln, ref_aa_miniprot_count)
        
    l_c_idx += 1
    m_c_idx += 1
    cdss = process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, fw)
    cds_list += cdss

    return cds_list