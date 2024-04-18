from lifton import get_id_fraction, lifton_class, logger
import math

def get_protein_boundary(cdss_aln_boundary, c_idx_last, c_idx, DEBUG):
    aa_start = cdss_aln_boundary[c_idx_last][0]
    aa_end = cdss_aln_boundary[c_idx-1][1]
    return aa_start, aa_end


def get_protein_reference_length_single(lifton_aln, c_idx, DEBUG):
    aa_start = 0
    aa_end = lifton_aln.cdss_protein_aln_boundaries[c_idx][1]
    aa_end = math.ceil(aa_end)
    ref_count = 0
    # logger.log(f"\t### protein chunk boundaries: {aa_start} - {aa_end}", debug=DEBUG)
    # logger.log("\tlifton_aln.ref_seq[aa_start:aa_start]: ", lifton_aln.ref_seq[aa_start:aa_end], debug=DEBUG)
    for i, letter in enumerate(lifton_aln.ref_seq[aa_start:aa_end]):
        if letter != "-":
            ref_count += 1
    return ref_count


def push_cds_idx(c_idx, lifton_aln, ref_aa_count, DEBUG):
    ref_aa_count = get_protein_reference_length_single(lifton_aln, c_idx, DEBUG)
    c_idx += 1
    return c_idx, ref_aa_count


def process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, chains, DEBUG):
    # logger.log(">> process_m_l_children!! ", debug=DEBUG)
    # logger.log("\tl_c_idx: ", l_c_idx, debug=DEBUG)
    # logger.log("\tm_c_idx: ", m_c_idx, debug=DEBUG)
    # logger.log("\tl_c_idx_last: ", l_c_idx_last, debug=DEBUG)
    # logger.log("\tm_c_idx_last: ", m_c_idx_last, debug=DEBUG)
    # logger.log("## m_cdss_aln_boundary: ", m_lifton_aln.cdss_protein_aln_boundaries, debug=DEBUG)
    # logger.log("## l_cdss_aln_boundary: ", l_lifton_aln.cdss_protein_aln_boundaries, debug=DEBUG)
    m_aa_start, m_aa_end = get_protein_boundary(m_lifton_aln.cdss_protein_aln_boundaries, m_c_idx_last, m_c_idx, DEBUG)
    l_aa_start, l_aa_end = get_protein_boundary(l_lifton_aln.cdss_protein_aln_boundaries, l_c_idx_last, l_c_idx, DEBUG)
    m_matches, m_length = get_id_fraction.get_partial_id_fraction(m_lifton_aln.ref_aln, m_lifton_aln.query_aln, math.floor(m_aa_start), math.ceil(m_aa_end))
    l_matches, l_length = get_id_fraction.get_partial_id_fraction(l_lifton_aln.ref_aln, l_lifton_aln.query_aln, math.floor(l_aa_start), math.ceil(l_aa_end))
    # logger.log(f"\tminiprot: {m_aa_start} - {m_aa_end}", debug=DEBUG)
    # logger.log(f"\tliftoff : {l_aa_start} - {l_aa_end}", debug=DEBUG)
    # logger.log(f"\tm_matches: {m_matches};  m_length: {m_length}", debug=DEBUG)
    # logger.log(f"\tl_matches: {l_matches};  l_length: {l_length}", debug=DEBUG)
    if (m_matches/m_length <= l_matches/l_length):
        cds_ls = create_lifton_entries(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, False)
        chains.append(f"liftoff[{l_aa_start}-{l_aa_end}]")
    elif (m_matches/m_length > l_matches/l_length):
        cds_ls = create_lifton_entries(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, True)
        chains.append(f"miniprot[{m_aa_start}-{m_aa_end}]")
    return cds_ls


def create_lifton_entries(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, miniprot_is_better):
    cds_list = []
    source_aln = m_lifton_aln if miniprot_is_better else l_lifton_aln
    for c_idx in range(m_c_idx_last, m_c_idx) if miniprot_is_better else range(l_c_idx_last, l_c_idx):
        if source_aln.db_entry.strand == "+":
            c_idx_fix = c_idx
        elif source_aln.db_entry.strand == "-":
            c_idx_fix = len(source_aln.cds_children) - c_idx - 1
        lifton_cds = source_aln.cds_children[c_idx_fix]
        lifton_cds.attributes = l_lifton_aln.cds_children[0].attributes
        cds_list.append(lifton_class.Lifton_CDS(lifton_cds))
    return cds_list


def chaining_algorithm(l_lifton_aln, m_lifton_aln, fai, DEBUG):
    l_children = l_lifton_aln.cds_children
    m_children = m_lifton_aln.cds_children
    m_c_idx = 0
    l_c_idx = 0
    m_c_idx_last = m_c_idx
    l_c_idx_last = l_c_idx
    cds_list = []
    ref_aa_liftoff_count = 0
    ref_aa_miniprot_count = 0
    chains = []
    while m_c_idx != (len(m_children)-1) or l_c_idx != (len(l_children)-1):
        if (m_c_idx == len(m_children)-1) and (l_c_idx < (len(l_children)-1)):
            l_c_idx, ref_aa_liftoff_count = push_cds_idx(l_c_idx, l_lifton_aln, ref_aa_liftoff_count, DEBUG)
            continue
        if (m_c_idx < len(m_children)-1) and (l_c_idx == (len(l_children)-1)):
            m_c_idx, ref_aa_miniprot_count = push_cds_idx(m_c_idx, m_lifton_aln, ref_aa_miniprot_count, DEBUG)
            continue
        # logger.log(f">> ref_aa_miniprot_count: {ref_aa_miniprot_count}", debug=DEBUG)
        # logger.log(f">> ref_aa_liftoff_count : {ref_aa_liftoff_count}", debug=DEBUG)
        # logger.log(f">> len(m_children)      : {len(m_children)}", debug=DEBUG)
        # logger.log(f">> len(l_children)      : {len(l_children)}", debug=DEBUG)
        if m_lifton_aln.db_entry.strand == "+":
            m_c = m_children[m_c_idx]
            l_c = l_children[l_c_idx]
            if ref_aa_liftoff_count < ref_aa_miniprot_count:
                l_c_idx, ref_aa_liftoff_count = push_cds_idx(l_c_idx, l_lifton_aln, ref_aa_liftoff_count, DEBUG)
            elif ref_aa_liftoff_count > ref_aa_miniprot_count:
                m_c_idx, ref_aa_miniprot_count = push_cds_idx(m_c_idx, m_lifton_aln, ref_aa_miniprot_count, DEBUG)
            elif ref_aa_liftoff_count == ref_aa_miniprot_count:
                # 1. Accumulated aa in the reference of Liftoff and miniprot are the same
                # 2. CDS ends at the same position.
                if l_c_idx > 0 and m_c_idx > 0 and m_c.end == l_c.end:
                    cdss = process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, chains, DEBUG)
                    cds_list += cdss
                    m_c_idx_last = m_c_idx
                    l_c_idx_last = l_c_idx
                l_c_idx, ref_aa_liftoff_count = push_cds_idx(l_c_idx, l_lifton_aln, ref_aa_liftoff_count, DEBUG)
                m_c_idx, ref_aa_miniprot_count = push_cds_idx(m_c_idx, m_lifton_aln, ref_aa_miniprot_count, DEBUG)
        elif m_lifton_aln.db_entry.strand == "-":
            m_c = m_children[len(m_children) - m_c_idx - 1]
            l_c = l_children[len(l_children) - l_c_idx - 1]
            if ref_aa_liftoff_count < ref_aa_miniprot_count:
                l_c_idx, ref_aa_liftoff_count = push_cds_idx(l_c_idx, l_lifton_aln, ref_aa_liftoff_count, DEBUG)
            elif ref_aa_liftoff_count > ref_aa_miniprot_count:
                m_c_idx, ref_aa_miniprot_count = push_cds_idx(m_c_idx, m_lifton_aln, ref_aa_miniprot_count, DEBUG)
            elif ref_aa_liftoff_count == ref_aa_miniprot_count:
                # 1. Accumulated aa in the reference of Liftoff and miniprot are the same
                # 2. CDS ends at the same position.
                if l_c_idx > 0 and m_c_idx > 0 and m_c.end == l_c.end:
                    cdss = process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, chains, DEBUG)
                    cds_list += cdss
                    m_c_idx_last = m_c_idx
                    l_c_idx_last = l_c_idx
                l_c_idx, ref_aa_liftoff_count = push_cds_idx(l_c_idx, l_lifton_aln, ref_aa_liftoff_count, DEBUG)
                m_c_idx, ref_aa_miniprot_count = push_cds_idx(m_c_idx, m_lifton_aln, ref_aa_miniprot_count, DEBUG)
    l_c_idx += 1
    m_c_idx += 1
    cdss = process_m_l_children(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, chains, DEBUG)
    cds_list += cdss
    return cds_list, chains