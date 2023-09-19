from lifton import write_lifton_entry, lifton_class

def create_lifton_entries(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, fw, miniprot_is_better):
    cds_list = []
    if miniprot_is_better:
        # print(">> miniprot is better!")
        for c_idx in range(m_c_idx_last, m_c_idx):
            lifton_cds = m_lifton_aln.cds_children[c_idx]
            lifton_cds.attributes = l_lifton_aln.cds_children[0].attributes

            cds = lifton_class.Lifton_CDS(lifton_cds)
            cds_list.append(cds)

            # print(lifton_cds)
            # print(lifton_cds.attributes)
            # write_lifton_entry.write_lifton_entry(fw, lifton_cds)
    else:
        # print(">> liftoff is better!")
        for c_idx in range(l_c_idx_last, l_c_idx):
            lifton_cds = l_lifton_aln.cds_children[c_idx]
            
            cds = lifton_class.Lifton_CDS(lifton_cds)
            cds_list.append(cds)

            # print(lifton_cds)
            # write_lifton_entry.write_lifton_entry(fw, lifton_cds)
            # for f in c:
            #     print("\t", f)

    return cds_list