from lifton import write_lifton_entry

def create_lifton_entries(m_c_idx, m_c_idx_last, m_lifton_aln, l_c_idx, l_c_idx_last, l_lifton_aln, fai, fw, miniprot_is_better):
    if miniprot_is_better:
        for c_idx in range(m_c_idx_last, m_c_idx):
            lifton_cds = m_lifton_aln.cds_children[c_idx]
            # lifton_cds = c.copy()
            lifton_cds.source = "Lifton"
            lifton_cds.attributes = l_lifton_aln.cds_children[0].attributes
            print(lifton_cds)
            print(lifton_cds.attributes)
    else:
        for c_idx in range(l_c_idx_last, l_c_idx):
            lifton_cds = l_lifton_aln.cds_children[c_idx]
            # lifton_cds = c.copy()
            lifton_cds.source = "Lifton"
            print(lifton_cds)
            print(lifton_cds.attributes)

            # for f in c:
            #     print("\t", f)

    write_lifton_entry.write_lifton_entry(fw, lifton_cds)



