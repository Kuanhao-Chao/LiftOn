
def adjust_cdss_protein_boundary(cdss_protein_aln_boundary, cigar_accum_len, length):
    
    cds_boundary_shift = 0
    # cdss_protein_aln_boundary = cdss_protein_boundary.copy()
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
        # print("### cdss_start, cdss_end: ", cdss_start, cdss_end)
        cdss_protein_aln_boundary[i] = (cdss_start, cdss_end)
    return cdss_protein_aln_boundary