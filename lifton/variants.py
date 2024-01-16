def has_stop_codon(ref_align, target_align):
    for i, letter in enumerate(target_align):
        if letter == "*" and ref_align[i] != "*":
            return True
    return False


def is_frameshift(s):
    # Initialize a variable to keep track of consecutive '-' characters.
    consecutive_count = 0
    for char in s:
        if char == '-':
            consecutive_count += 1
        else:
            # If we encounter a non-'-' character, check if the consecutive count is not divisible by three.
            if consecutive_count % 3 != 0:
                return True
            consecutive_count = 0
    # After the loop, check the last substring if it's not divisible by three.
    if consecutive_count % 3 != 0:
        return True    
    # If no non-divisible substrings are found, return False.
    return False


def find_variants(align_dna, align_protein, lifton_status, peps):
    # Mutation types:
    #   identical
    #   synonymous
    #   frameshift
    #   start_lost
    #   inframe_insertion
    #   inframe_deletion
    #   nonsynonymous
    #   stop_missing
    #   stop_codon_gain
    mutation_type = []
    if align_dna == None:
        mutation_type.append('full_transcript_loss')
        lifton_status.status = mutation_type
        return 
    if align_protein == None:
        mutation_type.append('full_protein_loss')
        lifton_status.status = mutation_type
        return 
    # 1. return cases
    if align_dna.identity == 1.0:
        mutation_type.append('identical')
        lifton_status.status = mutation_type
        return 
    # 2. return cases
    frameshift = False
    if align_protein.identity == 1.0:
        mutation_type.append('synonymous')
        lifton_status.status = mutation_type
        return
    # 3. return cases
    if is_frameshift(align_dna.query_aln):
        mutation_type.append('frameshift')
        frameshift = True
    if is_frameshift(align_dna.ref_aln) and frameshift == False:
        mutation_type.append('frameshift')
        frameshift = True
    if align_dna.query_aln[0:3] != align_dna.ref_aln[0:3] and \
        align_dna.query_aln[0:3] != 'ATG' and \
            align_protein.query_aln[0] != align_protein.ref_aln[0] and \
                align_protein.query_aln[0] != "M" :
        mutation_type.append("start_lost")

    if "-" in align_dna.ref_aln and not frameshift:
        mutation_type.append("inframe_insertion")

    if "-" in align_dna.query_aln and not frameshift:
        mutation_type.append("inframe_deletion")

    if len(peps) == 2 and str(peps[1]) == "":
        # This is a valid protein ends with stop codon *
        # But alignment is not 100% identical
        if len(mutation_type) == 0:
            mutation_type.append("nonsynonymous")
        
    elif len(peps) == 1:
        # This is a protein without stop codon
        mutation_type.append("stop_missing")
    
    else:
        mutation_type.append("stop_codon_gain")
        stop_codon_count = 0
        for idx, ele in enumerate(align_protein.query_seq):
            if ele == "*" and idx != len(align_protein.query_seq)-1:
                stop_codon_count += 1
        # self.entry.attributes["StopCodon"] = str(stop_codon_count)
        # lifton_status.status = "early_stop_codon"

        # print("extracted_parasail_res: ", extracted_parasail_res)
        # print("extracted_seq: ", extracted_seq)
        # print("reference_seq: ", reference_seq)
        # print("stop_codon_count: ", stop_codon_count)

        # self.__find_orfs(trans_seq, exon_lens, ref_protein_seq, lifton_aa_aln, lifton_status)
        # return lifton_aa_aln, False

    lifton_status.status = mutation_type