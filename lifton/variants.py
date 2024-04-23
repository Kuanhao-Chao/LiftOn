def has_stop_codon(ref_align, target_align):
    """
        This function checks if a string contains a stop codon.

        Parameters:
        - ref_align: reference sequence alignment
        - target_align: target sequence alignment

        Returns:
        True if a stop codon is found, False otherwise.
    """
    for i, letter in enumerate(target_align):
        if letter == "*" and ref_align[i] != "*":
            return True
    return False


def is_frameshift(s):
    """
        This function checks if a string contains a non-divisible substring of '-' characters.

        Parameters:
        - s: a string

        Returns:
        True if s contains a non-divisible substring of '-' characters, False otherwise.
    """
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


def find_variants(align_dna, align_protein, lifton_status, peps, is_non_coding):
    """
        This function finds the variants between two sequences.

        Parameters:
        - align_dna: DNA pairwise alignment
        - align_protein: Protein pairwise alignment
        - lifton_status: Lifton_Status object
        - peps: Protein sequence split by stop codon

        Returns:
        None
    """
    # Mutation types:
    #   (1) identical
    #   (2) synonymous
    #   (3) frameshift
    #   (4) start_lost
    #   (5) inframe_insertion
    #   (6) inframe_deletion
    #   (7) nonsynonymous
    #   (8) stop_missing
    #   (9) stop_codon_gain
    mutation_type = []
    if is_non_coding:
        mutation_type.append('non_coding')
        lifton_status.status = mutation_type
        return
    if align_dna == None:
        mutation_type.append('full_transcript_loss')
        lifton_status.status = mutation_type
        return 
    if align_protein == None:
        mutation_type.append('no_protein')
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
    lifton_status.status = mutation_type