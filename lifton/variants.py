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


def _coding_subalignment(align_dna, cds_start, cds_end):
    """Return (query_sub, ref_sub): the alignment columns whose QUERY (target)
    position falls within [cds_start, cds_end) of the spliced transcript.

    Used to scope frameshift detection to the CODING region (GitHub issue #46). The
    DNA alignment spans the FULL transcript (5'UTR + CDS + 3'UTR); a 1-2 bp indel in
    an UTR is not a coding frameshift, yet the unscoped check flagged it — extremely
    common on close pairs (e.g. human<->chimp), producing false 'frameshift' labels
    in score.txt and a needless ORF re-search that can perturb the emitted model.

    ``qpos`` never decreases, so the selected columns form a CONTIGUOUS run --
    only its two boundaries have to be located and the body is a slice. That
    turns a per-character Python loop that built two character lists (13.5 s,
    ~3% of Step 7 on the drosophila profile) into a boundary scan plus two
    C-level slices, with an ungapped fast path for close pairs, which is
    exactly where this check runs most often."""
    query = align_dna.query_aln
    ref = align_dna.ref_aln
    if cds_start >= cds_end:
        return "", ""
    if "-" not in query:
        # Ungapped query: column index == transcript position.
        return query[cds_start:cds_end], ref[cds_start:cds_end]
    length = len(query)
    start_col = end_col = length
    qpos = 0
    for col, qc in enumerate(query):
        if qpos == cds_start and start_col == length:
            start_col = col
        if qpos == cds_end:
            end_col = col
            break
        if qc != '-':
            qpos += 1
    return query[start_col:end_col], ref[start_col:end_col]


def find_variants(align_dna, align_protein, lifton_status, peps, is_non_coding, cds_span=None):
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
    # GH #46: scope the frameshift test to the CODING region. align_dna spans the
    # full transcript (incl. UTRs); a 1-2 bp UTR indel is not a coding frameshift, so
    # the unscoped check produced false 'frameshift' labels (+ a needless ORF re-search)
    # on close pairs. When the caller supplies the CDS span (query/target coords), test
    # only the coding sub-alignment; otherwise fall back to the full alignment (unchanged).
    if cds_span is not None:
        fs_query, fs_ref = _coding_subalignment(align_dna, cds_span[0], cds_span[1])
    else:
        fs_query, fs_ref = align_dna.query_aln, align_dna.ref_aln
    if is_frameshift(fs_query):
        mutation_type.append('frameshift')
        frameshift = True
    if is_frameshift(fs_ref) and frameshift == False:
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