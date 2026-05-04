def get_partial_id_fraction(reference, target, start, end):
    matches = 0
    # print("reference[start:end]: ", reference[start:end])
    # print("target[start:end]   : ", target[start:end])
    reference = reference.upper()
    target = target.upper()
    gaps_in_ref = 0
    for i, letter in enumerate(reference[start:end]):
        if letter == '-':
            gaps_in_ref += 1
        if letter == target[i+start]:
            matches += 1
        if target[i+start] == "*":
            break
    # Modify the region length by considering gaps in the reference (as long as it's a open reading frame)
    total_length = (end-start) - gaps_in_ref
    if total_length == 0:
        return matches, 1
    return matches, total_length


# Gap-collapsed protein sequence identity
def get_AA_id_fraction(reference, target):
    matches = 0
    # gap-compressed BLAST identity
    reference = reference.upper()
    target = target.upper()
    gaps_in_ref = 0
    for i, letter in enumerate(reference):
        if letter == '-':
            gaps_in_ref += 1
        if letter == target[i]:
            matches += 1
        if target[i] == "*":
            break        
    if max(len(reference), len(target)) == 0:  
        return matches, 1
    # Modify the region length by considering gaps in the reference (as long as it's a open reading frame)
    total_length = max(len(reference), len(target)) - gaps_in_ref    
    return matches, total_length


def get_DNA_id_fraction(reference, target):
    matches = 0
    # BLAST identity
    reference = reference.upper()
    target = target.upper()
    # V2.3 fix: guard length mismatch BEFORE the indexing loop.
    # parasail traceback strings are always equal-length, but a future
    # refactor (or a Pandas round-trip) could produce different-length
    # inputs and crash this loop with IndexError mid-iteration.
    if len(reference) != len(target):
        raise ValueError(
            f"get_DNA_id_fraction: reference length ({len(reference)}) "
            f"does not match target length ({len(target)}). The two "
            "sequences must be aligned to equal length before identity "
            "is computed."
        )
    for i, letter in enumerate(reference):
        if letter == target[i]:
            matches += 1
    if max(len(reference), len(target)) == 0:
        return matches, 1
    return matches, max(len(reference), len(target))