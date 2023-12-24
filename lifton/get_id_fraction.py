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
    if (end-start) == 0:
        return matches, 1
    
    # Modify the region length by considering gaps in the reference (as long as it's a open reading frame)
    total_length = (end-start) - gaps_in_ref
    return matches, total_length


# I collapsed gaps at the start and the end of the alignment.
def get_AA_id_fraction(reference, target):
    matches = 0
    reference = reference.upper()
    target = target.upper()
    gaps_in_ref = 0

    # # Count leading '-' characters
    # leading_gaps = 0
    # for i in range(min(len(reference), len(target))):
    #     if reference[i] == '-':
    #         leading_gaps += 1
    #     else:
    #         break

    # # Count trailing '-' characters
    # trailing_gaps = 0
    # for i in range(1, min(len(reference), len(target)) + 1):
    #     if reference[-i] == '-':
    #         trailing_gaps += 1
    #     else:
    #         break

    # # Only count 1 penalty for leading and trailing gaps
    # if leading_gaps > 0:
    #     leading_gaps -= 1
    # if trailing_gaps > 0:
    #     trailing_gaps -= 1
    
    # # Adjust total length by considering leading and trailing gaps
    # total_length = max(len(reference), len(target)) - leading_gaps - trailing_gaps
    
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
    for i, letter in enumerate(reference):
        if letter == target[i]:
            matches += 1
    if max(len(reference), len(target)) == 0:  
        return matches, 1
    return matches, max(len(reference), len(target))