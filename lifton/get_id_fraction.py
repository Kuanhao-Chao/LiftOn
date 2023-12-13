def get_id_fraction(reference, target, start, end):
    matches = 0
    # print("reference[start:end]: ", reference[start:end])
    # print("target[start:end]   : ", target[start:end])
    reference = reference.upper()
    target = target.upper()
    for i, letter in enumerate(reference[start:end]):
        if letter == target[i+start]:
            matches += 1
        if target[i] == "*":
            break
    if (end-start) == 0:
        return matches, 1
    return matches, end-start


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