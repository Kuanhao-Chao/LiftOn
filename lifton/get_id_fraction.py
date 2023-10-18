def get_id_fraction(reference, target, start, end):
    matches = 0

    # print("reference[start:end]: ", reference[start:end])
    # print("target[start:end]   : ", target[start:end])
    for i, letter in enumerate(reference[start:end]):
        if letter == target[i+start]:
            matches += 1
        if target[i] == "*":
            break
    return matches, end-start