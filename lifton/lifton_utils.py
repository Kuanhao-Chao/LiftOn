def segments_overlap(segment1, segment2):
    # Check if the segments have valid endpoints
    # print("Checking two segments overlapping.!")

    # print(segment1, segment2)

    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    
    # Sort the segments by their left endpoints
    segment1, segment2 = sorted([segment1, segment2], key=lambda x: x[0])


    # Check if the right endpoint of the first segment is greater than or equal to the left endpoint of the second segment
    # print(segment1[1] >= segment2[0])

    return segment1[1] >= segment2[0]


def custom_bisect_insert(sorted_list, element_to_insert):
    low = 0
    high = len(sorted_list)

    while low < high:
        mid = (low + high) // 2
        if sorted_list[mid].entry.end < element_to_insert.entry.end:
            low = mid + 1
        else:
            high = mid

    sorted_list.insert(low, element_to_insert)