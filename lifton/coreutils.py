"""Dependency-free core helpers (Iteration 16).

These three pure functions used to live in ``lifton.lifton_utils``. They were
moved here to break the ``lifton_utils`` <-> ``lifton_class`` package-level
import cycle: ``lifton_class`` needs only these helpers from ``lifton_utils``,
so importing this leaf module instead removes that edge entirely. This module
must import nothing from elsewhere in ``lifton`` so it stays a true leaf —
``lifton_utils`` re-exports these names for backward compatibility.
"""


def custom_bisect_insert(sorted_list, element_to_insert):
    """
        This function bisects the sorted list and inserts the element.

        Parameters:
        - sorted_list: sorted list
        - element_to_insert: element to insert

        Returns:
        None
    """
    low = 0
    high = len(sorted_list)

    while low < high:
        mid = (low + high) // 2
        if sorted_list[mid].entry.end < element_to_insert.entry.end:
            low = mid + 1
        else:
            high = mid
    sorted_list.insert(low, element_to_insert)


def get_ID_base(id, ref_features_dict=None):
    """
        This function gets the ID base by removing the last substring after "_" if it's a copy number.
        Only removes the suffix if it's confirmed to be a copy number added by LiftOn.

        The key insight: We can only safely remove a suffix if the ID without the suffix exists
        in the reference features dictionary. This confirms the suffix was added by LiftOn,
        not part of the original ID.

        Parameters:
        - id: ID
        - ref_features_dict: Optional reference features dictionary to verify if the base ID exists.
                            If provided, only removes suffix if base ID exists in the dictionary.

        Returns:
        ID base (original ID if suffix removal is not safe or cannot be verified)
    """
    splits = id.split("_")
    # Only try to remove suffix if there are at least 2 parts after splitting.
    # This also guards against the Phase 5 bug #5 case where a single-component
    # numeric id (e.g. "0") would otherwise be reduced to the empty string.
    if len(splits) < 2:
        return id

    # Check if the last part is a pure integer (potential copy number)
    try:
        last_part = splits[-1]
        copy_num = int(last_part)
        id_base = "_".join(splits[:-1])

        # If ref_features_dict is provided, verify that removing the suffix gives a valid ID.
        # This is the safest approach: only remove if we can confirm the base exists.
        if ref_features_dict is not None:
            # Only remove the suffix if the base ID exists in the reference dictionary
            # This confirms the suffix was added by LiftOn, not part of the original ID
            if id_base in ref_features_dict.keys():
                return id_base
            else:
                # The suffix is likely part of the original ID (e.g., FMUND_1), don't remove it
                return id
        else:
            # Without reference dict, we cannot safely determine if the suffix is a copy number.
            # Be conservative: don't remove the suffix to avoid breaking IDs that naturally end
            # with numbers (e.g. FMUND_1 in GCA_013396205.1-transcript_rna-gnl-WGS:JAAOAN-mrna.FMUND_1).
            return id
    except (ValueError, IndexError):
        # Last part is not an integer, so it's not a copy number suffix
        return id


def segments_overlap_length(segment1, segment2):
    """
        This function gets the length of the overlapping segments.

        Parameters:
        - segment1: segment 1 in tuple (start, end)
        - segment2: segment 2 in tuple (start, end)

        Returns:
        The length of the overlapping segments.
    """
    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    # Bug fix #4 (Phase 5): make the sort tie-break deterministic so the
    # function is symmetric. Previously sorting by start alone left the
    # tie-break to Python's stable-sort and the input order leaked into
    # the result. Computing the overlap length directly from the sorted
    # endpoints removes the ambiguity entirely.
    s1, e1 = segment1
    s2, e2 = segment2
    ovp_len = min(e1, e2) - max(s1, s2) + 1
    # V2.9 fix: disjoint segments produced negative overlap, which
    # callers treating the value as "bp shared" would silently divide
    # by. Clamp to 0 so the returned tuple always has ovp_len >= 0
    # AND ovp_len > 0 ⟺ ovp.
    if ovp_len < 0:
        ovp_len = 0
    ovp = ovp_len > 0
    return ovp_len, ovp
