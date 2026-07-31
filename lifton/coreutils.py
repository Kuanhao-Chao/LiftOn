"""Dependency-free core helpers (Iteration 16).

The first three pure functions used to live in ``lifton.lifton_utils``. They
were moved here to break the ``lifton_utils`` <-> ``lifton_class``
package-level import cycle: ``lifton_class`` needs only these helpers from
``lifton_utils``, so importing this leaf module instead removes that edge
entirely. This module must import nothing from elsewhere in ``lifton`` so it
stays a true leaf — ``lifton_utils`` re-exports these names for backward
compatibility.

``clone_feature`` / ``clone_attributes`` were added later for the same reason:
they are shared pure helpers with no LiftOn dependencies, so they belong in the
leaf rather than in whichever module happened to need them first.
"""

from itertools import groupby


def parse_cigar_items(cigar_string):
    """Yield ``(length, operation)`` pairs from a CIGAR string.

    This replaces the single call LiftOn made into the third-party ``cigar``
    package (``list(Cigar(s).items())`` in :mod:`lifton.align`), which was
    dropped as a dependency: ``cigar`` 0.1.3 (2015, unmaintained) ships only an
    sdist whose ``setup.py`` calls ``ez_setup.use_setuptools()``, downloading a
    setuptools egg from a URL that no longer serves one. Any ``pip install
    lifton`` without a cached wheel therefore died with
    ``tarfile.ReadError: not a gzip file``.

    The upstream generator also ran ``raise StopIteration`` on a ``*`` CIGAR,
    which PEP 479 turned into a ``RuntimeError`` on Python >= 3.7; this returns
    instead. Behaviour is otherwise identical, including yielding ``(0, None)``
    for the unavailable-alignment CIGAR ``*``.
    """
    if cigar_string == "*":
        yield (0, None)
        return
    groups = groupby(cigar_string, lambda character: character.isdigit())
    for _, digits in groups:
        yield int("".join(digits)), "".join(next(groups)[1])


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


def clone_attributes(attributes):
    """Copy a gffutils attribute mapping deeply enough to be independent.

    An attribute value is a list of strings, so ``list(v)`` copies everything
    that can be mutated; strings are immutable. ``copy.deepcopy`` reaches the
    same result by recursing into every string object, which is where a large
    part of the Step-7 copy cost went (43% of a Feature deepcopy, and more now
    that merged CDS rows carry the reference's full attribute set).

    The result keeps the source's class, so a ``gffutils.attributes.Attributes``
    stays one and a plain dict stays a dict.
    """
    if attributes is None:
        return None
    clone = attributes.__class__()
    for key, value in attributes.items():
        clone[key] = list(value) if isinstance(value, (list, tuple)) else value
    return clone


def _iter_slots(cls):
    """Every ``__slots__`` name declared anywhere in ``cls``'s MRO."""
    seen = []
    for klass in cls.__mro__:
        slots = klass.__dict__.get("__slots__", ())
        if isinstance(slots, str):
            slots = (slots,)
        for name in slots:
            if name not in ("__dict__", "__weakref__") and name not in seen:
                seen.append(name)
    return seen


def clone_feature(feature):
    """Return an independent copy of a ``gffutils.Feature``.

    Equivalent to ``copy.deepcopy`` for every field LiftOn reads or writes, and
    measured ~6.6x faster on a CDS carrying a full RefSeq attribute set
    (33.1 us -> 5.0 us). ``copy.deepcopy`` was the second-hottest entry in both
    Step-7 profiles, at 36 M calls.

    Everything on a Feature except ``attributes``, ``extra`` and ``dialect`` is
    an immutable scalar, so a ``__dict__`` copy plus fresh copies of those two
    mutable containers is a complete clone.

    ``dialect`` is shared ON PURPOSE, and that is a quarter of the win: it is a
    single parser-metadata dict shared by every feature of a database already,
    and nothing in LiftOn writes to it (the only readers are
    ``dialect['fmt']`` in the two vendored Liftoff writers). Should that ever
    change, this is the line to revisit.

    BOTH feature classes must work. ``gffutils.Feature`` keeps its state in
    ``__dict__``; the vendored ``lifton.gffbase.feature.Feature`` declares
    ``__slots__`` and therefore has no ``__dict__`` at all, which is the class
    the ``--inmemory-liftoff`` path uses.
    """
    cls = feature.__class__
    clone = cls.__new__(cls)
    state = getattr(feature, "__dict__", None)
    if state is not None:
        clone.__dict__.update(state)
    else:
        for name in _iter_slots(cls):
            try:
                setattr(clone, name, getattr(feature, name))
            except AttributeError:
                # An unset slot stays unset on the clone, as deepcopy leaves it.
                pass
    clone.attributes = clone_attributes(feature.attributes)
    extra = getattr(feature, "extra", None)
    if isinstance(extra, (list, tuple)):
        clone.extra = list(extra)
    return clone
