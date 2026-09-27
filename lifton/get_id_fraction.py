"""Sequence-identity counters over parasail traceback strings.

These run after EVERY alignment -- and the best-of-outcome merge plus candidate-3 plus
the ORF rescue mean roughly 6-12 identity scans per transcript -- so they sit directly in
the Step-7 hot loop. They used to be per-character Python `for` loops, which measured
8-25% of the cost of the alignment+scoring pair and, being pure Python, also held the GIL
and capped `--threads N` scaling.

The implementations below are the same computations expressed with C-level primitives
(`str.count`, `str.find`, `map(operator.eq, ...)`). Semantics are preserved EXACTLY,
including the subtle parts:

* the `*` (stop codon) scan stops AFTER processing the position it is found at, and only
  within the range the original loop visited;
* gaps are counted in the reference over that same truncated range;
* `total_length` uses the full requested span, not the truncated one;
* the `(matches, 1)` degenerate returns keep identity bounded in [0, 1].

`tests/test_get_id_fraction_equivalence.py` fuzzes these against the literal loop forms.

**Declared read-through.** A codon the annotation declares to translate through a stop
(`transl_except`: selenocysteine, pyrrolysine, stop readthrough) is `*` in both the
reference and the lifted protein, and the first-`*` cut-off above scored an identical
selenoprotein as residue/length (SEPHS2: 60/449). Callers that know the declared codons
pass their alignment columns as `readthrough_cols`: a target `*` in one of those columns
neither ends the scan nor costs a match. With no columns the original bodies run
unchanged.
"""

import operator


def readthrough_columns(ref_aln, query_aln, residues):
    """Alignment columns where the target reads through a declared stop.

    ``residues`` are 0-based reference-protein indices (``transl_except``
    read-through codons). A column qualifies when its reference residue is one
    of them and the target has ``*`` there -- except the target's last residue:
    a model that ENDS at a selenocysteine codon stops there, it does not read
    through it.
    """
    if not residues:
        return frozenset()
    last_query = max((i for i, char in enumerate(query_aln) if char != "-"),
                     default=-1)
    columns = set()
    residue = -1
    for column, ref_char in enumerate(ref_aln):
        if ref_char == "-":
            continue
        residue += 1
        if (residue in residues and query_aln[column] == "*"
                and column != last_query):
            columns.add(column)
    return frozenset(columns)


def readthrough_query_residues(query_aln, columns):
    """0-based target-protein indices of the read-through columns."""
    if not columns:
        return frozenset()
    residues = set()
    residue = -1
    for column, char in enumerate(query_aln):
        if char == "-":
            continue
        residue += 1
        if column in columns:
            residues.add(residue)
    return frozenset(residues)


def mask_readthrough(protein, query_aln, columns):
    """``protein`` with its read-through stops written as ``U``, so a split on
    ``*`` sees only the stops that end translation."""
    residues = readthrough_query_residues(query_aln, columns)
    if not residues:
        return protein
    return "".join("U" if i in residues and char == "*" else char
                   for i, char in enumerate(protein))


def _first_stop(target, start, end, readthrough_cols):
    stop = target.find("*", start, end)
    while stop >= 0 and stop in readthrough_cols:
        stop = target.find("*", stop + 1, end)
    return stop


def _matches(reference, target, start, limit, readthrough_cols):
    return sum(1 for i in range(start, start + limit)
               if reference[i] == target[i] or i in readthrough_cols)


def get_partial_id_fraction(reference, target, start, end, readthrough_cols=None):
    reference = reference.upper()
    target = target.upper()
    if readthrough_cols:
        window = len(reference[start:end])
        stop = _first_stop(target, start, start + window, readthrough_cols)
        limit = window if stop < 0 else (stop - start) + 1
        matches = _matches(reference, target, start, limit, readthrough_cols)
        gaps_in_ref = reference[start:end].count("-", 0, limit)
        total_length = (end - start) - gaps_in_ref
        if total_length == 0:
            return matches, 1
        return matches, total_length
    ref_segment = reference[start:end]
    # The original loop indexed target[i + start] and broke AFTER the position whose
    # target character is '*', so only stops inside the visited window count.
    stop = target.find("*", start, start + len(ref_segment))
    limit = len(ref_segment) if stop < 0 else (stop - start) + 1
    matches = sum(map(operator.eq,
                      ref_segment[:limit], target[start:start + limit]))
    gaps_in_ref = ref_segment.count("-", 0, limit)
    # Modify the region length by considering gaps in the reference (as long as it's a
    # open reading frame). NOTE: the span is the requested (end - start), not `limit`.
    total_length = (end - start) - gaps_in_ref
    if total_length == 0:
        return matches, 1
    return matches, total_length


# Gap-collapsed protein sequence identity
def get_AA_id_fraction(reference, target, readthrough_cols=None):
    # gap-compressed BLAST identity
    reference = reference.upper()
    target = target.upper()
    if readthrough_cols:
        stop = _first_stop(target, 0, len(reference), readthrough_cols)
        limit = len(reference) if stop < 0 else stop + 1
        matches = _matches(reference, target, 0, limit, readthrough_cols)
        gaps_in_ref = reference.count("-", 0, limit)
        span = max(len(reference), len(target))
        if span == 0:
            return matches, 1
        total_length = span - gaps_in_ref
        if total_length <= 0:
            return matches, 1
        return matches, total_length
    stop = target.find("*", 0, len(reference))
    limit = len(reference) if stop < 0 else stop + 1
    matches = sum(map(operator.eq, reference[:limit], target[:limit]))
    gaps_in_ref = reference.count("-", 0, limit)
    span = max(len(reference), len(target))
    if span == 0:
        return matches, 1
    # Modify the region length by considering gaps in the reference (as long as it's a
    # open reading frame).
    total_length = span - gaps_in_ref
    # V2.4 fix: when the reference is all gaps (or a stop-codon truncation
    # leaves total_length == 0), avoid producing a denominator that would
    # ZeroDivision in the caller. Identity is undefined in this case;
    # returning (matches, 1) keeps the value bounded in [0, 1].
    if total_length <= 0:
        return matches, 1
    return matches, total_length


def get_DNA_id_fraction(reference, target):
    # BLAST identity
    reference = reference.upper()
    target = target.upper()
    # V2.3 fix: guard length mismatch BEFORE the counting pass.
    # parasail traceback strings are always equal-length, but a future
    # refactor (or a Pandas round-trip) could produce different-length
    # inputs and silently truncate the comparison.
    if len(reference) != len(target):
        raise ValueError(
            f"get_DNA_id_fraction: reference length ({len(reference)}) "
            f"does not match target length ({len(target)}). The two "
            "sequences must be aligned to equal length before identity "
            "is computed."
        )
    matches = sum(map(operator.eq, reference, target))
    span = max(len(reference), len(target))
    if span == 0:
        return matches, 1
    return matches, span
