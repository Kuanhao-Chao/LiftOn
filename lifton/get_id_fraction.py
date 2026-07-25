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
"""

import operator


def get_partial_id_fraction(reference, target, start, end):
    reference = reference.upper()
    target = target.upper()
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
def get_AA_id_fraction(reference, target):
    # gap-compressed BLAST identity
    reference = reference.upper()
    target = target.upper()
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
