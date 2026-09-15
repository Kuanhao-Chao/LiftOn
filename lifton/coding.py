"""Coding semantics shared by extraction, scoring and model completion."""
from lifton.exceptions import LiftOnInputError


def initial_phase(features, strand):
    """Initial incomplete-codon bases, counted once in transcript orientation.

    Later CDS phases describe codons spanning splice junctions; trimming every
    segment would delete coding bases. Older interval-only callers have no phase.
    """
    cds = [f for f in features if getattr(f, 'featuretype', 'CDS') == 'CDS']
    if not cds:
        return 0
    first = max(cds, key=lambda f: f.end) if strand == '-' else min(cds, key=lambda f: f.start)
    phase = str(getattr(first, 'frame', '0'))
    if phase == '.':
        return 0
    if phase not in ('0', '1', '2'):
        raise LiftOnInputError(f'Invalid initial CDS phase {phase!r} for {getattr(first, "id", "CDS")}')
    return int(phase)


def phase_adjusted_lengths(lengths, phase):
    """Lengths contributing to translated sequence, preserving segment count."""
    adjusted = []
    for length in lengths:
        skipped = min(length, phase)
        adjusted.append(length - skipped)
        phase -= skipped
    return adjusted
