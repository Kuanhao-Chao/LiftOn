"""Give a miniprot-derived model its terminal stop codon.

miniprot reports a coding alignment, and its CDS stops at the last aligned
codon: the stop codon itself is not part of the hit. Everything else LiftOn
emits follows the reference convention, where a CDS includes its stop, and both
the reference protein (``extract_sequence.get_protein_sequence``) and a lifted
one (``Lifton_TRANS.translate_coding_seq``) are translated from the full CDS. So
a miniprot-derived model ends one codon short of the protein it is compared
against.

The ORF search cannot repair this. ``Lifton_TRANS.__find_orfs`` scans the
spliced *transcript* sequence, and a miniprot model has no UTR -- its exons are
its CDS -- so the stop codon sitting immediately downstream in the genome is
outside the sequence being searched. Measured on the v1.0.12 whole-genome
output, only half the rescued models on human to zebrafish end in a stop.

This module closes that gap and nothing else: when the codon immediately 3' of
the terminal CDS is a stop and the model does not already end in one, the
terminal CDS and its exon grow by exactly three bases. The encoded amino acids
cannot change, because the added codon terminates translation.

Deliberately NOT done here: the suppression interval a rescued gene contributes
is left at its pre-extension value, so which genes are placed, and where, is
exactly what it would have been. Only coordinates and the emitted protein move.
"""
import os

from lifton import coreutils

#: Default for the terminal-stop completion, in one place so a promotion or a
#: revert is a one-line change (the pattern the other rescue switches follow).
STOP_COMPLETION_DEFAULT = True

STOP_CODONS = frozenset(("TAA", "TAG", "TGA"))

#: Bases added when the downstream codon is a stop.
STOP_CODON_LENGTH = 3


def enabled(args):
    """Is terminal-stop completion on? ``LIFTON_ORF_STOP_COMPLETION`` wins over
    the resolved flag, as every other rescue switch does."""
    env = os.environ.get("LIFTON_ORF_STOP_COMPLETION")
    if env is not None:
        return env.strip().lower() not in ("", "0", "false", "no", "off")
    resolved = getattr(args, "orf_stop_completion", None)
    return STOP_COMPLETION_DEFAULT if resolved is None else bool(resolved)


def _codon(entry, fai, start, end):
    """The three bases at ``start..end`` read in ``entry``'s orientation."""
    probe = coreutils.clone_feature(entry)
    probe.start = start
    probe.end = end
    try:
        return str(probe.sequence(fai)).upper()
    except (KeyError, ValueError, IndexError):
        return ""


def _terminal_exon(lifton_trans):
    """The exon holding the last CDS in transcript orientation, with its CDS,
    or ``(None, None)`` when the model is not a candidate: the terminal CDS has
    to be flush with its exon boundary, so growing the CDS cannot overrun a UTR
    the model already has."""
    coding = [exon for exon in lifton_trans.exons if exon.cds is not None]
    if not coding:
        return None, None
    minus = lifton_trans.entry.strand == "-"
    exon = (min(coding, key=lambda e: e.cds.entry.start) if minus
            else max(coding, key=lambda e: e.cds.entry.end))
    flush = (exon.cds.entry.start == exon.entry.start if minus
             else exon.cds.entry.end == exon.entry.end)
    return (exon, exon.cds) if flush else (None, None)


def complete_terminal_stop(lifton_trans, fai):
    """Extend this transcript's terminal CDS and exon over a downstream stop
    codon. Returns True when it did.

    Refuses unless every precondition holds: the model has coding exons, its
    total CDS length is a whole number of codons, the terminal CDS ends flush
    with its exon, the model does not already end in a stop, the next codon in
    the genome is one, and it lies inside the sequence.
    """
    if not lifton_trans.exons or lifton_trans.entry.strand not in ("+", "-"):
        return False
    exon, cds = _terminal_exon(lifton_trans)
    if cds is None:
        return False
    coding_length = sum(e.cds.entry.end - e.cds.entry.start + 1
                        for e in lifton_trans.exons if e.cds is not None)
    if coding_length < STOP_CODON_LENGTH or coding_length % 3:
        return False
    minus = lifton_trans.entry.strand == "-"
    try:
        sequence_length = len(fai[cds.entry.seqid])
    except (KeyError, TypeError):
        return False
    if minus:
        last = (cds.entry.start, cds.entry.start + STOP_CODON_LENGTH - 1)
        nxt = (cds.entry.start - STOP_CODON_LENGTH, cds.entry.start - 1)
        if nxt[0] < 1:
            return False
    else:
        last = (cds.entry.end - STOP_CODON_LENGTH + 1, cds.entry.end)
        nxt = (cds.entry.end + 1, cds.entry.end + STOP_CODON_LENGTH)
        if nxt[1] > sequence_length:
            return False
    if _codon(cds.entry, fai, *last) in STOP_CODONS:
        return False
    if _codon(cds.entry, fai, *nxt) not in STOP_CODONS:
        return False
    if minus:
        cds.entry.start -= STOP_CODON_LENGTH
        exon.entry.start = min(exon.entry.start, cds.entry.start)
    else:
        cds.entry.end += STOP_CODON_LENGTH
        exon.entry.end = max(exon.entry.end, cds.entry.end)
    return True
