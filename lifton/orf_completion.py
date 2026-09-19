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

Three rules keep it from being able to make anything worse:

* It runs **after** the ORF search, so the ORF search sees exactly the sequence
  it sees today and cannot take a different path. Applying it before instead
  cost one transcript's identity on the drosophila-to-anopheles ladder cell,
  because completing the model suppressed the ``stop_missing`` mutation that had
  been triggering the search.
* It fires only when the reference protein itself ends in a stop. An annotation
  whose CDS excludes the stop produces a reference protein without one, and
  adding a stop the reference does not have would be a mismatch, not a match.
* The suppression interval a rescued gene contributes is left at its
  pre-extension value, so which genes are placed, and where, is exactly what it
  would have been. Only coordinates and the emitted protein move.

Only the protein identity is re-derived. ``dna_identity`` still describes the
model before the three bases were added, because re-deriving it means realigning
the whole transcript -- the most expensive part of scoring -- to move a value by
three bases in several hundred. The evaluator translates the emitted CDS itself,
so no measurement depends on the attribute.
"""
import os

from lifton import align, coreutils

#: Default for the terminal-stop completion, in one place so a promotion or a
#: revert is a one-line change (the pattern the other rescue switches follow).
STOP_COMPLETION_DEFAULT = True

#: Standard-code stop codons. A model whose annotation declares another code
#: is completed against that code's stops instead -- see `_stop_codons`.
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


def _stop_codons(lifton_trans):
    """The stop codons of the code this model declares."""
    resolve = getattr(lifton_trans, "transl_table", None)
    if resolve is None:
        return STOP_CODONS
    from lifton import coding
    return coding.stop_codons(resolve())


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


def reference_protein_has_stop(ref_proteins, ref_trans_id):
    """Does the reference protein for this transcript end in a stop? Only then
    does completing a model move it toward the reference rather than away."""
    try:
        protein = ref_proteins[ref_trans_id]
    except (KeyError, TypeError):
        return False
    protein = str(protein)
    return bool(protein) and protein.rstrip().endswith("*")


def complete_terminal_stop(lifton_trans, fai):
    """Extend this transcript's terminal CDS and exon over a downstream stop
    codon. Returns True when it did."""
    return apply_terminal_stop(lifton_trans, fai) is not None


def apply_terminal_stop(lifton_trans, fai):
    """Extend the terminal CDS and its exon over a downstream stop codon.

    Returns the ``(exon, cds, before)`` it changed so the caller can undo it,
    or None when it changed nothing. Refuses unless every precondition holds:
    the model has coding exons, its total CDS length is a whole number of
    codons, the terminal CDS ends flush with its exon, the model does not
    already end in a stop, the next codon in the genome is one, and it lies
    inside the sequence.
    """
    if not lifton_trans.exons or lifton_trans.entry.strand not in ("+", "-"):
        return None
    exon, cds = _terminal_exon(lifton_trans)
    if cds is None:
        return None
    coding_length = sum(e.cds.entry.end - e.cds.entry.start + 1
                        for e in lifton_trans.exons if e.cds is not None)
    from lifton.coding import initial_phase
    coding_length -= initial_phase([e.cds.entry for e in lifton_trans.exons if e.cds is not None],
                                    lifton_trans.entry.strand)
    if coding_length < STOP_CODON_LENGTH or coding_length % 3:
        return None
    minus = lifton_trans.entry.strand == "-"
    try:
        sequence_length = len(fai[cds.entry.seqid])
    except (KeyError, TypeError):
        return None
    if minus:
        last = (cds.entry.start, cds.entry.start + STOP_CODON_LENGTH - 1)
        nxt = (cds.entry.start - STOP_CODON_LENGTH, cds.entry.start - 1)
        if nxt[0] < 1:
            return None
    else:
        last = (cds.entry.end - STOP_CODON_LENGTH + 1, cds.entry.end)
        nxt = (cds.entry.end + 1, cds.entry.end + STOP_CODON_LENGTH)
        if nxt[1] > sequence_length:
            return None
    stops = _stop_codons(lifton_trans)
    if _codon(cds.entry, fai, *last) in stops:
        return None
    if _codon(cds.entry, fai, *nxt) not in stops:
        return None
    before = (exon.entry.start, exon.entry.end, cds.entry.start, cds.entry.end)
    if minus:
        cds.entry.start -= STOP_CODON_LENGTH
        exon.entry.start = min(exon.entry.start, cds.entry.start)
    else:
        cds.entry.end += STOP_CODON_LENGTH
        exon.entry.end = max(exon.entry.end, cds.entry.end)
    return exon, cds, before


def undo_terminal_stop(applied):
    """Put back the coordinates ``apply_terminal_stop`` changed."""
    exon, cds, (exon_start, exon_end, cds_start, cds_end) = applied
    exon.entry.start, exon.entry.end = exon_start, exon_end
    cds.entry.start, cds.entry.end = cds_start, cds_end


def complete_and_rescore(lifton_trans, m_entry, fai, ref_proteins, ref_trans_id,
                         lifton_status, args=None, enabled_override=None):
    """Complete this miniprot-derived model and refresh its identity.

    Call after the ORF search and before the status attributes are written, so
    the ORF search is unaffected and the recorded identity describes the model
    that is actually emitted. Returns True when the model was completed.
    """
    on = enabled_override if enabled_override is not None else enabled(args)
    if not on or not reference_protein_has_stop(ref_proteins, ref_trans_id):
        return False
    applied = apply_terminal_stop(lifton_trans, fai)
    if applied is None:
        return False
    # Appending a residue can shift a global alignment, so the stop the model
    # gains is usually a new match but occasionally costs one elsewhere: on
    # drosophila -> anopheles, 1 of 122 completed transcripts lost 0.0013
    # identity that way. Score it and keep the extension only when it does not
    # make the model worse, which makes completion non-regressing by
    # construction rather than by argument.
    before = lifton_status.lifton_aa
    alignment = align.lifton_parasail_align(lifton_trans, m_entry, fai,
                                            ref_proteins, ref_trans_id)
    identity = getattr(alignment, "identity", None)
    if identity is None or (before is not None and identity < before):
        undo_terminal_stop(applied)
        return False
    lifton_status.lifton_aa = identity
    lifton_trans.entry.attributes["orf_stop_completed"] = ["true"]
    return True
