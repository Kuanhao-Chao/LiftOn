"""Tier-1 audit fix — the chaining whole-side fallbacks must honour the CDS order
contract that ``Lifton_TRANS.update_cds_list`` relies on.

`create_lifton_entries` (the normal chaining path) walks a "-"-strand transcript with
``c_idx_fix = n - c_idx - 1``, so it emits CDS in **descending** genomic order, and
``update_cds_list`` applies an unconditional ``cds_list.reverse()`` for "-" to get
ascending order.

The two whole-side fallbacks in ``chaining_algorithm`` (miniprot produced nothing /
Liftoff produced nothing) used to hand back ``cds_children`` directly, which is
**ascending** (queried ``order_by='start'``). That same ``reverse()`` then INVERTED
them, so ``update_cds_list`` would rebuild the exon/CDS structure from a back-to-front
list. The bug was latent only because an empty ``cds_children`` implied a ``None``
alignment that crashed earlier — fixing that crash makes this path live, so both fixes
land together.
"""
from types import SimpleNamespace
from unittest.mock import MagicMock

from lifton import lifton_class, protein_maximization


def _cds(start, end, strand):
    f = MagicMock()
    f.seqid, f.start, f.end, f.strand = "chr1", start, end, strand
    f.featuretype, f.frame = "CDS", "0"
    f.attributes = {"Parent": ["tx1"], "ID": [f"cds_{start}"]}
    return f


def _aln(cds_children, strand):
    return lifton_class.Lifton_Alignment(
        extracted_identity=0.0, cds_children=cds_children,
        alignment_query="", alignment_comp="", alignment_ref="",
        cdss_protein_boundary={}, cdss_protein_aln_boundary=[],
        extracted_seq="", reference_seq="",
        db_entry=SimpleNamespace(strand=strand),
    )


def _spans(result):
    return [(c.entry.start, c.entry.end) for c in result]


ASCENDING = [(100, 199), (300, 399), (500, 599)]


def _run(empty_side, strand):
    """empty_side='miniprot' -> Liftoff fallback; 'liftoff' -> miniprot fallback."""
    children = [_cds(s, e, strand) for s, e in ASCENDING]
    if empty_side == "miniprot":
        l_aln, m_aln = _aln(children, strand), _aln([], strand)
    else:
        l_aln, m_aln = _aln([], strand), _aln(children, strand)
    result, _ = protein_maximization.chaining_algorithm(
        l_aln, m_aln, fai=None, DEBUG=False)
    return result


def test_plus_strand_fallbacks_stay_ascending():
    """"+" needs no reversal in update_cds_list, so the fallback must stay ascending."""
    for empty_side in ("miniprot", "liftoff"):
        assert _spans(_run(empty_side, "+")) == ASCENDING


def test_minus_strand_fallbacks_emit_descending_like_the_normal_path():
    """"-" is reversed by update_cds_list, so the fallback must emit DESCENDING —
    matching create_lifton_entries — otherwise the reverse() inverts it."""
    for empty_side in ("miniprot", "liftoff"):
        assert _spans(_run(empty_side, "-")) == list(reversed(ASCENDING))


def test_minus_strand_fallback_survives_update_cds_list_reverse():
    """End-to-end contract: after update_cds_list's reverse(), the CDS list is
    ascending — i.e. the exon/CDS structure is rebuilt front-to-back, not inverted."""
    cds_list = _run("miniprot", "-")
    cds_list.reverse()                      # what update_cds_list does for "-"
    assert _spans(cds_list) == ASCENDING
