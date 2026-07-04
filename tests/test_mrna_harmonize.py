"""GitHub issue #28 — harmonize transcript vs mRNA in the output.

LiftOn preserved the reference featuretype on the Liftoff / chaining path (often the
generic ``transcript``) while the miniprot path emitted ``mRNA``, so one output
labelled the same kind of feature two ways. A transcript that HAS a CDS is an mRNA
by SO definition, so ``feature_serializer.write_trans`` now canonicalizes coding
transcripts to ``mRNA`` (non-coding transcripts keep their type).
``LIFTON_NO_MRNA_HARMONIZE=1`` preserves the reference featuretype.
"""
import io

import gffutils

from lifton.lifton_class import Lifton_TRANS, Lifton_EXON, Lifton_CDS
from lifton.io import feature_serializer


def _ent(ftype, s, e, fid, parent):
    f = gffutils.Feature(seqid="chr1", source="LiftOn", featuretype=ftype, start=s,
                         end=e, strand="+", frame=".",
                         attributes={"ID": [fid], "Parent": [parent]})
    f.id = fid
    return f


def _trans(ftype, coding):
    trans = Lifton_TRANS.__new__(Lifton_TRANS)
    trans.entry = _ent(ftype, 101, 399, "tx1", "gene1")
    ex = Lifton_EXON(_ent("exon", 101, 199, "exon1", "tx1"))
    ex.cds = Lifton_CDS(_ent("CDS", 101, 199, "c", "tx1")) if coding else None
    trans.exons = [ex]
    return trans


def _transcript_line_type(trans):
    fw = io.StringIO()
    feature_serializer.write_trans(trans, fw)
    return fw.getvalue().splitlines()[0].split("\t")[2]


def test_coding_transcript_harmonized_to_mrna():
    assert _transcript_line_type(_trans("transcript", coding=True)) == "mRNA"


def test_already_mrna_unchanged():
    assert _transcript_line_type(_trans("mRNA", coding=True)) == "mRNA"


def test_noncoding_transcript_keeps_its_type():
    assert _transcript_line_type(_trans("ncRNA", coding=False)) == "ncRNA"


def test_harmonize_disabled_preserves_reference_type(monkeypatch):
    monkeypatch.setenv("LIFTON_NO_MRNA_HARMONIZE", "1")
    assert _transcript_line_type(_trans("transcript", coding=True)) == "transcript"
