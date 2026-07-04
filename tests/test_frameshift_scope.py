"""GitHub issue #46 — false 'frameshift' labels from UTR indels.

variants.find_variants ran is_frameshift on the FULL-transcript DNA alignment
(5'UTR + CDS + 3'UTR). A 1-2 bp indel in a UTR (ubiquitous on close pairs such as
human<->chimp) is a gap-run not divisible by three, so it was mislabelled
'frameshift' in score.txt AND triggered a needless ORF re-search. The fix scopes
the frameshift test to the coding sub-alignment (the CDS columns), via a CDS span
the caller derives from the contiguous CDS block in the spliced transcript.
"""
from types import SimpleNamespace

from lifton import variants, lifton_class


def _dna(query_aln, ref_aln, identity=0.9):
    return SimpleNamespace(query_aln=query_aln, ref_aln=ref_aln, identity=identity,
                           query_seq=query_aln.replace("-", ""))


def _prot(identity=0.9):
    return SimpleNamespace(identity=identity, query_aln="MK", ref_aln="MK",
                           query_seq="MK*")


def _status_types(dna, cds_span):
    status = lifton_class.Lifton_Status()
    variants.find_variants(dna, _prot(), status, ["MK", ""], is_non_coding=False,
                           cds_span=cds_span)
    return status.status


# ---------------- the scoping helper ----------------
def test_coding_subalignment_keeps_only_cds_columns():
    dna = _dna("A-AATGCGT", "AGAATGCGT")   # 1-bp query deletion in the 5'UTR (qpos 1)
    q, r = variants._coding_subalignment(dna, 3, 8)
    assert "-" not in q and len(q) == 5    # only the coding columns, gap-free


# ---------------- the fix ----------------
def test_utr_indel_not_labelled_frameshift():
    # deletion at query pos 1 (5'UTR); CDS is [3, 8) -> scoping excludes the UTR gap.
    types = _status_types(_dna("A-AATGCGT", "AGAATGCGT"), cds_span=(3, 8))
    assert "frameshift" not in types


def test_cds_indel_still_labelled_frameshift():
    # deletion at query pos 5 (INSIDE the CDS [3, 9)) -> still a frameshift.
    types = _status_types(_dna("AAATG-CGT", "AAATGACGT"), cds_span=(3, 9))
    assert "frameshift" in types


def test_without_cds_span_falls_back_to_full_alignment():
    # No CDS span -> unscoped (old behaviour): the UTR indel IS flagged. Proves the
    # change is backward-compatible when a caller does not supply the span.
    types = _status_types(_dna("A-AATGCGT", "AGAATGCGT"), cds_span=None)
    assert "frameshift" in types
