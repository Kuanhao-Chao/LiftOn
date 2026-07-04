"""Regression test for GitHub issue #39 — KeyError crash on an unknown parent.

On GRCz11 zebrafish, a copy/fragment query name from a gene family (the reporter's
`gene-trim`) converted via `convert_id_to_original` to an id that is not a key in
`feature_hierarchy.parents`. The three `feature_hierarchy.parents[...]` /
`children[...]` lookups in align_features indexed it unconditionally, so ONE such
alignment raised `KeyError: 'gene-trim'` and crashed the WHOLE chromosome (only
chr22 failed to lift while every other chromosome succeeded). The fix routes every
lookup through `_known_parent_key`, which returns None (warning once) for an
unknown id so that alignment is skipped and the rest of the chromosome still lifts.
"""
from types import SimpleNamespace

from lifton.liftoff import align_features


def test_known_parent_key_present_absent_and_warns_once(capsys):
    align_features._WARNED_MISSING_PARENT.clear()
    fh = SimpleNamespace(parents={"gene-real": 1})
    assert align_features._known_parent_key(fh, "gene-real_0") == "gene-real"
    assert align_features._known_parent_key(fh, "gene-trim_0") is None
    assert align_features._known_parent_key(fh, "gene-trim_0") is None   # no re-warn
    err = capsys.readouterr().err
    assert err.count("GitHub issue #39") == 1
    assert "gene-trim" in err


def test_get_aligned_blocks_skips_unknown_parent():
    """An aligned query whose converted id is not a known parent returns [] instead
    of raising KeyError (the #39 crash site, align_features.py get_aligned_blocks)."""
    align_features._WARNED_MISSING_PARENT.clear()
    fh = SimpleNamespace(parents={"gene-real": object()}, children={"gene-real": []})
    aln = SimpleNamespace(query_name="gene-trim_0", cigar=[])
    assert align_features.get_aligned_blocks(aln, 1, fh, "copies") == []


def test_remove_alignments_without_children_tolerates_unknown_parent():
    align_features._WARNED_MISSING_PARENT.clear()
    fh = SimpleNamespace(parents={"gene-real": "P"}, children={})
    unmapped = []
    blocks = {"gene-trim_0": [], "gene-real_0": []}
    out = align_features.remove_alignments_without_children(blocks, unmapped, fh)
    assert out == {}                 # both empty-block entries removed, no crash
    assert unmapped == ["P"]         # only the KNOWN parent recorded as unmapped
