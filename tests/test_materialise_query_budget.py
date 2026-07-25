"""Tier-3 perf fix — the scalar materialiser must not walk exon/CDS leaves.

`_walk_and_cache_features` (the gffutils path used by every `--threads N
--locus-pipeline` run) ran ALL FOUR `children(...)` query variants for EVERY feature and
then recursed into `children_l1` -- which for a transcript is its exons and CDS. Each
leaf therefore paid ~4 more SQLite round-trips and got a cache entry of its own, for
signatures the runtime never asks a leaf about. On a human RefSeq lift that is millions
of redundant queries, and it is the dominant cost of the materialise phase.

Its batched (gffbase) sibling `_walk_and_cache_features_batched` already keeps the
correct container/terminal split and never places leaves on the frontier. These tests
pin (a) that the two walkers now agree, and (b) the query budget, which is the
deterministic mechanism proof behind the speedup.
"""
from types import SimpleNamespace

from lifton.locus_pipeline import (
    MaterialisedLocus, _walk_and_cache_features,
)


def _feat(fid, ftype, start, end):
    return SimpleNamespace(id=fid, featuretype=ftype, start=start, end=end,
                           seqid="chr1", strand="+")


GENE = _feat("gene1", "gene", 100, 400)
TX = _feat("tx1", "mRNA", 100, 400)
EXONS = [_feat("e1", "exon", 100, 199), _feat("e2", "exon", 300, 400)]
CDSS = [_feat("c1", "CDS", 100, 199), _feat("c2", "CDS", 300, 400)]


class _CountingDB:
    """A realistic gene -> mRNA -> exon/CDS hierarchy that counts every query."""

    def __init__(self):
        self.calls = []

    def children(self, feature, featuretype=None, level=None, order_by=None):
        fid = getattr(feature, "id", feature)
        self.calls.append((fid, featuretype, level))
        if fid == "gene1":
            return iter([TX]) if (featuretype is None and level == 1) else iter([])
        if fid == "tx1":
            if featuretype == "exon":
                return iter(EXONS)
            if featuretype == ("CDS", "stop_codon"):
                return iter(CDSS)
            if featuretype is None and level == 1:
                return iter(EXONS + CDSS)
            return iter([])
        return iter([])          # exon/CDS leaves have no children


def _walk():
    db = _CountingDB()
    ctx = SimpleNamespace(l_feature_db=db)
    payload = MaterialisedLocus(submission_index=0, locus=GENE, locus_id="gene1")
    _walk_and_cache_features(GENE, ctx, payload)
    return db, payload


def test_leaves_are_never_queried_or_cached():
    db, payload = _walk()
    queried = {fid for fid, _ft, _lvl in db.calls}
    assert queried == {"gene1", "tx1"}, f"walked leaves: {queried}"
    assert set(payload.feature_cache) == {"gene1", "tx1"}


def test_query_budget_is_one_per_container_plus_two_per_transcript():
    db, _payload = _walk()
    # gene1 (container): exon@l1 (classify) + children@l1.
    # tx1  (terminal):   exon@l1 (classify) + CDS/stop.  No level-1 query, no recursion.
    assert len(db.calls) == 4, db.calls


def test_transcript_cache_still_serves_every_signature_the_runtime_uses():
    _db, payload = _walk()
    tx = payload.feature_cache["tx1"]
    assert [e.id for e in tx.exon_children_l1] == ["e1", "e2"]
    # Variant 3 (no-level exons) equals Variant 1 under the leaf-exon convention.
    assert [e.id for e in tx.exon_children_full] == ["e1", "e2"]
    assert [c.id for c in tx.cds_stop_children] == ["c1", "c2"]


def test_container_still_exposes_its_level1_children_for_the_recursion():
    _db, payload = _walk()
    gene = payload.feature_cache["gene1"]
    assert [c.id for c in gene.children_l1] == ["tx1"]
    assert gene.exon_children_l1 == []
