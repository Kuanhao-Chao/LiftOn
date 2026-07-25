"""Tier-3 perf fix — committed logical intervals must be pruned.

`Step7StateCoordinator._logical_intervals` makes an allocated-but-not-yet-committed
gene visible to `overlap()`. It was never pruned, so the list grew for the whole run and
every `overlap()` call -- two per (transcript, miniprot-candidate) pair -- rescanned it
linearly while holding the single coordinator lock. That is O(n^2) on a gene-dense
chromosome and silently caps the `--threads N` speedup.

Once `commit_locus_delta` has added the interval to the real `tree_dict`, the logical
copy is a pure duplicate (both built by the same `_make_interval`, and `overlap()` unions
into a set), so pruning it cannot change a query result. Entries whose locus FAILED to
serialize are never committed and must be retained.
"""
from types import SimpleNamespace

from lifton.locus_pipeline import LocusDelta, Step7StateCoordinator


def _coordinator():
    return Step7StateCoordinator(ref_features_dict={}, tree_dict={})


def _gene(gene_id, start, end, seqid="chr1"):
    return SimpleNamespace(id=gene_id, seqid=seqid, start=start, end=end,
                           featuretype="gene", attributes={})


def _allocate(coord, index, gene_id, start, end):
    delta = LocusDelta()
    coord.allocate_gene(index, gene_id, _gene(gene_id, start, end), delta)
    return delta


def test_commit_prunes_the_logical_placeholder():
    coord = _coordinator()
    delta = _allocate(coord, 0, "g1", 100, 200)
    assert len(coord._logical_intervals["chr1"]) == 1

    coord.commit(delta)

    # Pruned from the logical list...
    assert "chr1" not in coord._logical_intervals
    # ...but still found, now via the real tree.
    assert coord.overlap("chr1", 150, 160, through_index=0)


def test_overlap_result_is_unchanged_by_pruning():
    coord = _coordinator()
    delta = _allocate(coord, 0, "g1", 100, 200)
    before = coord.overlap("chr1", 150, 160, through_index=0)

    coord.commit(delta)
    after = coord.overlap("chr1", 150, 160, through_index=0)

    assert before == after, (before, after)


def test_uncommitted_locus_keeps_its_placeholder():
    """A locus that failed to serialize never commits, so its interval must stay."""
    coord = _coordinator()
    committed = _allocate(coord, 0, "g1", 100, 200)
    _failed = _allocate(coord, 1, "g2", 300, 400)

    coord.commit(committed)          # only the first one commits

    rows = coord._logical_intervals.get("chr1", [])
    assert [row[1].data for row in rows] == ["g2"]
    assert coord.overlap("chr1", 350, 360, through_index=1)


def test_logical_list_does_not_grow_across_many_commits():
    """The whole point: memory and per-query cost stay bounded."""
    coord = _coordinator()
    for index in range(50):
        start = 100 + index * 1000
        coord.commit(_allocate(coord, index, f"g{index}", start, start + 500))

    assert coord._logical_intervals.get("chr1", []) == []
    # every committed gene is still discoverable
    assert coord.overlap("chr1", 100, 600, through_index=49)
    assert coord.overlap("chr1", 49100, 49600, through_index=49)


def test_seqid_visibility_is_preserved_after_pruning():
    coord = _coordinator()
    delta = _allocate(coord, 0, "g1", 100, 200)
    coord.commit(delta)

    assert coord.has_seqid("chr1", through_index=0)
    assert "chr1" in coord.seqids(through_index=0)
