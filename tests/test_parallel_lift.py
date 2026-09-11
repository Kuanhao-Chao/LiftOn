"""Process-parallel Liftoff lift loop (v1.0.12 S2).

The real equivalence proof is a fresh Liftoff lift compared byte for byte
(``benchmarks``); the test suite cannot run minimap2. These tests pin the
orchestration that makes the parallel loop equal to the serial one: batching
by reference chromosome in serial order, neighbour lookups that see this
batch plus earlier passes, results merged in serial insertion order, unmapped
parents returned as the original objects, and the serial fallbacks.
"""
from __future__ import annotations

import types

import pytest

from lifton.liftoff import lift_features, parallel_lift


class _Parent:
    def __init__(self, identifier, seqid, start):
        self.id, self.seqid, self.start = identifier, seqid, start


class _Segment:
    def __init__(self, query_name):
        self.query_name = query_name


class _FeatureDb:
    def featuretypes(self):
        return iter(["gene", "mRNA", "exon", "CDS"])


# Genes on two reference chromosomes; g4 fails to lift.
PARENTS = {
    "g1": _Parent("g1", "chr1", 100), "g2": _Parent("g2", "chr1", 200),
    "g3": _Parent("g3", "chr2", 50), "g4": _Parent("g4", "chr2", 150),
    "g5": _Parent("g5", "chr2", 300), "g6": _Parent("g6", "chr3", 10),
}
ORDER = ["g1", "g2", "g3", "g4", "g5", "g6"]


def _alignments(names=ORDER):
    # Deliberately out of order: the loop must sort by reference position.
    return {name: [_Segment(f"{name}_0")] for name in reversed(names)}


def _fake_neighbour(ref_parents, alignment, lifted, ref_parent_order):
    """Upstream same-chromosome neighbour, as the real helper computes it."""
    name = alignment[0].query_name[:-2]
    index = ORDER.index(name)
    if index == 0 or PARENTS[ORDER[index - 1]].seqid != PARENTS[name].seqid:
        return 0, "", 0
    key = ORDER[index - 1] + "_0"
    if key in lifted:
        return lifted[key][0], "seen", 0
    return 0, "", 0


def _fake_lift(threshold, feature_order, features_to_lift, feature_hierarchy,
               prev_start, prev_ref_start, prev_seq, unmapped, alignment,
               seq_id_threshold, feature_locations, lifted_list, args):
    name = alignment[0].query_name
    parent = features_to_lift[name[:-2]]
    if parent.id == "g4":
        unmapped.append(parent)
        return [], name
    # Record what the neighbour hint saw, so any divergence shows in the output.
    return [f"{name}|after={prev_start}|{prev_seq}"], name


@pytest.fixture
def patched(monkeypatch):
    monkeypatch.setattr(lift_features, "find_neighbor_location", _fake_neighbour)
    monkeypatch.setattr(lift_features, "lift_single_feature", _fake_lift)
    monkeypatch.delenv("LIFTON_PARALLEL_LIFT", raising=False)
    monkeypatch.delenv("LIFTON_PARALLEL_LIFT_WORKERS", raising=False)


def _run(parallel, *, earlier=None, names=ORDER, threads=4, feature_locations=None):
    hierarchy = types.SimpleNamespace(parents=dict(PARENTS))
    lifted = dict(earlier or {})
    unmapped = []
    args = types.SimpleNamespace(threads=threads, parallel_lift=parallel)
    lift_features.lift_all_features(
        _alignments(names), 0.5, _FeatureDb(), hierarchy, unmapped, lifted, 0.5,
        feature_locations, args, ref_parent_order=None)
    return lifted, unmapped, hierarchy


def test_parallel_equals_serial(patched):
    serial, serial_unmapped, _ = _run(False)
    parallel, parallel_unmapped, hierarchy = _run(True)
    assert list(parallel.items()) == list(serial.items())
    assert [f.id for f in parallel_unmapped] == [f.id for f in serial_unmapped] == ["g4"]
    # The unmapped parent is the original object, not a pickled copy.
    assert parallel_unmapped[0] is hierarchy.parents["g4"]
    # Neighbour hints worked inside a chromosome and stopped at its boundary.
    assert serial["g2_0"] == ["g2_0|after=g1_0|after=0||seen"]
    assert serial["g3_0"] == ["g3_0|after=0|"]


def test_neighbours_from_an_earlier_pass_are_visible(patched):
    # A later pass (the unmapped-genes pass) re-lifts only some genes; its
    # neighbour hints must still see what the earlier pass lifted.
    earlier = {"g1_0": ["from-primary-pass"]}
    names = ["g2", "g3", "g5", "g6"]
    serial, _, _ = _run(False, earlier=earlier, names=names)
    parallel, _, _ = _run(True, earlier=earlier, names=names)
    assert list(parallel.items()) == list(serial.items())
    assert serial["g2_0"] == ["g2_0|after=from-primary-pass|seen"]


def test_batches_follow_reference_chromosomes():
    alignments = lift_features.sort_alignments(PARENTS, _alignments())
    batches = parallel_lift.group_by_reference_chromosome(alignments, PARENTS)
    assert [[a[0].query_name for a in batch] for batch in batches] == [
        ["g1_0", "g2_0"], ["g3_0", "g4_0", "g5_0"], ["g6_0"]]


@pytest.mark.parametrize("threads, feature_locations", [
    (1, None),                 # one worker: nothing to parallelize
    (4, object()),             # overlap intervals: the loop is not independent
])
def test_serial_fallbacks(patched, monkeypatch, threads, feature_locations):
    monkeypatch.setattr(parallel_lift.multiprocessing, "get_context",
                        lambda *_: pytest.fail("no pool expected"))
    assert parallel_lift.lift_all_features_parallel(
        _alignments(), 0.5, _FeatureDb(),
        types.SimpleNamespace(parents=dict(PARENTS)), [], {}, 0.5,
        feature_locations, types.SimpleNamespace(threads=threads), None) is False


def test_switch_resolution(monkeypatch):
    monkeypatch.delenv("LIFTON_PARALLEL_LIFT", raising=False)
    assert parallel_lift.enabled(types.SimpleNamespace(parallel_lift=True))
    assert not parallel_lift.enabled(types.SimpleNamespace(parallel_lift=False))
    assert parallel_lift.enabled(types.SimpleNamespace()) is \
        parallel_lift.PARALLEL_LIFT_DEFAULT
    monkeypatch.setenv("LIFTON_PARALLEL_LIFT", "1")
    assert parallel_lift.enabled(types.SimpleNamespace(parallel_lift=False))
