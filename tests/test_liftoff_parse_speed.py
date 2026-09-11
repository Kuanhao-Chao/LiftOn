"""Faster SAM-to-block parsing in the vendored Liftoff (v1.0.12).

``get_aligned_blocks`` now computes the children's parent-relative spans once
per alignment and tests block overlap with a bisect, building an aligned
segment only for blocks it keeps. These tests compare it against the original
algorithm (reconstructed here from the retained ``add_block`` and
``find_overlapping_children``) on thousands of random alignments, both
strands, with and without hard clips.
"""
from __future__ import annotations

import random
import types

import numpy as np
import pytest

from lifton.liftoff import align_features, liftoff_utils

OPS = align_features.get_cigar_operations()


def _original_blocks(alignment, aln_id, feature_hierarchy, search_type):
    """The pre-v1.0.12 ``get_aligned_blocks`` body."""
    cigar = alignment.cigar
    parent_key = liftoff_utils.convert_id_to_original(alignment.query_name)
    parent = feature_hierarchy.parents[parent_key]
    query_start, query_end = align_features.get_query_start_and_end(alignment, cigar, OPS)
    children = feature_hierarchy.children[parent_key]
    end_to_end = align_features.is_end_to_end_alignment(parent, query_start, query_end)
    if search_type == "copies" and end_to_end is False:
        return []
    reference_block_start = reference_block_pos = alignment.reference_start
    query_block_start = query_block_pos = query_start
    new_blocks, mismatches = [], []
    merged = liftoff_utils.merge_children_intervals(children)
    for operation, length in cigar:
        if align_features.base_is_aligned(operation, OPS):
            if operation == OPS["mismatch"]:
                for i in range(query_block_pos, query_block_pos + length):
                    mismatches.append(i)
            query_block_pos, reference_block_pos = align_features.adjust_position(
                operation, query_block_pos, reference_block_pos, length, OPS)
            if query_block_pos == query_end:
                align_features.add_block(query_block_pos, reference_block_pos, aln_id, alignment,
                                         query_block_start, reference_block_start, mismatches,
                                         new_blocks, merged, parent)
                break
        elif align_features.is_alignment_gap(operation, OPS):
            align_features.add_block(query_block_pos, reference_block_pos, aln_id, alignment,
                                     query_block_start, reference_block_start, mismatches,
                                     new_blocks, merged, parent)
            mismatches, query_block_start, reference_block_start, query_block_pos, reference_block_pos = \
                align_features.end_block_at_gap(operation, query_block_pos, reference_block_pos,
                                                length, OPS)
    return new_blocks


def _random_case(rng):
    parent_start = rng.randrange(1, 5000)
    parent_length = rng.randrange(50, 3000)
    parent = types.SimpleNamespace(id="g", start=parent_start,
                                   end=parent_start + parent_length - 1)
    children = []
    for _ in range(rng.randrange(0, 8)):
        start = rng.randrange(parent.start, parent.end + 1)
        end = min(parent.end, start + rng.randrange(0, 400))
        children.append(types.SimpleNamespace(start=start, end=end))
    cigar = []
    if rng.random() < 0.3:
        cigar.append((OPS["hard_clip"], rng.randrange(1, 30)))
    aligned = 0
    for _ in range(rng.randrange(1, 40)):
        operation = rng.choice([OPS["match"]] * 5 + [OPS["mismatch"]] * 2
                               + [OPS["insertion"], OPS["deletion"]])
        length = rng.randrange(1, 60)
        cigar.append((operation, length))
        if operation in (OPS["match"], OPS["mismatch"], OPS["insertion"]):
            aligned += length
    query_alignment_start = rng.randrange(0, 20)
    alignment = types.SimpleNamespace(
        query_name="g_0", reference_name="chr1", cigar=cigar,
        reference_start=rng.randrange(0, 10000), is_reverse=rng.random() < 0.5,
        query_alignment_start=query_alignment_start,
        query_alignment_end=query_alignment_start + aligned)
    hierarchy = types.SimpleNamespace(parents={"g": parent}, children={"g": children})
    return alignment, hierarchy


def _as_tuples(blocks):
    return [(b.aln_id, b.query_name, b.reference_name, b.query_block_start,
             b.query_block_end, b.reference_block_start, b.reference_block_end,
             b.is_reverse, b.mismatches.dtype, b.mismatches.tolist()) for b in blocks]


@pytest.mark.parametrize("search_type", ["chrm_by_chrm", "copies"])
def test_blocks_match_the_original_algorithm(search_type):
    rng = random.Random(2026)
    for case in range(3000):
        alignment, hierarchy = _random_case(rng)
        expected = _original_blocks(alignment, case, hierarchy, search_type)
        actual = align_features.get_aligned_blocks(alignment, case, hierarchy, search_type)
        assert _as_tuples(actual) == _as_tuples(expected), case


def test_overlap_test_matches_per_child_scan():
    rng = random.Random(11)
    for _ in range(5000):
        alignment, hierarchy = _random_case(rng)
        parent = hierarchy.parents["g"]
        merged = liftoff_utils.merge_children_intervals(hierarchy.children["g"])
        spans = align_features._relative_child_spans(parent, merged, alignment.is_reverse)
        query_start = rng.randrange(-50, 3100)
        query_end = query_start + rng.randrange(0, 500)
        block = types.SimpleNamespace(query_block_start=query_start,
                                      query_block_end=query_end,
                                      is_reverse=alignment.is_reverse)
        expected = align_features.find_overlapping_children(block, merged, parent) != []
        assert align_features._overlaps_any_child(spans, query_start, query_end) is expected


def test_mismatch_array_type_is_unchanged():
    alignment = types.SimpleNamespace(
        query_name="g_0", reference_name="chr1",
        cigar=[(OPS["match"], 5), (OPS["mismatch"], 2), (OPS["match"], 3)],
        reference_start=100, is_reverse=False,
        query_alignment_start=0, query_alignment_end=10)
    parent = types.SimpleNamespace(id="g", start=1, end=10)
    hierarchy = types.SimpleNamespace(
        parents={"g": parent},
        children={"g": [types.SimpleNamespace(start=1, end=10)]})
    (block,) = align_features.get_aligned_blocks(alignment, 1, hierarchy, "chrm_by_chrm")
    assert block.mismatches.tolist() == [5, 6]
    assert block.mismatches.dtype == np.array([1]).astype(int).dtype
