# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""Merge predicates — pure functions consumed by ``FeatureDB.merge``.

Signature: ``(acc: Feature, cur: Feature, components: list[Feature]) -> bool``.
``acc`` is the running accumulator; ``cur`` is the candidate to fold in;
``components`` is the list of already-folded features. All callables return
True if the candidate should be merged into the accumulator.

Mirrors the legacy `gffutils.merge_criteria` module.
"""

from __future__ import annotations


def seqid(acc, cur, components):
    return acc.seqid == cur.seqid


def strand(acc, cur, components):
    return acc.strand == cur.strand


def feature_type(acc, cur, components):
    return acc.featuretype == cur.featuretype


def exact_coordinates_only(acc, cur, components):
    return acc.start == cur.start and acc.end == cur.end


def overlap_end_inclusive(acc, cur, components):
    """True if ``cur.start`` falls within or immediately after ``acc``."""
    return acc.start <= cur.start <= acc.end + 1


def overlap_start_inclusive(acc, cur, components):
    return acc.start <= cur.end + 1 <= acc.end + 1


def overlap_any_inclusive(acc, cur, components):
    return overlap_end_inclusive(acc, cur, components) or overlap_start_inclusive(
        acc, cur, components
    )


def overlap_end_threshold(threshold: int):
    def predicate(acc, cur, components):
        return abs(acc.end - cur.start) <= threshold
    return predicate


def overlap_start_threshold(threshold: int):
    def predicate(acc, cur, components):
        return abs(acc.start - cur.end) <= threshold
    return predicate


def overlap_any_threshold(threshold: int):
    end_thr = overlap_end_threshold(threshold)
    start_thr = overlap_start_threshold(threshold)
    def predicate(acc, cur, components):
        return end_thr(acc, cur, components) or start_thr(acc, cur, components)
    return predicate
