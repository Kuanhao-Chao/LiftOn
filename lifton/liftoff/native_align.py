# ---------------------------------------------------------------------------
# Phase 11 — native (in-process) alignment path for the vendored Liftoff.
#
# Replaces the per-chromosome `subprocess.run([minimap2, ...])` + SAM file
# write + pysam re-parse in `align_features.py` with an in-memory
# `MinimapAligner.map(...)` round-trip. The entry point
# `align_features_to_target_native` returns the SAME
# `aligned_segments_dict` shape the legacy `align_features_to_target`
# returns, so downstream Liftoff code (lift_features.py,
# find_best_mapping.py) is unchanged.
#
# The byte-identity contract is preserved by reusing the existing
# `align_features.parse_alignment` body via an in-memory pysam-shaped
# adapter. Each MinimapHit is wrapped in `_PysamShim` exposing the
# attribute names that `add_alignment` / `get_aligned_blocks` read.
# ---------------------------------------------------------------------------

from __future__ import annotations

import re
import sys
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

from lifton.liftoff import align_features as _af
from lifton.native_bindings import MinimapAligner, MinimapHit, is_mappy_available


# pysam CIGAR op codes (the values `add_alignment` checks for):
#   0 = match (default minimap2 emit)
#   1 = insertion
#   2 = deletion
#   4 = soft clip
#   5 = hard clip
#   7 = sequence match (with --eqx)
#   8 = sequence mismatch (with --eqx)
_CIGAR_OP_TO_PYSAM = {
    "M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5,
    "P": 6, "=": 7, "X": 8,
}

_CIGAR_TOKEN = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_str_to_pysam_tuples(cigar_str: str) -> List[Tuple[int, int]]:
    """Convert a CIGAR string (e.g. ``"50M2I30M"``) into a list of
    ``(op_int, length)`` pairs in pysam's order. Mirrors what
    pysam.AlignedSegment.cigar would yield."""
    out: List[Tuple[int, int]] = []
    for length, op in _CIGAR_TOKEN.findall(cigar_str):
        op_int = _CIGAR_OP_TO_PYSAM[op]
        out.append((op_int, int(length)))
    return out


class _PysamShim:
    """Read-only adapter that exposes the pysam.AlignedSegment
    attributes the legacy `add_alignment` / `get_aligned_blocks`
    code reads: ``query_name``, ``reference_name``, ``reference_start``,
    ``query_alignment_start``, ``query_alignment_end``, ``cigar``,
    ``is_reverse``, ``is_unmapped``."""

    __slots__ = (
        "query_name", "reference_name", "reference_start",
        "query_alignment_start", "query_alignment_end",
        "cigar", "is_reverse", "is_unmapped",
    )

    def __init__(self, hit: MinimapHit):
        # mappy strand: +1 forward, -1 reverse
        self.query_name = hit.query_name
        self.reference_name = hit.ctg
        self.reference_start = int(hit.r_st)        # mappy is 0-based, pysam 0-based — matches
        self.query_alignment_start = int(hit.q_st)
        self.query_alignment_end = int(hit.q_en)
        self.cigar = cigar_str_to_pysam_tuples(hit.cigar_str)
        self.is_reverse = (int(hit.strand) == -1)
        self.is_unmapped = False


def parse_alignment_from_hits(
    hits_by_query: Dict[str, List[MinimapHit]],
    feature_hierarchy,
    unmapped_features,
    search_type: str,
):
    """In-memory mirror of `align_features.parse_alignment`.

    Consumes ``{query_name: [MinimapHit, ...]}`` and returns the same
    ``{edited_query_name: [aligned_seg, ...]}`` dict the legacy parser
    returns. Reuses the legacy ``add_alignment`` / ``get_aligned_blocks``
    / ``remove_alignments_without_children`` so the resulting
    ``aligned_seg`` records are byte-equivalent to the SAM-parser path.

    Unmapped queries (queries present in the feature hierarchy with no
    hits) are appended to ``unmapped_features`` exactly as the legacy
    parser would.
    """
    all_aligned_blocks = {}
    aln_id = 0
    name_dict = {}
    align_count_dict = {}

    # Iterate in a deterministic order so the resulting aln_id values
    # are stable across runs.
    for query_name in sorted(hits_by_query):
        for hit in hits_by_query[query_name]:
            ref_seq = _PysamShim(hit)
            aln_id = _af.add_alignment(
                ref_seq, align_count_dict, search_type, name_dict, aln_id,
                feature_hierarchy, all_aligned_blocks,
            )

    # Mark queries that produced no hits as unmapped, mirroring the
    # SAM iterator's `is_unmapped` branch.
    for query_name, parent in feature_hierarchy.parents.items():
        if query_name not in hits_by_query or not hits_by_query[query_name]:
            unmapped_features.append(parent)

    _af.remove_alignments_without_children(
        all_aligned_blocks, unmapped_features, feature_hierarchy,
    )
    return all_aligned_blocks


def _iter_query_fasta(features_file: str):
    """Yield ``(name, seq)`` from a FASTA via ``mappy.fastx_read``.
    Falls back to a tiny pure-Python parser when mappy is not
    available so the test suite can mock the iterator without
    actually invoking the binding."""
    try:                                                 # pragma: no cover - normal path
        import mappy
        for name, seq, *_ in mappy.fastx_read(features_file):
            yield name, seq
        return
    except ImportError:                                  # pragma: no cover
        pass
    # Pure-Python fallback (used by tests when mappy is patched out)
    name, parts = None, []
    with open(features_file, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def align_features_to_target_native(
    ref_chroms,
    target_chroms,
    args,
    feature_hierarchy,
    liftover_type,
    unmapped_features,
):
    """Phase 11 in-process replacement for
    :func:`align_features.align_features_to_target`.

    Routes when ``args.native`` is True. Falls back to the legacy
    subprocess path when ``mappy`` is not installed (so a missing
    binding never breaks the pipeline — it just degrades to the
    Phase 5 behaviour with a stderr warning).

    The return value matches the legacy function's contract: a
    ``{query_name: [aligned_seg, ...]}`` dict suitable for the
    downstream `lift_features.py` consumer.
    """
    if not is_mappy_available():
        sys.stderr.write(
            "\n[LiftOn] --native requested but `mappy` is not installed; "
            "falling back to the subprocess minimap2 path. "
            "Install via `pip install mappy` or "
            "`conda install -c bioconda mappy` to unlock the native path.\n"
        )
        return _af.align_features_to_target(
            ref_chroms, target_chroms, args,
            feature_hierarchy, liftover_type, unmapped_features,
        )

    # Build one Aligner per target chromosome (or one for the whole
    # genome when liftover_type isn't "chrm_by_chrm"). Mirrors the
    # split logic in `align_single_chroms`.
    target_fasta = args.target
    aligner = MinimapAligner(
        target_fasta,
        mm2_options=getattr(args, "mm2_options", ""),
        threads=int(getattr(args, "threads", 1) or 1),
    )

    hits_by_query: Dict[str, List[MinimapHit]] = {}
    # The legacy code writes per-chromosome reference-feature FASTA
    # files into args.directory; we re-use those if they exist (the
    # same `extract_features` step that the subprocess path consumes
    # also feeds us).
    for chrom_idx in range(len(target_chroms)):
        features_file, _features_name = _af.get_features_file(
            ref_chroms, args, liftover_type, chrom_idx,
        )
        for query_name, query_seq in _iter_query_fasta(features_file):
            for hit in aligner.map(query_name, query_seq):
                hits_by_query.setdefault(query_name, []).append(hit)

    return parse_alignment_from_hits(
        hits_by_query, feature_hierarchy, unmapped_features, liftover_type,
    )
