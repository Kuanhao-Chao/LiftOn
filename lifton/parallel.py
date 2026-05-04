# ---------------------------------------------------------------------------
# Phase 9 — deterministic parallel dispatcher for Step 7.
#
# `parallel_step7` runs the per-Liftoff-gene work either serially or
# through a ThreadPoolExecutor, ALWAYS emitting output in submission
# order via a heap-backed ordered-writer buffer. The contract is:
#
#   parallel_step7(..., threads=N)  ==  parallel_step7(..., threads=1)
#
# byte-for-byte. This is what makes `--threads 4` and `--threads 1`
# produce identical output GFF3 — the only externally visible change
# under parallel execution is wall-clock time.
#
# Thread-vs-process choice
# ────────────────────────
# We use ThreadPoolExecutor because both candidate FeatureDB backends
# (gffutils.FeatureDB over SQLite, gffbase.FeatureDB over DuckDB) hold
# connection objects that are not picklable across processes. parasail
# (the alignment hot loop) explicitly releases the GIL during
# `nw_trace_scan_sat`, so threads give real parallel speedup on the
# CPU-heavy alignment work without the IPC + DB-reopen complexity that
# processes would require.
# ---------------------------------------------------------------------------

from __future__ import annotations

import heapq
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Iterable, Iterator, Tuple

from lifton.locus_pipeline import LocusResult, StepContext, consume, process_locus


def _iter_loci(features: Iterable[str], l_feature_db) -> Iterator[Tuple[str, object]]:
    """Yield ``(feature_type, locus)`` pairs in the same order the
    legacy serial loop visits them.

    Iteration order matters: the submission index assigned by
    :func:`enumerate` here becomes the deterministic emit order in the
    parallel path, and the serial baseline is byte-equal to it by
    construction.
    """
    for feature in features:
        for locus in l_feature_db.features_of_type(feature):
            yield feature, locus


def _backend_supports_threads(*dbs, native: bool = False) -> bool:
    """Return True iff per-locus tasks can safely run on worker
    threads given the supplied FeatureDB inputs.

    Phase 11 contract:

      * ``native=True`` AND every supplied FeatureDB is gffbase
        (DuckDB-backed) — the dispatcher pre-materialises every
        locus's children in the parent thread via
        :func:`lifton.locus_pipeline.materialise_locus` and submits
        :class:`MaterialisedLocus` payloads to workers. DuckDB
        tolerates cross-thread reads against the same connection
        when the parent's iterators have been drained first, so
        the guard returns True and Step 7 fans out.
      * ``native=True`` BUT any supplied FeatureDB is gffutils
        (SQLite) — SQLite hard-binds connections to their creator
        thread regardless of pre-materialisation. The guard returns
        False and the dispatcher falls back to serial. Switch to
        the gffbase backend (``--stream --inmemory-liftoff`` or
        ``LIFTON_USE_GFFBASE=1``) to unlock parallelism.
      * ``native=False`` — workers still issue cold
        ``db.children(...)`` reads against the shared connection.
        Returns False (Phase 9 contract preserved).

    The escape hatch ``LIFTON_PARALLEL_FORCE`` overrides everything
    for advanced users / CI experiments that have already verified
    thread-safety on their stack.
    """
    import os as _os
    if _os.environ.get("LIFTON_PARALLEL_FORCE"):
        return True
    if not native:
        return False
    # Phase 11: native is on; unlock threads only when no SQLite
    # connection is in the worker hot path.
    for db in dbs:
        if db is None:
            continue
        module = type(db).__module__ or ""
        if module.startswith("gffutils"):
            return False
    return True


def parallel_step7(
    features: Iterable[str],
    l_feature_db,
    ctx: StepContext,
    fw,
    transcripts_stats_dict: dict,
    *,
    threads: int = 1,
    progress_every: int = 20,
) -> int:
    """Drive Step 7's per-Liftoff-gene work.

    Parameters
    ----------
    features
        Feature types (typically ``["gene"]``) to iterate from
        ``l_feature_db``.
    l_feature_db
        Liftoff FeatureDB object exposing ``features_of_type(...)``.
    ctx
        Immutable context bundle for :func:`process_locus`.
    fw, transcripts_stats_dict
        The output GFF3 file handle and the running stats dict
        (gene → counter); both owned by the parent and only mutated
        from the parent thread by :func:`consume`.
    threads
        Worker count. ``threads <= 1`` falls back to a serial
        for-loop (Phase 5 contract preserved byte-for-byte).
    progress_every
        How often to emit the ``>> LiftOn processed: N features.``
        progress line (matches the legacy loop).

    Returns
    -------
    int
        The number of features that were processed (whether
        emittable or not). Mirrors the legacy ``processed_features``
        counter so the surrounding stats/log code is unchanged.
    """
    submission = enumerate(_iter_loci(features, l_feature_db))
    processed = 0

    # Phase 9/10 thread-safety guard. Workers issue
    # `db.children(...)` reads against the shared FeatureDB
    # connection; both gffutils (SQLite) and gffbase (DuckDB) fail
    # under concurrent iteration unless --native is active (Phase 10
    # bindings sidestep the shared-cursor bottleneck) or
    # LIFTON_PARALLEL_FORCE=1 is set.
    native_active = bool(getattr(ctx.args, "native", False))
    backend_safe = _backend_supports_threads(
        l_feature_db, ctx.ref_db, ctx.m_feature_db,
        native=native_active,
    )
    if not backend_safe and threads is not None and threads > 1:
        sys.stderr.write(
            "\n[LiftOn] --locus-pipeline requested with threads={} but "
            "the active FeatureDB backend cannot serve concurrent "
            "worker reads (no --native, no LIFTON_PARALLEL_FORCE). "
            "Falling back to serial execution. Use --native to unlock "
            "parallelism via mappy/pyminiprot bindings.\n".format(threads)
        )
        threads = 1

    if threads is None or threads <= 1:
        # ── Serial path ──────────────────────────────────────────────
        # Produces output identical to the pre-Phase-9 inline loop in
        # lifton.py (the new path goes through `consume` instead of an
        # inline `if/write_entry`, but the resulting bytes are equal).
        for idx, (_feature, locus) in submission:
            result = process_locus(idx, locus, ctx=ctx)
            consume(result, fw, transcripts_stats_dict)
            processed = idx + 1
            if processed % progress_every == 0:
                sys.stdout.write(
                    f"\r>> LiftOn processed: {processed} features."
                )
        return processed

    # ── Parallel path with ordered-writer buffer ─────────────────────
    # Materialise the locus list in the parent thread BEFORE any
    # worker starts. Both backend iterators (gffutils SQLite cursors,
    # gffbase DuckDB result sets) are not safe to keep open while
    # worker threads issue concurrent reads on the same connection.
    materialised = list(submission)
    processed = len(materialised)

    # Phase 11: when --native is on, additionally pre-materialise
    # each locus's children/ref-attrs in the parent thread so worker
    # tasks consume :class:`MaterialisedLocus` payloads. Workers still
    # delegate to the legacy `run_liftoff.process_liftoff` for the
    # actual gene assembly (the byte-identity contract of Phase 10
    # carries through), but the parent-thread pre-fetch primes the
    # DB caches and removes the cold-iteration aliasing window that
    # gated the Phase 9 thread-safety guard.
    payloads = None
    if native_active:
        from lifton.locus_pipeline import materialise_locus, process_locus_native
        payloads = [
            materialise_locus(idx, locus, ctx)
            for idx, (_feature, locus) in materialised
        ]

    next_to_emit = 0
    pending: list = []  # min-heap keyed by submission index

    with ThreadPoolExecutor(max_workers=int(threads)) as ex:
        if payloads is not None:
            futures = {
                ex.submit(process_locus_native, p, ctx): p.submission_index
                for p in payloads
            }
        else:
            futures = {
                ex.submit(process_locus, idx, locus, ctx=ctx): idx
                for idx, (_feature, locus) in materialised
            }

        for fut in as_completed(futures):
            result: LocusResult = fut.result()
            heapq.heappush(pending, (result.index, result))
            # Drain prefix that is now contiguous from `next_to_emit`.
            while pending and pending[0][0] == next_to_emit:
                _, drained = heapq.heappop(pending)
                consume(drained, fw, transcripts_stats_dict)
                next_to_emit += 1
                if next_to_emit % progress_every == 0:
                    sys.stdout.write(
                        f"\r>> LiftOn processed: {next_to_emit} features."
                    )

    # Defensive: drain anything left (should be empty by construction
    # because `as_completed` enumerates every future).
    while pending:
        _, drained = heapq.heappop(pending)
        consume(drained, fw, transcripts_stats_dict)
        next_to_emit += 1

    return processed
