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
import os
import pickle
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Iterable, Iterator, Optional, Tuple

from lifton.locus_pipeline import LocusResult, StepContext, consume, process_locus


# ---------------------------------------------------------------------------
# Phase 15d — bounded ordered-writer with optional on-disk spill
# ---------------------------------------------------------------------------

class _OrderedWriter:
    """Heap-backed ordered emitter for :class:`LocusResult`. Keeps the
    in-RAM pending set bounded at ``max_pending``; oldest pending
    entries spill to ``spill_dir`` (one pickle file per index) and are
    re-loaded when ``next_to_emit`` advances to their position.

    The output sequence is identical to the unbounded heap path —
    ``offer(...)`` plus a final ``drain()`` reproduces submission
    order regardless of how many entries spilled.
    """

    def __init__(self, *, spill_dir: str, max_pending: int,
                 consume_fn, fw, stats: dict,
                 progress_every: int = 0):
        self.spill_dir = spill_dir
        self.max_pending = max(1, int(max_pending))
        self.consume_fn = consume_fn
        self.fw = fw
        self.stats = stats
        self.progress_every = progress_every
        self.pending: list = []           # min-heap of (idx, LocusResult)
        self.spilled: dict[int, str] = {} # idx -> pickle path
        self.next_to_emit = 0
        self.spill_count = 0
        self._spill_dir_made = False

    def _ensure_spill_dir(self) -> None:
        if not self._spill_dir_made:
            os.makedirs(self.spill_dir, exist_ok=True)
            self._spill_dir_made = True

    def _spill_one(self) -> None:
        """Spill the highest-index entry (least likely to be needed
        soon). Removes it from the in-RAM heap and writes a pickle
        side-car to ``spill_dir``."""
        # heapq is a min-heap; max element is at heap[-1] only when
        # heap is small — to find max we scan, which is O(n) but n is
        # bounded by max_pending so this is cheap.
        if not self.pending:
            return
        max_pos = max(range(len(self.pending)),
                      key=lambda i: self.pending[i][0])
        idx, result = self.pending[max_pos]
        # Remove via swap-with-last + heapify.
        last = self.pending.pop()
        if max_pos < len(self.pending):
            self.pending[max_pos] = last
            heapq.heapify(self.pending)
        self._ensure_spill_dir()
        path = os.path.join(self.spill_dir, f"locus_{idx}.pkl")
        with open(path, "wb") as fh:
            pickle.dump(result, fh, protocol=pickle.HIGHEST_PROTOCOL)
        self.spilled[idx] = path
        self.spill_count += 1

    def _restore(self, idx: int) -> Optional[LocusResult]:
        """Re-load a spilled result from disk; deletes the side-car."""
        path = self.spilled.pop(idx, None)
        if path is None:
            return None
        try:
            with open(path, "rb") as fh:
                return pickle.load(fh)
        finally:
            try:
                os.unlink(path)
            except OSError:
                pass

    def _emit_one(self, result: LocusResult) -> None:
        self.consume_fn(result, self.fw, self.stats)
        self.next_to_emit += 1
        if self.progress_every and \
                self.next_to_emit % self.progress_every == 0:
            sys.stdout.write(
                f"\r>> LiftOn processed: {self.next_to_emit} features."
            )

    def offer(self, result: LocusResult) -> None:
        heapq.heappush(self.pending, (result.index, result))
        # Drain any prefix that just became contiguous.
        self._drain_ready()
        # Bound heap depth.
        while len(self.pending) > self.max_pending:
            self._spill_one()

    def _drain_ready(self) -> None:
        """Pop from heap and (if needed) restore from spill while the
        next-expected index is available."""
        while True:
            if self.pending and self.pending[0][0] == self.next_to_emit:
                _, result = heapq.heappop(self.pending)
                self._emit_one(result)
                continue
            if self.next_to_emit in self.spilled:
                restored = self._restore(self.next_to_emit)
                if restored is not None:
                    self._emit_one(restored)
                    continue
            return

    def drain(self) -> None:
        """Emit every remaining buffered/spilled result in submission
        order. Safe to call when there's nothing left."""
        self._drain_ready()
        # If there's still pending material that hasn't reached
        # next_to_emit (because we spilled future indices), push it
        # all back through the heap and drain.
        while self.pending or self.spilled:
            # Pull the lowest-index from heap or spill, whichever is smaller.
            heap_min = self.pending[0][0] if self.pending else None
            spill_min = min(self.spilled) if self.spilled else None
            candidates = [c for c in (heap_min, spill_min) if c is not None]
            if not candidates:
                return
            target = min(candidates)
            if heap_min is not None and target == heap_min:
                _, result = heapq.heappop(self.pending)
            else:
                result = self._restore(target)
                if result is None:
                    return
            self.next_to_emit = target
            self._emit_one(result)
            # After advancing, see if more contiguous work is unlocked.
            self._drain_ready()


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

    # Phase 15d: bounded ordered-writer with on-disk spill, opt-in via
    # ctx.args.writer_max_pending (default = unbounded ⇒ legacy heap
    # path is byte-identical). When set, the writer caps in-RAM
    # pending entries and spills the rest to a temp directory.
    max_pending = int(getattr(ctx.args, "writer_max_pending", 0) or 0)
    if max_pending > 0:
        spill_dir = tempfile.mkdtemp(prefix="lifton_writer_spill_")
        writer = _OrderedWriter(
            spill_dir=spill_dir,
            max_pending=max_pending,
            consume_fn=consume,
            fw=fw,
            stats=transcripts_stats_dict,
            progress_every=progress_every,
        )
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
                writer.offer(fut.result())
        writer.drain()
        try:
            os.rmdir(spill_dir)
        except OSError:
            pass
        return processed

    # Legacy unbounded heap path — byte-identical to Phase 9 behaviour.
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
