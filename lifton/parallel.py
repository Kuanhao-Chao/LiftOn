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

    Iteration 8 — the materialisation + proxy-DB path is now universal.
    Parallel Step 7 ALWAYS routes through
    :func:`lifton.locus_pipeline.process_locus_native`: each locus is
    pre-materialised on the parent thread into a
    :class:`MaterialisedLocus`, and workers read it through the
    read-only proxies built by
    :func:`lifton.locus_pipeline._build_proxied_ctx` (``_RefDbProxy``,
    ``_LFeatureDbProxy``, ``_MFeatureDbProxy``). Because workers are
    fully decoupled from the shared connection, SQLite's
    hard-thread-binding no longer applies — gffutils (SQLite) and
    gffbase (DuckDB) backends are equally safe, regardless of
    ``--native``.

    Pre-Iteration-8 this was gated to ``native=True`` only (Phase 17b
    wired the proxy path behind ``--native``); the ``native`` parameter
    is retained for signature/back-compat but **no longer gates** —
    the function returns True for any backend.

    Escape hatches:

      * ``LIFTON_PARALLEL_FORCE`` — overrides everything (force True).
      * ``LIFTON_PARALLEL_BLOCK_GFFUTILS`` — restores the
        pre-Iteration-8 strict rejection for gffutils backends, the
        only remaining way back to serial-on-gffutils (useful if a
        future regression re-introduces a worker-side DB read).
    """
    import os as _os
    if _os.environ.get("LIFTON_PARALLEL_FORCE"):
        return True
    # Iteration 8: workers go through the MaterialisedLocus payload +
    # proxy DBs (process_locus_native); the real FeatureDB connection is
    # never touched from a worker thread, so any backend is safe — this
    # is decoupled from --native. ``native`` is accepted for back-compat.
    if _os.environ.get("LIFTON_PARALLEL_BLOCK_GFFUTILS"):
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

    # Thread-safety guard. Iteration 8: workers never issue reads
    # against the shared FeatureDB connection — they go through the
    # MaterialisedLocus + proxy-DB path (`process_locus_native`), so
    # any backend (gffutils SQLite, gffbase DuckDB) is safe regardless
    # of --native. `native_active` is still computed and passed for
    # signature/back-compat, but no longer gates (only
    # LIFTON_PARALLEL_BLOCK_GFFUTILS can force serial-on-gffutils now).
    native_active = bool(getattr(ctx.args, "native", False))
    backend_safe = _backend_supports_threads(
        l_feature_db, ctx.ref_db, ctx.m_feature_db,
        native=native_active,
    )
    if not backend_safe and threads is not None and threads > 1:
        # Iteration 8: parallel Step 7 routes through the
        # MaterialisedLocus + proxy-DB path on ANY backend, so this
        # branch is no longer reached on a normal run. It now fires
        # ONLY when LIFTON_PARALLEL_BLOCK_GFFUTILS is set — the opt-out
        # that restores the pre-Iteration-8 strict gffutils rejection.
        sys.stderr.write(
            "\n[LiftOn] --locus-pipeline requested with threads={} but "
            "LIFTON_PARALLEL_BLOCK_GFFUTILS is set, so the gffutils "
            "backend is held to serial execution. Falling back to "
            "serial. Unset LIFTON_PARALLEL_BLOCK_GFFUTILS to use the "
            "materialised-payload + proxy-DB parallel path.\n".format(threads)
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

    # Iteration 8: parallel Step 7 ALWAYS pre-materialises each locus
    # into a :class:`MaterialisedLocus` payload so workers consume it
    # through the read-only proxy DBs (`process_locus_native`) and never
    # touch the shared FeatureDB connection. This is what makes the
    # parallel path backend-agnostic AND thread-safe WITHOUT --native
    # (pre-Iteration-8 the materialise step was gated behind
    # `native_active`). Workers still delegate to the legacy
    # `run_liftoff.process_liftoff` via the proxy ctx for gene assembly,
    # so byte-identity carries through (pinned by the 24-cell matrix,
    # whose native=False threaded cells now exercise exactly this path).
    #
    # Phase 17c Item 2b: when the DBs are on disk, the materialise step
    # runs in a small prefetcher pool of 2-4 threads, each holding a
    # thread-local re-opened FeatureDB connection. This shrinks the
    # parent-thread serial bottleneck (~33 s on bee, ~84 s projected
    # on rice, ~5 min on human at 110K transcripts). When the factory
    # cannot reopen thread-local DBs (in-memory blob backends), fall
    # back to the serial parent-thread materialise loop.
    from lifton.locus_pipeline import (
        materialise_locus, process_locus_native,
        _ThreadLocalCtxFactory, materialise_locus_with_factory,
    )
    factory = _ThreadLocalCtxFactory(ctx)
    if factory.viable and len(materialised) > 0:
        # Phase 17c parallel prefetcher pool. Cap at 4 prefetchers
        # since the marginal gain plateaus (~50-60 % reduction at
        # N=4 per the Phase 17 exploration sizing); larger pool
        # increases SQLite contention without commensurate speedup.
        prefetch_workers = min(4, int(threads),
                               max(1, len(materialised)))
        payloads = [None] * len(materialised)
        with ThreadPoolExecutor(
                max_workers=prefetch_workers,
                thread_name_prefix="lifton-prefetch") as _pf_pool:
            _pf_futures = {
                _pf_pool.submit(
                    materialise_locus_with_factory,
                    idx, locus, factory,
                ): idx
                for idx, (_feature, locus) in materialised
            }
            for _fut in as_completed(_pf_futures):
                _p = _fut.result()
                payloads[_p.submission_index] = _p
    else:
        # In-memory backends or non-extractable dbfn — keep the
        # serial parent-thread materialise loop (correct, just no
        # prefetcher speedup).
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
