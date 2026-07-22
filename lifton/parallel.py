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
import threading
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from typing import Callable, Iterable, Iterator, Optional, Tuple

from lifton.locus_pipeline import (
    LocusResult,
    Step7StateCoordinator,
    StepContext,
    consume,
    make_worker_context,
    process_locus,
)


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

    def __init__(self, *, spill_dir: str, max_pending: Optional[int],
                 consume_fn, fw, stats: dict,
                 progress_every: int = 0, progress_stream=None):
        self.spill_dir = spill_dir
        self.max_pending = (max(1, int(max_pending))
                            if max_pending is not None else None)
        self.consume_fn = consume_fn
        self.fw = fw
        self.stats = stats
        self.progress_every = progress_every
        self.progress_stream = (sys.stdout if progress_stream is None
                                else progress_stream)
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
            self.progress_stream.write(
                f"\r>> LiftOn processed: {self.next_to_emit} features."
            )

    def offer(self, result: LocusResult) -> None:
        heapq.heappush(self.pending, (result.index, result))
        # Drain any prefix that just became contiguous.
        self._drain_ready()
        # Bound heap depth.
        while (self.max_pending is not None and
               len(self.pending) > self.max_pending):
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


def _iter_loci_serial_safe(
    features: Iterable[str], l_feature_db,
) -> Iterator[Tuple[str, object]]:
    """Root scan safe to interleave with hierarchy reads on one handle."""
    safe_scan = getattr(l_feature_db, "features_of_type_stream_safe", None)
    for feature in features:
        loci = (safe_scan(feature) if callable(safe_scan)
                else l_feature_db.features_of_type(feature))
        for locus in loci:
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


def _step7_inflight_limit(args, threads: int) -> int:
    """Return the maximum number of submitted-but-not-emitted loci.

    The ordered commit boundary, not merely future completion, drives this
    limit.  Consequently a slow early locus cannot let later completed results
    accumulate without bound.  ``LIFTON_STEP7_MAX_INFLIGHT`` is primarily an
    operational escape hatch; embedders may set ``args.step7_max_inflight``
    directly without requiring a CLI-only option.
    """
    for raw in (getattr(args, "step7_max_inflight", None),
                os.environ.get("LIFTON_STEP7_MAX_INFLIGHT")):
        if not isinstance(raw, (int, str)) or isinstance(raw, bool):
            continue
        try:
            configured = int(raw)
        except (TypeError, ValueError):
            continue
        if configured > 0:
            return configured
    return max(1, 2 * int(threads))


def _ordered_bounded_map(
    executor,
    items,
    submit_fn,
    max_inflight: int,
    *,
    high_water_callback: Optional[Callable[[int], None]] = None,
):
    """Yield submitted work in input order while bounding the whole window.

    ``Executor.map`` eagerly consumes its input on older Python releases and
    limiting only unfinished futures still permits completed later results to
    accumulate behind a slow first item.  This primitive counts both states
    against one limit.  It is shared by Step 7's materialisation prefetch and
    locus-processing stages so their memory contracts cannot drift apart.
    """
    source = iter(items)
    futures = {}
    buffered = {}
    submitted = 0
    next_to_yield = 0
    exhausted = False

    def report_high_water():
        if high_water_callback is not None:
            high_water_callback(len(futures) + len(buffered))

    def fill_window():
        nonlocal submitted, exhausted
        while not exhausted and len(futures) + len(buffered) < max_inflight:
            try:
                item = next(source)
            except StopIteration:
                exhausted = True
                break
            future = submit_fn(executor, item)
            futures[future] = submitted
            submitted += 1
            report_high_water()

    fill_window()
    while futures or buffered:
        while next_to_yield in buffered:
            result = buffered.pop(next_to_yield)
            next_to_yield += 1
            yield result
            fill_window()
        if not futures:
            continue
        done, _ = wait(tuple(futures), return_when=FIRST_COMPLETED)
        for future in done:
            ordinal = futures.pop(future)
            buffered[ordinal] = future.result()
            report_high_water()


def _iter_bounded_ordered(executor, items, submit_fn,
                          max_inflight: int):
    """Backward-compatible alias for the shared bounded-map primitive."""
    yield from _ordered_bounded_map(
        executor, items, submit_fn, max_inflight,
    )


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
    processed = 0
    progress_stream = (sys.stderr
                       if getattr(ctx.args, "output", None) == "stdout"
                       else sys.stdout)
    state_coordinator = Step7StateCoordinator(
        ctx.ref_features_dict, ctx.tree_dict,
    )
    # Callers occasionally reuse an argparse namespace across focused runs.
    # Clear prior runtime diagnostics before recording this dispatch.
    for marker in (
        "_step7_parallel_fallback_reason",
        "_step7_parallel_strategy",
        "_step7_materialisation_constraint",
    ):
        try:
            delattr(ctx.args, marker)
        except (AttributeError, TypeError):
            pass

    def _consume_ordered(result, output_handle, stats):
        return consume(
            result, output_handle, stats, ctx=ctx,
            state_coordinator=state_coordinator,
        )

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
        try:
            setattr(
                ctx.args,
                "_step7_parallel_fallback_reason",
                "LIFTON_PARALLEL_BLOCK_GFFUTILS requested serial execution",
            )
        except Exception:
            pass
        threads = 1

    if threads is None or threads <= 1:
        # ── Serial path ──────────────────────────────────────────────
        # Produces output identical to the pre-Phase-9 inline loop in
        # lifton.py (the new path goes through `consume` instead of an
        # inline `if/write_entry`, but the resulting bytes are equal).
        for idx, (_feature, locus) in enumerate(
                _iter_loci_serial_safe(features, l_feature_db)):
            worker_ctx = make_worker_context(ctx, idx, state_coordinator)
            result = process_locus(idx, locus, ctx=worker_ctx)
            _consume_ordered(result, fw, transcripts_stats_dict)
            processed = idx + 1
            if progress_every and processed % progress_every == 0:
                progress_stream.write(
                    f"\r>> LiftOn processed: {processed} features."
                )
        try:
            setattr(
                ctx.args,
                "_step7_max_inflight_observed",
                1 if processed else 0,
            )
            setattr(ctx.args, "_step7_parallel_strategy", "serial")
        except Exception:
            pass
        return processed

    # ── Parallel path with ordered-writer buffer ─────────────────────
    # Roots come from a dedicated reopened database, so worker hierarchy
    # queries cannot invalidate the scan cursor. Roots and full payloads are
    # consumed lazily within the same bounded dispatch window.
    # Iteration 10: parallel Step 7 FUSES materialise + process into one
    # pipelined pool. Each worker materialises ITS OWN locus (reads from a
    # thread-local re-opened DB via `_ThreadLocalCtxFactory`) and then
    # immediately processes it (`process_locus_native`, which reads only the
    # materialised payload's in-process caches + the parent ctx's non-DB
    # fields). The SQLite-bound materialise on one locus overlaps the
    # GIL-released parasail on others, so wall collapses from the
    # pre-Iteration-10 `T_materialise + T_process` (a 4-thread prefetcher pool
    # drained FULLY before the worker pool started — a hard barrier) toward
    # `max(T_materialise, T_process)`.
    #
    # Byte-identity carries through unchanged: the fused worker builds the
    # SAME `MaterialisedLocus` payload (same locus, same on-disk DB,
    # `order_by='start'`) the prefetcher pool built, and runs the SAME
    # `process_locus_native`; only the SCHEDULING changes. Pinned by the
    # 24-cell matrix (its on-disk gffutils cells take the fused path) +
    # tests/test_fresh_parallel_step7.py.
    #
    # The fuse applies only when `factory.viable` (DBs on-disk, re-openable
    # per thread). In-memory/blob backends keep every real DB access on the
    # parent thread, then fan fully materialised payloads out to workers.
    # `LIFTON_FUSE_STEP7=0` restores the
    # pre-Iteration-10 two-phase prefetcher-pool path (escape hatch + A/B
    # baseline). `LIFTON_FUSE_MAT_CONCURRENCY=k` (default unset) caps how many
    # workers may be inside the materialise half at once (a Semaphore) — the
    # lever if SQLite contention ever regresses the fused win.
    from lifton.locus_pipeline import (
        materialise_locus, process_locus_native,
        _ThreadLocalCtxFactory, materialise_locus_with_factory,
    )
    factory = _ThreadLocalCtxFactory(ctx)
    open_root_scan = getattr(factory, "open_root_scan", None)
    submission = None
    _fuse_enabled = os.environ.get("LIFTON_FUSE_STEP7", "1") != "0"
    fused = _fuse_enabled and factory.viable
    inflight_limit = _step7_inflight_limit(ctx.args, int(threads))
    if fused:
        strategy = "parallel_fused_worker_materialise"
    elif factory.viable:
        strategy = "parallel_prefetch_worker_materialise"
    else:
        strategy = "parallel_parent_materialise_worker_process"
    try:
        setattr(ctx.args, "_step7_parallel_strategy", strategy)
        if not factory.viable:
            setattr(
                ctx.args,
                "_step7_materialisation_constraint",
                getattr(
                    factory, "fallback_reason",
                    "materialisation databases are not reopenable",
                ),
            )
    except Exception:
        pass

    _mat_sema = None
    if fused:
        _k = os.environ.get("LIFTON_FUSE_MAT_CONCURRENCY")
        if _k:
            try:
                _mat_sema = threading.Semaphore(max(1, int(_k)))
            except ValueError:
                _mat_sema = None

    if fused:
        def _mat_and_process(idx, locus):
            # Materialise + process on ONE worker thread. Wrap the
            # materialise half so a failure becomes a per-locus error
            # LocusResult (mirrors process_locus_native) instead of
            # crashing the pool; Exception (not BaseException) so
            # KeyboardInterrupt/SystemExit still propagate.
            worker_ctx = make_worker_context(ctx, idx, state_coordinator)
            try:
                if _mat_sema is not None:
                    with _mat_sema:
                        payload = materialise_locus_with_factory(
                            idx, locus, factory)
                else:
                    payload = materialise_locus_with_factory(
                        idx, locus, factory)
            except Exception as exc:
                import traceback as _tb
                return LocusResult(
                    index=idx,
                    locus_id=getattr(locus, "id", "<unknown>"),
                    lifton_gene=None,
                    error=exc,
                    error_tb=_tb.format_exc(),
                    delta=worker_ctx.state_journal.finish(),
                )
            except BaseException:
                # Preserve KeyboardInterrupt/SystemExit propagation without
                # leaving later workers blocked at the ordered copy gate.
                worker_ctx.state_journal.finish()
                raise
            return process_locus_native(payload, worker_ctx)

        task_items = None

        def _submit_one(executor, item):
            idx, (_feature, locus) = item
            return executor.submit(_mat_and_process, idx, locus)
    else:
        # Non-fused: materialisation remains a distinct stage, but it now
        # feeds the processing executor through a bounded iterator instead of
        # retaining one full payload and one future per locus.  The on-disk
        # path still uses a separate prefetch pool (the two-pool A/B contract);
        # in-memory/blob backends materialise lazily on the parent thread and
        # run only fully detached payloads on the processing workers.
        if factory.viable:
            # Phase 17c parallel prefetcher pool. Cap at 4 prefetchers
            # since the marginal gain plateaus (~50-60 % reduction at
            # N=4 per the Phase 17 exploration sizing); larger pool
            # increases SQLite contention without commensurate speedup.
            prefetch_workers = min(4, int(threads))

            def _payload_items():
                def _submit_prefetch(executor, item):
                    idx, (_feature, locus) = item
                    return executor.submit(
                        materialise_locus_with_factory,
                        idx, locus, factory,
                    )

                with ThreadPoolExecutor(
                        max_workers=prefetch_workers,
                        thread_name_prefix="lifton-prefetch") as pool:
                    yield from _iter_bounded_ordered(
                        pool, submission, _submit_prefetch,
                        inflight_limit,
                    )
        else:
            # In-memory backends or non-extractable dbfn: the owning parent
            # thread is the only code allowed to touch a real DB connection.
            # Materialisation remains lazy; the bounded map has one vacant
            # slot while next(source) builds a payload, so materialising +
            # submitted + completed-but-not-emitted loci never exceed the
            # configured dispatch window.
            def _payload_items():
                for idx, (_feature, locus) in submission:
                    try:
                        yield materialise_locus(idx, locus, ctx)
                    except Exception as exc:
                        # Match fused-path failure isolation. Resolving the
                        # journal here also prevents later workers from
                        # blocking forever at the ordered state gate.
                        import traceback as _tb
                        worker_ctx = make_worker_context(
                            ctx, idx, state_coordinator,
                        )
                        yield LocusResult(
                            index=idx,
                            locus_id=getattr(locus, "id", "<unknown>"),
                            lifton_gene=None,
                            error=exc,
                            error_tb=_tb.format_exc(),
                            delta=worker_ctx.state_journal.finish(),
                        )

        task_items = _payload_items()

        def _process_materialised(payload):
            if isinstance(payload, LocusResult):
                return payload
            worker_ctx = make_worker_context(
                ctx, payload.submission_index, state_coordinator,
            )
            return process_locus_native(payload, worker_ctx)

        def _submit_one(executor, payload):
            return executor.submit(
                _process_materialised, payload,
            )

    # Every submitted-but-not-emitted locus counts against the window.  This
    # is stronger than limiting executor futures alone: if locus 0 is slow,
    # completed loci 1..N cannot trigger unlimited replenishment.
    raw_max_pending = getattr(ctx.args, "writer_max_pending", 0)
    if not isinstance(raw_max_pending, (int, str)) or \
            isinstance(raw_max_pending, bool):
        raw_max_pending = 0
    try:
        max_pending = max(0, int(raw_max_pending))
    except (TypeError, ValueError):
        max_pending = 0

    spill_tmp = (tempfile.TemporaryDirectory(
        prefix="lifton_writer_spill_") if max_pending > 0 else None)
    writer = _OrderedWriter(
        spill_dir=(spill_tmp.name if spill_tmp is not None else ""),
        max_pending=(max_pending if max_pending > 0 else None),
        consume_fn=_consume_ordered,
        fw=fw,
        stats=transcripts_stats_dict,
        progress_every=progress_every,
        progress_stream=progress_stream,
    )
    root_scan_db = None
    task_source = None
    observed_high_water = 0

    def _record_high_water(value):
        nonlocal observed_high_water
        observed_high_water = max(observed_high_water, value)

    try:
        if factory.viable:
            root_scan_db = (
                open_root_scan()
                if open_root_scan is not None else l_feature_db
            )
            root_items = _iter_loci(features, root_scan_db)
        else:
            # Hierarchy reads and the root scan share the sole parent-owned
            # connection. DuckDB's stream-safe keyset scan prevents child
            # queries from invalidating its active cursor.
            root_scan_db = l_feature_db
            root_items = _iter_loci_serial_safe(features, root_scan_db)
        submission = enumerate(root_items)
        task_source = iter(submission if fused else task_items)
        with ThreadPoolExecutor(max_workers=int(threads)) as executor:
            for result in _ordered_bounded_map(
                    executor, task_source, _submit_one, inflight_limit,
                    high_water_callback=_record_high_water):
                writer.offer(result)
        writer.drain()
        processed = writer.next_to_emit
    finally:
        close_source = getattr(task_source, "close", None)
        if close_source is not None:
            close_source()
        close_factory = getattr(factory, "close", None)
        if close_factory is not None:
            close_factory()
        if root_scan_db is not None and root_scan_db is not l_feature_db:
            close_scan = getattr(factory, "close_root_scan", None)
            if close_scan is not None:
                close_scan(root_scan_db)
        if spill_tmp is not None:
            spill_tmp.cleanup()
        try:
            setattr(
                ctx.args, "_step7_max_inflight_observed",
                observed_high_water,
            )
        except Exception:
            pass

    return processed
