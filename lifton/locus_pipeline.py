# ---------------------------------------------------------------------------
# Phase 9 — locus-major task wrapper.
#
# Provides a single function `process_locus(submission_index, locus, ctx)`
# that calls into the existing `run_liftoff.process_liftoff` exactly once
# for a single Liftoff gene locus and packages the result + the
# submission index into a `LocusResult`. The submission index is what
# the ordered-writer in `lifton.parallel` uses to keep output
# deterministic regardless of worker completion order.
#
# Crucially this module does NOT call `write_entry` on the resulting
# Lifton_GENE. The parent process owns the file handle + stats dict
# and emits in submission order. That separation is what guarantees
# byte-identical output across `--threads ∈ {1,2,4,...}`.
# ---------------------------------------------------------------------------

from __future__ import annotations

import copy as _copy
import io as _io
import os as _os
import threading as _threading
import traceback as _traceback
from collections import OrderedDict as _OrderedDict
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import Any, Dict, List, Optional, Tuple


def safe_exception_text(error, traceback_text=None, *, prefer_traceback=False):
    """Render an exception without trusting its ``__str__`` method.

    Third-party exceptions occasionally implement a broken ``__str__`` that
    raises while LiftOn is already handling the original failure.  A captured
    traceback is a reliable fallback; ``repr`` and finally the type name keep
    reporting non-fatal even for deliberately hostile exception objects.
    """
    if (prefer_traceback and isinstance(traceback_text, str)
            and traceback_text):
        return traceback_text
    try:
        return str(error)
    except Exception:
        if isinstance(traceback_text, str) and traceback_text:
            return traceback_text
        try:
            return repr(error)
        except Exception:
            return f"<unprintable {type(error).__name__}>"


@dataclass
class StepContext:
    """Bundle of every immutable-ish input that ``process_liftoff``
    needs. Built once in the parent process, shared by reference to
    every worker thread (only the alignment helpers mutate the
    objects they receive, and they do so via private members the
    test suite already covers)."""

    ref_db: Any                     # gffutils.FeatureDB or gffbase.FeatureDB
    l_feature_db: Any               # liftoff DB
    m_feature_db: Optional[Any]     # miniprot DB or None
    ref_id_2_m_id_trans_dict: dict
    tree_dict: dict
    tgt_fai: Any                    # pyfaidx.Fasta of the target genome
    ref_proteins: dict
    ref_trans: dict
    ref_features_dict: dict
    fw_score: Any                   # text file handle (writes are short)
    fw_chain: Optional[Any]         # text file handle or None
    args: Any                       # argparse.Namespace (carries --legacy-merge etc.)
    # Set only on a per-task context constructed by ``make_worker_context``.
    # The public/default context keeps this as None, preserving direct callers.
    state_journal: Optional[Any] = None
    # Ordered parent-thread ledger. ``parallel_step7`` keeps its integer return
    # API while callers can inspect this list for manifest/failure policy data.
    failure_records: List[dict] = field(default_factory=list)
    # Canonical reference genes whose GFF3 hierarchy was actually emitted.
    # Only ordered parent-thread ``consume`` mutates this set.
    emitted_ref_gene_ids: set = field(default_factory=set)
    emitted_intervals: List[Tuple[str, int, int, str]] = field(
        default_factory=list,
    )


@dataclass
class LocusDelta:
    """Buffered side effects produced while processing one locus.

    Worker threads may build arbitrarily complex ``Lifton_GENE`` objects, but
    they must not mutate Step 7's shared counters/interval tree or write the
    shared score/chain streams.  Those effects are captured here and committed
    by :func:`consume` in submission order.
    """

    copy_increments: Dict[str, int] = field(default_factory=dict)
    tree_intervals: List[Tuple[str, int, int, str]] = field(default_factory=list)
    score_text: str = ""
    chain_text: str = ""
    committed: bool = False


def commit_locus_delta(delta: LocusDelta, ref_features_dict: dict,
                       tree_dict: dict) -> None:
    """Commit a construction delta exactly once.

    This helper is intentionally independent of the Step 7 coordinator so
    serial evaluate-then-commit paths (notably miniprot rescue) can construct a
    gene with :class:`DeferredStateJournal`, serialize it privately, and only
    publish its copy/tree effects after successful serialization.
    """
    if delta.committed:
        return
    for ref_gene_id, increment in delta.copy_increments.items():
        feature = ref_features_dict.get(ref_gene_id)
        if feature is not None:
            feature.copy_num += increment
    if delta.tree_intervals:
        from intervaltree import IntervalTree
        from lifton.intervals import _make_interval
        for seqid, start, end, gene_id in delta.tree_intervals:
            if seqid not in tree_dict:
                tree_dict[seqid] = IntervalTree()
            tree_dict[seqid].add(_make_interval(start, end, gene_id))
    delta.committed = True


class DeferredStateJournal:
    """Journal for a serial candidate whose state commits only on success.

    Unlike the Step 7 journal, this class has no ordered allocation gate.  It
    reads the current committed copy counter, making it suitable for serial
    candidate evaluation: discarding the delta leaves the next candidate free
    to reuse the unconsumed copy number.
    """

    def __init__(self, ref_features_dict: dict, *, buffer_score: bool = False,
                 buffer_chain: bool = False):
        self.ref_features_dict = ref_features_dict
        self.delta = LocusDelta()
        self.score_handle = _io.StringIO() if buffer_score else None
        self.chain_handle = _io.StringIO() if buffer_chain else None
        self._allocated = False
        self._finished = False

    def allocate_gene_copy(self, ref_gene_id, entry, _ref_features_dict):
        if self._allocated:
            raise RuntimeError("candidate attempted more than one root gene allocation")
        attrs = getattr(entry, "attributes", {}) or {}
        explicit = (int(attrs["extra_copy_number"][0])
                    if "extra_copy_number" in attrs else None)
        if ref_gene_id in self.ref_features_dict:
            copy_num = (explicit if explicit is not None else int(getattr(
                self.ref_features_dict[ref_gene_id], "copy_num", 0,
            )))
            self.delta.copy_increments[ref_gene_id] = 1
        else:
            copy_num = explicit if explicit is not None else 0
        gene_id = (f"{ref_gene_id}_{copy_num}"
                   if copy_num > 0 else ref_gene_id)
        if getattr(entry, "featuretype", None) == "gene":
            self.delta.tree_intervals.append((
                entry.seqid, int(entry.start), int(entry.end), gene_id,
            ))
        self._allocated = True
        return copy_num

    def finish(self) -> LocusDelta:
        if self._finished:
            return self.delta
        if self.score_handle is not None:
            self.delta.score_text = self.score_handle.getvalue()
        if self.chain_handle is not None:
            self.delta.chain_text = self.chain_handle.getvalue()
        self._finished = True
        return self.delta


class _JournalIntervalTree:
    """Read-only interval-tree view used by one worker.

    Reads are serialized with parent-thread commits and augmented with the
    deterministic logical deltas registered for this locus and its predecessors.
    This reproduces the serial visibility rules without exposing a mutating
    ``IntervalTree`` to worker threads.
    """

    __slots__ = ("_coordinator", "_seqid", "_index")

    def __init__(self, coordinator, seqid: str, index: int):
        self._coordinator = coordinator
        self._seqid = seqid
        self._index = index

    def overlap(self, begin: int, end: int):
        return self._coordinator.overlap(
            self._seqid, begin, end, through_index=self._index,
        )


class _JournalTreeDict(Mapping):
    """Mapping facade whose values expose the read-only ``overlap`` API."""

    __slots__ = ("_coordinator", "_index")

    def __init__(self, coordinator, index: int):
        self._coordinator = coordinator
        self._index = index

    def __getitem__(self, seqid: str):
        if not self._coordinator.has_seqid(seqid, self._index):
            raise KeyError(seqid)
        return _JournalIntervalTree(self._coordinator, seqid, self._index)

    def __iter__(self):
        return iter(self._coordinator.seqids(self._index))

    def __len__(self):
        return len(self._coordinator.seqids(self._index))


class Step7StateCoordinator:
    """Allocate copy IDs deterministically and commit worker deltas safely.

    Copy allocation is a very short ordered gate near the start of each locus;
    expensive sequence/alignment work remains concurrent.  Global counters,
    interval trees, and real output handles are mutated only by ordered
    :func:`consume` calls on the parent thread.
    """

    def __init__(self, ref_features_dict: dict, tree_dict: dict):
        self.ref_features_dict = ref_features_dict
        self.tree_dict = tree_dict
        self._condition = _threading.Condition(_threading.RLock())
        self._next_index = 0
        self._resolved_indices = set()
        self._virtual_copy_counts: Dict[str, int] = {}
        self._logical_intervals: Dict[str, List[Tuple[int, Any]]] = {}

    def allocate_gene(self, index: int, ref_gene_id: Optional[str], entry,
                      delta: LocusDelta) -> int:
        """Return the copy number this locus receives in serial order."""
        with self._condition:
            while index != self._next_index:
                self._condition.wait()

            explicit = None
            attrs = getattr(entry, "attributes", {}) or {}
            if "extra_copy_number" in attrs:
                explicit = int(attrs["extra_copy_number"][0])

            if ref_gene_id in self.ref_features_dict:
                if ref_gene_id not in self._virtual_copy_counts:
                    self._virtual_copy_counts[ref_gene_id] = int(
                        getattr(self.ref_features_dict[ref_gene_id],
                                "copy_num", 0)
                    )
                copy_num = (explicit if explicit is not None else
                            self._virtual_copy_counts[ref_gene_id])
                self._virtual_copy_counts[ref_gene_id] += 1
                delta.copy_increments[ref_gene_id] = (
                    delta.copy_increments.get(ref_gene_id, 0) + 1
                )
            else:
                copy_num = explicit if explicit is not None else 0

            gene_id = (f"{ref_gene_id}_{copy_num}"
                       if copy_num > 0 else ref_gene_id)
            if getattr(entry, "featuretype", None) == "gene":
                seqid = entry.seqid
                interval_tuple = (seqid, int(entry.start), int(entry.end),
                                  gene_id)
                delta.tree_intervals.append(interval_tuple)
                from lifton.intervals import _make_interval
                interval = _make_interval(entry.start, entry.end, gene_id)
                self._logical_intervals.setdefault(seqid, []).append(
                    (index, interval)
                )

            self._resolved_indices.add(index)
            self._next_index += 1
            self._condition.notify_all()
            return copy_num

    def finish_without_allocation(self, index: int) -> None:
        """Release the ordered gate when a locus never built a gene."""
        with self._condition:
            if index in self._resolved_indices:
                return
            while index != self._next_index:
                self._condition.wait()
            self._resolved_indices.add(index)
            self._next_index += 1
            self._condition.notify_all()

    def seqids(self, through_index: int):
        with self._condition:
            result = set(self.tree_dict.keys())
            result.update(
                seqid for seqid, rows in self._logical_intervals.items()
                if any(index <= through_index for index, _ in rows)
            )
            return tuple(result)

    def has_seqid(self, seqid: str, through_index: int) -> bool:
        with self._condition:
            if seqid in self.tree_dict:
                return True
            return any(
                index <= through_index
                for index, _ in self._logical_intervals.get(seqid, ())
            )

    def overlap(self, seqid: str, begin: int, end: int,
                *, through_index: int):
        with self._condition:
            base_tree = self.tree_dict.get(seqid)
            overlaps = (set(base_tree.overlap(begin, end))
                        if base_tree is not None else set())
            for index, interval in self._logical_intervals.get(seqid, ()):
                if index <= through_index and interval.overlaps(begin, end):
                    overlaps.add(interval)
            return overlaps

    def commit(self, delta: LocusDelta) -> None:
        """Apply one delta exactly once; caller guarantees index order."""
        with self._condition:
            commit_locus_delta(
                delta, self.ref_features_dict, self.tree_dict,
            )


class LocusStateJournal:
    """Per-worker journal backing a :class:`LocusDelta`."""

    def __init__(self, index: int, coordinator: Step7StateCoordinator,
                 *, buffer_score: bool, buffer_chain: bool):
        self.index = index
        self.coordinator = coordinator
        self.delta = LocusDelta()
        self.score_handle = _io.StringIO() if buffer_score else None
        self.chain_handle = _io.StringIO() if buffer_chain else None
        self._allocated = False
        self._finished = False

    @property
    def tree_view(self):
        return _JournalTreeDict(self.coordinator, self.index)

    def allocate_gene_copy(self, ref_gene_id, entry, _ref_features_dict):
        if self._allocated:
            raise RuntimeError(
                f"locus {self.index} attempted more than one root gene allocation"
            )
        copy_num = self.coordinator.allocate_gene(
            self.index, ref_gene_id, entry, self.delta,
        )
        self._allocated = True
        return copy_num

    def finish(self) -> LocusDelta:
        if self._finished:
            return self.delta
        if not self._allocated:
            self.coordinator.finish_without_allocation(self.index)
        if self.score_handle is not None:
            self.delta.score_text = self.score_handle.getvalue()
        if self.chain_handle is not None:
            self.delta.chain_text = self.chain_handle.getvalue()
        self._finished = True
        return self.delta


def make_worker_context(ctx: StepContext, index: int,
                        coordinator: Step7StateCoordinator) -> StepContext:
    """Return a task-local context with journaled state and sidecars."""
    journal = LocusStateJournal(
        index, coordinator,
        buffer_score=ctx.fw_score is not None,
        buffer_chain=ctx.fw_chain is not None,
    )
    return replace(
        ctx,
        tree_dict=journal.tree_view,
        fw_score=journal.score_handle,
        fw_chain=journal.chain_handle,
        state_journal=journal,
    )


@dataclass
class LocusResult:
    """Result of one per-locus task. ``error`` is populated when
    `process_liftoff` raises so the parent process can log the
    failure and keep the rest of the pipeline running."""

    index: int                       # submission index (0-based)
    locus_id: str                    # gffutils-feature id, for logging
    lifton_gene: Optional[Any] = None
    error: Optional[BaseException] = None
    error_tb: Optional[str] = None   # formatted traceback captured at the raise
                                     # site, so the parent can log the real cause
                                     # even when the exception's __str__ is broken
    delta: LocusDelta = field(default_factory=LocusDelta)

    @property
    def emittable(self) -> bool:
        """True iff the result has a valid Lifton_GENE that should be
        passed to ``write_entry``."""
        return (
            self.lifton_gene is not None
            and getattr(self.lifton_gene, "ref_gene_id", None) is not None
            and self.error is None
        )


def process_locus(submission_index: int, locus, *, ctx: StepContext) -> LocusResult:
    """Run the Step 7 per-locus body for a single gene.

    This is the call that workers in `lifton.parallel` schedule. It is
    intentionally thin: it dispatches to the existing
    `run_liftoff.process_liftoff` so that the byte-output is byte-equal
    to the legacy serial loop.

    Exceptions are caught and packaged in :class:`LocusResult` so a
    single bad locus cannot bring down sibling workers.
    """
    # Local import keeps the module-load graph free of cycles and
    # matches the lazy-import pattern used elsewhere in lifton/.
    from lifton import run_liftoff
    journal = getattr(ctx, "state_journal", None)

    # The merge mode (default best-of-outcome vs --legacy-merge) is read from
    # ctx.args inside process_liftoff, so nothing extra is threaded here.
    try:
        gene = run_liftoff.process_liftoff(
            None, locus,
            ctx.ref_db,
            ctx.l_feature_db,
            ctx.ref_id_2_m_id_trans_dict,
            ctx.m_feature_db,
            ctx.tree_dict,
            ctx.tgt_fai,
            ctx.ref_proteins,
            ctx.ref_trans,
            ctx.ref_features_dict,
            ctx.fw_score,
            ctx.fw_chain,
            ctx.args,
            ENTRY_FEATURE=True,
            state_journal=journal,
        )
    except Exception as exc:        # V1.4 fix: narrowed from BaseException
        # KeyboardInterrupt, SystemExit, GeneratorExit are BaseException
        # but NOT Exception, so they propagate and let Ctrl-C kill the
        # whole pool instead of being silently packaged into LocusResult.
        return LocusResult(
            index=submission_index,
            locus_id=getattr(locus, "id", "<unknown>"),
            lifton_gene=None,
            error=exc,
            error_tb=_traceback.format_exc(),
            delta=(journal.finish() if journal is not None else LocusDelta()),
        )
    finally:
        # This also releases the coordinator's ordered allocation gate when
        # the locus returned before constructing a root gene. BaseException
        # still propagates, but cannot strand sibling workers behind the gate.
        if journal is not None:
            journal.finish()
    return LocusResult(
        index=submission_index,
        locus_id=getattr(locus, "id", "<unknown>"),
        lifton_gene=gene,
        delta=(journal.delta if journal is not None else LocusDelta()),
    )


# ---------------------------------------------------------------------------
# Phase 11 — pre-materialised per-locus payload + native task wrapper.
# ---------------------------------------------------------------------------

@dataclass
class _FeaturePreFetch:
    """Pre-fetched ``l_feature_db`` data for a single feature visited
    during the recursive descent of ``process_liftoff``. Populated by
    :func:`_walk_and_cache_features` on the parent thread; read by
    :class:`_LFeatureDbProxy` on workers.

    Each list captures the exact shape of children the per-locus body
    queries via ``l_feature_db.children(feature, ...)``:

    - ``children_l1``         — ``children(level=1)``
    - ``exon_children_l1``    — ``children(featuretype='exon', level=1, order_by='start')``
    - ``exon_children_full``  — ``children(featuretype='exon', order_by='start')`` (no level)
    - ``cds_stop_children``   — ``children(featuretype=('CDS','stop_codon'), order_by='start')``
    """
    feature: Any
    children_l1: list = field(default_factory=list)
    exon_children_l1: list = field(default_factory=list)
    exon_children_full: list = field(default_factory=list)
    cds_stop_children: list = field(default_factory=list)


@dataclass
class _MiniprotPreFetch:
    """Pre-fetched ``m_feature_db`` data for one miniprot ``m_id``.

    The per-locus body at ``lifton_utils.LiftOn_miniprot_alignment``
    accesses ``m_feature_db[m_id]`` (returns the mRNA Feature) and
    ``m_feature_db.children(m_entry, featuretype=('CDS','stop_codon'), order_by='start')``.
    Both are pre-fetched here so workers never touch the shared
    miniprot DB connection.
    """
    feature: Any                            # the mRNA Feature
    cds_stop_children: list = field(default_factory=list)


class _RefDbProxy:
    """Read-only proxy that satisfies the ``ref_db[id]`` access pattern
    used by :mod:`lifton.run_liftoff` from a pre-materialised
    ``{id -> attributes}`` cache.

    The legacy ``run_liftoff`` reads only ``ref_db[id].attributes``
    (lines 112, 131) plus a ``ref_db[id]`` existence check (line 252)
    that's caught with ``(KeyError, FeatureNotFoundError)``. A minimal
    proxy that returns an attribute-bearing stub on hit and raises
    ``KeyError`` on miss is sufficient.
    """
    __slots__ = ("_attrs_cache",)

    def __init__(self, attrs_cache: "Dict[str, dict]"):
        self._attrs_cache = attrs_cache

    def __getitem__(self, feature_id: str) -> "_RefFeatureStub":
        if feature_id not in self._attrs_cache:
            raise KeyError(feature_id)
        return _RefFeatureStub(feature_id, self._attrs_cache[feature_id])


class _RefFeatureStub:
    """Minimal Feature-shaped object exposing ``.id`` and
    ``.attributes`` — the only attributes the per-locus body reads
    from ``ref_db[id]``."""
    __slots__ = ("id", "attributes")

    def __init__(self, feature_id: str, attributes: dict):
        self.id = feature_id
        self.attributes = attributes


class _LFeatureDbProxy:
    """Read-only proxy for ``ctx.l_feature_db``. Satisfies the four
    ``children(...)`` signatures used by ``run_liftoff.process_liftoff``
    (and ``lifton_add_trans_exon_cds``) plus the ``[id]`` bracket
    lookup, all from a ``{feature_id -> _FeaturePreFetch}`` cache.

    Because the legacy code uses ``order_by='start'`` consistently and
    the cache is built with the same ordering, output is byte-equal.
    """
    __slots__ = ("_cache",)

    def __init__(self, cache: "Dict[str, _FeaturePreFetch]"):
        self._cache = cache

    def __getitem__(self, feature_id: str):
        if feature_id not in self._cache:
            raise KeyError(feature_id)
        return self._cache[feature_id].feature

    def children(self, feature, featuretype=None, level=None,
                 order_by=None):
        # Normalise the lookup key. The cache stores by feature.id;
        # the call site passes either a Feature object or a raw id.
        feature_id = getattr(feature, "id", feature)
        entry = self._cache.get(feature_id)
        if entry is None:
            return iter([])
        # Dispatch on (featuretype, level). The cache pre-fetches the
        # exact four signatures the per-locus body uses.
        if featuretype == "exon" and level == 1:
            return iter(entry.exon_children_l1)
        if featuretype == "exon" and level is None:
            return iter(entry.exon_children_full)
        if featuretype == ("CDS", "stop_codon"):
            return iter(entry.cds_stop_children)
        if featuretype is None and level == 1:
            return iter(entry.children_l1)
        # Any other call signature is a programmer error: the cache
        # was not pre-fetched for this shape. Fail loudly so a future
        # refactor that adds a new children() variant in the per-locus
        # body is caught immediately.
        raise NotImplementedError(
            f"_LFeatureDbProxy.children: un-cached signature "
            f"(featuretype={featuretype!r}, level={level!r}); "
            f"add the call signature to _walk_and_cache_features."
        )


class _MFeatureDbProxy:
    """Read-only proxy for ``ctx.m_feature_db``. Workers need
    ``m_feature_db[m_id]`` (returns the mRNA Feature) and
    ``m_feature_db.children(m_entry, featuretype=('CDS','stop_codon'),
    order_by='start')`` — both served from the
    ``{m_id -> _MiniprotPreFetch}`` cache.
    """
    __slots__ = ("_cache",)

    def __init__(self, cache: "Dict[str, _MiniprotPreFetch]"):
        self._cache = cache

    def __getitem__(self, m_id: str):
        if m_id not in self._cache:
            raise KeyError(m_id)
        return self._cache[m_id].feature

    def children(self, feature, featuretype=None, level=None,
                 order_by=None):
        m_id = getattr(feature, "id", feature)
        entry = self._cache.get(m_id)
        if entry is None:
            return iter([])
        if featuretype == ("CDS", "stop_codon"):
            return iter(entry.cds_stop_children)
        raise NotImplementedError(
            f"_MFeatureDbProxy.children: un-cached signature "
            f"(featuretype={featuretype!r}); the only signature the "
            f"per-locus body uses is featuretype=('CDS','stop_codon')."
        )


@dataclass
class MaterialisedLocus:
    """Frozen-by-convention snapshot of every input ``process_liftoff``
    needs for a single Liftoff gene. All fields are populated in the
    parent thread via :func:`materialise_locus` BEFORE any worker
    starts; workers consume the payload read-only via the proxy DBs
    in :func:`process_locus_native`.

    The point: workers never call ``l_feature_db.children(...)`` or
    any other shared-FeatureDB method directly — they go through
    :class:`_LFeatureDbProxy`, :class:`_RefDbProxy`,
    :class:`_MFeatureDbProxy`, all backed by the caches below. With
    workers fully decoupled from the real DBs, ``ThreadPoolExecutor``
    is safe under any backend (gffutils SQLite, gffbase DuckDB).

    Phase 17 (option b — completed Phase 11/12 materialisation):
    extended with ``feature_cache``, ``ref_attrs_cache``,
    ``miniprot_cache`` so the existing legacy flat fields
    (``children_l1`` etc.) become a thin compatibility surface and
    the proxies become the canonical worker access path.
    """
    submission_index: int
    locus: Any                              # gffutils.Feature snapshot
    locus_id: str

    # Pre-fetched DB reads (parent-thread only):
    children_l1: list = field(default_factory=list)
    exon_children: list = field(default_factory=list)
    cds_children: list = field(default_factory=list)
    cds_stop_children: list = field(default_factory=list)

    # Reference-side snapshots:
    ref_gene_id: Optional[str] = None
    ref_trans_id: Optional[str] = None
    ref_gene_attrs: dict = field(default_factory=dict)
    ref_trans_attrs: dict = field(default_factory=dict)

    # Phase 17 (option b) — proxy-backed caches. Workers read the per-
    # locus DB-equivalents from these dicts via the proxy classes
    # above. ``feature_cache`` follows the same branch descent as
    # ``run_liftoff.process_liftoff``: containers recurse through
    # level-1 children, while features with direct exons are terminal;
    # ``ref_attrs_cache`` carries the ref-side ``[id].attributes``
    # lookups; ``miniprot_cache`` carries every ``m_feature_db[m_id]``
    # the chaining algorithm might consult.
    feature_cache: Dict[str, _FeaturePreFetch] = field(default_factory=dict)
    ref_attrs_cache: Dict[str, dict] = field(default_factory=dict)
    miniprot_cache: Dict[str, _MiniprotPreFetch] = field(default_factory=dict)


class HierarchyBatchLoader:
    """Backend-neutral full-feature hierarchy loader.

    gffbase serves each anchor batch with one set-based query; gffutils and
    lightweight embedders retain scalar generator semantics.  The result shape
    is identical in both modes, which keeps materialisation independent of the
    selected annotation backend.
    """

    def __init__(self, db, *, slots: int = 1, batch_size: Optional[int] = None):
        self.db = db
        raw = (_os.environ.get("LIFTON_HIERARCHY_BATCH_SIZE")
               if batch_size is None else batch_size)
        if raw is None:
            resolved = min(128, max(1, int(slots or 1)))
        else:
            try:
                resolved = max(1, int(raw))
            except (TypeError, ValueError):
                resolved = min(128, max(1, int(slots or 1)))
        self.batch_size = resolved
        self.batched = callable(getattr(db, "children_batched_features", None))

    @staticmethod
    def _ids(anchors):
        return [getattr(anchor, "id", anchor) for anchor in anchors]

    @staticmethod
    def _is_missing_feature_error(exc):
        """Recognize lookup misses from supported annotation backends."""
        if isinstance(exc, KeyError):
            return True
        error_type = type(exc)
        module = str(getattr(error_type, "__module__", ""))
        return (
            error_type.__name__ == "FeatureNotFoundError"
            and module.startswith(("gffutils", "gffbase", "lifton.gffbase"))
        )

    def children_many(
        self, anchors, *, featuretype=None, level=None, order_by=None,
    ):
        ids = self._ids(anchors)
        result = _OrderedDict((anchor, ()) for anchor in ids)
        if self.batched:
            for offset in range(0, len(ids), self.batch_size):
                chunk = ids[offset:offset + self.batch_size]
                result.update(self.db.children_batched_features(
                    chunk, featuretype=featuretype, level=level,
                    order_by=order_by,
                ))
            return result
        for anchor in ids:
            result[anchor] = tuple(self.db.children(
                anchor, featuretype=featuretype, level=level,
                order_by=order_by,
            ))
        return result

    def features_many(self, feature_ids):
        ids = list(_OrderedDict.fromkeys(feature_ids))
        batched = getattr(self.db, "features_batched", None)
        if callable(batched):
            result = _OrderedDict()
            for offset in range(0, len(ids), self.batch_size):
                result.update(batched(ids[offset:offset + self.batch_size]))
            return result
        result = _OrderedDict()
        for feature_id in ids:
            try:
                result[feature_id] = self.db[feature_id]
            except Exception as exc:
                if not self._is_missing_feature_error(exc):
                    raise
        return result


def _walk_and_cache_features_batched(
    feature, ctx: StepContext, payload: MaterialisedLocus,
    loader: HierarchyBatchLoader, *, max_depth: int = 8,
) -> None:
    """Snapshot exactly the hierarchy branches ``process_liftoff`` reads.

    A feature with direct level-1 exons is transcript-shaped: the runtime
    consumes its exons and CDS/stop-codon descendants, then stops.  A feature
    without direct exons is a container: the runtime reads its level-1
    children and recurses.  Keeping that split here avoids querying every
    signature for every node and avoids ever placing ordinary exon/CDS leaves
    on the next frontier.
    """
    frontier = [feature]
    depth = 0
    while frontier and depth <= max_depth:
        unique = []
        frontier_ids = set()
        for current in frontier:
            feature_id = getattr(current, "id", None)
            if (feature_id is not None
                    and feature_id not in payload.feature_cache
                    and feature_id not in frontier_ids):
                unique.append(current)
                frontier_ids.add(feature_id)
        if not unique:
            return

        exon_l1 = loader.children_many(
            unique, featuretype="exon", level=1, order_by="start",
        )
        containers = [
            current for current in unique
            if not exon_l1.get(getattr(current, "id", current), ())
        ]
        terminals = [
            current for current in unique
            if exon_l1.get(getattr(current, "id", current), ())
        ]
        children_l1 = (
            loader.children_many(containers, level=1) if containers else {}
        )
        cds_stop = (
            loader.children_many(
                terminals,
                featuretype=("CDS", "stop_codon"),
                order_by="start",
            ) if terminals else {}
        )

        next_frontier = []
        for current in unique:
            feature_id = current.id
            direct_exons = list(exon_l1.get(feature_id, ()))
            direct_children = (
                [] if direct_exons else list(children_l1.get(feature_id, ()))
            )
            entry = _FeaturePreFetch(
                feature=current,
                children_l1=direct_children,
                exon_children_l1=direct_exons,
                # The runtime only asks for the no-level exon form after the
                # direct-exon branch has succeeded. Standard GFF3 exons are
                # leaves, so that result is the same ordered list.
                exon_children_full=list(direct_exons),
                cds_stop_children=list(cds_stop.get(feature_id, ())),
            )
            payload.feature_cache[feature_id] = entry
            next_frontier.extend(direct_children)
        frontier = next_frontier
        depth += 1


def _walk_and_cache_features(feature, ctx: StepContext, payload: MaterialisedLocus,
                             *, depth: int = 0, max_depth: int = 8) -> None:
    """Recursively pre-fetch every ``l_feature_db`` access the legacy
    ``run_liftoff.process_liftoff`` will issue under ``feature``.

    The legacy code (``lifton/run_liftoff.py:226-260``) does:

    1. ``children(feature, featuretype='exon', level=1, order_by='start')``
    2. If empty, ``children(feature, level=1)`` and recurse into each
       returned child.
    3. If non-empty, ``children(feature, featuretype='exon', order_by='start')``
       (no level) and ``children(feature, featuretype=('CDS','stop_codon'),
       order_by='start')`` via ``lifton_add_trans_exon_cds``.

    All four signatures get cached for every visited feature. The
    ``max_depth`` guard mirrors the V5.2 cycle-detection guard at
    ``run_liftoff.py:212-224`` — if the GFF3 hierarchy exceeds 8 levels
    we bail; deeper nesting is pathological.
    """
    from lifton import logger
    feature_id = getattr(feature, "id", None)
    if feature_id is None or feature_id in payload.feature_cache:
        return
    if depth > max_depth:
        logger.log_warning(
            f"_walk_and_cache_features({feature_id}): exceeded max depth "
            f"{max_depth}; deeper hierarchy not pre-fetched (recursion "
            f"will surface as KeyError on the proxy)."
        )
        return

    entry = _FeaturePreFetch(feature=feature)
    # Variant 1: featuretype='exon', level=1. This also CLASSIFIES the feature, exactly
    # as the batched sibling `_walk_and_cache_features_batched` does:
    #   * direct level-1 exons  -> TRANSCRIPT-shaped ("terminal"): the runtime consumes
    #     its exons plus CDS/stop-codon descendants and stops.
    #   * no direct exons       -> CONTAINER: the runtime reads its level-1 children and
    #     recurses into them.
    # The scalar walker used to run all four query variants for EVERY feature and then
    # recurse into `children_l1` -- which for a transcript is its exons and CDS. Each of
    # those leaves then paid ~4 more SQLite round-trips and got a cache entry of its own,
    # for signatures the runtime never asks a leaf about (the proxy answers an un-cached
    # feature with an empty iterator either way). On a human RefSeq lift that is millions
    # of redundant queries. Keeping the container/terminal split here removes them.
    try:
        entry.exon_children_l1 = list(
            ctx.l_feature_db.children(
                feature, featuretype="exon", level=1, order_by="start",
            )
        )
    except Exception as exc:
        logger.log_warning(
            f"_walk_and_cache_features({feature_id}): exon@l1 failed: {exc}"
        )

    if entry.exon_children_l1:
        # Terminal. Variant 3 (no-level exons) equals Variant 1 under the leaf-exon
        # convention (Phase 17c Win 1), and Variant 2 is only ever asked of containers.
        entry.exon_children_full = list(entry.exon_children_l1)
        # Variant 4: featuretype=('CDS','stop_codon') (no level)
        try:
            entry.cds_stop_children = list(
                ctx.l_feature_db.children(
                    feature, featuretype=("CDS", "stop_codon"), order_by="start",
                )
            )
        except Exception as exc:
            logger.log_warning(
                f"_walk_and_cache_features({feature_id}): CDS+stop_codon failed: {exc}"
            )
        payload.feature_cache[feature_id] = entry
        return

    # Container. Variant 2: level=1 (any type)
    try:
        entry.children_l1 = list(
            ctx.l_feature_db.children(feature, level=1)
        )
    except Exception as exc:
        logger.log_warning(
            f"_walk_and_cache_features({feature_id}): children@l1 failed: {exc}"
        )
    payload.feature_cache[feature_id] = entry

    # Recurse into every level-1 child so the proxy can satisfy the
    # process_liftoff recursion at run_liftoff.py:239.
    for child in entry.children_l1:
        _walk_and_cache_features(child, ctx, payload,
                                 depth=depth + 1, max_depth=max_depth)


def _maybe_cache_ref_attrs(payload: MaterialisedLocus, ref_db,
                           ref_id: Optional[str]) -> None:
    """Pre-fetch ``ref_db[ref_id].attributes`` into the payload cache.
    No-op if ``ref_id`` is None or already cached. Errors are
    swallowed and tracked via the cache miss (KeyError raised by the
    proxy on access)."""
    if ref_id is None or ref_id in payload.ref_attrs_cache:
        return
    try:
        attrs = _copy.deepcopy(ref_db[ref_id].attributes)
        payload.ref_attrs_cache[ref_id] = attrs
    except Exception:
        # Leave it absent → proxy will raise KeyError, which the
        # legacy code at run_liftoff.py:253 catches.
        pass


def _populate_ref_attrs_for_descent(payload: MaterialisedLocus,
                                    ctx: StepContext) -> None:
    """For every feature in ``payload.feature_cache`` that the
    per-locus body would treat as a transcript (i.e. has level-1 exon
    children), look up its ``ref_gene_id`` + ``ref_trans_id`` via
    ``lifton_utils.get_ref_ids_liftoff`` and pre-fetch the
    corresponding ``ref_db[id].attributes``.
    """
    from lifton import lifton_utils as _lu
    # The top-level locus is the gene; lookup is keyed by (gene_id, None).
    locus_id = payload.locus_id
    try:
        ref_gene_id, ref_trans_id_for_gene = _lu.get_ref_ids_liftoff(
            ctx.ref_features_dict, locus_id, None,
        )
    except Exception:
        ref_gene_id, ref_trans_id_for_gene = None, None
    payload.ref_gene_id = ref_gene_id
    payload.ref_trans_id = ref_trans_id_for_gene
    ref_ids = [
        ref_id for ref_id in (ref_gene_id, ref_trans_id_for_gene)
        if ref_id is not None
    ]

    # For each cached transcript-shaped feature (one with level-1 exon
    # children) under the locus, look up the (ref_gene_id, ref_trans_id)
    # pair the recursive call site at run_liftoff.py:244 will compute.
    for fid, entry in payload.feature_cache.items():
        if not entry.exon_children_l1:
            continue
        if fid == locus_id:
            # Already handled above (the locus IS the transcript path).
            continue
        try:
            r_gene_id, r_trans_id = _lu.get_ref_ids_liftoff(
                ctx.ref_features_dict, locus_id, fid,
            )
        except Exception:
            continue
        for ref_id in (r_gene_id, r_trans_id):
            if ref_id is not None and ref_id not in ref_ids:
                ref_ids.append(ref_id)

    loader = HierarchyBatchLoader(
        ctx.ref_db, slots=getattr(ctx.args, "threads", 1),
    )
    for ref_id, feature in loader.features_many(ref_ids).items():
        payload.ref_attrs_cache[ref_id] = _copy.deepcopy(feature.attributes)


def _populate_miniprot_cache(payload: MaterialisedLocus,
                             ctx: StepContext) -> None:
    """Pre-fetch every ``m_feature_db[m_id]`` + its CDS/stop_codon
    children that ``LiftOn_miniprot_alignment`` might consult for
    the ref_trans_ids the per-locus body will encounter.
    """
    if ctx.m_feature_db is None:
        return
    from lifton import logger
    # Collect every ref_trans_id we cached (the only ones the legacy
    # code looks up via ref_id_2_m_id_trans_dict).
    m_ids = []
    for r_trans_id in payload.ref_attrs_cache:
        for m_id in ctx.ref_id_2_m_id_trans_dict.get(r_trans_id) or []:
            if m_id not in m_ids:
                m_ids.append(m_id)
    loader = HierarchyBatchLoader(
        ctx.m_feature_db, slots=getattr(ctx.args, "threads", 1),
    )
    entries = loader.features_many(m_ids)
    try:
        descendants = loader.children_many(
            list(entries.values()),
            featuretype=("CDS", "stop_codon"), order_by="start",
        )
    except Exception as exc:
        logger.log_warning(
            f"_populate_miniprot_cache: batched children failed: {exc}"
        )
        descendants = {}
    for m_id in m_ids:
        m_entry = entries.get(m_id)
        if m_entry is None:
            logger.log_warning(
                f"_populate_miniprot_cache({m_id}): lookup failed"
            )
            continue
        payload.miniprot_cache[m_id] = _MiniprotPreFetch(
            feature=m_entry,
            cds_stop_children=list(descendants.get(m_id, ())),
        )


def materialise_locus(submission_index: int, locus,
                      ctx: StepContext) -> MaterialisedLocus:
    """Parent-thread-only function. Pre-fetches **everything**
    ``run_liftoff.process_liftoff`` would have read from
    ``ctx.l_feature_db``, ``ctx.ref_db``, and ``ctx.m_feature_db`` so
    the worker thread is purely CPU + payload reads.

    Phase 17 (option b): walks the locus's runtime-reachable descent and
    populates ``payload.feature_cache`` (l_feature_db),
    ``payload.ref_attrs_cache`` (ref_db), and
    ``payload.miniprot_cache`` (m_feature_db). The legacy flat fields
    (``children_l1``, ``exon_children``, ``ref_gene_attrs``, etc.)
    are still populated to keep existing test invariants stable.
    """
    from lifton import logger

    locus_id = getattr(locus, "id", "<unknown>")
    payload = MaterialisedLocus(
        submission_index=submission_index,
        locus=locus,
        locus_id=locus_id,
    )

    hierarchy_loader = HierarchyBatchLoader(
        ctx.l_feature_db,
        slots=getattr(ctx.args, "threads", 1),
    )

    # ── Phase 17b: cache the runtime-reachable hierarchy under locus.
    try:
        if hierarchy_loader.batched:
            _walk_and_cache_features_batched(
                locus, ctx, payload, hierarchy_loader,
            )
        else:
            _walk_and_cache_features(locus, ctx, payload)
    except Exception as exc:
        logger.log_warning(
            f"materialise_locus({locus_id}): feature-tree walk failed: {exc}"
        )

    # Populate the legacy flat fields from the cache for back-compat.
    locus_entry = payload.feature_cache.get(locus_id)
    if locus_entry is not None:
        payload.children_l1 = list(locus_entry.children_l1)
        payload.exon_children = list(locus_entry.exon_children_l1)
        payload.cds_stop_children = list(locus_entry.cds_stop_children)
        if not payload.cds_stop_children:
            # A container locus no longer performs a redundant all-descendant
            # CDS query. Reconstruct the legacy flat view from terminal cache
            # entries; workers use ``feature_cache`` directly.
            seen_descendants = set()
            for entry in payload.feature_cache.values():
                if not entry.exon_children_l1:
                    continue
                for descendant in entry.cds_stop_children:
                    descendant_id = getattr(descendant, "id", None)
                    identity = (
                        ("id", descendant_id)
                        if descendant_id is not None
                        else ("object", id(descendant))
                    )
                    if identity not in seen_descendants:
                        seen_descendants.add(identity)
                        payload.cds_stop_children.append(descendant)
        # Retain the historical field without its separate database query.
        payload.cds_children = [
            child for child in payload.cds_stop_children
            if getattr(child, "featuretype", None) == "CDS"
        ]

    # ── Phase 17b: populate the ref-side cache.
    _populate_ref_attrs_for_descent(payload, ctx)

    # Mirror the ref attrs into the legacy flat fields.
    if payload.ref_gene_id is not None:
        payload.ref_gene_attrs = payload.ref_attrs_cache.get(
            payload.ref_gene_id, {},
        )
    if payload.ref_trans_id is not None:
        payload.ref_trans_attrs = payload.ref_attrs_cache.get(
            payload.ref_trans_id, {},
        )

    # ── Phase 17b: populate the miniprot cache.
    _populate_miniprot_cache(payload, ctx)

    return payload


# ---------------------------------------------------------------------------
# Phase 17c Item 2b — thread-local DB factory for parallel materialise.
# ---------------------------------------------------------------------------

class _ThreadLocalCtxFactory:
    """Build per-thread :class:`StepContext` instances whose FeatureDB
    fields point to **fresh** ``gffutils.FeatureDB`` (or
    ``lifton.gffbase.FeatureDB``) connections opened on the calling
    thread.

    The materialise-side bottleneck (Phase 17b discovery, ~33 s on bee
    / ~84 s on rice / ~5 min on human) is SQLite query latency on a
    single shared connection from the parent thread. With each thread
    holding its own DB, the materialise step parallelises and the
    parent-thread serial bottleneck dissolves.

    On-disk databases are reopened by path. DuckDB-backed in-memory gffbase
    databases are duplicated with ``DuckDBPyConnection.cursor()``; DuckDB's
    duplicate connection shares the same in-memory catalog while giving each
    worker its own query context. Other in-memory/database-shaped objects are
    not shared: callers keep materialisation on the owning parent thread, then
    safely fan fully detached payloads out to processing workers.

    All non-DB ``StepContext`` fields (``ref_features_dict``,
    ``ref_id_2_m_id_trans_dict``, ``tree_dict``, ``tgt_fai``,
    ``ref_proteins``, ``ref_trans``, file handles, args) are copied by
    reference — they are read-only on the worker hot path or
    single-threaded mutations (file handles, but writes are short and
    GIL-protected).
    """
    __slots__ = (
        "_parent_ctx", "_ref_dbfn", "_l_dbfn", "_m_dbfn", "_local",
        "_contexts", "_contexts_lock", "_fallback_reason",
    )

    def __init__(self, parent_ctx: StepContext):
        self._parent_ctx = parent_ctx
        self._ref_dbfn = self._extract_dbfn(parent_ctx.ref_db)
        self._l_dbfn = self._extract_dbfn(parent_ctx.l_feature_db)
        self._m_dbfn = self._extract_dbfn(parent_ctx.m_feature_db)
        import threading as _t
        self._local = _t.local()
        self._contexts = []
        self._contexts_lock = _t.Lock()
        missing = []
        for label, db, dbfn in (
            ("reference", parent_ctx.ref_db, self._ref_dbfn),
            ("liftoff", parent_ctx.l_feature_db, self._l_dbfn),
            ("miniprot", parent_ctx.m_feature_db, self._m_dbfn),
        ):
            if (db is not None and not self._safe_to_share(db)
                    and dbfn is None
                    and not self._can_clone_connection(db)):
                missing.append(label)
        self._fallback_reason = (
            "non-reopenable materialisation database(s): " + ", ".join(missing)
            if missing else None
        )

    @staticmethod
    def _extract_dbfn(db) -> Optional[str]:
        """Return the on-disk DB path for ``db`` if reopen-able, else
        None (in-memory, blob, or non-FeatureDB inputs)."""
        if db is None:
            return None
        dbfn = getattr(db, "dbfn", None)
        if not isinstance(dbfn, str):
            return None
        if dbfn == ":memory:" or not dbfn:
            return None
        import os as _os
        if not _os.path.exists(dbfn):
            return None
        return dbfn

    @staticmethod
    def _safe_to_share(db) -> bool:
        """Return whether ``db`` is an immutable in-process value, not a DB.

        Mappings are used by focused embedders/tests for reference attributes
        and are safe to share.  FeatureDB-shaped objects are never returned as
        safe merely because reopening failed.
        """
        from collections.abc import Mapping as _Mapping
        return (
            isinstance(db, _Mapping)
            or (not hasattr(db, "conn") and not hasattr(db, "dbfn"))
        )

    @staticmethod
    def _can_clone_connection(db) -> bool:
        """Return whether ``db`` can make an independent thread cursor.

        ``duckdb.connect(":memory:")`` cannot be reopened by name: a new
        connection would point at a different empty catalog.  DuckDB's
        ``cursor()`` operation instead duplicates the connection while sharing
        that catalog, which is the supported shape for concurrent Python
        readers.  Limit this path to gffbase objects so a similarly named
        method on an arbitrary backend is never treated as a safety promise.
        """
        module = type(db).__module__ or ""
        conn = getattr(db, "conn", None)
        return (
            module.startswith("lifton.gffbase")
            and conn is not None
            and callable(getattr(conn, "cursor", None))
        )

    @property
    def viable(self) -> bool:
        """True iff every materialisation DB can reopen or clone per thread."""
        liftoff_reopenable = (
            self._l_dbfn is not None
            or self._safe_to_share(self._parent_ctx.l_feature_db)
            or self._can_clone_connection(self._parent_ctx.l_feature_db)
        )
        return liftoff_reopenable and self._fallback_reason is None

    @property
    def fallback_reason(self) -> Optional[str]:
        return self._fallback_reason

    @staticmethod
    def _open_thread_db(parent_db, dbfn: Optional[str]):
        """Open or clone a DB handle owned by the calling thread.

        Immutable test mappings may be shared. Database-shaped objects must
        reopen by path or, for gffbase, duplicate their DuckDB connection;
        the owning parent connection is never handed to a worker.
        """
        if parent_db is None:
            return None
        if dbfn is None:
            if _ThreadLocalCtxFactory._safe_to_share(parent_db):
                return parent_db
            if _ThreadLocalCtxFactory._can_clone_connection(parent_db):
                duplicate = None
                try:
                    duplicate = parent_db.conn.cursor()
                    return type(parent_db)(duplicate)
                except Exception as exc:
                    if duplicate is not None:
                        try:
                            duplicate.close()
                        except Exception:
                            pass
                    raise RuntimeError(
                        "could not clone in-memory gffbase connection"
                    ) from exc
            raise RuntimeError(
                "materialisation database is not reopenable; refusing to "
                "share its parent connection with a worker"
            )
        # Match the parent's class (gffutils vs gffbase).
        cls = type(parent_db)
        try:
            return cls(dbfn)
        except Exception as exc:
            raise RuntimeError(
                f"could not reopen materialisation database {dbfn!r}"
            ) from exc

    def open_root_scan(self):
        """Open the dedicated root iterator connection.

        It never services hierarchy reads, so worker queries cannot invalidate
        its cursor.  SQLite is additionally placed in query-only mode when the
        backend exposes a compatible connection.
        """
        if not self.viable:
            raise RuntimeError(self.fallback_reason or "root DB is not reopenable")
        try:
            db = self._open_thread_db(
                self._parent_ctx.l_feature_db, self._l_dbfn,
            )
        except RuntimeError:
            # Connection-less test doubles/embedders may advertise a dbfn to
            # exercise fused scheduling without owning an actual DB handle.
            # Real FeatureDB objects always expose ``conn`` and must reopen.
            if hasattr(self._parent_ctx.l_feature_db, "conn"):
                raise
            return self._parent_ctx.l_feature_db
        conn = getattr(db, "conn", None)
        if conn is not None and type(db).__module__.startswith("gffutils"):
            try:
                conn.execute("PRAGMA query_only=ON")
            except Exception:
                pass
        return db

    @staticmethod
    def close_root_scan(db) -> None:
        _ThreadLocalCtxFactory._close_db(db)

    @staticmethod
    def _close_db(db) -> None:
        if db is None or _ThreadLocalCtxFactory._safe_to_share(db):
            return
        conn = getattr(db, "conn", None)
        if conn is not None:
            try:
                conn.close()
            except Exception:
                pass

    def get(self) -> StepContext:
        """Return the calling thread's :class:`StepContext`. First
        call on a new thread opens fresh DB connections; subsequent
        calls reuse them."""
        ctx = getattr(self._local, "ctx", None)
        if ctx is not None:
            return ctx
        opened = []
        try:
            ref_db = self._open_thread_db(
                self._parent_ctx.ref_db, self._ref_dbfn,
            )
            opened.append(ref_db)
            l_feature_db = self._open_thread_db(
                self._parent_ctx.l_feature_db, self._l_dbfn,
            )
            opened.append(l_feature_db)
            m_feature_db = self._open_thread_db(
                self._parent_ctx.m_feature_db, self._m_dbfn,
            )
            opened.append(m_feature_db)
            thread_local_ctx = StepContext(
                ref_db=ref_db,
                l_feature_db=l_feature_db,
                m_feature_db=m_feature_db,
                ref_id_2_m_id_trans_dict=(
                    self._parent_ctx.ref_id_2_m_id_trans_dict
                ),
                tree_dict=self._parent_ctx.tree_dict,
                tgt_fai=self._parent_ctx.tgt_fai,
                ref_proteins=self._parent_ctx.ref_proteins,
                ref_trans=self._parent_ctx.ref_trans,
                ref_features_dict=self._parent_ctx.ref_features_dict,
                fw_score=self._parent_ctx.fw_score,
                fw_chain=self._parent_ctx.fw_chain,
                args=self._parent_ctx.args,
                state_journal=None,
                failure_records=self._parent_ctx.failure_records,
            )
            with self._contexts_lock:
                self._contexts.append(thread_local_ctx)
            self._local.ctx = thread_local_ctx
            return thread_local_ctx
        except BaseException:
            for db in reversed(opened):
                self._close_db(db)
            raise

    def close(self) -> None:
        """Close every worker-owned DB handle after its executor joins."""
        with self._contexts_lock:
            contexts, self._contexts = self._contexts, []
        parent_dbs = {
            id(self._parent_ctx.ref_db), id(self._parent_ctx.l_feature_db),
            id(self._parent_ctx.m_feature_db),
        }
        seen = set()
        for ctx in contexts:
            for db in (ctx.ref_db, ctx.l_feature_db, ctx.m_feature_db):
                if db is None or id(db) in parent_dbs or id(db) in seen:
                    continue
                seen.add(id(db))
                self._close_db(db)


def materialise_locus_with_factory(
        submission_index: int, locus,
        factory: _ThreadLocalCtxFactory) -> MaterialisedLocus:
    """Worker-side materialise: gets a thread-local ctx from the
    factory, then calls the existing :func:`materialise_locus`
    against the per-thread DB connections. The output payload is
    byte-equivalent to the parent-thread variant by construction —
    the same ``locus`` queried against the same on-disk database or shared
    in-memory DuckDB catalog gives the same children + ref attributes."""
    ctx = factory.get()
    return materialise_locus(submission_index, locus, ctx)


# ---------------------------------------------------------------------------
# Phase 11 — proxy ctx builder.
# ---------------------------------------------------------------------------

def _build_proxied_ctx(payload: MaterialisedLocus,
                       ctx: StepContext) -> StepContext:
    """Construct a worker-local :class:`StepContext` whose
    ``ref_db``, ``l_feature_db``, and ``m_feature_db`` are read-only
    proxies backed by the materialised payload caches.

    Every other field passes through unchanged: ``tree_dict`` is a
    read-only ``IntervalTree`` lookup (immutable after Step 6),
    ``ref_proteins`` / ``ref_trans`` are ``pyfaidx.Fasta`` mmaps
    (thread-safe), ``ref_features_dict`` / ``ref_id_2_m_id_trans_dict``
    are plain Python dicts (read-only on the worker hot path),
    ``fw_score`` / ``fw_chain`` are short-write text file handles
    (CPython buffered I/O is GIL-protected), and ``args`` is an
    immutable namespace.
    """
    return StepContext(
        ref_db=_RefDbProxy(payload.ref_attrs_cache),
        l_feature_db=_LFeatureDbProxy(payload.feature_cache),
        m_feature_db=(_MFeatureDbProxy(payload.miniprot_cache)
                      if ctx.m_feature_db is not None else None),
        ref_id_2_m_id_trans_dict=ctx.ref_id_2_m_id_trans_dict,
        tree_dict=ctx.tree_dict,
        tgt_fai=ctx.tgt_fai,
        ref_proteins=ctx.ref_proteins,
        ref_trans=ctx.ref_trans,
        ref_features_dict=ctx.ref_features_dict,
        fw_score=ctx.fw_score,
        fw_chain=ctx.fw_chain,
        args=ctx.args,
        state_journal=getattr(ctx, "state_journal", None),
        failure_records=getattr(ctx, "failure_records", []),
    )


def process_locus_native(payload: MaterialisedLocus,
                         ctx: StepContext) -> "LocusResult":
    """Worker-thread-safe entry point.

    Phase 17 (option b — completed Phase 12 materialisation): workers
    get a **proxied StepContext** whose ``ref_db``, ``l_feature_db``,
    and ``m_feature_db`` are read-only views of the materialised
    payload caches (built on the parent thread). The worker then
    calls the legacy ``run_liftoff.process_liftoff`` unchanged — the
    function works exactly as before, but reads come from in-process
    Python dicts instead of the shared DB connections. SQLite's
    hard-thread-binding no longer applies; the parallelism guard at
    ``lifton/parallel.py:_backend_supports_threads`` can therefore
    return True for ``--native`` regardless of backend type.

    Byte-identity is preserved by construction: the proxies return
    the same ``Feature`` objects (or attribute-bearing stubs for
    ``ref_db``), in the same order (cache built with
    ``order_by='start'``), as the real DBs would have. The legacy
    code's algorithmic path is unchanged.
    """
    from lifton import run_liftoff
    journal = getattr(ctx, "state_journal", None)
    proxy_ctx = _build_proxied_ctx(payload, ctx)
    try:
        gene = run_liftoff.process_liftoff(
            None, payload.locus,
            proxy_ctx.ref_db, proxy_ctx.l_feature_db,
            proxy_ctx.ref_id_2_m_id_trans_dict, proxy_ctx.m_feature_db,
            proxy_ctx.tree_dict,
            proxy_ctx.tgt_fai, proxy_ctx.ref_proteins, proxy_ctx.ref_trans,
            proxy_ctx.ref_features_dict,
            proxy_ctx.fw_score, proxy_ctx.fw_chain, proxy_ctx.args,
            ENTRY_FEATURE=True,
            state_journal=journal,
        )
    except Exception as exc:        # V1.4 fix: narrowed from BaseException
        return LocusResult(
            index=payload.submission_index,
            locus_id=payload.locus_id,
            lifton_gene=None,
            error=exc,
            error_tb=_traceback.format_exc(),
            delta=(journal.finish() if journal is not None else LocusDelta()),
        )
    finally:
        if journal is not None:
            journal.finish()
    return LocusResult(
        index=payload.submission_index,
        locus_id=payload.locus_id,
        lifton_gene=gene,
        delta=(journal.delta if journal is not None else LocusDelta()),
    )


def consume(result: LocusResult, fw, transcripts_stats_dict: dict, *,
            ctx: Optional[StepContext] = None,
            state_coordinator: Optional[Step7StateCoordinator] = None) -> bool:
    """Apply a :class:`LocusResult` to the output file handle and the
    stats dict. Mirrors the existing inline body in `lifton.py`'s
    Step 7 loop; centralising it here means the serial path and the
    parallel path share one implementation, removing yet another
    source of drift.

    Returns
    -------
    bool
        True iff the gene was written. Errors are logged to stderr;
        absent or invalid genes are silently skipped (Phase 5
        contract).
    """
    failure_records = (getattr(ctx, "failure_records", None)
                       if ctx is not None else None)
    if result.error is not None:
        # Mirror the inline `try/except` in lifton.py:425-426 — log
        # to stderr, swallow the error, keep going.
        #
        # Format DEFENSIVELY: a worker can package an exception whose own
        # __str__ returns a non-string (seen on full eudicot->monocot lifts:
        # `TypeError: __str__ returned non-string (type NoneType)`), and the
        # old `f"{result.error}"` then raised *inside the error handler*,
        # crashing the whole run instead of skipping one locus. Prefer the
        # traceback captured at the raise site (always a plain str); the
        # shared formatter safely falls back through str(), repr(), and the
        # exception type name.
        from lifton import logger
        detail = safe_exception_text(
            result.error, result.error_tb, prefer_traceback=True,
        )
        logger.log_error(
            f"Error during Liftoff gene processing ({result.locus_id}): "
            f"{detail}"
        )
        if failure_records is not None:
            failure_records.append({
                "stage": "step7",
                "kind": "processing_error",
                "locus_id": result.locus_id,
                "detail": detail,
            })
        return False
    if not result.emittable:
        if failure_records is not None:
            failure_records.append({
                "stage": "step7",
                "kind": "not_emittable",
                "locus_id": result.locus_id,
            })
        return False
    cds_id_allocator = (
        getattr(getattr(ctx, "args", None), "_cds_id_allocator", None)
        if ctx is not None else None
    )
    if cds_id_allocator is None:
        write_result = result.lifton_gene.write_entry(
            fw, transcripts_stats_dict
        )
    else:
        write_result = result.lifton_gene.write_entry(
            fw,
            transcripts_stats_dict,
            cds_id_allocator=cds_id_allocator,
        )
    serialization_failures = getattr(
        result.lifton_gene, "_serialization_failures", [],
    )
    if failure_records is not None:
        for feature_id, detail in serialization_failures:
            failure_records.append({
                "stage": "step7",
                "kind": "serialization_error",
                "locus_id": result.locus_id,
                "feature_id": feature_id,
                "detail": detail,
            })
        if write_result is False and not serialization_failures:
            failure_records.append({
                "stage": "step7",
                "kind": "serialization_error",
                "locus_id": result.locus_id,
            })
    # ``None`` is the historical success return.  The atomic serializer now
    # returns an explicit boolean, so only False means the model was rejected.
    if write_result is False:
        return False

    # Copy counters, suppression intervals, and sidecars describe published
    # models. Commit them only after the complete hierarchy reached the staged
    # GFF3 stream; a failed attempt must not inflate completeness statistics or
    # leave score/chain rows with no corresponding feature.
    if state_coordinator is not None:
        state_coordinator.commit(result.delta)
    if ctx is not None:
        if result.delta.score_text and ctx.fw_score is not None:
            ctx.fw_score.write(result.delta.score_text)
        if result.delta.chain_text and ctx.fw_chain is not None:
            ctx.fw_chain.write(result.delta.chain_text)
        emitted = getattr(ctx, "emitted_ref_gene_ids", None)
        ref_gene_id = getattr(result.lifton_gene, "ref_gene_id", None)
        if emitted is not None and ref_gene_id is not None:
            emitted.add(ref_gene_id)
        emitted_intervals = getattr(ctx, "emitted_intervals", None)
        if emitted_intervals is not None:
            entry = getattr(result.lifton_gene, "entry", None)
            # Step 8's overlap suppression applies only to successfully
            # published protein-coding gene loci.  Lightweight test doubles and
            # non-gene top-level features need not expose genomic coordinates.
            if entry is not None and getattr(entry, "featuretype", None) == "gene":
                emitted_intervals.append((
                    entry.seqid, int(entry.start), int(entry.end), entry.id,
                ))
    return True
