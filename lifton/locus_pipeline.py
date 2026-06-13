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
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


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
    args: Any                       # argparse.Namespace
    optimize: bool = False          # --optimize v2 lane (RESERVED; no-op today)


@dataclass
class LocusResult:
    """Result of one per-locus task. ``error`` is populated when
    `process_liftoff` raises so the parent process can log the
    failure and keep the rest of the pipeline running."""

    index: int                       # submission index (0-based)
    locus_id: str                    # gffutils-feature id, for logging
    lifton_gene: Optional[Any] = None
    error: Optional[BaseException] = None

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

    # --optimize v2 lane (RESERVED): future output-mutating accuracy/speed work
    # (smarter Liftoff/miniprot merge, banded/seeded alignment) branches from
    # ctx.optimize here. Today it is a pure pass-through, byte-identical to the
    # default path. TODO(v2): wire the alternate per-locus refinement.
    _optimize = bool(getattr(ctx, "optimize", False))  # noqa: F841

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
        )
    return LocusResult(
        index=submission_index,
        locus_id=getattr(locus, "id", "<unknown>"),
        lifton_gene=gene,
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
    # above. ``feature_cache`` covers the entire transitive descent
    # of ``locus`` through ``l_feature_db`` (so the recursion at
    # ``run_liftoff.process_liftoff:239`` is satisfied locally);
    # ``ref_attrs_cache`` carries the ref-side ``[id].attributes``
    # lookups; ``miniprot_cache`` carries every ``m_feature_db[m_id]``
    # the chaining algorithm might consult.
    feature_cache: Dict[str, _FeaturePreFetch] = field(default_factory=dict)
    ref_attrs_cache: Dict[str, dict] = field(default_factory=dict)
    miniprot_cache: Dict[str, _MiniprotPreFetch] = field(default_factory=dict)


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
    # Variant 1: featuretype='exon', level=1
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
    # Variant 2: level=1 (any type)
    try:
        entry.children_l1 = list(
            ctx.l_feature_db.children(feature, level=1)
        )
    except Exception as exc:
        logger.log_warning(
            f"_walk_and_cache_features({feature_id}): children@l1 failed: {exc}"
        )
    # Variant 3: featuretype='exon' (no level, used by lifton_add_trans_exon_cds).
    #
    # Phase 17c Win 1: dedup against Variant 1 when the feature has any
    # level-1 exon children. In standard GFF3 (gene → mRNA → exon, with
    # exons as leaves), the per-locus body only invokes Variant 3 on
    # transcript-shaped features (i.e. those reached via the
    # ``len(exon_children) > 0`` branch at run_liftoff.py:240) — and on
    # those, all exons are at level=1 by the leaf-exon convention. So
    # Variant 3 ⊇ Variant 1 holds with equality, and the second query
    # is redundant. The fallback path (run a real query when there are
    # no level-1 exons) is preserved for gene-level features where the
    # proxy might still be asked for Variant 3 defensively.
    if entry.exon_children_l1:
        entry.exon_children_full = list(entry.exon_children_l1)
    else:
        try:
            entry.exon_children_full = list(
                ctx.l_feature_db.children(
                    feature, featuretype="exon", order_by="start",
                )
            )
        except Exception as exc:
            logger.log_warning(
                f"_walk_and_cache_features({feature_id}): exon (no-level) failed: {exc}"
            )
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
    _maybe_cache_ref_attrs(payload, ctx.ref_db, ref_gene_id)
    _maybe_cache_ref_attrs(payload, ctx.ref_db, ref_trans_id_for_gene)

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
        _maybe_cache_ref_attrs(payload, ctx.ref_db, r_gene_id)
        _maybe_cache_ref_attrs(payload, ctx.ref_db, r_trans_id)


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
    candidate_ref_trans_ids = set(payload.ref_attrs_cache.keys())
    for r_trans_id in candidate_ref_trans_ids:
        m_ids = ctx.ref_id_2_m_id_trans_dict.get(r_trans_id) or []
        for m_id in m_ids:
            if m_id in payload.miniprot_cache:
                continue
            try:
                m_entry = ctx.m_feature_db[m_id]
            except Exception as exc:
                logger.log_warning(
                    f"_populate_miniprot_cache({m_id}): lookup failed: {exc}"
                )
                continue
            try:
                cds_stop = list(ctx.m_feature_db.children(
                    m_entry, featuretype=("CDS", "stop_codon"),
                    order_by="start",
                ))
            except Exception as exc:
                logger.log_warning(
                    f"_populate_miniprot_cache({m_id}): children failed: {exc}"
                )
                cds_stop = []
            payload.miniprot_cache[m_id] = _MiniprotPreFetch(
                feature=m_entry, cds_stop_children=cds_stop,
            )


def materialise_locus(submission_index: int, locus,
                      ctx: StepContext) -> MaterialisedLocus:
    """Parent-thread-only function. Pre-fetches **everything**
    ``run_liftoff.process_liftoff`` would have read from
    ``ctx.l_feature_db``, ``ctx.ref_db``, and ``ctx.m_feature_db`` so
    the worker thread is purely CPU + payload reads.

    Phase 17 (option b): walks the locus's transitive descent and
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

    # ── Phase 17b: walk + cache the entire l_feature_db tree under locus.
    try:
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
        # ``cds_children`` (CDS-only, no stop_codon) is a separate
        # historical query the test suite asserts on; pre-fetch it
        # explicitly so legacy tests keep passing.
        try:
            payload.cds_children = list(
                ctx.l_feature_db.children(
                    locus, featuretype="CDS", order_by="start",
                )
            )
        except Exception as exc:
            logger.log_warning(
                f"materialise_locus({locus_id}): cds_children "
                f"(legacy) failed: {exc}"
            )

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

    The factory only constructs a non-trivial path when the DBs in the
    parent ``ctx`` advertise an on-disk ``dbfn``. For in-memory DBs
    (gffbase blob path) or DBs without a path, it returns ``None``
    (caller falls back to the legacy parent-thread materialise loop —
    correctness preserved, just no parallelism speedup).

    All non-DB ``StepContext`` fields (``ref_features_dict``,
    ``ref_id_2_m_id_trans_dict``, ``tree_dict``, ``tgt_fai``,
    ``ref_proteins``, ``ref_trans``, file handles, args) are copied by
    reference — they are read-only on the worker hot path or
    single-threaded mutations (file handles, but writes are short and
    GIL-protected).
    """
    __slots__ = ("_parent_ctx", "_ref_dbfn", "_l_dbfn", "_m_dbfn",
                 "_local")

    def __init__(self, parent_ctx: StepContext):
        self._parent_ctx = parent_ctx
        self._ref_dbfn = self._extract_dbfn(parent_ctx.ref_db)
        self._l_dbfn = self._extract_dbfn(parent_ctx.l_feature_db)
        self._m_dbfn = self._extract_dbfn(parent_ctx.m_feature_db)
        import threading as _t
        self._local = _t.local()

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

    @property
    def viable(self) -> bool:
        """True iff at least the Liftoff FeatureDB has a reopen-able
        path. ``ref_db`` and ``m_feature_db`` may be None / in-memory
        and we still benefit from parallelising the l_feature_db
        queries (the dominant cost)."""
        return self._l_dbfn is not None

    @staticmethod
    def _open_thread_db(parent_db, dbfn: Optional[str]):
        """Open a fresh DB connection for the calling thread. Returns
        the parent ``parent_db`` unchanged when no path is available,
        which preserves correctness for in-memory backends (their
        cross-thread reads are safe — only on-disk SQLite hard-binds
        connections)."""
        if parent_db is None:
            return None
        if dbfn is None:
            # In-memory backends or non-FeatureDB inputs — fall back
            # to the parent's instance. Safe under our usage (the only
            # backends we share are dict-like or blob-built FeatureDBs
            # that don't hold thread-bound cursors).
            return parent_db
        # Match the parent's class (gffutils vs gffbase).
        cls = type(parent_db)
        try:
            return cls(dbfn)
        except Exception:
            # Fall back to the parent's instance if we can't reopen
            # — still correct (legacy single-thread materialise).
            return parent_db

    def get(self) -> StepContext:
        """Return the calling thread's :class:`StepContext`. First
        call on a new thread opens fresh DB connections; subsequent
        calls reuse them."""
        ctx = getattr(self._local, "ctx", None)
        if ctx is not None:
            return ctx
        thread_local_ctx = StepContext(
            ref_db=self._open_thread_db(
                self._parent_ctx.ref_db, self._ref_dbfn),
            l_feature_db=self._open_thread_db(
                self._parent_ctx.l_feature_db, self._l_dbfn),
            m_feature_db=self._open_thread_db(
                self._parent_ctx.m_feature_db, self._m_dbfn),
            ref_id_2_m_id_trans_dict=self._parent_ctx.ref_id_2_m_id_trans_dict,
            tree_dict=self._parent_ctx.tree_dict,
            tgt_fai=self._parent_ctx.tgt_fai,
            ref_proteins=self._parent_ctx.ref_proteins,
            ref_trans=self._parent_ctx.ref_trans,
            ref_features_dict=self._parent_ctx.ref_features_dict,
            fw_score=self._parent_ctx.fw_score,
            fw_chain=self._parent_ctx.fw_chain,
            args=self._parent_ctx.args,
        )
        self._local.ctx = thread_local_ctx
        return thread_local_ctx


def materialise_locus_with_factory(
        submission_index: int, locus,
        factory: _ThreadLocalCtxFactory) -> MaterialisedLocus:
    """Worker-side materialise: gets a thread-local ctx from the
    factory, then calls the existing :func:`materialise_locus`
    against the per-thread DB connections. The output payload is
    byte-equivalent to the parent-thread variant by construction —
    the same ``locus`` queried against the same on-disk SQLite gives
    the same children + ref attributes."""
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
        )
    except Exception as exc:        # V1.4 fix: narrowed from BaseException
        return LocusResult(
            index=payload.submission_index,
            locus_id=payload.locus_id,
            lifton_gene=None,
            error=exc,
        )
    return LocusResult(
        index=payload.submission_index,
        locus_id=payload.locus_id,
        lifton_gene=gene,
    )


def consume(result: LocusResult, fw, transcripts_stats_dict: dict) -> bool:
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
    if result.error is not None:
        # Mirror the inline `try/except` in lifton.py:425-426 — log
        # to stderr, swallow the error, keep going.
        from lifton import logger
        logger.log_error(
            f"Error during Liftoff gene processing ({result.locus_id}): "
            f"{result.error}"
        )
        return False
    if not result.emittable:
        return False
    result.lifton_gene.write_entry(fw, transcripts_stats_dict)
    return True
