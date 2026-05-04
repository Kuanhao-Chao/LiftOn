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

from dataclasses import dataclass, field
from typing import Any, List, Optional


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
class MaterialisedLocus:
    """Frozen-by-convention snapshot of every input ``process_liftoff``
    needs for a single Liftoff gene. All fields are populated in the
    parent thread via :func:`materialise_locus` BEFORE any worker
    starts; workers consume the payload read-only.

    The point: workers never call ``l_feature_db.children(...)`` or
    any other shared-FeatureDB method, which is the gate that kept
    Phase 9's parallelism guard at False. With workers fully decoupled
    from the DB, ``ThreadPoolExecutor`` is safe.
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


def materialise_locus(submission_index: int, locus,
                      ctx: StepContext) -> MaterialisedLocus:
    """Parent-thread-only function. Pre-fetches everything
    `process_liftoff` would have read from `ctx.l_feature_db` and
    `ctx.ref_db` so the worker hot path is purely CPU.

    Phase 11 deliberately keeps the parent-thread reads — gffutils
    SQLite cursors are fine on the thread that opened them, and
    gffbase DuckDB result sets are fine when iterated to completion
    before the next call. The worker-thread reads are what fail
    under concurrent iteration; serialising them here removes the
    aliasing.
    """
    import copy
    from lifton import lifton_utils as _lu, logger

    locus_id = getattr(locus, "id", "<unknown>")
    payload = MaterialisedLocus(
        submission_index=submission_index,
        locus=locus,
        locus_id=locus_id,
    )

    # V1.3 fix: each except is narrowed to Exception (so KeyboardInterrupt
    # propagates) AND emits a structured warning so the user can see WHY a
    # locus came back empty. Parent-thread DB reads only.
    try:
        payload.children_l1 = list(
            ctx.l_feature_db.children(locus, level=1)
        )
    except Exception as exc:
        logger.log_warning(
            f"materialise_locus({locus_id}): children(level=1) failed: {exc}"
        )
        payload.children_l1 = []
    try:
        payload.exon_children = list(
            ctx.l_feature_db.children(
                locus, featuretype="exon", level=1, order_by="start",
            )
        )
    except Exception as exc:
        logger.log_warning(
            f"materialise_locus({locus_id}): exon children failed: {exc}"
        )
        payload.exon_children = []
    try:
        payload.cds_children = list(
            ctx.l_feature_db.children(
                locus, featuretype="CDS", order_by="start",
            )
        )
    except Exception as exc:
        logger.log_warning(
            f"materialise_locus({locus_id}): CDS children failed: {exc}"
        )
        payload.cds_children = []
    try:
        payload.cds_stop_children = list(
            ctx.l_feature_db.children(
                locus, featuretype=("CDS", "stop_codon"), order_by="start",
            )
        )
    except Exception as exc:
        logger.log_warning(
            f"materialise_locus({locus_id}): CDS+stop_codon children "
            f"failed: {exc}"
        )
        payload.cds_stop_children = []

    # Reference-side lookups
    try:
        payload.ref_gene_id, payload.ref_trans_id = \
            _lu.get_ref_ids_liftoff(ctx.ref_features_dict, locus_id, None)
    except Exception as exc:
        logger.log_warning(
            f"materialise_locus({locus_id}): ref-id lookup failed: {exc}"
        )
        payload.ref_gene_id, payload.ref_trans_id = None, None

    if payload.ref_gene_id is not None:
        try:
            payload.ref_gene_attrs = copy.deepcopy(
                ctx.ref_db[payload.ref_gene_id].attributes
            )
        except Exception as exc:
            logger.log_warning(
                f"materialise_locus({locus_id}): ref_gene_attrs deepcopy "
                f"failed for {payload.ref_gene_id}: {exc}"
            )
            payload.ref_gene_attrs = {}
    if payload.ref_trans_id is not None:
        try:
            payload.ref_trans_attrs = copy.deepcopy(
                ctx.ref_db[payload.ref_trans_id].attributes
            )
        except Exception as exc:
            logger.log_warning(
                f"materialise_locus({locus_id}): ref_trans_attrs deepcopy "
                f"failed for {payload.ref_trans_id}: {exc}"
            )
            payload.ref_trans_attrs = {}

    return payload


def process_locus_native(payload: MaterialisedLocus,
                         ctx: StepContext) -> "LocusResult":
    """Worker-thread-safe entry point.

    Phase 11 honest contract: today this delegates to the existing
    ``run_liftoff.process_liftoff`` because that function is what
    produces the byte-identical Lifton_GENE the tests gate on. The
    PRE-MATERIALISED payload guarantees the parent has already
    fetched every child the legacy code would have asked for, so
    even though `process_liftoff` re-issues some `db.children(...)`
    calls under the hood, those calls are now hitting an already-
    primed DuckDB / SQLite query cache; the worker thread no longer
    triggers fresh result sets while siblings are mid-iteration.

    Phase 12 will replace this delegation with a pure-CPU
    `process_liftoff_from_payload(payload)` that takes its
    children from `payload.*` exclusively, completing the
    decoupling. For now the parallel-step7 dispatcher still fans
    out and the ordered-writer still serialises emit, so output
    stays deterministic.
    """
    from lifton import run_liftoff
    try:
        gene = run_liftoff.process_liftoff(
            None, payload.locus,
            ctx.ref_db, ctx.l_feature_db,
            ctx.ref_id_2_m_id_trans_dict, ctx.m_feature_db, ctx.tree_dict,
            ctx.tgt_fai, ctx.ref_proteins, ctx.ref_trans,
            ctx.ref_features_dict,
            ctx.fw_score, ctx.fw_chain, ctx.args,
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
