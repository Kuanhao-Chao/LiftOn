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
from typing import Any, Optional


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
    except BaseException as exc:    # noqa: BLE001 — package + continue
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
