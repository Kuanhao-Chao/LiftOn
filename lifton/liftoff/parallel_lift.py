"""Process-parallel lift loop for the vendored Liftoff (v1.0.12).

``lift_features.lift_all_features`` is GIL-bound Python (the alignment DAG in
``find_best_mapping``, coordinate conversion, feature merging), so threads
cannot speed it up. Processes can, because without overlap intervals
(``feature_locations is None``: the primary pass and the unmapped-genes pass)
the loop has exactly one dependency between features. That dependency is the
neighbour hint from ``find_neighbor_location``, and the hint never crosses a
reference chromosome: ``liftoff_utils.find_nonoverlapping_upstream_neighbor``
returns ``None`` whenever the upstream neighbour lies on another sequence.

Each reference chromosome's alignments are therefore lifted, in the serial
order, by one forked worker. The worker looks neighbours up first in what it
has lifted so far and then in what earlier passes lifted, which is exactly
what the serial loop can see. Results merge back in chromosome order, the
serial insertion order, and unmapped features are mapped back to the original
parent objects. The output is identical to the serial loop by construction.

Fork is used so the read-only inputs (feature hierarchy, parent order, the
earlier lifted features) are shared copy-on-write rather than pickled; only
each batch of alignments goes to a worker and only its results come back.
Workers never touch the reference database: the parent computes the feature
order before forking.

On by default whenever ``--threads`` is greater than 1; ``--no-parallel-lift``
or ``LIFTON_PARALLEL_LIFT=0`` runs the serial loop. Byte-identical on the
drosophila and dog -> cat whole genomes (aligner phase -33 % and -32 % at
``-t 8``).
"""
from __future__ import annotations

import multiprocessing
import os
from collections import ChainMap

from lifton.liftoff import lift_features, liftoff_utils

PARALLEL_LIFT_DEFAULT = True

# Read-only state for the forked workers, set in the parent just before the
# pool is created and cleared right after.
_SHARED: dict = {}


def enabled(args) -> bool:
    env = os.environ.get("LIFTON_PARALLEL_LIFT")
    if env is not None:
        return env.strip().lower() not in ("0", "", "false", "no")
    value = getattr(args, "parallel_lift", None)
    return PARALLEL_LIFT_DEFAULT if value is None else bool(value)


def _worker_budget(args) -> int:
    override = os.environ.get("LIFTON_PARALLEL_LIFT_WORKERS")
    if override:
        return max(1, int(override))
    return max(1, int(getattr(args, "threads", 1) or 1))


def group_by_reference_chromosome(alignments, features_to_lift):
    """Split alignments (already in reference ``(seqid, start)`` order) into
    contiguous same-reference-chromosome batches."""
    batches, current, current_seqid = [], [], None
    for alignment in alignments:
        parent = features_to_lift[
            liftoff_utils.convert_id_to_original(alignment[0].query_name)]
        if parent.seqid != current_seqid and current:
            batches.append(current)
            current = []
        current_seqid = parent.seqid
        current.append(alignment)
    if current:
        batches.append(current)
    return batches


def _reset_inherited_signal_handlers():
    """A forked worker must not run the parent's SIGTERM/SIGHUP handlers
    (``Pool.terminate`` signals its workers); restore the defaults."""
    import signal
    for signum in (signal.SIGTERM, getattr(signal, "SIGHUP", None)):
        if signum is not None:
            signal.signal(signum, signal.SIG_DFL)


def _lift_batch(batch):
    shared = _SHARED
    hierarchy = shared["feature_hierarchy"]
    lifted = {}
    lookup = ChainMap(lifted, shared["earlier"])
    unmapped = []
    for alignment in batch:
        start, seqid, ref_start = lift_features.find_neighbor_location(
            hierarchy.parents, alignment, lookup, shared["ref_parent_order"])
        features, name = lift_features.lift_single_feature(
            shared["threshold"], shared["feature_order"], hierarchy.parents,
            hierarchy, start, ref_start, seqid, unmapped, alignment,
            shared["seq_id_threshold"], None, lookup, shared["args"])
        if features != []:
            lifted[name] = features
    return lifted, unmapped


def lift_all_features_parallel(alns, threshold, feature_db, feature_hierarchy,
                               unmapped_features, lifted_feature_list,
                               seq_id_threshold, feature_locations, args,
                               ref_parent_order):
    """Parallel twin of ``lift_features.lift_all_features``.

    Returns False without doing anything when the serial loop must run:
    overlap intervals are in play, or there is at most one batch or worker.
    """
    if feature_locations is not None or _worker_budget(args) <= 1:
        return False
    features_to_lift = feature_hierarchy.parents
    alignments = lift_features.sort_alignments(features_to_lift, alns)
    batches = group_by_reference_chromosome(alignments, features_to_lift)
    workers = min(_worker_budget(args), len(batches))
    if workers <= 1:
        return False
    _SHARED.update(
        feature_hierarchy=feature_hierarchy,
        feature_order=lift_features.get_feature_order(feature_db),
        threshold=threshold, seq_id_threshold=seq_id_threshold,
        ref_parent_order=ref_parent_order, args=args,
        earlier=dict(lifted_feature_list),
    )
    try:
        try:
            with multiprocessing.get_context("fork").Pool(
                    workers, initializer=_reset_inherited_signal_handlers) as pool:
                results = pool.map(_lift_batch, batches, chunksize=1)
        except OSError as error:
            # fork() can fail even with plenty of RAM free. Under strict overcommit
            # accounting (vm.overcommit_memory=2, no swap) the kernel reserves the
            # parent's whole address space for each child with no copy-on-write
            # credit, so a large parent forking many workers exceeds CommitLimit:
            # observed here as ENOMEM from a ~35 GB parent at 32 workers while
            # ~930 GB was physically free. The serial path below produces the same
            # result, so degrade to it instead of losing a whole-genome run.
            from lifton import logger
            logger.log_warning(
                f"Parallel lift could not start {workers} worker(s) ({error}); "
                f"falling back to the serial lift. Set LIFTON_PARALLEL_LIFT_WORKERS "
                f"to a smaller value to keep the parallel path."
            )
            return False
    finally:
        _SHARED.clear()
    for lifted, unmapped in results:
        lifted_feature_list.update(lifted)
        unmapped_features.extend(
            features_to_lift.get(feature.id, feature) for feature in unmapped)
    return True
