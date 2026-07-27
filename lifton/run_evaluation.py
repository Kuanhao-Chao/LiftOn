"""Evaluation-mode helpers.

Evaluation consumes an existing target annotation and writes per-transcript
identity/status rows.  It must not apply the copy-number or interval-tree
mutations used by the annotation-producing pipeline.
"""

import copy
import io
import traceback
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any, Iterable, Iterator

import gffutils

from lifton import lifton_class, lifton_utils


_MISSING_FEATURE_ERRORS = (KeyError, gffutils.exceptions.FeatureNotFoundError)


@dataclass(frozen=True)
class EvaluationPayload:
    """Self-contained target/reference hierarchy for one evaluation root."""

    index: int
    locus: Any
    target_features: dict[str, Any]
    target_children: dict[str, tuple[Any, ...]]
    reference_features: dict[str, Any]


@dataclass
class EvaluationResult:
    """One ordered evaluation result with buffered score output."""

    index: int
    locus_id: str
    score_text: str = ""
    evaluated: bool = False
    error: BaseException | None = None
    error_tb: str | None = None


class _FeatureLookupProxy:
    """Minimal read-only ``FeatureDB[id]`` proxy backed by copied features."""

    def __init__(self, features):
        self._features = features

    def __getitem__(self, feature_id):
        try:
            return self._features[str(feature_id)]
        except KeyError as exc:
            raise KeyError(feature_id) from exc


class _TargetHierarchyProxy(_FeatureLookupProxy):
    """Read-only target DB view supporting evaluation's children queries."""

    def __init__(self, features, children):
        super().__init__(features)
        self._children = children

    @staticmethod
    def _feature_types(value):
        if value is None:
            return None
        if isinstance(value, str):
            return {value}
        return set(value)

    def children(self, parent, featuretype=None, level=None, order_by=None):
        parent_id = getattr(parent, "id", parent)
        wanted = self._feature_types(featuretype)
        direct = list(self._children.get(str(parent_id), ()))
        if level == 1:
            candidates = direct
        else:
            candidates = []
            pending = list(reversed(direct))
            while pending:
                child = pending.pop()
                candidates.append(child)
                pending.extend(reversed(self._children.get(str(child.id), ())))
        if wanted is not None:
            candidates = [item for item in candidates if item.featuretype in wanted]
        if order_by == "start":
            candidates.sort(key=_feature_sort_key)
        return iter(candidates)


def _chm13_mode(args):
    return bool(getattr(args, "evaluation_liftoff_chm13", False))


def _logical_id(feature_id, prefix):
    """Return the unprefixed ID used by evaluation's sequence dictionaries."""
    if prefix and feature_id.startswith(prefix):
        return feature_id[len(prefix):]
    return feature_id


def _copy_base(feature_id):
    """Return a possible LiftOn/Liftoff copy base, or ``None``.

    The base is only a candidate.  Lookup always tries the exact ID first, so
    natural reference IDs such as ``TX_1`` cannot be mistaken for copies.
    """
    base, separator, suffix = feature_id.rpartition("_")
    if separator and base and suffix.isdigit():
        return base
    return None


def _lookup_reference_feature(ref_db, feature_id, prefix="", allow_copy_suffix=True):
    """Resolve a reference feature without swallowing unexpected failures.

    CHM13 evaluation presents IDs without RefSeq's ``gene-``/``rna-``
    prefixes.  Copied target annotations may additionally append ``_<n>``.
    Exact lookup precedes copy-suffix fallback to preserve genuine numeric
    suffixes in reference IDs.
    """
    if feature_id is None:
        return None, None

    logical_id = _logical_id(str(feature_id), prefix)
    candidates = [logical_id]
    copy_base = _copy_base(logical_id) if allow_copy_suffix else None
    if copy_base is not None:
        candidates.append(copy_base)

    for candidate in candidates:
        key = f"{prefix}{candidate}" if prefix else candidate
        try:
            return candidate, ref_db[key]
        except _MISSING_FEATURE_ERRORS:
            continue
    return None, None


def _isolated_ref_features(ref_features_dict, ref_gene_id):
    """Clone only mutable copy-counter state needed by ``Lifton_GENE``."""
    isolated = dict(ref_features_dict)
    if ref_gene_id in isolated:
        isolated[ref_gene_id] = copy.copy(isolated[ref_gene_id])
    return isolated


def _feature_sort_key(feature):
    return (
        str(getattr(feature, "seqid", "")),
        int(getattr(feature, "start", 0)),
        int(getattr(feature, "end", 0)),
        str(getattr(feature, "featuretype", "")),
        str(getattr(feature, "id", "")),
    )


def _children(tgt_db, locus, *, featuretype=None, level=None):
    kwargs = {"order_by": "start"}
    if featuretype is not None:
        kwargs["featuretype"] = featuretype
    if level is not None:
        kwargs["level"] = level
    return sorted(list(tgt_db.children(locus, **kwargs)), key=_feature_sort_key)


def _copy_target_hierarchy(tgt_db, root):
    """Copy one hierarchy while the owning database connection is local."""
    features = {}
    children = {}
    pending = [root]
    while pending:
        feature = pending.pop()
        feature_id = str(feature.id)
        if feature_id in features:
            continue
        feature_copy = copy.deepcopy(feature)
        features[feature_id] = feature_copy
        direct = list(tgt_db.children(feature, level=1, order_by="start"))
        direct.sort(key=_feature_sort_key)
        copied_children = tuple(copy.deepcopy(child) for child in direct)
        children[feature_id] = copied_children
        pending.extend(reversed(direct))
    return features, children


def _reference_keys_for_payload(root, target_features, ref_features_dict, args):
    """Return the small set of reference IDs evaluation may look up."""
    keys = set()
    ref_gene_id, _ = lifton_utils.get_ref_ids_liftoff(
        ref_features_dict, root.id, None,
    )
    if ref_gene_id is not None:
        prefix = "gene-" if _chm13_mode(args) else ""
        logical = _logical_id(str(ref_gene_id), prefix)
        keys.add(f"{prefix}{logical}")
        base = _copy_base(logical)
        if base:
            keys.add(f"{prefix}{base}")

    transcript_prefix = "rna-" if _chm13_mode(args) else ""
    for feature in target_features.values():
        if feature.featuretype in {"gene", "exon", "CDS", "stop_codon"}:
            continue
        logical = _logical_id(str(feature.id), transcript_prefix)
        keys.add(f"{transcript_prefix}{logical}")
        base = _copy_base(logical)
        if base:
            keys.add(f"{transcript_prefix}{base}")
    return keys


def materialize_evaluation_payload(index, locus, ref_db, tgt_db,
                                   ref_features_dict, args):
    """Create a worker-safe evaluation payload on the DB-owning thread."""
    target_features, target_children = _copy_target_hierarchy(tgt_db, locus)
    reference_features = {}
    for key in sorted(_reference_keys_for_payload(
            locus, target_features, ref_features_dict, args)):
        try:
            reference_features[key] = copy.deepcopy(ref_db[key])
        except _MISSING_FEATURE_ERRORS:
            continue
    root_copy = target_features[str(locus.id)]
    return EvaluationPayload(
        index=index,
        locus=root_copy,
        target_features=target_features,
        target_children=target_children,
        reference_features=reference_features,
    )


def evaluate_payload(payload, *, tree_dict, tgt_fai, ref_proteins,
                     ref_trans, ref_features_dict, args):
    """Evaluate one materialized hierarchy without touching shared DB/output."""
    score = io.StringIO()
    try:
        gene = evaluation(
            None,
            payload.locus,
            _FeatureLookupProxy(payload.reference_features),
            _TargetHierarchyProxy(
                payload.target_features, payload.target_children,
            ),
            tree_dict,
            tgt_fai,
            ref_proteins,
            ref_trans,
            ref_features_dict,
            score,
            args,
            ENTRY_FEATURE=True,
        )
        return EvaluationResult(
            index=payload.index,
            locus_id=str(payload.locus.id),
            score_text=score.getvalue(),
            evaluated=gene is not None,
        )
    except Exception as exc:
        return EvaluationResult(
            index=payload.index,
            locus_id=str(getattr(payload.locus, "id", "<unknown>")),
            error=exc,
            error_tb=traceback.format_exc(),
        )


def parallel_evaluate_loci(loci: Iterable[Any], ref_db, tgt_db, tree_dict,
                           tgt_fai, ref_proteins, ref_trans,
                           ref_features_dict, fw_score, args, *, threads=1,
                           max_inflight=None) -> tuple[int, list[EvaluationResult]]:
    """Evaluate root loci with bounded workers and ordered score publication."""
    threads = max(1, int(threads or 1))
    max_inflight = max(1, int(max_inflight or (2 * threads)))

    def payloads() -> Iterator[EvaluationPayload]:
        for index, locus in enumerate(loci):
            yield materialize_evaluation_payload(
                index, locus, ref_db, tgt_db, ref_features_dict, args,
            )

    worker_kwargs = {
        "tree_dict": tree_dict,
        "tgt_fai": tgt_fai,
        "ref_proteins": ref_proteins,
        "ref_trans": ref_trans,
        "ref_features_dict": ref_features_dict,
        "args": args,
    }

    def worker(payload):
        return evaluate_payload(payload, **worker_kwargs)

    failures = []
    processed = 0
    observed_high_water = 0

    def record_high_water(value):
        nonlocal observed_high_water
        observed_high_water = max(observed_high_water, value)

    if threads <= 1:
        results = map(worker, payloads())
        executor = None
        observed_high_water = 1
    else:
        from lifton.parallel import _ordered_bounded_map

        executor = ThreadPoolExecutor(
            max_workers=threads, thread_name_prefix="lifton-evaluation",
        )
        results = _ordered_bounded_map(
            executor,
            payloads(),
            lambda pool, payload: pool.submit(worker, payload),
            max_inflight,
            high_water_callback=record_high_water,
        )
    try:
        for result in results:
            processed += 1
            if result.error is not None:
                failures.append(result)
                continue
            fw_score.write(result.score_text)
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    try:
        setattr(args, "_evaluation_max_inflight_observed",
                observed_high_water if processed else 0)
    except Exception:
        pass
    return processed, failures


def initialize_lifton_gene_eval(locus, ref_db, tree_dict, ref_features_dict,
                                args, with_exons=False):
    """Initialize an evaluation-only gene without mutating pipeline state."""
    ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(
        ref_features_dict, locus.id, None,
    )
    if ref_gene_id is None:
        return None, None, None

    prefix = "gene-" if _chm13_mode(args) else ""
    ref_gene_id, ref_feature = _lookup_reference_feature(
        ref_db, ref_gene_id, prefix=prefix,
    )
    if ref_feature is None:
        return None, None, None

    # Lifton_GENE normally increments a global copy counter and inserts into a
    # shared interval tree.  Neither side effect belongs in read-only
    # evaluation, so give it isolated state for both.
    isolated_features = _isolated_ref_features(ref_features_dict, ref_gene_id)
    lifton_gene = lifton_class.Lifton_GENE(
        ref_gene_id,
        copy.deepcopy(locus),
        copy.deepcopy(ref_feature.attributes),
        {},
        isolated_features,
        args,
        tmp=with_exons,
    )
    return lifton_gene, ref_gene_id, ref_trans_id


def lifton_add_trans_exon_cds_eval(lifton_gene, locus, ref_db, tgt_db,
                                   ref_trans_id, args):
    """Add one target transcript hierarchy for evaluation.

    The target database features are deep-copied because LiftOn constructors
    rewrite IDs and Parent attributes.  If a valid target has CDS records but
    omits explicit exons, CDS spans are used as evaluation-only exon spans.
    """
    prefix = "rna-" if _chm13_mode(args) else ""
    ref_trans_id, ref_feature = _lookup_reference_feature(
        ref_db, ref_trans_id, prefix=prefix,
    )
    if ref_feature is None:
        return None, 0

    lifton_trans = lifton_gene.add_transcript(
        ref_trans_id,
        copy.deepcopy(locus),
        copy.deepcopy(ref_feature.attributes),
    )

    exons = _children(tgt_db, locus, featuretype="exon")
    cds_features = _children(
        tgt_db, locus, featuretype=("CDS", "stop_codon"),
    )
    if not exons:
        # GFF3 permits CDS children without separate exon records.  LiftOn's
        # sequence extraction operates on exon objects, so create transient
        # exon views from CDS segments without altering the target database.
        for cds in _children(tgt_db, locus, featuretype="CDS"):
            exon = copy.deepcopy(cds)
            exon.featuretype = "exon"
            exon.frame = "."
            lifton_gene.add_exon(lifton_trans.entry.id, exon)
    else:
        for exon in exons:
            lifton_gene.add_exon(lifton_trans.entry.id, copy.deepcopy(exon))

    for cds in cds_features:
        lifton_gene.add_cds(lifton_trans.entry.id, copy.deepcopy(cds))
    return lifton_trans, len(cds_features)


def evaluation(lifton_gene, locus, ref_db, tgt_db, tree_dict, tgt_fai,
               ref_proteins, ref_trans, ref_features_dict, fw_score, args,
               ENTRY_FEATURE=False):
    """Evaluate every transcript below ``locus`` and return its LiftOn gene.

    Unmapped target genes/transcripts are skipped.  Unexpected database or
    alignment errors still propagate so real input/tool failures remain
    visible to callers.
    """
    direct_exons = _children(tgt_db, locus, featuretype="exon", level=1)
    direct_cdss = _children(tgt_db, locus, featuretype="CDS", level=1)
    has_transcript_structure = bool(direct_exons or direct_cdss)

    ref_gene_id = None
    if lifton_gene is None:
        if not ENTRY_FEATURE:
            return None
        lifton_gene, ref_gene_id, _ = initialize_lifton_gene_eval(
            locus,
            ref_db,
            tree_dict,
            ref_features_dict,
            args,
            with_exons=has_transcript_structure,
        )
        if lifton_gene is None or lifton_gene.ref_gene_id is None:
            return None

    if not has_transcript_structure:
        for feature in _children(tgt_db, locus, level=1):
            if feature.featuretype in {"exon", "CDS", "stop_codon"}:
                continue
            evaluation(
                lifton_gene,
                feature,
                ref_db,
                tgt_db,
                tree_dict,
                tgt_fai,
                ref_proteins,
                ref_trans,
                ref_features_dict,
                fw_score,
                args,
            )
        return lifton_gene

    if ENTRY_FEATURE:
        ref_trans_id = ref_gene_id or lifton_gene.ref_gene_id
    else:
        _, ref_trans_id = lifton_utils.get_ref_ids_liftoff(
            ref_features_dict, lifton_gene.ref_gene_id, locus.id,
        )

    lifton_status = lifton_class.Lifton_Status()
    lifton_trans, _ = lifton_add_trans_exon_cds_eval(
        lifton_gene, locus, ref_db, tgt_db, ref_trans_id, args,
    )
    if lifton_trans is None:
        return lifton_gene

    canonical_ref_trans_id = lifton_trans.ref_tran_id
    lifton_gene.orf_search_protein(
        lifton_trans.entry.id,
        canonical_ref_trans_id,
        tgt_fai,
        ref_proteins,
        ref_trans,
        lifton_status,
        eval_only=True,
        eval_liftoff_chm13=_chm13_mode(args),
    )
    lifton_utils.print_lifton_status(
        lifton_trans.entry.id, locus, lifton_status,
        DEBUG=bool(getattr(args, "debug", False)),
    )
    lifton_utils.write_lifton_eval_status(
        fw_score, lifton_trans.entry.id, locus, lifton_status,
    )
    return lifton_gene
