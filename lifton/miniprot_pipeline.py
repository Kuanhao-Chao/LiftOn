"""Bounded deterministic Step-8 candidate analysis and publication.

Workers receive copied miniprot/reference records and isolated copy state.  The
parent rechecks overlap, assigns the real copy number, serializes, and commits
state in input order.  This preserves serial LiftOn semantics while allowing
the alignment/ORF work to overlap across candidates.
"""

from __future__ import annotations

import copy
import io
import traceback
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Any, Iterable, Iterator

from intervaltree import Interval
import gffutils

from lifton import lifton_utils, run_miniprot
from lifton.locus_pipeline import (
    DeferredStateJournal,
    commit_locus_delta,
    safe_exception_text,
)


@dataclass(frozen=True)
class MiniprotPayload:
    index: int
    transcript: Any
    children: tuple[Any, ...]
    reference_features: dict[str, Any]
    error: BaseException | None = None
    error_tb: str | None = None


@dataclass
class MiniprotAnalysis:
    index: int
    transcript: Any
    gene: Any | None = None
    score_text: str = ""
    explicit_copy_number: int | None = None
    error: BaseException | None = None
    error_tb: str | None = None


@dataclass
class Step8Outcome:
    processed: int = 0
    emitted: int = 0
    emitted_ref_gene_ids: set[str] | None = None
    failures: list[dict] | None = None
    max_inflight_observed: int = 0

    def __post_init__(self):
        if self.emitted_ref_gene_ids is None:
            self.emitted_ref_gene_ids = set()
        if self.failures is None:
            self.failures = []


class _MiniprotDbProxy:
    def __init__(self, transcript, children):
        self._transcript = transcript
        self._children = children

    def __getitem__(self, feature_id):
        logical_id = self._transcript.attributes.get(
            "ID", [self._transcript.id],
        )[0]
        if str(feature_id) not in {str(self._transcript.id), str(logical_id)}:
            raise KeyError(feature_id)
        return self._transcript

    def children(self, parent, featuretype=None, order_by=None, **_kwargs):
        if str(getattr(parent, "id", parent)) != str(self._transcript.id):
            return iter(())
        if featuretype is None:
            wanted = None
        elif isinstance(featuretype, str):
            wanted = {featuretype}
        else:
            wanted = set(featuretype)
        rows = [
            child for child in self._children
            if wanted is None or child.featuretype in wanted
        ]
        if order_by == "start":
            rows.sort(key=lambda row: (
                int(row.start), int(row.end), str(row.id),
            ))
        return iter(rows)


class _ReferenceDbProxy:
    def __init__(self, features):
        self._features = features

    def __getitem__(self, feature_id):
        try:
            return self._features[str(feature_id)]
        except KeyError as exc:
            raise KeyError(feature_id) from exc


class _AnnotationProxy:
    def __init__(self, features):
        self.db_connection = _ReferenceDbProxy(features)


def materialize_miniprot_payload(index, transcript, ref_db, m_feature_db,
                                 m_id_2_ref_id_trans_dict,
                                 ref_features_reverse_dict):
    """Copy the records a Step-8 worker needs on the DB-owning thread."""
    try:
        transcript_copy = copy.deepcopy(transcript)
    except Exception as exc:
        return MiniprotPayload(
            index=index,
            transcript=transcript,
            children=(),
            reference_features={},
            error=exc,
            error_tb=traceback.format_exc(),
        )
    try:
        children = tuple(
            copy.deepcopy(child) for child in m_feature_db.children(
                transcript,
                featuretype=("CDS", "stop_codon"),
                order_by="start",
            )
        )
        ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(
            ref_features_reverse_dict,
            transcript.attributes["ID"][0],
            m_id_2_ref_id_trans_dict,
        )
        reference_features = {}
        for feature_id in (ref_gene_id, ref_trans_id):
            if feature_id is None:
                continue
            try:
                reference_features[feature_id] = copy.deepcopy(
                    ref_db.db_connection[feature_id]
                )
            except (KeyError, gffutils.exceptions.FeatureNotFoundError):
                continue
        return MiniprotPayload(
            index=index,
            transcript=transcript_copy,
            children=children,
            reference_features=reference_features,
        )
    except Exception as exc:
        return MiniprotPayload(
            index=index,
            transcript=transcript_copy,
            children=(),
            reference_features={},
            error=exc,
            error_tb=traceback.format_exc(),
        )


def analyze_miniprot_candidate(payload, *, tgt_fai, ref_proteins, ref_trans,
                               ref_features_dict,
                               m_id_2_ref_id_trans_dict,
                               ref_features_len_dict,
                               ref_trans_exon_num_dict,
                               ref_features_reverse_dict, args):
    """Run expensive candidate work against copied records and private state."""
    if payload.error is not None:
        return MiniprotAnalysis(
            index=payload.index,
            transcript=payload.transcript,
            error=payload.error,
            error_tb=payload.error_tb,
        )
    transcript_id = payload.transcript.attributes["ID"][0]
    ref_gene_id, _ = lifton_utils.get_ref_ids_miniprot(
        ref_features_reverse_dict,
        transcript_id,
        m_id_2_ref_id_trans_dict,
    )
    isolated_features = {}
    if ref_gene_id in ref_features_dict:
        clone = copy.copy(ref_features_dict[ref_gene_id])
        clone.copy_num = 0
        isolated_features[ref_gene_id] = clone
    explicit_copy_number = None
    raw_explicit = payload.transcript.attributes.get("extra_copy_number")
    if raw_explicit:
        try:
            explicit_copy_number = int(raw_explicit[0])
        except (TypeError, ValueError, IndexError):
            explicit_copy_number = None
    journal = DeferredStateJournal(isolated_features, buffer_score=True)
    try:
        gene = run_miniprot.process_miniprot(
            payload.transcript,
            _AnnotationProxy(payload.reference_features),
            _MiniprotDbProxy(payload.transcript, payload.children),
            {},
            tgt_fai,
            ref_proteins,
            ref_trans,
            isolated_features,
            journal.score_handle,
            m_id_2_ref_id_trans_dict,
            ref_features_len_dict,
            ref_trans_exon_num_dict,
            ref_features_reverse_dict,
            args,
            state_journal=journal,
        )
        delta = journal.finish()
        return MiniprotAnalysis(
            index=payload.index,
            transcript=payload.transcript,
            gene=gene,
            score_text=delta.score_text,
            explicit_copy_number=explicit_copy_number,
        )
    except Exception as exc:
        journal.finish()
        return MiniprotAnalysis(
            index=payload.index,
            transcript=payload.transcript,
            error=exc,
            error_tb=traceback.format_exc(),
        )


def _rebase_gene_copy(gene, copy_num):
    """Change a privately analyzed candidate to its ordered copy identity."""
    old_gene_id = str(gene.entry.id)
    new_gene_id = (
        f"{gene.ref_gene_id}_{copy_num}" if copy_num > 0 else gene.ref_gene_id
    )
    gene.copy_num = copy_num
    gene.entry.id = new_gene_id
    gene.entry.attributes["ID"] = [new_gene_id]
    if copy_num > 0:
        gene.entry.attributes["extra_copy_number"] = [str(copy_num)]
    else:
        gene.entry.attributes.pop("extra_copy_number", None)

    transcript_ids = {}
    rebased = {}
    for old_key, transcript in gene.transcripts.items():
        if not hasattr(transcript, "ref_tran_id"):
            rebased[old_key] = transcript
            continue
        old_transcript_id = str(transcript.entry.id)
        new_transcript_id = (
            f"{transcript.ref_tran_id}_{copy_num}"
            if copy_num > 0 else transcript.ref_tran_id
        )
        transcript.entry.id = new_transcript_id
        transcript.entry.attributes["ID"] = [new_transcript_id]
        if "Parent" in transcript.entry.attributes:
            transcript.entry.attributes["Parent"] = [new_gene_id]
        if "transcript_id" in transcript.entry.attributes:
            transcript.entry.attributes["transcript_id"] = [new_transcript_id]
        for exon in transcript.exons:
            exon.entry.attributes["Parent"] = [new_transcript_id]
            if exon.cds is not None:
                exon.cds.entry.attributes["Parent"] = [new_transcript_id]
        transcript_ids[old_transcript_id] = new_transcript_id
        rebased[new_transcript_id] = transcript
    gene.transcripts = rebased
    return old_gene_id, new_gene_id, transcript_ids


def _rebase_score_text(score_text, transcript_ids):
    lines = []
    for line in score_text.splitlines(keepends=True):
        body = line[:-1] if line.endswith("\n") else line
        newline = "\n" if line.endswith("\n") else ""
        feature_id, separator, remainder = body.partition("\t")
        feature_id = transcript_ids.get(feature_id, feature_id)
        lines.append(feature_id + separator + remainder + newline)
    return "".join(lines)


def parallel_step8(transcripts: Iterable[Any], ref_db, m_feature_db, tree_dict,
                   tgt_fai, ref_proteins, ref_trans, ref_features_dict,
                   fw, fw_score, transcripts_stats_dict,
                   m_id_2_ref_id_trans_dict, ref_features_len_dict,
                   ref_trans_exon_num_dict, ref_features_reverse_dict, args,
                   *, threads=1, max_inflight=None):
    """Analyze Step-8 candidates concurrently and finalize them in order."""
    threads = max(1, int(threads or 1))
    limit = max(1, int(max_inflight or (2 * threads)))
    outcome = Step8Outcome()

    def payloads() -> Iterator[MiniprotPayload]:
        for index, transcript in enumerate(transcripts):
            yield materialize_miniprot_payload(
                index, transcript, ref_db, m_feature_db,
                m_id_2_ref_id_trans_dict, ref_features_reverse_dict,
            )

    def worker(payload):
        return analyze_miniprot_candidate(
            payload,
            tgt_fai=tgt_fai,
            ref_proteins=ref_proteins,
            ref_trans=ref_trans,
            ref_features_dict=ref_features_dict,
            m_id_2_ref_id_trans_dict=m_id_2_ref_id_trans_dict,
            ref_features_len_dict=ref_features_len_dict,
            ref_trans_exon_num_dict=ref_trans_exon_num_dict,
            ref_features_reverse_dict=ref_features_reverse_dict,
            args=args,
        )

    high_water = [0]
    if threads <= 1:
        analyses = map(worker, payloads())
        executor = None
        high_water[0] = 1
    else:
        from lifton.parallel import _ordered_bounded_map

        executor = ThreadPoolExecutor(
            max_workers=threads, thread_name_prefix="lifton-step8",
        )
        analyses = _ordered_bounded_map(
            executor,
            payloads(),
            lambda pool, payload: pool.submit(worker, payload),
            limit,
            high_water_callback=lambda value: high_water.__setitem__(
                0, max(high_water[0], value)
            ),
        )

    try:
        for analysis in analyses:
            outcome.processed += 1
            transcript_id = str(getattr(analysis.transcript, "id", "<unknown>"))
            interval = Interval(
                analysis.transcript.start,
                analysis.transcript.end,
                transcript_id,
            )
            if lifton_utils.check_ovps_ratio(
                    analysis.transcript, interval, args.overlap, tree_dict):
                continue
            if analysis.error is not None:
                outcome.failures.append({
                    "mRNA": transcript_id,
                    "kind": "processing_error",
                    "message": safe_exception_text(
                        analysis.error, analysis.error_tb,
                    ),
                    "traceback": analysis.error_tb,
                })
                continue
            gene = analysis.gene
            if gene is None or gene.ref_gene_id is None:
                continue

            journal = DeferredStateJournal(ref_features_dict)
            allocation_entry = copy.copy(gene.entry)
            allocation_entry.attributes = copy.deepcopy(gene.entry.attributes)
            if analysis.explicit_copy_number is None:
                allocation_entry.attributes.pop("extra_copy_number", None)
            else:
                allocation_entry.attributes["extra_copy_number"] = [
                    str(analysis.explicit_copy_number)
                ]
            copy_num = journal.allocate_gene_copy(
                gene.ref_gene_id, allocation_entry, ref_features_dict,
            )
            _, _, transcript_ids = _rebase_gene_copy(gene, copy_num)
            block = io.StringIO()
            cds_id_allocator = getattr(
                args, "_cds_id_allocator", None
            )
            write_result = (
                gene.write_entry(block, transcripts_stats_dict)
                if cds_id_allocator is None
                else gene.write_entry(
                    block,
                    transcripts_stats_dict,
                    cds_id_allocator=cds_id_allocator,
                )
            )
            if write_result is False:
                outcome.failures.append({
                    "mRNA": transcript_id,
                    "kind": "serialization_error",
                    "message": "miniprot hierarchy was not serializable",
                })
                continue
            fw.write(block.getvalue())
            delta = journal.finish()
            commit_locus_delta(delta, ref_features_dict, tree_dict)
            fw_score.write(_rebase_score_text(
                analysis.score_text, transcript_ids,
            ))
            outcome.emitted += 1
            outcome.emitted_ref_gene_ids.add(gene.ref_gene_id)
            for feature_id, detail in getattr(
                    gene, "_serialization_failures", []):
                outcome.failures.append({
                    "mRNA": transcript_id,
                    "feature_id": feature_id,
                    "kind": "serialization_error",
                    "message": detail,
                })
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    outcome.max_inflight_observed = high_water[0]
    return outcome
