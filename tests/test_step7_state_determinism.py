"""Regression tests for Step 7's ordered state/output commit boundary."""

from __future__ import annotations

import copy
import io
import random
import threading
import time
from types import SimpleNamespace

import gffutils
import pytest

from lifton import lifton_class, parallel
from lifton.locus_pipeline import (
    DeferredStateJournal,
    LocusDelta,
    LocusResult,
    StepContext,
    Step7StateCoordinator,
    commit_locus_delta,
    consume,
)


def _feature(feature_id, start, end, featuretype="gene", parent=None):
    attrs = {"ID": [feature_id]}
    if parent is not None:
        attrs["Parent"] = [parent]
    return gffutils.Feature(
        seqid="chr1", source="Liftoff", featuretype=featuretype,
        start=start, end=end, score=".", strand="+", frame=".",
        attributes=attrs,
    )


class _FeatureDB:
    def __init__(self, loci):
        self._loci = loci

    def features_of_type(self, _featuretype):
        return iter(self._loci)

    def children(self, _feature, **_kwargs):
        return iter(())


class _ParentOnlyIO(io.StringIO):
    """Record which threads touched a real output handle."""

    def __init__(self):
        super().__init__()
        self.writer_threads = set()

    def write(self, text):
        self.writer_threads.add(threading.current_thread().name)
        return super().write(text)


def _run_copy_rich(monkeypatch, threads, seed):
    from lifton import run_liftoff

    loci = [
        _feature("refA" if index == 0 else f"refA_{index}",
                 100 + index * 20, 110 + index * 20)
        for index in range(12)
    ]
    locus_index = {locus.id: index for index, locus in enumerate(loci)}
    delays = list(range(len(loci)))
    random.Random(seed).shuffle(delays)

    ref_feature = lifton_class.Lifton_feature("refA")
    ref_features = {"refA": ref_feature}
    ref_attrs = {
        "ID": ["refA"],
        "gene_biotype": ["protein_coding"],
    }
    ref_db = {"refA": SimpleNamespace(attributes=ref_attrs)}
    tree_dict = {}
    gff = _ParentOnlyIO()
    score = _ParentOnlyIO()
    chain = _ParentOnlyIO()
    stats = {"coding": {}, "non-coding": {}, "other": {}}
    args = SimpleNamespace(
        annotation_database="RefSeq", native=False, write_chains=True,
        writer_max_pending=0,
    )
    ctx = StepContext(
        ref_db=ref_db,
        l_feature_db=_FeatureDB(loci),
        m_feature_db=None,
        ref_id_2_m_id_trans_dict={},
        tree_dict=tree_dict,
        tgt_fai=None,
        ref_proteins={},
        ref_trans={},
        ref_features_dict=ref_features,
        fw_score=score,
        fw_chain=chain,
        args=args,
    )

    def fake_process(*call_args, ENTRY_FEATURE=False,
                     state_journal=None, **_kwargs):
        assert ENTRY_FEATURE is True
        locus = call_args[1]
        index = locus_index[locus.id]
        # Delay before copy allocation so workers reach the constructor in a
        # deliberately non-submission order.
        time.sleep(delays[index] * 0.0005)
        gene = lifton_class.Lifton_GENE(
            "refA", copy.deepcopy(locus), copy.deepcopy(ref_attrs),
            call_args[6], call_args[10], call_args[13],
            state_journal=state_journal,
        )
        transcript = _feature(
            f"lifted_tx_{index}", locus.start, locus.end,
            featuretype="mRNA", parent=locus.id,
        )
        gene.add_transcript(
            "txA", transcript,
            {"ID": ["txA"], "Parent": ["refA"]},
        )
        call_args[11].write(f"score\t{gene.entry.id}\n")
        call_args[12].write(f"chain\t{gene.entry.id}\n")
        return gene

    monkeypatch.setattr(run_liftoff, "process_liftoff", fake_process)
    parallel.parallel_step7(
        ["gene"], ctx.l_feature_db, ctx, gff, stats,
        threads=threads, progress_every=0,
    )

    expected_ids = {"refA"} | {f"refA_{index}" for index in range(1, 12)}
    assert ref_feature.copy_num == 12
    assert {interval.data for interval in tree_dict["chr1"]} == expected_ids
    # The real handles belong exclusively to ordered parent-thread consume.
    assert gff.writer_threads == {"MainThread"}
    assert score.writer_threads == {"MainThread"}
    assert chain.writer_threads == {"MainThread"}
    return gff.getvalue(), score.getvalue(), chain.getvalue(), stats


@pytest.mark.parametrize("seed", [3, 29, 101])
def test_copy_score_chain_and_stats_are_deterministic(monkeypatch, seed):
    baseline = _run_copy_rich(monkeypatch, 1, seed)
    for threads in (2, 4, 8):
        assert _run_copy_rich(monkeypatch, threads, seed) == baseline
    assert baseline[3] == {
        "coding": {"txA": 12}, "non-coding": {}, "other": {},
    }


def test_deferred_journal_only_consumes_copy_on_commit():
    ref_feature = lifton_class.Lifton_feature("refA")
    ref_features = {"refA": ref_feature}
    args = SimpleNamespace(annotation_database="RefSeq")
    attrs = {"ID": ["refA"], "gene_biotype": ["protein_coding"]}

    rejected = DeferredStateJournal(ref_features)
    rejected_gene = lifton_class.Lifton_GENE(
        "refA", _feature("candidate_1", 1, 10), copy.deepcopy(attrs),
        {}, ref_features, args, state_journal=rejected,
    )
    assert rejected_gene.entry.id == "refA"
    assert ref_feature.copy_num == 0

    accepted = DeferredStateJournal(ref_features)
    accepted_gene = lifton_class.Lifton_GENE(
        "refA", _feature("candidate_2", 20, 30), copy.deepcopy(attrs),
        {}, ref_features, args, state_journal=accepted,
    )
    tree_dict = {}
    commit_locus_delta(accepted.finish(), ref_features, tree_dict)
    assert accepted_gene.entry.id == "refA"
    assert ref_feature.copy_num == 1
    assert [interval.data for interval in tree_dict["chr1"]] == ["refA"]


def test_consume_records_processing_and_serialization_failures():
    stats = {"coding": {}, "non-coding": {}, "other": {}}
    ctx = SimpleNamespace(
        failure_records=[], fw_score=None, fw_chain=None,
    )
    processing = LocusResult(
        index=0, locus_id="bad_process", error=RuntimeError("boom"),
        error_tb="synthetic traceback",
    )
    assert consume(processing, io.StringIO(), stats, ctx=ctx) is False

    class RejectedGene:
        ref_gene_id = "refA"
        _serialization_failures = [("tx_bad", "invalid coordinates")]

        @staticmethod
        def write_entry(_fw, _stats):
            return False

    rejected = LocusResult(
        index=1, locus_id="bad_write", lifton_gene=RejectedGene(),
    )
    assert consume(rejected, io.StringIO(), stats, ctx=ctx) is False
    assert [record["kind"] for record in ctx.failure_records] == [
        "processing_error", "serialization_error",
    ]


def test_failed_locus_does_not_commit_state_or_sidecars():
    ref_feature = lifton_class.Lifton_feature("refA")
    ref_features = {"refA": ref_feature}
    tree_dict = {}
    score = io.StringIO()
    chain = io.StringIO()
    ctx = SimpleNamespace(
        failure_records=[], fw_score=score, fw_chain=chain,
        emitted_ref_gene_ids=set(), emitted_intervals=[],
    )
    delta = LocusDelta(
        copy_increments={"refA": 1},
        tree_intervals=[("chr1", 1, 10, "refA")],
        score_text="orphan score\n",
        chain_text="orphan chain\n",
    )
    result = LocusResult(
        index=0, locus_id="refA", error=RuntimeError("failed"),
        error_tb="synthetic traceback", delta=delta,
    )
    coordinator = Step7StateCoordinator(ref_features, tree_dict)

    assert consume(
        result, io.StringIO(),
        {"coding": {}, "non-coding": {}, "other": {}},
        ctx=ctx, state_coordinator=coordinator,
    ) is False
    assert ref_feature.copy_num == 0
    assert tree_dict == {}
    assert score.getvalue() == ""
    assert chain.getvalue() == ""
    assert ctx.emitted_ref_gene_ids == set()
    assert ctx.emitted_intervals == []
