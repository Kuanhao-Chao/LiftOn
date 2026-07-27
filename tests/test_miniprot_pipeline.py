from types import SimpleNamespace
from io import StringIO
import time

import gffutils
from intervaltree import Interval, IntervalTree

from lifton import lifton_class, miniprot_pipeline
from lifton import run_miniprot
from lifton.gffbase import create_db as create_gffbase_db
from lifton.locus_pipeline import DeferredStateJournal, commit_locus_delta


def _transcript(feature_id, start):
    return SimpleNamespace(
        id=feature_id,
        seqid="chr1",
        start=start,
        end=start + 20,
        attributes={"ID": [feature_id]},
    )


class _Gene:
    def __init__(self, gene_id, transcript_id, *, writable=True):
        self.ref_gene_id = gene_id
        self.copy_num = 0
        self.entry = SimpleNamespace(
            id=gene_id,
            seqid="chr1",
            start=1,
            end=20,
            featuretype="gene",
            attributes={"ID": [gene_id]},
        )
        trans_entry = SimpleNamespace(
            id=transcript_id,
            attributes={
                "ID": [transcript_id],
                "Parent": [gene_id],
                "transcript_id": [transcript_id],
            },
        )
        self.transcripts = {
            transcript_id: SimpleNamespace(
                ref_tran_id=transcript_id,
                entry=trans_entry,
                exons=[],
            )
        }
        self._writable = writable
        self._serialization_failures = []

    def write_entry(self, handle, _stats):
        if not self._writable:
            return False
        transcript_id = next(iter(self.transcripts))
        handle.write(f"{self.entry.id}\t{transcript_id}\n")
        return True


def _run(monkeypatch, analyses, *, threads=2, max_inflight=2, tree=None):
    transcripts = [analysis.transcript for analysis in analyses]

    def fake_materialize(index, transcript, *_args):
        return miniprot_pipeline.MiniprotPayload(index, transcript, (), {})

    by_id = {analysis.transcript.id: analysis for analysis in analyses}

    def fake_analyze(payload, **_kwargs):
        if payload.index == 0 and threads > 1:
            time.sleep(0.03)
        return by_id[payload.transcript.id]

    monkeypatch.setattr(
        miniprot_pipeline, "materialize_miniprot_payload", fake_materialize,
    )
    monkeypatch.setattr(
        miniprot_pipeline, "analyze_miniprot_candidate", fake_analyze,
    )
    output = StringIO()
    score = StringIO()
    ref_features = {"gene": lifton_class.Lifton_feature("gene")}
    outcome = miniprot_pipeline.parallel_step8(
        transcripts,
        object(),
        object(),
        {} if tree is None else tree,
        object(),
        {},
        {},
        ref_features,
        output,
        score,
        {"coding": {}, "non-coding": {}, "other": {}},
        {},
        {},
        {},
        {},
        SimpleNamespace(overlap=0.1),
        threads=threads,
        max_inflight=max_inflight,
    )
    return outcome, output.getvalue(), score.getvalue(), ref_features


def test_step8_parallel_analysis_publishes_in_order_and_rebases_copy_ids(monkeypatch):
    first = _transcript("mp1", 1)
    second = _transcript("mp2", 101)
    analyses = [
        miniprot_pipeline.MiniprotAnalysis(
            0, first, _Gene("gene", "tx"), "tx\tfirst\n",
        ),
        miniprot_pipeline.MiniprotAnalysis(
            1, second, _Gene("gene", "tx"), "tx\tsecond\n",
        ),
    ]

    outcome, output, score, ref_features = _run(monkeypatch, analyses)

    assert output == "gene\ttx\ngene_1\ttx_1\n"
    assert score == "tx\tfirst\ntx_1\tsecond\n"
    assert outcome.processed == 2
    assert outcome.emitted == 2
    assert outcome.max_inflight_observed <= 2
    assert ref_features["gene"].copy_num == 2


def test_failed_serialization_does_not_consume_copy_number(monkeypatch):
    first = _transcript("mp1", 1)
    second = _transcript("mp2", 101)
    analyses = [
        miniprot_pipeline.MiniprotAnalysis(
            0, first, _Gene("gene", "tx", writable=False), "tx\tbad\n",
        ),
        miniprot_pipeline.MiniprotAnalysis(
            1, second, _Gene("gene", "tx"), "tx\tgood\n",
        ),
    ]

    outcome, output, score, ref_features = _run(monkeypatch, analyses)

    assert output == "gene\ttx\n"
    assert score == "tx\tgood\n"
    assert outcome.emitted == 1
    assert outcome.failures[0]["kind"] == "serialization_error"
    assert ref_features["gene"].copy_num == 1


def test_worker_failure_is_ordered_and_does_not_stop_following_candidates(monkeypatch):
    first = _transcript("mp1", 1)
    second = _transcript("mp2", 101)
    analyses = [
        miniprot_pipeline.MiniprotAnalysis(
            0, first, error=ValueError("broken"), error_tb="trace",
        ),
        miniprot_pipeline.MiniprotAnalysis(
            1, second, _Gene("gene", "tx"), "tx\tgood\n",
        ),
    ]

    outcome, output, score, _ = _run(monkeypatch, analyses)

    assert output == "gene\ttx\n"
    assert score == "tx\tgood\n"
    assert outcome.failures[0]["mRNA"] == "mp1"
    assert outcome.failures[0]["kind"] == "processing_error"


def test_worker_failure_with_broken_str_is_still_reported(monkeypatch):
    class BrokenStrError(RuntimeError):
        def __str__(self):
            raise TypeError("broken exception renderer")

    transcript = _transcript("mp1", 1)
    analyses = [miniprot_pipeline.MiniprotAnalysis(
        0, transcript, error=BrokenStrError("broken"), error_tb="trace",
    )]

    outcome, output, score, _ = _run(monkeypatch, analyses)

    assert output == ""
    assert score == ""
    assert outcome.failures[0]["message"] == "trace"


def test_overlapped_worker_failure_matches_serial_skip(monkeypatch):
    transcript = _transcript("redundant", 1)
    analyses = [miniprot_pipeline.MiniprotAnalysis(
        0, transcript, error=ValueError("malformed child"), error_tb="trace",
    )]
    tree = {"chr1": IntervalTree([Interval(1, 21, "existing")])}

    outcome, output, score, _ = _run(
        monkeypatch, analyses, tree=tree,
    )

    assert outcome.processed == 1
    assert outcome.emitted == 0
    assert outcome.failures == []
    assert outcome.preoverlap_elided == 1
    assert outcome.analysis_submitted == 0
    assert output == ""
    assert score == ""


def test_static_overlap_elides_materialization_and_analysis(monkeypatch):
    transcript = _transcript("redundant", 1)
    tree = {"chr1": IntervalTree([Interval(1, 21, "existing")])}

    def forbidden(*_args, **_kwargs):
        raise AssertionError("static overlap reached candidate work")

    monkeypatch.setattr(
        miniprot_pipeline, "materialize_miniprot_payload", forbidden,
    )
    monkeypatch.setattr(
        miniprot_pipeline, "analyze_miniprot_candidate", forbidden,
    )
    outcome = miniprot_pipeline.parallel_step8(
        [transcript],
        object(),
        object(),
        tree,
        object(),
        {},
        {},
        {},
        StringIO(),
        StringIO(),
        {"coding": {}, "non-coding": {}, "other": {}},
        {},
        {},
        {},
        {},
        SimpleNamespace(overlap=0.1),
        threads=2,
        max_inflight=2,
    )

    assert outcome.processed == 1
    assert outcome.emitted == 0
    assert outcome.failures == []
    assert outcome.preoverlap_elided == 1
    assert outcome.analysis_submitted == 0
    assert outcome.processed == (
        outcome.preoverlap_elided + outcome.analysis_submitted
    )
    assert outcome.child_batch_calls == 0
    assert outcome.child_scalar_materializations == 0


def test_parent_rechecks_residual_after_earlier_step8_commit(monkeypatch):
    first = _transcript("mp1", 1)
    second = _transcript("mp2", 1)
    analyses = [
        miniprot_pipeline.MiniprotAnalysis(
            0, first, _Gene("gene", "tx"), "tx\tfirst\n",
        ),
        miniprot_pipeline.MiniprotAnalysis(
            1, second, _Gene("gene", "tx"), "tx\tsecond\n",
        ),
    ]

    outcome, output, score, ref_features = _run(monkeypatch, analyses)

    assert output == "gene\ttx\n"
    assert score == "tx\tfirst\n"
    assert outcome.processed == 2
    assert outcome.emitted == 1
    assert outcome.preoverlap_elided == 0
    assert outcome.analysis_submitted == 2
    assert outcome.processed == (
        outcome.preoverlap_elided + outcome.analysis_submitted
    )
    assert ref_features["gene"].copy_num == 1


def test_residual_children_are_batched_within_inflight_window(monkeypatch):
    transcripts = [
        _transcript(f"mp{index}", index * 100 + 1)
        for index in range(5)
    ]

    class BatchedChildren:
        def __init__(self):
            self.batch_calls = []
            self.scalar_calls = []

        def children_batched_features(
                self, anchors, *, featuretype=None, level=None,
                order_by=None):
            anchors = tuple(anchors)
            self.batch_calls.append(
                (anchors, featuretype, level, order_by)
            )
            return {
                anchor: (
                    SimpleNamespace(
                        id=f"{anchor}-cds",
                        featuretype="CDS",
                        start=index * 100 + 1,
                        end=index * 100 + 20,
                        attributes={"ID": [f"{anchor}-cds"]},
                    ),
                )
                for index, anchor in enumerate(anchors)
            }

        def children(self, transcript, **_kwargs):
            self.scalar_calls.append(transcript.id)
            return iter(())

    database = BatchedChildren()
    observed_children = {}

    def fake_analyze(payload, **_kwargs):
        observed_children[payload.index] = tuple(
            child.id for child in payload.children
        )
        return miniprot_pipeline.MiniprotAnalysis(
            payload.index, payload.transcript,
        )

    monkeypatch.setattr(
        miniprot_pipeline, "analyze_miniprot_candidate", fake_analyze,
    )
    tree = {
        "chr1": IntervalTree([
            Interval(101, 121, "existing-1"),
            Interval(301, 321, "existing-3"),
        ])
    }
    outcome = miniprot_pipeline.parallel_step8(
        transcripts,
        SimpleNamespace(db_connection={}),
        database,
        tree,
        object(),
        {},
        {},
        {},
        StringIO(),
        StringIO(),
        {"coding": {}, "non-coding": {}, "other": {}},
        {},
        {},
        {},
        {},
        SimpleNamespace(overlap=0.1),
        threads=2,
        max_inflight=2,
    )

    assert [
        call[0] for call in database.batch_calls
    ] == [("mp0",), ("mp2",), ("mp4",)]
    assert all(
        featuretype == ("CDS", "stop_codon")
        and level is None
        and order_by == "start"
        for _anchors, featuretype, level, order_by in database.batch_calls
    )
    assert database.scalar_calls == []
    assert observed_children == {
        index: (f"mp{index}-cds",) for index in (0, 2, 4)
    }
    assert outcome.processed == 5
    assert outcome.preoverlap_elided == 2
    assert outcome.analysis_submitted == 3
    assert outcome.processed == (
        outcome.preoverlap_elided + outcome.analysis_submitted
    )
    assert outcome.child_batch_calls == 3
    assert outcome.child_batch_fallbacks == 0
    assert outcome.child_scalar_materializations == 0
    assert outcome.max_inflight_observed <= 2


def test_batch_failure_falls_back_and_isolates_candidate_error():
    transcripts = [
        _transcript("good", 1),
        _transcript("bad", 101),
        _transcript("later", 201),
    ]

    class BrokenBatch:
        def __init__(self):
            self.batch_calls = 0
            self.scalar_calls = []

        def children_batched_features(self, *_args, **_kwargs):
            self.batch_calls += 1
            raise RuntimeError("batch unavailable")

        def children(self, transcript, **_kwargs):
            self.scalar_calls.append(transcript.id)
            if transcript.id == "bad":
                raise RuntimeError("malformed hierarchy")
            return iter(())

    database = BrokenBatch()
    outcome = miniprot_pipeline.parallel_step8(
        transcripts,
        SimpleNamespace(db_connection={}),
        database,
        {},
        object(),
        {},
        {},
        {},
        StringIO(),
        StringIO(),
        {"coding": {}, "non-coding": {}, "other": {}},
        {},
        {},
        {},
        {},
        SimpleNamespace(overlap=0.1),
        threads=2,
        max_inflight=3,
    )

    assert database.batch_calls == 1
    assert database.scalar_calls == ["good", "bad", "later"]
    assert outcome.processed == 3
    assert outcome.analysis_submitted == 3
    assert outcome.processed == (
        outcome.preoverlap_elided + outcome.analysis_submitted
    )
    assert outcome.child_batch_calls == 1
    assert outcome.child_batch_fallbacks == 1
    assert outcome.child_scalar_materializations == 3
    assert len(outcome.failures) == 1
    assert outcome.failures[0]["mRNA"] == "bad"
    assert "malformed hierarchy" in outcome.failures[0]["message"]


def test_materialization_failure_becomes_ordered_candidate_error():
    transcript = _transcript("broken", 1)

    class BrokenChildren:
        def children(self, *_args, **_kwargs):
            raise RuntimeError("malformed hierarchy")

    payload = miniprot_pipeline.materialize_miniprot_payload(
        0,
        transcript,
        SimpleNamespace(db_connection={}),
        BrokenChildren(),
        {},
        {},
    )

    assert isinstance(payload.error, RuntimeError)
    assert "malformed hierarchy" in str(payload.error)


def _db(rows):
    return gffutils.create_db(
        "##gff-version 3\n" + "\n".join(rows) + "\n",
        dbfn=":memory:",
        from_string=True,
        force=True,
        keep_order=True,
        merge_strategy="merge",
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )


def _gffbase_db(rows):
    return create_gffbase_db(
        "##gff-version 3\n" + "\n".join(rows) + "\n",
        ":memory:",
        from_string=True,
        force=True,
        build_rtree=False,
    )


def test_real_candidate_proxy_matches_serial_finalize(monkeypatch):
    ref_db = _db([
        "chr1\tref\tgene\t1\t90\t.\t+\t.\tID=gene;gene_biotype=protein_coding",
        "chr1\tref\tmRNA\t1\t90\t.\t+\t.\tID=tx;Parent=gene",
    ])
    miniprot_rows = [
        "chr1\tminiprot\tmRNA\t101\t190\t.\t+\t.\tID=mp1",
        "chr1\tminiprot\tCDS\t101\t190\t.\t+\t0\tID=cds1;Parent=mp1",
        "chr1\tminiprot\tmRNA\t301\t390\t.\t+\t.\tID=mp2",
        "chr1\tminiprot\tCDS\t301\t390\t.\t+\t0\tID=cds2;Parent=mp2",
    ]
    miniprot_db = _db(miniprot_rows)
    parallel_miniprot_db = _gffbase_db(miniprot_rows)
    args = SimpleNamespace(
        overlap=0.1,
        min_miniprot=0.5,
        max_miniprot=2.0,
        annotation_database="RefSeq",
        debug=False,
    )
    mapping = {"mp1": "tx", "mp2": "tx"}
    reverse = {"tx": "gene"}
    lengths = {"gene": 90}
    exon_counts = {"tx": 1}

    monkeypatch.setattr(
        run_miniprot.align,
        "lifton_parasail_align",
        lambda *_args, **_kwargs: SimpleNamespace(identity=0.875),
    )

    def fake_orf(self, *_args, **_kwargs):
        return None, True

    monkeypatch.setattr(lifton_class.Lifton_GENE, "orf_search_protein", fake_orf)

    def run_serial():
        ref_features = {"gene": lifton_class.Lifton_feature("gene")}
        tree = {}
        output = StringIO()
        score = StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        for transcript in miniprot_db.features_of_type("mRNA"):
            journal = DeferredStateJournal(ref_features, buffer_score=True)
            gene = run_miniprot.process_miniprot(
                transcript,
                SimpleNamespace(db_connection=ref_db),
                miniprot_db,
                tree,
                object(),
                {"tx": "M"},
                {"tx": "ATG"},
                ref_features,
                journal.score_handle,
                mapping,
                lengths,
                exon_counts,
                reverse,
                args,
                state_journal=journal,
            )
            assert gene.write_entry(output, stats) is True
            delta = journal.finish()
            commit_locus_delta(delta, ref_features, tree)
            score.write(delta.score_text)
        return output.getvalue(), score.getvalue(), stats

    serial_output, serial_score, serial_stats = run_serial()
    parallel_output = StringIO()
    parallel_score = StringIO()
    parallel_stats = {"coding": {}, "non-coding": {}, "other": {}}
    ref_features = {"gene": lifton_class.Lifton_feature("gene")}
    outcome = miniprot_pipeline.parallel_step8(
        parallel_miniprot_db.features_of_type("mRNA"),
        SimpleNamespace(db_connection=ref_db),
        parallel_miniprot_db,
        {},
        object(),
        {"tx": "M"},
        {"tx": "ATG"},
        ref_features,
        parallel_output,
        parallel_score,
        parallel_stats,
        mapping,
        lengths,
        exon_counts,
        reverse,
        args,
        threads=2,
        max_inflight=2,
    )

    assert outcome.failures == []
    assert outcome.emitted == 2
    assert outcome.processed == (
        outcome.preoverlap_elided + outcome.analysis_submitted
    )
    assert outcome.child_batch_calls == 1
    assert outcome.child_batch_fallbacks == 0
    assert outcome.child_scalar_materializations == 0
    assert parallel_output.getvalue() == serial_output
    assert parallel_score.getvalue() == serial_score
    assert parallel_stats == serial_stats
