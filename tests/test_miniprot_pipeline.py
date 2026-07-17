from types import SimpleNamespace
from io import StringIO
import time

import gffutils
from intervaltree import Interval, IntervalTree

from lifton import lifton_class, miniprot_pipeline
from lifton import run_miniprot
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
    assert output == ""
    assert score == ""


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


def test_real_candidate_proxy_matches_serial_finalize(monkeypatch):
    ref_db = _db([
        "chr1\tref\tgene\t1\t90\t.\t+\t.\tID=gene;gene_biotype=protein_coding",
        "chr1\tref\tmRNA\t1\t90\t.\t+\t.\tID=tx;Parent=gene",
    ])
    miniprot_db = _db([
        "chr1\tminiprot\tmRNA\t101\t190\t.\t+\t.\tID=mp1",
        "chr1\tminiprot\tCDS\t101\t190\t.\t+\t0\tID=cds1;Parent=mp1",
        "chr1\tminiprot\tmRNA\t301\t390\t.\t+\t.\tID=mp2",
        "chr1\tminiprot\tCDS\t301\t390\t.\t+\t0\tID=cds2;Parent=mp2",
    ])
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
        miniprot_db.features_of_type("mRNA"),
        SimpleNamespace(db_connection=ref_db),
        miniprot_db,
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
    assert parallel_output.getvalue() == serial_output
    assert parallel_score.getvalue() == serial_score
    assert parallel_stats == serial_stats
