"""Focused contracts for read-only annotation evaluation."""

from io import StringIO
from types import SimpleNamespace

import gffutils
import pytest

from lifton import lifton_class, run_evaluation


def _db(rows):
    text = "##gff-version 3\n" + "\n".join(rows) + "\n"
    return gffutils.create_db(
        text,
        dbfn=":memory:",
        from_string=True,
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )


def _row(featuretype, start, end, attributes, source="test", frame="."):
    return (
        f"chr1\t{source}\t{featuretype}\t{start}\t{end}\t.\t+\t{frame}\t"
        f"{attributes}"
    )


def _args(chm13=False):
    return SimpleNamespace(
        annotation_database="RefSeq",
        evaluation_liftoff_chm13=chm13,
        debug=False,
    )


def _ref_features(gene_id):
    return {gene_id: lifton_class.Lifton_feature(gene_id)}


def _capture_orf(monkeypatch):
    calls = []

    def fake_orf(self, trans_id, ref_trans_id, _fai, _proteins, _transcripts,
                 status, eval_only=False, eval_liftoff_chm13=False):
        calls.append({
            "trans_id": trans_id,
            "ref_trans_id": ref_trans_id,
            "eval_only": eval_only,
            "chm13": eval_liftoff_chm13,
            "exons": len(self.transcripts[trans_id].exons),
        })
        status.lifton_dna = 0.75
        status.lifton_aa = 0.5
        status.status = ["nonsynonymous"]
        return None, None

    monkeypatch.setattr(lifton_class.Lifton_GENE, "orf_search_protein", fake_orf)
    return calls


def _evaluate(ref_db, tgt_db, ref_features, args):
    score = StringIO()
    genes = list(tgt_db.features_of_type("gene"))
    assert len(genes) == 1
    result = run_evaluation.evaluation(
        None,
        genes[0],
        ref_db,
        tgt_db,
        {},
        object(),
        {},
        {},
        ref_features,
        score,
        args,
        ENTRY_FEATURE=True,
    )
    return result, score.getvalue()


@pytest.mark.parametrize(
    ("chm13", "ref_gene", "ref_trans", "target_gene", "target_trans"),
    [
        (False, "gene1", "tx1", "gene1", "tx1"),
        (True, "gene-G1", "rna-T1", "G1", "T1"),
    ],
)
def test_normal_and_chm13_ids_use_their_own_reference_namespaces(
        monkeypatch, chm13, ref_gene, ref_trans, target_gene, target_trans):
    ref_db = _db([
        _row("gene", 1, 90, f"ID={ref_gene};gene_biotype=protein_coding"),
        _row("mRNA", 1, 90, f"ID={ref_trans};Parent={ref_gene}"),
    ])
    tgt_db = _db([
        _row("gene", 1, 90, f"ID={target_gene}"),
        _row("mRNA", 1, 90, f"ID={target_trans};Parent={target_gene}"),
        _row("exon", 1, 90, f"ID=ex1;Parent={target_trans}"),
        _row("CDS", 1, 90, f"ID=cds1;Parent={target_trans}", frame="0"),
    ])
    feature_key = "G1" if chm13 else "gene1"
    calls = _capture_orf(monkeypatch)

    _, score = _evaluate(ref_db, tgt_db, _ref_features(feature_key), _args(chm13))

    logical_transcript = "T1" if chm13 else "tx1"
    assert calls == [{
        "trans_id": logical_transcript,
        "ref_trans_id": logical_transcript,
        "eval_only": True,
        "chm13": chm13,
        "exons": 1,
    }]
    assert score == (
        f"{logical_transcript}\t0.75\t0.5\tnonsynonymous\tchr1:1-90\n"
    )


def test_copy_suffix_falls_back_to_reference_base(monkeypatch):
    ref_db = _db([
        _row("gene", 1, 90, "ID=gene1;gene_biotype=protein_coding"),
        _row("mRNA", 1, 90, "ID=tx1;Parent=gene1"),
    ])
    tgt_db = _db([
        _row("gene", 1, 90, "ID=gene1_1;extra_copy_number=1"),
        _row("mRNA", 1, 90, "ID=tx1_1;Parent=gene1_1"),
        _row("exon", 1, 90, "ID=ex1_1;Parent=tx1_1"),
        _row("CDS", 1, 90, "ID=cds1_1;Parent=tx1_1", frame="0"),
    ])
    calls = _capture_orf(monkeypatch)

    _, score = _evaluate(ref_db, tgt_db, _ref_features("gene1"), _args())

    assert calls[0]["ref_trans_id"] == "tx1"
    assert calls[0]["trans_id"] == "tx1_1"
    assert score.startswith("tx1_1\t")


def test_exact_numeric_suffix_wins_over_copy_fallback(monkeypatch):
    ref_db = _db([
        _row("gene", 1, 90, "ID=gene1;gene_biotype=protein_coding"),
        _row("mRNA", 1, 90, "ID=tx;Parent=gene1"),
        _row("mRNA", 1, 90, "ID=tx_1;Parent=gene1"),
    ])
    tgt_db = _db([
        _row("gene", 1, 90, "ID=gene1"),
        _row("mRNA", 1, 90, "ID=tx_1;Parent=gene1"),
        _row("exon", 1, 90, "ID=ex1;Parent=tx_1"),
        _row("CDS", 1, 90, "ID=cds1;Parent=tx_1", frame="0"),
    ])
    calls = _capture_orf(monkeypatch)

    _evaluate(ref_db, tgt_db, _ref_features("gene1"), _args())

    assert calls[0]["ref_trans_id"] == "tx_1"
    assert calls[0]["trans_id"] == "tx_1"


def test_missing_reference_ids_skip_only_unmapped_models(monkeypatch):
    ref_db = _db([
        _row("gene", 1, 200, "ID=gene1;gene_biotype=protein_coding"),
        _row("mRNA", 101, 200, "ID=valid;Parent=gene1"),
    ])
    tgt_db = _db([
        _row("gene", 1, 200, "ID=gene1"),
        _row("mRNA", 1, 90, "ID=missing;Parent=gene1"),
        _row("exon", 1, 90, "ID=missing-ex;Parent=missing"),
        _row("mRNA", 101, 200, "ID=valid;Parent=gene1"),
        _row("exon", 101, 200, "ID=valid-ex;Parent=valid"),
    ])
    calls = _capture_orf(monkeypatch)

    _, score = _evaluate(ref_db, tgt_db, _ref_features("gene1"), _args())

    assert [call["ref_trans_id"] for call in calls] == ["valid"]
    assert score.startswith("valid\t")
    assert "missing" not in score

    absent_gene_db = _db([
        _row("gene", 1, 90, "ID=other;gene_biotype=protein_coding"),
    ])
    result, absent_score = _evaluate(
        absent_gene_db,
        _db([_row("gene", 1, 90, "ID=gene1")]),
        _ref_features("gene1"),
        _args(),
    )
    assert result is None
    assert absent_score == ""


def test_evaluation_does_not_mutate_copy_counter_or_interval_tree():
    ref_db = _db([
        _row("gene", 1, 90, "ID=gene1;gene_biotype=protein_coding"),
    ])
    tgt_db = _db([_row("gene", 1, 90, "ID=gene1")])
    ref_features = _ref_features("gene1")
    tree = {"sentinel": object()}

    gene, ref_gene_id, ref_trans_id = run_evaluation.initialize_lifton_gene_eval(
        next(tgt_db.features_of_type("gene")),
        ref_db,
        tree,
        ref_features,
        _args(),
    )

    assert gene is not None
    assert (ref_gene_id, ref_trans_id) == ("gene1", None)
    assert ref_features["gene1"].copy_num == 0
    assert list(tree) == ["sentinel"]


def test_cds_only_and_nested_target_hierarchies_are_evaluated(monkeypatch):
    ref_db = _db([
        _row("gene", 1, 90, "ID=gene1;gene_biotype=protein_coding"),
        _row("mRNA", 1, 90, "ID=tx1;Parent=gene1"),
    ])
    tgt_db = _db([
        _row("gene", 1, 90, "ID=gene1"),
        _row("primary_transcript", 1, 90, "ID=container;Parent=gene1"),
        _row("mRNA", 1, 90, "ID=tx1;Parent=container"),
        _row("CDS", 1, 90, "ID=cds1;Parent=tx1", frame="0"),
    ])
    calls = _capture_orf(monkeypatch)

    _, score = _evaluate(ref_db, tgt_db, _ref_features("gene1"), _args())

    assert calls[0]["ref_trans_id"] == "tx1"
    assert calls[0]["exons"] == 1
    assert score.startswith("tx1\t")


def test_transcript_score_order_is_coordinate_deterministic(monkeypatch):
    ref_db = _db([
        _row("gene", 1, 300, "ID=gene1;gene_biotype=protein_coding"),
        _row("mRNA", 200, 300, "ID=txB;Parent=gene1"),
        _row("mRNA", 1, 100, "ID=txA;Parent=gene1"),
    ])
    tgt_db = _db([
        _row("gene", 1, 300, "ID=gene1"),
        _row("mRNA", 200, 300, "ID=txB;Parent=gene1"),
        _row("exon", 200, 300, "ID=exB;Parent=txB"),
        _row("mRNA", 1, 100, "ID=txA;Parent=gene1"),
        _row("exon", 1, 100, "ID=exA;Parent=txA"),
    ])
    _capture_orf(monkeypatch)

    _, score = _evaluate(ref_db, tgt_db, _ref_features("gene1"), _args())

    assert [line.split("\t", 1)[0] for line in score.splitlines()] == [
        "txA", "txB",
    ]


def test_unexpected_reference_lookup_interrupt_propagates():
    class InterruptedDB:
        def __getitem__(self, _key):
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        run_evaluation._lookup_reference_feature(InterruptedDB(), "tx1")
