"""Every class of dropped feature is counted and reported.

The Liftoff ``-copies`` bug dropped a transcript and all its exons and CDS
~4,400 times across the benchmark corpus, in every release up to v1.0.11, while
printing a warning every single time. A warning per event in a stream of
hundreds of thousands of lines is indistinguishable from noise; one number at
the end is not. ``locus_pipeline._record_childless_gene`` added that for one
class, and this generalises it.
"""
from __future__ import annotations

import textwrap
import threading
import types

import pytest

from lifton import drop_ledger, lifton
from tests.test_integration_pipeline import hermetic_pipeline  # noqa: F401


@pytest.fixture(autouse=True)
def _clean_ledger():
    drop_ledger.reset()
    yield
    drop_ledger.reset()


class TestLedger:
    def test_a_clean_run_reports_nothing(self, capsys):
        counts = {}
        manifest = types.SimpleNamespace(
            record_count=lambda name, value: counts.__setitem__(name, value))
        drop_ledger.report(manifest)
        assert drop_ledger.total() == 0
        # Every class is recorded as zero, so the manifest shape is stable...
        assert counts["dropped_features_total"] == 0
        assert all(counts[f"dropped_{kind}"] == 0 for kind in drop_ledger.CLASSES)
        # ...but a clean run does not grow a section saying nothing happened.
        assert "dropped" not in capsys.readouterr().err

    def test_counts_examples_and_the_summary(self, capsys):
        for index in range(3):
            drop_ledger.record("unresolvable_transcript", f"rna-{index}")
        drop_ledger.record("reference_protein_sequence", "rna-P")
        assert drop_ledger.counts() == {"unresolvable_transcript": 3,
                                        "reference_protein_sequence": 1}
        assert drop_ledger.total() == 4
        assert drop_ledger.examples("unresolvable_transcript") == \
            ["rna-0", "rna-1", "rna-2"]
        drop_ledger.report(None)
        message = capsys.readouterr().err
        assert "4 reference feature(s) were dropped" in message
        assert "3 x" in message and "rna-0" in message

    def test_examples_are_capped_but_counts_are_not(self):
        for index in range(drop_ledger.EXAMPLE_CAP + 25):
            drop_ledger.record("unresolvable_gene", f"gene-{index}")
        assert drop_ledger.counts()["unresolvable_gene"] == \
            drop_ledger.EXAMPLE_CAP + 25
        assert len(drop_ledger.examples("unresolvable_gene")) == \
            drop_ledger.EXAMPLE_CAP

    def test_an_undeclared_class_is_refused(self):
        # A counter nobody can interpret is barely better than no counter.
        with pytest.raises(KeyError, match="drop_ledger.CLASSES"):
            drop_ledger.record("something_new", "x")

    def test_every_class_has_an_explanation(self):
        for kind, sentence in drop_ledger.CLASSES.items():
            assert sentence and not sentence.endswith(".")
            assert kind.islower()

    def test_recording_is_safe_from_several_threads(self):
        # Step 7 dispatches to a thread pool, so these sites are reached
        # concurrently; a lost update would under-report the very number the
        # ledger exists to make trustworthy.
        def worker():
            for _ in range(500):
                drop_ledger.record("unresolvable_transcript")
        threads = [threading.Thread(target=worker) for _ in range(8)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        assert drop_ledger.counts()["unresolvable_transcript"] == 4000

    def test_reset_clears_between_runs(self):
        drop_ledger.record("unresolvable_gene", "g")
        drop_ledger.reset()
        assert drop_ledger.counts() == {} and drop_ledger.total() == 0


def _wrap(sequence):
    return "\n".join(textwrap.wrap(sequence, 60)) + "\n"


class TestEndToEnd:
    """A reference naming a transcript it does not declare."""

    def _workspace(self, work, *, break_the_reference):
        work.mkdir(parents=True, exist_ok=True)
        chromosome = ["A"] * 600
        cds = "ATG" + "GCT" * 31 + "TAA"
        chromosome[100:100 + len(cds)] = cds
        sequence = ">chr1\n" + _wrap("".join(chromosome))
        (work / "ref.fa").write_text(sequence)
        (work / "tgt.fa").write_text(sequence)
        rows = (
            "chr1\t{src}\tgene\t101\t199\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
            "chr1\t{src}\tmRNA\t101\t199\t.\t+\t.\tID=tx1;Parent=gene1\n"
            "chr1\t{src}\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
            "chr1\t{src}\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n")
        (work / "ref.gff3").write_text("##gff-version 3\n" + rows.format(src="test"))
        # The lifted transcript claims a reference id the reference lacks, which
        # is what the -copies bug produced ~4,400 times.
        lifted = rows.format(src="Liftoff")
        if break_the_reference:
            lifted = lifted.replace("ID=tx1;", "ID=tx_absent;").replace(
                "Parent=tx1", "Parent=tx_absent")
        (work / "liftoff.gff3").write_text("##gff-version 3\n" + lifted)
        # A real hit: an empty miniprot file cannot be loaded as a database.
        (work / "miniprot.gff3").write_text(
            "##gff-version 3\n"
            "chr1\tminiprot\tmRNA\t101\t199\t.\t+\t.\tID=MP1;Rank=1;"
            "Identity=1.0000;Target=tx1 1 32\n"
            "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.c1;Parent=MP1\n")
        (work / "out").mkdir()
        return work

    def _run(self, work):
        output = work / "out" / "lifton.gff3"
        lifton.run_all_lifton_steps(lifton.parse_args([
            str(work / "tgt.fa"), str(work / "ref.fa"),
            "-g", str(work / "ref.gff3"), "-L", str(work / "liftoff.gff3"),
            "-M", str(work / "miniprot.gff3"), "-o", str(output),
            "-ad", "RefSeq", "--force"]))
        return output.read_text()

    def test_a_clean_reference_drops_nothing(self, tmp_path, hermetic_pipeline,
                                             capsys):
        self._run(self._workspace(tmp_path / "clean", break_the_reference=False))
        assert drop_ledger.total() == 0
        assert "were dropped" not in capsys.readouterr().err

    def test_an_unresolvable_transcript_is_counted_and_reported(
            self, tmp_path, hermetic_pipeline, capsys):
        self._run(self._workspace(tmp_path / "broken", break_the_reference=True))
        assert drop_ledger.counts().get("unresolvable_transcript", 0) >= 1
        message = capsys.readouterr().err
        assert "were dropped" in message
        assert "exons and CDS" in message

    def test_the_manifest_carries_the_totals(self, tmp_path, hermetic_pipeline):
        import json
        work = self._workspace(tmp_path / "manifest", break_the_reference=True)
        self._run(work)
        manifest = json.loads(
            (work / "out" / "lifton_output" / "run_manifest.json").read_text())
        counts = manifest.get("counts") or {}
        assert counts["dropped_features_total"] >= 1
        assert counts["dropped_unresolvable_transcript"] >= 1
