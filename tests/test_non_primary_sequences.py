"""Tell the user their reference carries alternate-locus / patch sequences.

14.6 % of GRCh38 RefSeq coding genes sit on _alt/_fix contigs as copies of a
gene on the primary assembly, competing for the same target locus. LiftOn does
not override which one wins (notes/primary_assembly_id_preference_nogo.md), but
it should not leave the situation to be discovered by measurement either.
"""
from __future__ import annotations

import textwrap
import types

import pytest

from lifton import coreutils, lifton
from tests.test_integration_pipeline import hermetic_pipeline  # noqa: F401


class TestClassification:
    @pytest.mark.parametrize("seqid, expected", [
        ("chr1", None),
        ("chr17_KI270857v1_alt", "alt"),
        ("chr1_KI270762v1_fix", "fix"),
        # Unlocalized and unplaced scaffolds ARE the primary assembly; they are
        # not copies of anything and must not be reported.
        ("chr1_KI270706v1_random", None),
        ("chrUn_GL000195v1", None),
        # An accession-style name cannot be classified from the name alone.
        ("NC_000001.11", None),
    ])
    def test_only_alt_and_fix_are_non_primary(self, seqid, expected):
        assert coreutils.non_primary_seqid_class(seqid) == expected


class TestReport:
    def _ref_db(self, rows):
        connection = types.SimpleNamespace(execute=lambda _sql: iter(rows))
        return types.SimpleNamespace(db_connection=connection)

    def _args(self):
        counts = {}
        manifest = types.SimpleNamespace(
            record_count=lambda name, value: counts.__setitem__(name, value))
        args = types.SimpleNamespace(_run_manifest=manifest)
        return args, counts

    def test_a_primary_only_reference_says_nothing(self, capsys):
        args, counts = self._args()
        lifton._report_non_primary_sequences(
            self._ref_db([("chr1", 500), ("chr2", 400)]), args)
        assert counts == {}
        assert "alternate-locus" not in capsys.readouterr().err

    def test_alt_and_fix_genes_are_counted_and_reported(self, capsys):
        args, counts = self._args()
        lifton._report_non_primary_sequences(self._ref_db([
            ("chr1", 500), ("chr17_KI270857v1_alt", 7),
            ("chr1_KI270762v1_fix", 3), ("chrUn_GL000195v1", 2),
        ]), args)
        assert counts == {"reference_non_primary_sequences": 2,
                          "reference_genes_on_non_primary_sequences": 10}
        message = capsys.readouterr().err
        assert "2 alternate-locus" in message and "10 gene(s)" in message

    def test_a_backend_without_the_table_is_not_a_failure(self):
        def _raise(_sql):
            raise RuntimeError("no such table: features")
        args, counts = self._args()
        lifton._report_non_primary_sequences(
            types.SimpleNamespace(
                db_connection=types.SimpleNamespace(execute=_raise)), args)
        assert counts == {}

    def test_a_run_without_a_manifest_still_reports(self, capsys):
        args = types.SimpleNamespace()
        lifton._report_non_primary_sequences(
            self._ref_db([("chr9_KI270717v1_alt", 4)]), args)
        assert "alternate-locus" in capsys.readouterr().err
