"""Focused tests for the consolidated annotation metadata scan."""

from __future__ import annotations

import hashlib
from contextlib import contextmanager

import pytest

from lifton import annotation
from lifton.annotation_validator import (
    AnnotationScanResult,
    _attribute_values,
    scan_annotation,
    validate_annotation_file,
)
from lifton.exceptions import LiftOnInputError


def test_identity_attribute_parser_preserves_compatible_edge_cases():
    assert _attribute_values(
        "Dbxref=ID=decoy; ID = gene1 ;Parent=tx1, tx2;Note=Parent=fake"
    ) == {"ID": ("gene1",), "Parent": ("tx1", "tx2")}
    assert _attribute_values("ID=;Parent=") == {"ID": (), "Parent": ()}


def test_scan_rejects_empty_id_without_crashing(tmp_path):
    source = tmp_path / "empty-id.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chr1\tRefSeq\tgene\t1\t90\t.\t+\t.\tID=\n"
    )

    scan = scan_annotation(str(source))

    assert scan.data_lines == 1
    assert scan.duplicate_ids == ()
    assert any(finding.rule == "empty_id" for finding in scan.ncbi_findings)


def test_scan_is_stable_and_classifies_only_true_id_collisions(tmp_path):
    source = tmp_path / "annotation.gff3"
    source.write_text(
        "##gff-version 3\n"
        "#!genome-build test\n"
        "chr1\tRefSeq\tgene\t1\t90\t.\t+\t.\tID=g1\n"
        "chr1\tRefSeq\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tRefSeq\tCDS\t1\t30\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\tRefSeq\tCDS\t61\t90\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\tRefSeq\tgene\t100\t120\t.\t+\t.\tID=g1\n"
    )

    scan = scan_annotation(str(source), target_seqids={"chr1"})

    assert scan.sha256 == hashlib.sha256(source.read_bytes()).hexdigest()
    assert scan.size == source.stat().st_size
    assert scan.directives == ("##gff-version 3", "#!genome-build test")
    assert scan.discontinuous_cds_ids == ("c1",)
    assert scan.duplicate_id_map == {"g1": [3, 7]}
    assert scan.missing_target_seqids == ()
    validation = validate_annotation_file(str(source), scan_result=scan)
    assert not validation.is_valid
    assert validation.duplicate_ids == {"g1": [3, 7]}


def test_cds_repeat_with_changed_source_parent_or_strand_is_collision(tmp_path):
    variants = {
        "source": "chr1\tother\tCDS\t31\t60\t.\t+\t0\tID=c1;Parent=t1",
        "parent": "chr1\tRefSeq\tCDS\t31\t60\t.\t+\t0\tID=c1;Parent=t2",
        "strand": "chr1\tRefSeq\tCDS\t31\t60\t.\t-\t0\tID=c1;Parent=t1",
    }
    for name, second in variants.items():
        source = tmp_path / f"{name}.gff3"
        source.write_text(
            "##gff-version 3\n"
            "chr1\tRefSeq\tgene\t1\t90\t.\t+\t.\tID=g1\n"
            "chr1\tRefSeq\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
            "chr1\tRefSeq\tCDS\t1\t30\t.\t+\t0\tID=c1;Parent=t1\n"
            f"{second}\n"
        )
        assert "c1" in scan_annotation(str(source)).duplicate_id_map


def test_annotation_reuses_supplied_scan_without_rescanning(gff_standard, monkeypatch):
    scan = scan_annotation(str(gff_standard))

    def unexpected_scan(*args, **kwargs):
        raise AssertionError("matching supplied scan must be reused")

    monkeypatch.setattr(annotation, "scan_annotation", unexpected_scan)
    loaded = annotation.Annotation(
        str(gff_standard),
        False,
        False,
        force=True,
        scan_result=scan,
    )
    loaded.db_connection.conn.close()


def test_orphan_and_target_seqid_findings_are_reused(tmp_path):
    source = tmp_path / "orphans.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chrX\tRefSeq\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=missing\n"
    )

    scan = scan_annotation(str(source), target_seqids={"chr1"})

    assert scan.orphan_parents == ("missing",)
    assert scan.missing_target_seqids == ("chrX",)
    assert any(finding.rule == "dangling_parent" for finding in scan.ncbi_findings)
    assert any(finding.rule == "unknown_seqid" for finding in scan.ncbi_findings)


def test_database_suffixes_row_keys_without_rewriting_logical_cds_id(tmp_path):
    source = tmp_path / "cds.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chr1\tRefSeq\tgene\t1\t90\t.\t+\t.\tID=g1\n"
        "chr1\tRefSeq\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tRefSeq\tCDS\t1\t30\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\tRefSeq\tCDS\t61\t90\t.\t+\t0\tID=c1;Parent=t1\n"
    )

    loaded = annotation.Annotation(str(source), False, False, force=True)
    cds_rows = list(loaded.db_connection.features_of_type("CDS"))

    assert len({row.id for row in cds_rows}) == 2
    assert [row.attributes["ID"] for row in cds_rows] == [["c1"], ["c1"]]
    loaded.db_connection.conn.close()


def test_gffbase_suffixes_row_keys_without_rewriting_logical_cds_id(
        tmp_path, monkeypatch):
    monkeypatch.setenv("LIFTON_DISABLE_RTREE", "1")
    source = tmp_path / "cds.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chr1\tRefSeq\tgene\t1\t90\t.\t+\t.\tID=g1\n"
        "chr1\tRefSeq\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tRefSeq\tCDS\t1\t30\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\tRefSeq\tCDS\t61\t90\t.\t+\t0\tID=c1;Parent=t1\n"
    )

    loaded = annotation.Annotation(
        str(source), False, False, force=True, backend="gffbase",
    )
    cds_rows = list(loaded.db_connection.features_of_type("CDS"))

    assert len({row.id for row in cds_rows}) == 2
    assert [row.attributes["ID"] for row in cds_rows] == [["c1"], ["c1"]]
    loaded.db_connection.conn.close()


def test_every_scan_error_stops_database_build():
    loaded = annotation.Annotation.__new__(annotation.Annotation)
    loaded.file_name = "unstable.gff3"
    loaded.scan_result = AnnotationScanResult(
        file_path=loaded.file_name,
        errors=("Annotation changed while it was being scanned",),
    )

    with pytest.raises(SystemExit):
        loaded._run_preflight_validation()


def test_source_replacement_while_waiting_for_cache_lock_fails_closed(
        tmp_path, monkeypatch):
    source = tmp_path / "annotation.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chr1\tRefSeq\tgene\t1\t90\t.\t+\t.\tID=g1\n"
        "chr1\tRefSeq\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tRefSeq\texon\t1\t90\t.\t+\t.\tID=e1;Parent=t1\n"
    )
    replacement = tmp_path / "replacement.gtf"
    replacement.write_text(
        "chr1\tother\ttranscript\t1\t90\t.\t+\t.\t"
        'gene_id "changed"; transcript_id "changed.t1";\n'
    )

    @contextmanager
    def replacing_lock(_db_path):
        replacement.replace(source)
        yield

    def unexpected_build(*args, **kwargs):
        raise AssertionError("unvalidated replacement must not be ingested")

    monkeypatch.setattr(annotation.annotation_cache, "cache_lock", replacing_lock)
    monkeypatch.setattr(annotation.gffutils, "create_db", unexpected_build)

    with pytest.raises(LiftOnInputError, match="changed after preflight"):
        annotation.Annotation(str(source), False, False, force=True)
