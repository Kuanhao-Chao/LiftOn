"""Focused tests for atomic GFF3 output publication."""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from lifton.output_transaction import (
    OutputTransaction,
    OutputTransactionState,
    OutputTransactionStateError,
)


def test_file_destination_is_unchanged_until_commit(tmp_path):
    destination = tmp_path / "lifted.gff3"
    destination.write_text("old output\n")
    transaction = OutputTransaction(destination)

    assert transaction.working_path.parent == destination.parent
    assert transaction.working_path != destination
    transaction.open().write("##gff-version 3\nchr1\tLiftOn\tgene\t1\t3\t.\t+\t.\tID=g1\n")
    transaction.close()

    assert destination.read_text() == "old output\n"
    assert transaction.state is OutputTransactionState.CLOSED

    assert transaction.commit() == destination
    assert destination.read_text().startswith("##gff-version 3\n")
    assert not transaction.working_path.exists()
    assert transaction.state is OutputTransactionState.COMMITTED


def test_commit_publishes_file_with_os_replace(tmp_path, monkeypatch):
    destination = tmp_path / "result.gff3"
    transaction = OutputTransaction(destination)
    transaction.open().write("complete\n")
    replacements: list[tuple[Path, Path]] = []

    real_replace = __import__("os").replace

    def recording_replace(source, target):
        replacements.append((Path(source), Path(target)))
        real_replace(source, target)

    monkeypatch.setattr("lifton.output_transaction.os.replace", recording_replace)
    transaction.commit()

    assert replacements == [(transaction.working_path, destination)]
    assert destination.read_text() == "complete\n"


def test_abort_preserves_predictable_partial_artifact(tmp_path):
    destination = tmp_path / "result.gff3"
    destination.write_text("previous good output\n")
    transaction = OutputTransaction(destination)
    transaction.open().write("incomplete output\n")

    partial = transaction.abort()

    assert partial == tmp_path / "result.partial.gff3"
    assert partial.read_text() == "incomplete output\n"
    assert destination.read_text() == "previous good output\n"
    assert not transaction.working_path.exists()
    assert transaction.state is OutputTransactionState.ABORTED


def test_stdout_receives_nothing_before_commit(tmp_path):
    stdout = io.StringIO()
    transaction = OutputTransaction("stdout", stdout=stdout, work_dir=tmp_path)
    transaction.open().write("buffered GFF3\n")
    transaction.close()

    assert stdout.getvalue() == ""

    assert transaction.commit() is None
    assert stdout.getvalue() == "buffered GFF3\n"
    assert not transaction.working_path.exists()


def test_stdout_defaults_to_current_sys_stdout(tmp_path, monkeypatch):
    stdout = io.StringIO()
    monkeypatch.setattr("sys.stdout", stdout)
    transaction = OutputTransaction("stdout", work_dir=tmp_path)
    transaction.open().write("default stream\n")

    transaction.commit()

    assert stdout.getvalue() == "default stream\n"


def test_stdout_abort_preserves_artifact_without_writing_stream(tmp_path):
    stdout = io.StringIO()
    transaction = OutputTransaction("stdout", stdout=stdout, work_dir=tmp_path)
    transaction.open().write("partial stdout output\n")

    partial = transaction.abort()

    assert stdout.getvalue() == ""
    assert partial == tmp_path / "lifton.partial.gff3"
    assert partial.read_text() == "partial stdout output\n"


def test_terminal_operations_are_idempotent_and_conflicts_are_errors(tmp_path):
    committed = OutputTransaction(tmp_path / "committed.gff3")
    committed.open().write("done\n")
    committed.close()
    committed.close()
    assert committed.commit() == tmp_path / "committed.gff3"
    assert committed.commit() == tmp_path / "committed.gff3"
    committed.close()
    with pytest.raises(OutputTransactionStateError, match="cannot abort a committed"):
        committed.abort()

    aborted = OutputTransaction(tmp_path / "aborted.gff3")
    aborted.open().write("not done\n")
    expected_partial = aborted.abort()
    assert aborted.abort() == expected_partial
    aborted.close()
    with pytest.raises(OutputTransactionStateError, match="cannot commit an aborted"):
        aborted.commit()


def test_closed_transaction_cannot_be_reopened(tmp_path):
    transaction = OutputTransaction(tmp_path / "closed.gff3")
    transaction.close()

    with pytest.raises(OutputTransactionStateError, match="state 'closed'"):
        transaction.open()

    transaction.abort()


def test_context_manager_commits_or_preserves_partial_on_error(tmp_path):
    destination = tmp_path / "success.gff3"
    with OutputTransaction(destination) as stream:
        stream.write("success\n")
    assert destination.read_text() == "success\n"

    failed_destination = tmp_path / "failed.gff3"
    transaction = OutputTransaction(failed_destination)
    with pytest.raises(ValueError, match="synthetic failure"):
        with transaction as stream:
            stream.write("recoverable\n")
            raise ValueError("synthetic failure")

    assert not failed_destination.exists()
    assert transaction.partial_path.read_text() == "recoverable\n"


def test_missing_staging_directory_is_reported(tmp_path):
    missing = tmp_path / "missing"

    with pytest.raises(FileNotFoundError, match="staging directory"):
        OutputTransaction(missing / "out.gff3")
