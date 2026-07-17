"""Focused tests for dependency-free run provenance collection."""

from __future__ import annotations

import hashlib
import json
import threading
from pathlib import Path

import pytest

from lifton import run_manifest


def test_sanitizers_redact_credentials_and_filter_environment() -> None:
    argv = run_manifest.sanitize_argv(
        [
            "lifton",
            "--threads=8",
            "--api-token=abc123",
            "--password",
            "hunter2",
            "https://user:secret@example.org/input.fa",
        ]
    )
    assert argv == [
        "lifton",
        "--threads=8",
        "--api-token=<redacted>",
        "--password",
        "<redacted>",
        "https://<redacted>@example.org/input.fa",
    ]

    options = run_manifest.sanitize_options(
        {
            "threads": 8,
            "api_key": "do-not-store",
            "nested": {"auth_token": "also-secret", "mode": "stream"},
            "output": Path("result.gff3"),
        }
    )
    assert options == {
        "threads": 8,
        "api_key": "<redacted>",
        "nested": {"auth_token": "<redacted>", "mode": "stream"},
        "output": "result.gff3",
    }

    environment = run_manifest.collect_environment(
        {
            "PATH": "/private/bin",
            "HOME": "/private/home",
            "LIFTON_BACKEND": "gffbase",
            "LIFTON_AUTH_TOKEN": "do-not-store",
        }
    )
    assert environment == {"LIFTON_AUTH_TOKEN": "<redacted>", "LIFTON_BACKEND": "gffbase"}


def test_fingerprint_input_hashes_files_and_tolerates_missing(tmp_path: Path) -> None:
    source = tmp_path / "reference.fa"
    source.write_bytes(b">chr1\nACGT\n")

    fingerprint = run_manifest.fingerprint_input(source, block_size=3)
    assert fingerprint["exists"] is True
    assert fingerprint["is_file"] is True
    assert fingerprint["size_bytes"] == len(b">chr1\nACGT\n")
    assert fingerprint["sha256"] == hashlib.sha256(b">chr1\nACGT\n").hexdigest()
    assert "error" not in fingerprint

    missing = run_manifest.fingerprint_input(tmp_path / "missing.gff3")
    assert missing["exists"] is False
    assert missing["sha256"] is None
    assert missing["error"].startswith("FileNotFoundError:")


def test_missing_dependencies_tools_and_nonrepository_are_nonfatal(tmp_path: Path) -> None:
    distributions = run_manifest.collect_dependency_versions(("package-that-cannot-exist-lifton-test",))
    assert distributions == {"package-that-cannot-exist-lifton-test": None}

    tool = run_manifest.probe_tool_version((str(tmp_path / "missing-tool"), "--version"))
    assert tool["available"] is False
    assert tool["version"] is None
    assert tool["error"] == "not found"

    git = run_manifest.collect_git_metadata(tmp_path)
    assert git["repository"] is False
    assert git["commit"] is None
    assert git["branch"] is None


def test_manifest_records_timing_counts_validation_choices_and_json(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("LIFTON_TEST_MODE", "1")
    monkeypatch.setenv("UNRELATED_SECRET", "must-not-appear")
    source = tmp_path / "reference.gff3"
    source.write_text("##gff-version 3\n", encoding="utf-8")

    manifest = run_manifest.RunManifest(
        argv=["lifton", "--token", "hidden", str(source)],
        options={"threads": 4, "password": "hidden"},
        inputs={"annotation": source, "missing": tmp_path / "not-there.fa"},
        backend={"annotation": "gffbase"},
        cache={"minimap2_index": "miss"},
        dependency_names=(),
        tool_commands={},
        collect_git=False,
    )
    manifest.set_backend_choice("alignment", "subprocess")
    manifest.set_cache_choice("miniprot", "hit")
    with manifest.phase("alignment", {"workers": 4}):
        manifest.record_count("genes", 2)
        assert manifest.increment_count("genes") == 3
    manifest.record_failure("writer", "token=leaked-value", fatal=False)
    manifest.record_validation("gff3", passed=True, errors=0, warnings=1)
    manifest.finish("success")

    output = tmp_path / "nested" / "run.json"
    manifest.write(output)
    stored = json.loads(output.read_text(encoding="utf-8"))

    assert stored["schema_version"] == 2
    assert stored["run"]["status"] == "success"
    assert stored["run"]["argv"][2] == "<redacted>"
    assert stored["run"]["options"]["password"] == "<redacted>"
    assert stored["run"]["environment"] == {"LIFTON_TEST_MODE": "1"}
    assert stored["run"]["backend"] == {"alignment": "subprocess", "annotation": "gffbase"}
    assert stored["run"]["cache"] == {"minimap2_index": "miss", "miniprot": "hit"}
    assert stored["inputs"]["annotation"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert stored["inputs"]["missing"]["exists"] is False
    fingerprint_phase = stored["phases"]["fingerprint_inputs"]
    assert fingerprint_phase["status"] == "success"
    assert fingerprint_phase["details"] == {
        "algorithm": "sha256",
        "block_size_bytes": 1024 * 1024,
        "execution": "background-thread",
    }
    assert fingerprint_phase["duration_seconds"] >= 0
    assert fingerprint_phase["metrics"]["input_count"] == 2
    assert fingerprint_phase["metrics"]["hashed_files"] == 1
    assert fingerprint_phase["metrics"]["unavailable_inputs"] == 1
    assert fingerprint_phase["metrics"]["bytes_hashed"] == source.stat().st_size
    assert fingerprint_phase["metrics"]["changed_inputs"] == 0
    assert fingerprint_phase["metrics"]["join_wait_seconds"] >= 0
    assert stored["phases"]["alignment"]["status"] == "success"
    assert stored["phases"]["alignment"]["duration_seconds"] >= 0
    assert stored["counts"] == {"genes": 3}
    assert stored["failures"][0]["message"] == "token=<redacted>"
    assert stored["validation"]["gff3"] == {
        "details": {},
        "errors": 0,
        "passed": True,
        "warnings": 1,
    }
    assert stored["software"]["lifton"]
    assert stored["software"]["python"]["version"]
    assert stored["resources"]["start"]["peak_bytes"] is None or stored["resources"]["start"]["peak_bytes"] > 0

    snapshot = manifest.to_dict()
    snapshot["counts"]["genes"] = 999
    assert manifest.to_dict()["counts"]["genes"] == 3


def test_manifest_hashes_in_background_and_joins_before_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "large.fa"
    source.write_bytes(b">chr1\n" + b"ACGT" * 1024)
    worker_started = threading.Event()
    release_worker = threading.Event()
    publisher_started = threading.Event()
    original = run_manifest.fingerprint_inputs

    def delayed_fingerprints(inputs):
        worker_started.set()
        if not release_worker.wait(timeout=5):
            raise RuntimeError("test fingerprint worker timed out")
        return original(inputs)

    monkeypatch.setattr(run_manifest, "fingerprint_inputs", delayed_fingerprints)
    manifest = run_manifest.RunManifest(
        inputs={"target": source},
        dependency_names=(),
        tool_commands={},
        collect_git=False,
    )
    assert worker_started.wait(timeout=1)
    assert manifest.to_dict()["inputs"]["target"]["fingerprint_status"] == "pending"
    assert manifest.to_dict()["phases"]["fingerprint_inputs"]["status"] == "running"

    output = tmp_path / "run.json"
    publisher_errors = []

    def publish():
        publisher_started.set()
        try:
            manifest.finish("success")
            manifest.write(output)
        except BaseException as exc:  # pragma: no cover - asserted below
            publisher_errors.append(exc)

    publisher = threading.Thread(target=publish)
    publisher.start()
    assert publisher_started.wait(timeout=1)
    assert publisher.is_alive()
    assert not output.exists()

    release_worker.set()
    publisher.join(timeout=5)
    assert not publisher.is_alive()
    assert publisher_errors == []
    stored = json.loads(output.read_text(encoding="utf-8"))
    assert stored["inputs"]["target"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert stored["inputs"]["target"]["fingerprint_status"] == "complete"
    assert stored["phases"]["fingerprint_inputs"]["status"] == "success"


def test_manifest_records_worker_failure_without_pending_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "target.fa"
    source.write_bytes(b">chr1\nACGT\n")

    def failed_fingerprints(_inputs):
        raise RuntimeError("worker boom")

    monkeypatch.setattr(run_manifest, "fingerprint_inputs", failed_fingerprints)
    manifest = run_manifest.RunManifest(
        inputs={"target": source},
        dependency_names=(),
        tool_commands={},
        collect_git=False,
    )
    manifest.finish("success")
    snapshot = manifest.to_dict()
    assert snapshot["run"]["status"] == "success"
    assert snapshot["inputs"]["target"]["fingerprint_status"] == "unavailable"
    assert snapshot["inputs"]["target"]["sha256"] is None
    assert "worker boom" in snapshot["inputs"]["target"]["error"]
    phase = snapshot["phases"]["fingerprint_inputs"]
    assert phase["status"] == "failed"
    assert phase["metrics"]["unavailable_inputs"] == 1
    assert snapshot["failures"][0]["phase"] == "fingerprint_inputs"
    assert snapshot["failures"][0]["message"] == "worker boom"


def test_fingerprint_rejects_file_changed_while_hashing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "changing.fa"
    source.write_bytes(b">chr1\nACGT\n")
    real_sha256 = hashlib.sha256
    mutated = False

    class MutatingDigest:
        def __init__(self):
            self.digest = real_sha256()

        def update(self, block):
            nonlocal mutated
            self.digest.update(block)
            if not mutated:
                mutated = True
                with source.open("ab") as handle:
                    handle.write(b"N")

        def hexdigest(self):
            return self.digest.hexdigest()

    monkeypatch.setattr(run_manifest.hashlib, "sha256", MutatingDigest)
    fingerprint = run_manifest.fingerprint_input(source, block_size=3)
    assert fingerprint["changed_during_hash"] is True
    assert fingerprint["fingerprint_status"] == "changed"
    assert fingerprint["sha256"] is None
    assert fingerprint["error"] == "input changed while computing SHA-256"
    assert fingerprint["stat_after"]["size_bytes"] > fingerprint["size_bytes"]

    source.write_bytes(b">chr1\nTGCA\n")
    mutated = False
    manifest = run_manifest.RunManifest(
        inputs={"target": source},
        dependency_names=(),
        tool_commands={},
        collect_git=False,
    )
    manifest.finish("success")
    snapshot = manifest.to_dict()
    assert snapshot["inputs"]["target"]["sha256"] is None
    assert snapshot["phases"]["fingerprint_inputs"]["status"] == "failed"
    assert snapshot["phases"]["fingerprint_inputs"]["metrics"]["changed_inputs"] == 1
    assert snapshot["failures"][0]["phase"] == "fingerprint_inputs"


def test_fingerprint_classifies_deleted_path_as_changed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "deleted-during-hash.fa"
    source.write_bytes(b">chr1\nACGT\n")
    real_sha256 = hashlib.sha256
    deleted = False

    class DeletingDigest:
        def __init__(self):
            self.digest = real_sha256()

        def update(self, block):
            nonlocal deleted
            self.digest.update(block)
            if not deleted:
                deleted = True
                source.unlink()

        def hexdigest(self):
            return self.digest.hexdigest()

    monkeypatch.setattr(run_manifest.hashlib, "sha256", DeletingDigest)
    fingerprint = run_manifest.fingerprint_input(source, block_size=3)
    assert fingerprint["changed_during_hash"] is True
    assert fingerprint["fingerprint_status"] == "changed"
    assert fingerprint["sha256"] is None
    assert fingerprint["error"].startswith(
        "input path changed while computing SHA-256: FileNotFoundError:"
    )


def test_fingerprint_classifies_delete_between_stat_and_open_as_changed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "deleted-before-open.fa"
    source.write_bytes(b">chr1\nACGT\n")
    original_open = Path.open
    deleted = False

    def deleting_open(path, *args, **kwargs):
        nonlocal deleted
        if path == source and not deleted:
            deleted = True
            source.unlink()
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", deleting_open)
    fingerprint = run_manifest.fingerprint_input(source)

    assert fingerprint["exists"] is True
    assert fingerprint["is_file"] is True
    assert fingerprint["changed_during_hash"] is True
    assert fingerprint["fingerprint_status"] == "changed"
    assert fingerprint["sha256"] is None
    assert fingerprint["error"].startswith(
        "input path changed or became unreadable while computing SHA-256: "
        "FileNotFoundError:"
    )


def test_fingerprint_keeps_stable_unreadable_file_unavailable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "unreadable.fa"
    source.write_bytes(b">chr1\nACGT\n")
    original_open = Path.open

    def denied_open(path, *args, **kwargs):
        if path == source:
            raise PermissionError("permission denied")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", denied_open)
    fingerprint = run_manifest.fingerprint_input(source)

    assert fingerprint["exists"] is True
    assert fingerprint["changed_during_hash"] is False
    assert fingerprint["fingerprint_status"] == "unavailable"
    assert fingerprint["sha256"] is None
    assert fingerprint["error"] == "PermissionError: permission denied"


def test_phase_context_records_failure_and_reraises() -> None:
    manifest = run_manifest.RunManifest(
        dependency_names=(),
        tool_commands={},
        collect_git=False,
    )

    with pytest.raises(RuntimeError, match="boom"):
        with manifest.phase("failing"):
            raise RuntimeError("boom password=hidden")

    snapshot = manifest.to_dict()
    assert snapshot["phases"]["failing"]["status"] == "failed"
    assert snapshot["failures"][0]["type"] == "RuntimeError"
    assert snapshot["failures"][0]["message"] == "boom password=<redacted>"


def test_atomic_write_keeps_previous_file_and_cleans_temp_on_failure(tmp_path: Path) -> None:
    destination = tmp_path / "run.json"
    destination.write_text('{"old": true}\n', encoding="utf-8")

    with pytest.raises(TypeError):
        run_manifest.atomic_write_json(destination, {"bad": object()})

    assert destination.read_text(encoding="utf-8") == '{"old": true}\n'
    assert list(tmp_path.glob(".run.json.*.tmp")) == []
