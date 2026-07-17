"""Focused tests for dependency-free run provenance collection."""

from __future__ import annotations

import hashlib
import json
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
