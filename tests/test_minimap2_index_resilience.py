"""Regression coverage for the unresolved Git LFS minimap2 index in GH #57.

These tests intentionally use tiny synthetic files and fake subprocesses.  They
exercise LiftOn's validation and failure boundaries without requiring minimap2
or a hydrated copy of the 152 MB chr22 index.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from lifton.exceptions import LiftOnAlignmentError
from lifton import run_liftoff as run_liftoff_module
from lifton.liftoff import align_features, liftoff_main


LFS_POINTER = (
    b"version https://git-lfs.github.com/spec/v1\n"
    b"oid sha256:a6c9eab78d65520a3628bdc8ccf6c3ac5ea073f0a681b85d306fbae83386e6f7\n"
    b"size 151911282\n"
)


def _args(run_dir: Path) -> SimpleNamespace:
    return SimpleNamespace(
        directory=str(run_dir),
        mm2_options="-a --eqx -N 50 -p 0.5",
    )


def _target_fasta(tmp_path: Path) -> Path:
    target = tmp_path / "target.fa"
    target.write_text(">chr1\nACGTACGT\n")
    return target


def _successful_index_builder(calls):
    def run(command, **kwargs):
        calls.append((command, kwargs))
        output = Path(command[command.index("-d") + 1])
        output.write_bytes(align_features.MMI_MAGIC + b"synthetic-index")
        return subprocess.CompletedProcess(command, 0)

    return run


@pytest.mark.parametrize(
    ("payload", "expected"),
    [
        (None, "missing"),
        (b"MMI\x02synthetic-index", "valid"),
        (LFS_POINTER, "lfs_pointer"),
        (b"not a minimap2 index", "invalid"),
        (b"MMI", "invalid"),
    ],
)
def test_classify_minimap2_index(tmp_path, payload, expected):
    index = tmp_path / "target.fa.mmi"
    if payload is not None:
        index.write_bytes(payload)

    assert align_features.classify_minimap2_index(index) == expected


def test_valid_sidecar_index_is_reused_without_subprocess(tmp_path, monkeypatch):
    target = _target_fasta(tmp_path)
    sidecar = Path(str(target) + ".mmi")
    sidecar.write_bytes(align_features.MMI_MAGIC + b"existing-index")

    def unexpected_run(*args, **kwargs):  # pragma: no cover - assertion path
        raise AssertionError("a valid cached minimap2 index must be reused")

    monkeypatch.setattr(align_features.subprocess, "run", unexpected_run)
    result = align_features.build_minimap2_index(
        str(target), _args(tmp_path / "run"), "2", "minimap2",
        target_prefix="target_all",
    )

    assert Path(result) == sidecar


def test_missing_index_is_built_in_run_directory(tmp_path, monkeypatch):
    target = _target_fasta(tmp_path)
    run_dir = tmp_path / "run"
    calls = []
    monkeypatch.setattr(
        align_features.subprocess, "run", _successful_index_builder(calls)
    )

    result = Path(
        align_features.build_minimap2_index(
            str(target), _args(run_dir), "2", "minimap2",
            target_prefix="target_all",
        )
    )

    assert run_dir in result.parents
    assert result.read_bytes().startswith(align_features.MMI_MAGIC)
    assert calls and calls[0][1]["check"] is True
    assert str(target) in calls[0][0]


@pytest.mark.parametrize("bad_index", [LFS_POINTER, b"corrupt-or-truncated"])
def test_bad_sidecar_is_preserved_and_rebuilt_locally(
    tmp_path, monkeypatch, bad_index
):
    target = _target_fasta(tmp_path)
    sidecar = Path(str(target) + ".mmi")
    sidecar.write_bytes(bad_index)
    run_dir = tmp_path / "run"
    calls = []
    monkeypatch.setattr(
        align_features.subprocess, "run", _successful_index_builder(calls)
    )

    result = Path(
        align_features.build_minimap2_index(
            str(target), _args(run_dir), "1", "minimap2",
            target_prefix="target_all",
        )
    )

    assert result != sidecar
    assert run_dir in result.parents
    assert result.read_bytes().startswith(align_features.MMI_MAGIC)
    assert sidecar.read_bytes() == bad_index
    assert len(calls) == 1


def test_index_build_subprocess_failure_is_actionable_and_atomic(
    tmp_path, monkeypatch
):
    target = _target_fasta(tmp_path)
    run_dir = tmp_path / "run"

    def fail(command, **kwargs):
        raise subprocess.CalledProcessError(17, command, stderr="bad index")

    monkeypatch.setattr(align_features.subprocess, "run", fail)

    with pytest.raises(LiftOnAlignmentError) as exc:
        align_features.build_minimap2_index(
            str(target), _args(run_dir), "1", "minimap2",
            target_prefix="target_all",
        )

    message = str(exc.value).lower()
    assert "minimap2" in message
    assert "index" in message
    assert str(target) in str(exc.value)
    assert not list(run_dir.glob("*.tmp*"))
    assert not any(
        p.read_bytes().startswith(align_features.MMI_MAGIC)
        for p in run_dir.glob("*.mmi")
    )


def test_successful_builder_that_writes_invalid_index_is_rejected(
    tmp_path, monkeypatch
):
    target = _target_fasta(tmp_path)
    run_dir = tmp_path / "run"

    def write_invalid(command, **kwargs):
        Path(command[command.index("-d") + 1]).write_bytes(b"not-an-index")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(align_features.subprocess, "run", write_invalid)

    with pytest.raises(LiftOnAlignmentError, match="(?i)index"):
        align_features.build_minimap2_index(
            str(target), _args(run_dir), "1", "minimap2",
            target_prefix="target_all",
        )

    assert not list(run_dir.glob("*.tmp*"))


def test_stale_sidecar_is_ignored_and_rebuilt_locally(tmp_path, monkeypatch):
    target = _target_fasta(tmp_path)
    sidecar = Path(str(target) + ".mmi")
    sidecar.write_bytes(align_features.MMI_MAGIC + b"stale-index")
    # A target modified after its sidecar invalidates that cache even when its
    # MMI magic still looks correct.
    target.touch()
    target_mtime = target.stat().st_mtime_ns
    sidecar_mtime = max(0, target_mtime - 1_000_000)
    sidecar.touch()
    import os
    os.utime(sidecar, ns=(sidecar_mtime, sidecar_mtime))

    calls = []
    run_dir = tmp_path / "run"
    monkeypatch.setattr(
        align_features.subprocess, "run", _successful_index_builder(calls)
    )

    result = Path(
        align_features.build_minimap2_index(
            str(target), _args(run_dir), "1", "minimap2",
            target_prefix="target_all",
        )
    )

    assert result.parent == run_dir
    assert result != sidecar
    assert len(calls) == 1


def test_sam_without_sequence_dictionary_is_rejected(tmp_path):
    target = _target_fasta(tmp_path)
    sam = tmp_path / "empty-header.sam"
    sam.write_text("@HD\tVN:1.6\tSO:unknown\n")

    with pytest.raises(LiftOnAlignmentError) as exc:
        align_features.validate_sam_target(
            str(sam), str(target), index_file=str(target) + ".mmi"
        )

    message = str(exc.value)
    assert "@SQ" in message or "sequence" in message.lower()
    assert str(target) in message


@pytest.mark.parametrize(
    "sq_line",
    [
        "@SQ\tSN:chr2\tLN:8\n",
        "@SQ\tSN:chr1\tLN:9\n",
    ],
)
def test_sam_dictionary_must_match_target_names_and_lengths(tmp_path, sq_line):
    target = _target_fasta(tmp_path)
    sam = tmp_path / "wrong-target.sam"
    sam.write_text("@HD\tVN:1.6\tSO:unknown\n" + sq_line)

    with pytest.raises(LiftOnAlignmentError) as exc:
        align_features.validate_sam_target(
            str(sam), str(target), index_file=str(target) + ".mmi"
        )

    message = str(exc.value).lower()
    assert "target" in message
    assert "match" in message or "sequence" in message


def test_valid_sam_with_only_unmapped_queries_is_not_mislabeled_corrupt(tmp_path):
    """Zero biological hits are valid SAM, not evidence of a broken index."""
    target = _target_fasta(tmp_path)
    sam = tmp_path / "unmapped.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:unknown\n"
        "@SQ\tSN:chr1\tLN:8\n"
        "gene1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
    )

    assert align_features.validate_sam_target(str(sam), str(target)) is None


def test_lifton_owned_sam_output_wins_over_user_mm2_output(tmp_path, monkeypatch):
    target = _target_fasta(tmp_path)
    final_sam = tmp_path / "final.sam"
    user_sam = tmp_path / "user.sam"
    observed = {}

    def fake_run(command, **kwargs):
        observed["command"] = command
        output_positions = [i for i, token in enumerate(command) if token == "-o"]
        actual_output = Path(command[output_positions[-1] + 1])
        actual_output.write_text("@SQ\tSN:chr1\tLN:8\n")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(align_features.subprocess, "run", fake_run)
    align_features._run_minimap2_to_sam(
        ["minimap2", "-a", "-o", str(user_sam), "-t", "1"],
        ["target.mmi", "queries.fa"], str(final_sam), str(target),
        {"chr1": 8}, index_file="target.mmi",
    )

    output_positions = [
        i for i, token in enumerate(observed["command"]) if token == "-o"
    ]
    assert Path(observed["command"][output_positions[-1] + 1]) != user_sam
    assert final_sam.exists()
    assert not user_sam.exists()


def test_failed_atomic_sam_keeps_previous_output_and_removes_temp(
    tmp_path, monkeypatch
):
    target = _target_fasta(tmp_path)
    final_sam = tmp_path / "final.sam"
    final_sam.write_text("previous-good-output\n")

    def fail(command, **kwargs):
        raise subprocess.CalledProcessError(9, command)

    monkeypatch.setattr(align_features.subprocess, "run", fail)
    with pytest.raises(LiftOnAlignmentError, match="alignment"):
        align_features._run_minimap2_to_sam(
            ["minimap2", "-a", "-t", "1"],
            ["target.mmi", "queries.fa"], str(final_sam), str(target),
            {"chr1": 8}, index_file="target.mmi",
        )

    assert final_sam.read_text() == "previous-good-output\n"
    assert not list(tmp_path.glob(".*final.sam*.tmp"))


@pytest.mark.parametrize("failure_type", ["header", "command"])
def test_reused_index_failure_rebuilds_and_retries_once(
    tmp_path, monkeypatch, failure_type
):
    target = _target_fasta(tmp_path)
    sidecar = Path(str(target) + ".mmi")
    sidecar.write_bytes(align_features.MMI_MAGIC + b"cached")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "reference_all_genes.fa").write_text(">gene1\nACGT\n")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGTACGT\n")
    rebuilt = run_dir / "target_all.mmi"
    build_calls = []
    align_calls = []

    def fake_build(*args, **kwargs):
        build_calls.append(kwargs)
        if kwargs.get("force_rebuild"):
            rebuilt.write_bytes(align_features.MMI_MAGIC + b"rebuilt")
            return str(rebuilt)
        return str(sidecar)

    def fake_align(command, positional_inputs, output_file, *args, **kwargs):
        align_calls.append((command, positional_inputs))
        if len(align_calls) == 1:
            if failure_type == "header":
                raise align_features._SamTargetMismatch("stale header")
            raise align_features._Minimap2CommandError("incompatible index")
        Path(output_file).write_text("@SQ\tSN:chr1\tLN:8\n")

    monkeypatch.setattr(align_features, "build_minimap2_index", fake_build)
    monkeypatch.setattr(align_features, "_run_minimap2_to_sam", fake_align)
    args = SimpleNamespace(
        directory=str(run_dir), reference=str(reference), target=str(target),
        m=None, mm2_options="-a",
    )

    result = align_features.align_single_chroms(
        [str(reference)], [str(target)], 1, args, 8,
        "chrm_by_chrm", 0, target_lengths={"chr1": 8},
    )

    assert Path(result).name == "reference_all_to_target_all.sam"
    assert len(align_calls) == 2
    assert len(build_calls) == 2
    assert build_calls[1]["force_rebuild"] is True
    assert align_calls[1][1][0] == str(rebuilt)


def test_worker_alignment_error_terminates_and_joins_pool(tmp_path, monkeypatch):
    observed = {"close": 0, "terminate": 0, "join": 0}

    class FailingPool:
        def __init__(self, *args, **kwargs):
            pass

        def imap_unordered(self, *args, **kwargs):
            raise LiftOnAlignmentError("synthetic worker failure")

        def close(self):
            observed["close"] += 1

        def terminate(self):
            observed["terminate"] += 1

        def join(self):
            observed["join"] += 1

    monkeypatch.setattr(align_features, "Pool", FailingPool)
    monkeypatch.setattr(
        align_features, "split_target_sequence", lambda *args, **kwargs: {"chr1": "ACGT"}
    )
    args = SimpleNamespace(
        native=False, subcommand=None, target=str(tmp_path / "target.fa"),
        threads=1, directory=str(tmp_path), mm2_options="-a",
    )

    with pytest.raises(LiftOnAlignmentError, match="worker"):
        align_features.align_features_to_target(
            ["reference.fa"], [args.target], args, object(), "chrm_by_chrm", [],
        )

    assert observed == {"close": 0, "terminate": 1, "join": 1}


def test_empty_liftoff_result_fails_before_gff_is_written(tmp_path, monkeypatch):
    args = SimpleNamespace(
        chroms=None,
        reference=str(tmp_path / "reference.fa"),
        target=str(tmp_path / "target.fa"),
        features=None,
        u=str(tmp_path / "unmapped.txt"),
        cds=False,
        polish=False,
        output=str(tmp_path / "output.gff3"),
    )
    feature_db = object()
    feature_hierarchy = object()

    monkeypatch.setattr(
        liftoff_main.liftover_types,
        "lift_original_annotation",
        lambda *a, **k: (feature_db, feature_hierarchy, []),
    )
    monkeypatch.setattr(
        liftoff_main,
        "map_unmapped_features",
        lambda *a, **k: [],
    )
    monkeypatch.setattr(
        liftoff_main, "map_features_from_unplaced_seq", lambda *a, **k: None
    )
    monkeypatch.setattr(
        liftoff_main, "write_unmapped_features_file", lambda *a, **k: None
    )
    monkeypatch.setattr(liftoff_main, "map_extra_copies", lambda *a, **k: None)

    writes = []
    monkeypatch.setattr(
        liftoff_main.write_new_gff,
        "write_new_gff",
        lambda *a, **k: writes.append((a, k)),
    )

    with pytest.raises(LiftOnAlignmentError) as exc:
        liftoff_main.run_all_liftoff_steps(args, ref_db=object())

    message = str(exc.value).lower()
    assert "zero" in message or "no features" in message
    assert "alignment" in message or "lift" in message
    assert writes == []

    with pytest.raises(LiftOnAlignmentError):
        liftoff_main.run_all_liftoff_steps_inmemory(args, ref_db=object())


def test_run_liftoff_reports_alignment_failure_without_annotation_blame(
    tmp_path, monkeypatch, capsys
):
    args = SimpleNamespace(
        polish=False, inmemory_liftoff=False,
        output=str(tmp_path / "output.gff3"),
    )

    def fail_alignment(*args, **kwargs):
        raise LiftOnAlignmentError("SAM contains no @SQ target records")

    monkeypatch.setattr(
        run_liftoff_module.liftoff_main, "run_all_liftoff_steps", fail_alignment
    )
    with pytest.raises(SystemExit) as exc:
        run_liftoff_module.run_liftoff(str(tmp_path) + "/", None, args)

    assert exc.value.code == 1
    stderr = capsys.readouterr().err
    assert "Liftoff alignment failed" in stderr
    assert "no @SQ" in stderr
    assert "ANNOTATION FILE ERROR" not in stderr
