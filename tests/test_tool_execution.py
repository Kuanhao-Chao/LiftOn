"""Resource policy and factual native-process diagnostics for GH #71."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest

from lifton import tool_execution


@pytest.mark.parametrize(
    ("target_bases", "expected"),
    [
        (4_000_000_000, "parallel"),
        (4_000_000_001, "serial"),
    ],
)
def test_default_schedule_boundary(target_bases, expected):
    schedule = tool_execution.resolve_aligner_schedule(
        target_bases, needs_liftoff=True, needs_miniprot=True,
    )
    assert schedule.mode == expected
    assert schedule.parallel is (expected == "parallel")


def test_explicit_schedule_overrides_and_single_tool_cases():
    large_parallel = tool_execution.resolve_aligner_schedule(
        20_000_000_000, needs_liftoff=True, needs_miniprot=True,
        force_parallel=True,
    )
    small_serial = tool_execution.resolve_aligner_schedule(
        10, needs_liftoff=True, needs_miniprot=True, force_serial=True,
    )
    single = tool_execution.resolve_aligner_schedule(
        20_000_000_000, needs_liftoff=False, needs_miniprot=True,
    )
    cached = tool_execution.resolve_aligner_schedule(
        20_000_000_000, needs_liftoff=False, needs_miniprot=False,
    )

    assert (large_parallel.mode, large_parallel.forced) == ("parallel", True)
    assert (small_serial.mode, small_serial.forced) == ("serial", True)
    assert single.mode == "single"
    assert single.active_tools == ("miniprot",)
    assert cached.mode == "cached"


def test_conflicting_schedule_overrides_are_rejected_in_pure_resolver():
    with pytest.raises(ValueError, match="exclusive"):
        tool_execution.resolve_aligner_schedule(
            1, needs_liftoff=True, needs_miniprot=True,
            force_serial=True, force_parallel=True,
        )


def test_parser_rejects_both_schedule_flags():
    from lifton import lifton

    with pytest.raises(SystemExit) as error:
        lifton.parse_args([
            "target.fa", "reference.fa", "-g", "reference.gff3",
            "--serial-aligners", "--parallel-aligners",
        ])
    assert error.value.code == 2


def test_target_statistics_use_index_lengths_without_loading_sequences():
    class ExplodingFasta:
        faidx = SimpleNamespace(index={
            "a": SimpleNamespace(rlen=10),
            "b": SimpleNamespace(rlen=7),
        })

        def __getitem__(self, _name):  # pragma: no cover - assertion path
            raise AssertionError("target sequence was loaded")

    assert tool_execution.target_fasta_statistics(ExplodingFasta()) == {
        "sequence_count": 2,
        "total_bases": 17,
        "maximum_sequence_bases": 10,
    }


@pytest.mark.parametrize(
    ("returncode", "rendered", "signal_name"),
    [
        (-11, "terminated by signal 11 (SIGSEGV)", "SIGSEGV"),
        (17, "exit status 17", None),
        (0, "exit status 0", None),
    ],
)
def test_return_codes_distinguish_signals_from_exit_statuses(
    returncode, rendered, signal_name,
):
    assert tool_execution.describe_returncode(returncode) == rendered
    evidence = tool_execution.returncode_evidence(returncode)
    assert evidence["termination"] == rendered
    assert (evidence["signal"] or {}).get("name") == signal_name
    assert "out of memory" not in str(evidence).lower()


def test_miniprot_stage_matches_issue_71_failure_boundary():
    stderr = """
[M::mp_ntseq_read@1.0*1.0] read 20029007188 bases in 657 contigs
[M::mp_idx_build@2.0*1.0] 156477268 blocks
[M::mp_idx_build@3.0*1.0] collected syncmers
"""
    assert tool_execution.detect_miniprot_stage(stderr) == "syncmers_collected"


@pytest.mark.parametrize(
    ("stderr", "expected"),
    [
        ("[M::mm_idx_gen] collected minimizers\n", "target_minimizers_collected"),
        ("[M::mm_idx_gen] sorted minimizers\n", "target_minimizers_sorted"),
        ("[M::main] loaded/built the index for 657 target sequence(s)\n",
         "target_index_ready"),
        ("[M::worker_pipeline] mapped 74405 sequences\n",
         "query_mapping_complete"),
        ("allocation failed\n", "split-index alignment"),
    ],
)
def test_minimap2_stage_uses_latest_native_milestone(stderr, expected):
    assert tool_execution.detect_minimap2_stage(
        stderr, "split-index alignment",
    ) == expected


def test_execution_record_bounds_stderr_and_preserves_signal():
    record = tool_execution.execution_record(
        "miniprot", command=["miniprot", "target.fa", "proteins.fa"],
        stage="syncmers_collected", status="failed", returncode=-11,
        stderr_tail="x" * (tool_execution.STDERR_TAIL_BYTES + 100),
    )
    assert len(record["stderr_tail"].encode()) == tool_execution.STDERR_TAIL_BYTES
    assert record["signal"]["name"] == "SIGSEGV"


def test_subprocess_runner_drains_and_bounds_stderr():
    result = tool_execution.run_with_bounded_stderr(
        [
            sys.executable, "-c",
            "import sys; sys.stderr.write('ERROR\\n' + 'x' * 70000); "
            "raise SystemExit(17)",
        ],
        tail_limit=4096,
        relay_stderr=False,
    )

    assert result.returncode == 17
    assert result.stderr_tail == "x" * 4096
    assert result.stderr_error_seen is True


@pytest.mark.parametrize(
    ("banner", "expected"),
    [
        ("0.13-r248", (0, 13)),
        ("miniprot 0.18-r300", (0, 18)),
        ("development build", None),
        (None, None),
    ],
)
def test_parse_miniprot_version(banner, expected):
    assert tool_execution.parse_miniprot_version(banner) == expected


@pytest.mark.parametrize(
    ("maximum_bases", "banner", "expected"),
    [
        ((1 << 31) - 1, "0.13-r248", "not_applicable"),
        (1 << 31, "0.13-r248", "incompatible"),
        ((1 << 31) + 1, "0.13-r248", "incompatible"),
        ((1 << 31) + 1, "0.14-r250", "compatible"),
        ((1 << 31) + 1, "development build", "unknown"),
    ],
)
def test_miniprot_long_sequence_compatibility(
    maximum_bases, banner, expected,
):
    assert tool_execution.miniprot_long_sequence_compatibility(
        maximum_bases, banner,
    ) == expected


def test_worker_execution_events_are_collected_into_parent_manifest(tmp_path):
    class RecordingManifest:
        def __init__(self):
            self.records = []

        def record_aligner_execution(self, record):
            self.records.append(record)

    worker_args = SimpleNamespace(aligner_diagnostics_dir=str(tmp_path))
    record = tool_execution.execution_record(
        "minimap2", command=["minimap2", "-d", "target.mmi"],
        stage="index build", status="failed", returncode=-11,
        stderr_tail="indexing\n",
    )
    tool_execution.record_execution(worker_args, record)

    manifest = RecordingManifest()
    parent_args = SimpleNamespace(
        aligner_diagnostics_dir=str(tmp_path), _run_manifest=manifest,
    )
    collected = tool_execution.collect_execution_events(parent_args)

    assert collected == manifest.records
    assert manifest.records[0]["signal"]["name"] == "SIGSEGV"
    assert not list(tmp_path.glob("execution-*.json"))
