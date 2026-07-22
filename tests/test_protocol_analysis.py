from __future__ import annotations

import copy

import pytest

from benchmarks.compare import protocol_analysis


def _records(kind: str):
    records = []
    for case in protocol_analysis.protocol_cases(kind):
        row = dict(case)
        threads = case["threads"]
        row.update({
            "candidate_sha": "1" * 40,
            "profile": {
                "wall_clock_seconds": 100.0 / threads,
                "peak_rss_mb": 500.0 + threads,
                "user_cpu_seconds": 90.0,
                "sys_cpu_seconds": 5.0,
                "filesystem_inputs": 10,
                "filesystem_outputs": 20,
                "major_page_faults": 0,
                "minor_page_faults": 30,
            },
            "fingerprints": {
                "byte_sha256": "2" * 64,
                "semantic_sha256": "3" * 64,
            },
        })
        records.append(row)
    return records


def test_thread_schedule_is_complete_and_balanced():
    cases = protocol_analysis.thread_scaling_cases()
    assert len(cases) == 36
    for repetition in range(1, 7):
        assert {
            case["threads"] for case in cases
            if case["repetition"] == repetition
        } == set(protocol_analysis.THREAD_LEVELS)
    for position in range(1, 7):
        assert {
            case["threads"] for case in cases
            if case["position"] == position
        } == set(protocol_analysis.THREAD_LEVELS)


def test_io_schedule_has_four_copies_of_each_mode():
    cases = protocol_analysis.io_mode_cases()
    assert len(cases) == 32
    assert {
        mode: sum(case["mode"] == mode for case in cases)
        for mode in protocol_analysis.IO_MODE_ORDER
    } == {mode: 4 for mode in protocol_analysis.IO_MODE_ORDER}


def test_thread_summary_reports_speedup_and_efficiency():
    result = protocol_analysis.summarize_protocol(
        _records("thread_scaling"),
        kind="thread_scaling",
    )
    assert result["n_arms"] == 36
    assert result["outputs_identical"] is True
    by_threads = {row["threads"]: row for row in result["summaries"]}
    assert by_threads[1]["speedup"] == pytest.approx(1.0)
    assert by_threads[8]["speedup"] == pytest.approx(8.0)
    assert by_threads[8]["parallel_efficiency"] == pytest.approx(1.0)
    assert by_threads[8]["median_filesystem_outputs"] == 20.0


def test_io_summary_preserves_mode_order():
    result = protocol_analysis.summarize_protocol(
        _records("io_modes"),
        kind="io_modes",
    )
    assert [row["mode"] for row in result["summaries"]] == list(
        protocol_analysis.IO_MODE_ORDER
    )


def test_summary_rejects_missing_case():
    records = _records("thread_scaling")
    with pytest.raises(ValueError, match="incomplete"):
        protocol_analysis.summarize_protocol(
            records[:-1],
            kind="thread_scaling",
        )


def test_summary_rejects_output_drift():
    records = _records("io_modes")
    changed = copy.deepcopy(records)
    changed[-1]["fingerprints"]["byte_sha256"] = "4" * 64
    with pytest.raises(ValueError, match="not byte- and semantic-identical"):
        protocol_analysis.summarize_protocol(changed, kind="io_modes")


def test_summary_rejects_malformed_optional_counter():
    records = _records("io_modes")
    records[0]["profile"]["filesystem_inputs"] = -1
    with pytest.raises(ValueError, match="non-negative integer or null"):
        protocol_analysis.summarize_protocol(records, kind="io_modes")
