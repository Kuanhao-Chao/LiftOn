"""Unit tests for the Phase 16 benchmark-harness time-log parsers.

The original `_parse_gnu_time` cut on the first colon in the matched
prefix, which silently corrupted the wall-clock value because GNU
`/usr/bin/time -v` emits the elapsed-time line as
`"Elapsed (wall clock) time (h:mm:ss or m:ss): 2:19.75"` — the colon
inside the prefix was being treated as the value separator. This file
guards against that regression and against the call-site fragility that
let a single parse failure collapse the entire benchmark report.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# `benchmarks/run_benchmarks.py` is a script, not part of the `lifton`
# package, so we put the benchmarks dir on sys.path and import by name.
_HARNESS_DIR = Path(__file__).resolve().parent.parent / "benchmarks"
sys.path.insert(0, str(_HARNESS_DIR))

import run_benchmarks as harness  # noqa: E402


GNU_TIME_SAMPLE = """\
Command exited with non-zero status 1
\tCommand being timed: "lifton --stream ..."
\tUser time (seconds): 127.05
\tSystem time (seconds): 6.47
\tPercent of CPU this job got: 95%
\tElapsed (wall clock) time (h:mm:ss or m:ss): 2:19.75
\tAverage shared text size (kbytes): 0
\tMaximum resident set size (kbytes): 5472328
\tExit status: 1
"""


class TestParseGnuTime:
    def test_extracts_wall_clock_str_with_colon_in_prefix(self):
        parsed = harness._parse_gnu_time(GNU_TIME_SAMPLE)
        # Pre-fix this used to be "mm:ss or m:ss): 2:19.75".
        assert parsed["wall_clock_str"] == "2:19.75"

    def test_extracts_user_cpu_seconds(self):
        parsed = harness._parse_gnu_time(GNU_TIME_SAMPLE)
        assert parsed["user_cpu_seconds"] == "127.05"

    def test_extracts_sys_cpu_seconds(self):
        parsed = harness._parse_gnu_time(GNU_TIME_SAMPLE)
        assert parsed["sys_cpu_seconds"] == "6.47"

    def test_extracts_max_rss_kb(self):
        parsed = harness._parse_gnu_time(GNU_TIME_SAMPLE)
        assert parsed["max_rss_kb"] == "5472328"

    def test_real_fixture_from_bee_run(self):
        # The bee/logs/lift.time.log on disk is a *runtime artefact* that
        # gets overwritten — and sometimes truncated to zero bytes — by
        # every benchmark invocation. Skip if the fixture is missing or
        # empty (e.g., a previous run was killed before the GNU time
        # finalised the file). When present and non-empty, assert the
        # parser produces clean colon-separated numeric values (the
        # pre-fix bug produced the corruption "mm:ss or m:ss): 2:19.75"
        # instead).
        bee_log = (Path(__file__).resolve().parent.parent
                   / "benchmarks" / "results" / "bee" / "logs"
                   / "lift.time.log")
        if not bee_log.exists() or bee_log.stat().st_size == 0:
            return
        parsed = harness._parse_gnu_time(bee_log.read_text())
        if "wall_clock_str" not in parsed:
            # Time-log was written but doesn't contain the
            # "Elapsed (wall clock) time" line — partial / truncated.
            return
        # Acceptable shapes: "m:ss[.ss]", "mm:ss[.ss]", "h:mm:ss[.ss]".
        assert re.match(
            r"^\d+(?::\d+){1,2}(?:\.\d+)?$", parsed["wall_clock_str"],
        ), f"unexpected wall_clock_str shape: {parsed['wall_clock_str']!r}"
        # max_rss_kb is a plain integer; the parser fix didn't touch
        # this key but it's still a useful invariant.
        assert parsed.get("max_rss_kb", "").isdigit(), parsed


class TestWallClockStrToSeconds:
    def test_two_part_m_ss(self):
        assert harness._wall_clock_str_to_seconds("2:19.75") == 139.75

    def test_two_part_mm_ss(self):
        # human dataset shape (21 minutes 06.51 seconds).
        assert abs(harness._wall_clock_str_to_seconds("21:06.51")
                   - (21 * 60 + 6.51)) < 1e-9

    def test_three_part_h_mm_ss(self):
        assert abs(harness._wall_clock_str_to_seconds("1:23:45.6")
                   - (3600 + 23 * 60 + 45.6)) < 1e-9


class TestSafeFloat:
    def test_returns_default_on_garbage_string(self):
        # The exact pre-fix corruption that crashed the harness.
        assert harness._safe_float("mm:ss or m:ss): 21:06.51", 0.0) == 0.0

    def test_returns_default_on_none(self):
        assert harness._safe_float(None, 7.5) == 7.5

    def test_returns_default_on_empty_string(self):
        assert harness._safe_float("", 1.25) == 1.25

    def test_passes_through_valid_float_string(self):
        assert harness._safe_float("127.05", 0.0) == 127.05

    def test_passes_through_real_float(self):
        assert harness._safe_float(3.14, 0.0) == 3.14
