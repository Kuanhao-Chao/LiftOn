"""Re-export the audited GNU/BSD ``/usr/bin/time`` profiling helpers from the
Phase-16 harness (``benchmarks/run_benchmarks.py``) so the comparison pipeline
reuses them verbatim rather than reimplementing wall-clock / peak-RSS parsing.
"""
import sys
from pathlib import Path

# Make ``benchmarks`` importable regardless of cwd (repo root on sys.path).
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from benchmarks.run_benchmarks import (  # noqa: E402
    run_profiled,
    ProfileResult,
    _parse_gnu_time,
    _wall_clock_str_to_seconds,
)

__all__ = [
    "run_profiled",
    "ProfileResult",
    "_parse_gnu_time",
    "_wall_clock_str_to_seconds",
]
