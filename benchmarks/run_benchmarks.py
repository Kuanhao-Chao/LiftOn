#!/usr/bin/env python3
"""Phase 16 — biological validation benchmark harness for LiftOn.

A single Python driver that:
  1. Downloads (resumably) the FASTA + GFF inputs from the JHU CCB
     mirror for any subset of the five published same-species datasets.
  2. Runs the full LiftOn pipeline with the Phase 11/15 fastest
     parallel configuration (`--stream --inmemory-liftoff
     --locus-pipeline -t 8 --native`).
  3. Captures Peak RSS + wall-clock via the platform's `time` utility
     (auto-detects GNU `/usr/bin/time -v` on Linux, BSD `/usr/bin/time
     -l` on macOS) plus an in-process `resource.getrusage` fallback.
  4. Optionally re-invokes LiftOn in `--evaluation` (-E) mode on the newly
     generated candidate annotation to extract identity metrics.  A published
     target GFF remains an optional comparison asset; it is never substituted
     for the candidate being measured.
  5. Writes per-dataset JSON + a roll-up summary table to stdout.

Usage examples:

    # Run everything (download + lift + evaluate) for all 5 datasets
    python benchmarks/run_benchmarks.py --all

    # Only download (useful on a Slurm login node)
    python benchmarks/run_benchmarks.py --download-only --datasets bee arabidopsis

    # Run a single dataset and emit JSON to a known location
    python benchmarks/run_benchmarks.py --datasets bee \\
        --output benchmarks/results/run_$(date +%Y%m%d).json

    # Skip evaluation (lift-only, fastest)
    python benchmarks/run_benchmarks.py --datasets human --no-evaluation

The script is hermetic about state: every dataset writes under
`benchmarks/data/<id>/` (inputs) and `benchmarks/results/<id>/`
(outputs + logs). Re-running the same dataset is idempotent: if
inputs already exist on disk and pass a sentinel size check, they
are not redownloaded; if a `.lifton.done` flag already exists, the
LiftOn invocation is skipped (use `--force` to override).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata as importlib_metadata
import importlib.util
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
DEFAULT_REGISTRY = HERE / "datasets.json"
DEFAULT_DATA_DIR = HERE / "data"
DEFAULT_RESULTS_DIR = HERE / "results"
LIFT_CACHE_DISTRIBUTIONS = (
    "duckdb",
    "gffutils",
    "mappy",
    "numpy",
    "parasail",
    "pyarrow",
    "pyfaidx",
    "pysam",
)


# ---------------------------------------------------------------------------
# Dataset registry
# ---------------------------------------------------------------------------

@dataclass
class Dataset:
    id: str
    species: str
    reference_fa: str
    target_fa: str
    reference_gff: str
    target_gff: Optional[str] = None
    approx_size_gb: float = 0.0
    cross_species: bool = False
    annotation_database: str = "RefSeq"
    truth_gff: Optional[str] = None
    ortholog_map: Optional[str] = None
    truth_id_policy: str = "ortholog-map"


@dataclass
class Registry:
    datasets: list[Dataset]
    lifton_flags: list[str] = field(default_factory=list)
    evaluation_flags: list[str] = field(default_factory=list)


def load_registry(path: Path) -> Registry:
    raw = json.loads(path.read_text())
    # Strip free-form comment keys (anything starting with "_") both from
    # comment-only entries and from real dataset entries.  Several registry
    # rows carry an inline ``_note``; passing it to Dataset(**entry) used to
    # raise TypeError outside the build-controller's private sanitizer.
    ds_entries = []
    for entry in raw["datasets"]:
        if not isinstance(entry, dict):
            continue
        sanitized = {
            key: value for key, value in entry.items()
            if not str(key).startswith("_")
        }
        if sanitized:
            ds_entries.append(sanitized)
    datasets = [Dataset(**entry) for entry in ds_entries]
    return Registry(
        datasets=datasets,
        lifton_flags=list(raw.get("lifton_flags", [])),
        evaluation_flags=list(raw.get("evaluation_flags", [])),
    )


# ---------------------------------------------------------------------------
# Resilient HTTP/FTP download
# ---------------------------------------------------------------------------

def _filename_for(url: str) -> str:
    return url.rsplit("/", 1)[-1]


def download(url: str, dest: Path, *, min_bytes: int = 1024,
             retries: int = 3, log=print) -> Path:
    """Download ``url`` into ``dest``. Idempotent: if the file already
    exists and is at least ``min_bytes`` in size, skip. Retries on
    failure with exponential backoff."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        log(f"  ✓ already on disk: {dest.name} "
            f"({dest.stat().st_size / 1e6:.1f} MB)")
        return dest
    last_err: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        log(f"  ↓ {url}  (attempt {attempt}/{retries})")
        try:
            tmp = dest.with_suffix(dest.suffix + ".part")
            with urllib.request.urlopen(url, timeout=120) as resp, \
                    open(tmp, "wb") as fh:
                shutil.copyfileobj(resp, fh, length=1 << 20)
            os.replace(tmp, dest)
            return dest
        except Exception as exc:
            last_err = exc
            log(f"    ! download failed: {exc}")
            time.sleep(2 ** attempt)
    raise RuntimeError(
        f"download failed after {retries} attempts: {url} ({last_err})"
    )


# ---------------------------------------------------------------------------
# Profiling — peak RSS + wall-clock
# ---------------------------------------------------------------------------

@dataclass
class ProfileResult:
    wall_clock_seconds: float
    peak_rss_mb: float
    user_cpu_seconds: float
    sys_cpu_seconds: float
    exit_code: int
    stdout_path: str
    stderr_path: str
    time_log_path: Optional[str]
    filesystem_inputs: Optional[int] = None
    filesystem_outputs: Optional[int] = None
    major_page_faults: Optional[int] = None
    minor_page_faults: Optional[int] = None
    voluntary_context_switches: Optional[int] = None
    involuntary_context_switches: Optional[int] = None


def _platform_time_argv() -> tuple[Optional[list[str]], str]:
    """Return ``(argv_prefix, kind)`` where:
      * ``argv_prefix`` is the list of args to prepend to the actual
        command (e.g. ``["/usr/bin/time", "-v"]`` on Linux,
        ``["/usr/bin/time", "-l"]`` on macOS), or None when no system
        ``time`` is suitable.
      * ``kind`` is one of "gnu", "bsd", "rusage" — used by the
        log parser to know which output format to expect.
    """
    p = "/usr/bin/time"
    if not Path(p).exists():
        return None, "rusage"
    sysname = platform.system()
    if sysname == "Linux":
        return [p, "-v"], "gnu"
    if sysname == "Darwin":
        return [p, "-l"], "bsd"
    return None, "rusage"


_GNU_KEYS = {
    "Maximum resident set size (kbytes)": "max_rss_kb",
    "User time (seconds)": "user_cpu_seconds",
    "System time (seconds)": "sys_cpu_seconds",
    "Elapsed (wall clock) time (h:mm:ss or m:ss)": "wall_clock_str",
    "File system inputs": "filesystem_inputs",
    "File system outputs": "filesystem_outputs",
    "Major (requiring I/O) page faults": "major_page_faults",
    "Minor (reclaiming a frame) page faults": "minor_page_faults",
    "Voluntary context switches": "voluntary_context_switches",
    "Involuntary context switches": "involuntary_context_switches",
}


def _parse_gnu_time(log: str) -> dict:
    out: dict[str, Any] = {}
    for line in log.splitlines():
        line = line.strip()
        for prefix, key in _GNU_KEYS.items():
            if line.startswith(prefix):
                # Slice past the matched prefix rather than splitting on
                # the first colon — the wall-clock prefix itself
                # contains colons (`(h:mm:ss or m:ss)`), so a naive
                # split(":", 1) would cut inside the prefix and capture
                # garbage like "mm:ss or m:ss): 2:19.75" instead of the
                # numeric value.
                out[key] = line[len(prefix):].lstrip(": ").strip()
                break
    return out


def _safe_float(x: Any, default: float) -> float:
    """Convert ``x`` to float, returning ``default`` on failure or when
    ``x`` is None / empty. Used at the time-log parsing call site so
    one malformed field never collapses the whole benchmark report."""
    if x is None or x == "":
        return default
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def _optional_int(value: Any) -> Optional[int]:
    """Parse an optional integer emitted by ``time``.

    Resource counters are absent on unsupported platforms and in historical
    result JSON.  Returning ``None`` preserves that distinction instead of
    incorrectly reporting an unavailable counter as zero.
    """

    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _parse_bsd_time(log: str) -> dict:
    """macOS BSD `/usr/bin/time -l` puts everything in one trailing
    block. Parse the lines we care about by suffix."""
    out: dict[str, Any] = {}
    for line in log.splitlines():
        m = re.match(r"\s*(\d+(?:\.\d+)?)\s+(.*)$", line)
        if not m:
            continue
        val, label = m.group(1), m.group(2).strip()
        if "maximum resident set size" in label:
            # macOS reports bytes, not kbytes.
            out["max_rss_kb"] = float(val) / 1024.0
        elif label == "real":
            out["wall_clock_seconds"] = float(val)
        elif label == "user":
            out["user_cpu_seconds"] = float(val)
        elif label == "sys":
            out["sys_cpu_seconds"] = float(val)
    return out


def _wall_clock_str_to_seconds(s: str) -> float:
    """`/usr/bin/time -v` prints "h:mm:ss" or "m:ss.ss"."""
    parts = s.split(":")
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    return float(s)


def run_profiled(argv: list[str], *, label: str, log_dir: Path,
                 env: Optional[dict] = None,
                 cwd: Optional[Path] = None,
                 log=print) -> ProfileResult:
    """Run ``argv`` under the platform `time` utility and capture
    stdout/stderr + a peak-RSS report. Falls back to in-process
    ``resource.getrusage`` when no system time is available."""
    log_dir.mkdir(parents=True, exist_ok=True)
    out_path = log_dir / f"{label}.stdout.log"
    err_path = log_dir / f"{label}.stderr.log"
    time_path = log_dir / f"{label}.time.log"

    prefix, kind = _platform_time_argv()
    full_argv: list[str]
    if prefix is None:
        full_argv = list(argv)
    else:
        # /usr/bin/time writes its summary to stderr in both GNU and BSD
        # modes; we tee stderr -> file AND -> the time log via the
        # `-o` switch.
        full_argv = list(prefix) + ["-o", str(time_path)] + list(argv)
        kind_log = kind

    log(f"\n[bench] {label} — invoking:")
    log(f"  {' '.join(argv)}")
    t0 = time.time()
    try:
        with open(out_path, "wb") as out_fh, open(err_path, "wb") as err_fh:
            proc = subprocess.run(
                full_argv,
                stdout=out_fh, stderr=err_fh,
                env={**os.environ, **(env or {})},
                cwd=str(cwd) if cwd else None,
                check=False,
            )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"could not exec {full_argv[0]!r}: {exc}"
        ) from exc
    wall = time.time() - t0
    rc = proc.returncode

    parsed: dict[str, Any] = {}
    user = 0.0
    sys_t = 0.0
    rss_kb = 0.0
    if prefix is not None and time_path.exists():
        log_text = time_path.read_text(errors="replace")
        if kind == "gnu":
            parsed = _parse_gnu_time(log_text)
            if "wall_clock_str" in parsed:
                try:
                    wall = _wall_clock_str_to_seconds(parsed["wall_clock_str"])
                except (ValueError, IndexError) as exc:
                    log(f"  ! could not parse wall-clock from "
                        f"{time_path} ({exc!s}); falling back to "
                        f"in-process timer ({wall:.2f}s).")
            user = _safe_float(parsed.get("user_cpu_seconds"), 0.0)
            sys_t = _safe_float(parsed.get("sys_cpu_seconds"), 0.0)
            rss_kb = _safe_float(parsed.get("max_rss_kb"), 0.0)
        else:
            parsed = _parse_bsd_time(log_text)
            wall = _safe_float(parsed.get("wall_clock_seconds"), wall)
            user = _safe_float(parsed.get("user_cpu_seconds"), 0.0)
            sys_t = _safe_float(parsed.get("sys_cpu_seconds"), 0.0)
            rss_kb = _safe_float(parsed.get("max_rss_kb"), 0.0)
    else:
        # Fallback to in-process rusage. Note: only sees current proc,
        # NOT the LiftOn child — so it's a lower bound. Use only when
        # /usr/bin/time is unavailable.
        import resource
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        rss_kb = float(usage.ru_maxrss)
        if platform.system() == "Darwin":
            rss_kb /= 1024.0
        user = usage.ru_utime
        sys_t = usage.ru_stime
        time_path = None  # type: ignore[assignment]

    return ProfileResult(
        wall_clock_seconds=wall,
        peak_rss_mb=rss_kb / 1024.0,
        user_cpu_seconds=user,
        sys_cpu_seconds=sys_t,
        exit_code=rc,
        stdout_path=str(out_path),
        stderr_path=str(err_path),
        time_log_path=str(time_path) if time_path else None,
        filesystem_inputs=_optional_int(parsed.get("filesystem_inputs")),
        filesystem_outputs=_optional_int(parsed.get("filesystem_outputs")),
        major_page_faults=_optional_int(parsed.get("major_page_faults")),
        minor_page_faults=_optional_int(parsed.get("minor_page_faults")),
        voluntary_context_switches=_optional_int(
            parsed.get("voluntary_context_switches")
        ),
        involuntary_context_switches=_optional_int(
            parsed.get("involuntary_context_switches")
        ),
    )


# ---------------------------------------------------------------------------
# Evaluation summary parsing
# ---------------------------------------------------------------------------

@dataclass
class EvalSummary:
    mapped: Optional[int] = None
    lost: Optional[int] = None
    extra_copies: Optional[int] = None
    avg_identity: Optional[float] = None
    score_file: Optional[str] = None
    format: Optional[str] = None
    records: int = 0
    coding_records: int = 0
    noncoding_records: int = 0
    malformed_records: int = 0


_NONCODING_STATUSES = frozenset({"non_coding"})


def _identity(value: str) -> Optional[float]:
    """Parse a finite identity in the closed interval [0, 1]."""
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(parsed) or parsed < 0.0 or parsed > 1.0:
        return None
    return parsed


def _is_noncoding_status(value: str) -> bool:
    tokens = {
        token.strip().lower()
        for token in str(value).split(";")
        if token.strip()
    }
    return bool(tokens & _NONCODING_STATUSES)


def parse_score_txt(path: Path) -> EvalSummary:
    """Parse the LiftOn score.txt file into mapped/lost/identity stats.

    The eight-column format is emitted by
    ``lifton_utils.write_lifton_status`` — one row per emitted transcript:

        transcript_id, liftoff_PI, miniprot_PI, LiftOn_DNA_identity,
        LiftOn_protein_identity, annotation_source, status, target_interval

    Older harness code documented a different seven-column layout and read the
    protein identity and status from the wrong columns.  It consequently
    reported ``mapped=0`` and ``avg_identity=None`` for valid production runs.
    """
    if not path.exists():
        return EvalSummary(score_file=str(path), format="lifton_score_v1")
    records = 0
    coding_records = 0
    noncoding_records = 0
    malformed_records = 0
    identities: list[float] = []
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            if len(parts) != 8:
                malformed_records += 1
                continue
            aa_identity = _identity(parts[4])
            if aa_identity is None:
                malformed_records += 1
                continue
            records += 1
            if _is_noncoding_status(parts[6]):
                noncoding_records += 1
            else:
                coding_records += 1
                identities.append(aa_identity)
    avg_id = sum(identities) / len(identities) if identities else None
    return EvalSummary(
        mapped=records,
        lost=None,
        avg_identity=avg_id,
        score_file=str(path),
        format="lifton_score_v1",
        records=records,
        coding_records=coding_records,
        noncoding_records=noncoding_records,
        malformed_records=malformed_records,
    )


def parse_eval_txt(path: Path) -> EvalSummary:
    """Parse LiftOn's five-column evaluation-mode output.

    ``lifton_utils.write_lifton_eval_status`` emits:

        transcript_id, DNA_identity, protein_identity, status, target_interval

    Evaluation output contains evaluated target records, not a reference-left
    join, so it does not encode a meaningful ``lost`` count.  Feature-level
    mapped/lost counts are returned separately in ``biological_summary``.
    """
    if not path.exists():
        return EvalSummary(score_file=str(path), format="lifton_eval_v1")
    records = 0
    coding_records = 0
    noncoding_records = 0
    malformed_records = 0
    identities: list[float] = []
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            if len(parts) != 5:
                malformed_records += 1
                continue
            dna_identity = _identity(parts[1])
            protein_identity = _identity(parts[2])
            if dna_identity is None or protein_identity is None:
                malformed_records += 1
                continue
            records += 1
            if _is_noncoding_status(parts[3]):
                noncoding_records += 1
            else:
                coding_records += 1
                identities.append(protein_identity)
    avg_id = sum(identities) / len(identities) if identities else None
    return EvalSummary(
        mapped=records,
        lost=None,
        avg_identity=avg_id,
        score_file=str(path),
        format="lifton_eval_v1",
        records=records,
        coding_records=coding_records,
        noncoding_records=noncoding_records,
        malformed_records=malformed_records,
    )


def count_feature_stats(stats_dir: Path) -> dict[str, int]:
    """LiftOn writes plaintext lists at:
        stats/mapped_feature.txt      — one reference feature per line
        stats/unmapped_features.txt   — one feature id per line
        stats/extra_copy_features.txt — same shape

    These are all feature-level quantities, so mapped and lost retain a common
    denominator instead of being mixed with transcript rows from score.txt.
    """
    def _count(name: str) -> int:
        p = stats_dir / name
        if not p.exists():
            return 0
        with p.open() as handle:
            return sum(1 for line in handle if line.strip())

    return {
        "mapped_features": _count("mapped_feature.txt"),
        "lost_features": _count("unmapped_features.txt"),
        "extra_copy_features": _count("extra_copy_features.txt"),
        "mapped_transcripts": _count("mapped_transcript.txt"),
    }


def count_unmapped_extra(stats_dir: Path) -> tuple[int, int]:
    """Backward-compatible two-count wrapper used by older callers."""
    feature_stats = count_feature_stats(stats_dir)
    return (
        feature_stats["lost_features"],
        feature_stats["extra_copy_features"],
    )


def validate_biological_result(
    result: dict[str, Any],
    *,
    require_evaluation: bool = False,
    require_identity: bool = False,
) -> list[str]:
    """Validate the strict E2E biology contract without trusting legacy totals."""

    errors: list[str] = []
    biological = result.get("biological_summary")
    if not isinstance(biological, dict):
        return ["biological_summary is missing or not an object"]
    if biological.get("schema_version") != 1:
        errors.append("biological_summary schema_version is not 1")
    integer_fields = [
        "reference_features",
        "mapped_features",
        "lost_features",
        "extra_copy_features",
        "emitted_transcript_records",
        "mapped_transcripts_reported",
    ]
    if require_evaluation:
        integer_fields.extend([
            "evaluated_transcript_records",
            "evaluated_coding_records",
        ])
    for field_name in integer_fields:
        value = biological.get(field_name)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            errors.append(
                f"biological_summary.{field_name} is not a non-negative integer"
            )

    reference_features = biological.get("reference_features")
    mapped_features = biological.get("mapped_features")
    lost_features = biological.get("lost_features")
    if isinstance(reference_features, int) and not isinstance(reference_features, bool):
        if reference_features <= 0:
            errors.append("biological_summary has no reference features")
        if (
            isinstance(mapped_features, int)
            and not isinstance(mapped_features, bool)
            and isinstance(lost_features, int)
            and not isinstance(lost_features, bool)
            and mapped_features + lost_features != reference_features
        ):
            errors.append(
                "biological_summary mapped_features + lost_features does not "
                "equal reference_features"
            )
        completeness = biological.get("feature_completeness")
        if (
            not isinstance(completeness, (int, float))
            or isinstance(completeness, bool)
            or not math.isfinite(float(completeness))
            or not 0.0 <= float(completeness) <= 1.0
        ):
            errors.append(
                "biological_summary.feature_completeness is not finite in [0, 1]"
            )
        elif (
            reference_features > 0
            and isinstance(mapped_features, int)
            and not math.isclose(
                float(completeness),
                mapped_features / reference_features,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        ):
            errors.append("biological_summary feature completeness is inconsistent")

    emitted = biological.get("emitted_transcript_records")
    if isinstance(emitted, int) and not isinstance(emitted, bool) and emitted <= 0:
        errors.append("biological_summary has no emitted transcript records")
    if require_evaluation:
        evaluated = biological.get("evaluated_transcript_records")
        evaluated_coding = biological.get("evaluated_coding_records")
        if (
            isinstance(evaluated, int)
            and not isinstance(evaluated, bool)
            and evaluated <= 0
        ):
            errors.append("biological_summary has no evaluated transcript records")
        if (
            isinstance(evaluated, int)
            and not isinstance(evaluated, bool)
            and isinstance(evaluated_coding, int)
            and not isinstance(evaluated_coding, bool)
            and evaluated_coding > evaluated
        ):
            errors.append(
                "biological_summary.evaluated_coding_records exceeds "
                "evaluated_transcript_records"
            )

    summaries = [("score_summary", "lifton_score_v1")]
    if require_evaluation:
        summaries.append(("evaluation_summary", "lifton_eval_v1"))
    for field_name, expected_format in summaries:
        summary = result.get(field_name)
        if not isinstance(summary, dict):
            errors.append(f"{field_name} is missing or not an object")
            continue
        if summary.get("format") != expected_format:
            errors.append(f"{field_name}.format is not {expected_format}")
        records = summary.get("records")
        coding_records = summary.get("coding_records")
        noncoding_records = summary.get("noncoding_records")
        malformed = summary.get("malformed_records")
        if not isinstance(records, int) or isinstance(records, bool) or records <= 0:
            errors.append(f"{field_name}.records is not a positive integer")
        for count_name, count in (
            ("coding_records", coding_records),
            ("noncoding_records", noncoding_records),
        ):
            if (
                not isinstance(count, int)
                or isinstance(count, bool)
                or count < 0
            ):
                errors.append(
                    f"{field_name}.{count_name} is not a non-negative integer"
                )
        if (
            isinstance(records, int)
            and isinstance(coding_records, int)
            and isinstance(noncoding_records, int)
            and coding_records + noncoding_records != records
        ):
            errors.append(
                f"{field_name} coding + noncoding records does not equal records"
            )
        if malformed != 0:
            errors.append(f"{field_name}.malformed_records is not zero")
        corresponding = (
            "evaluated_transcript_records"
            if field_name == "evaluation_summary"
            else "emitted_transcript_records"
        )
        if (
            isinstance(records, int)
            and records != biological.get(corresponding)
        ):
            errors.append(
                f"{field_name}.records disagrees with "
                f"biological_summary.{corresponding}"
            )
        if (
            field_name == "evaluation_summary"
            and isinstance(coding_records, int)
            and coding_records != biological.get("evaluated_coding_records")
        ):
            errors.append(
                "evaluation_summary.coding_records disagrees with "
                "biological_summary.evaluated_coding_records"
            )

    identity = biological.get("mean_protein_identity")
    score = result.get("score_summary")
    evaluation = result.get("evaluation_summary")
    coding_records = (
        (evaluation or {}).get("coding_records")
        if require_evaluation
        else (score or {}).get("coding_records")
    )
    identity_required = require_identity or (
        isinstance(coding_records, int) and coding_records > 0
    )
    if identity_required and (
        not isinstance(identity, (int, float))
        or isinstance(identity, bool)
        or not math.isfinite(float(identity))
        or not 0.0 <= float(identity) <= 1.0
    ):
        errors.append(
            "biological_summary.mean_protein_identity is not finite in [0, 1]"
        )
    return errors


# ---------------------------------------------------------------------------
# Per-dataset orchestration
# ---------------------------------------------------------------------------

def _which(name: str) -> Optional[str]:
    return shutil.which(name)


def _ensure_runtime() -> dict:
    """Verify the host has the binaries LiftOn needs at runtime."""
    info = {}
    for tool in ("lifton", "minimap2", "miniprot"):
        info[tool] = _which(tool) or None
    info["python"] = sys.executable
    info["platform"] = platform.platform()
    return info


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_fingerprint(path: Path) -> dict[str, Any]:
    resolved = Path(path).resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size": stat.st_size,
        "sha256": _sha256_file(resolved),
    }


def _lift_tool_provenance() -> dict[str, Any]:
    """Fingerprint the executable stack and importable LiftOn source tree."""

    tools = {}
    for name in ("lifton", "minimap2", "miniprot"):
        resolved = _which(name)
        tools[name] = (
            _file_fingerprint(Path(resolved))
            if resolved is not None else {"path": None}
        )
    tools["python"] = {
        **_file_fingerprint(Path(sys.executable)),
        "version": platform.python_version(),
    }
    spec = importlib.util.find_spec("lifton")
    origin = (
        Path(spec.origin).resolve()
        if spec is not None and spec.origin is not None else None
    )
    source_files = []
    if origin is not None:
        source_root = origin.parent
        for path in sorted(source_root.rglob("*.py")):
            record = _file_fingerprint(path)
            record["relative_path"] = str(path.relative_to(source_root))
            source_files.append(record)
        source_payload = json.dumps(
            source_files,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
        source = {
            "root": str(source_root),
            "n_files": len(source_files),
            "sha256": hashlib.sha256(source_payload).hexdigest(),
        }
    else:
        source = {"root": None, "n_files": 0, "sha256": None}
    distributions = {}
    for requested in LIFT_CACHE_DISTRIBUTIONS:
        try:
            distribution = importlib_metadata.distribution(requested)
        except importlib_metadata.PackageNotFoundError as exc:
            raise RuntimeError(
                f"required LiftOn distribution is unavailable: {requested}"
            ) from exc
        metadata = {}
        for filename in ("METADATA", "WHEEL", "INSTALLER", "direct_url.json"):
            try:
                text = distribution.read_text(filename)
            except (OSError, UnicodeError) as exc:
                raise RuntimeError(
                    f"cannot fingerprint {requested}/{filename}: {exc}"
                ) from exc
            if text is None:
                continue
            encoded = text.encode()
            metadata[filename] = {
                "size": len(encoded),
                "sha256": hashlib.sha256(encoded).hexdigest(),
            }
        version = str(distribution.version or "").strip()
        if not version:
            raise RuntimeError(
                f"required LiftOn distribution has no version: {requested}"
            )
        distributions[requested] = {
            "name": str(distribution.metadata.get("Name") or requested),
            "version": version,
            "metadata": metadata,
        }
    return {
        "tools": tools,
        "lifton_source": source,
        "distributions": distributions,
        "environment": {
            name: value
            for name, value in sorted(os.environ.items())
            if name.startswith("LIFTON_") or name == "PYTHONPATH"
        },
    }


def _lift_cache_request(
    argv: list[str],
    paths: dict[str, Path],
) -> dict[str, Any]:
    inputs = {
        name: _file_fingerprint(paths[name])
        for name in ("reference_fa", "target_fa", "reference_gff")
    }
    return {
        "argv": list(argv),
        "inputs": inputs,
        "provenance": _lift_tool_provenance(),
    }


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _valid_lift_profile(profile: ProfileResult) -> bool:
    numeric_profile = (
        profile.wall_clock_seconds,
        profile.peak_rss_mb,
        profile.user_cpu_seconds,
        profile.sys_cpu_seconds,
    )
    return (
        profile.exit_code == 0
        and all(
            isinstance(value, (int, float))
            and not isinstance(value, bool)
            and math.isfinite(float(value))
            for value in numeric_profile
        )
        and float(profile.wall_clock_seconds) > 0
        and float(profile.peak_rss_mb) > 0
        and float(profile.user_cpu_seconds) >= 0
        and float(profile.sys_cpu_seconds) >= 0
    )


def _cached_lift_profile(
    done_flag: Path,
    *,
    request: dict[str, Any],
    out_gff: Path,
) -> Optional[ProfileResult]:
    try:
        document = json.loads(done_flag.read_text())
    except (OSError, TypeError, ValueError):
        return None
    if (
        not isinstance(document, dict)
        or document.get("schema_version") != 1
        or document.get("kind") != "lifton_lift_cache"
        or document.get("request") != request
        or document.get("key") != _canonical_sha256(request)
    ):
        return None
    artifact = document.get("artifact")
    try:
        observed = _file_fingerprint(out_gff)
    except OSError:
        return None
    if artifact != observed:
        return None
    raw_profile = document.get("profile")
    if not isinstance(raw_profile, dict):
        return None
    try:
        profile = ProfileResult(**raw_profile)
    except (TypeError, ValueError):
        return None
    if not _valid_lift_profile(profile):
        return None
    return profile


def _write_lift_cache(
    done_flag: Path,
    *,
    request: dict[str, Any],
    profile: ProfileResult,
    out_gff: Path,
) -> None:
    if not _valid_lift_profile(profile):
        raise ValueError(
            "cannot cache a failed, zero, or non-finite LiftOn profile"
        )
    document = {
        "schema_version": 1,
        "kind": "lifton_lift_cache",
        "key": _canonical_sha256(request),
        "request": request,
        "profile": profile.__dict__,
        "artifact": _file_fingerprint(out_gff),
    }
    temporary = done_flag.with_name(done_flag.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, done_flag)
    directory_fd = os.open(done_flag.parent, os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def _prepare_evaluation_target(out_gff: Path, evaluation_dir: Path) -> Path:
    """Expose the candidate GFF under an isolated evaluation directory.

    LiftOn derives ``lifton_output`` and gffutils sidecars from the directory
    containing ``-o``.  Passing the downloaded comparison annotation to ``-o``
    therefore evaluated the wrong file and wrote into the input tree; passing
    the original candidate path would overwrite the lift run's manifest and
    intermediates.  A same-filesystem hard link is cheap even for mammalian
    GFFs and gives evaluation its own artifact namespace.  Copy only when hard
    links are unavailable.
    """
    evaluation_dir.mkdir(parents=True, exist_ok=True)
    target = evaluation_dir / "candidate.gff3"
    eval_outdir = evaluation_dir / "lifton_output"
    for stale in (
        target,
        Path(str(target) + "_db"),
        Path(str(target) + ".eval_db"),
        eval_outdir / "eval.txt",
        eval_outdir / "eval.partial.txt",
        eval_outdir / "run_manifest.json",
    ):
        try:
            stale.unlink()
        except FileNotFoundError:
            pass
    try:
        os.link(out_gff, target)
    except OSError:
        shutil.copy2(out_gff, target)
    return target


def run_dataset(
    ds: Dataset, *,
    data_dir: Path, results_dir: Path,
    lifton_flags: list[str],
    evaluation_flags: list[str],
    do_download: bool = True,
    do_lift: bool = True,
    do_evaluate: bool = True,
    force: bool = False,
    log=print,
) -> dict:
    log(f"\n========================================================")
    log(f" {ds.id}: {ds.species}")
    log(f"========================================================")
    ds_data = data_dir / ds.id
    ds_results = results_dir / ds.id
    log_dir = ds_results / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: download ────────────────────────────────────────────────
    paths = {
        "reference_fa":  ds_data / _filename_for(ds.reference_fa),
        "target_fa":     ds_data / _filename_for(ds.target_fa),
        "reference_gff": ds_data / _filename_for(ds.reference_gff),
    }
    if ds.target_gff:
        paths["target_gff"] = ds_data / _filename_for(ds.target_gff)

    if do_download:
        log(f"\n[bench] {ds.id} — downloading inputs into {ds_data}")
        for url_key in ("reference_fa", "target_fa", "reference_gff"):
            url = getattr(ds, url_key)
            download(url, paths[url_key], log=log)
        if ds.target_gff:
            try:
                download(ds.target_gff, paths["target_gff"], log=log)
            except Exception as exc:
                log(f"  ! published comparison download failed (candidate "
                    f"evaluation will continue): {exc}")
                paths.pop("target_gff", None)

    if not do_lift:
        return {
            "dataset": ds.id, "species": ds.species,
            "downloaded_only": True,
            "input_paths": {k: str(v) for k, v in paths.items()},
        }

    # ── Step 2: run LiftOn ──────────────────────────────────────────────
    out_gff = ds_results / "lifton.gff3"
    done_flag = ds_results / ".lifton.done"
    argv = ["lifton"]
    argv.extend(lifton_flags)
    argv.extend(["-g", str(paths["reference_gff"])])
    argv.extend(["-o", str(out_gff)])
    # Positional: <target.fa> <ref.fa>
    argv.extend([str(paths["target_fa"]), str(paths["reference_fa"])])
    try:
        cache_request = _lift_cache_request(argv, paths)
    except (OSError, RuntimeError) as exc:
        return {
            "dataset": ds.id,
            "species": ds.species,
            "error": f"lift cache fingerprinting failed: {exc}",
        }
    lift_profile = (
        _cached_lift_profile(
            done_flag,
            request=cache_request,
            out_gff=out_gff,
        )
        if done_flag.exists() and not force else None
    )
    if lift_profile is not None:
        log(f"[bench] {ds.id} — lift result already on disk "
            f"with matching provenance ({done_flag})")
    else:
        try:
            done_flag.unlink()
        except FileNotFoundError:
            pass
        prof = run_profiled(argv, label="lift", log_dir=log_dir, log=log)
        if prof.exit_code != 0:
            log(f"  ! lifton failed with exit code {prof.exit_code}; "
                f"see {prof.stderr_path}")
            return {
                "dataset": ds.id, "species": ds.species,
                "lift_profile": prof.__dict__,
                "error": "lifton non-zero exit",
            }
        try:
            _write_lift_cache(
                done_flag,
                request=cache_request,
                profile=prof,
                out_gff=out_gff,
            )
        except (OSError, ValueError) as exc:
            return {
                "dataset": ds.id,
                "species": ds.species,
                "lift_profile": prof.__dict__,
                "error": f"lift cache publication failed: {exc}",
            }
        lift_profile = prof

    # ── Step 3: gather post-lift stats from LiftOn's own output dir ─────
    lifton_outdir = ds_results / "lifton_output"
    score_file = lifton_outdir / "score.txt"
    stats_subdir = lifton_outdir / "stats"
    score_summary = parse_score_txt(score_file)
    feature_summary = count_feature_stats(stats_subdir)

    # ── Step 4: optional --evaluation pass against truth set ────────────
    eval_profile = None
    evaluation_summary = None
    evaluation_target = None
    evaluation_manifest = None
    if do_evaluate:
        # Evaluate the annotation produced above, not the downloaded comparison
        # GFF.  Keep the evaluation copy in its own directory so its gffutils
        # database, run manifest, eval.txt, and intermediate files cannot
        # overwrite either input data or the profiled lift's artifacts.
        evaluation_dir = ds_results / "evaluation"
        evaluation_target = _prepare_evaluation_target(out_gff, evaluation_dir)
        eval_argv = ["lifton"]
        eval_argv.extend(evaluation_flags)
        eval_argv.extend(["-g", str(paths["reference_gff"])])
        eval_argv.extend(["-o", str(evaluation_target)])
        eval_argv.extend([str(paths["target_fa"]),
                          str(paths["reference_fa"])])
        eval_profile = run_profiled(
            eval_argv, label="evaluation", log_dir=log_dir, log=log,
        )
        eval_outdir = evaluation_dir / "lifton_output"
        eval_score = eval_outdir / "eval.txt"
        evaluation_manifest = eval_outdir / "run_manifest.json"
        evaluation_summary = parse_eval_txt(eval_score)

    identity_summary = evaluation_summary or score_summary
    mapped_features = feature_summary["mapped_features"]
    lost_features = feature_summary["lost_features"]
    reference_features = mapped_features + lost_features
    biological_summary = {
        "schema_version": 1,
        "reference_features": reference_features,
        "mapped_features": mapped_features,
        "lost_features": lost_features,
        "extra_copy_features": feature_summary["extra_copy_features"],
        "feature_completeness": (
            mapped_features / reference_features
            if reference_features else None
        ),
        "emitted_transcript_records": score_summary.records,
        "mapped_transcripts_reported": feature_summary["mapped_transcripts"],
        "evaluated_transcript_records": (
            evaluation_summary.records if evaluation_summary else None
        ),
        "evaluated_coding_records": (
            evaluation_summary.coding_records if evaluation_summary else None
        ),
        "mean_protein_identity": identity_summary.avg_identity,
    }
    # Keep the original five-key aggregate for downstream consumers, but make
    # mapped/lost explicitly feature-level and source identity from the correctly
    # parsed candidate evaluation.  New validation should use biological_summary.
    eval_summary = {
        "mapped": mapped_features,
        "lost": lost_features,
        "extra_copies": feature_summary["extra_copy_features"],
        "avg_identity": identity_summary.avg_identity,
        "score_file": identity_summary.score_file,
    }

    result = {
        "dataset": ds.id,
        "species": ds.species,
        "approx_size_gb": ds.approx_size_gb,
        "lift_profile": lift_profile.__dict__,
        "eval_profile": eval_profile.__dict__ if eval_profile else None,
        "eval_summary": eval_summary,
        "score_summary": score_summary.__dict__,
        "evaluation_summary": (
            evaluation_summary.__dict__ if evaluation_summary else None
        ),
        "biological_summary": biological_summary,
        "evaluation_target": (
            str(evaluation_target) if evaluation_target else None
        ),
        "evaluation_manifest": (
            str(evaluation_manifest) if evaluation_manifest else None
        ),
        "out_gff": str(out_gff),
        "input_paths": {k: str(v) for k, v in paths.items()},
    }
    biology_errors = validate_biological_result(
        result,
        require_evaluation=do_evaluate,
    )
    if eval_profile is not None and eval_profile.exit_code != 0:
        result["error"] = "evaluation non-zero exit"
    elif evaluation_summary is not None and evaluation_summary.records <= 0:
        result["error"] = "evaluation produced no valid biological records"
    elif (evaluation_summary is not None
          and evaluation_summary.malformed_records):
        result["error"] = "evaluation produced malformed biological records"
    elif biology_errors:
        result["error"] = "invalid biological summary: " + "; ".join(biology_errors)
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="run_benchmarks.py",
        description="LiftOn Phase 16 biological validation harness.",
    )
    p.add_argument("--registry", default=str(DEFAULT_REGISTRY),
                   help="Dataset registry JSON "
                        f"(default: {DEFAULT_REGISTRY.relative_to(REPO_ROOT)})")
    p.add_argument("--datasets", nargs="*", default=None,
                   help="Dataset ids to run (default: every entry in "
                        "the registry). Use --all as a synonym for "
                        "explicit selection.")
    p.add_argument("--all", action="store_true",
                   help="Equivalent to listing every dataset id.")
    p.add_argument("--data-dir", default=str(DEFAULT_DATA_DIR),
                   help=f"Where to store inputs "
                        f"(default: {DEFAULT_DATA_DIR.relative_to(REPO_ROOT)})")
    p.add_argument("--results-dir", default=str(DEFAULT_RESULTS_DIR),
                   help=f"Where to write outputs "
                        f"(default: {DEFAULT_RESULTS_DIR.relative_to(REPO_ROOT)})")
    p.add_argument("--output", default=None,
                   help="Roll-up JSON path. Default: "
                        "<results-dir>/summary_<UTC>.json")
    p.add_argument("--download-only", action="store_true",
                   help="Fetch inputs only — skip lift + evaluation.")
    p.add_argument("--no-evaluation", action="store_true",
                   help="Skip the --evaluation pass (lift only).")
    p.add_argument("--force", action="store_true",
                   help="Re-run lift even when .lifton.done exists.")
    args = p.parse_args(argv)

    registry = load_registry(Path(args.registry))
    selected_ids = (
        [d.id for d in registry.datasets]
        if (args.all or not args.datasets)
        else args.datasets
    )
    selected = [d for d in registry.datasets if d.id in selected_ids]
    missing = set(selected_ids) - {d.id for d in selected}
    if missing:
        p.error(f"unknown dataset id(s): {sorted(missing)}; "
                f"available: {[d.id for d in registry.datasets]}")

    data_dir = Path(args.data_dir).resolve()
    results_dir = Path(args.results_dir).resolve()
    data_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    runtime = _ensure_runtime()
    print("[bench] runtime:")
    for k, v in runtime.items():
        print(f"  {k:>10s} = {v}")
    if runtime["lifton"] is None:
        print("[bench] WARN: `lifton` is not on PATH; the lift step will "
              "fail unless you `pip install -e .` first.", file=sys.stderr)

    rows: list[dict] = []
    for ds in selected:
        try:
            row = run_dataset(
                ds,
                data_dir=data_dir, results_dir=results_dir,
                lifton_flags=registry.lifton_flags,
                evaluation_flags=registry.evaluation_flags,
                do_download=True,
                do_lift=not args.download_only,
                do_evaluate=not args.no_evaluation,
                force=args.force,
            )
        except Exception as exc:
            row = {
                "dataset": ds.id, "species": ds.species,
                "error": str(exc),
            }
            print(f"[bench] {ds.id} — FAILED: {exc}", file=sys.stderr)
        rows.append(row)

    # ── Summary table to stdout ──────────────────────────────────────────
    print("\n" + "=" * 78)
    print(f"{'dataset':<14} {'wall(s)':>10} {'peakRSS(MB)':>12} "
          f"{'mapped':>8} {'lost':>8} {'extra':>8} {'meanID':>8}")
    print("-" * 78)
    for row in rows:
        if "error" in row:
            print(f"{row['dataset']:<14} ERROR: {row['error']}")
            continue
        prof = row.get("lift_profile") or {}
        es = row.get("eval_summary") or {}
        wall = float(prof.get("wall_clock_seconds") or 0.0)
        rss = float(prof.get("peak_rss_mb") or 0.0)
        mapped = es.get("mapped") or 0
        lost = es.get("lost") or 0
        extra = es.get("extra_copies") or 0
        mid = es.get("avg_identity")
        mid_s = f"{mid:.4f}" if mid is not None else "n/a"
        print(f"{row['dataset']:<14} {wall:>10.1f} {rss:>12.1f} "
              f"{mapped:>8} {lost:>8} {extra:>8} {mid_s:>8}")
    print("=" * 78)

    out_json = (
        Path(args.output)
        if args.output
        else results_dir / f"summary_{time.strftime('%Y%m%dT%H%M%SZ')}.json"
    )
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(
        {"runtime": runtime, "rows": rows}, indent=2, default=str,
    ))
    print(f"\n[bench] wrote {out_json}")
    bad = [r for r in rows if "error" in r]
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
