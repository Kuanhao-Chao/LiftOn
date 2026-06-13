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
  4. Optionally re-invokes LiftOn in `--evaluation` (-E) mode against
     the dataset's pre-computed truth set to extract mapped /
     lost / extra-copy / identity metrics.
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
import json
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


@dataclass
class Registry:
    datasets: list[Dataset]
    lifton_flags: list[str] = field(default_factory=list)
    evaluation_flags: list[str] = field(default_factory=list)


def load_registry(path: Path) -> Registry:
    raw = json.loads(path.read_text())
    # Strip free-form comment keys (anything starting with "_").
    ds_entries = [d for d in raw["datasets"]
                  if not (isinstance(d, dict) and
                          all(k.startswith("_") for k in d.keys()))]
    datasets = [Dataset(**d) for d in ds_entries]
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


def parse_score_txt(path: Path) -> EvalSummary:
    """Parse the LiftOn score.txt file into mapped/lost/identity stats.

    The file format is a tab-separated table emitted by
    `lifton_utils.write_lifton_status` — one row per transcript:

        <transcript_id>\\t<chr>\\t<start>\\t<end>\\t<status>\\t<dna_id>\\t<aa_id>\\t...

    We only need aggregate counts for the harness summary.
    """
    if not path.exists():
        return EvalSummary(score_file=str(path))
    mapped = 0
    lost = 0
    identities: list[float] = []
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            status = parts[4].strip().lower()
            try:
                aa_id = float(parts[6])
            except ValueError:
                aa_id = 0.0
            if status in ("mapped", "liftoff", "miniprot", "lifton"):
                mapped += 1
                identities.append(aa_id)
            elif status in ("lost", "unmapped"):
                lost += 1
    avg_id = sum(identities) / len(identities) if identities else None
    return EvalSummary(
        mapped=mapped,
        lost=lost,
        avg_identity=avg_id,
        score_file=str(path),
    )


def count_unmapped_extra(stats_dir: Path) -> tuple[int, int]:
    """LiftOn writes plaintext lists at:
        stats/unmapped_features.txt   — one feature id per line
        stats/extra_copy_features.txt — same shape

    Both are produced regardless of --evaluation; line-count gives
    the per-run aggregates we need."""
    def _count(name: str) -> int:
        p = stats_dir / name
        if not p.exists():
            return 0
        return sum(1 for _ in p.open())
    return _count("unmapped_features.txt"), _count("extra_copy_features.txt")


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
                log(f"  ! truth-set download failed (evaluation will be "
                    f"skipped): {exc}")
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
    if done_flag.exists() and not force:
        log(f"[bench] {ds.id} — lift result already on disk "
            f"(remove {done_flag} or pass --force to redo)")
    else:
        argv = ["lifton"]
        argv.extend(lifton_flags)
        argv.extend(["-g", str(paths["reference_gff"])])
        argv.extend(["-o", str(out_gff)])
        # Positional: <target.fa> <ref.fa>
        argv.extend([str(paths["target_fa"]), str(paths["reference_fa"])])
        prof = run_profiled(argv, label="lift", log_dir=log_dir, log=log)
        if prof.exit_code != 0:
            log(f"  ! lifton failed with exit code {prof.exit_code}; "
                f"see {prof.stderr_path}")
            return {
                "dataset": ds.id, "species": ds.species,
                "lift_profile": prof.__dict__,
                "error": "lifton non-zero exit",
            }
        done_flag.write_text(time.strftime("%Y-%m-%dT%H:%M:%SZ"))
        lift_profile = prof
    if "lift_profile" not in dir():
        # Re-stat the existing run without a profile snapshot.
        lift_profile = ProfileResult(
            wall_clock_seconds=0.0,
            peak_rss_mb=0.0,
            user_cpu_seconds=0.0,
            sys_cpu_seconds=0.0,
            exit_code=0,
            stdout_path="<cached>",
            stderr_path="<cached>",
            time_log_path=None,
        )

    # ── Step 3: gather post-lift stats from LiftOn's own output dir ─────
    lifton_outdir = ds_results / "lifton_output"
    score_file = lifton_outdir / "score.txt"
    stats_subdir = lifton_outdir / "stats"
    eval_summary = parse_score_txt(score_file)
    unmapped_n, extra_n = count_unmapped_extra(stats_subdir)
    eval_summary.lost = (eval_summary.lost or 0) + unmapped_n
    eval_summary.extra_copies = extra_n

    # ── Step 4: optional --evaluation pass against truth set ────────────
    eval_profile = None
    if do_evaluate and ds.target_gff and "target_gff" in paths:
        eval_argv = ["lifton"]
        eval_argv.extend(evaluation_flags)
        eval_argv.extend(["-g", str(paths["reference_gff"])])
        eval_argv.extend(["-o", str(paths["target_gff"])])
        eval_argv.extend([str(paths["target_fa"]),
                          str(paths["reference_fa"])])
        eval_profile = run_profiled(
            eval_argv, label="evaluation", log_dir=log_dir, log=log,
        )
        eval_score = lifton_outdir / "eval.txt"
        if eval_score.exists():
            eval_summary = parse_score_txt(eval_score)

    return {
        "dataset": ds.id,
        "species": ds.species,
        "approx_size_gb": ds.approx_size_gb,
        "lift_profile": lift_profile.__dict__,
        "eval_profile": eval_profile.__dict__ if eval_profile else None,
        "eval_summary": eval_summary.__dict__,
        "out_gff": str(out_gff),
        "input_paths": {k: str(v) for k, v in paths.items()},
    }


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
