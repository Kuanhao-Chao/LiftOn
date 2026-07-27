"""Run provenance and metrics helpers for LiftOn.

The module deliberately has no third-party requirements.  Pipeline code can
create a :class:`RunManifest`, update it while work progresses, and atomically
write the final JSON document without making provenance collection a new
failure mode.

Typical integration::

    manifest = RunManifest(
        argv=sys.argv,
        options=vars(args),
        inputs={"reference": args.reference, "target": args.target},
        backend={"annotation": "gffbase"},
        cache={"minimap2_index": "miss"},
    )
    with manifest.phase("liftoff"):
        run_liftoff()
    manifest.record_count("genes_written", gene_count)
    manifest.record_validation("gff3", passed=True, errors=0)
    manifest.finish("success")
    manifest.write("lifton.run.json")
"""

from __future__ import annotations

import copy
import concurrent.futures
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import stat
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator, Mapping, Sequence

from . import __version__


SCHEMA_VERSION = 2
REDACTED = "<redacted>"
DEFAULT_DEPENDENCIES = (
    "lifton",
    "numpy",
    "biopython",
    "parasail",
    "intervaltree",
    "interlap",
    "networkx",
    "pyfaidx",
    "pysam",
    "gffutils",
    "ujson",
    "duckdb",
    "pyarrow",
    "mappy",
)
DEFAULT_TOOL_COMMANDS = {
    "minimap2": ("minimap2", "--version"),
    "miniprot": ("miniprot", "--version"),
}

_SENSITIVE_KEY = re.compile(
    r"(?:^|[_-])(?:password|passwd|passphrase|secret|token|api[_-]?key|"
    r"credential|authorization|cookie)(?:$|[_-])",
    re.IGNORECASE,
)
_SENSITIVE_ASSIGNMENT = re.compile(
    r"\b(password|passwd|passphrase|secret|token|api[-_]?key|credential|"
    r"authorization|cookie)(\s*[:=]\s*)([^\s,;]+)",
    re.IGNORECASE,
)
_URL_USERINFO = re.compile(r"([a-z][a-z0-9+.-]*://)[^/@\s]+@", re.IGNORECASE)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _is_sensitive_key(key: object) -> bool:
    return bool(_SENSITIVE_KEY.search(str(key)))


def sanitize_text(value: object) -> str:
    """Return text with common credential assignments and URL userinfo removed."""
    text = str(value)
    text = _URL_USERINFO.sub(r"\1<redacted>@", text)
    return _SENSITIVE_ASSIGNMENT.sub(lambda match: match.group(1) + match.group(2) + REDACTED, text)


def sanitize_data(value: Any, *, _seen: set[int] | None = None) -> Any:
    """Convert arbitrary option data to JSON-safe values and redact secrets."""
    if value is None or isinstance(value, (bool, int, float)):
        return value
    if isinstance(value, str):
        return sanitize_text(value)
    if isinstance(value, os.PathLike):
        return sanitize_text(os.fspath(value))

    seen = set() if _seen is None else _seen
    value_id = id(value)
    if value_id in seen:
        return "<recursive>"

    if isinstance(value, Mapping):
        seen.add(value_id)
        sanitized = {}
        for key, item in value.items():
            string_key = str(key)
            sanitized[string_key] = REDACTED if _is_sensitive_key(key) else sanitize_data(item, _seen=seen)
        seen.remove(value_id)
        return sanitized

    if isinstance(value, (list, tuple, set, frozenset)):
        seen.add(value_id)
        items = value
        if isinstance(value, (set, frozenset)):
            items = sorted(value, key=repr)
        sanitized = [sanitize_data(item, _seen=seen) for item in items]
        seen.remove(value_id)
        return sanitized

    if hasattr(value, "__dict__"):
        return sanitize_data(vars(value), _seen=seen)
    return sanitize_text(value)


def sanitize_argv(argv: Sequence[object] | None) -> list[str]:
    """Sanitize a command line, redacting values of credential-like flags."""
    if argv is None:
        return []

    result: list[str] = []
    redact_next = False
    for raw_token in argv:
        token = str(raw_token)
        if redact_next:
            result.append(REDACTED)
            redact_next = False
            continue

        if token.startswith("--") and "=" in token:
            flag, value = token.split("=", 1)
            result.append(f"{flag}={REDACTED}" if _is_sensitive_key(flag) else sanitize_text(f"{flag}={value}"))
            continue

        result.append(sanitize_text(token))
        if token.startswith("--") and _is_sensitive_key(token):
            redact_next = True
    return result


def sanitize_options(options: Mapping[str, Any] | object | None) -> dict[str, Any]:
    """Return an argparse namespace or mapping as sanitized JSON data."""
    if options is None:
        return {}
    raw = options if isinstance(options, Mapping) else vars(options)
    sanitized = sanitize_data(raw)
    return sanitized if isinstance(sanitized, dict) else {"value": sanitized}


def collect_environment(environ: Mapping[str, str] | None = None) -> dict[str, str]:
    """Collect only explicitly LiftOn-scoped environment variables."""
    source = os.environ if environ is None else environ
    result = {}
    for key in sorted(source):
        if key.startswith("LIFTON_"):
            result[key] = REDACTED if _is_sensitive_key(key) else sanitize_text(source[key])
    return result


def _stat_signature(value: os.stat_result) -> tuple[int, int, int, int, int]:
    """Return the file identity and mutation-sensitive metadata."""
    return (
        int(value.st_dev),
        int(value.st_ino),
        int(value.st_size),
        int(value.st_mtime_ns),
        int(value.st_ctime_ns),
    )


def _stat_evidence(value: os.stat_result) -> dict[str, int]:
    return {
        "device": int(value.st_dev),
        "inode": int(value.st_ino),
        "size_bytes": int(value.st_size),
        "mtime_ns": int(value.st_mtime_ns),
        "ctime_ns": int(value.st_ctime_ns),
    }


def fingerprint_input(path: str | os.PathLike[str] | None, block_size: int = 1024 * 1024) -> dict[str, Any]:
    """Return a stable SHA-256 fingerprint, recording errors instead of raising."""
    display_path = None if path is None else sanitize_text(os.fspath(path))
    result: dict[str, Any] = {
        "path": display_path,
        "exists": False,
        "is_file": False,
        "size_bytes": None,
        "mtime_ns": None,
        "sha256": None,
        "bytes_read": 0,
        "changed_during_hash": False,
        "fingerprint_status": "pending",
    }
    if path is None:
        result["error"] = "no path supplied"
        result["fingerprint_status"] = "unavailable"
        return result

    initial_path_stat = None
    try:
        file_path = Path(path)
        initial_path_stat = file_path.stat()
        result.update(
            exists=True,
            is_file=stat.S_ISREG(initial_path_stat.st_mode),
            size_bytes=initial_path_stat.st_size,
            mtime_ns=initial_path_stat.st_mtime_ns,
            identity={
                "device": int(initial_path_stat.st_dev),
                "inode": int(initial_path_stat.st_ino),
            },
        )
        if not result["is_file"]:
            result["error"] = "path is not a regular file"
            result["fingerprint_status"] = "unavailable"
            return result

        digest = hashlib.sha256()
        with file_path.open("rb") as handle:
            opened_stat = os.fstat(handle.fileno())
            while True:
                block = handle.read(block_size)
                if not block:
                    break
                digest.update(block)
                result["bytes_read"] += len(block)
            final_handle_stat = os.fstat(handle.fileno())
        try:
            final_path_stat = file_path.stat()
        except OSError as exc:
            result.update(
                changed_during_hash=True,
                fingerprint_status="changed",
                error=(
                    "input path changed while computing SHA-256: "
                    f"{type(exc).__name__}: {sanitize_text(exc)}"
                ),
            )
            return result
        signatures = {
            _stat_signature(initial_path_stat),
            _stat_signature(opened_stat),
            _stat_signature(final_handle_stat),
            _stat_signature(final_path_stat),
        }
        if len(signatures) != 1:
            result.update(
                changed_during_hash=True,
                fingerprint_status="changed",
                stat_after=_stat_evidence(final_path_stat),
                error="input changed while computing SHA-256",
            )
            return result
        result["sha256"] = digest.hexdigest()
        result["fingerprint_status"] = "complete"
    except OSError as exc:
        if initial_path_stat is not None:
            try:
                error_path_stat = file_path.stat()
            except OSError:
                error_path_stat = None
            changed = (
                error_path_stat is None
                or _stat_signature(error_path_stat)
                != _stat_signature(initial_path_stat)
            )
            if changed:
                result.update(
                    changed_during_hash=True,
                    fingerprint_status="changed",
                    error=(
                        "input path changed or became unreadable while "
                        f"computing SHA-256: {type(exc).__name__}: "
                        f"{sanitize_text(exc)}"
                    ),
                )
                if error_path_stat is not None:
                    result["stat_after"] = _stat_evidence(error_path_stat)
            else:
                result["error"] = (
                    f"{type(exc).__name__}: {sanitize_text(exc)}"
                )
                result["fingerprint_status"] = "unavailable"
        else:
            result["error"] = f"{type(exc).__name__}: {sanitize_text(exc)}"
            result["fingerprint_status"] = "unavailable"
    except (ValueError, TypeError) as exc:
        result["error"] = f"{type(exc).__name__}: {sanitize_text(exc)}"
        result["fingerprint_status"] = "unavailable"
    return result


def fingerprint_inputs(inputs: Mapping[str, str | os.PathLike[str] | None] | None) -> dict[str, dict[str, Any]]:
    """Fingerprint named inputs in deterministic key order."""
    if not inputs:
        return {}
    return {str(name): fingerprint_input(inputs[name]) for name in sorted(inputs, key=str)}


def _pending_fingerprints(
    inputs: Mapping[str, str | os.PathLike[str] | None],
) -> dict[str, dict[str, Any]]:
    pending = {}
    for name in sorted(inputs, key=str):
        path = inputs[name]
        pending[str(name)] = {
            "path": None if path is None else sanitize_text(os.fspath(path)),
            "exists": None,
            "is_file": None,
            "size_bytes": None,
            "mtime_ns": None,
            "sha256": None,
            "bytes_read": 0,
            "changed_during_hash": False,
            "fingerprint_status": "pending",
        }
    return pending


def _fingerprint_inputs_task(
    inputs: Mapping[str, str | os.PathLike[str] | None],
) -> dict[str, Any]:
    started_at = _utc_now()
    started = time.perf_counter()
    fingerprints = fingerprint_inputs(inputs)
    duration = time.perf_counter() - started
    changed_names = sorted(
        name for name, value in fingerprints.items()
        if value.get("changed_during_hash")
    )
    metrics = {
        "input_count": len(fingerprints),
        "regular_files": sum(
            1 for value in fingerprints.values() if value.get("is_file")
        ),
        "hashed_files": sum(
            1 for value in fingerprints.values()
            if value.get("fingerprint_status") == "complete"
        ),
        "bytes_hashed": sum(
            int(value.get("bytes_read") or 0)
            for value in fingerprints.values()
        ),
        "unavailable_inputs": sum(
            1 for value in fingerprints.values()
            if value.get("fingerprint_status") == "unavailable"
        ),
        "changed_inputs": len(changed_names),
        "changed_input_names": changed_names,
    }
    return {
        "inputs": fingerprints,
        "started_at": started_at,
        "finished_at": _utc_now(),
        "duration_seconds": duration,
        "metrics": metrics,
    }


def collect_dependency_versions(names: Sequence[str] = DEFAULT_DEPENDENCIES) -> dict[str, str | None]:
    """Return installed distribution versions; missing packages map to ``None``."""
    versions: dict[str, str | None] = {}
    for name in names:
        try:
            versions[str(name)] = importlib.metadata.version(name)
        except (importlib.metadata.PackageNotFoundError, ValueError, TypeError):
            versions[str(name)] = None
    return versions


def probe_tool_version(command: Sequence[str] | str, timeout: float = 5.0) -> dict[str, Any]:
    """Run a version command without a shell and return a non-throwing result."""
    argv = [command] if isinstance(command, str) else [str(part) for part in command]
    result: dict[str, Any] = {
        "command": sanitize_argv(argv),
        "available": False,
        "version": None,
    }
    if not argv:
        result["error"] = "empty command"
        return result

    try:
        completed = subprocess.run(
            argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
            timeout=timeout,
        )
        output = completed.stdout or ""
        first_line = next((line.strip() for line in output.splitlines() if line.strip()), None)
        result.update(
            available=True,
            version=sanitize_text(first_line[:500]) if first_line else None,
            returncode=completed.returncode,
        )
        if completed.returncode != 0:
            result["error"] = f"exit status {completed.returncode}"
    except FileNotFoundError:
        result["error"] = "not found"
    except subprocess.TimeoutExpired:
        result["error"] = f"timed out after {timeout:g}s"
    except (OSError, ValueError) as exc:
        result["error"] = f"{type(exc).__name__}: {sanitize_text(exc)}"
    return result


def collect_tool_versions(
    commands: Mapping[str, Sequence[str] | str] = DEFAULT_TOOL_COMMANDS,
    timeout: float = 5.0,
) -> dict[str, dict[str, Any]]:
    """Probe named external tools independently so one failure cannot abort collection."""
    return {str(name): probe_tool_version(commands[name], timeout) for name in sorted(commands, key=str)}


def _git_value(arguments: Sequence[str], cwd: str | os.PathLike[str] | None, timeout: float) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *arguments],
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
            text=True,
            timeout=timeout,
        )
    except (FileNotFoundError, OSError, ValueError, subprocess.TimeoutExpired):
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip() or None


def collect_git_metadata(cwd: str | os.PathLike[str] | None = None, timeout: float = 5.0) -> dict[str, Any]:
    """Collect Git executable/repository data without recording remotes."""
    tool = probe_tool_version(("git", "--version"), timeout)
    result: dict[str, Any] = {
        "available": tool["available"],
        "version": tool["version"],
        "repository": False,
        "commit": None,
        "branch": None,
        "dirty": None,
    }
    if not tool["available"]:
        result["error"] = tool.get("error")
        return result

    commit = _git_value(("rev-parse", "HEAD"), cwd, timeout)
    if commit is None:
        result["error"] = "not a Git work tree"
        return result

    branch = _git_value(("symbolic-ref", "--quiet", "--short", "HEAD"), cwd, timeout)
    status = _git_value(("status", "--porcelain", "--untracked-files=no"), cwd, timeout)
    result.update(repository=True, commit=commit, branch=branch, dirty=bool(status))
    return result


def process_rss() -> dict[str, Any]:
    """Return current and peak resident memory where the platform exposes it."""
    current_bytes = None
    peak_bytes = None
    sources: list[str] = []

    try:
        with open("/proc/self/status", encoding="ascii") as handle:
            for line in handle:
                if line.startswith("VmRSS:"):
                    current_bytes = int(line.split()[1]) * 1024
                elif line.startswith("VmHWM:"):
                    peak_bytes = int(line.split()[1]) * 1024
        if current_bytes is not None or peak_bytes is not None:
            sources.append("/proc/self/status")
    except (OSError, ValueError, IndexError):
        pass

    try:
        import resource

        maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        resource_peak = int(maximum if sys.platform == "darwin" else maximum * 1024)
        peak_bytes = max(peak_bytes or 0, resource_peak)
        sources.append("resource.getrusage")
    except (ImportError, OSError, ValueError, AttributeError):
        pass

    return {
        "current_bytes": current_bytes,
        "peak_bytes": peak_bytes,
        "source": "+".join(sources) if sources else None,
    }


def atomic_write_json(path: str | os.PathLike[str], document: Any, mode: int = 0o600) -> None:
    """Atomically replace *path* with a formatted JSON document.

    The temporary file is created in the destination directory, flushed, and
    fsynced before :func:`os.replace`.  On serialization or I/O failure, the
    old destination remains untouched and the temporary file is removed.
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=target.parent)
    replaced = False
    try:
        os.chmod(temporary_name, mode)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            descriptor = -1
            json.dump(document, handle, indent=2, sort_keys=True, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, target)
        replaced = True

        try:
            directory_fd = os.open(target.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        except OSError:
            pass
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if not replaced:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


class RunManifest:
    """Mutable builder for one LiftOn run manifest."""

    def __init__(
        self,
        *,
        argv: Sequence[object] | None = None,
        options: Mapping[str, Any] | object | None = None,
        inputs: Mapping[str, str | os.PathLike[str] | None] | None = None,
        backend: Mapping[str, Any] | None = None,
        cache: Mapping[str, Any] | None = None,
        dependency_names: Sequence[str] = DEFAULT_DEPENDENCIES,
        tool_commands: Mapping[str, Sequence[str] | str] = DEFAULT_TOOL_COMMANDS,
        git_cwd: str | os.PathLike[str] | None = None,
        collect_git: bool = True,
    ) -> None:
        self._run_started = time.perf_counter()
        self._phase_starts: dict[str, float] = {}
        self._finished = False
        self._fingerprints_complete = False
        self._fingerprint_executor: concurrent.futures.ThreadPoolExecutor | None = None
        self._fingerprint_future: concurrent.futures.Future[dict[str, Any]] | None = None

        normalized_inputs = {
            str(name): inputs[name] for name in sorted(inputs or {}, key=str)
        }
        fingerprint_queued_at = None
        if normalized_inputs:
            fingerprint_queued_at = _utc_now()
            self._fingerprint_executor = concurrent.futures.ThreadPoolExecutor(
                max_workers=1,
                thread_name_prefix="lifton-input-fingerprint",
            )
            self._fingerprint_future = self._fingerprint_executor.submit(
                _fingerprint_inputs_task, normalized_inputs,
            )
        else:
            self._fingerprints_complete = True

        git = collect_git_metadata(git_cwd) if collect_git else None
        phases = {}
        if normalized_inputs:
            phases["fingerprint_inputs"] = {
                "started_at": fingerprint_queued_at,
                "finished_at": None,
                "duration_seconds": None,
                "status": "running",
                "details": {
                    "algorithm": "sha256",
                    "execution": "background-thread",
                    "block_size_bytes": 1024 * 1024,
                },
                "metrics": {},
            }
        self._document: dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "run": {
                "started_at": _utc_now(),
                "finished_at": None,
                "duration_seconds": None,
                "status": "running",
                "argv": sanitize_argv(argv),
                "options": sanitize_options(options),
                "environment": collect_environment(),
                "backend": sanitize_data(backend or {}),
                "cache": sanitize_data(cache or {}),
            },
            "software": {
                "lifton": __version__,
                "python": {
                    "version": platform.python_version(),
                    "implementation": platform.python_implementation(),
                    "executable": sys.executable,
                },
                "platform": {
                    "system": platform.system(),
                    "release": platform.release(),
                    "machine": platform.machine(),
                },
                "git": git,
                "dependencies": collect_dependency_versions(dependency_names),
                "tools": collect_tool_versions(tool_commands),
            },
            "inputs": _pending_fingerprints(normalized_inputs),
            "phases": phases,
            "counts": {},
            "failures": [],
            "validation": {},
            "resources": {"start": process_rss(), "end": None, "last_write": None},
        }

    def _complete_input_fingerprints(self) -> None:
        """Join the background hash worker and publish only stable evidence."""
        if self._fingerprints_complete:
            return

        future = self._fingerprint_future
        executor = self._fingerprint_executor
        join_started = time.perf_counter()
        try:
            if future is None:
                raise RuntimeError("input fingerprint worker was not started")
            result = future.result()
            join_wait = time.perf_counter() - join_started
            self._document["inputs"] = result["inputs"]
            metrics = dict(result["metrics"])
            metrics["join_wait_seconds"] = join_wait
            phase = self._document["phases"]["fingerprint_inputs"]
            phase.update(
                started_at=result["started_at"],
                finished_at=result["finished_at"],
                duration_seconds=result["duration_seconds"],
                status=(
                    "failed" if metrics["changed_inputs"] else "success"
                ),
                metrics=metrics,
            )
            if metrics["changed_inputs"]:
                self.record_failure(
                    "fingerprint_inputs",
                    "input changed while computing SHA-256",
                    details={
                        "inputs": metrics["changed_input_names"],
                    },
                )
        except Exception as exc:
            error = (
                "input fingerprint worker failed: "
                f"{type(exc).__name__}: {sanitize_text(exc)}"
            )
            for fingerprint in self._document["inputs"].values():
                if fingerprint.get("fingerprint_status") == "pending":
                    fingerprint.update(
                        fingerprint_status="unavailable",
                        error=error,
                    )
            join_wait = time.perf_counter() - join_started
            phase = self._document["phases"]["fingerprint_inputs"]
            phase.update(
                finished_at=_utc_now(),
                duration_seconds=join_wait,
                status="failed",
                metrics={
                    "input_count": len(self._document["inputs"]),
                    "regular_files": 0,
                    "hashed_files": 0,
                    "bytes_hashed": 0,
                    "unavailable_inputs": len(self._document["inputs"]),
                    "changed_inputs": 0,
                    "changed_input_names": [],
                    "join_wait_seconds": join_wait,
                },
            )
            self.record_failure("fingerprint_inputs", exc)
        finally:
            if executor is not None:
                executor.shutdown(wait=True)
            self._fingerprint_future = None
            self._fingerprint_executor = None
            self._fingerprints_complete = True

    def set_backend_choice(self, name: str, value: Any) -> None:
        sanitized = sanitize_data({str(name): value})
        self._document["run"]["backend"].update(sanitized)

    def set_cache_choice(self, name: str, value: Any) -> None:
        sanitized = sanitize_data({str(name): value})
        self._document["run"]["cache"].update(sanitized)

    def start_phase(self, name: str, details: Mapping[str, Any] | None = None) -> None:
        """Start a uniquely named phase using a monotonic timer."""
        name = str(name)
        if name in self._phase_starts:
            raise ValueError(f"phase already active: {name}")
        if name in self._document["phases"]:
            raise ValueError(f"phase already recorded: {name}")
        self._phase_starts[name] = time.perf_counter()
        self._document["phases"][name] = {
            "started_at": _utc_now(),
            "finished_at": None,
            "duration_seconds": None,
            "status": "running",
            "details": sanitize_data(details or {}),
        }

    def end_phase(
        self,
        name: str,
        *,
        status: str = "success",
        metrics: Mapping[str, Any] | None = None,
    ) -> float:
        """End an active phase and return its elapsed seconds."""
        name = str(name)
        try:
            started = self._phase_starts.pop(name)
        except KeyError as exc:
            raise ValueError(f"phase is not active: {name}") from exc
        elapsed = time.perf_counter() - started
        phase = self._document["phases"][name]
        phase.update(
            finished_at=_utc_now(),
            duration_seconds=elapsed,
            status=sanitize_text(status),
            metrics=sanitize_data(metrics or {}),
        )
        return elapsed

    @contextmanager
    def phase(self, name: str, details: Mapping[str, Any] | None = None) -> Iterator[None]:
        """Context manager that records success or an exception-backed failure."""
        self.start_phase(name, details)
        try:
            yield
        except BaseException as exc:
            self.end_phase(name, status="failed")
            self.record_failure(name, exc)
            raise
        else:
            self.end_phase(name)

    def record_count(self, name: str, value: int) -> None:
        self._document["counts"][str(name)] = int(value)

    def increment_count(self, name: str, amount: int = 1) -> int:
        name = str(name)
        value = int(self._document["counts"].get(name, 0)) + int(amount)
        self._document["counts"][name] = value
        return value

    def record_failure(
        self,
        phase: str,
        error: BaseException | str,
        *,
        fatal: bool = False,
        details: Mapping[str, Any] | None = None,
    ) -> None:
        error_type = type(error).__name__ if isinstance(error, BaseException) else "Error"
        self._document["failures"].append(
            {
                "time": _utc_now(),
                "phase": str(phase),
                "type": error_type,
                "message": sanitize_text(error),
                "fatal": bool(fatal),
                "details": sanitize_data(details or {}),
            }
        )

    def record_validation(
        self,
        name: str,
        *,
        passed: bool,
        errors: int = 0,
        warnings: int = 0,
        details: Mapping[str, Any] | None = None,
    ) -> None:
        self._document["validation"][str(name)] = {
            "passed": bool(passed),
            "errors": int(errors),
            "warnings": int(warnings),
            "details": sanitize_data(details or {}),
        }

    def finish(self, status: str = "success") -> None:
        """Finalize the run, marking any still-active phases incomplete."""
        if self._finished:
            return
        for phase_name in list(self._phase_starts):
            self.end_phase(phase_name, status="incomplete")
        self._complete_input_fingerprints()
        self._document["run"].update(
            finished_at=_utc_now(),
            duration_seconds=time.perf_counter() - self._run_started,
            status=sanitize_text(status),
        )
        self._document["resources"]["end"] = process_rss()
        self._finished = True

    def to_dict(self) -> dict[str, Any]:
        """Return an independent snapshot safe for caller modification."""
        return copy.deepcopy(self._document)

    def write(self, path: str | os.PathLike[str]) -> None:
        """Atomically write the current snapshot without implicitly finalizing it."""
        self._complete_input_fingerprints()
        self._document["resources"]["last_write"] = process_rss()
        atomic_write_json(path, self._document)
