"""Resource-aware scheduling and bounded external-tool diagnostics.

The native aligners are independent programs.  LiftOn can decide when to
overlap them, but it cannot turn a target index that is larger than available
memory into a successful one.  This module therefore keeps policy separate
from biology: it chooses a schedule from target size and explicit overrides,
and records factual process evidence without guessing that a signal was OOM.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile
import threading
import time
from typing import Any, Iterable, Mapping, Sequence


MAX_CONCURRENT_TARGET_BASES = 4_000_000_000
# miniprot <=0.13 returned a signed 32-bit sequence length from its FASTA
# reader; 2^31 is therefore the first unrepresentable positive length.
MINIPROT_LONG_SEQUENCE_BASES = 1 << 31
MINIPROT_LONG_SEQUENCE_MIN_VERSION = (0, 14)
STDERR_TAIL_BYTES = 64 << 10
EXTERNAL_ALIGNER_INSTALL_HELP = (
    "pip does not install the minimap2/miniprot executables. Install them with "
    "`conda install --override-channels -c conda-forge -c bioconda "
    "--strict-channel-priority minimap2 miniprot` and ensure they are on PATH. "
    "See https://khchao.com/LiftOn/content/installation.html."
)


@dataclass(frozen=True)
class AlignerSchedule:
    """Resolved execution policy for the external DNA/protein aligners."""

    mode: str
    reason: str
    active_tools: tuple[str, ...]
    forced: bool
    boundary_bases: int = MAX_CONCURRENT_TARGET_BASES

    @property
    def parallel(self) -> bool:
        return self.mode == "parallel"

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["active_tools"] = list(self.active_tools)
        value["parallel"] = self.parallel
        return value


@dataclass(frozen=True)
class CommandResult:
    """A process result with a bounded tail and stream-wide ERROR sentinel."""

    returncode: int
    stderr_tail: str
    stderr_error_seen: bool = False


class BoundedStderr:
    """Thread-safe byte tail and ERROR sentinel for relayed stderr."""

    def __init__(self, limit: int = STDERR_TAIL_BYTES):
        self.limit = max(1, int(limit))
        self._tail = bytearray()
        self._error_carry = b""
        self.error_seen = False
        self._lock = threading.Lock()

    def feed(self, chunk: bytes | str) -> None:
        payload = chunk.encode("utf-8", errors="replace") \
            if isinstance(chunk, str) else bytes(chunk)
        with self._lock:
            upper = self._error_carry + payload.upper()
            if b"ERROR" in upper:
                self.error_seen = True
            self._error_carry = upper[-4:]
            self._tail.extend(payload)
            if len(self._tail) > self.limit:
                del self._tail[:-self.limit]

    @property
    def text(self) -> str:
        with self._lock:
            payload = bytes(self._tail)
        return payload.decode("utf-8", errors="replace")


def target_fasta_statistics(fasta: Any) -> dict[str, int]:
    """Summarize an opened :class:`pyfaidx.Fasta` without reading sequences."""

    index = getattr(getattr(fasta, "faidx", None), "index", None)
    if index is None:
        lengths: Iterable[int] = (len(fasta[name]) for name in fasta.keys())
    else:
        lengths = (int(record.rlen) for record in index.values())

    sequence_count = 0
    total_bases = 0
    maximum_sequence_bases = 0
    for length in lengths:
        length = int(length)
        sequence_count += 1
        total_bases += length
        maximum_sequence_bases = max(maximum_sequence_bases, length)
    return {
        "sequence_count": sequence_count,
        "total_bases": total_bases,
        "maximum_sequence_bases": maximum_sequence_bases,
    }


def resolve_aligner_schedule(
    target_bases: int,
    *,
    needs_liftoff: bool,
    needs_miniprot: bool,
    force_serial: bool = False,
    force_parallel: bool = False,
) -> AlignerSchedule:
    """Resolve one deterministic schedule from work present and target size."""

    if force_serial and force_parallel:
        raise ValueError("serial and parallel aligner overrides are exclusive")
    active = tuple(
        tool for tool, needed in (
            ("liftoff", needs_liftoff),
            ("miniprot", needs_miniprot),
        ) if needed
    )
    if not active:
        return AlignerSchedule("cached", "all_aligner_outputs_precomputed", active, False)
    if len(active) == 1:
        return AlignerSchedule("single", "only_one_aligner_required", active, False)
    if force_serial:
        return AlignerSchedule("serial", "user_forced_serial", active, True)
    if force_parallel:
        return AlignerSchedule("parallel", "user_forced_parallel", active, True)
    if int(target_bases) > MAX_CONCURRENT_TARGET_BASES:
        return AlignerSchedule(
            "serial", "target_exceeds_concurrent_boundary", active, False,
        )
    return AlignerSchedule(
        "parallel", "target_within_concurrent_boundary", active, False,
    )


def describe_returncode(returncode: int | None) -> str:
    """Render POSIX signal return codes without asserting an OOM cause."""

    if returncode is None:
        return "process did not start"
    value = int(returncode)
    if value >= 0:
        return f"exit status {value}"
    signum = -value
    try:
        name = signal.Signals(signum).name
    except ValueError:
        name = f"SIGNAL_{signum}"
    return f"terminated by signal {signum} ({name})"


def returncode_evidence(returncode: int | None) -> dict[str, Any]:
    """Return machine-readable status and signal facts."""

    result: dict[str, Any] = {
        "returncode": None if returncode is None else int(returncode),
        "termination": describe_returncode(returncode),
        "signal": None,
    }
    if returncode is not None and int(returncode) < 0:
        signum = -int(returncode)
        try:
            name = signal.Signals(signum).name
        except ValueError:
            name = None
        try:
            description = signal.strsignal(signum)
        except (AttributeError, ValueError):
            description = None
        result["signal"] = {
            "number": signum,
            "name": name,
            "description": description,
        }
    return result


def _relay_stderr(chunk: bytes) -> None:
    stream = sys.stderr
    try:
        binary = getattr(stream, "buffer", None)
        if binary is not None:
            binary.write(chunk)
            binary.flush()
        else:
            stream.write(chunk.decode("utf-8", errors="replace"))
            stream.flush()
    except (AttributeError, OSError, ValueError):
        pass


def run_with_bounded_stderr(
    command: Sequence[str],
    *,
    stdout: Any = None,
    tail_limit: int = STDERR_TAIL_BYTES,
    relay_stderr: bool = True,
) -> CommandResult:
    """Run a command, relay stderr live, and retain only its bounded tail."""

    proc = subprocess.Popen(
        [str(part) for part in command], stdout=stdout, stderr=subprocess.PIPE,
    )
    state = BoundedStderr(tail_limit)

    def _drain() -> None:
        if proc.stderr is None:
            return
        try:
            while True:
                chunk = proc.stderr.read(1 << 16)
                if not chunk:
                    return
                state.feed(chunk)
                if relay_stderr:
                    _relay_stderr(chunk)
        finally:
            try:
                proc.stderr.close()
            except Exception:
                pass

    thread = threading.Thread(
        target=_drain, name="lifton-tool-stderr", daemon=True,
    )
    thread.start()
    try:
        returncode = proc.wait()
    except BaseException:
        try:
            proc.terminate()
            proc.wait(timeout=5)
        except Exception:
            try:
                proc.kill()
                proc.wait(timeout=5)
            except Exception:
                pass
        thread.join(timeout=5)
        raise
    thread.join()
    return CommandResult(
        returncode=returncode,
        stderr_tail=state.text,
        stderr_error_seen=state.error_seen,
    )


def bounded_text_tail(text: bytes | str | None, limit: int = STDERR_TAIL_BYTES) -> str:
    state = BoundedStderr(limit)
    if text:
        state.feed(text)
    return state.text


def detect_miniprot_stage(stderr_text: str) -> str:
    """Identify the last miniprot indexing milestone visible in stderr."""

    text = stderr_text.lower()
    if "kmer-block pairs" in text:
        return "kmer_block_pairs_built"
    if "collected syncmers" in text:
        return "syncmers_collected"
    if re.search(r"\b\d+\s+blocks\b", text):
        return "target_blocks_built"
    if re.search(r"read\s+\d+\s+bases", text):
        return "target_loaded"
    return "process_started"


def detect_minimap2_stage(stderr_text: str, operation: str) -> str:
    """Identify the last minimap2 milestone, retaining an operation fallback."""

    text = stderr_text.lower()
    if re.search(r"\bmapped\s+\d+\s+sequences?\b", text):
        return "query_mapping_complete"
    if "loaded/built the index" in text:
        return "target_index_ready"
    if "sorted minimizers" in text:
        return "target_minimizers_sorted"
    if "collected minimizers" in text:
        return "target_minimizers_collected"
    return str(operation)


def parse_miniprot_version(banner: str | None) -> tuple[int, int] | None:
    if not banner:
        return None
    match = re.search(r"(?<!\d)(\d+)\.(\d+)(?!\d)", str(banner))
    if match is None:
        return None
    return int(match.group(1)), int(match.group(2))


def miniprot_long_sequence_compatibility(
    maximum_sequence_bases: int, banner: str | None,
) -> str:
    """Classify the known pre-0.14 long-sequence compatibility boundary."""

    if int(maximum_sequence_bases) < MINIPROT_LONG_SEQUENCE_BASES:
        return "not_applicable"
    version = parse_miniprot_version(banner)
    if version is None:
        return "unknown"
    if version < MINIPROT_LONG_SEQUENCE_MIN_VERSION:
        return "incompatible"
    return "compatible"


def execution_record(
    tool: str,
    *,
    command: Sequence[str] | None,
    stage: str,
    status: str,
    returncode: int | None,
    stderr_tail: str = "",
    details: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    record = {
        "tool": str(tool),
        "event_time_ns": time.time_ns(),
        "process_id": os.getpid(),
        "command": [str(part) for part in command] if command else None,
        "stage": str(stage),
        "status": str(status),
        "stderr_tail": bounded_text_tail(stderr_tail),
        **returncode_evidence(returncode),
    }
    if details:
        record["details"] = dict(details)
    return record


def record_execution(args: Any, record: Mapping[str, Any]) -> None:
    """Record directly, or use one per-process JSON event across a worker."""

    manifest = getattr(args, "_run_manifest", None)
    if manifest is not None:
        manifest.record_aligner_execution(record)
        return
    directory = getattr(args, "aligner_diagnostics_dir", None)
    if not directory:
        return
    os.makedirs(directory, exist_ok=True)
    descriptor, path = tempfile.mkstemp(
        prefix="execution-", suffix=".json", dir=directory,
    )
    complete = False
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            descriptor = -1
            json.dump(dict(record), handle, sort_keys=True)
            handle.write("\n")
            complete = True
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if not complete:
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass


def collect_execution_events(args: Any) -> list[dict[str, Any]]:
    """Move completed worker events into the parent manifest exactly once."""

    directory = getattr(args, "aligner_diagnostics_dir", None)
    if not directory:
        return []
    records = []
    for path in sorted(Path(directory).glob("execution-*.json")):
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(value, dict):
                records.append(value)
        except (OSError, UnicodeError, json.JSONDecodeError):
            continue
        finally:
            try:
                path.unlink()
            except FileNotFoundError:
                pass
    records.sort(key=lambda value: (
        int(value.get("event_time_ns", 0)), int(value.get("process_id", 0)),
    ))
    manifest = getattr(args, "_run_manifest", None)
    if manifest is not None:
        for record in records:
            manifest.record_aligner_execution(record)
    else:
        args._aligner_execution_events = (
            list(getattr(args, "_aligner_execution_events", [])) + records
        )
    return records
