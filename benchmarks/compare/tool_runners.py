"""Profiled standalone runners for Liftoff, miniprot, and LiftOn.

Each runner writes its output under ``work/<id>/tools/<tool>/`` and a profile
JSON (wall-clock + peak RSS) under ``tools/logs/``.  Reuse is fail-closed: a
versioned cache manifest binds the exact request and hashes both the output and
profile.  The legacy ``.done/<tool>.done`` sentinel is retained for compatibility
but is never sufficient cache evidence by itself.  LiftOn reuses the standalone
Liftoff + miniprot outputs via ``-L`` / ``-M``.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Mapping

from .profiling import run_profiled


TOOL_CACHE_SCHEMA_VERSION = 1
TOOL_CACHE_BUILDER_VERSION = "content-addressed-tool-runner-v1"
TOOL_CACHE_MANIFEST_NAME = "cache.manifest.json"

# These variables can alter executable discovery, Python imports, native
# libraries, locale-sensitive output, threading, or temporary-file behaviour.
# Job/scheduler metadata is intentionally excluded: it cannot change a tool's
# result and would make an otherwise identical cache unusable on every job.
_RELEVANT_ENVIRONMENT_NAMES = frozenset({
    "CONDA_PREFIX",
    "HOME",
    "LANG",
    "LC_ALL",
    "LD_LIBRARY_PATH",
    "MPLCONFIGDIR",
    "PATH",
    "PYTHONHASHSEED",
    "PYTHONDONTWRITEBYTECODE",
    "PYTHONNOUSERSITE",
    "PYTHONPATH",
    "TMPDIR",
    "VIRTUAL_ENV",
})
_RELEVANT_ENVIRONMENT_PREFIXES = (
    "LIFTON_",
    "LIFTOFF_",
    "MINIMAP2_",
    "MINIPROT_",
    "MKL_",
    "NUMEXPR_",
    "OMP_",
    "OPENBLAS_",
)


def _compose_env(tools: dict) -> dict:
    """PATH that lets every subprocess find minimap2/miniprot/samtools, plus
    PYTHONNOUSERSITE to avoid user-site bleed across conda envs."""
    extra = tools.get("extra_path", "")
    env = {
        "PATH": f"{extra}:{os.environ.get('PATH', '')}" if extra else os.environ.get("PATH", ""),
        "PYTHONNOUSERSITE": "1",
    }
    return env


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(document: Any) -> str:
    payload = json.dumps(
        document,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _absolute_path(path: str | Path) -> Path:
    return Path(os.path.abspath(Path(path).expanduser()))


def _file_record(path: str | Path, *, require_nonempty: bool = True) -> dict[str, Any]:
    """Return a stable content fingerprint, rejecting a concurrently changed file."""

    absolute = _absolute_path(path)
    try:
        before = absolute.stat()
    except OSError as exc:
        raise RuntimeError(f"cache file is missing or unreadable: {absolute}") from exc
    if not absolute.is_file() or (require_nonempty and before.st_size <= 0):
        qualifier = "missing or empty" if require_nonempty else "missing"
        raise RuntimeError(f"cache file is {qualifier}: {absolute}")
    digest = _sha256_file(absolute)
    try:
        after = absolute.stat()
    except OSError as exc:
        raise RuntimeError(f"cache file disappeared while hashing: {absolute}") from exc
    stable_fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns")
    if any(getattr(before, field) != getattr(after, field) for field in stable_fields):
        raise RuntimeError(f"cache file changed while hashing: {absolute}")
    return {
        "path": str(absolute),
        "size": after.st_size,
        "sha256": digest,
    }


def _normalise_path_list(value: str) -> str:
    entries = []
    for entry in value.split(os.pathsep):
        entries.append(str(_absolute_path(entry or ".")))
    return os.pathsep.join(entries)


def _normalised_relevant_environment(overlay: Mapping[str, str]) -> dict[str, str]:
    effective = {**os.environ, **{str(key): str(value) for key, value in overlay.items()}}
    keys = {
        key for key in effective
        if key in _RELEVANT_ENVIRONMENT_NAMES
        or key.startswith(_RELEVANT_ENVIRONMENT_PREFIXES)
    }
    result = {key: effective[key] for key in sorted(keys)}
    if "PATH" in result:
        result["PATH"] = _normalise_path_list(result["PATH"])
    if "PYTHONPATH" in result:
        result["PYTHONPATH"] = _normalise_path_list(result["PYTHONPATH"])
    return result


def _resolve_executable(executable: str | Path, env: Mapping[str, str]) -> Path:
    configured = str(executable)
    has_separator = os.sep in configured or (
        os.altsep is not None and os.altsep in configured
    )
    candidate = configured if has_separator else shutil.which(
        configured, path=env.get("PATH", os.environ.get("PATH", "")),
    )
    if candidate is None:
        raise RuntimeError(f"required executable is not on PATH: {configured}")
    resolved = Path(candidate).expanduser().resolve()
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        raise RuntimeError(f"required executable is not executable: {resolved}")
    return resolved


def _normalise_version_output(stdout: str, stderr: str) -> str:
    selected = stdout.strip() or stderr.strip()
    return "\n".join(line.rstrip() for line in selected.splitlines()).strip()


def _executable_record(
    executable: str | Path,
    env: Mapping[str, str],
) -> dict[str, Any]:
    resolved = _resolve_executable(executable, env)
    version_argv = [str(resolved), "--version"]
    try:
        result = subprocess.run(
            version_argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=10,
            check=False,
            env={**os.environ, **env},
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"could not probe executable version: {resolved}: {exc}") from exc
    version = _normalise_version_output(result.stdout or "", result.stderr or "")
    if result.returncode != 0 or not version:
        raise RuntimeError(
            f"executable version probe failed for {resolved} "
            f"(exit {result.returncode})"
        )
    record = _file_record(resolved)
    return {
        "configured": str(executable),
        "path": record["path"],
        "size": record["size"],
        "sha256": record["sha256"],
        "version": version,
        "version_argv": version_argv,
        "version_exit_code": result.returncode,
    }


def _tool_request(
    *,
    tool_name: str,
    executable: str | Path,
    argv: list[str],
    inputs: Mapping[str, str | Path],
    env: Mapping[str, str],
    threads: int,
    mode: str,
    options: Mapping[str, Any],
) -> dict[str, Any]:
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise RuntimeError("tool threads must be a positive integer")
    request = {
        "schema_version": TOOL_CACHE_SCHEMA_VERSION,
        "builder": {
            "id": TOOL_CACHE_BUILDER_VERSION,
            "source_sha256": _sha256_file(Path(__file__).resolve()),
        },
        "tool": tool_name,
        "argv": [str(value) for value in argv],
        "environment": _normalised_relevant_environment(env),
        "inputs": {
            label: _file_record(path)
            for label, path in sorted(inputs.items())
        },
        "executable": _executable_record(executable, env),
        "configuration": {
            "threads": threads,
            "mode": mode,
            "options": dict(sorted(options.items())),
        },
    }
    return request


def _atomic_write_json(path: Path, document: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(document, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _cache_paths(
    work_dir: Path,
    tool_name: str,
    mode: str,
) -> tuple[Path, Path, Path, Path]:
    if tool_name not in {"liftoff", "miniprot", "lifton", "lifton_optimize"}:
        raise ValueError(f"unsupported standalone tool cache: {tool_name!r}")
    if mode not in {"genes", "allfeat"}:
        raise ValueError(f"unsupported standalone tool mode: {mode!r}")
    if tool_name == "miniprot" and mode != "genes":
        raise ValueError("miniprot has only the shared genes-mode cache")
    tools_root, suffix = _mode_dirs(Path(work_dir), mode)
    out_dir = tools_root / tool_name
    profile_path = tools_root / "logs" / f"{tool_name}.profile.json"
    done = Path(work_dir) / ".done" / f"{tool_name}{suffix}.done"
    return out_dir / TOOL_CACHE_MANIFEST_NAME, done, out_dir, profile_path


def _record_matches(path: Path, record: Any) -> bool:
    if not isinstance(record, dict) or set(record) != {"path", "size", "sha256"}:
        return False
    if record.get("path") != str(_absolute_path(path)):
        return False
    try:
        return _file_record(path) == record
    except (OSError, RuntimeError):
        return False


def _request_is_live(request: Any) -> bool:
    if not isinstance(request, dict) or set(request) != {
        "schema_version", "builder", "tool", "argv", "environment",
        "inputs", "executable", "configuration",
    }:
        return False
    if request.get("schema_version") != TOOL_CACHE_SCHEMA_VERSION:
        return False
    builder = request.get("builder")
    if builder != {
        "id": TOOL_CACHE_BUILDER_VERSION,
        "source_sha256": _sha256_file(Path(__file__).resolve()),
    }:
        return False
    argv = request.get("argv")
    if not isinstance(argv, list) or not argv or any(not isinstance(x, str) for x in argv):
        return False
    environment = request.get("environment")
    if not isinstance(environment, dict) or any(
        not isinstance(key, str) or not isinstance(value, str)
        for key, value in environment.items()
    ):
        return False
    inputs = request.get("inputs")
    if not isinstance(inputs, dict) or not inputs:
        return False
    for record in inputs.values():
        if not isinstance(record, dict) or set(record) != {"path", "size", "sha256"}:
            return False
        try:
            if _file_record(record["path"]) != record:
                return False
        except (KeyError, OSError, RuntimeError, TypeError, ValueError):
            return False
    executable = request.get("executable")
    if not isinstance(executable, dict):
        return False
    try:
        live_executable = _executable_record(executable["configured"], environment)
    except (KeyError, OSError, RuntimeError, TypeError, ValueError):
        return False
    if live_executable != executable or argv[0] != executable["configured"]:
        return False
    configuration = request.get("configuration")
    return (
        isinstance(configuration, dict)
        and set(configuration) == {"threads", "mode", "options"}
        and isinstance(configuration["threads"], int)
        and not isinstance(configuration["threads"], bool)
        and configuration["threads"] >= 1
        and isinstance(configuration["mode"], str)
        and isinstance(configuration["options"], dict)
    )


def _verified_cached_profile(
    manifest_path: Path,
    done: Path,
    output_path: Path,
    profile_path: Path,
    *,
    expected_request: Mapping[str, Any] | None = None,
) -> dict[str, Any] | None:
    if not done.is_file() or not manifest_path.is_file():
        return None
    try:
        document = json.loads(manifest_path.read_text(encoding="utf-8"))
        request = document["request"]
        request_sha256 = _canonical_sha256(request)
        if (
            not isinstance(document, dict)
            or set(document) != {
                "schema_version", "builder_version", "request_sha256",
                "request", "artifacts",
            }
            or document["schema_version"] != TOOL_CACHE_SCHEMA_VERSION
            or document["builder_version"] != TOOL_CACHE_BUILDER_VERSION
            or document["request_sha256"] != request_sha256
            or (expected_request is not None and request != expected_request)
            or not _request_is_live(request)
        ):
            return None
        artifacts = document["artifacts"]
        if not isinstance(artifacts, dict) or set(artifacts) != {"output", "profile"}:
            return None
        if not _record_matches(output_path, artifacts["output"]):
            return None
        if not _record_matches(profile_path, artifacts["profile"]):
            return None
        profile = json.loads(profile_path.read_text(encoding="utf-8"))
        if not isinstance(profile, dict) or profile.get("exit_code") != 0:
            return None
        return profile
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None


def verify_cached_tool_run(
    work_dir: str | Path,
    tool_name: str,
    mode: str = "genes",
) -> bool:
    """Verify a standalone cache against all live content it fingerprints.

    This is intentionally usable without the benchmark registry so acquisition
    preflight can reject tampered or incomplete caches.  A runner additionally
    compares the manifest with a newly constructed request before reuse, which
    catches changed command-line options and environment.
    """

    manifest_path, done, out_dir, profile_path = _cache_paths(
        Path(work_dir), tool_name, mode,
    )
    output_name = "lifton.gff3" if tool_name.startswith("lifton") else f"{tool_name}.gff3"
    return _verified_cached_profile(
        manifest_path,
        done,
        out_dir / output_name,
        profile_path,
    ) is not None


def _prepare_cache_rebuild(manifest_path: Path, done: Path) -> None:
    # Removal happens before launching the external process.  A failed rebuild
    # can therefore never leave a previous successful manifest/sentinel that
    # appears to bless a partial new output.
    manifest_path.unlink(missing_ok=True)
    done.unlink(missing_ok=True)


def _write_cache_manifest(
    manifest_path: Path,
    done: Path,
    request: Mapping[str, Any],
    output_path: Path,
    profile_path: Path,
) -> None:
    document = {
        "schema_version": TOOL_CACHE_SCHEMA_VERSION,
        "builder_version": TOOL_CACHE_BUILDER_VERSION,
        "request_sha256": _canonical_sha256(request),
        "request": request,
        "artifacts": {
            "output": _file_record(output_path),
            "profile": _file_record(profile_path),
        },
    }
    _atomic_write_json(manifest_path, document)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("ok\n", encoding="ascii")


def _save_profile(pr, path: Path) -> dict:
    d = dataclasses.asdict(pr)
    path.write_text(json.dumps(d, indent=2))
    return d


def _clean_input_dbs(*gffs) -> None:
    """Unlink stale gffutils '<gff>_db' SQLite siblings so the next lifton run
    rebuilds the input DBs fresh. Re-running lifton on cached -L/-M inputs
    otherwise hits a stale/partial sibling ("table features already exists") or,
    worse, a db built from a now-superseded GFF (e.g. an old protein-space
    miniprot that silently suppresses the merge). Deterministic — the inputs are
    unchanged, lifton rebuilds. All paths live under work/<id>/subset|tools/
    (salz3), never the read-only salz2 originals."""
    for g in gffs:
        dbf = Path(str(g) + "_db")
        if dbf.exists():
            dbf.unlink()


def _mode_dirs(work_dir: Path, mode: str):
    """Return (tools_root, done_suffix) for a feature mode. "genes" keeps the
    legacy paths (tools/, *.done) so already-run results stay byte-identical;
    other modes use parallel namespaces (tools_<mode>/, *__<mode>.done)."""
    if mode == "genes":
        return work_dir / "tools", ""
    return work_dir / f"tools_{mode}", f"__{mode}"


def run_liftoff(bench: dict, manifest: dict, work_dir: Path, tools: dict,
                threads: int, force: bool, log=print, mode="genes",
                types_file=None) -> dict:
    tools_root, suf = _mode_dirs(work_dir, mode)
    out_dir = tools_root / "liftoff"
    log_dir = tools_root / "logs"
    done = work_dir / ".done" / f"liftoff{suf}.done"
    out_gff = out_dir / "liftoff.gff3"
    prof_path = log_dir / "liftoff.profile.json"
    cache_manifest = out_dir / TOOL_CACHE_MANIFEST_NAME
    p = manifest["paths"]
    inter = out_dir / "intermediate"
    raw_out = out_dir / "liftoff_raw.gff3"
    argv = [tools["liftoff_bin"]]
    if mode != "genes":
        if types_file is None:
            raise RuntimeError(f"liftoff mode {mode!r} requires a feature-types file")
        argv += ["-f", str(types_file)]   # lift all listed top-level types
    argv += [
        "-g", p["ref_gff"],
        "-o", str(raw_out),
        "-u", str(out_dir / "unmapped_features.txt"),
        "-dir", str(inter),
        "-p", str(threads),
        "-copies", "-polish",
        p["tgt_fa"], p["ref_fa"],
    ]
    env = _compose_env(tools)
    inputs = {
        "reference_annotation": p["ref_gff"],
        "reference_fasta": p["ref_fa"],
        "target_fasta": p["tgt_fa"],
    }
    if types_file is not None:
        inputs["feature_types"] = types_file
    request = _tool_request(
        tool_name="liftoff",
        executable=tools["liftoff_bin"],
        argv=argv,
        inputs=inputs,
        env=env,
        threads=threads,
        mode=mode,
        options={
            "copies": True,
            "polish": True,
            "feature_types": mode != "genes",
        },
    )
    if not force:
        cached = _verified_cached_profile(
            cache_manifest, done, out_gff, prof_path,
            expected_request=request,
        )
        if cached is not None:
            log(f"  [liftoff:{mode}] verified cache")
            return cached
    _prepare_cache_rebuild(cache_manifest, done)
    out_dir.mkdir(parents=True, exist_ok=True)
    polished = Path(str(raw_out) + "_polished")
    for stale in (raw_out, polished, out_gff):
        stale.unlink(missing_ok=True)
    pr = run_profiled(argv, label="liftoff", log_dir=log_dir,
                      env=env, log=log)
    if pr.exit_code != 0:
        raise RuntimeError(f"liftoff failed (exit {pr.exit_code}); see {pr.stderr_path}")
    # -polish makes Liftoff emit <out>_polished; normalize to liftoff.gff3.
    src = polished if polished.exists() else raw_out
    if not src.exists():
        raise RuntimeError(f"liftoff produced no output ({raw_out} / {polished})")
    shutil.copyfile(src, out_gff)
    prof = _save_profile(pr, prof_path)
    _write_cache_manifest(cache_manifest, done, request, out_gff, prof_path)
    log(f"  [liftoff:{mode}] {pr.wall_clock_seconds:.1f}s, {pr.peak_rss_mb:.0f} MB")
    return prof


def run_miniprot(bench: dict, manifest: dict, work_dir: Path, tools: dict,
                 threads: int, force: bool, log=print) -> dict:
    out_dir = work_dir / "tools" / "miniprot"
    log_dir = work_dir / "tools" / "logs"
    done = work_dir / ".done" / "miniprot.done"
    out_gff = out_dir / "miniprot.gff3"
    prof_path = log_dir / "miniprot.profile.json"
    cache_manifest = out_dir / TOOL_CACHE_MANIFEST_NAME
    p = manifest["paths"]
    argv = [
        tools["miniprot_bin"], "-t", str(threads), "--gff-only",
        p["tgt_fa"], p["ref_faa"],
    ]
    env = _compose_env(tools)
    request = _tool_request(
        tool_name="miniprot",
        executable=tools["miniprot_bin"],
        argv=argv,
        inputs={
            "reference_proteins": p["ref_faa"],
            "target_fasta": p["tgt_fa"],
        },
        env=env,
        threads=threads,
        mode="genes",
        options={"gff_only": True},
    )
    if not force:
        cached = _verified_cached_profile(
            cache_manifest, done, out_gff, prof_path,
            expected_request=request,
        )
        if cached is not None:
            log("  [miniprot] verified cache")
            return cached
    _prepare_cache_rebuild(cache_manifest, done)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_gff.unlink(missing_ok=True)
    # run_profiled redirects stdout to <label>.stdout.log; miniprot --gff-only
    # writes GFF to stdout, so the GFF is that stdout file.
    pr = run_profiled(argv, label="miniprot", log_dir=log_dir,
                      env=env, log=log)
    if pr.exit_code != 0:
        raise RuntimeError(f"miniprot failed (exit {pr.exit_code}); see {pr.stderr_path}")
    err = Path(pr.stderr_path).read_text(errors="replace") if pr.stderr_path else ""
    if "ERROR" in err.upper():
        raise RuntimeError(f"miniprot reported ERROR on stderr; see {pr.stderr_path}")
    shutil.copyfile(pr.stdout_path, out_gff)
    if out_gff.stat().st_size == 0:
        raise RuntimeError("miniprot produced empty output")
    prof = _save_profile(pr, prof_path)
    _write_cache_manifest(cache_manifest, done, request, out_gff, prof_path)
    log(f"  [miniprot] {pr.wall_clock_seconds:.1f}s, {pr.peak_rss_mb:.0f} MB")
    return prof


def run_lifton(bench: dict, manifest: dict, work_dir: Path, tools: dict,
               threads: int, force: bool, log=print, mode="genes",
               types_file=None, optimize=False) -> dict:
    # The optimize variant is a parallel, separately-cached run (own dir +
    # sentinel + profile) so the byte-frozen default lifton output is never
    # disturbed. It is genes-mode only (the merge is a coding-protein refinement;
    # allfeat adds only non-coding top-level types it never touches).
    name = "lifton_optimize" if optimize else "lifton"
    tools_root, suf = _mode_dirs(work_dir, mode)
    out_dir = tools_root / name
    log_dir = tools_root / "logs"
    done = work_dir / ".done" / f"{name}{suf}.done"
    out_gff = out_dir / "lifton.gff3"   # keep the filename; the dir disambiguates
    prof_path = log_dir / f"{name}.profile.json"
    cache_manifest = out_dir / TOOL_CACHE_MANIFEST_NAME
    p = manifest["paths"]
    liftoff_gff = tools_root / "liftoff" / "liftoff.gff3"          # mode-specific
    miniprot_gff = work_dir / "tools" / "miniprot" / "miniprot.gff3"  # shared, mode-independent
    if not verify_cached_tool_run(work_dir, "liftoff", mode):
        raise RuntimeError(
            f"{name}:{mode} needs a verified Liftoff cache; "
            f"run liftoff mode={mode} first"
        )
    if not verify_cached_tool_run(work_dir, "miniprot"):
        raise RuntimeError(
            f"{name}:{mode} needs the verified shared genes-mode miniprot cache; "
            "run genes mode first"
        )
    argv = [
        tools["lifton_bin"], "-t", str(threads), "-copies",
        "-ad", bench.get("annotation_database", "RefSeq"),
        "-g", p["ref_gff"],
        "-L", str(liftoff_gff),
        "-M", str(miniprot_gff),
        "-o", str(out_gff),
    ]
    # Parity with liftoff (-p threads) and miniprot (-t threads): unlock LiftOn's
    # threaded Step 7. --native routes per-locus workers through the Phase-17b
    # materialised proxy DBs (gffutils-safe), and --locus-pipeline enables the
    # ThreadPoolExecutor. The 24-cell matrix proves this is byte-identical to the
    # serial default; -L/-M are still honoured (the native aligners never run).
    if threads and threads > 1:
        argv += ["--native", "--locus-pipeline"]
    if optimize:
        argv += ["--optimize"]   # opt-in best-of-outcome merge + Case-1 fix
    if mode != "genes":
        if types_file is None:
            raise RuntimeError(f"lifton mode {mode!r} requires a feature-types file")
        argv += ["-f", str(types_file)]   # process all listed top-level types
    argv += [p["tgt_fa"], p["ref_fa"]]
    env = _compose_env(tools)
    inputs = {
        "liftoff_gff": liftoff_gff,
        "miniprot_gff": miniprot_gff,
        "reference_annotation": p["ref_gff"],
        "reference_fasta": p["ref_fa"],
        "target_fasta": p["tgt_fa"],
    }
    if types_file is not None:
        inputs["feature_types"] = types_file
    request = _tool_request(
        tool_name=name,
        executable=tools["lifton_bin"],
        argv=argv,
        inputs=inputs,
        env=env,
        threads=threads,
        mode=mode,
        options={
            "annotation_database": bench.get("annotation_database", "RefSeq"),
            "copies": True,
            "feature_types": mode != "genes",
            "locus_pipeline": threads > 1,
            "native": threads > 1,
            "optimize": bool(optimize),
        },
    )
    if not force:
        cached = _verified_cached_profile(
            cache_manifest, done, out_gff, prof_path,
            expected_request=request,
        )
        if cached is not None:
            log(f"  [{name}:{mode}] verified cache")
            return cached
    _prepare_cache_rebuild(cache_manifest, done)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_gff.unlink(missing_ok=True)
    # Any run that proceeds here is non-cached and reuses the verified -L/-M
    # inputs. Clear stale gffutils '<gff>_db' siblings so lifton rebuilds from
    # the current GFFs. run_all_tools invokes variants sequentially, so the
    # shared input DBs are never raced.
    _clean_input_dbs(p["ref_gff"], liftoff_gff, miniprot_gff)
    pr = run_profiled(argv, label=name, log_dir=log_dir,
                      env=env, log=log)
    if pr.exit_code != 0:
        raise RuntimeError(f"{name} failed (exit {pr.exit_code}); see {pr.stderr_path}")
    if not out_gff.exists() or out_gff.stat().st_size == 0:
        raise RuntimeError(f"{name} produced no/empty output {out_gff}")
    prof = _save_profile(pr, prof_path)
    _write_cache_manifest(cache_manifest, done, request, out_gff, prof_path)
    log(f"  [{name}:{mode}] {pr.wall_clock_seconds:.1f}s, {pr.peak_rss_mb:.0f} MB")
    return prof


def run_all_tools(bench, manifest, work_dir, tools, threads, force, log=print,
                  mode="genes") -> dict:
    profiles = {}
    types_file = None
    if mode != "genes":
        from . import feature_types
        types_file = feature_types.write_types_file(manifest, work_dir, log)
    profiles["liftoff"] = run_liftoff(bench, manifest, work_dir, tools, threads,
                                      force, log, mode=mode, types_file=types_file)
    if mode == "genes":
        profiles["miniprot"] = run_miniprot(bench, manifest, work_dir, tools,
                                            threads, force, log)
    else:
        # miniprot is protein-only / mode-independent: reuse the genes-mode run
        # (do NOT re-run it); load its profile so the allfeat eval/report has its timing.
        if not verify_cached_tool_run(work_dir, "miniprot"):
            raise RuntimeError(
                f"{bench['id']}: mode {mode!r} needs the verified shared "
                "genes-mode miniprot — run genes mode first"
            )
        mp_prof = work_dir / "tools" / "logs" / "miniprot.profile.json"
        profiles["miniprot"] = json.loads(mp_prof.read_text())
    profiles["lifton"] = run_lifton(bench, manifest, work_dir, tools, threads,
                                    force, log, mode=mode, types_file=types_file)
    # 4th variant: LiftOn with --optimize (best-of-outcome merge). Genes-mode
    # only; runs after the default lifton (sequential → no shared-DB race).
    if mode == "genes":
        profiles["lifton_optimize"] = run_lifton(
            bench, manifest, work_dir, tools, threads, force, log,
            mode=mode, types_file=types_file, optimize=True)
    return profiles
