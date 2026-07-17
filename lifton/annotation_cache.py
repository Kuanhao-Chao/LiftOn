"""Content-aware cache manifests for annotation databases.

Annotation databases are derived artifacts.  A database file by itself does
not say which source bytes or parser options produced it, so it must never be
reused without its matching manifest.
"""

from __future__ import annotations

import hashlib
import json
import os
import uuid
from pathlib import Path
from typing import Any, Mapping


MANIFEST_FORMAT_VERSION = 1
PARSER_CONTRACT_VERSION = 1


def manifest_path_for(db_path: str) -> str:
    """Return the JSON sidecar path for an annotation database."""
    return os.fspath(db_path) + ".manifest.json"


def _json_value(value: Any) -> Any:
    """Convert parser settings into a stable, JSON-compatible value."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if isinstance(value, Mapping):
        return {
            str(key): _json_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, (set, frozenset)):
        normalized = [_json_value(item) for item in value]
        return sorted(normalized, key=lambda item: json.dumps(item, sort_keys=True))
    if callable(value):
        module = getattr(value, "__module__", type(value).__module__)
        name = getattr(value, "__qualname__", getattr(value, "__name__", type(value).__qualname__))
        return {"callable": f"{module}.{name}"}
    return {"type": f"{type(value).__module__}.{type(value).__qualname__}", "repr": repr(value)}


def source_fingerprint(source_path: str) -> dict[str, Any]:
    """Hash a stable snapshot of *source_path* and return its size.

    The before/after metadata check prevents recording a digest assembled
    while another process was replacing or modifying the annotation.
    """
    path = os.fspath(source_path)
    before = os.stat(path)
    digest = hashlib.sha256()
    bytes_read = 0
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            bytes_read += len(chunk)
    after = os.stat(path)
    before_key = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    after_key = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if before_key != after_key or bytes_read != after.st_size:
        raise OSError(f"Annotation changed while it was being fingerprinted: {path}")
    return {"sha256": digest.hexdigest(), "size": bytes_read}


def schema_digest(schema: str) -> str:
    """Return a compact version identifier for a textual database schema."""
    return hashlib.sha256(schema.encode("utf-8")).hexdigest()


def build_manifest(
    source_path: str,
    *,
    backend: str,
    backend_version: str,
    schema_version: str,
    tool_version: str,
    settings: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the complete expected manifest for a derived database."""
    return {
        "manifest_format_version": MANIFEST_FORMAT_VERSION,
        "source": source_fingerprint(source_path),
        "backend": {
            "name": backend,
            "version": str(backend_version),
            "schema_version": str(schema_version),
        },
        "tool": {"name": "LiftOn", "version": str(tool_version)},
        "settings": _json_value(settings),
    }


def load_manifest(db_path: str) -> dict[str, Any] | None:
    """Load a cache manifest, returning ``None`` for missing/corrupt data."""
    try:
        with open(manifest_path_for(db_path), "r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def cache_matches(db_path: str, expected: Mapping[str, Any]) -> bool:
    """Return whether a non-empty database and exact manifest both exist."""
    try:
        if not os.path.isfile(db_path) or os.path.getsize(db_path) == 0:
            return False
    except OSError:
        return False
    return load_manifest(db_path) == dict(expected)


def write_manifest_atomic(db_path: str, manifest: Mapping[str, Any]) -> None:
    """Publish *manifest* using an atomic same-directory rename."""
    target = manifest_path_for(db_path)
    parent = os.path.dirname(os.path.abspath(target))
    Path(parent).mkdir(parents=True, exist_ok=True)
    temporary = f"{target}.tmp.{os.getpid()}.{uuid.uuid4().hex}"
    try:
        with open(temporary, "x", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, target)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def temporary_db_path(db_path: str) -> str:
    """Return an unused same-directory path suitable for atomic publish."""
    target = os.path.abspath(os.fspath(db_path))
    parent = os.path.dirname(target)
    name = os.path.basename(target)
    return os.path.join(parent, f".{name}.tmp.{os.getpid()}.{uuid.uuid4().hex}")
