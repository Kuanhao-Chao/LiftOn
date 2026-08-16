"""Recompute release-evaluation metrics without mutating sealed run roots."""
from __future__ import annotations

import argparse
import copy
import datetime as dt
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from . import build_controller, evaluator, id_mapping, release_evaluation
from . import release_report


OVERLAY_SCHEMA_VERSION = 1
OVERLAY_KIND = "lifton_release_evaluation_overlay"
RECORD_KIND = "lifton_release_evaluation_overlay_record"
RUNTIME_PACKAGES = (
    "biopython",
    "gffutils",
    "intervaltree",
    "numpy",
    "parasail",
    "pyfaidx",
    "pysam",
)


def _utc_now() -> str:
    return (
        dt.datetime.now(dt.timezone.utc)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _artifact_record(path: Path) -> dict[str, Any]:
    resolved = Path(path).resolve()
    stat = resolved.stat()
    if not resolved.is_file() or stat.st_size <= 0:
        raise ValueError(f"artifact is missing or empty: {resolved}")
    digest = release_evaluation.sha256_file(resolved)
    after = resolved.stat()
    if stat.st_size != after.st_size or stat.st_mtime_ns != after.st_mtime_ns:
        raise ValueError(f"artifact changed while it was hashed: {resolved}")
    return {
        "path": str(resolved),
        "size": after.st_size,
        "sha256": digest,
    }


def _verify_artifact_record(record: Any, *, label: str) -> dict[str, Any]:
    if not isinstance(record, Mapping):
        raise ValueError(f"{label} artifact record is malformed")
    path = record.get("path")
    if not isinstance(path, str) or not Path(path).is_absolute():
        raise ValueError(f"{label} artifact path must be absolute")
    observed = _artifact_record(Path(path))
    if dict(record) != observed:
        raise ValueError(f"{label} artifact no longer matches its SHA-256 record")
    return observed


def _tooling_records() -> dict[str, dict[str, Any]]:
    repo_root = Path(build_controller.REPO_ROOT).resolve()
    files = {
        "build_controller": Path(build_controller.__file__).resolve(),
        "release_rescore": Path(__file__).resolve(),
        "release_evaluation": Path(release_evaluation.__file__).resolve(),
        "release_report": Path(release_report.__file__).resolve(),
        "evaluator": Path(evaluator.__file__).resolve(),
        "id_mapping": Path(id_mapping.__file__).resolve(),
    }
    package_root = repo_root / "lifton"
    files.update({
        f"lifton/{path.relative_to(package_root).as_posix()}": path
        for path in sorted(package_root.rglob("*.py"))
        if path.stat().st_size > 0
    })
    return {
        name: _artifact_record(path)
        for name, path in sorted(files.items())
    }


def _runtime_record() -> dict[str, Any]:
    packages = {}
    for name in RUNTIME_PACKAGES:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "python": sys.version,
        "python_executable": str(Path(sys.executable).resolve()),
        "platform": platform.platform(),
        "packages": packages,
    }


class _SourceVerifier:
    def __init__(self) -> None:
        self._records: dict[Path, dict[str, Any]] = {}

    def verify(self, record: Any, *, label: str) -> dict[str, Any]:
        if not isinstance(record, Mapping):
            raise ValueError(f"{label} fingerprint is malformed")
        raw_path = record.get("path")
        if not isinstance(raw_path, str) or not Path(raw_path).is_absolute():
            raise ValueError(f"{label} fingerprint path must be absolute")
        path = Path(raw_path).resolve()
        observed = self._records.get(path)
        if observed is None:
            observed = _artifact_record(path)
            self._records[path] = observed
        expected = {
            "path": str(path),
            "size": record.get("size"),
            "sha256": record.get("sha256"),
        }
        if expected != observed:
            raise ValueError(f"{label} no longer matches sealed input evidence")
        return dict(observed)


def _safe_component(value: str) -> str:
    normalized = "".join(
        character if character.isalnum() or character in "-_" else "_"
        for character in value
    ).strip("_")
    if not normalized:
        raise ValueError(f"unsafe empty path component derived from {value!r}")
    return normalized[:120]


def _pair_key(pair: Mapping[str, Any]) -> tuple[str, str, int]:
    panel = pair.get("panel")
    benchmark = pair.get("benchmark")
    repetition = pair.get("repetition")
    if (
        panel not in {"subset", "full", "e2e"}
        or not isinstance(benchmark, str)
        or not benchmark
        or not isinstance(repetition, int)
        or isinstance(repetition, bool)
        or repetition <= 0
    ):
        raise ValueError("paired result identity is malformed")
    return panel, benchmark, repetition


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _isolated_link(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.symlink_to(Path(source).resolve())
    return destination


def _isolated_fasta(source: Path, destination: Path) -> Path:
    linked = _isolated_link(source, destination)
    for canonical, isolated in (
        (Path(str(source) + ".fai"), Path(str(destination) + ".fai")),
        (Path(str(source) + ".gzi"), Path(str(destination) + ".gzi")),
        (source.with_suffix(".dict"), destination.with_suffix(".dict")),
    ):
        if canonical.is_file():
            shutil.copy2(canonical, isolated)
    return linked


def _group_signature(pair: Mapping[str, Any]) -> str:
    inputs = pair.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ValueError("paired result lacks input fingerprints")
    selected = {
        name: inputs.get(name)
        for name in ("ref_gff", "ref_fa", "tgt_fa")
    }
    if any(not isinstance(value, Mapping) for value in selected.values()):
        raise ValueError("paired result lacks reference/target scoring inputs")
    return _sha256_text(_canonical_json(selected))


def _summary_manifest(pair: Mapping[str, Any]) -> dict[str, Any]:
    candidate = (pair.get("versions") or {}).get("candidate")
    summary = candidate.get("summary") if isinstance(candidate, Mapping) else None
    if not isinstance(summary, Mapping):
        raise ValueError("candidate evaluator summary is missing")
    return {
        "id": pair["benchmark"],
        "species": str(summary.get("species", pair["benchmark"])),
        "cross_species": bool(summary.get("cross_species", False)),
        "miniprot_target_space": "transcript",
        "protein_acc_to_mrna": {},
    }


def _output_record(
    pair_path: Path,
    pair: Mapping[str, Any],
    label: str,
) -> dict[str, Any]:
    version = (pair.get("versions") or {}).get(label)
    if not isinstance(version, Mapping):
        raise ValueError(f"{pair_path}: {label} version evidence is malformed")
    raw_path = version.get("output_gff")
    if not isinstance(raw_path, str) or not Path(raw_path).is_absolute():
        raise ValueError(f"{pair_path}: {label} output path is malformed")
    path = Path(raw_path).resolve()
    fingerprints = release_evaluation.gff3_fingerprints(path)
    if fingerprints != version.get("fingerprints"):
        raise ValueError(
            f"{pair_path}: {label} output changed after sealed evaluation"
        )
    return {
        "path": str(path),
        "size": path.stat().st_size,
        "fingerprints": fingerprints,
    }


def _score_group(
    rows: Sequence[tuple[Path, Mapping[str, Any], dict[str, Any]]],
    *,
    staging: Path,
    tooling_sha256: str,
    threads: int,
) -> list[dict[str, Any]]:
    _first_path, first_pair, first_source = rows[0]
    panel, benchmark, _ = _pair_key(first_pair)
    signature = _group_signature(first_pair)
    scratch = staging / ".scratch" / f"{panel}-{_safe_component(benchmark)}-{signature[:12]}"
    scratch.mkdir(parents=True, exist_ok=False)
    inputs = first_source["inputs"]
    ref_gff = _isolated_link(
        Path(inputs["ref_gff"]["path"]), scratch / "reference.gff3",
    )
    ref_fa = _isolated_fasta(
        Path(inputs["ref_fa"]["path"]), scratch / "reference.fa",
    )
    tgt_fa = _isolated_fasta(
        Path(inputs["tgt_fa"]["path"]), scratch / "target.fa",
    )
    reference, reference_index = evaluator.build_reference(
        str(ref_gff),
        str(ref_fa),
        log=lambda message: print(message, flush=True),
    )
    reference_id_index = release_evaluation._declared_id_index(ref_gff)
    manifest = _summary_manifest(first_pair)
    manifest_rows = []
    for pair_path, pair, source in rows:
        row_panel, row_benchmark, repetition = _pair_key(pair)
        if (row_panel, row_benchmark) != (panel, benchmark):
            raise AssertionError("rescore group identity changed")
        pair_sha = source["pair_result"]["sha256"]
        relative = Path("records") / panel / _safe_component(benchmark) / (
            f"repetition-{repetition:02d}-{pair_sha[:12]}"
        )
        record_dir = staging / relative
        evaluation_dir = record_dir / "evaluation"
        evaluation_dir.mkdir(parents=True, exist_ok=False)
        versions = {}
        for label in ("candidate", "reference"):
            output_source = Path(source["outputs"][label]["path"])
            output = _isolated_link(
                output_source,
                scratch / f"repetition-{repetition:02d}-{label}.gff3",
            )
            original_version = pair["versions"][label]
            summary = dict(evaluator.evaluate_tool(
                label,
                str(output),
                str(tgt_fa),
                reference,
                manifest,
                evaluation_dir,
                original_version.get("profile"),
                log=lambda message: print(message, flush=True),
                ref_index=reference_index,
                threads=threads,
            ))
            transcript_path = (
                evaluation_dir / f"{label}.transcripts.tsv"
            ).resolve()
            if (
                Path(str(summary.get("transcripts_tsv", ""))).resolve()
                != transcript_path
                or not transcript_path.is_file()
                or transcript_path.stat().st_size <= 0
            ):
                raise RuntimeError(
                    f"{pair_path}: {label} rescore did not publish its TSV"
                )
            summary["stable_id_preservation"] = (
                release_evaluation.stable_id_preservation(
                    ref_gff,
                    output,
                    reference_index=reference_id_index,
                )
            )
            relative_tsv = str(transcript_path.relative_to(record_dir))
            summary["transcripts_tsv"] = relative_tsv
            versions[label] = {
                "summary": summary,
                "transcripts_tsv": {
                    "path": relative_tsv,
                    "size": transcript_path.stat().st_size,
                    "sha256": release_evaluation.sha256_file(transcript_path),
                },
            }
        record = {
            "schema_version": OVERLAY_SCHEMA_VERSION,
            "kind": RECORD_KIND,
            "key": {
                "panel": panel,
                "benchmark": benchmark,
                "repetition": repetition,
            },
            "source": source,
            "tooling_sha256": tooling_sha256,
            "versions": versions,
        }
        record_path = record_dir / "overlay.json"
        _write_json(record_path, record)
        manifest_rows.append({
            "key": dict(record["key"]),
            "source_pair": dict(source["pair_result"]),
            "record": str(record_path.relative_to(staging)),
            "record_sha256": release_evaluation.sha256_file(record_path),
        })
    return manifest_rows


def _write_overlay_once(
    roots: Sequence[Path],
    staging: Path,
    *,
    threads: int,
    allow_incomplete: bool,
) -> dict[str, Any]:
    pairs = release_report.load_pairs(roots, allow_incomplete=allow_incomplete)
    if not pairs:
        raise RuntimeError("no successful pair_result.json files found")
    seen_keys = {}
    verifier = _SourceVerifier()
    prepared = []
    for pair_path, pair in sorted(pairs, key=lambda item: _pair_key(item[1])):
        key = _pair_key(pair)
        if key in seen_keys:
            raise ValueError(
                f"duplicate paired repetition {key}: {seen_keys[key]} and {pair_path}"
            )
        seen_keys[key] = pair_path
        pair_record = _artifact_record(pair_path)
        inputs = pair.get("inputs")
        if not isinstance(inputs, Mapping) or not inputs:
            raise ValueError(f"{pair_path}: input evidence is missing")
        input_records = {
            name: verifier.verify(record, label=f"{pair_path}:{name}")
            for name, record in sorted(inputs.items())
        }
        outputs = {
            label: _output_record(pair_path, pair, label)
            for label in ("candidate", "reference")
        }
        prepared.append((pair_path, pair, {
            "pair_result": pair_record,
            "inputs": input_records,
            "outputs": outputs,
        }))

    tooling = _tooling_records()
    tooling_sha256 = _sha256_text(_canonical_json(tooling))
    grouped = defaultdict(list)
    for row in prepared:
        pair = row[1]
        grouped[(pair["panel"], pair["benchmark"], _group_signature(pair))].append(row)
    records = []
    for key in sorted(grouped):
        records.extend(_score_group(
            grouped[key],
            staging=staging,
            tooling_sha256=tooling_sha256,
            threads=threads,
        ))
    shutil.rmtree(staging / ".scratch", ignore_errors=True)
    manifest = {
        "schema_version": OVERLAY_SCHEMA_VERSION,
        "kind": OVERLAY_KIND,
        "created_at": _utc_now(),
        "run_roots": [str(Path(root).resolve()) for root in roots],
        "source_pair_count": len(prepared),
        "tooling": tooling,
        "tooling_sha256": tooling_sha256,
        "runtime": _runtime_record(),
        "git_state": build_controller.collect_git_state(
            build_controller.REPO_ROOT
        ),
        "records": sorted(
            records,
            key=lambda row: (
                row["key"]["panel"],
                row["key"]["benchmark"],
                row["key"]["repetition"],
            ),
        ),
    }
    if str(staging) in _canonical_json(manifest):
        raise ValueError("overlay manifest leaked its staging path")
    _write_json(staging / "manifest.json", manifest)
    return manifest


def _validate_destination(path: Path) -> None:
    if not (path.exists() or path.is_symlink()):
        return
    if path.is_symlink() or not path.is_dir():
        raise ValueError(f"refusing to replace non-directory overlay: {path}")
    manifest_path = path / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError(f"refusing to replace unrecognized overlay: {exc}") from exc
    if (
        not isinstance(manifest, Mapping)
        or manifest.get("schema_version") != OVERLAY_SCHEMA_VERSION
        or manifest.get("kind") != OVERLAY_KIND
    ):
        raise ValueError(f"refusing to replace unrecognized overlay: {path}")


def write_overlay(
    roots: Sequence[Path],
    output_dir: Path,
    *,
    threads: int = 8,
    allow_incomplete: bool = False,
) -> dict[str, Any]:
    if threads <= 0:
        raise ValueError("threads must be positive")
    roots = tuple(Path(root).resolve() for root in roots)
    output_dir = Path(output_dir).resolve()
    if any(
        output_dir == root
        or output_dir.is_relative_to(root)
        or root.is_relative_to(output_dir)
        for root in roots
    ):
        raise ValueError("overlay output may not overlap a sealed run root")
    _validate_destination(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(
        prefix=f".{output_dir.name}.tmp-", dir=output_dir.parent,
    ))
    quarantined = None
    try:
        manifest = _write_overlay_once(
            roots,
            staging,
            threads=threads,
            allow_incomplete=allow_incomplete,
        )
        if output_dir.exists():
            quarantined = Path(tempfile.mkdtemp(
                prefix=f".{output_dir.name}.stale-", dir=output_dir.parent,
            ))
            quarantined.rmdir()
            os.replace(output_dir, quarantined)
        os.replace(staging, output_dir)
        if quarantined is not None:
            shutil.rmtree(quarantined)
        return manifest
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        if quarantined is not None and quarantined.exists() and not output_dir.exists():
            os.replace(quarantined, output_dir)
        raise


def _relative_artifact(
    record_path: Path,
    raw: Any,
    *,
    label: str,
) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ValueError(f"{record_path}: {label} artifact is malformed")
    relative = raw.get("path")
    if not isinstance(relative, str) or Path(relative).is_absolute():
        raise ValueError(f"{record_path}: {label} path must be relative")
    path = (record_path.parent / relative).resolve()
    if not path.is_relative_to(record_path.parent.resolve()):
        raise ValueError(f"{record_path}: {label} escapes its overlay record")
    expected = {
        "path": relative,
        "size": raw.get("size"),
        "sha256": raw.get("sha256"),
    }
    observed = _artifact_record(path)
    comparable = {**observed, "path": relative}
    if expected != comparable:
        raise ValueError(f"{record_path}: {label} artifact was modified")
    return observed


def _load_overlay_root(root: Path) -> tuple[dict[Path, dict[str, Any]], dict[str, Any]]:
    root = Path(root).resolve()
    manifest_path = root / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError(f"{manifest_path}: overlay manifest is unreadable: {exc}") from exc
    required = {
        "schema_version", "kind", "created_at", "run_roots",
        "source_pair_count", "tooling", "tooling_sha256", "runtime",
        "git_state", "records",
    }
    if (
        not isinstance(manifest, Mapping)
        or set(manifest) != required
        or manifest.get("schema_version") != OVERLAY_SCHEMA_VERSION
        or manifest.get("kind") != OVERLAY_KIND
        or not isinstance(manifest.get("records"), list)
        or manifest.get("source_pair_count") != len(manifest["records"])
    ):
        raise ValueError(f"{manifest_path}: overlay manifest is malformed")
    tooling = manifest.get("tooling")
    if not isinstance(tooling, Mapping) or not tooling:
        raise ValueError(f"{manifest_path}: overlay tooling is missing")
    verified_tooling = {
        name: _verify_artifact_record(record, label=f"overlay tooling {name}")
        for name, record in sorted(tooling.items())
    }
    tooling_sha256 = _sha256_text(_canonical_json(verified_tooling))
    if tooling_sha256 != manifest.get("tooling_sha256"):
        raise ValueError(f"{manifest_path}: overlay tooling digest is inconsistent")
    index = {}
    seen_keys = set()
    for entry in manifest["records"]:
        if not isinstance(entry, Mapping):
            raise ValueError(f"{manifest_path}: overlay record entry is malformed")
        relative = entry.get("record")
        if not isinstance(relative, str) or Path(relative).is_absolute():
            raise ValueError(f"{manifest_path}: overlay record path must be relative")
        record_path = (root / relative).resolve()
        if not record_path.is_relative_to(root):
            raise ValueError(f"{manifest_path}: overlay record escapes its root")
        record_sha = release_evaluation.sha256_file(record_path)
        if record_sha != entry.get("record_sha256"):
            raise ValueError(f"{record_path}: overlay record was modified")
        record = json.loads(record_path.read_text())
        if (
            not isinstance(record, Mapping)
            or record.get("schema_version") != OVERLAY_SCHEMA_VERSION
            or record.get("kind") != RECORD_KIND
            or record.get("tooling_sha256") != tooling_sha256
        ):
            raise ValueError(f"{record_path}: overlay record is malformed")
        key_row = record.get("key")
        if not isinstance(key_row, Mapping):
            raise ValueError(f"{record_path}: overlay key is malformed")
        key = (
            key_row.get("panel"),
            key_row.get("benchmark"),
            key_row.get("repetition"),
        )
        if key in seen_keys or dict(key_row) != entry.get("key"):
            raise ValueError(f"{record_path}: duplicate or inconsistent overlay key")
        seen_keys.add(key)
        source = record.get("source")
        if not isinstance(source, Mapping):
            raise ValueError(f"{record_path}: overlay source is malformed")
        source_pair = _verify_artifact_record(
            source.get("pair_result"), label=f"{record_path} source pair",
        )
        if source_pair != entry.get("source_pair"):
            raise ValueError(f"{record_path}: source-pair evidence is inconsistent")
        versions = record.get("versions")
        if not isinstance(versions, Mapping) or set(versions) != {
            "candidate", "reference",
        }:
            raise ValueError(f"{record_path}: overlay versions are malformed")
        transcript_paths = {}
        for label in ("candidate", "reference"):
            version = versions[label]
            if not isinstance(version, Mapping) or not isinstance(
                version.get("summary"), Mapping
            ):
                raise ValueError(f"{record_path}: {label} overlay is malformed")
            transcript = _relative_artifact(
                record_path,
                version.get("transcripts_tsv"),
                label=f"{label} transcripts TSV",
            )
            if version["summary"].get("transcripts_tsv") != (
                version["transcripts_tsv"].get("path")
            ):
                raise ValueError(
                    f"{record_path}: {label} summary and TSV path disagree"
                )
            transcript_paths[label] = transcript
        pair_path = Path(source_pair["path"]).resolve()
        if pair_path in index:
            raise ValueError(f"duplicate overlay source pair: {pair_path}")
        index[pair_path] = {
            "root": root,
            "record_path": record_path,
            "record_sha256": record_sha,
            "record": record,
            "transcript_paths": transcript_paths,
        }
    evidence = {
        "root": str(root),
        "manifest": _artifact_record(manifest_path),
        "tooling_sha256": tooling_sha256,
        "evaluator_sha256": verified_tooling["evaluator"]["sha256"],
        "records": len(index),
        "created_at": manifest["created_at"],
        "runtime": manifest["runtime"],
        "git_state": manifest["git_state"],
    }
    return index, evidence


def apply_overlays(
    pairs: Sequence[tuple[Path, Mapping[str, Any]]],
    roots: Iterable[Path],
) -> tuple[list[tuple[Path, dict[str, Any]]], dict[str, Any]]:
    combined = {}
    evidence = []
    tooling_sha256 = None
    for root in roots:
        index, root_evidence = _load_overlay_root(Path(root))
        if tooling_sha256 is None:
            tooling_sha256 = root_evidence["tooling_sha256"]
        elif root_evidence["tooling_sha256"] != tooling_sha256:
            raise ValueError("evaluation overlays use mixed tooling snapshots")
        duplicated = set(combined) & set(index)
        if duplicated:
            raise ValueError(f"duplicate evaluation overlay pairs: {sorted(duplicated)}")
        combined.update(index)
        evidence.append(root_evidence)
    expected = {Path(path).resolve() for path, _ in pairs}
    observed = set(combined)
    if expected != observed:
        raise ValueError(
            "evaluation overlay coverage does not match selected pairs; "
            f"missing={sorted(expected - observed)}, extra={sorted(observed - expected)}"
        )
    transformed = []
    for pair_path, raw_pair in pairs:
        pair_path = Path(pair_path).resolve()
        overlay = combined[pair_path]
        record = overlay["record"]
        pair = copy.deepcopy(raw_pair)
        if _pair_key(pair) != (
            record["key"]["panel"],
            record["key"]["benchmark"],
            record["key"]["repetition"],
        ):
            raise ValueError(f"{pair_path}: evaluation overlay identity mismatch")
        if record["source"].get("inputs") != pair.get("inputs"):
            raise ValueError(f"{pair_path}: evaluation overlay input mismatch")
        for label in ("candidate", "reference"):
            output = record["source"]["outputs"][label]
            version = pair["versions"][label]
            output_path = Path(str(version.get("output_gff", ""))).resolve()
            if (
                output.get("path") != str(output_path)
                or output.get("size") != output_path.stat().st_size
                or output.get("fingerprints") != version.get("fingerprints")
            ):
                raise ValueError(f"{pair_path}: {label} overlay output mismatch")
            overlay_version = record["versions"][label]
            summary = copy.deepcopy(overlay_version["summary"])
            original_summary = version.get("summary")
            if isinstance(original_summary, Mapping) and "target_truth" in original_summary:
                summary["target_truth"] = copy.deepcopy(
                    original_summary["target_truth"]
                )
            transcript = overlay["transcript_paths"][label]
            summary["transcripts_tsv"] = transcript["path"]
            artifacts = copy.deepcopy(version.get("evaluation_artifacts") or {})
            artifacts["transcripts_tsv"] = dict(transcript)
            version["summary"] = summary
            version["evaluation_artifacts"] = artifacts
            version["evaluation_overlay"] = {
                "root": str(overlay["root"]),
                "record": str(overlay["record_path"]),
                "record_sha256": overlay["record_sha256"],
                "tooling_sha256": record["tooling_sha256"],
            }
            if pair["panel"] == "e2e":
                version["evaluation_summary"] = (
                    release_evaluation._neutral_e2e_summary(summary)
                )
                version["biological_summary"] = (
                    release_evaluation._neutral_e2e_biology(
                        summary, version.get("score_summary"),
                    )
                )
        transformed.append((pair_path, pair))
    return transformed, {
        "validated": True,
        "tooling_sha256": tooling_sha256,
        "pair_count": len(transformed),
        "roots": evidence,
    }


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-root", action="append", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--threads", type=_positive_int, default=8)
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help=(
            "select successful pair documents from diagnostic controllers "
            "while skipping failed, pending, or active cells"
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    write_overlay(
        [Path(root) for root in args.runs_root],
        Path(args.output_dir),
        threads=args.threads,
        allow_incomplete=args.allow_incomplete,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
