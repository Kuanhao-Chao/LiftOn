"""Content-addressed cache coverage for standalone benchmark tools."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.compare import tool_runners
from benchmarks.compare.profiling import ProfileResult


def _write_executable(path: Path, version: str) -> None:
    path.write_text(
        "#!/bin/sh\n"
        "if [ \"$1\" = \"--version\" ]; then\n"
        f"  echo '{version}'\n"
        "  exit 0\n"
        "fi\n"
        "exit 97\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


@pytest.fixture
def runner_case(tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    paths = {
        "ref_gff": inputs / "reference.gff3",
        "ref_fa": inputs / "reference.fa",
        "tgt_fa": inputs / "target.fa",
        "ref_faa": inputs / "reference.faa",
    }
    paths["ref_gff"].write_text(
        "##gff-version 3\nchr1\ttest\tgene\t1\t4\t.\t+\t.\tID=g1\n",
        encoding="utf-8",
    )
    paths["ref_fa"].write_text(">chr1\nACGT\n", encoding="utf-8")
    paths["tgt_fa"].write_text(">chr1\nACGT\n", encoding="utf-8")
    paths["ref_faa"].write_text(">p1\nM\n", encoding="utf-8")
    types_file = inputs / "types.txt"
    types_file.write_text("gene\n", encoding="utf-8")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    executables = {}
    for name in ("liftoff", "miniprot", "lifton"):
        executable = bin_dir / name
        _write_executable(executable, f"{name} 1.0")
        executables[name] = executable

    tools = {
        "liftoff_bin": str(executables["liftoff"]),
        "miniprot_bin": str(executables["miniprot"]),
        "lifton_bin": str(executables["lifton"]),
        "extra_path": str(bin_dir),
    }
    manifest = {"paths": {key: str(value) for key, value in paths.items()}}
    bench = {"id": "demo", "annotation_database": "RefSeq"}
    work = tmp_path / "work" / "demo"
    calls = []

    def fake_run_profiled(argv, *, label, log_dir, env=None, cwd=None, log=print):
        del env, cwd, log
        calls.append([str(value) for value in argv])
        log_dir.mkdir(parents=True, exist_ok=True)
        stdout = log_dir / f"{label}.stdout.log"
        stderr = log_dir / f"{label}.stderr.log"
        time_log = log_dir / f"{label}.time.log"
        stdout.write_text(
            "##gff-version 3\nchr1\ttest\tmRNA\t1\t4\t.\t+\t.\tID=t1\n",
            encoding="utf-8",
        )
        stderr.write_text("clean stderr\n", encoding="utf-8")
        time_log.write_text("profile clock\n", encoding="utf-8")
        if label == "liftoff":
            output = Path(argv[argv.index("-o") + 1])
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(
                "##gff-version 3\nchr1\ttest\tgene\t1\t4\t.\t+\t.\tID=g1\n",
                encoding="utf-8",
            )
        elif label.startswith("lifton"):
            output = Path(argv[argv.index("-o") + 1])
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(
                "##gff-version 3\nchr1\ttest\tgene\t1\t4\t.\t+\t.\tID=g1\n",
                encoding="utf-8",
            )
        return ProfileResult(
            wall_clock_seconds=1.25,
            peak_rss_mb=16.0,
            user_cpu_seconds=1.0,
            sys_cpu_seconds=0.1,
            exit_code=0,
            stdout_path=str(stdout),
            stderr_path=str(stderr),
            time_log_path=str(time_log),
        )

    monkeypatch.setattr(tool_runners, "run_profiled", fake_run_profiled)
    return {
        "bench": bench,
        "manifest": manifest,
        "work": work,
        "tools": tools,
        "paths": paths,
        "types_file": types_file,
        "executables": executables,
        "calls": calls,
    }


def _run_miniprot(case, *, threads=4):
    return tool_runners.run_miniprot(
        case["bench"], case["manifest"], case["work"], case["tools"],
        threads=threads, force=False, log=lambda _message: None,
    )


def test_miniprot_manifest_binds_complete_request_and_reuses(runner_case):
    case = runner_case
    first = _run_miniprot(case)
    manifest_path = (
        case["work"] / "tools" / "miniprot" /
        tool_runners.TOOL_CACHE_MANIFEST_NAME
    )
    document = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert document["schema_version"] == tool_runners.TOOL_CACHE_SCHEMA_VERSION
    assert document["builder_version"] == tool_runners.TOOL_CACHE_BUILDER_VERSION
    assert document["request_sha256"] == tool_runners._canonical_sha256(
        document["request"],
    )
    request = document["request"]
    assert request["argv"] == [
        case["tools"]["miniprot_bin"], "-t", "4", "--gff-only",
        case["manifest"]["paths"]["tgt_fa"],
        case["manifest"]["paths"]["ref_faa"],
    ]
    assert request["configuration"] == {
        "threads": 4,
        "mode": "genes",
        "options": {"gff_only": True},
    }
    assert set(request["inputs"]) == {"reference_proteins", "target_fasta"}
    assert set(request["executable"]) == {
        "configured", "path", "size", "sha256", "version",
        "version_argv", "version_exit_code",
    }
    assert request["executable"]["version"] == "miniprot 1.0"
    assert request["builder"]["source_sha256"] == tool_runners._sha256_file(
        Path(tool_runners.__file__),
    )
    assert "PATH" in request["environment"]
    assert request["environment"]["PYTHONNOUSERSITE"] == "1"
    assert set(document["artifacts"]) == {"output", "profile"}
    assert tool_runners.verify_cached_tool_run(case["work"], "miniprot")

    second = _run_miniprot(case)
    assert second == first
    assert len(case["calls"]) == 1


def test_done_and_outputs_without_manifest_never_permit_reuse(runner_case):
    case = runner_case
    _run_miniprot(case)
    manifest_path = (
        case["work"] / "tools" / "miniprot" /
        tool_runners.TOOL_CACHE_MANIFEST_NAME
    )
    manifest_path.unlink()

    assert (case["work"] / ".done" / "miniprot.done").is_file()
    assert not tool_runners.verify_cached_tool_run(case["work"], "miniprot")
    _run_miniprot(case)
    assert len(case["calls"]) == 2
    assert tool_runners.verify_cached_tool_run(case["work"], "miniprot")


@pytest.mark.parametrize("tamper", ["output", "profile", "input", "executable"])
def test_tampered_or_stale_miniprot_cache_is_rebuilt(runner_case, tamper):
    case = runner_case
    _run_miniprot(case)
    if tamper == "output":
        path = case["work"] / "tools" / "miniprot" / "miniprot.gff3"
        path.write_text(path.read_text(encoding="utf-8") + "# tampered\n", encoding="utf-8")
    elif tamper == "profile":
        path = case["work"] / "tools" / "logs" / "miniprot.profile.json"
        path.write_text("{}\n", encoding="utf-8")
    elif tamper == "input":
        case["paths"]["tgt_fa"].write_text(">chr1\nTGCA\n", encoding="utf-8")
    else:
        _write_executable(case["executables"]["miniprot"], "miniprot 2.0")

    assert not tool_runners.verify_cached_tool_run(case["work"], "miniprot")
    _run_miniprot(case)
    assert len(case["calls"]) == 2
    assert tool_runners.verify_cached_tool_run(case["work"], "miniprot")


def test_threads_change_invalidates_exact_request(runner_case):
    case = runner_case
    _run_miniprot(case, threads=4)
    _run_miniprot(case, threads=2)

    assert len(case["calls"]) == 2
    manifest_path = (
        case["work"] / "tools" / "miniprot" /
        tool_runners.TOOL_CACHE_MANIFEST_NAME
    )
    request = json.loads(manifest_path.read_text(encoding="utf-8"))["request"]
    assert request["configuration"]["threads"] == 2
    assert request["argv"][2] == "2"


def test_liftoff_and_lifton_variants_require_verified_dependency_caches(runner_case):
    case = runner_case
    run_liftoff = lambda **kwargs: tool_runners.run_liftoff(
        case["bench"], case["manifest"], case["work"], case["tools"],
        threads=4, force=False, log=lambda _message: None, **kwargs,
    )
    run_lifton = lambda **kwargs: tool_runners.run_lifton(
        case["bench"], case["manifest"], case["work"], case["tools"],
        threads=4, force=False, log=lambda _message: None, **kwargs,
    )

    run_liftoff()
    _run_miniprot(case)
    run_lifton()
    run_lifton(optimize=True)
    run_liftoff(mode="allfeat", types_file=case["types_file"])
    run_lifton(mode="allfeat", types_file=case["types_file"])

    assert tool_runners.verify_cached_tool_run(case["work"], "liftoff")
    assert tool_runners.verify_cached_tool_run(case["work"], "miniprot")
    assert tool_runners.verify_cached_tool_run(case["work"], "lifton")
    assert tool_runners.verify_cached_tool_run(case["work"], "lifton_optimize")
    assert tool_runners.verify_cached_tool_run(case["work"], "liftoff", "allfeat")
    assert tool_runners.verify_cached_tool_run(case["work"], "lifton", "allfeat")

    liftoff_output = case["work"] / "tools" / "liftoff" / "liftoff.gff3"
    liftoff_output.write_text("tampered dependency\n", encoding="utf-8")
    assert not tool_runners.verify_cached_tool_run(case["work"], "liftoff")
    with pytest.raises(RuntimeError, match="verified Liftoff cache"):
        tool_runners.run_lifton(
            case["bench"], case["manifest"], case["work"], case["tools"],
            threads=2, force=False, log=lambda _message: None,
        )
