"""Profiled standalone runners for Liftoff, miniprot, and LiftOn.

Each runner writes its output under ``work/<id>/tools/<tool>/`` and a profile
JSON (wall-clock + peak RSS) under ``tools/logs/``; each is guarded by a
``.done/<tool>.done`` sentinel for resumability. LiftOn reuses the standalone
Liftoff + miniprot outputs via ``-L`` / ``-M``.
"""
from __future__ import annotations

import dataclasses
import json
import os
import shutil
from pathlib import Path

from .profiling import run_profiled


def _compose_env(tools: dict) -> dict:
    """PATH that lets every subprocess find minimap2/miniprot/samtools, plus
    PYTHONNOUSERSITE to avoid user-site bleed across conda envs."""
    extra = tools.get("extra_path", "")
    env = {
        "PATH": f"{extra}:{os.environ.get('PATH', '')}" if extra else os.environ.get("PATH", ""),
        "PYTHONNOUSERSITE": "1",
    }
    return env


def _save_profile(pr, path: Path) -> dict:
    d = dataclasses.asdict(pr)
    path.write_text(json.dumps(d, indent=2))
    return d


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
    if done.exists() and not force and out_gff.exists():
        log(f"  [liftoff:{mode}] cached")
        return json.loads(prof_path.read_text()) if prof_path.exists() else {}
    out_dir.mkdir(parents=True, exist_ok=True)
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
    pr = run_profiled(argv, label="liftoff", log_dir=log_dir,
                      env=_compose_env(tools), log=log)
    if pr.exit_code != 0:
        raise RuntimeError(f"liftoff failed (exit {pr.exit_code}); see {pr.stderr_path}")
    # -polish makes Liftoff emit <out>_polished; normalize to liftoff.gff3.
    polished = Path(str(raw_out) + "_polished")
    src = polished if polished.exists() else raw_out
    if not src.exists():
        raise RuntimeError(f"liftoff produced no output ({raw_out} / {polished})")
    shutil.copyfile(src, out_gff)
    prof = _save_profile(pr, prof_path)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("ok\n")
    log(f"  [liftoff:{mode}] {pr.wall_clock_seconds:.1f}s, {pr.peak_rss_mb:.0f} MB")
    return prof


def run_miniprot(bench: dict, manifest: dict, work_dir: Path, tools: dict,
                 threads: int, force: bool, log=print) -> dict:
    out_dir = work_dir / "tools" / "miniprot"
    log_dir = work_dir / "tools" / "logs"
    done = work_dir / ".done" / "miniprot.done"
    out_gff = out_dir / "miniprot.gff3"
    prof_path = log_dir / "miniprot.profile.json"
    if done.exists() and not force and out_gff.exists():
        log("  [miniprot] cached")
        return json.loads(prof_path.read_text()) if prof_path.exists() else {}
    out_dir.mkdir(parents=True, exist_ok=True)
    p = manifest["paths"]
    argv = [
        tools["miniprot_bin"], "-t", str(threads), "--gff-only",
        p["tgt_fa"], p["ref_faa"],
    ]
    # run_profiled redirects stdout to <label>.stdout.log; miniprot --gff-only
    # writes GFF to stdout, so the GFF is that stdout file.
    pr = run_profiled(argv, label="miniprot", log_dir=log_dir,
                      env=_compose_env(tools), log=log)
    if pr.exit_code != 0:
        raise RuntimeError(f"miniprot failed (exit {pr.exit_code}); see {pr.stderr_path}")
    err = Path(pr.stderr_path).read_text(errors="replace") if pr.stderr_path else ""
    if "ERROR" in err.upper():
        raise RuntimeError(f"miniprot reported ERROR on stderr; see {pr.stderr_path}")
    shutil.copyfile(pr.stdout_path, out_gff)
    if out_gff.stat().st_size == 0:
        raise RuntimeError("miniprot produced empty output")
    prof = _save_profile(pr, prof_path)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("ok\n")
    log(f"  [miniprot] {pr.wall_clock_seconds:.1f}s, {pr.peak_rss_mb:.0f} MB")
    return prof


def run_lifton(bench: dict, manifest: dict, work_dir: Path, tools: dict,
               threads: int, force: bool, log=print, mode="genes",
               types_file=None) -> dict:
    tools_root, suf = _mode_dirs(work_dir, mode)
    out_dir = tools_root / "lifton"
    log_dir = tools_root / "logs"
    done = work_dir / ".done" / f"lifton{suf}.done"
    out_gff = out_dir / "lifton.gff3"
    prof_path = log_dir / "lifton.profile.json"
    if done.exists() and not force and out_gff.exists():
        log(f"  [lifton:{mode}] cached")
        return json.loads(prof_path.read_text()) if prof_path.exists() else {}
    out_dir.mkdir(parents=True, exist_ok=True)
    p = manifest["paths"]
    liftoff_gff = tools_root / "liftoff" / "liftoff.gff3"          # mode-specific
    miniprot_gff = work_dir / "tools" / "miniprot" / "miniprot.gff3"  # shared, mode-independent
    if not liftoff_gff.exists():
        raise RuntimeError(f"lifton:{mode} needs {liftoff_gff} (run liftoff mode={mode} first)")
    if not miniprot_gff.exists():
        raise RuntimeError(f"lifton:{mode} needs the shared miniprot {miniprot_gff} "
                           f"(run genes mode first)")
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
    if mode != "genes":
        if types_file is None:
            raise RuntimeError(f"lifton mode {mode!r} requires a feature-types file")
        argv += ["-f", str(types_file)]   # process all listed top-level types
    argv += [p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label="lifton", log_dir=log_dir,
                      env=_compose_env(tools), log=log)
    if pr.exit_code != 0:
        raise RuntimeError(f"lifton failed (exit {pr.exit_code}); see {pr.stderr_path}")
    if not out_gff.exists() or out_gff.stat().st_size == 0:
        raise RuntimeError(f"lifton produced no/empty output {out_gff}")
    prof = _save_profile(pr, prof_path)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("ok\n")
    log(f"  [lifton:{mode}] {pr.wall_clock_seconds:.1f}s, {pr.peak_rss_mb:.0f} MB")
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
        if not (work_dir / "tools" / "miniprot" / "miniprot.gff3").exists():
            raise RuntimeError(f"{bench['id']}: mode {mode!r} needs the shared genes-mode "
                               f"miniprot — run genes mode first")
        mp_prof = work_dir / "tools" / "logs" / "miniprot.profile.json"
        profiles["miniprot"] = json.loads(mp_prof.read_text()) if mp_prof.exists() else {}
    profiles["lifton"] = run_lifton(bench, manifest, work_dir, tools, threads,
                                    force, log, mode=mode, types_file=types_file)
    return profiles
