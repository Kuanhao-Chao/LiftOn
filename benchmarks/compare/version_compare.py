#!/usr/bin/env python
"""Comprehensive A/B between the previous stable release (LiftOn v1.0.8) and the
current ``devel`` HEAD.

Runs the REAL v1.0.8 binary (built from the ``v1.0.8`` tag into the
``lifton_stable`` conda env) and the ``devel`` binary on the same inputs, then
scores BOTH outputs with the SAME version-agnostic neutral re-scorer
(``benchmarks.compare.evaluator``, which runs under ``lifton_devel`` and aligns
each predicted protein to the reference protein with LiftOn's own parasail
kernel). So accuracy is impartial — it never trusts either version's
self-reported ``protein_identity`` (those are parsed too, only as a cross-check).

Three arms (run independently via ``--arm``):

  controlled  cached -L/-M, ``-t 1`` for BOTH versions. The fair apples-to-apples
              path for per-transcript accuracy + the in-process LiftOn-logic
              wall/RSS. (devel ``-t 1`` is byte-identical to its threaded path
              per the 24-cell matrix, and v1.0.8 has no threaded Step 7.)
  fresh       no -L/-M — each version runs its OWN Liftoff+miniprot from scratch.
              The only fair source for COMPLETENESS (devel auto-detects gene-like
              top-level types; v1.0.8 lifts ``gene`` only). ``-t 1`` + HASHSEED=0
              for determinism (Liftoff ``-copies`` is otherwise non-deterministic).
  full        whole-genome (benchmarks/data/<id>), each version end-to-end. devel
              at its documented best config (``-t 8`` + stream/inmemory/locus),
              v1.0.8 at ``-t 1`` (its only option). Real-world scale.

A 4th synthetic version ``devel_legacy`` (devel + ``--legacy-merge --full-dp-align
--gene-only``) is run on the controlled arm only: it emulates v1.0.8's ALGORITHM
defaults on the devel binary, so the residual (devel_legacy vs real stable)
isolates the non-flag-gated differences (the Phase-4 gene-ID fix + the
Phase-13.5C canonical writer) from the flag-gated algorithm promotions.

Provenance is pinned on each install's git SHA (``lifton -V`` is unreliable — the
devel ``__version__`` still reads ``v1.0.8``).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.version_compare --arm controlled [IDS...]
    python -m benchmarks.compare.version_compare --arm fresh drosophila arabidopsis
    python -m benchmarks.compare.version_compare --arm full bee rice arabidopsis
Results accumulate into benchmarks/compare/version_compare.results.json.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from . import evaluator, id_mapping
from .profiling import run_profiled
from .tool_runners import _clean_input_dbs

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
DATA = HERE.parent / "data"           # benchmarks/data/<id> full-genome inputs
REG = json.loads((HERE / "benchmarks.json").read_text())
RESULTS = HERE / "version_compare.results.json"

DEVEL_ENV = "/home/kh.chao/miniconda3/envs/lifton_devel"
STABLE_ENV = "/home/kh.chao/miniconda3/envs/lifton_stable"
STABLE_WT = "/ccb/salz3/kh.chao/lifton_v1_0_8_worktree"
DEVEL_WT = "/ccb/salz3/kh.chao/LiftOn"

VERSIONS = {
    "stable": {                       # the real v1.0.8 release
        "bin": f"{STABLE_ENV}/bin/lifton",
        "validate": f"{STABLE_ENV}/bin/gff3-validate",
        "worktree": STABLE_WT,
        "expect_sha": "e503643d",
        "must_lack_native": True,
        "extra_flags": [],
    },
    "devel": {                        # current devel HEAD (Iteration 19)
        "bin": f"{DEVEL_ENV}/bin/lifton",
        "validate": f"{DEVEL_ENV}/bin/gff3-validate",
        "worktree": DEVEL_WT,
        "expect_sha": None,           # resolved at runtime (HEAD)
        "must_lack_native": False,
        "extra_flags": [],
    },
    "devel_legacy": {                 # devel emulating v1.0.8's algorithm defaults
        "bin": f"{DEVEL_ENV}/bin/lifton",
        "validate": f"{DEVEL_ENV}/bin/gff3-validate",
        "worktree": DEVEL_WT,
        "expect_sha": None,
        "must_lack_native": False,
        "extra_flags": ["--legacy-merge", "--full-dp-align", "--gene-only"],
    },
}

# Always-on devel validator as the consistent validity yardstick for both outputs.
DEVEL_VALIDATE = f"{DEVEL_ENV}/bin/gff3-validate"

# Full-genome inputs (salz3-cached, writable, .fai present).
FULL_INPUTS = {
    "bee": ("HAv3.1_genomic.fna", "HAv3.1_genomic.gff", "ASM1932182v1_genomic.fna"),
    "arabidopsis": ("TAIR10.fna", "TAIR10.gff", "ASM2311539v1_genomic.fna"),
    "rice": ("IRGSP_genomic.fna", "IRGSP_genomic.gff", "ASM3414082v1_genomic.fna"),
}

ENV_EXTRA_PATH = "/ccb/sw/bin:/home/kh.chao/bin"


def _env():
    import os
    return {
        "PATH": f"{ENV_EXTRA_PATH}:{os.environ.get('PATH', '')}",
        "PYTHONNOUSERSITE": "1",
        "PYTHONHASHSEED": "0",
    }


def _bench(bid):
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b
    raise KeyError(bid)


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------

def _git_head(worktree):
    return subprocess.run(["git", "-C", worktree, "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()


def provenance_gate(versions, log=print):
    """Pin each install on its git SHA + lifton.__file__ pointing into the
    worktree (lifton -V is unreliable). stable must additionally lack --native."""
    info = {}
    for v in versions:
        spec = VERSIONS[v]
        head = _git_head(spec["worktree"])
        exp = spec["expect_sha"] or head
        ok_sha = head.startswith(exp) or exp == head[:len(exp)]
        if spec["expect_sha"] and not head.startswith(spec["expect_sha"]):
            raise RuntimeError(f"{v}: worktree HEAD {head} != expect {spec['expect_sha']}")
        # confirm the binary's package resolves into the worktree (derive the
        # interpreter from the bin DIR, not a substring replace — the env path
        # itself contains "lifton").
        py = spec["bin"].rsplit("/", 1)[0] + "/python"
        # run from a NEUTRAL cwd (/tmp) — running from the devel repo root would
        # let cwd's ./lifton shadow the editable install and mis-report. The real
        # lifts run with cwd=statedir (no local lifton/), matching this.
        pkg = subprocess.run([py, "-c",
                              "import lifton,os;print(os.path.dirname(lifton.__file__))"],
                             capture_output=True, text=True, env=_env_os(),
                             cwd="/tmp").stdout.strip()
        if not pkg.startswith(spec["worktree"]):
            raise RuntimeError(f"{v}: lifton package resolves to {pkg}, "
                               f"not the expected worktree {spec['worktree']}")
        h = subprocess.run([spec["bin"], "-h"], capture_output=True, text=True)
        has_native = "--native" in ((h.stdout or "") + (h.stderr or ""))
        if spec["must_lack_native"] and has_native:
            raise RuntimeError(f"{v}: expected NO --native but help has it (wrong binary?)")
        info[v] = {"head": head, "pkg_dir": pkg, "has_native": has_native}
        log(f"  [prov] {v}: HEAD={head[:10]} pkg={pkg} has_native={has_native}")
    return info


def _env_os():
    import os
    e = dict(os.environ)
    e.update(_env())
    return e


# ---------------------------------------------------------------------------
# running a single lift
# ---------------------------------------------------------------------------

def _build_argv(version, paths, anndb, threads, L, M, out_gff, devel_fast):
    spec = VERSIONS[version]
    argv = [spec["bin"], "-t", str(threads), "-copies", "-ad", anndb,
            "-g", str(paths["ref_gff"])]
    if L is not None:
        argv += ["-L", str(L)]
    if M is not None:
        argv += ["-M", str(M)]
    if devel_fast and version.startswith("devel"):
        # BUG #2 (2026-06-18): --stream / --inmemory-liftoff UNDER-RECOVER (and
        # can segfault) on large distant genomes — the in-memory miniprot drain
        # drops Step-8 rescues at scale (full eudicot->monocot: 2730 vs the
        # default path's 7856 mRNA). This is a byte-identity-at-scale defect the
        # synthetic 24-cell fixture cannot see; tracked separately. Use ONLY
        # --locus-pipeline here: parallel Step 7 with default (disk) I/O is the
        # path the 24-cell matrix pins byte-identical at the stream=off /
        # inmemory=off / -t4 cell, so full-genome devel stays both fast and
        # correct. (Restore the in-memory flags once bug #2 is fixed.)
        argv += ["--locus-pipeline"]
    argv += list(spec["extra_flags"])
    if version.startswith("devel"):
        # Iteration 23: the miniprot-only rescue is now default-ON. Pin the FROZEN
        # fourway baselines explicitly OFF so the committed fourway_results.json
        # stays apples-to-apples (the rescue-ON recall gain is reported as a
        # separate arm via benchmarks/compare/miniprot_rescue_ab.py, NOT folded
        # into / re-baselined against the 4-way devel column). stable (v1.0.8)
        # has no such flag, so this is scoped to devel/devel_legacy.
        argv += ["--no-miniprot-rescue"]
    argv += ["-o", str(out_gff), str(paths["tgt_fa"]), str(paths["ref_fa"])]
    return argv


def run_lift(version, paths, anndb, threads, L, M, statedir, devel_fast=False,
             log=print):
    statedir.mkdir(parents=True, exist_ok=True)
    out_gff = statedir / f"{version}.gff3"
    argv = _build_argv(version, paths, anndb, threads, L, M, out_gff, devel_fast)
    pr = run_profiled(argv, label=version, log_dir=statedir / "logs",
                      env=_env(), cwd=statedir, log=log)
    if pr.exit_code != 0 or not out_gff.exists() or out_gff.stat().st_size == 0:
        raise RuntimeError(f"{version} lift failed (exit {pr.exit_code}); see {pr.stderr_path}")
    return out_gff, pr


# ---------------------------------------------------------------------------
# validity (consistent yardstick: devel gff3-validate on BOTH outputs)
# ---------------------------------------------------------------------------

def validate_gff(gff, log=print):
    """Run devel gff3-validate and parse the boxed summary's authoritative
    ``Errors : N`` / ``Warnings : N`` counts + the VALID/INVALID verdict."""
    import re
    r = subprocess.run([DEVEL_VALIDATE, str(gff)], capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    m_err = re.search(r"Errors\s*:\s*(\d+)", out)
    m_warn = re.search(r"Warnings\s*:\s*(\d+)", out)
    valid = ("VALID" in out) and ("INVALID" not in out)
    return {"exit": r.returncode,
            "n_errors": int(m_err.group(1)) if m_err else None,
            "n_warnings": int(m_warn.group(1)) if m_warn else None,
            "valid": valid}


# ---------------------------------------------------------------------------
# self-reported protein_identity (cross-check only)
# ---------------------------------------------------------------------------

def self_pi(gff, ref_ids):
    out = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "mRNA":
                continue
            attrs = {}
            for kv in c[8].split(";"):
                if "=" in kv:
                    k, val = kv.split("=", 1)
                    attrs[k] = val
            pi = attrs.get("protein_identity")
            if pi in (None, "", "None"):
                continue
            fid = attrs.get("ID", "")
            base, _ = id_mapping.strip_copy_suffix(fid, ref_ids)
            try:
                out[base] = float(pi)
            except ValueError:
                pass
    return out


# ---------------------------------------------------------------------------
# metric helpers over an evaluator transcripts.tsv
# ---------------------------------------------------------------------------

def _neutral_pi(tsv):
    """ref_mrna_id -> neutral protein_identity for coding, recovered rows."""
    import csv
    out = {}
    with open(tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("recovered") != "1" or row.get("is_coding") != "1":
                continue
            pi = row.get("protein_identity")
            if pi in (None, "", "None"):
                continue
            try:
                out[row["ref_mrna_id"]] = float(pi)
            except ValueError:
                pass
    return out


def _recovered_ids(tsv):
    import csv
    out = set()
    with open(tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("recovered") == "1":
                out.add(row["ref_mrna_id"])
    return out


def _mean(vals):
    vals = [v for v in vals if v is not None]
    return round(sum(vals) / len(vals), 5) if vals else None


def _delta(base, new):
    common = set(base) & set(new)
    improved = sum(1 for k in common if new[k] > base[k] + 1e-9)
    regressed = sum(1 for k in common if new[k] < base[k] - 1e-9)
    net = round(sum(new[k] - base[k] for k in common) / len(common), 6) if common else 0.0
    return {"n_common": len(common), "improved": improved,
            "regressed": regressed, "net_per_transcript": net}


def _self_vs_neutral(self_d, neutral_d):
    common = set(self_d) & set(neutral_d)
    if not common:
        return None
    return round(sum(self_d[k] - neutral_d[k] for k in common) / len(common), 6)


# ---------------------------------------------------------------------------
# arms
# ---------------------------------------------------------------------------

def _resolve(bid, arm):
    """Return (manifest, paths, anndb, L, M, devel_fast, versions)."""
    b = _bench(bid)
    anndb = b.get("annotation_database", "RefSeq")
    if arm in ("controlled", "fresh"):
        man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
        paths = man["paths"]
        if arm == "controlled":
            L = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
            M = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
            versions = ["stable", "devel", "devel_legacy"]
            return man, paths, anndb, L, M, False, versions
        return man, paths, anndb, None, None, False, ["stable", "devel"]
    # full genome
    ref_fa, ref_gff, tgt_fa = FULL_INPUTS[bid]
    d = DATA / bid
    paths = {"ref_fa": d / ref_fa, "ref_gff": d / ref_gff, "tgt_fa": d / tgt_fa}
    man = {"id": bid, "species": b["species"], "cross_species": b["cross_species"],
           "miniprot_target_space": b.get("miniprot_target_space", "transcript"),
           "protein_acc_to_mrna": {}, "paths": {k: str(v) for k, v in paths.items()}}
    return man, paths, anndb, None, None, True, ["stable", "devel"]


def run_arm(bid, arm, log=print):
    man, paths, anndb, L, M, devel_fast, versions = _resolve(bid, arm)
    provenance_gate(versions, log=log)
    ab_root = WORK / bid / "_version_ab" / arm
    ab_root.mkdir(parents=True, exist_ok=True)

    # clean cached -L/-M input DBs once (controlled arm) so the first run rebuilds.
    if arm == "controlled":
        _clean_input_dbs(paths["ref_gff"], L, M)

    gffs, profs, crashed = {}, {}, {}
    for v in versions:
        threads = 1
        if arm == "full" and v == "devel":
            threads = 8
        # Clean stale gffutils '_db' siblings before EACH version, on every arm:
        # the fresh/full arms share one (subset or full) ref_gff across both
        # versions, so without this the 2nd version collides with the 1st's ref_db
        # ("table features already exists"). For controlled this also drops -L/-M dbs.
        _clean_input_dbs(paths["ref_gff"], L, M)
        log(f"--- {bid} [{arm}] version {v} (-t {threads}) ---")
        try:
            g, pr = run_lift(v, paths, anndb, threads, L, M, ab_root / v,
                             devel_fast=devel_fast, log=log)
            gffs[v], profs[v] = g, pr
        except Exception as e:   # one version crashing must NOT lose the others
            sd = ab_root / v / "logs" / f"{v}.stderr.log"
            crashed[v] = str(sd)
            log(f"  !! {v} CRASHED ({e}); recording as crashed, continuing — see {sd}")
    versions = [v for v in versions if v in gffs]   # keep only what succeeded
    if not versions:
        raise RuntimeError(f"{bid}/{arm}: ALL versions failed: {list(crashed)}")

    # score every version's output with the SAME neutral evaluator
    ref, ref_index = evaluator.build_reference(str(paths["ref_gff"]),
                                               str(paths["ref_fa"]), log=log)
    ref_ids = set(ref.keys())
    eval_dir = ab_root / "eval"
    summaries = {}
    for v in versions:
        prof = {"wall_clock_seconds": profs[v].wall_clock_seconds,
                "peak_rss_mb": profs[v].peak_rss_mb}
        summaries[v] = evaluator.evaluate_tool(
            v, str(gffs[v]), str(paths["tgt_fa"]), ref, man, eval_dir, prof,
            log=log, ref_index=ref_index, threads=8)

    # assemble the comparison record
    neutral = {v: _neutral_pi(eval_dir / f"{v}.transcripts.tsv") for v in versions}
    selfpi = {v: self_pi(gffs[v], ref_ids) for v in versions}
    recovered = {v: _recovered_ids(eval_dir / f"{v}.transcripts.tsv") for v in versions}
    valid = {v: validate_gff(gffs[v], log=log) for v in versions}

    rec = {
        "benchmark": bid, "arm": arm,
        "n_reference_coding": summaries[versions[0]]["n_reference_coding"],
        "n_reference_total": summaries[versions[0]]["n_reference_total"],
        "versions": versions, "crashed": crashed,
        "neutral_mean_pi": {v: summaries[v]["protein_identity"]["mean"] for v in versions},
        "neutral_median_pi": {v: summaries[v]["protein_identity"]["median"] for v in versions},
        "neutral_pct_identical": {v: summaries[v]["protein_identity"]["pct_identical"] for v in versions},
        "self_mean_pi": {v: _mean(selfpi[v].values()) for v in versions},
        "self_vs_neutral_bias": {v: _self_vs_neutral(selfpi[v], neutral[v]) for v in versions},
        "completeness_coding": {v: summaries[v]["completeness_coding"] for v in versions},
        "completeness_feature_total": {v: summaries[v]["completeness_feature_total"] for v in versions},
        "n_recovered_coding": {v: summaries[v]["n_recovered_coding"] for v in versions},
        "n_recovered_any": {v: summaries[v]["n_recovered_any"] for v in versions},
        "wall_s": {v: round(profs[v].wall_clock_seconds, 2) for v in versions},
        "peak_rss_mb": {v: round(profs[v].peak_rss_mb, 1) for v in versions},
        "cpu_s": {v: round(profs[v].user_cpu_seconds + profs[v].sys_cpu_seconds, 1) for v in versions},
        "validity": valid,
        "feature_census": {v: summaries[v].get("completeness_by_type", {}) for v in versions},
    }
    # headline deltas vs stable
    if "devel" in versions and "stable" in versions:
        rec["delta_devel_vs_stable"] = _delta(neutral["stable"], neutral["devel"])
        s, d = rec["neutral_mean_pi"]["stable"], rec["neutral_mean_pi"]["devel"]
        rec["neutral_mean_pi_delta"] = (round(d - s, 6) if (s is not None and d is not None) else None)
        rec["recovered_gained_by_devel"] = len(recovered["devel"] - recovered["stable"])
        rec["recovered_lost_by_devel"] = len(recovered["stable"] - recovered["devel"])
    if "devel_legacy" in versions and "stable" in versions:
        # residual: devel_legacy (algorithm-emulated stable) vs real stable
        rec["delta_devel_legacy_vs_stable"] = _delta(neutral["stable"], neutral["devel_legacy"])
        rec["delta_devel_vs_devel_legacy"] = _delta(neutral["devel_legacy"], neutral["devel"])

    _save(rec)
    _print_rec(rec, log)
    return rec


def _save(rec):
    db = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}
    db[f"{rec['arm']}:{rec['benchmark']}"] = rec
    RESULTS.write_text(json.dumps(db, indent=2))


def _print_rec(rec, log):
    nm = rec["neutral_mean_pi"]
    log(f"  [{rec['benchmark']}/{rec['arm']}] neutral mean PI: " +
        " ".join(f"{v}={nm[v]}" for v in rec["versions"]))
    if "neutral_mean_pi_delta" in rec:
        d = rec["delta_devel_vs_stable"]
        log(f"      devel vs stable: Δmean={rec['neutral_mean_pi_delta']:+} "
            f"({d['improved']} improved / {d['regressed']} regressed); "
            f"recovered +{rec.get('recovered_gained_by_devel')}/-{rec.get('recovered_lost_by_devel')}; "
            f"wall {rec['wall_s']['stable']}->{rec['wall_s']['devel']}s")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=["controlled", "fresh", "full"])
    ap.add_argument("ids", nargs="*")
    a = ap.parse_args(argv)
    defaults = {
        "controlled": ["drosophila", "mouse_to_rat", "rice", "arabidopsis", "bee"],
        "fresh": ["drosophila", "arabidopsis"],
        "full": ["bee", "rice", "arabidopsis"],
    }
    ids = a.ids or defaults[a.arm]
    for bid in ids:
        print(f"\n=== {bid} [{a.arm}] : v1.0.8 stable vs devel ===", flush=True)
        try:
            run_arm(bid, a.arm, log=print)
        except Exception as e:
            import traceback
            print(f"!! {bid}/{a.arm} FAILED: {e}\n{traceback.format_exc()}", flush=True)
    print(f"\nDONE [{a.arm}] -> {RESULTS}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
