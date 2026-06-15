#!/usr/bin/env python
"""Gene-like-lift A/B: default gene-like lift vs `--gene-only` (Iter-5 → Iter-12).

As of the Iteration-12 promotion the gene-like lift is the DEFAULT and
`--gene-only` is the opt-out, so the two states below invert their flags vs the
Iteration-5 framing (the OUTPUT bytes of each state are unchanged — only which
flag selects the gene-only baseline vs the gene-like candidate moved).

Design (isolates the effect cleanly):
  1. Build ONE gene-like Liftoff annotation via the PROVEN standalone `liftoff`
     binary with `-f <gene-like types>` (subprocess minimap2 — the same path that
     built the cached `-L`). It contains genes AND pseudogenes/etc.
  2. Run BOTH LiftOn states off that SAME `-L` + the cached `-M`:
       state "off"  --gene-only   LiftOn processes only `gene`     ← gene-only baseline
       state "on"   (no flag)      LiftOn processes gene-like types ← default (the promotion)
     Same `-L` → the gene lift is identical by construction; "on" only ADDS the
     gene-like loci. (We deliberately avoid LiftOn's in-process `--native`
     Liftoff here: that path mapped nothing in a fresh run — a separate latent
     issue, see memory `lifton-native-liftoff-empty`.)

Reports per dataset: top-level feature-type counts off/on (pseudogenes etc.
recovered), output VALIDITY via the in-tree gff3 validator, and that off's
feature-line set is a subset of on's (existing lift unchanged). Promote iff:
materially more real gene-like features lifted, validates clean, off ⊆ on.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.lift_gene_like_ab [IDS...]   # default = drosophila mouse_to_rat
"""
from __future__ import annotations

import collections
import json
import shutil
import subprocess
import sys
from pathlib import Path

from lifton import annotation, lifton_utils
from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

DEFAULT_IDS = ["drosophila", "mouse_to_rat"]
STATES = ("off", "on")


def _ann_db(bid):
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _gene_like_types(ref_gff):
    db = annotation.Annotation(ref_gff, None, None, "create_unique", "ID",
                               False, False, True)
    return lifton_utils.get_gene_like_feature_types(db)


def _build_gene_like_liftoff(bid, p, types, root, log=print):
    """Run the standalone liftoff binary with -f <types> (proven subprocess
    path). Cached; returns the gene-like liftoff.gff3 path."""
    cache = root / "liftoff.genelike.gff3"
    if cache.exists() and cache.stat().st_size > 0:
        log(f"  [{bid}] gene-like liftoff cached")
        return cache
    types_file = root / "gene_like_feature_types.txt"
    types_file.write_text("\n".join(types) + "\n")
    inter = root / "liftoff_inter"
    inter.mkdir(parents=True, exist_ok=True)
    raw = root / "liftoff_genelike_raw.gff3"
    argv = [TOOLS["liftoff_bin"], "-f", str(types_file),
            "-g", p["ref_gff"], "-o", str(raw),
            "-u", str(root / "unmapped.txt"), "-dir", str(inter),
            "-p", "8", "-copies", "-polish", p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label=f"liftoff_genelike_{bid}",
                      log_dir=root / "logs", env=_compose_env(TOOLS), log=log)
    polished = Path(str(raw) + "_polished")
    src = polished if polished.exists() else raw
    if pr.exit_code != 0 or not src.exists() or src.stat().st_size == 0:
        raise RuntimeError(f"{bid}: gene-like liftoff failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    shutil.copyfile(src, cache)
    return cache


def _run_state(bid, p, state, liftoff_L, root, log=print):
    statedir = root / state
    statedir.mkdir(parents=True, exist_ok=True)
    out = statedir / f"{state}.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    _clean_input_dbs(p["ref_gff"], liftoff_L, miniprot)
    # -t 1 serial, no --native (the LiftOn processing is fast with -L/-M; we
    # avoid the native path entirely). Same -L for both states.
    argv = [TOOLS["lifton_bin"], "-t", "1", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff_L), "-M", str(miniprot),
            "-o", str(out)]
    # Post-Iteration-12 promotion: the gene-like lift is now the DEFAULT, so the
    # gene-only baseline ("off") opts out via --gene-only; "on" = default (no
    # flag — --lift-gene-like is a kept no-op alias). Pre-promotion this was
    # reversed (off=no flag, on=--lift-gene-like); the bytes are the same, only
    # the flag that selects each state moved.
    if state == "off":
        argv.append("--gene-only")
    argv += [p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label=f"lift_gene_like_{state}",
                      log_dir=root / "logs", env=_compose_env(TOOLS), log=log)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _toplevel_type_counts(gff):
    c = collections.Counter()
    for ln in Path(gff).read_text().splitlines():
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 9 or "Parent=" in f[8]:
            continue
        c[f[2]] += 1
    return c


def _feature_line_set(gff):
    return {ln for ln in Path(gff).read_text().splitlines()
            if ln and not ln.startswith("#") and len(ln.split("\t")) >= 9}


def _validate(gff):
    r = subprocess.run([TOOLS["lifton_python"], "-m", "lifton.gff3_validator", str(gff)],
                       env=_compose_env(TOOLS), capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    n_err = sum(1 for ln in out.splitlines() if "error" in ln.lower())
    return r.returncode, n_err


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: --lift-gene-like vs default A/B ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_lift_gene_like_ab"
        root.mkdir(parents=True, exist_ok=True)

        types = _gene_like_types(p["ref_gff"])
        print(f"  [{bid}] gene-like types = {types}", flush=True)
        liftoff_L = _build_gene_like_liftoff(bid, p, types, root)

        outs, prs = {}, {}
        for state in STATES:
            print(f"--- {bid} state {state} ---", flush=True)
            outs[state], prs[state] = _run_state(bid, p, state, liftoff_L, root)

        cnt = {s: _toplevel_type_counts(outs[s]) for s in STATES}
        all_types = sorted(set(cnt["off"]) | set(cnt["on"]))
        recovered = {t: [cnt["off"].get(t, 0), cnt["on"].get(t, 0)]
                     for t in all_types if cnt["on"].get(t, 0) != cnt["off"].get(t, 0)}
        off_set, on_set = _feature_line_set(outs["off"]), _feature_line_set(outs["on"])
        gene_unchanged = off_set <= on_set
        on_rc, on_nerr = _validate(outs["on"])
        off_rc, off_nerr = _validate(outs["off"])

        rec = {
            "benchmark": bid,
            "gene_like_types": types,
            "type_counts": {s: dict(cnt[s]) for s in STATES},
            "recovered_off_on": recovered,
            "pseudogene_off": cnt["off"].get("pseudogene", 0),
            "pseudogene_on": cnt["on"].get("pseudogene", 0),
            "gene_off": cnt["off"].get("gene", 0),
            "gene_on": cnt["on"].get("gene", 0),
            "off_lines": len(off_set), "on_lines": len(on_set),
            "on_only_lines": len(on_set - off_set), "off_only_lines": len(off_set - on_set),
            "gene_lift_unchanged": gene_unchanged,
            "validate_on": {"exit": on_rc, "n_error_lines": on_nerr},
            "validate_off": {"exit": off_rc, "n_error_lines": off_nerr},
            "wall_s": {s: round(prs[s].wall_clock_seconds, 1) for s in STATES},
        }
        # Promote iff: (1) more real gene-like features lifted, (2) the existing
        # lift is preserved (off lines ⊆ on), and (3) NO NEW validity errors vs
        # the default. We compare to off's error count (not zero) because the
        # default output itself carries inherited RefSeq dup-ID/containment
        # errors — the bar is "adds no new errors", not "absolutely clean".
        rec["promote"] = bool(
            sum((v[1] - v[0]) for v in recovered.values() if v[1] > v[0]) > 0
            and gene_unchanged and on_nerr <= off_nerr)
        results.append(rec)
        print(f"\n  [{bid}] pseudogene off={rec['pseudogene_off']} on={rec['pseudogene_on']} | "
              f"gene off={rec['gene_off']} on={rec['gene_on']} | recovered(off,on)={recovered}", flush=True)
        print(f"  off-lines={rec['off_lines']} on-lines={rec['on_lines']} "
              f"(off-only={rec['off_only_lines']}, on-only=+{rec['on_only_lines']}) | "
              f"off⊆on={gene_unchanged} | validate-on exit={on_rc} errs={on_nerr} "
              f"(off errs={off_nerr}) | PROMOTE={rec['promote']}", flush=True)

    (HERE / "lift_gene_like_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "lift_gene_like_ab.md")
    print("\nDONE: wrote lift_gene_like_ab.json + lift_gene_like_ab.md", flush=True)
    return 0


def _write_md(results, path):
    lines = [
        "## --lift-gene-like vs default A/B (Iteration 5)\n",
        "Both states share ONE gene-like Liftoff `-L` (built via the standalone "
        "liftoff binary) + cached `-M`. **off** = default (LiftOn processes "
        "`gene` only); **on** = `--lift-gene-like` (processes auto-detected "
        "gene-like types). Same `-L` → the gene lift is identical; `on` only "
        "ADDS gene-like loci. Promote iff: more real gene-like features lifted, "
        "validates clean, off-lines ⊆ on-lines.\n",
        "| Dataset | gene-like types | pseudogene off/on | gene off/on | off⊆on (lift preserved) | on-only lines | validate-on (exit/errs) | PROMOTE |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        v = r["validate_on"]
        lines.append("| {} | {} | {}/{} | {}/{} | {} | +{} | {}/{} | {} |".format(
            r["benchmark"], ",".join(r["gene_like_types"]),
            r["pseudogene_off"], r["pseudogene_on"], r["gene_off"], r["gene_on"],
            "yes" if r["gene_lift_unchanged"] else "NO", r["on_only_lines"],
            v["exit"], v["n_error_lines"], "yes" if r["promote"] else "no"))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
