#!/usr/bin/env python
"""Iteration-4 feasibility diagnostic — measure the CHAINING merge's headroom.

Read-only. Runs LiftOn ONCE per benchmark (default best-of-outcome merge, with
`--write_chains`) on the cached `-L`/`-M`, then parses the per-transcript chain
log (`lifton_output/chain.txt`) to quantify where the protein-maximization
chaining (`lifton/protein_maximization.py:chaining_algorithm`) leaves identity
on the table. No code under `lifton/` is modified or monkeypatched — the chain
log already records each chunk's winner and AA window.

Each chain.txt line is `transcript_id\t<label>;<label>;...` where a label is
`miniprot[a-b]`, `liftoff[a-b]`, or `empty[a-b]` (a chunk both sources scored
0.0 on, whose CDS is currently DROPPED — protein_maximization.py:104-108). The
[a-b] is the chunk's protein-coordinate window, so `b-a` is its AA span (chunk
coarseness). chain.txt has one line per merge-FIRED transcript (the chains are
written inside the `has_valid_miniprot` branch), which is exactly the population
the chaining acts on.

Signals (the go/no-go for Iteration 4):
  - empty-drop frequency  -> headroom for refinement #2 (both-zero -> Liftoff
    fallback): if ~0, #2 is a no-op and we pivot or stop.
  - chunk coarseness (single-chunk fraction, max AA span) -> headroom for
    refinement #1 (relax the exact-genomic-end sync at :259): a high
    single-chunk fraction means sync points rarely fire and the merge is
    all-or-nothing over large spans.
  - merge-kept fraction (score.txt col 6) -> how often the merge candidate
    actually wins the best-of-outcome valve.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.chaining_headroom [IDS...]   # default = drosophila mouse_to_rat
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

DEFAULT_IDS = ["drosophila", "mouse_to_rat"]
_LABEL_RE = re.compile(r"(miniprot|liftoff|empty)\[([0-9.]+)-([0-9.]+)\]")


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run_default(bid: str, p: dict, root: Path, log=print):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    out = root / "default.gff3"
    argv = [TOOLS["lifton_bin"], "-t", "8", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), "--native", "--locus-pipeline",
            p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label=f"chain_headroom_{bid}",
                      log_dir=root / "logs", env=_compose_env(TOOLS), log=log)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: lifton failed (exit {pr.exit_code}); see {pr.stderr_path}")
    return root / "lifton_output" / "chain.txt", root / "lifton_output" / "score.txt"


def _parse_chains(chain_txt: Path):
    """transcript_id -> per-chunk stats parsed from chain.txt."""
    stats = {}
    if not chain_txt.exists():
        return stats
    for line in chain_txt.read_text().splitlines():
        parts = line.split("\t")
        if not parts or not parts[0]:
            continue
        tid = parts[0]
        labels = _LABEL_RE.findall(parts[1]) if len(parts) > 1 else []
        spans = [float(b) - float(a) for _, a, b in labels]
        stats[tid] = {
            "n_chunks": len(labels),
            "n_empty": sum(1 for k, _, _ in labels if k == "empty"),
            "n_miniprot": sum(1 for k, _, _ in labels if k == "miniprot"),
            "n_liftoff": sum(1 for k, _, _ in labels if k == "liftoff"),
            "max_span": round(max(spans), 1) if spans else 0.0,
            "total_span": round(sum(spans), 1) if spans else 0.0,
        }
    return stats


def _score_annotation_counts(score_txt: Path) -> dict:
    counts = {"LiftOn_chaining_algorithm": 0, "Liftoff": 0, "other": 0, "total": 0}
    if not score_txt.exists():
        return counts
    for line in score_txt.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) < 6:
            continue
        ann = cols[5]
        counts["total"] += 1
        counts[ann if ann in ("LiftOn_chaining_algorithm", "Liftoff") else "other"] += 1
    return counts


def _pct(n, d):
    return round(100.0 * n / d, 1) if d else 0.0


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: chaining headroom diagnostic ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_chain_headroom"
        root.mkdir(parents=True, exist_ok=True)
        _clean_input_dbs(p["ref_gff"],
                         W / "tools" / "liftoff" / "liftoff.gff3",
                         W / "tools" / "miniprot" / "miniprot.gff3")

        chain_txt, score_txt = _run_default(bid, p, root)
        chains = _parse_chains(chain_txt)
        fb = _score_annotation_counts(score_txt)

        n_fired = len(chains)
        n_with_empty = sum(1 for s in chains.values() if s["n_empty"] > 0)
        total_empty = sum(s["n_empty"] for s in chains.values())
        n_single_chunk = sum(1 for s in chains.values() if s["n_chunks"] == 1)
        n_multi = n_fired - n_single_chunk
        mean_chunks = round(sum(s["n_chunks"] for s in chains.values()) / n_fired, 2) if n_fired else 0.0
        # coarseness: among multi-CDS transcripts, how often is everything one chunk?
        spans = sorted((s["max_span"] for s in chains.values()), reverse=True)
        top_coarse = [(t, s["n_chunks"], s["max_span"])
                      for t, s in sorted(chains.items(), key=lambda kv: -kv[1]["max_span"])[:10]]

        rec = {
            "benchmark": bid,
            "n_merge_fired": n_fired,
            "merge_kept_valve": fb["LiftOn_chaining_algorithm"],
            "liftoff_kept_valve": fb["Liftoff"],
            # --- refinement #2 (both-zero -> Liftoff fallback) headroom ---
            "empty_drop": {
                "n_transcripts_with_empty": n_with_empty,
                "pct_of_fired": _pct(n_with_empty, n_fired),
                "total_empty_chunks": total_empty,
            },
            # --- refinement #1 (relax exact-end sync) headroom ---
            "coarseness": {
                "n_single_chunk": n_single_chunk,
                "pct_single_chunk": _pct(n_single_chunk, n_fired),
                "n_multi_chunk": n_multi,
                "mean_chunks_per_transcript": mean_chunks,
                "top_max_span_aa": spans[:5],
            },
            "top_coarse_transcripts": top_coarse,
        }
        results.append(rec)
        e, c = rec["empty_drop"], rec["coarseness"]
        print(f"  [{bid}] merge-fired={n_fired} (kept={rec['merge_kept_valve']} "
              f"liftoff={rec['liftoff_kept_valve']})", flush=True)
        print(f"    #2 empty-drop: {e['n_transcripts_with_empty']} transcripts "
              f"({e['pct_of_fired']}%), {e['total_empty_chunks']} dropped chunks", flush=True)
        print(f"    #1 coarseness: {c['n_single_chunk']} single-chunk "
              f"({c['pct_single_chunk']}%), mean {c['mean_chunks_per_transcript']} chunks/transcript",
              flush=True)

    (HERE / "chaining_headroom.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "chaining_headroom.md")
    print("\nDONE: wrote chaining_headroom.json + chaining_headroom.md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Chaining headroom diagnostic (Iteration 4, read-only)\n",
        "Default best-of-outcome run on cached `-L`/`-M`; per-transcript chain "
        "log parsed. **#2 empty-drop** = transcripts with a both-zero chunk whose "
        "CDS is currently dropped (headroom for the Liftoff-fallback refinement). "
        "**#1 coarseness** = single-chunk fraction (headroom for relaxing the "
        "exact-genomic-end sync). Decision: material #2 -> implement #2; mostly "
        "coarse -> consider #1 (higher risk); neither -> chaining near ceiling.\n",
        "| Dataset | merge-fired | kept/liftoff (valve) | #2 empty-drop txpts (%) | dropped chunks | #1 single-chunk (%) | mean chunks |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in results:
        e, c = r["empty_drop"], r["coarseness"]
        lines.append("| {} | {} | {}/{} | {} ({}%) | {} | {} ({}%) | {} |".format(
            r["benchmark"], r["n_merge_fired"], r["merge_kept_valve"],
            r["liftoff_kept_valve"], e["n_transcripts_with_empty"], e["pct_of_fired"],
            e["total_empty_chunks"], c["n_single_chunk"], c["pct_single_chunk"],
            c["mean_chunks_per_transcript"]))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
