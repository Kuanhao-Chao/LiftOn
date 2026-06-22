#!/usr/bin/env python
"""A/B for the giant-gene windowed-alignment fix (Iteration 2).

Runs LiftOn on a benchmark TWICE on the same cached `-L`/`-M`, differing only by
the alignment kernel:
  state "windowed"  default                          anchor-windowed above threshold
  state "fulldp"    LIFTON_ALIGN_WINDOW_{AA,NT}=1e9  forced full-DP (the baseline)

Both are single-threaded (`-t 1`) so at most one giant alignment is in flight
(bounded peak RSS) and the comparison isolates the windowing effect. It reports
wall-clock + peak RSS for each state, asserts that every transcript whose output
differs is a "giant" (ref protein > AA threshold or ref transcript > NT
threshold) — i.e. normal genes are byte-identical — and reports the per-giant
protein_identity delta (windowed vs full-DP; should be ~0).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.align_window_ab [IDS...]   # default = mouse
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

DEFAULT_IDS = ["mouse"]
AA_THRESH = 8000        # must match lifton/align.py _ALIGN_WINDOW_AA default
NT_THRESH = 25000       # must match lifton/align.py _ALIGN_WINDOW_NT default

# NOTE: as of Iteration 3, "band everything" is the DEFAULT, so the giant-only
# windowing this Iter-2 A/B measures is now reached via LIFTON_FULL_DP_ALIGN=1
# (the --full-dp-align escape hatch). "fulldp" forces pure full DP incl. giants.
STATES = {
    "windowed": {"LIFTON_FULL_DP_ALIGN": "1"},
    "fulldp": {"LIFTON_ALIGN_WINDOW_AA": "1000000000",
               "LIFTON_ALIGN_WINDOW_NT": "1000000000"},
}


def _ann_db(bid):
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _giant_mrna_ids(bid, p):
    """mRNA IDs whose ref protein > AA_THRESH aa OR ref transcript > NT_THRESH nt
    — the transcripts that take the windowed path (output may legitimately
    differ from full-DP)."""
    man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
    acc2mrna = man.get("protein_acc_to_mrna", {})
    giant = set()
    # protein-length giants
    prot = WORK / bid / "subset" / "ref.proteins.subset.faa"
    if prot.exists():
        acc = None
        ln = 0
        def flush():
            if acc and ln > AA_THRESH and acc in acc2mrna:
                giant.add(acc2mrna[acc])
        for line in prot.read_text().splitlines():
            if line.startswith(">"):
                flush(); acc = line[1:].split()[0]; ln = 0
            else:
                ln += len(line.strip())
        flush()
    # transcript-length giants (sum of exon lengths per mRNA in the ref GFF)
    exon_len = {}
    for line in Path(p["ref_gff"]).read_text().splitlines():
        if line.startswith("#"):
            continue
        c = line.split("\t")
        if len(c) < 9 or c[2] != "exon":
            continue
        attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
        par = attrs.get("Parent", "")
        exon_len[par] = exon_len.get(par, 0) + (int(c[4]) - int(c[3]) + 1)
    giant |= {mid for mid, L in exon_len.items() if L > NT_THRESH}
    return giant


_TRANSCRIPT_TYPES = {"mRNA", "lnc_RNA", "transcript", "ncRNA", "snRNA", "snoRNA",
                     "miRNA", "rRNA", "tRNA", "primary_transcript", "antisense_RNA",
                     "scRNA", "guide_RNA", "RNase_P_RNA", "telomerase_RNA",
                     "SRP_RNA", "RNase_MRP_RNA", "Y_RNA", "vault_RNA"}


def _mrna_blocks(gff_path):
    """Map transcript-id -> list of its output lines, and transcript-id ->
    protein_identity. A line is keyed to the TRANSCRIPT it belongs to: a
    transcript-level feature (mRNA/lnc_RNA/... whose Parent is a gene) by its
    own ID; a child (exon/CDS) by its Parent; a gene line by its own ID. This
    keeps coding AND non-coding transcripts correctly attributed (a long
    lnc_RNA that takes the windowed DNA path is keyed by its own rna-* id, not
    its parent gene)."""
    lines = [ln for ln in Path(gff_path).read_text().splitlines()
             if ln and not ln.startswith("#") and len(ln.split("\t")) >= 9]
    gene_ids = set()
    for ln in lines:
        c = ln.split("\t")
        if c[2] == "gene" or c[2].endswith("_gene") or c[2] == "pseudogene":
            a = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            if "ID" in a:
                gene_ids.add(a["ID"])
    blocks, pid = {}, {}
    for ln in lines:
        c = ln.split("\t")
        attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
        par = attrs.get("Parent", "")
        if par and par in gene_ids:          # transcript-level
            tid = attrs.get("ID")
        elif par:                            # child (exon/CDS) of a transcript
            tid = par
        else:                                # gene / top-level
            tid = attrs.get("ID")
        if not tid:
            continue
        blocks.setdefault(tid, []).append(ln)
        if c[2] in _TRANSCRIPT_TYPES and "protein_identity" in attrs:
            pid[attrs["ID"]] = float(attrs["protein_identity"])
    return blocks, pid


def _run(bid, state, env_extra, p, root):
    out_dir = root / state
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{state}.gff3"
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    argv = [TOOLS["lifton_bin"], "-t", "1", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), p["tgt_fa"], p["ref_fa"]]
    env = _compose_env(TOOLS)
    env.update(env_extra)
    pr = run_profiled(argv, label=f"align_window_{state}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: windowed vs full-DP alignment A/B ===", flush=True)
        man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = WORK / bid / "_align_window_ab"
        root.mkdir(parents=True, exist_ok=True)
        _clean_input_dbs(p["ref_gff"],
                         WORK / bid / "tools" / "liftoff" / "liftoff.gff3",
                         WORK / bid / "tools" / "miniprot" / "miniprot.gff3")

        outs, prs = {}, {}
        for state, env_extra in STATES.items():
            print(f"--- {bid} state {state} ---", flush=True)
            outs[state], prs[state] = _run(bid, state, env_extra, p, root)

        giant = _giant_mrna_ids(bid, p)
        wb, wpid = _mrna_blocks(outs["windowed"])
        fb, fpid = _mrna_blocks(outs["fulldp"])
        all_tids = set(wb) | set(fb)
        differing = [t for t in all_tids if wb.get(t) != fb.get(t)]
        unexpected = [t for t in differing if t not in giant]
        giant_pid_delta = sorted(
            ((t, round(wpid[t] - fpid[t], 6)) for t in giant if t in wpid and t in fpid),
            key=lambda kv: -abs(kv[1]))

        rec = {
            "benchmark": bid,
            "wall_s": {s: round(prs[s].wall_clock_seconds, 1) for s in STATES},
            "peak_rss_mb": {s: round(getattr(prs[s], "peak_rss_mb", 0) or 0, 0) for s in STATES},
            "n_giant_transcripts": len(giant),
            "n_differing_transcripts": len(differing),
            "n_unexpected_diffs_nongiant": len(unexpected),
            "unexpected_diff_ids": unexpected[:20],
            "max_giant_pid_delta": (giant_pid_delta[0] if giant_pid_delta else None),
            "giant_pid_deltas_top": giant_pid_delta[:10],
            "speedup": (round(prs["fulldp"].wall_clock_seconds /
                              max(prs["windowed"].wall_clock_seconds, 1e-9), 1)),
            "rss_reduction": (round((getattr(prs["fulldp"], "peak_rss_mb", 0) or 0) /
                                    max(getattr(prs["windowed"], "peak_rss_mb", 0) or 1, 1), 1)),
        }
        results.append(rec)
        print(f"\n  [{bid}] wall windowed={rec['wall_s']['windowed']}s fulldp={rec['wall_s']['fulldp']}s "
              f"({rec['speedup']}x) | RSS windowed={rec['peak_rss_mb']['windowed']}MB "
              f"fulldp={rec['peak_rss_mb']['fulldp']}MB ({rec['rss_reduction']}x)", flush=True)
        print(f"  giants={rec['n_giant_transcripts']} differing={rec['n_differing_transcripts']} "
              f"unexpected_nongiant_diffs={rec['n_unexpected_diffs_nongiant']} "
              f"max_giant_pid_delta={rec['max_giant_pid_delta']}", flush=True)

    (HERE / "align_window_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "align_window_ab.md")
    print("\nDONE: wrote align_window_ab.json + align_window_ab.md", flush=True)
    # non-zero exit if any non-giant transcript changed (a real regression)
    return 0 if all(r["n_unexpected_diffs_nongiant"] == 0 for r in results) else 1


def _write_md(results, path):
    lines = [
        "## Giant-gene windowed alignment: windowed (default) vs forced full-DP A/B\n",
        "Same cached `-L`/`-M`, single-threaded, differing only by "
        "`LIFTON_ALIGN_WINDOW_{AA,NT}`. Acceptance: large speed + RSS reduction; "
        "**0 non-giant transcripts differ** (normal genes byte-identical); "
        "per-giant protein_identity delta ~0.\n",
        "| Dataset | wall windowed/fulldp (×) | peak RSS windowed/fulldp (×) | giants | non-giant diffs | max giant PI Δ |",
        "|---|---|---|---|---|---|",
    ]
    for r in results:
        w, rss = r["wall_s"], r["peak_rss_mb"]
        lines.append(
            "| {} | {}s / {}s ({}×) | {}MB / {}MB ({}×) | {} | {} | {} |".format(
                r["benchmark"], w["windowed"], w["fulldp"], r["speedup"],
                rss["windowed"], rss["fulldp"], r["rss_reduction"],
                r["n_giant_transcripts"], r["n_unexpected_diffs_nongiant"],
                r["max_giant_pid_delta"]))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
