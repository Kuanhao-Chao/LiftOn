#!/usr/bin/env python
"""Iteration-15 feasibility diagnostic — Step-8 weak-Liftoff rescue headroom (Gap #2).

AUDIT RECORD (Iteration-15 NO-GO). This headroom looked promising (drosophila 44
"recoverable" at mean mp_aa 0.925 / weak local Liftoff PI 0.68; mouse_to_rat 8),
so the feature was BUILT (`--miniprot-rescue-weak-liftoff`) and A/B'd
(`miniprot_only_rescue_ab.py`). The A/B then REVEALED what this upside-only
projection could not: the recoverable population is **overwhelmingly DUPLICATES**
of ref genes Liftoff already lifted weakly — only **5 (drosophila) / 0
(mouse_to_rat)** are genuinely-new genes; the rest emit redundant overlapping
gene models (44 / 8 extra). The feature + the byte-neutral `LIFTON_MINIPROT_RESCUE_DIAG`
instrumentation this driver relies on were **REVERTED** (NO-GO). This script +
its committed `miniprot_only_rescue_headroom.{json,md}` are the audit trail;
re-running requires restoring the reverted instrumentation. See the CLAUDE.md
Iteration-15 note + memory.

Read-only / BYTE-NEUTRAL. Step 8 (`run_miniprot.process_miniprot`) drops a
miniprot-only mRNA whenever it overlaps a Liftoff locus (`is_overlapped`),
REGARDLESS of how poor that Liftoff lift is. This is distinct from the Iter-13
NO-GO (which only relaxed the overlap/length THRESHOLDS): the lever here is a
Liftoff-*quality*-aware gate. This driver runs LiftOn ONCE per dataset with
`LIFTON_MINIPROT_RESCUE_DIAG=<path>` (DEFAULT thresholds, output byte-frozen,
cached `-L`/`-M`), which makes `process_miniprot` emit one TSV line per
SUPPRESSED (overlapped) candidate:

  mtid  reason  mp_aa  length_ratio  max_overlap_ratio  ref_gene  ref_trans
        n_cds  n_overlap_genes  locus

`reason` ∈ {PASS, no_ref_gene, no_ref_protein, processed_pseudogene,
length_ratio}; only PASS candidates pass every emit-path filter. `mp_aa` is what
the suppressed mRNA WOULD score (built+scored on a throwaway tree, real state
untouched). We coordinate-join the candidate `locus` against the run's
`score.txt` to recover the **local Liftoff protein identity** (max emitted PI
over Liftoff transcripts whose locus overlaps the candidate). A candidate is
**recoverable** iff it is a clean rescue at a weak Liftoff locus:
    mp_aa >= MP_FLOOR (0.70)  AND  mp_aa >= local_liftoff_pi + MARGIN (0.05)

GO/NO-GO: GO iff n_recoverable >= MIN_RECOVERABLE (10) with mean recoverable
mp_aa >= 0.75 on any dataset; else NO-GO (negligible / well-tuned, document).
Recoverable is the UPPER BOUND; the A/B (`miniprot_only_rescue_ab.py`) confirms
on emitted output, including completeness (off ⊆ on) + no new validity errors.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.miniprot_only_rescue_headroom [IDS...]   # default below
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

DEFAULT_IDS = ["drosophila", "mouse_to_rat"]
MP_FLOOR = 0.70      # below this miniprot is garbage (the Iter-13 garbage floor)
MARGIN = 0.05        # miniprot must beat local Liftoff by this to be a real win
MIN_RECOVERABLE = 10
MEAN_MP_BAR = 0.75


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run_diag(bid: str, p: dict, root: Path, log=print):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    out = root / "default.gff3"
    diag = root / "mp_rescue_diag.tsv"
    if diag.exists():
        diag.unlink()
    _clean_input_dbs(p["ref_gff"], liftoff, miniprot)
    env = _compose_env(TOOLS)
    env["LIFTON_MINIPROT_RESCUE_DIAG"] = str(diag)
    argv = [TOOLS["lifton_bin"], "-t", "8", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), "--native", "--locus-pipeline",
            p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label=f"mp_rescue_headroom_{bid}",
                      log_dir=root / "logs", env=env, log=log)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: lifton failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return diag, root / "lifton_output" / "score.txt"


def _parse_diag(diag: Path):
    """Dedupe by mtid (defensive); keep the first occurrence."""
    rows = {}
    if not diag.exists():
        return []
    for line in diag.read_text().splitlines():
        c = line.split("\t")
        if len(c) < 10:
            continue
        mtid = c[0]
        if mtid in rows:
            continue

        def _f(x):
            try:
                return float(x)
            except (ValueError, TypeError):
                return None
        rows[mtid] = {"mtid": mtid, "reason": c[1], "mp_aa": _f(c[2]),
                      "len_ratio": _f(c[3]), "max_ovp": _f(c[4]),
                      "ref_gene": c[5], "ref_trans": c[6], "n_cds": c[7],
                      "n_ovp_genes": _f(c[8]), "locus": c[9]}
    return list(rows.values())


def _parse_locus(locus: str):
    """`seqid:start-end` -> (seqid, start, end). Returns None on malformed."""
    try:
        seqid, rest = locus.rsplit(":", 1)
        s, e = rest.split("-")
        return seqid, int(s), int(e)
    except (ValueError, AttributeError):
        return None


def _build_liftoff_pi_index(score_txt: Path):
    """Per-seqid list of (start, end, lifton_aa) from score.txt for the local-PI
    coordinate join (col4 = lifton_aa, col7 = seqid:start-end)."""
    idx = {}
    if not score_txt.exists():
        return idx
    for ln in score_txt.read_text().splitlines():
        c = ln.rstrip("\n").split("\t")
        if len(c) < 8:
            continue
        try:
            pi = float(c[4])
        except (ValueError, IndexError):
            continue
        loc = _parse_locus(c[7])
        if loc is None:
            continue
        seqid, s, e = loc
        idx.setdefault(seqid, []).append((s, e, pi))
    return idx


def _local_liftoff_pi(idx, locus):
    """Max emitted Liftoff PI over score.txt loci overlapping the candidate."""
    loc = _parse_locus(locus)
    if loc is None:
        return None
    seqid, s, e = loc
    best = None
    for (ls, le, pi) in idx.get(seqid, []):
        if min(le, e) - max(ls, s) + 1 > 0:  # overlap
            best = pi if best is None else max(best, pi)
    return best


def _mean(xs):
    return round(sum(xs) / len(xs), 5) if xs else 0.0


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: weak-Liftoff miniprot-rescue headroom (Gap #2) ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_mp_rescue_headroom"
        root.mkdir(parents=True, exist_ok=True)

        diag, score_txt = _run_diag(bid, p, root)
        rows = _parse_diag(diag)
        pi_idx = _build_liftoff_pi_index(score_txt)

        reasons = ["PASS", "no_ref_gene", "no_ref_protein",
                   "processed_pseudogene", "length_ratio"]
        by_reason = {r: sum(1 for x in rows if x["reason"] == r) for r in reasons}
        passed = [r for r in rows if r["reason"] == "PASS" and r["mp_aa"] is not None]

        # mp_aa distribution over PASS candidates.
        mp_buckets = {
            "<0.5": sum(1 for r in passed if r["mp_aa"] < 0.5),
            "[0.5,0.7)": sum(1 for r in passed if 0.5 <= r["mp_aa"] < 0.7),
            "[0.7,0.85)": sum(1 for r in passed if 0.7 <= r["mp_aa"] < 0.85),
            ">=0.85": sum(1 for r in passed if r["mp_aa"] >= 0.85),
        }

        recoverable = []
        for r in passed:
            local_pi = _local_liftoff_pi(pi_idx, r["locus"])
            r["local_liftoff_pi"] = local_pi
            if (r["mp_aa"] >= MP_FLOOR
                    and (local_pi is None or r["mp_aa"] >= local_pi + MARGIN)):
                recoverable.append(r)

        rec = {
            "benchmark": bid,
            "n_suppressed_candidates": len(rows),
            "by_reason": by_reason,
            "n_pass_filters": len(passed),
            "mp_aa_buckets": mp_buckets,
            "n_recoverable": len(recoverable),
            "mean_mp_aa_recoverable": _mean([r["mp_aa"] for r in recoverable]),
            "mean_local_liftoff_pi_recoverable": _mean(
                [r["local_liftoff_pi"] for r in recoverable
                 if r["local_liftoff_pi"] is not None]),
        }
        verdict = ("GO (worth an A/B)"
                   if (rec["n_recoverable"] >= MIN_RECOVERABLE
                       and rec["mean_mp_aa_recoverable"] >= MEAN_MP_BAR)
                   else "NO-GO / negligible")
        rec["verdict"] = verdict
        results.append(rec)
        print(f"  [{bid}] suppressed={len(rows)} pass-filters={len(passed)} "
              f"by_reason={by_reason}", flush=True)
        print(f"    mp_aa buckets (PASS) = {mp_buckets}", flush=True)
        print(f"    RECOVERABLE (mp>={MP_FLOOR} & mp>=localLiftoff+{MARGIN}) = "
              f"{rec['n_recoverable']}; mean mp_aa = "
              f"{rec['mean_mp_aa_recoverable']}; mean local Liftoff PI = "
              f"{rec['mean_local_liftoff_pi_recoverable']}", flush=True)
        print(f"    --> {verdict}", flush=True)

    (HERE / "miniprot_only_rescue_headroom.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "miniprot_only_rescue_headroom.md")
    any_go = any(r["n_recoverable"] >= MIN_RECOVERABLE
                 and r["mean_mp_aa_recoverable"] >= MEAN_MP_BAR for r in results)
    print(f"\nOVERALL: {'GO' if any_go else 'NO-GO'} "
          f"(GO if n_recoverable >= {MIN_RECOVERABLE} & mean mp_aa >= "
          f"{MEAN_MP_BAR} on any dataset).", flush=True)
    print("DONE: wrote miniprot_only_rescue_headroom.json + "
          "miniprot_only_rescue_headroom.md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Weak-Liftoff miniprot-rescue headroom (Iteration 15, Gap #2, read-only)\n",
        "Default run on cached `-L`/`-M`, `LIFTON_MINIPROT_RESCUE_DIAG` set "
        "(byte-frozen output). Per SUPPRESSED (overlapped) miniprot mRNA: the "
        "drop reason + what it WOULD score (`mp_aa`), coordinate-joined against "
        "`score.txt` for the **local Liftoff PI**. **recoverable** = a clean "
        f"rescue at a weak locus: mp_aa >= {MP_FLOOR} AND mp_aa >= local Liftoff "
        f"PI + {MARGIN}. Distinct from the Iter-13 threshold lever (this is "
        f"Liftoff-quality-aware). GO if recoverable >= {MIN_RECOVERABLE} & mean "
        f"mp_aa >= {MEAN_MP_BAR} on any dataset; else NO-GO.\n",
        "| Dataset | suppressed | pass-filters | mp_aa buckets | recoverable | mean mp_aa (rec) | mean local Liftoff PI (rec) | verdict |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {} |".format(
            r["benchmark"], r["n_suppressed_candidates"], r["n_pass_filters"],
            r["mp_aa_buckets"], r["n_recoverable"], r["mean_mp_aa_recoverable"],
            r["mean_local_liftoff_pi_recoverable"], r["verdict"]))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
