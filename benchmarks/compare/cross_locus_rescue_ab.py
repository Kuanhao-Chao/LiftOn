#!/usr/bin/env python
"""Phase D A/B — duplicate-safe CROSS-LOCUS miniprot rescue.

Tests `lifton.cross_locus_rescue` on a cross-locus-dominant very-distant pair
(default human->zebrafish; ~86% of the Figure-4 regressed transcripts are
cross-locus paralog disagreement candidate-3 cannot reach). Step-7-only via the
-T/-P cached-FASTA bypass + cached -L/-M; the cached ref `.gff_db` is opened
directly (no rebuild).

Both arms are -t8 --locus-pipeline, NO -copies (deterministic), default
miniprot-rescue (Iter-23) ON in BOTH so the ONLY difference is the cross-locus
pass:
  off : default                         (cross_locus OFF)
  on  : --miniprot-cross-locus-rescue   (cross_locus ON)

Cross-locus REPLACE improves a weak gene's best protein identity but drops that
gene's other (weak) variants -> a per-TRANSCRIPT recall cost (the synteny /
variant-multiplicity tradeoff). So unlike candidate-3 (strictly additive), this
is NOT expected to be a clean win; the harness reports BOTH axes plus the
load-bearing safety check (ZERO duplicate gene ids) so the promote/keep-opt-in
decision is informed.

Run:  python -m benchmarks.compare.cross_locus_rescue_ab [PAIR]
"""
from __future__ import annotations
import csv, json, sys
from collections import Counter
from pathlib import Path
from . import evaluator
from . import version_compare as vc
from .profiling import run_profiled
from .tool_runners import _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
TOOLS = json.loads((HERE / "benchmarks.json").read_text())["tools"]
_ALLB = json.loads((HERE / "benchmarks.json").read_text())["benchmarks"]

_ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
FORCE = "--force" in sys.argv[1:]
PYTHONPATH = next((a.split("=", 1)[1] for a in sys.argv[1:]
                   if a.startswith("--pythonpath=")), None)
PAIR = _ARGS[0] if _ARGS else "human_to_zebrafish"
_CFG = next(c for c in _ALLB if c.get("id") == PAIR)
REF_GFF = Path(_CFG["ref_gff"])
REF_FA = Path(_CFG["ref_genome"])
TGT_FA = Path(_CFG["tgt_genome"])
FULL = WORK / PAIR / "_fourway_full" / "devel" / "lifton_output"
L = FULL / "liftoff" / "liftoff.gff3"
M = FULL / "miniprot" / "miniprot.gff3"
T = FULL / "intermediate_files" / "transcripts.fa"
P = FULL / "intermediate_files" / "proteins.fa"
MINI_TSV = WORK / PAIR / "_fourway_full" / "eval" / "miniprot.transcripts.tsv"


def _run(arm, root):
    out = root / arm / f"{arm}.gff3"
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.stat().st_size > 1_000_000 and not FORCE:
        print(f"[bench] reusing existing {arm}.gff3 ({out.stat().st_size} bytes)", flush=True)
        return out
    argv = [TOOLS["lifton_bin"], "-t", "8", "--locus-pipeline", "-ad", "RefSeq",
            "-g", str(REF_GFF), "-L", str(L), "-M", str(M),
            "-T", str(T), "-P", str(P), "-o", str(out),
            str(TGT_FA), str(REF_FA)]
    if arm == "on":
        argv.insert(1, "--miniprot-cross-locus-rescue")
    env = dict(_compose_env(TOOLS))
    if PYTHONPATH:
        env["PYTHONPATH"] = PYTHONPATH
    pr = run_profiled(argv, label=f"xlocus_ab_{arm}", log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"arm {arm} failed (exit {pr.exit_code}); see {pr.stderr_path}")
    return out


def _gene_ids_and_tags(gff: Path):
    """Return (Counter of top-level gene-line IDs, n cross_locus-tagged mRNA)."""
    ids = Counter()
    n_tag = 0
    with gff.open() as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.split("\t")
            if len(c) < 9:
                continue
            attrs = c[8]
            if c[2] == "gene" and "Parent=" not in attrs:
                for kv in attrs.split(";"):
                    if kv.startswith("ID="):
                        ids[kv[3:].strip()] += 1
                        break
            if c[2] == "mRNA" and "lifton_rescue=cross_locus" in attrs:
                n_tag += 1
    return ids, n_tag


def _pi_by_ref(tsv: Path) -> dict:
    out = {}
    if not tsv.exists():
        return out
    for r in csv.DictReader(tsv.open(), delimiter="\t"):
        if str(r.get("is_coding")).strip().lower() not in ("1", "true", "yes"):
            continue
        if str(r.get("recovered")).strip().lower() not in ("1", "true", "yes"):
            continue
        pi = r.get("protein_identity")
        if pi in (None, "", "None"):
            continue
        try:
            v = float(pi)
        except ValueError:
            continue
        rid = r["ref_mrna_id"]
        out[rid] = v if rid not in out else max(out[rid], v)
    return out


def _tag_counts(gff: Path):
    """(cross-locus transcripts, of which attached isoforms)."""
    replaced = isoforms = 0
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "mRNA":
                continue
            if "lifton_rescue=cross_locus" in cols[8]:
                replaced += 1
                isoforms += "rescue_isoform=true" in cols[8]
    return replaced, isoforms


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 5) if xs else None


def main():
    root = WORK / PAIR / "_xlocus_ab"
    root.mkdir(parents=True, exist_ok=True)
    for f in (REF_GFF, REF_FA, TGT_FA, L, M, T, P):
        if not Path(f).exists():
            sys.exit(f"missing input: {f}")
    man = {"id": PAIR, "species": _CFG.get("species", PAIR),
           "cross_species": _CFG.get("cross_species", True),
           "paths": {"ref_gff": str(REF_GFF), "ref_fa": str(REF_FA), "tgt_fa": str(TGT_FA)},
           "miniprot_target_space": "transcript", "protein_acc_to_mrna": {}}
    print(f"=== build reference ({PAIR}) ===", flush=True)
    ref, ref_index = evaluator.build_reference(str(REF_GFF), str(REF_FA), log=print)
    eval_dir = root / "eval"; eval_dir.mkdir(parents=True, exist_ok=True)
    pis = {}
    gene_dups = {}
    n_tags = {}
    for arm in ("off", "on"):
        print(f"=== run arm {arm} ===", flush=True)
        out = _run(arm, root)
        ids, n_tag = _gene_ids_and_tags(out)
        gene_dups[arm] = [g for g, n in ids.items() if n > 1]
        n_tags[arm] = n_tag
        evaluator.evaluate_tool(arm, str(out), str(TGT_FA), ref, man, eval_dir, None,
                                log=print, ref_index=ref_index, threads=8)
        pis[arm] = _pi_by_ref(eval_dir / f"{arm}.transcripts.tsv")
    mini = _pi_by_ref(MINI_TSV)

    off, on = pis["off"], pis["on"]
    common = set(off) & set(on)
    improved = sum(1 for k in common if on[k] > off[k] + 1e-9)
    regressed = sum(1 for k in common if on[k] < off[k] - 1e-9)
    lost = set(off) - set(on)        # dropped variants (the recall cost)
    gained = set(on) - set(off)
    cm_off = set(off) & set(mini); cm_on = set(on) & set(mini)
    a2a_off = _mean([off[k] - mini[k] for k in cm_off])
    a2a_on = _mean([on[k] - mini[k] for k in cm_on])
    rec = {
        "pair": PAIR,
        "n_replaced(cross_locus_tag_on)": n_tags["on"],
        "DUPLICATE_GENE_IDS_on": len(gene_dups["on"]),       # MUST be 0 (safety)
        "DUPLICATE_GENE_IDS_off": len(gene_dups["off"]),
        "dup_examples_on": gene_dups["on"][:5],
        # recovery (per recovered coding transcript) -- the recall axis
        "recovered_off": len(off), "recovered_on": len(on),
        "recovered_delta": len(on) - len(off),
        "n_transcripts_lost(off-minus-on)": len(lost),
        "n_transcripts_gained(on-minus-off)": len(gained),
        # per-transcript identity on the common set
        "n_improved_common": improved, "n_regressed_common": regressed,
        "mean_pi_off": _mean(list(off.values())), "mean_pi_on": _mean(list(on.values())),
        # the Figure-4 apples-to-apples deficit vs miniprot
        "apples_to_apples_vs_miniprot_off": a2a_off,
        "apples_to_apples_vs_miniprot_on": a2a_on,
        "deficit_shrunk": (a2a_on is not None and a2a_off is not None and a2a_on > a2a_off),
        # SAFETY gate (duplicate-free) is load-bearing; the rest is the tradeoff
        "DUPLICATE_SAFE": len(gene_dups["on"]) == len(gene_dups["off"]) == 0,
    }
    replaced_on, isoforms_on = _tag_counts(root / "on" / "on.gff3")
    validity = {arm: vc.validate_gff(root / arm / f"{arm}.gff3", log=print)
                for arm in ("off", "on")}
    errors = {arm: validity[arm].get("n_errors") or 0 for arm in ("off", "on")}
    # A replace can only be promoted if it stops costing transcripts: the whole
    # reason it stayed opt-in was -401 net on this pair.
    gate = {
        "duplicate_safe": rec["DUPLICATE_SAFE"],
        "no_transcript_cost": len(on) >= len(off),
        "identity_improved": (rec["mean_pi_on"] or 0) > (rec["mean_pi_off"] or 0),
        "no_regression_on_common": regressed == 0,
        "validity_not_worse": errors["on"] <= errors["off"],
    }
    rec["n_isoforms_attached"] = isoforms_on
    rec["n_cross_locus_transcripts"] = replaced_on
    rec["validity"] = validity
    rec["gate"] = gate
    rec["gate_pass"] = all(gate.values())
    (HERE / f"cross_locus_rescue_ab.{PAIR}.json").write_text(json.dumps(rec, indent=2))
    print(f"\n=== cross-locus rescue A/B ({PAIR}) ===")
    print(json.dumps(rec, indent=2))
    print("done:gate " + ("PASS" if rec["gate_pass"] else "FAIL"))


if __name__ == "__main__":
    main()
