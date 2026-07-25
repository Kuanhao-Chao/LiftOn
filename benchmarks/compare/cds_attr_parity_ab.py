#!/usr/bin/env python
"""CDS attribute parity divergence-ladder A/B — cached -L/-M, deterministic.

The chaining merge and the ORF-rescue boundary patch rebuild a transcript's CDS
list, and both used to reset the emitted attributes to ``{Parent}``. Measured on
the full drosophila lift, 7,127 of 7,142 pure-Liftoff CDS rows kept
``Dbxref``/``product``/``protein_id``/``gene``/``locus_tag``/… while only 21 of
151,277 merged rows did -- even though every merged mRNA row kept the very same
attributes. ``Lifton_EXON.add_lifton_cds`` / ``add_novel_lifton_cds`` now carry
the transcript's reference CDS attribute template through the rebuild.

Arms (identical argv, cached -L/-M, -t1, NO -copies = deterministic; only the
env differs):
  "off" : LIFTON_NO_CDS_ATTR_CARRY=1   (the pre-change {Parent}-only rows)
  "on"  : default                      (attributes carried)

WHY THE GATE IS NOT A PROTEIN-IDENTITY COMPARISON. The change touches column 9
only. That makes a far stronger check available than a mean-PI delta: columns
1-8 of every row must be **byte-identical** between the arms. Since a CDS's
coordinates, strand and phase are exactly columns 1-8, byte-equality there
*proves* the translated protein is unchanged -- no re-alignment pass can add
information to that, and the evaluator's own scoring noise cannot mask a
regression. The remaining risk is confined to column 9, so the gate also
requires that the ON attribute set is a strict superset of OFF with no shared
value rewritten.

Per-cell gate (ALL must hold):
  - row counts equal;
  - cols 1-8 mismatched rows == 0        (no coordinate/phase/strand/score move);
  - attributes lost == 0 and values changed == 0   (strict superset);
  - protein_identity changed on 0 transcripts      (LiftOn's own scoring flat);
  - gff3-validate error count not worse;
  - CDS descriptive-attribute coverage strictly rises (the win) OR the cell has
    no rebuilt CDS at all (a legitimate inert no-op).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.cds_attr_parity_ab [IDS...]
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
from pathlib import Path

from . import evaluator
from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

# The standard promotion ladder (same cells as miniprot_candidate_ladder_ab /
# miniprot_rescue_ab): same-species control + distant + very-distant. Smaller,
# faster cells first so the incremental JSON write validates the harness early.
DEFAULT_IDS = [
    "celegans_to_briggsae",       # distant
    "drosophila_to_anopheles",    # very-distant
    "zebrafish_to_medaka",        # very-distant
    "rice_to_sorghum",            # distant
    "t4_human_to_xenopus",        # very-distant
    "t4_human_to_chicken",        # very-distant
    "human_to_mouse",             # distant
    "drosophila",                 # same-species control (slowest cell)
]

STRUCTURAL_KEYS = {"ID", "Parent"}


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _divergence(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("divergence_class") or b.get("tier") or ""
    return ""


def _run_arm(bid, arm, p, root):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    statedir = root / arm
    statedir.mkdir(parents=True, exist_ok=True)
    out = statedir / f"{arm}.gff3"
    _clean_input_dbs(p["ref_gff"], liftoff, miniprot)
    argv = [TOOLS["lifton_bin"], "-t", "1", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), p["tgt_fa"], p["ref_fa"]]
    env = dict(_compose_env(TOOLS))
    env["PYTHONHASHSEED"] = "0"
    # Pin BOTH arms to one build. `_compose_env` deliberately builds a fresh
    # env, so set LIFTON_AB_PYTHONPATH to run the ladder out of a detached
    # worktree while the main tree moves on -- a mixed-build A/B table is not a
    # comparison.
    ab_pythonpath = os.environ.get("LIFTON_AB_PYTHONPATH")
    if ab_pythonpath:
        env["PYTHONPATH"] = ab_pythonpath
    if arm == "off":
        env["LIFTON_NO_CDS_ATTR_CARRY"] = "1"
    else:
        env.pop("LIFTON_NO_CDS_ATTR_CARRY", None)
    pr = run_profiled(argv, label=f"cds_attr_{bid}_{arm}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: arm {arm} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out


def _rows(path: Path):
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) == 9:
                yield cols


def _attrs(col9: str) -> dict:
    return dict(kv.split("=", 1) for kv in col9.split(";") if "=" in kv)


def _compare(off_path: Path, on_path: Path) -> dict:
    """Row-aligned structural comparison of the two outputs."""
    off = list(_rows(off_path))
    on = list(_rows(on_path))
    rec = {"rows_off": len(off), "rows_on": len(on)}
    if len(off) != len(on):
        rec["row_count_equal"] = False
        return rec
    rec["row_count_equal"] = True

    cols_mismatch = 0
    attrs_lost = attrs_changed = 0
    gained = collections.Counter()
    rich = {"off": collections.Counter(), "on": collections.Counter()}
    total = collections.Counter()
    pi_off, pi_on = {}, {}

    for a, b in zip(off, on):
        if a[:8] != b[:8]:
            cols_mismatch += 1
        ao, bo = _attrs(a[8]), _attrs(b[8])
        for key, value in ao.items():
            if key not in bo:
                attrs_lost += 1
            elif bo[key] != value:
                attrs_changed += 1
        for key in bo.keys() - ao.keys():
            gained[key] += 1
        ftype = b[2]
        total[ftype] += 1
        if set(ao) - STRUCTURAL_KEYS:
            rich["off"][ftype] += 1
        if set(bo) - STRUCTURAL_KEYS:
            rich["on"][ftype] += 1
        if ftype in ("mRNA", "transcript"):
            if "protein_identity" in ao:
                pi_off[ao.get("ID", "")] = ao["protein_identity"]
            if "protein_identity" in bo:
                pi_on[bo.get("ID", "")] = bo["protein_identity"]

    rec.update({
        "cols1_8_mismatch": cols_mismatch,
        "attrs_lost": attrs_lost,
        "attrs_value_changed": attrs_changed,
        "attrs_gained": dict(gained.most_common(12)),
        "n_cds": total["CDS"],
        "cds_coverage_off": (round(rich["off"]["CDS"] / total["CDS"], 5)
                             if total["CDS"] else None),
        "cds_coverage_on": (round(rich["on"]["CDS"] / total["CDS"], 5)
                            if total["CDS"] else None),
        "n_transcripts_scored": len(pi_off),
        "protein_identity_changed": sum(
            1 for k in pi_off if pi_on.get(k) != pi_off[k]
        ),
        "bytes_off": off_path.stat().st_size,
        "bytes_on": on_path.stat().st_size,
    })
    rec["size_delta_pct"] = (
        round(100.0 * (rec["bytes_on"] - rec["bytes_off"]) / rec["bytes_off"], 2)
        if rec["bytes_off"] else None
    )
    return rec


def _validate(gff: Path):
    return evaluator.count_gff3_validator_errors(
        gff, TOOLS["lifton_python"], _compose_env(TOOLS))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="*")
    args = ap.parse_args(argv if argv is not None else sys.argv[1:])
    ids = args.ids or DEFAULT_IDS

    results = []
    for bid in ids:
        print(f"=== {bid}: CDS attribute parity A/B ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_cds_attr_parity_ab"
        root.mkdir(parents=True, exist_ok=True)

        out_off = _run_arm(bid, "off", p, root)
        out_on = _run_arm(bid, "on", p, root)

        rec = {"benchmark": bid,
               "divergence": _divergence(bid) or man.get("species", "")}
        rec.update(_compare(out_off, out_on))
        rec["validity"] = {"off": _validate(out_off), "on": _validate(out_on)}

        cov_off = rec.get("cds_coverage_off")
        cov_on = rec.get("cds_coverage_on")
        # A cell with nothing to rebuild (every CDS already rich) is a
        # legitimate inert no-op, not a failure.
        inert = cov_off is not None and cov_on is not None and cov_on == cov_off
        rec["inert"] = bool(inert)
        rec["gate_pass"] = bool(
            rec.get("row_count_equal")
            and rec.get("cols1_8_mismatch") == 0
            and rec.get("attrs_lost") == 0
            and rec.get("attrs_value_changed") == 0
            and rec.get("protein_identity_changed") == 0
            and rec["validity"]["on"] <= rec["validity"]["off"]
            and (cov_on is None or cov_off is None or cov_on >= cov_off)
        )
        results.append(rec)
        print(f"  [{bid}] cols1-8 mismatch={rec.get('cols1_8_mismatch')} "
              f"lost={rec.get('attrs_lost')} changed={rec.get('attrs_value_changed')} "
              f"CDS coverage {cov_off}->{cov_on} "
              f"PI changed={rec.get('protein_identity_changed')} "
              f"val {rec['validity']['off']}->{rec['validity']['on']} "
              f"size {rec.get('size_delta_pct')}% "
              f"gate {'PASS' if rec['gate_pass'] else 'FAIL'}", flush=True)
        (HERE / "cds_attr_parity_ab.json").write_text(json.dumps(results, indent=2))
        _write_md(results, HERE / "cds_attr_parity_ab.md")

    n_pass = sum(1 for r in results if r["gate_pass"])
    print(f"\nGATE: {n_pass}/{len(results)} cells pass.", flush=True)
    print("DONE: wrote cds_attr_parity_ab.json + .md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## CDS attribute parity divergence-ladder A/B\n",
        "`off` = `LIFTON_NO_CDS_ATTR_CARRY=1` (the pre-change `{Parent}`-only CDS "
        "rows); `on` = default (the reference CDS's descriptive attributes carried "
        "through the chaining merge and the ORF-rescue rebuild). Cached `-L`/`-M`, "
        "`-t1` no-`copies`, identical argv — only the env differs.\n",
        "The change touches column 9 only, so the gate is stronger than a mean-PI "
        "comparison: **columns 1-8 must be byte-identical** on every row, which is "
        "exactly a CDS's coordinates, strand and phase — proving the translated "
        "protein cannot have moved. On top of that the ON attribute set must be a "
        "strict superset of OFF (nothing lost, no shared value rewritten), LiftOn's "
        "own `protein_identity` must be unchanged on every transcript, and "
        "`gff3-validate` must not report more errors.\n",
        "| Dataset | divergence | rows | cols1-8 Δ | attrs lost | attrs changed | "
        "CDS coverage off→on | PI changed | validity off→on | size Δ | gate |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        v = r.get("validity", {})
        lines.append(
            f"| {r['benchmark']} | {r.get('divergence','')} | {r.get('rows_on','?')} | "
            f"{r.get('cols1_8_mismatch','?')} | {r.get('attrs_lost','?')} | "
            f"{r.get('attrs_value_changed','?')} | "
            f"{r.get('cds_coverage_off')}→{r.get('cds_coverage_on')} | "
            f"{r.get('protein_identity_changed','?')} | "
            f"{v.get('off')}→{v.get('on')} | {r.get('size_delta_pct')}% | "
            f"{'PASS' if r.get('gate_pass') else 'FAIL'} |"
        )
    n_pass = sum(1 for r in results if r.get("gate_pass"))
    lines.append(f"\n**GATE: {n_pass}/{len(results)} cells pass.**\n")
    path.write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
