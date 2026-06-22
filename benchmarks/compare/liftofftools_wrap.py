"""LiftoffTools ``variants`` cross-check (Liftoff & LiftOn only).

miniprot's GFF lacks the gene->mRNA->CDS hierarchy ``liftofftools variants``
expects, so it is excluded. This is a confirming cross-check, not the source of
truth (the custom evaluator is). Runs with the working env binary (the PATH
``liftofftools`` is broken on a numpy ``intc`` issue).
"""
from __future__ import annotations

import json
import os
import subprocess
from collections import Counter
from pathlib import Path


def _find_variant_effects(out_dir: Path) -> Path | None:
    for name in ("variant_effects", "variant_effects.tsv"):
        p = out_dir / name
        if p.exists():
            return p
    hits = list(out_dir.glob("variant*"))
    return hits[0] if hits else None


def parse_variant_effects(path: Path) -> dict:
    """Parse the variant_effects TSV → class counts + identity columns if present."""
    classes: Counter = Counter()
    dna_ids, aa_ids = [], []
    n = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            n += 1
            # heuristic: last column = variant effect class; cols 3/4 = dna/aa id
            eff = parts[-1].strip()
            classes[eff] += 1
            for idx in (2, 3):
                if len(parts) > idx:
                    try:
                        val = float(parts[idx])
                        (dna_ids if idx == 2 else aa_ids).append(val)
                    except ValueError:
                        pass
    out = {"n": n, "classes": dict(classes)}
    if dna_ids:
        out["dna_id_mean"] = round(sum(dna_ids) / len(dna_ids), 5)
    if aa_ids:
        out["aa_id_mean"] = round(sum(aa_ids) / len(aa_ids), 5)
    out["n_identical"] = classes.get("identical", 0)
    return out


def run_variants(tool: str, manifest: dict, work_dir: Path, tools: dict,
                 force: bool, log=print) -> dict | None:
    p = manifest["paths"]
    tool_gff = work_dir / "tools" / tool / f"{tool}.gff3"
    if not tool_gff.exists():
        return None
    out_dir = work_dir / "eval" / "liftofftools" / tool
    done = work_dir / ".done" / f"liftofftools_{tool}.done"
    res_path = out_dir / "parsed.json"
    if done.exists() and not force and res_path.exists():
        log(f"  [liftofftools:{tool}] cached")
        return json.loads(res_path.read_text())
    out_dir.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "PATH": f"{tools.get('extra_path','')}:{os.environ.get('PATH','')}",
           "PYTHONNOUSERSITE": "1"}
    argv = [
        tools["liftofftools_bin"], "variants",
        "-r", p["ref_fa"], "-t", p["tgt_fa"],
        "-rg", p["ref_gff"], "-tg", str(tool_gff),
        "-dir", str(out_dir), "-force",
    ]
    log(f"  [liftofftools:{tool}] running variants ...")
    proc = subprocess.run(argv, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True)
    (out_dir / "liftofftools.stdout.log").write_text(proc.stdout or "")
    (out_dir / "liftofftools.stderr.log").write_text(proc.stderr or "")
    ve = _find_variant_effects(out_dir)
    if proc.returncode != 0 or ve is None:
        log(f"  [liftofftools:{tool}] unavailable (rc={proc.returncode}); see logs")
        res = {"tool": tool, "available": False, "returncode": proc.returncode}
    else:
        res = {"tool": tool, "available": True, **parse_variant_effects(ve)}
    res_path.write_text(json.dumps(res, indent=2))
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("ok\n")
    return res


def run_crosscheck(manifest: dict, work_dir: Path, tools: dict, force: bool,
                   log=print) -> dict:
    out = {}
    for tool in ("liftoff", "lifton"):
        r = run_variants(tool, manifest, work_dir, tools, force, log)
        if r is not None:
            out[tool] = r
    (work_dir / "eval" / "liftofftools" / "crosscheck.json").write_text(json.dumps(out, indent=2))
    return out
