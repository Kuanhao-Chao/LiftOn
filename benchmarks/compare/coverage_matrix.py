#!/usr/bin/env python
"""Audit coverage of the 8 examples on https://ccb.jhu.edu/lifton/ against the two
benchmark registries (subset `benchmarks.json` + full-genome `datasets.json`) and
the on-disk inputs. Emits `coverage_matrix.md`.

Usage (repo root): python -m benchmarks.compare.coverage_matrix
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
DATA = HERE.parent / "data"
BENCH = json.loads((HERE / "benchmarks.json").read_text())
DATASETS = json.loads((HERE.parent / "datasets.json").read_text())
OUT = HERE / "coverage_matrix.md"

# canonical 8 website examples → (subset benchmarks.json id, full datasets.json id)
WEBSITE = [
    ("Human GRCh38 → T2T-CHM13", "same-species", "human_mane", "human",
     "content/same_species_liftover/liftover_GRCh38_2_T2TCHM13.html"),
    ("Mouse GRCm39 → NOD/SCID", "same-species", "mouse", "mouse",
     "content/same_species_liftover/liftover_mouse.html"),
    ("Honey bee HAv3.1 → ASM1932182v1", "same-species", "bee", "bee",
     "content/same_species_liftover/liftover_bee_insect.html"),
    ("Arabidopsis TAIR10 → ASM2311539v1", "same-species", "arabidopsis", "arabidopsis",
     "content/same_species_liftover/liftover_arabidopsis_plant.html"),
    ("Rice IRGSP → ASM3414082v1", "same-species", "rice", "rice",
     "content/same_species_liftover/liftover_rice_plant.html"),
    ("Human GRCh38 → Chimpanzee", "close cross-species", "human_to_chimp", "human_to_chimp",
     "content/close_species_liftover/liftover_GRCh38_2_chimp.html"),
    ("D. melanogaster → D. erecta", "distant cross-species", "drosophila", "drosophila",
     "content/distant_species_liftover/liftover_drosophila_erecta.html"),
    ("Mouse GRCm39 → Rat mRatBN7.2", "distant cross-species", "mouse_to_rat", "mouse_to_rat",
     "content/distant_species_liftover/liftover_mouse_2_rat.html"),
]


def _bench(bid):
    return next((b for b in BENCH["benchmarks"] if b["id"] == bid), None)


def _dataset(did):
    return next((d for d in DATASETS.get("datasets", []) if d.get("id") == did), None)


def _readable(p):
    try:
        return p and Path(p).exists() and Path(p).stat().st_size > 0
    except OSError:
        return False


def _subset_cached(bid):
    w = WORK / bid
    return {
        "subset": _readable(w / "subset" / "ref.chrom.gff3"),
        "liftoff": _readable(w / "tools" / "liftoff" / "liftoff.gff3"),
        "miniprot": _readable(w / "tools" / "miniprot" / "miniprot.gff3"),
    }


def _yn(b):
    return "✅" if b else "—"


def main():
    L = ["# LiftOn website-example coverage matrix\n",
         "Audit of the 8 examples on https://ccb.jhu.edu/lifton/ against the subset registry "
         "(`benchmarks.json`), the full-genome FTP registry (`datasets.json`), and on-disk inputs.\n"]

    L.append("| # | website example | class | subset id | in benchmarks.json | full id | in datasets.json | ref inputs on disk | subset built (S/L/M) |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    all_ok = True
    for i, (name, cls, bid, did, _url) in enumerate(WEBSITE, 1):
        b = _bench(bid)
        d = _dataset(did)
        inputs_ok = bool(b) and all(_readable(b.get(k)) for k in ("ref_genome", "ref_gff", "tgt_genome"))
        cache = _subset_cached(bid) if b else {"subset": False, "liftoff": False, "miniprot": False}
        slm = f"{_yn(cache['subset'])}{_yn(cache['liftoff'])}{_yn(cache['miniprot'])}"
        if not (b and inputs_ok):
            all_ok = False
        L.append("| {} | {} | {} | `{}` | {} | `{}` | {} | {} | {} |".format(
            i, name, cls, bid, _yn(bool(b)), did, _yn(bool(d)), _yn(inputs_ok), slm))
    L.append("")
    L.append(f"**All 8 website examples covered by a subset benchmark with on-disk inputs:** "
             f"{'YES' if all_ok else 'NO — see gaps above'}.\n")

    # datasets.json membership detail (the registry-completeness fix)
    L.append("## Full-genome FTP registry (`datasets.json`) membership\n")
    have = {d.get("id") for d in DATASETS.get("datasets", [])}
    L.append("| website full id | present in datasets.json |")
    L.append("|---|---|")
    for _n, _c, _bid, did, _u in WEBSITE:
        L.append(f"| `{did}` | {_yn(did in have)} |")
    L.append("")
    missing = [did for _n, _c, _b, did, _u in WEBSITE if did not in have]
    L.append("Missing from `datasets.json`: " +
             (", ".join(f"`{m}`" for m in missing) if missing else "none — registry complete.") + "\n")

    OUT.write_text("\n".join(L) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
