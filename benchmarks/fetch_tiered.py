#!/usr/bin/env python
"""Fetch the tiered benchmark-pair inputs into a clean tier/pair tree.

Reads ``benchmarks/tiers.json`` and, for each pair, builds::

    <data_root>/<tier>/<id>/{ref.fna, ref.gff, ref.faa?, tgt.fna}

with canonical symlink names so ``benchmarks/compare/benchmarks.json`` can point
at fixed paths regardless of the messy NCBI filenames (same convention as
``benchmarks/compare/fetch_new_pairs.py``, whose download helpers this script
reuses).

Robustness:
  * A reference/target already on disk (a pair's ``ref.local_fna`` /
    ``tgt.local_fna`` hint) is **symlinked** instead of re-downloaded; an absent
    or unreadable hint falls back to the accession.
  * Each unique accession is fetched **once** into ``<data_root>/_cache/<acc>/``
    and symlinked into every pair that uses it (so a shared genome — e.g. the
    tomato Micro-Tom reference used by two pairs — downloads a single time).
  * ``_ncbi_download`` (reused) already retries 3x; a **per-pair** failure is
    logged and skipped (not fatal), so one bad accession never aborts the run.
  * A ``manifest.json`` at the data_root records source / byte-size / status for
    every canonical file.
  * Idempotent — re-running re-links from the (cached) downloads without
    re-fetching.

Usage (repo root, lifton_devel env)::

    python -m benchmarks.fetch_tiered                 # all 12 pairs
    python -m benchmarks.fetch_tiered t2_human_to_gorilla
    python -m benchmarks.fetch_tiered close_cross_species   # a whole tier
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

# Reuse the proven downloader helpers from the sibling compare/ script.
from benchmarks.compare.fetch_new_pairs import (
    _glob1,
    _link,
    _ncbi_download,
)

REGISTRY = Path(__file__).resolve().parent / "tiers.json"


def load_registry(path: Path = REGISTRY) -> dict:
    with open(path) as fh:
        return json.load(fh)


def _resolved_size(p: Path) -> int:
    """Byte-size of a (possibly symlinked) file, following the link."""
    try:
        return os.path.getsize(p)  # stat() follows symlinks
    except OSError:
        return -1


def _cache_dir(data_root: Path) -> Path:
    d = data_root / "_cache"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download_acc(cache: Path, acc: str, include: str) -> Path:
    """Download ``acc`` once into the per-accession cache; return its data dir.

    No accession in tiers.json is used as both a reference (needs gff3+protein)
    and a target (genome only), so a single cache entry per accession is safe.
    """
    return _ncbi_download(cache, acc, acc, include)


def _link_optional(src, dst) -> bool:
    """Symlink ``dst`` -> ``src`` if ``src`` resolves; return True on success."""
    try:
        _link(src, dst)
        return True
    except FileNotFoundError:
        return False


def fetch_pair(pair: dict, data_root: Path) -> dict:
    """Populate ``<data_root>/<tier>/<id>/`` for one pair; return its manifest."""
    pid, tier = pair["id"], pair["tier"]
    dest = data_root / tier / pid
    dest.mkdir(parents=True, exist_ok=True)
    cache = _cache_dir(data_root)
    files: dict = {}

    # --- reference: genome + annotation (+ optional proteins) ---
    ref = pair["ref"]
    if ref.get("local_fna") and Path(ref["local_fna"]).exists():
        _link(ref["local_fna"], dest / "ref.fna")
        _link(ref["local_gff"], dest / "ref.gff")
        files["ref.fna"] = {"source": "local", "target": ref["local_fna"]}
        files["ref.gff"] = {"source": "local", "target": ref["local_gff"]}
        if ref.get("local_faa") and _link_optional(ref["local_faa"], dest / "ref.faa"):
            files["ref.faa"] = {"source": "local", "target": ref["local_faa"]}
    else:
        d = _download_acc(cache, ref["acc"], "genome,gff3,protein")
        _link(_glob1(str(d / "*_genomic.fna")), dest / "ref.fna")
        _link(_glob1(str(d / "genomic.gff")), dest / "ref.gff")
        files["ref.fna"] = {"source": f"ncbi:{ref['acc']}"}
        files["ref.gff"] = {"source": f"ncbi:{ref['acc']}"}
        try:
            _link(_glob1(str(d / "protein.faa")), dest / "ref.faa")
            files["ref.faa"] = {"source": f"ncbi:{ref['acc']}"}
        except FileNotFoundError:
            print("  (no protein.faa — transcript-space miniprot translates anyway)",
                  flush=True)

    # --- target: genome only ---
    tgt = pair["tgt"]
    if tgt.get("local_fna") and Path(tgt["local_fna"]).exists():
        _link(tgt["local_fna"], dest / "tgt.fna")
        files["tgt.fna"] = {"source": "local", "target": tgt["local_fna"]}
    else:
        d = _download_acc(cache, tgt["acc"], "genome")
        _link(_glob1(str(d / "*_genomic.fna")), dest / "tgt.fna")
        files["tgt.fna"] = {"source": f"ncbi:{tgt['acc']}"}

    # --- record resolved sizes; required files must be non-empty ---
    required = ["ref.fna", "ref.gff", "tgt.fna"]
    ok = True
    for name, meta in files.items():
        size = _resolved_size(dest / name)
        meta["bytes"] = size
        if name in required and size <= 0:
            ok = False
    for name in required:
        if name not in files:
            ok = False

    return {"tier": tier, "dir": str(dest), "status": "ok" if ok else "incomplete",
            "files": files}


def main(argv=None) -> int:
    argv = list(argv) if argv is not None else []
    reg = load_registry()
    data_root = Path(reg["data_root"])
    data_root.mkdir(parents=True, exist_ok=True)
    pairs = reg["pairs"]

    if argv:
        sel = set(argv)
        pairs = [p for p in pairs if p["id"] in sel or p["tier"] in sel]
        if not pairs:
            print(f"no pairs matched {argv}; known ids: {[p['id'] for p in reg['pairs']]}",
                  flush=True)
            return 2

    manifest: dict = {}
    failed: list = []
    for pair in pairs:
        pid = pair["id"]
        print(f"\n=== fetch {pid}  [{pair['tier']}] ===", flush=True)
        try:
            entry = fetch_pair(pair, data_root)
            manifest[pid] = entry
            if entry["status"] != "ok":
                failed.append(pid)
                print(f"  !! {pid} incomplete (a required file is missing/empty)",
                      flush=True)
            else:
                print(f"  ok -> {entry['dir']}", flush=True)
        except Exception as e:   # noqa: BLE001
            import traceback
            manifest[pid] = {"tier": pair["tier"], "status": "failed", "error": str(e)}
            failed.append(pid)
            print(f"  !! {pid} FETCH FAILED: {e}\n{traceback.format_exc()}", flush=True)

    # Merge into any existing manifest so a partial re-fetch (one pair / one
    # tier) updates just its entries instead of clobbering the whole record.
    manifest_path = data_root / "manifest.json"
    merged: dict = {}
    if manifest_path.exists():
        try:
            merged = json.load(open(manifest_path)).get("pairs", {})
        except Exception:   # noqa: BLE001
            merged = {}
    merged.update(manifest)
    all_failed = [pid for pid, e in merged.items() if e.get("status") != "ok"]
    with open(manifest_path, "w") as fh:
        json.dump({"pairs": merged,
                   "summary": {"total": len(merged),
                               "ok": len(merged) - len(all_failed),
                               "failed": all_failed}}, fh, indent=2)
    print(f"\nfetch done (this run): ok={len(pairs) - len(failed)}/{len(pairs)} "
          f"failed={failed}\nmanifest ({len(merged)} pairs total) -> {manifest_path}",
          flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or None))
