#!/usr/bin/env python
"""Download + canonicalize the NEW benchmark-pair inputs into
``/ccb/salz3/kh.chao/lifton_benchmark_data/<id>/`` with STABLE symlink names
(``ref.fna`` / ``ref.gff`` / ``ref.faa`` / ``tgt.fna``; ``ref.gff3`` for the
Ensembl pair) so ``benchmarks.json`` can point at fixed paths regardless of the
messy NCBI/Ensembl filenames. Idempotent — skips a role whose canonical link
already resolves. Used by ``fourway_orchestrate.sh`` Phase A.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.fetch_new_pairs [ids...]   # default: all new pairs
"""
from __future__ import annotations

import glob
import subprocess
import sys
import zipfile
from pathlib import Path

BASE = Path("/ccb/salz3/kh.chao/lifton_benchmark_data")
DATASETS_BIN = "/ccb/sw/bin/datasets"
GFFREAD = "/ccb/sw/bin/gffread"

# NCBI pairs: which roles to download (ref reused from on-disk where absent).
NCBI = {
    "fly_mel_to_pseudoobscura": {"tgt": "GCF_009870125.1"},          # ref = on-disk d.melanogaster
    "yeast_cerevisiae_to_paradoxus": {"ref": "GCF_000146045.2", "tgt": "GCF_002079055.1"},
    "candida_albicans_to_dubliniensis": {"ref": "GCF_000182965.3", "tgt": "GCF_000026945.1"},
    "celegans_to_briggsae": {"ref": "GCF_000002985.6", "tgt": "GCF_000004555.2"},
    "chicken_to_quail": {"ref": "GCF_016699485.2", "tgt": "GCF_001577835.2"},
    # --- distant cross-species ladders (2026-06-18): plant / fish / insect / fungi ---
    "arabidopsis_to_lyrata": {"tgt": "GCF_000004255.2"},             # ref = on-disk TAIR10
    "rice_to_sorghum": {"tgt": "GCF_000003195.3"},                   # ref = on-disk rice IRGSP
    "zebrafish_to_medaka": {"ref": "GCF_000002035.6", "tgt": "GCF_002234675.1"},  # zebrafish ref also = human_to_zebrafish target
    "drosophila_to_anopheles": {"tgt": "GCF_000005575.2"},           # ref = on-disk d.melanogaster
    "cerevisiae_to_pombe": {"tgt": "GCF_000002945.2"},               # ref = on-disk yeast_cerevisiae ref
}
ENSEMBL = {
    "human_ensembl_gtf": {
        "fasta": "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/"
                 "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        "gtf": "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/"
               "Homo_sapiens.GRCh38.112.gtf.gz",
    },
}
# reuse-only pairs (no download) — listed so the driver can no-op them cleanly
# arabidopsis_to_rice: both genomes on-disk (TAIR10 ref + rice IRGSP target).
# human_to_zebrafish: on-disk human ref + the zebrafish genome fetched for
# zebrafish_to_medaka (so fetch zebrafish_to_medaka before building this pair).
REUSE = {"human_to_mouse", "human_pseudogene_stress",
         "arabidopsis_to_rice", "human_to_zebrafish"}


def _run(argv, **kw):
    print("  $ " + " ".join(map(str, argv)), flush=True)
    subprocess.run(argv, check=True, **kw)


def _link(src, dst):
    src, dst = Path(src), Path(dst)
    if not src.exists():
        raise FileNotFoundError(src)
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    dst.symlink_to(src)
    print(f"  link {dst.name} -> {src}", flush=True)


def _glob1(pat):
    hits = sorted(glob.glob(pat))
    if not hits:
        raise FileNotFoundError(pat)
    return hits[0]


def _ncbi_download(dest: Path, role: str, acc: str, include: str) -> Path:
    z = dest / f"{role}.zip"
    unz = dest / f"{role}_unz"
    data = unz / "ncbi_dataset" / "data" / acc
    if data.exists() and any(data.iterdir()):
        return data
    last = None
    for attempt in range(1, 4):
        try:
            _run([DATASETS_BIN, "download", "genome", "accession", acc,
                  "--include", include, "--filename", str(z)])
            with zipfile.ZipFile(z) as zf:
                zf.extractall(unz)
            return data
        except Exception as e:   # noqa: BLE001
            last = e
            print(f"  ! {acc} attempt {attempt}/3 failed: {e}", flush=True)
    raise RuntimeError(f"datasets download failed for {acc}: {last}")


def fetch_ncbi(pid: str, spec: dict):
    dest = BASE / pid
    dest.mkdir(parents=True, exist_ok=True)
    if "ref" in spec:
        d = _ncbi_download(dest, "ref", spec["ref"], "genome,gff3,protein")
        _link(_glob1(str(d / "*_genomic.fna")), dest / "ref.fna")
        _link(_glob1(str(d / "genomic.gff")), dest / "ref.gff")
        try:
            _link(_glob1(str(d / "protein.faa")), dest / "ref.faa")
        except FileNotFoundError:
            print("  (no protein.faa — transcript-space miniprot translates anyway)", flush=True)
    if "tgt" in spec:
        d = _ncbi_download(dest, "tgt", spec["tgt"], "genome")
        _link(_glob1(str(d / "*_genomic.fna")), dest / "tgt.fna")


def fetch_ensembl(pid: str, spec: dict):
    dest = BASE / pid
    dest.mkdir(parents=True, exist_ok=True)
    fa, gff3 = dest / "ref.fna", dest / "ref.gff3"
    if not fa.exists():
        gz = dest / "ref.fa.gz"
        _run(["wget", "-q", "-O", str(gz), spec["fasta"]])
        _run(["bash", "-c", f"gunzip -c {gz} > {fa}"])
    if not gff3.exists():
        gz, gtf = dest / "ref.gtf.gz", dest / "ref.gtf"
        _run(["wget", "-q", "-O", str(gz), spec["gtf"]])
        _run(["bash", "-c", f"gunzip -c {gz} > {gtf}"])
        # GTF -> GFF3, keeping gene records + all attributes (biotype etc.).
        _run([GFFREAD, str(gtf), "-F", "--keep-genes", "-o", str(gff3)])


def main(argv=None):
    ids = argv if argv else (list(NCBI) + list(ENSEMBL))
    failed = []
    for pid in ids:
        print(f"\n=== fetch {pid} ===", flush=True)
        try:
            if pid in NCBI:
                fetch_ncbi(pid, NCBI[pid])
            elif pid in ENSEMBL:
                fetch_ensembl(pid, ENSEMBL[pid])
            elif pid in REUSE:
                print("  reuse-only (on-disk) — nothing to download", flush=True)
            else:
                print(f"  ! unknown new pair {pid}", flush=True)
                failed.append(pid)
        except Exception as e:   # noqa: BLE001
            import traceback
            print(f"  !! {pid} FETCH FAILED: {e}\n{traceback.format_exc()}", flush=True)
            failed.append(pid)
    print(f"\nfetch done; failed={failed}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or None))
