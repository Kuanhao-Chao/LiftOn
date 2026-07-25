"""Verification-gap closer — serial vs threaded byte-identity for Step 8.

`miniprot_pipeline.parallel_step8` is a full threaded dispatcher with its own proxy-DB
materialisation layer, architecturally the twin of `parallel_step7`. But the 24-cell
byte-identity matrix's fixtures yield ~0 Step-8 candidates, so the load-bearing gate
never exercised it: `tests/test_miniprot_pipeline.py` covers the dispatcher with mocks,
and nothing drove the REAL pipeline through it and compared serial against threaded
output.

This fixture makes Step 8 actually fire: gene2 is present in the reference but absent
from the Liftoff lift, and its miniprot model is INSIDE the default
`-min_miniprot`/`-max_miniprot` band, so default Step 8 emits it as a miniprot-only
model (unlike `_build_rescue_workspace`, which deliberately puts it out of band so only
the separate rescue pass can recover it).
"""
from __future__ import annotations

import pytest

from tests.test_integration_pipeline import hermetic_pipeline  # noqa: F401


def _wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width)) + "\n"


@pytest.fixture
def step8_workspace(tmp_path):
    """Reference with two coding genes; Liftoff lifts only gene1, so gene2's
    in-band miniprot model is a genuine Step-8 candidate."""
    work = tmp_path / "ws"
    work.mkdir(parents=True, exist_ok=True)

    chrom = ["A"] * 900
    exon1 = "ATG" + "GCT" * 32              # 99 nt
    exon2 = ("GCT" * 32) + "TAA"            # 99 nt
    g2cds = "ATG" + "GCT" * 31 + "TAA"      # 99 nt, a clean ORF
    for i, ch in enumerate(exon1):
        chrom[100 + i] = ch
    for i, ch in enumerate(exon2):
        chrom[300 + i] = ch
    for i, ch in enumerate(g2cds):
        chrom[600 + i] = ch                 # gene2 CDS at 601-699
    seq = "".join(chrom)

    ref_fa = work / "ref.fa"; ref_fa.write_text(">chr1\n" + _wrap(seq))
    tgt_fa = work / "tgt.fa"; tgt_fa.write_text(">chr1\n" + _wrap(seq))

    gene1 = (
        "chr1\t{src}\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\t{src}\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\t{src}\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\t{src}\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\t{src}\tCDS\t301\t399\t.\t+\t0\tID=cds2;Parent=tx1\n"
    )
    gene2 = (
        "chr1\t{src}\tgene\t601\t699\t.\t+\t.\tID=gene2;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t601\t699\t.\t+\t.\tID=tx2;Parent=gene2\n"
        "chr1\t{src}\texon\t601\t699\t.\t+\t.\tID=exon3;Parent=tx2\n"
        "chr1\t{src}\tCDS\t601\t699\t.\t+\t0\tID=cds3;Parent=tx2\n"
    )
    (work / "ref.gff3").write_text(
        "##gff-version 3\n" + gene1.format(src="test") + gene2.format(src="test"))
    # Liftoff MISSES gene2 entirely.
    (work / "liftoff.gff3").write_text(
        "##gff-version 3\n" + gene1.format(src="Liftoff"))
    # MP2 spans exactly the reference gene2 (ratio 1.0 -> inside [0.9, 1.5]).
    (work / "miniprot.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t399\t.\t+\t.\tID=MP1;Target=tx1 1 66\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
        "chr1\tminiprot\tCDS\t301\t399\t.\t+\t0\tID=MP1.cds2;Parent=MP1\n"
        "chr1\tminiprot\tmRNA\t601\t699\t.\t+\t.\tID=MP2;Target=tx2 1 32\n"
        "chr1\tminiprot\tCDS\t601\t699\t.\t+\t0\tID=MP2.cds1;Parent=MP2\n"
    )
    out_dir = work / "out"; out_dir.mkdir()
    return {"ref_fa": ref_fa, "tgt_fa": tgt_fa, "ref_gff": work / "ref.gff3",
            "liftoff": work / "liftoff.gff3", "miniprot": work / "miniprot.gff3",
            "out": out_dir}


def _drive(ws, *, threads: int, suffix: str) -> bytes:
    from lifton import lifton as lifton_main

    out_gff = ws["out"] / f"lifton_{suffix}.gff3"
    argv = [
        str(ws["tgt_fa"]), str(ws["ref_fa"]),
        "-g", str(ws["ref_gff"]),
        "-L", str(ws["liftoff"]), "-M", str(ws["miniprot"]),
        "-o", str(out_gff), "-ad", "RefSeq", "--force",
        "--no-miniprot-rescue",       # isolate Step 8 from the separate rescue pass
    ]
    if threads > 1:
        argv += ["-t", str(threads), "--locus-pipeline"]
    lifton_main.run_all_lifton_steps(lifton_main.parse_args(argv))
    return out_gff.read_bytes()


def test_step8_actually_fires(step8_workspace, hermetic_pipeline):
    """Guard the fixture itself: without a real Step-8 emission this matrix is vacuous."""
    body = _drive(step8_workspace, threads=1, suffix="probe").decode()
    assert "gene2" in body, body
    assert "miniprot" in body


@pytest.mark.parametrize("threads", [2, 4])
def test_parallel_step8_is_byte_identical_to_serial(
        step8_workspace, hermetic_pipeline, threads):
    serial = _drive(step8_workspace, threads=1, suffix="serial")
    parallel = _drive(step8_workspace, threads=threads, suffix=f"t{threads}")
    assert parallel == serial, (
        f"parallel Step 8 at -t{threads} diverged from serial\n"
        f"--- serial ---\n{serial.decode()}\n--- parallel ---\n{parallel.decode()}"
    )


def test_threaded_runs_agree_with_each_other(step8_workspace, hermetic_pipeline):
    assert _drive(step8_workspace, threads=2, suffix="a") == \
           _drive(step8_workspace, threads=4, suffix="b")
