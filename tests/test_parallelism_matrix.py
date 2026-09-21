"""Phase 9 — full 12-cell parallelism + I/O matrix.

The Phase 5 zero-bug contract says the output GFF3 must be byte-identical
across every supported flag combination. Phase 9 expands that gate to
include a third axis: ``--threads ∈ {1, 2, 4}``.

Cells: ``--stream={off,on}`` × ``--inmemory-liftoff={off,on}`` ×
``--threads ∈ {1,2,4}`` (with the locus-pipeline flag on for any
threads > 1) = **12 cells**, all expected byte-identical.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from tests.test_integration_pipeline import (  # noqa: F401
    integration_workspace,
    hermetic_pipeline,
)


def _drive(workspace, *, stream: bool, inmem: bool, threads: int,
           suffix: str) -> bytes:
    from lifton import lifton as lifton_main

    out_gff = workspace["out"] / f"lifton_{suffix}.gff3"
    argv = [
        str(workspace["tgt_fa"]),
        str(workspace["ref_fa"]),
        "-g", str(workspace["ref_gff"]),
        "-L", str(workspace["liftoff"]),
        "-M", str(workspace["miniprot"]),
        "-o", str(out_gff),
        "-ad", "RefSeq",
        "--force",
    ]
    if stream:
        argv.append("--stream")
    if inmem:
        argv.append("--inmemory-liftoff")
    if threads > 1:
        argv += ["-t", str(threads), "--locus-pipeline"]
    args = lifton_main.parse_args(argv)
    lifton_main.run_all_lifton_steps(args)
    return out_gff.read_bytes()


# ---------------------------------------------------------------------------
# Threads-only axis (with stream + inmemory both off)
# ---------------------------------------------------------------------------

class TestThreadsAxis:
    @pytest.mark.parametrize("threads", [1, 2, 4])
    def test_thread_count_yields_valid_output(self, integration_workspace,
                                              hermetic_pipeline, threads):
        b = _drive(integration_workspace,
                   stream=False, inmem=False, threads=threads,
                   suffix=f"t{threads}")
        assert len(b) > 0

    def test_threads_1_2_4_all_byte_identical(self, integration_workspace,
                                              hermetic_pipeline):
        outputs = {
            t: _drive(integration_workspace,
                      stream=False, inmem=False, threads=t,
                      suffix=f"axis_t{t}")
            for t in (1, 2, 4)
        }
        assert outputs[1] == outputs[2] == outputs[4], (
            "Output diverged across thread counts — determinism gate failed"
        )


# ---------------------------------------------------------------------------
# Determinism with stream + inmemory both on
# ---------------------------------------------------------------------------

class TestThreadsAxisStreamingInmemory:
    @pytest.mark.parametrize("threads", [1, 2, 4])
    def test_under_streaming_path(self, integration_workspace,
                                  hermetic_pipeline, threads):
        b = _drive(integration_workspace,
                   stream=True, inmem=True, threads=threads,
                   suffix=f"si_t{threads}")
        assert len(b) > 0

    def test_streaming_threads_byte_identical(self, integration_workspace,
                                              hermetic_pipeline):
        outputs = {
            t: _drive(integration_workspace,
                      stream=True, inmem=True, threads=t,
                      suffix=f"si_axis_t{t}")
            for t in (1, 2, 4)
        }
        assert outputs[1] == outputs[2] == outputs[4]


# ---------------------------------------------------------------------------
# Full 12-cell matrix
# ---------------------------------------------------------------------------

class TestFull12CellMatrix:
    def test_all_twelve_combinations_byte_identical(
            self, integration_workspace, hermetic_pipeline):
        """The Phase 9 golden output gate: every cell of the
        2 × 2 × 3 matrix produces identical bytes."""
        outputs = {}
        for stream in (False, True):
            for inmem in (False, True):
                for threads in (1, 2, 4):
                    key = f"s{int(stream)}_i{int(inmem)}_t{threads}"
                    outputs[key] = _drive(
                        integration_workspace,
                        stream=stream, inmem=inmem, threads=threads,
                        suffix=key,
                    )
        baseline = outputs["s0_i0_t1"]
        diverged = [k for k, v in outputs.items() if v != baseline]
        assert not diverged, (
            f"Output diverged for: {diverged}. "
            f"Lengths: {[(k, len(outputs[k])) for k in diverged]}"
        )
        # Sanity check: every output non-trivial
        assert all(len(v) > 0 for v in outputs.values())


# ---------------------------------------------------------------------------
# Locus-pipeline path requires threads>1 to fan out
# ---------------------------------------------------------------------------

class TestLocusPipelineSemantics:
    def test_locus_pipeline_with_threads_1_uses_serial(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """When --locus-pipeline is set but --threads is 1 (the
        default), the dispatcher should still take the serial path
        (no ThreadPoolExecutor created)."""
        from lifton import parallel
        original = parallel.ThreadPoolExecutor
        constructed = {"n": 0}

        class TPESpy(original):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        from lifton import lifton as lifton_main
        out_gff = integration_workspace["out"] / "lifton_lp_t1.gff3"
        argv = [
            str(integration_workspace["tgt_fa"]),
            str(integration_workspace["ref_fa"]),
            "-g", str(integration_workspace["ref_gff"]),
            "-L", str(integration_workspace["liftoff"]),
            "-M", str(integration_workspace["miniprot"]),
            "-o", str(out_gff),
            "-ad", "RefSeq", "--force",
            "--locus-pipeline",
        ]
        args = lifton_main.parse_args(argv)
        lifton_main.run_all_lifton_steps(args)
        assert out_gff.exists()
        assert constructed["n"] == 0

    def test_locus_pipeline_with_threads_4_creates_pool(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """With LIFTON_PARALLEL_FORCE + a gffbase backend the dispatcher
        creates at least one ThreadPoolExecutor.

        ``>= 1`` (not ``== 1``) keeps this robust to the dispatch shape.
        Iteration 10 fused the materialise + process phases into ONE pool
        on the on-disk path (so a viable run now builds a single fused pool,
        not a prefetcher pool + a worker pool); the exact-count contracts are
        pinned in ``tests/test_fuse_step7.py`` (fused → 1, ``LIFTON_FUSE_STEP7=0``
        two-phase → 2).
        """
        monkeypatch.setenv("LIFTON_USE_GFFBASE", "1")
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")

        from lifton import parallel
        original = parallel.ThreadPoolExecutor
        constructed = {"n": 0}

        class TPESpy(original):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        from lifton import lifton as lifton_main
        out_gff = integration_workspace["out"] / "lifton_lp_t4.gff3"
        argv = [
            str(integration_workspace["tgt_fa"]),
            str(integration_workspace["ref_fa"]),
            "-g", str(integration_workspace["ref_gff"]),
            "-L", str(integration_workspace["liftoff"]),
            "-M", str(integration_workspace["miniprot"]),
            "-o", str(out_gff),
            "-ad", "RefSeq", "--force",
            "--locus-pipeline", "-t", "4",
        ]
        args = lifton_main.parse_args(argv)
        lifton_main.run_all_lifton_steps(args)
        assert out_gff.exists()
        assert constructed["n"] >= 1

    def test_locus_pipeline_parallel_on_sqlite_by_default(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """Iteration 8: --locus-pipeline + threads>1 on the DEFAULT
        gffutils (SQLite) backend now runs PARALLEL (via the
        materialised-payload + proxy-DB path) WITHOUT --native — it no
        longer silently downgrades to serial. A ThreadPoolExecutor is
        created and the output stays byte-identical to a serial run.
        """
        from lifton import parallel
        original = parallel.ThreadPoolExecutor
        constructed = {"n": 0}

        class TPESpy(original):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        from lifton import lifton as lifton_main
        out_serial = integration_workspace["out"] / "default_serial.gff3"
        out_parallel = integration_workspace["out"] / "default_parallel.gff3"

        base = [
            str(integration_workspace["tgt_fa"]),
            str(integration_workspace["ref_fa"]),
            "-g", str(integration_workspace["ref_gff"]),
            "-L", str(integration_workspace["liftoff"]),
            "-M", str(integration_workspace["miniprot"]),
            "-ad", "RefSeq", "--force",
        ]
        # Serial reference (no pool expected)
        lifton_main.run_all_lifton_steps(
            lifton_main.parse_args(base + ["-o", str(out_serial)]))
        assert constructed["n"] == 0, "serial run unexpectedly created a pool"

        # Parallel on the default gffutils backend — a pool IS created.
        lifton_main.run_all_lifton_steps(lifton_main.parse_args(
            base + ["-o", str(out_parallel), "--locus-pipeline", "-t", "4"]))
        assert constructed["n"] >= 1, (
            "parallel Step 7 did not create a ThreadPoolExecutor on the "
            "default gffutils backend (silently fell back to serial?)"
        )
        # Byte-identical to serial — the whole point of the contract.
        assert out_serial.read_bytes() == out_parallel.read_bytes()

    def test_block_gffutils_restores_serial_fallback(
            self, integration_workspace, hermetic_pipeline, monkeypatch, capsys):
        """The LIFTON_PARALLEL_BLOCK_GFFUTILS opt-out restores the
        pre-Iteration-8 strict serial fallback on gffutils: no pool is
        created, a warning is emitted, and output stays byte-identical.
        """
        monkeypatch.setenv("LIFTON_PARALLEL_BLOCK_GFFUTILS", "1")

        from lifton import parallel
        original = parallel.ThreadPoolExecutor
        constructed = {"n": 0}

        class TPESpy(original):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        from lifton import lifton as lifton_main
        out_serial = integration_workspace["out"] / "block_serial.gff3"
        out_attempted = integration_workspace["out"] / "block_attempted.gff3"

        base = [
            str(integration_workspace["tgt_fa"]),
            str(integration_workspace["ref_fa"]),
            "-g", str(integration_workspace["ref_gff"]),
            "-L", str(integration_workspace["liftoff"]),
            "-M", str(integration_workspace["miniprot"]),
            "-ad", "RefSeq", "--force",
        ]
        lifton_main.run_all_lifton_steps(
            lifton_main.parse_args(base + ["-o", str(out_serial)]))

        lifton_main.run_all_lifton_steps(lifton_main.parse_args(
            base + ["-o", str(out_attempted), "--locus-pipeline", "-t", "4"]))

        # Blocked → fell back to serial → no pool created, byte-identical.
        assert constructed["n"] == 0
        assert out_serial.read_bytes() == out_attempted.read_bytes()
        assert "LIFTON_PARALLEL_BLOCK_GFFUTILS" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# A 3-level hierarchy — the shape every other cell in this file is blind to
# ---------------------------------------------------------------------------

@pytest.fixture
def nested_exon_workspace(tmp_path):
    """A reference with a microRNA precursor, as RefSeq actually encodes one.

    Every other workspace in the suite is gene -> mRNA -> leaf exons, where a
    level-1 exon query and a recursive one are equal by construction. That is
    why a real `-t 1` vs `-t N` divergence sat in the tree while the 24-cell
    matrix, this matrix and `test_fresh_parallel_step7` were all green: the
    serial runtime asked for exons recursively and absorbed a nested
    transcript's exon, and the Step-7 proxy served level-1 only.

    ``test/GRCh38_chr22.gff3`` contains 46 features of this shape and the human
    RefSeq reference 1,915 -- so this fixture is a miniature of ordinary input,
    not a contrived one.
    """
    from tests.test_integration_pipeline import _wrap

    work = tmp_path / "work"
    work.mkdir()

    chrom = ["A"] * 600
    for i, ch in enumerate("ATG" + "GCT" * 32):
        chrom[100 + i] = ch
    for i, ch in enumerate(("GCT" * 32) + "TAA"):
        chrom[300 + i] = ch
    for i, ch in enumerate("GTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCA"):
        chrom[449 + i] = ch
    seq = "".join(chrom)

    ref_fa = work / "ref.fa"
    ref_fa.write_text(">chr1\n" + _wrap(seq))
    tgt_fa = work / "tgt.fa"
    tgt_fa.write_text(">chr1\n" + _wrap(seq))

    # gene2 is the point: the precursor owns exon-pt2-1 and NOTHING else.
    # exon-mir2-1 belongs to the miRNA nested under it.
    hierarchy = (
        "chr1\t{src}\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\t{src}\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\t{src}\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\t{src}\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\t{src}\tCDS\t301\t399\t.\t+\t0\tID=cds2;Parent=tx1\n"
        "chr1\t{src}\tgene\t450\t520\t.\t+\t.\tID=gene2;gene_biotype=miRNA\n"
        "chr1\t{src}\tprimary_transcript\t450\t520\t.\t+\t.\tID=pt2;Parent=gene2\n"
        "chr1\t{src}\texon\t450\t520\t.\t+\t.\tID=exon-pt2-1;Parent=pt2\n"
        "chr1\t{src}\tmiRNA\t460\t481\t.\t+\t.\tID=mir2;Parent=pt2\n"
        "chr1\t{src}\texon\t460\t481\t.\t+\t.\tID=exon-mir2-1;Parent=mir2\n"
    )

    ref_gff = work / "ref.gff3"
    ref_gff.write_text("##gff-version 3\n" + hierarchy.format(src="test"))
    liftoff_gff = work / "liftoff.gff3"
    liftoff_gff.write_text("##gff-version 3\n" + hierarchy.format(src="Liftoff"))

    miniprot_gff = work / "miniprot.gff3"
    miniprot_gff.write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t399\t.\t+\t.\tID=MP1;Target=tx1 1 66\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
        "chr1\tminiprot\tCDS\t301\t399\t.\t+\t0\tID=MP1.cds2;Parent=MP1\n"
    )

    out_dir = work / "out"
    out_dir.mkdir()
    return {
        "ref_fa": ref_fa, "tgt_fa": tgt_fa, "ref_gff": ref_gff,
        "liftoff": liftoff_gff, "miniprot": miniprot_gff,
        "out": out_dir, "work": work,
    }


class TestNestedHierarchyThreadAgreement:
    """`--threads N` must equal `--threads 1` on a 3-level hierarchy too."""

    def test_serial_and_threaded_are_byte_identical(self, nested_exon_workspace,
                                                    hermetic_pipeline):
        serial = _drive(nested_exon_workspace, stream=False, inmem=False,
                        threads=1, suffix="nested_t1")
        threaded = _drive(nested_exon_workspace, stream=False, inmem=False,
                          threads=4, suffix="nested_t4")
        assert serial == threaded, (
            "-t 1 and -t 4 disagree on a nested-exon locus; the serial runtime "
            "and the Step-7 proxy are asking the database different questions"
        )

    def test_the_precursor_emits_only_its_own_exon(self, nested_exon_workspace,
                                                   hermetic_pipeline):
        """Pin the biology, so agreeing on the wrong answer cannot pass.

        The reference gives the precursor one exon (450-520). The miRNA's exon
        (460-481) belongs to the miRNA. Emitting both under the precursor also
        produces two exons where one contains the other -- an overlapping-exon
        pair the reference does not have.
        """
        out = _drive(nested_exon_workspace, stream=False, inmem=False,
                     threads=1, suffix="nested_shape").decode()
        exons = [ln.split("\t") for ln in out.splitlines()
                 if not ln.startswith("#") and ln.split("\t")[2] == "exon"
                 and "Parent=pt2" in ln.split("\t")[8]]
        assert [(e[3], e[4]) for e in exons] == [("450", "520")], (
            f"the precursor absorbed a nested transcript's exon: {exons}"
        )
