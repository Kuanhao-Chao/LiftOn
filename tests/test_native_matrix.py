"""Phase 10 — extend the byte-identity gate to cover the --native axis.

Phase 9 established the 12-cell 2 × 2 × 3 matrix
(stream × inmemory × threads). Phase 10 adds a fourth axis ``--native
∈ {off, on}`` so the full matrix is 24 cells. Output must remain
byte-identical across all 24 — `--native` changes alignment plumbing,
not algorithms.

Hermeticity: the integration_workspace fixture provides pre-baked
liftoff and miniprot GFF3s via -L / -M, so the actual minimap2 /
miniprot binaries are never invoked. The native-bindings facade is
exercised directly by `tests/test_native_bindings.py`; this matrix
test confirms that toggling the `--native` flag does not perturb the
output shape end-to-end.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.test_integration_pipeline import (  # noqa: F401
    integration_workspace,
    merge_firing_workspace,
    hermetic_pipeline,
)


def _drive(workspace, *, stream: bool, inmem: bool, threads: int,
           native: bool, suffix: str, optimize: bool = False,
           legacy_merge: bool = False) -> bytes:
    from lifton import lifton as lifton_main

    out_gff = workspace["out"] / f"native_{suffix}.gff3"
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
    if native:
        argv.append("--native")
    if optimize:
        argv.append("--optimize")
    if legacy_merge:
        argv.append("--legacy-merge")
    args = lifton_main.parse_args(argv)
    lifton_main.run_all_lifton_steps(args)
    return out_gff.read_bytes()


# ---------------------------------------------------------------------------
# Native flag has zero effect when nothing else is on
# ---------------------------------------------------------------------------

class TestNativeFlagPlain:
    def test_native_only_byte_identical_to_baseline(
            self, integration_workspace, hermetic_pipeline):
        baseline = _drive(integration_workspace,
                          stream=False, inmem=False, threads=1,
                          native=False, suffix="off")
        with_native = _drive(integration_workspace,
                             stream=False, inmem=False, threads=1,
                             native=True, suffix="on")
        assert baseline == with_native


# ---------------------------------------------------------------------------
# Native + threads — parallelism is genuinely unlocked
# ---------------------------------------------------------------------------

class TestNativeUnlocksParallelism:
    @pytest.mark.parametrize("threads", [1, 2, 4])
    def test_native_threads_byte_identical(
            self, integration_workspace, hermetic_pipeline, threads):
        baseline = _drive(integration_workspace,
                          stream=False, inmem=False, threads=1,
                          native=False, suffix=f"native_baseline")
        result = _drive(integration_workspace,
                        stream=True, inmem=True, threads=threads,
                        native=True, suffix=f"native_t{threads}")
        assert result == baseline

    def test_inmemory_parallel_manifest_reports_cloned_fused_materialisation(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        from lifton import lifton_utils

        # The fixture normally supplies -L as a path, which intentionally
        # bypasses Liftoff's in-memory emitter. Feed the equivalent bytes here
        # so this regression exercises the cursor-cloned in-memory gffbase DB.
        liftoff_bytes = integration_workspace["liftoff"].read_bytes()
        monkeypatch.setattr(
            lifton_utils,
            "exec_liftoff",
            lambda _outdir, _ref_db, _args: liftoff_bytes,
        )
        _drive(
            integration_workspace,
            stream=True,
            inmem=True,
            threads=4,
            native=True,
            suffix="inmemory_manifest",
        )
        manifest = json.loads(
            (
                integration_workspace["out"]
                / "lifton_output"
                / "run_manifest.json"
            ).read_text()
        )
        backend = manifest["run"]["backend"]
        assert backend["step7_dispatch"] == (
            "parallel_fused_worker_materialise"
        )
        assert "step7_materialisation_constraint" not in backend
        assert "step7_parallel_fallback" not in backend

    def test_genuine_fast_path_sources_match_file_backed_baseline(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """Exercise real blob/direct-DB inputs, not merely their CLI flags.

        The broad matrix supplies ``-L`` and ``-M`` paths so the external-tool
        short circuit intentionally wins before either fast-path producer can
        run.  Here the mocked producer boundary returns the same fixture bytes
        as a Liftoff blob and a gffbase database as streamed miniprot would.
        This reaches annotation ingestion plus both serial and parallel Step 7
        while keeping the test hermetic and biologically input-identical.
        """
        from lifton import lifton_utils
        from lifton.gffbase import create_db

        stream_db_path = (
            integration_workspace["out"] / "miniprot-stream-fixture.duckdb"
        )
        stream_db = create_db(
            str(integration_workspace["miniprot"]),
            dbfn=str(stream_db_path),
            force=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            build_rtree=False,
        )
        stream_db.conn.execute("CHECKPOINT")
        stream_db.conn.close()

        liftoff_bytes = integration_workspace["liftoff"].read_bytes()
        observed_sources = []

        def liftoff_source(_outdir, _ref_db, args):
            in_memory = bool(args.inmemory_liftoff)
            observed_sources.append(("liftoff", in_memory))
            return (
                liftoff_bytes
                if in_memory
                else str(integration_workspace["liftoff"])
            )

        def miniprot_source(_outdir, args, _target, _proteins):
            streaming = bool(args.stream)
            observed_sources.append(("miniprot", streaming))
            return (
                str(stream_db_path)
                if streaming
                else str(integration_workspace["miniprot"])
            )

        monkeypatch.setattr(lifton_utils, "exec_liftoff", liftoff_source)
        monkeypatch.setattr(lifton_utils, "exec_miniprot", miniprot_source)

        baseline = _drive(
            integration_workspace,
            stream=False,
            inmem=False,
            threads=1,
            native=False,
            suffix="genuine_file_baseline",
        )

        cases = [
            # Liftoff blob: serial semantics at t1 and cursor-cloned fused
            # materialise/process workers at t4, with/without native mode.
            (False, True, 1, False, "inmemory_serial"),
            (False, True, 4, False, "inmemory_parallel"),
            (False, True, 4, True, "inmemory_native_parallel"),
            # Streamed miniprot's production return shape is a direct DuckDB.
            (True, False, 1, False, "stream_serial"),
            (True, False, 4, True, "stream_native_parallel"),
            # Exercise the combined fast path on both dispatch strategies.
            (True, True, 1, True, "combined_native_serial"),
            (True, True, 4, True, "combined_native_parallel"),
        ]
        for stream, inmem, threads, native, suffix in cases:
            result = _drive(
                integration_workspace,
                stream=stream,
                inmem=inmem,
                threads=threads,
                native=native,
                suffix=f"genuine_{suffix}",
            )
            assert result == baseline, suffix

            manifest = json.loads(
                (
                    integration_workspace["out"]
                    / "lifton_output"
                    / "run_manifest.json"
                ).read_text()
            )
            backend = manifest["run"]["backend"]
            counts = manifest["counts"]
            assert backend["liftoff_annotation"] == (
                "gffbase" if inmem else "gffutils"
            )
            assert backend["miniprot_annotation"] == (
                "gffbase" if stream else "gffutils"
            )
            if threads == 1:
                assert backend["step7_dispatch"] == "serial"
                assert backend["step8_dispatch"] == "serial"
                assert backend["step8_overlap_filter"] == "in_process"
            else:
                assert backend["step7_dispatch"] == (
                    "parallel_fused_worker_materialise"
                )
                assert "step7_materialisation_constraint" not in backend
                assert backend["step8_dispatch"] == (
                    "parallel_batched_parent_materialise"
                    if stream
                    else "parallel_scalar_parent_materialise"
                )
                assert backend["step8_overlap_filter"] == (
                    "monotonic_precheck_ordered_recheck"
                )
            assert counts["miniprot_candidates_processed"] == (
                counts["miniprot_candidates_pre_overlap_elided"]
                + counts["miniprot_candidates_submitted_for_analysis"]
            )

        assert set(observed_sources) == {
            ("liftoff", False),
            ("liftoff", True),
            ("miniprot", False),
            ("miniprot", True),
        }


# ---------------------------------------------------------------------------
# Full 24-cell golden gate
# ---------------------------------------------------------------------------

class TestFullNativeMatrix:
    def test_all_24_combinations_byte_identical(
            self, integration_workspace, hermetic_pipeline):
        """Phase 10 contract: every cell of the 2 × 2 × 3 × 2 matrix
        produces identical bytes. ``--native`` must not change the
        output."""
        outputs = {}
        for stream in (False, True):
            for inmem in (False, True):
                for threads in (1, 2, 4):
                    for native in (False, True):
                        key = (f"s{int(stream)}_i{int(inmem)}_"
                               f"t{threads}_n{int(native)}")
                        outputs[key] = _drive(
                            integration_workspace,
                            stream=stream, inmem=inmem,
                            threads=threads, native=native,
                            suffix=key,
                        )
        baseline = outputs["s0_i0_t1_n0"]
        diverged = [k for k, v in outputs.items() if v != baseline]
        assert not diverged, (
            f"Output diverged for: {diverged}. "
            f"Lengths: {[(k, len(outputs[k])) for k in diverged]}"
        )
        assert all(len(v) > 0 for v in outputs.values())


# ---------------------------------------------------------------------------
# Merge-promotion guard — the best-of-outcome merge is now the DEFAULT path
# ---------------------------------------------------------------------------

class TestMergePromotion:
    """Guards the promoted **best-of-outcome** Liftoff/miniprot merge, which is
    now the unconditional default (the pre-promotion unconditional merge is
    reachable only via ``--legacy-merge``; ``--optimize`` is a kept no-op alias).

    This runs on ``merge_firing_workspace`` — a fixture where the lifted CDS
    region carries a lesion so ``liftoff_aln.identity < 1`` AND a pre-baked
    miniprot supplies a higher-identity chunk, so the chaining/merge branch in
    ``run_liftoff.process_liftoff_with_protein`` is actually entered (the
    default ``integration_workspace`` has a perfect ORF → merge never fires →
    such a test would be vacuous). The merge protein identity here is < 1.0, so
    the full per-candidate best-of-outcome compare branch runs (not the
    perfect-protein fast-path skip).

    Pinned invariants (all reliably reproducible synthetically):
      1. the merge **fires** — emitted ``status == LiftOn_chaining_algorithm``
         (the changed code path is exercised; the test is not vacuous);
      2. ``--optimize`` is a true **no-op alias** — its output bytes are
         identical to the default path;
      3. ``--legacy-merge`` (the pre-promotion unconditional merge) runs and
         the promoted default **never lowers** protein identity below it
         (monotonic floor) nor changes the set of emitted transcripts.

    The *strict* revert — default protein identity **>** legacy when the merge
    corrupts a transcript, which is the headline win of best-of-outcome — is
    proven on real divergent data by the benchmark gate
    (``benchmarks/compare``: drosophila 119/4 and mouse_to_rat 294/2
    improved/regressed). It is intentionally NOT asserted here: ORF-rescue
    re-scans the whole spliced exon union, so a single-CDS synthetic merge is
    always equalized back to the Liftoff candidate (a tie), and a multi-exon
    frameshift construction that defeats rescue is too brittle to pin.
    """

    @staticmethod
    def _attr_by_mrna(gff_bytes, attr):
        """Map mRNA ID -> attribute value parsed from a GFF3 blob."""
        out = {}
        for line in gff_bytes.decode().splitlines():
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "mRNA":
                continue
            attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
            if "ID" in attrs and attr in attrs:
                out[attrs["ID"]] = attrs[attr]
        return out

    @classmethod
    def _protein_identity_by_mrna(cls, gff_bytes):
        return {k: float(v) for k, v in cls._attr_by_mrna(gff_bytes, "protein_identity").items()}

    def test_merge_actually_fires(self, merge_firing_workspace, hermetic_pipeline,
                                  monkeypatch):
        # Isolate to the 2-way best-of-outcome merge this class pins. candidate-3
        # (a separate, later, additive 3rd best-of-outcome candidate — a
        # standalone miniprot-only model, default-on, covered by its own
        # tests/test_miniprot_candidate.py) can ALSO win on this fixture once its
        # clean-rebuild fix makes miniprot's model genuinely better here than the
        # merge -- which is correct behaviour for that layer, but would make this
        # test vacuous for what it's actually pinning (the merge branch itself).
        monkeypatch.setenv("LIFTON_MINIPROT_CANDIDATE", "0")
        default = _drive(merge_firing_workspace,
                         stream=False, inmem=False, threads=1,
                         native=False, suffix="merge_default")
        status = self._attr_by_mrna(default, "status")
        assert status, "no mRNA with a status attribute was emitted"
        assert status.get("tx1") == "LiftOn_chaining_algorithm", (
            f"merge branch not exercised; tx1 status = {status.get('tx1')!r}"
        )

    def test_candidate3_fires_on_merge_fixture(self, merge_firing_workspace,
                                               hermetic_pipeline, monkeypatch):
        # Companion to test_merge_actually_fires: with candidate-3 DEFAULT-ON
        # (no env override), miniprot's standalone model is STRICTLY better than
        # the 2-way merge on this fixture, so the emitted tx1 status flips from
        # the merge branch ("LiftOn_chaining_algorithm", pinned above under
        # LIFTON_MINIPROT_CANDIDATE=0) to "LiftOn_miniprot". Pins that the
        # default pipeline actually wires the 3rd candidate end-to-end (the unit
        # coverage is tests/test_miniprot_candidate.py).
        monkeypatch.delenv("LIFTON_MINIPROT_CANDIDATE", raising=False)
        default = _drive(merge_firing_workspace,
                         stream=False, inmem=False, threads=1,
                         native=False, suffix="cand3_default")
        status = self._attr_by_mrna(default, "status")
        assert status.get("tx1") == "LiftOn_miniprot", (
            f"candidate-3 did not fire end-to-end; tx1 status = {status.get('tx1')!r}"
        )

    def test_optimize_is_noop_alias(self, merge_firing_workspace, hermetic_pipeline):
        default = _drive(merge_firing_workspace,
                         stream=False, inmem=False, threads=1,
                         native=False, suffix="alias_default", optimize=False)
        optimized = _drive(merge_firing_workspace,
                           stream=False, inmem=False, threads=1,
                           native=False, suffix="alias_on", optimize=True)
        assert len(default) > 0 and len(optimized) > 0
        assert optimized == default, "--optimize is no longer a byte-for-byte no-op alias"

    def test_default_never_below_legacy_merge(
            self, merge_firing_workspace, hermetic_pipeline):
        default = _drive(merge_firing_workspace,
                         stream=False, inmem=False, threads=1,
                         native=False, suffix="floor_default", legacy_merge=False)
        legacy = _drive(merge_firing_workspace,
                        stream=False, inmem=False, threads=1,
                        native=False, suffix="floor_legacy", legacy_merge=True)
        assert len(default) > 0 and len(legacy) > 0
        d = self._protein_identity_by_mrna(default)
        lg = self._protein_identity_by_mrna(legacy)
        # The merge only re-selects CDS; it never drops or adds transcripts.
        assert set(d) == set(lg), "default vs --legacy-merge emit different transcripts"
        assert d, "no protein_identity attributes emitted"
        # Monotonic floor: the promoted default is never worse than legacy.
        for mrna_id, lg_id in lg.items():
            assert d[mrna_id] >= lg_id - 1e-6, (
                f"promoted default regressed protein identity vs --legacy-merge "
                f"for {mrna_id}: {d[mrna_id]:.6f} < {lg_id:.6f}"
            )
