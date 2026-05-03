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

from pathlib import Path

import pytest

from tests.test_integration_pipeline import (  # noqa: F401
    integration_workspace,
    hermetic_pipeline,
)


def _drive(workspace, *, stream: bool, inmem: bool, threads: int,
           native: bool, suffix: str) -> bytes:
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
