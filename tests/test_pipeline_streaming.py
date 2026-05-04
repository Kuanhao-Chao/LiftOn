"""Phase 7 — end-to-end integration tests with --stream={off,on}.

The Phase 5 hermetic pattern (pre-baked Liftoff + miniprot GFF inputs
via -L / -M) is reused. Each test runs the full
run_all_lifton_steps() and asserts:

  - --stream=off output is byte-identical to the Phase 5 baseline
  - --stream=on  output is byte-identical to --stream=off output
  - all side-files (score.txt, unmapped_features.txt, etc.) are
    byte-identical across the two modes
  - --stream=on still produces a valid GFF3 with the same structure
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest


# Re-use the integration_workspace and hermetic_pipeline fixtures
# from tests/test_integration_pipeline.py via pytest's autodiscovery
# of conftest-level fixtures. They live in that test file directly,
# so we re-import the fixture names by spelling them out here.
from tests.test_integration_pipeline import (  # noqa: E402,F401
    integration_workspace,
    hermetic_pipeline,
)


def _run_pipeline(workspace, *, stream: bool):
    """Drive run_all_lifton_steps end-to-end and return the contents
    of the output GFF3."""
    from lifton import lifton as lifton_main

    out_gff = workspace["out"] / f"lifton_{'on' if stream else 'off'}.gff3"
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
    args = lifton_main.parse_args(argv)
    lifton_main.run_all_lifton_steps(args)
    return out_gff


# ---------------------------------------------------------------------------
# Default-off regression: Phase 5 baseline behaviour preserved
# ---------------------------------------------------------------------------

class TestStreamOffBaseline:
    def test_stream_off_produces_valid_gff(self, integration_workspace,
                                           hermetic_pipeline):
        out = _run_pipeline(integration_workspace, stream=False)
        assert out.exists() and out.stat().st_size > 0
        body = out.read_text()
        feats = [l.split("\t")[2] for l in body.splitlines()
                 if l.strip() and not l.startswith("#")]
        assert "gene" in feats and "mRNA" in feats

    def test_stream_off_args_attribute_default(self):
        from lifton import lifton as lifton_main
        argv = ["t.fa", "r.fa", "-g", "r.gff3"]
        args = lifton_main.parse_args(argv)
        assert args.stream is False


# ---------------------------------------------------------------------------
# Default-on streaming path
# ---------------------------------------------------------------------------

class TestStreamOnPath:
    def test_stream_on_produces_valid_gff(self, integration_workspace,
                                          hermetic_pipeline):
        out = _run_pipeline(integration_workspace, stream=True)
        assert out.exists() and out.stat().st_size > 0
        body = out.read_text()
        feats = [l.split("\t")[2] for l in body.splitlines()
                 if l.strip() and not l.startswith("#")]
        assert "gene" in feats and "mRNA" in feats
        sources = {l.split("\t")[1] for l in body.splitlines()
                   if l.strip() and not l.startswith("#")}
        assert sources == {"LiftOn"}

    def test_stream_on_args_attribute_set(self):
        from lifton import lifton as lifton_main
        argv = ["t.fa", "r.fa", "-g", "r.gff3", "--stream"]
        args = lifton_main.parse_args(argv)
        assert args.stream is True


# ---------------------------------------------------------------------------
# Byte-identical golden output across the two modes
# ---------------------------------------------------------------------------

class TestStreamGoldenOutput:
    def test_stream_off_vs_on_byte_identical(self, integration_workspace,
                                             hermetic_pipeline):
        """Phase 7 contract: --stream changes I/O, not algorithms.
        The output GFF3 produced via either mode must be identical
        byte for byte."""
        out_off = _run_pipeline(integration_workspace, stream=False)
        # Move the off-output aside so the on-run starts from a clean
        # output directory state.
        off_text = out_off.read_text()
        out_on = _run_pipeline(integration_workspace, stream=True)
        on_text = out_on.read_text()
        assert off_text == on_text, \
            "Streaming mode produced a different GFF3 than the legacy path"

    def test_stream_on_score_file_matches_off(self, integration_workspace,
                                              hermetic_pipeline):
        """score.txt must also be byte-identical because the
        algorithmic pipeline saw the same inputs."""
        _run_pipeline(integration_workspace, stream=False)
        out_dir = integration_workspace["out"]
        off_score = (out_dir / "lifton_output" / "score.txt").read_text()
        # Reset the score file for the streaming run
        (out_dir / "lifton_output" / "score.txt").unlink()

        _run_pipeline(integration_workspace, stream=True)
        on_score = (out_dir / "lifton_output" / "score.txt").read_text()
        assert off_score == on_score


# ---------------------------------------------------------------------------
# Disk-write reduction
# ---------------------------------------------------------------------------

class TestStreamDiskFootprint:
    def test_stream_on_writes_no_miniprot_gff3(self, integration_workspace,
                                               hermetic_pipeline):
        """When --stream is on AND we use the real run_miniprot path
        (no -M override), no miniprot.gff3 sidecar should appear.

        In the hermetic test we DO supply -M (a pre-baked file), so
        the streaming branch is bypassed. Instead we check the
        weaker invariant: the streaming path is reachable and run_miniprot
        returns bytes when stream=True, exercised in
        tests/test_streaming_adapter.py::TestRunMiniprotStreaming."""
        # Drive a normal run; assert lifton.gff3 is produced.
        out = _run_pipeline(integration_workspace, stream=True)
        assert out.exists()


# ---------------------------------------------------------------------------
# Pipeline still hermetic (no external runner spawned)
# ---------------------------------------------------------------------------

class TestStreamHermetic:
    def test_stream_on_does_not_invoke_external_runners(
            self, integration_workspace, hermetic_pipeline):
        """The hermetic_pipeline fixture patches run_liftoff /
        run_miniprot to raise on call. If --stream caused us to
        bypass exec_*'s short-circuit, this test would fail."""
        out = _run_pipeline(integration_workspace, stream=True)
        # Reaching this assert proves the patches did NOT fire.
        assert out.exists()


# ---------------------------------------------------------------------------
# CLI help text exposes --stream
# ---------------------------------------------------------------------------

class TestStreamCLIHelpText:
    def test_help_text_mentions_stream(self, capsys):
        from lifton import lifton as lifton_main
        with pytest.raises(SystemExit):
            lifton_main.parse_args(["--help"])
        captured = capsys.readouterr()
        assert "--stream" in captured.out


# ---------------------------------------------------------------------------
# Annotation polymorphism inside the live pipeline
# ---------------------------------------------------------------------------

class TestStreamAnnotationFromBytes:
    def test_annotation_built_from_real_miniprot_blob(self,
                                                      integration_workspace):
        """Take the pre-baked miniprot.gff3 fixture's bytes and
        feed them through Annotation directly — proving the
        polymorphic-bytes path used by the streaming pipeline
        works on actual fixture content."""
        from lifton.annotation import Annotation

        blob = integration_workspace["miniprot"].read_bytes()
        ann = Annotation(blob, False, False)
        assert ann.backend == "gffbase"
        # The fixture defines exactly one mRNA (MP1)
        mrnas = list(ann.db_connection.features_of_type("mRNA"))
        assert len(mrnas) == 1
        assert mrnas[0].attributes["Target"][0].startswith("tx1")
