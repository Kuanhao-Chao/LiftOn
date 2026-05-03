"""Phase 8 — vendored-Liftoff in-memory refactor.

Three test layers:

1. **Unit tests** for the pure serialiser
   :func:`lifton.liftoff.inmemory_emitter.lifted_features_to_gff3_bytes`
   and the new
   :func:`lifton.liftoff.liftoff_main.run_all_liftoff_steps_inmemory`
   entrypoint surface (without invoking minimap2).

2. **Differential parity tests** asserting that the in-memory bytes
   produced by ``inmemory_emitter`` are byte-identical to the disk
   bytes produced by the legacy ``write_new_gff.write_new_gff`` for
   the same ``lifted_feature_list``.

3. **Integration tests** driving the full
   :func:`lifton.lifton.run_all_lifton_steps` across all four
   ``--stream={off,on}`` × ``--inmemory-liftoff={off,on}``
   combinations and asserting byte-identical output GFF3.
"""

from __future__ import annotations

import io
import os
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

from tests.test_integration_pipeline import (  # noqa: F401
    integration_workspace,
    hermetic_pipeline,
)


# ---------------------------------------------------------------------------
# Unit — emitter & entrypoint surface
# ---------------------------------------------------------------------------

class TestEmitterModuleSurface:
    def test_lifted_features_to_gff3_bytes_imports(self):
        from lifton.liftoff.inmemory_emitter import lifted_features_to_gff3_bytes
        assert callable(lifted_features_to_gff3_bytes)

    def test_run_all_liftoff_steps_inmemory_exists(self):
        from lifton.liftoff.liftoff_main import run_all_liftoff_steps_inmemory
        assert callable(run_all_liftoff_steps_inmemory)

    def test_run_all_liftoff_steps_legacy_still_exists(self):
        from lifton.liftoff.liftoff_main import run_all_liftoff_steps
        assert callable(run_all_liftoff_steps)


class TestEmitterHeader:
    def test_writes_gff_version_directive_for_gff3_dialect(self):
        """Mirrors write_new_gff.write_header for out_type='gff3'."""
        from lifton.liftoff.inmemory_emitter import _write_header_to
        buf = io.StringIO()
        _write_header_to(buf, "gff3", sys_argv=["lifton", "x.fa", "y.fa"])
        text = buf.getvalue()
        assert text.startswith("##gff-version 3\n")
        assert "# LiftOn v" in text
        assert "# lifton x.fa y.fa\n" in text

    def test_no_version_directive_for_gtf_dialect(self):
        from lifton.liftoff.inmemory_emitter import _write_header_to
        buf = io.StringIO()
        _write_header_to(buf, "gtf", sys_argv=["lifton"])
        text = buf.getvalue()
        assert "##gff-version 3" not in text
        assert "# LiftOn v" in text


# ---------------------------------------------------------------------------
# Unit — emitter on a hand-built lifted_feature_list
# ---------------------------------------------------------------------------

def _make_synthetic_lifted_features():
    """Build a minimal ``lifted_feature_list`` shaped exactly like
    Liftoff would produce so we can drive the emitter end-to-end."""
    parent = SimpleNamespace(
        seqid="chr1", source="Liftoff", featuretype="gene",
        start=100, end=200, strand="+", frame=".",
        id="gene1",
        attributes={
            "ID": ["gene1"],
            "copy_id": ["gene1"],
            "copy_num_ID": ["gene1_0"],
            "coverage": ["0.99"],
            "sequence_ID": ["0.98"],
            "extra_copy_number": ["0"],
        },
    )
    child = SimpleNamespace(
        seqid="chr1", source="Liftoff", featuretype="mRNA",
        start=100, end=200, strand="+", frame=".",
        id="tx1",
        attributes={
            "ID": ["tx1"],
            "Parent": ["gene1"],
            "extra_copy_number": ["0"],
        },
    )
    return {parent.attributes["copy_id"][0]: [parent, child]}, parent, child


def _fake_args():
    return SimpleNamespace(a=0.5, s=0.5)


def _fake_feature_db(fmt="gff3"):
    return SimpleNamespace(dialect={"fmt": fmt})


class TestEmitterRoundTrip:
    def test_emits_well_formed_gff3(self):
        from lifton.liftoff.inmemory_emitter import lifted_features_to_gff3_bytes
        lifted, _parent, _child = _make_synthetic_lifted_features()
        blob = lifted_features_to_gff3_bytes(
            lifted, _fake_args(), _fake_feature_db("gff3"),
            sys_argv=["lifton", "demo"],
        )
        assert isinstance(blob, bytes)
        text = blob.decode("utf-8")
        # Header
        assert text.startswith("##gff-version 3\n")
        assert "# LiftOn v" in text
        # Body — gene + mRNA feature lines
        body_lines = [l for l in text.splitlines() if not l.startswith("#") and l]
        feat_types = [l.split("\t")[2] for l in body_lines]
        assert "gene" in feat_types
        assert "mRNA" in feat_types

    def test_emitter_matches_disk_write_byte_for_byte(self, tmp_path):
        """Differential gate: the bytes from the emitter must be the
        same bytes that ``write_new_gff.write_new_gff`` would have
        written to a file."""
        from lifton.liftoff import write_new_gff
        from lifton.liftoff.inmemory_emitter import lifted_features_to_gff3_bytes

        lifted, _, _ = _make_synthetic_lifted_features()
        args = _fake_args()
        feat_db = _fake_feature_db("gff3")

        # Path 1 — disk
        disk_args = SimpleNamespace(**vars(args), output=str(tmp_path / "x.gff3"))
        write_new_gff.write_new_gff(
            _make_synthetic_lifted_features()[0], disk_args, feat_db,
        )
        disk_bytes = (tmp_path / "x.gff3").read_bytes()

        # Path 2 — bytes
        # Both calls share the SAME sys.argv so the provenance comment
        # matches.
        ram_bytes = lifted_features_to_gff3_bytes(
            _make_synthetic_lifted_features()[0], args, feat_db,
            sys_argv=sys.argv,
        )

        assert ram_bytes == disk_bytes, (
            "In-memory emitter diverged from disk write_new_gff."
        )


# ---------------------------------------------------------------------------
# Unit — run_all_liftoff_steps_inmemory tuple shape
# ---------------------------------------------------------------------------

class TestRunAllLiftoffStepsInmemoryShape:
    def test_returns_four_tuple(self, monkeypatch):
        """We don't actually invoke minimap2; we monkey-patch the
        pipeline body so the function under test only runs glue."""
        from lifton.liftoff import liftoff_main

        sentinel_lifted = {"copy1": []}
        sentinel_db = SimpleNamespace(dialect={"fmt": "gff3"})
        sentinel_order = []
        sentinel_unmapped = []

        def fake_pipeline(args, ref_db, *, polish_intermediate_write):
            assert polish_intermediate_write is False, (
                "in-memory entrypoint must not write the polish "
                "intermediate file"
            )
            return (sentinel_lifted, sentinel_db, sentinel_order,
                    sentinel_unmapped)

        monkeypatch.setattr(liftoff_main, "_run_liftoff_pipeline", fake_pipeline)

        result = liftoff_main.run_all_liftoff_steps_inmemory(
            SimpleNamespace(), SimpleNamespace(),
        )
        assert result == (sentinel_lifted, sentinel_db,
                          sentinel_order, sentinel_unmapped)

    def test_legacy_entrypoint_still_writes_disk(self, monkeypatch):
        """Reverse contract: legacy callers must still get the
        post-pipeline write_new_gff call."""
        from lifton.liftoff import liftoff_main

        sentinel_lifted = {"copy1": []}
        sentinel_db = SimpleNamespace(dialect={"fmt": "gff3"})

        def fake_pipeline(args, ref_db, *, polish_intermediate_write):
            return (sentinel_lifted, sentinel_db, [], [])

        write_call = {"n": 0, "args": None}

        def fake_write(lifted, args, feature_db):
            write_call["n"] += 1
            write_call["args"] = (lifted, args, feature_db)

        monkeypatch.setattr(liftoff_main, "_run_liftoff_pipeline", fake_pipeline)
        monkeypatch.setattr(liftoff_main.write_new_gff, "write_new_gff", fake_write)

        liftoff_main.run_all_liftoff_steps(SimpleNamespace(), SimpleNamespace())
        assert write_call["n"] == 1


# ---------------------------------------------------------------------------
# Integration — 2x2 flag matrix end-to-end
# ---------------------------------------------------------------------------

def _drive_pipeline(workspace, *, stream, inmemory_liftoff,
                    suffix):
    """Run the full LiftOn pipeline through the integration fixture
    and return the bytes of the output GFF3."""
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
    if inmemory_liftoff:
        argv.append("--inmemory-liftoff")
    args = lifton_main.parse_args(argv)
    lifton_main.run_all_lifton_steps(args)
    return out_gff.read_bytes()


class TestInmemoryLiftoffMatrix:
    """Phase 8 contract: across all four flag combinations the output
    GFF3 must be byte-identical."""

    def test_off_off(self, integration_workspace, hermetic_pipeline):
        b = _drive_pipeline(integration_workspace, stream=False,
                            inmemory_liftoff=False, suffix="off_off")
        assert len(b) > 0

    def test_off_on(self, integration_workspace, hermetic_pipeline):
        b = _drive_pipeline(integration_workspace, stream=False,
                            inmemory_liftoff=True, suffix="off_on")
        assert len(b) > 0

    def test_on_off(self, integration_workspace, hermetic_pipeline):
        b = _drive_pipeline(integration_workspace, stream=True,
                            inmemory_liftoff=False, suffix="on_off")
        assert len(b) > 0

    def test_on_on(self, integration_workspace, hermetic_pipeline):
        b = _drive_pipeline(integration_workspace, stream=True,
                            inmemory_liftoff=True, suffix="on_on")
        assert len(b) > 0

    def test_all_four_combinations_byte_identical(
            self, integration_workspace, hermetic_pipeline):
        """The Phase 8 golden output gate: every cell in the 2×2
        flag matrix must produce the same bytes."""
        outputs = {}
        for stream in (False, True):
            for inmem in (False, True):
                key = f"s{int(stream)}_i{int(inmem)}"
                outputs[key] = _drive_pipeline(
                    integration_workspace,
                    stream=stream, inmemory_liftoff=inmem,
                    suffix=key,
                )
        baseline = outputs["s0_i0"]
        for key, blob in outputs.items():
            assert blob == baseline, (
                f"Output for {key} diverged from s0_i0 baseline "
                f"({len(blob)} bytes vs {len(baseline)} bytes)"
            )


# ---------------------------------------------------------------------------
# Integration — CLI flag plumbing
# ---------------------------------------------------------------------------

class TestInmemoryLiftoffCLI:
    def test_flag_default_is_false(self):
        from lifton import lifton as lifton_main
        argv = ["t.fa", "r.fa", "-g", "r.gff3"]
        args = lifton_main.parse_args(argv)
        assert args.inmemory_liftoff is False

    def test_flag_set_to_true(self):
        from lifton import lifton as lifton_main
        argv = ["t.fa", "r.fa", "-g", "r.gff3", "--inmemory-liftoff"]
        args = lifton_main.parse_args(argv)
        assert args.inmemory_liftoff is True

    def test_help_text_mentions_flag(self, capsys):
        from lifton import lifton as lifton_main
        with pytest.raises(SystemExit):
            lifton_main.parse_args(["--help"])
        out = capsys.readouterr().out
        assert "--inmemory-liftoff" in out


# ---------------------------------------------------------------------------
# Differential — run_liftoff returns bytes vs path under flag toggle
# ---------------------------------------------------------------------------

class TestRunLiftoffReturnShape:
    """The ``run_liftoff`` wrapper in :mod:`lifton.run_liftoff`
    returns a path string by default and a bytes blob when
    ``--inmemory-liftoff`` is on. Verify via mocked liftoff body."""

    def _mock_liftoff_body(self, monkeypatch):
        """Replace the heavy liftoff pipeline functions with no-op
        sentinels; we only test the wrapper's branching logic."""
        from lifton import run_liftoff
        sentinel_lifted = {"copyA": []}
        sentinel_db = SimpleNamespace(dialect={"fmt": "gff3"})

        def fake_run_all(args, ref_db):
            return None  # legacy path; writes via internal write_new_gff

        def fake_run_inmemory(args, ref_db):
            return (sentinel_lifted, sentinel_db, [], [])

        def fake_write(lifted, args, feature_db):
            # Legacy path uses this to materialise the GFF3 — touch
            # the path so existence checks elsewhere pass if needed.
            with open(args.output, "w") as fh:
                fh.write("##gff-version 3\n")

        def fake_emit(lifted, args, feature_db, *, sys_argv=None):
            return b"##gff-version 3\n# emitted-by-fake\n"

        monkeypatch.setattr(run_liftoff.liftoff_main,
                            "run_all_liftoff_steps", fake_run_all)
        monkeypatch.setattr(run_liftoff.liftoff_main,
                            "run_all_liftoff_steps_inmemory",
                            fake_run_inmemory)
        # Ensure the legacy disk path actually creates a file
        monkeypatch.setattr(
            "lifton.liftoff.write_new_gff.write_new_gff", fake_write,
        )
        monkeypatch.setattr(
            "lifton.liftoff.inmemory_emitter.lifted_features_to_gff3_bytes",
            fake_emit,
        )

    def test_default_returns_path(self, tmp_path, monkeypatch):
        from lifton import run_liftoff as run_liftoff_mod
        self._mock_liftoff_body(monkeypatch)
        args = SimpleNamespace(polish=False, inmemory_liftoff=False,
                               chroms=None, features=None,
                               reference="r", target="t", output="ignored",
                               u="u", a=0.5, s=0.5)
        result = run_liftoff_mod.run_liftoff(str(tmp_path) + "/", None, args)
        assert isinstance(result, str)
        assert result.endswith("liftoff.gff3")

    def test_inmemory_returns_bytes(self, tmp_path, monkeypatch):
        from lifton import run_liftoff as run_liftoff_mod
        self._mock_liftoff_body(monkeypatch)
        args = SimpleNamespace(polish=False, inmemory_liftoff=True,
                               chroms=None, features=None,
                               reference="r", target="t", output="ignored",
                               u="u", a=0.5, s=0.5)
        result = run_liftoff_mod.run_liftoff(str(tmp_path) + "/", None, args)
        assert isinstance(result, (bytes, bytearray))
        assert b"##gff-version 3" in result
