"""Intermediate Liftoff/miniprot outputs live under ``lifton_output/``.

v1.0.10 and v1.0.11 wrote them to ``<out>/lifton_outputliftoff/`` and
``<out>/lifton_outputminiprot/``: the caller passes ``<out>/lifton_output``
without a trailing separator and the runners concatenated ``"liftoff/"`` and
``"miniprot/"`` onto it. The external tools are stubbed; only paths are
checked.
"""
from __future__ import annotations

import os
import types

from lifton import lifton, run_liftoff, run_miniprot


def test_liftoff_writes_under_lifton_output(tmp_path, monkeypatch):
    seen = {}

    def fake_liftoff(liftoff_args, ref_db):
        seen["output"], seen["unmapped"] = liftoff_args.output, liftoff_args.u
        with open(liftoff_args.output, "w") as handle:
            handle.write("##gff-version 3\n")

    monkeypatch.setattr(run_liftoff.liftoff_main, "run_all_liftoff_steps", fake_liftoff)
    _, lifton_outdir = lifton._run_outdirs(str(tmp_path / "out.gff3"))
    assert not lifton_outdir.endswith(os.sep)   # the shape that exposed the bug
    args = types.SimpleNamespace(inmemory_liftoff=False, polish=False)
    result = run_liftoff.run_liftoff(lifton_outdir, None, args)
    expected_dir = tmp_path / "lifton_output" / "liftoff"
    assert result == str(expected_dir / "liftoff.gff3")
    assert seen["output"] == result
    assert seen["unmapped"] == str(expected_dir / "unmapped_features.txt")
    assert not (tmp_path / "lifton_outputliftoff").exists()


def test_miniprot_writes_under_lifton_output(tmp_path, monkeypatch):
    def fake_run(command, stdout=None, **kwargs):
        stdout.write("##gff-version 3\nchr1\tminiprot\tmRNA\t1\t9\t.\t+\t.\tID=MP1\n")
        return types.SimpleNamespace(returncode=0, stderr_tail="",
                                     stderr_error_seen=False)

    monkeypatch.setattr(run_miniprot, "run_with_bounded_stderr", fake_run)
    monkeypatch.setattr(run_miniprot, "record_execution", lambda *a, **k: None)
    _, lifton_outdir = lifton._run_outdirs(str(tmp_path / "out.gff3"))
    args = types.SimpleNamespace(mp_options="", threads=1, stream=False)
    result = run_miniprot.run_miniprot(lifton_outdir, args, "tgt.fa", "prot.fa")
    assert result == str(tmp_path / "lifton_output" / "miniprot" / "miniprot.gff3")
    assert os.path.getsize(result) > 0
    assert not (tmp_path / "lifton_outputminiprot").exists()
