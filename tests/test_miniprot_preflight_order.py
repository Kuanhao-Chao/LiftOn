"""Tier-2 audit fix — the miniprot preflight must not break the cached -L/-M workflow.

`lifton_utils.exec_miniprot` called `check_miniprot_installed()` (which `sys.exit(1)`s)
BEFORE the `-M` `os.path.exists` short-circuit, so the documented workflow

    lifton -g ref.gff3 -L prior.gff3 -M prior_miniprot.gff3 target.fa reference.fa

exited 1 on a machine without miniprot installed even though the binary is never
invoked. `main()` already gates its own preflight on whether the aligners are needed;
`exec_miniprot` did not — the classic guard-on-one-of-two-paths shape.

`check_miniprot_installed` is also hardened to match its documented twin
`run_liftoff.check_minimap2_installed`: suppress the child's output (under `-o stdout`
the version banner went to the inherited fd 1, i.e. straight into the GFF3 stream) and
require a successful exit so a present-but-broken binary is not reported as installed.
"""
from types import SimpleNamespace

import pytest

from lifton import lifton_utils, run_miniprot


def test_cached_miniprot_gff_skips_the_preflight(tmp_path, monkeypatch):
    """With a valid -M, miniprot is never run, so its absence must not matter."""
    cached = tmp_path / "prior_miniprot.gff3"
    cached.write_text("##gff-version 3\n")

    def _boom():
        raise AssertionError("preflight must not run when -M is reused")

    monkeypatch.setattr(lifton_utils, "check_miniprot_installed", _boom)
    monkeypatch.setattr(
        lifton_utils.run_miniprot, "run_miniprot",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("miniprot must not be executed when -M is reused")),
    )

    args = SimpleNamespace(miniprot=str(cached))
    assert lifton_utils.exec_miniprot(str(tmp_path), args, "tgt.fa", "ref.faa") \
        == str(cached)


def test_missing_miniprot_gff_still_preflights(tmp_path, monkeypatch):
    """Control: when miniprot really has to run, the preflight still guards it."""
    calls = []
    monkeypatch.setattr(lifton_utils, "check_miniprot_installed",
                        lambda: calls.append("checked"))
    monkeypatch.setattr(lifton_utils.run_miniprot, "run_miniprot",
                        lambda *a, **k: "generated.gff3")

    args = SimpleNamespace(miniprot=None)
    assert lifton_utils.exec_miniprot(str(tmp_path), args, "tgt.fa", "ref.faa") \
        == "generated.gff3"
    assert calls == ["checked"]


def test_check_miniprot_installed_suppresses_output_and_checks_returncode(monkeypatch):
    seen = {}

    def _fake_run(command, **kwargs):
        seen.update(kwargs)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(run_miniprot.subprocess, "run", _fake_run)
    assert run_miniprot.check_miniprot_installed() is True
    # The version banner must never reach the inherited fd 1 (it would corrupt
    # `-o stdout` GFF3 output).
    assert seen.get("stdout") is run_miniprot.subprocess.DEVNULL
    assert seen.get("stderr") is run_miniprot.subprocess.DEVNULL


def test_broken_binary_is_not_reported_as_installed(monkeypatch):
    monkeypatch.setattr(run_miniprot.subprocess, "run",
                        lambda *a, **k: SimpleNamespace(returncode=127))
    assert run_miniprot.check_miniprot_installed() is False


@pytest.mark.parametrize("exc", [FileNotFoundError, PermissionError,
                                 NotADirectoryError])
def test_missing_binary_is_not_installed(monkeypatch, exc):
    def _raise(*a, **k):
        raise exc()
    monkeypatch.setattr(run_miniprot.subprocess, "run", _raise)
    assert run_miniprot.check_miniprot_installed() is False
