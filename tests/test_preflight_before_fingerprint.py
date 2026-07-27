"""Tier-2 audit fix — a failed preflight must report immediately.

`main()` created the run manifest before checking that minimap2/miniprot are runnable.
Constructing `RunManifest` submits a full SHA-256 of the target genome, reference genome
and reference annotation to a `ThreadPoolExecutor`, and `concurrent.futures` registers
an interpreter-exit hook that JOINS that worker. So a typo'd path or a missing binary
appeared to hang for the duration of a multi-GB hash before printing its error, with no
way to opt out.

The manifest is now created only after the preflight passes; `run_all_lifton_steps`
calls `_ensure_run_manifest` itself, so provenance is unaffected for real runs.
"""
from types import SimpleNamespace

import pytest

from lifton import lifton as lifton_main


def _args(tmp_path, **over):
    base = dict(
        evaluation=False, liftoff=None, miniprot=None, native=False, m=None,
        output=str(tmp_path / "out.gff3"), target="tgt.fa", reference="ref.fa",
        reference_annotation="ref.gff3",
    )
    base.update(over)
    return SimpleNamespace(**base)


def test_missing_miniprot_does_not_fingerprint_first(tmp_path, monkeypatch):
    monkeypatch.setattr(lifton_main, "parse_args", lambda a: _args(tmp_path))
    monkeypatch.setattr(lifton_main.run_miniprot, "check_miniprot_installed",
                        lambda: False)
    monkeypatch.setattr(
        lifton_main, "_ensure_run_manifest",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("input fingerprinting started before the preflight failed")),
    )
    monkeypatch.setattr(lifton_main, "run_all_lifton_steps", lambda a: None)

    with pytest.raises(SystemExit) as exc:
        lifton_main.main([])
    assert "miniprot is not installed" in str(exc.value)


def test_missing_minimap2_does_not_fingerprint_first(tmp_path, monkeypatch):
    monkeypatch.setattr(lifton_main, "parse_args", lambda a: _args(tmp_path))
    monkeypatch.setattr(lifton_main.run_miniprot, "check_miniprot_installed",
                        lambda: True)
    monkeypatch.setattr(lifton_main.run_liftoff, "check_minimap2_installed",
                        lambda *a: False)
    monkeypatch.setattr(
        lifton_main, "_ensure_run_manifest",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("input fingerprinting started before the preflight failed")),
    )
    monkeypatch.setattr(lifton_main, "run_all_lifton_steps", lambda a: None)

    with pytest.raises(SystemExit) as exc:
        lifton_main.main([])
    assert "minimap2 is not installed" in str(exc.value)


def test_manifest_is_still_created_when_the_preflight_passes(tmp_path, monkeypatch):
    """Control: a healthy run still gets its manifest (and its fingerprints)."""
    monkeypatch.setattr(lifton_main, "parse_args", lambda a: _args(tmp_path))
    monkeypatch.setattr(lifton_main.run_miniprot, "check_miniprot_installed",
                        lambda: True)
    monkeypatch.setattr(lifton_main.run_liftoff, "check_minimap2_installed",
                        lambda *a: True)
    created = []
    monkeypatch.setattr(lifton_main, "_ensure_run_manifest",
                        lambda *a, **k: created.append(a))
    monkeypatch.setattr(lifton_main, "run_all_lifton_steps", lambda a: None)

    lifton_main.main([])
    assert created, "manifest must still be created for a healthy run"
