"""GitHub issue #43 — preflight-check minimap2 (and document it as a requirement).

minimap2 is a hard runtime dependency of the default (vendored Liftoff) DNA-lift,
but — unlike miniprot — it was neither listed in the requirements nor preflight-
checked, so a missing binary surfaced only as a deep FileNotFoundError mid-run.
run_liftoff.check_minimap2_installed mirrors run_miniprot.check_miniprot_installed,
and lifton.main preflight-checks it and exits fast with a clear message.
"""
import pytest

from lifton import run_liftoff
from lifton import lifton as lifton_main


def test_check_minimap2_installed_true(monkeypatch):
    monkeypatch.setattr(run_liftoff.subprocess, "run", lambda *a, **k: None)
    assert run_liftoff.check_minimap2_installed() is True


@pytest.mark.parametrize("exc", [FileNotFoundError, PermissionError,
                                 NotADirectoryError])
def test_check_minimap2_installed_false_when_unrunnable(monkeypatch, exc):
    def _raise(*a, **k):
        raise exc()
    monkeypatch.setattr(run_liftoff.subprocess, "run", _raise)
    assert run_liftoff.check_minimap2_installed() is False


def test_main_exits_fast_when_minimap2_missing(monkeypatch):
    monkeypatch.setattr(lifton_main.run_miniprot, "check_miniprot_installed", lambda: True)
    monkeypatch.setattr(lifton_main.run_liftoff, "check_minimap2_installed", lambda: False)
    monkeypatch.setattr(lifton_main, "parse_args", lambda a: object())
    reached = []
    monkeypatch.setattr(lifton_main, "run_all_lifton_steps", lambda a: reached.append(1))
    with pytest.raises(SystemExit) as exc:
        lifton_main.main([])
    assert "minimap2 is not installed" in str(exc.value)
    assert reached == []          # never reached the pipeline


def test_main_proceeds_when_both_installed(monkeypatch):
    monkeypatch.setattr(lifton_main.run_miniprot, "check_miniprot_installed", lambda: True)
    monkeypatch.setattr(lifton_main.run_liftoff, "check_minimap2_installed", lambda: True)
    monkeypatch.setattr(lifton_main, "parse_args", lambda a: object())
    reached = []
    monkeypatch.setattr(lifton_main, "run_all_lifton_steps", lambda a: reached.append(1))
    lifton_main.main([])
    assert reached == [1]         # both present -> pipeline runs
