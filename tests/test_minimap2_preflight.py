"""GitHub issue #43 — preflight-check minimap2 (and document it as a requirement).

minimap2 is a hard runtime dependency of the default (vendored Liftoff) DNA-lift,
but — unlike miniprot — it was neither listed in the requirements nor preflight-
checked, so a missing binary surfaced only as a deep FileNotFoundError mid-run.
run_liftoff.check_minimap2_installed mirrors run_miniprot.check_miniprot_installed,
and lifton.main preflight-checks it and exits fast with a clear message.
"""
import pytest
from types import SimpleNamespace

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


def test_main_skips_aligners_for_evaluation(monkeypatch):
    args = SimpleNamespace(evaluation=True, liftoff=None, miniprot=None, m=None)
    monkeypatch.setattr(lifton_main, "parse_args", lambda a: args)
    monkeypatch.setattr(
        lifton_main.run_miniprot, "check_miniprot_installed",
        lambda: (_ for _ in ()).throw(AssertionError("unexpected miniprot check")),
    )
    monkeypatch.setattr(
        lifton_main.run_liftoff, "check_minimap2_installed",
        lambda *a: (_ for _ in ()).throw(AssertionError("unexpected minimap2 check")),
    )
    reached = []
    monkeypatch.setattr(lifton_main, "run_all_lifton_steps", reached.append)

    lifton_main.main([])

    assert reached == [args]


def test_main_skips_aligners_for_valid_precomputed_inputs(tmp_path, monkeypatch):
    liftoff = tmp_path / "liftoff.gff3"
    miniprot = tmp_path / "miniprot.gff3"
    liftoff.write_text("##gff-version 3\n")
    miniprot.write_text("##gff-version 3\n")
    args = SimpleNamespace(
        evaluation=False, liftoff=str(liftoff), miniprot=str(miniprot),
        native=False, m=None,
    )
    monkeypatch.setattr(lifton_main, "parse_args", lambda a: args)
    monkeypatch.setattr(
        lifton_main.run_miniprot, "check_miniprot_installed",
        lambda: (_ for _ in ()).throw(AssertionError("unexpected miniprot check")),
    )
    monkeypatch.setattr(
        lifton_main.run_liftoff, "check_minimap2_installed",
        lambda *a: (_ for _ in ()).throw(AssertionError("unexpected minimap2 check")),
    )
    reached = []
    monkeypatch.setattr(lifton_main, "run_all_lifton_steps", reached.append)

    lifton_main.main([])

    assert reached == [args]


def test_main_checks_custom_minimap2_path(monkeypatch):
    custom = "/opt/custom/minimap2"
    args = SimpleNamespace(
        evaluation=False, liftoff=None, miniprot="cached.gff3",
        native=False, m=custom,
    )
    monkeypatch.setattr(lifton_main.os.path, "exists", lambda p: p == "cached.gff3")
    monkeypatch.setattr(lifton_main, "parse_args", lambda a: args)
    monkeypatch.setattr(lifton_main.run_miniprot, "check_miniprot_installed", lambda: True)
    observed = []
    monkeypatch.setattr(
        lifton_main.run_liftoff, "check_minimap2_installed",
        lambda path: observed.append(path) or True,
    )
    monkeypatch.setattr(lifton_main, "run_all_lifton_steps", lambda a: None)

    lifton_main.main([])

    assert observed == [custom]


def test_el_implies_evaluation(tmp_path):
    args = lifton_main.parse_args([
        "target.fa", "reference.fa", "-g", "reference.gff3", "-EL",
    ])
    assert args.evaluation is True
