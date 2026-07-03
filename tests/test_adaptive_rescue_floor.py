"""Tests for the divergence-adaptive miniprot-only rescue PI floor (Phase 2).

The miniprot-only rescue (Iteration 23, default-ON) used a FIXED protein-identity
floor of 0.5. The adaptive floor (PROMOTED to default ON here) lowers that floor
toward 0.30 as DNA-lift gene recall drops, admitting more genuinely-missing genes
on divergent pairs while leaving same/close-species (high-recall) behaviour
unchanged. Mirrors lifton2's shipped design.

Covers the pure mapping (``_adaptive_floor``), the toggle (``_adaptive_floor_on``:
env wins, args fallback, default ON, explicit/env opt-out), the divergence proxy
(``_dna_lift_recall``), and the CLI resolver (``resolve_miniprot_rescue_args``
stashing ``args.adaptive_rescue_floor`` honouring the env). The end-to-end
apply-level behaviour (lower floor -> more rescues on divergent pairs, inert on
high-recall) is proven by benchmarks/compare/adaptive_floor_ab.py.
"""
from types import SimpleNamespace

from lifton import miniprot_rescue as mr


# --------------------------- pure mapping ---------------------------------- #

def test_adaptive_floor_low_recall_returns_min():
    assert mr._adaptive_floor(0.05, 0.5) == 0.30
    assert mr._adaptive_floor(0.0, 0.5) == 0.30


def test_adaptive_floor_high_recall_returns_base():
    assert mr._adaptive_floor(0.6, 0.5) == 0.5
    assert mr._adaptive_floor(0.9, 0.5) == 0.5


def test_adaptive_floor_linear_midpoint():
    # recall 0.3 is halfway between R_LOW(0.1) and R_HIGH(0.5) -> floor halfway
    # between FLOOR_MIN(0.3) and base(0.5) = 0.4.
    assert abs(mr._adaptive_floor(0.3, 0.5) - 0.4) < 1e-9


def test_adaptive_floor_never_above_base():
    # base below FLOOR_MIN: lo = min(0.30, base) so the floor never rises above base.
    assert mr._adaptive_floor(0.0, 0.2) == 0.2
    assert mr._adaptive_floor(0.9, 0.2) == 0.2


def test_adaptive_floor_boundaries_exact():
    assert mr._adaptive_floor(mr.ADAPT_R_LOW, 0.5) == min(mr.ADAPT_FLOOR_MIN, 0.5)
    assert mr._adaptive_floor(mr.ADAPT_R_HIGH, 0.5) == 0.5


# --------------------------- toggle ---------------------------------------- #

def test_adaptive_floor_on_env_wins(monkeypatch):
    monkeypatch.setenv("LIFTON_RESCUE_ADAPTIVE_FLOOR", "1")
    assert mr._adaptive_floor_on(
        SimpleNamespace(adaptive_rescue_floor=False)) is True
    monkeypatch.setenv("LIFTON_RESCUE_ADAPTIVE_FLOOR", "0")
    assert mr._adaptive_floor_on(
        SimpleNamespace(adaptive_rescue_floor=True)) is False


def test_adaptive_floor_on_args_fallback(monkeypatch):
    monkeypatch.delenv("LIFTON_RESCUE_ADAPTIVE_FLOOR", raising=False)
    assert mr._adaptive_floor_on(
        SimpleNamespace(adaptive_rescue_floor=True)) is True
    assert mr._adaptive_floor_on(SimpleNamespace()) is True   # PROMOTED: default ON
    assert mr._adaptive_floor_on(                             # explicit opt-out
        SimpleNamespace(adaptive_rescue_floor=False)) is False


# --------------------------- divergence proxy ------------------------------ #

def test_dna_lift_recall_empty_universe_returns_one():
    # No miniprot-found genes -> rescue is marginal -> base floor (recall 1.0).
    assert mr._dna_lift_recall(set(), set()) == 1.0
    assert mr._dna_lift_recall(set(), {"g1"}) == 1.0


def test_dna_lift_recall_none_lifted_is_zero():
    assert mr._dna_lift_recall({"g1", "g2"}, set()) == 0.0


def test_dna_lift_recall_all_lifted_is_one():
    assert mr._dna_lift_recall({"g1", "g2"}, {"g1", "g2", "g3"}) == 1.0


def test_dna_lift_recall_half():
    assert mr._dna_lift_recall({"g1", "g2"}, {"g1"}) == 0.5


# --------------------------- CLI resolver ---------------------------------- #

class TestResolveAdaptiveFloorArg:
    def test_default_on(self, monkeypatch):
        monkeypatch.delenv("LIFTON_RESCUE_ADAPTIVE_FLOOR", raising=False)
        from lifton import lifton as lifton_main
        args = SimpleNamespace(miniprot_rescue=True, adaptive_rescue_floor=True)
        lifton_main.resolve_miniprot_rescue_args(args)
        assert args.adaptive_rescue_floor is True

    def test_no_flag_opts_out(self, monkeypatch):
        monkeypatch.delenv("LIFTON_RESCUE_ADAPTIVE_FLOOR", raising=False)
        from lifton import lifton as lifton_main
        args = SimpleNamespace(miniprot_rescue=True, adaptive_rescue_floor=False)
        lifton_main.resolve_miniprot_rescue_args(args)
        assert args.adaptive_rescue_floor is False

    def test_env_forces_off(self, monkeypatch):
        monkeypatch.setenv("LIFTON_RESCUE_ADAPTIVE_FLOOR", "0")
        from lifton import lifton as lifton_main
        args = SimpleNamespace(miniprot_rescue=True, adaptive_rescue_floor=True)
        lifton_main.resolve_miniprot_rescue_args(args)
        assert args.adaptive_rescue_floor is False

    def test_env_forces_on(self, monkeypatch):
        monkeypatch.setenv("LIFTON_RESCUE_ADAPTIVE_FLOOR", "1")
        from lifton import lifton as lifton_main
        args = SimpleNamespace(miniprot_rescue=True, adaptive_rescue_floor=False)
        lifton_main.resolve_miniprot_rescue_args(args)
        assert args.adaptive_rescue_floor is True

    def test_missing_attr_defaults_on(self, monkeypatch):
        monkeypatch.delenv("LIFTON_RESCUE_ADAPTIVE_FLOOR", raising=False)
        from lifton import lifton as lifton_main
        args = SimpleNamespace(miniprot_rescue=True)  # no adaptive_rescue_floor
        lifton_main.resolve_miniprot_rescue_args(args)
        assert args.adaptive_rescue_floor is True
