"""Miniprot-only rescue seam (Iteration 22 build; Iteration-23 promotion ON by default).

These tests pin the *seam* of the feature without invoking the full
14-parameter ``process_miniprot`` (which needs a real miniprot+ref FeatureDB):

  * the rescue parses ON by default / ``--no-miniprot-rescue`` opts OUT
    (``--miniprot-rescue`` kept as a no-op alias);
  * ``resolve_miniprot_rescue_args`` honours the env escape hatches and falls
    back safely on garbage;
  * ``_miniprot_rescue_band_ok`` implements the wider sanity band correctly.

Byte-identity of the default (flag-OFF) path is pinned by the 24-cell matrix
(``test_native_matrix.py``) + an end-to-end full-GFF3 ``cmp``; the real-data
rescue *behaviour* is proven by ``benchmarks/compare/miniprot_rescue_ab.py``.
"""
import types
import pytest

from lifton import lifton, run_miniprot


# --------------------------------------------------------------------------- #
# Flag parsing
# --------------------------------------------------------------------------- #
class TestMiniprotRescueFlag:
    def _parse(self, extra):
        return lifton.parse_args(["t.fa", "r.fa", "-g", "r.gff3", "-o", "o.gff3", *extra])

    def test_default_on(self):
        # Iteration-23 promotion: the miniprot-only rescue is now ON by default.
        assert self._parse([]).miniprot_rescue is True

    def test_no_flag_opts_out(self):
        # --no-miniprot-rescue restores the pre-Iteration-23 lift (rescue OFF).
        assert self._parse(["--no-miniprot-rescue"]).miniprot_rescue is False

    def test_flag_on(self):
        # --miniprot-rescue is now a no-op alias; rescue stays ON (the default).
        assert self._parse(["--miniprot-rescue"]).miniprot_rescue is True


# --------------------------------------------------------------------------- #
# Env resolution
# --------------------------------------------------------------------------- #
class TestResolveMiniprotRescueArgs:
    def _args(self, **kw):
        a = types.SimpleNamespace(miniprot_rescue=kw.pop("flag", False))
        return a

    def test_defaults_when_no_env(self, monkeypatch):
        for k in ("LIFTON_MINIPROT_RESCUE", "LIFTON_MINIPROT_RESCUE_MIN_ID", "LIFTON_MINIPROT_RESCUE_LEN"):
            monkeypatch.delenv(k, raising=False)
        a = lifton.resolve_miniprot_rescue_args(self._args())
        assert a.miniprot_rescue is False
        assert a.miniprot_rescue_min_id == 0.5
        assert a.miniprot_rescue_len == (0.5, 2.0)

    def test_env_enables_flag(self, monkeypatch):
        monkeypatch.setenv("LIFTON_MINIPROT_RESCUE", "1")
        a = lifton.resolve_miniprot_rescue_args(self._args(flag=False))
        assert a.miniprot_rescue is True

    def test_explicit_flag_stays_on_without_env(self, monkeypatch):
        monkeypatch.delenv("LIFTON_MINIPROT_RESCUE", raising=False)
        a = lifton.resolve_miniprot_rescue_args(self._args(flag=True))
        assert a.miniprot_rescue is True

    def test_env_floor_and_band(self, monkeypatch):
        monkeypatch.setenv("LIFTON_MINIPROT_RESCUE", "yes")
        monkeypatch.setenv("LIFTON_MINIPROT_RESCUE_MIN_ID", "0.7")
        monkeypatch.setenv("LIFTON_MINIPROT_RESCUE_LEN", "0.3,3.0")
        a = lifton.resolve_miniprot_rescue_args(self._args())
        assert a.miniprot_rescue is True
        assert a.miniprot_rescue_min_id == 0.7
        assert a.miniprot_rescue_len == (0.3, 3.0)

    def test_garbage_env_falls_back(self, monkeypatch):
        monkeypatch.setenv("LIFTON_MINIPROT_RESCUE_MIN_ID", "notafloat")
        monkeypatch.setenv("LIFTON_MINIPROT_RESCUE_LEN", "bogus")
        a = lifton.resolve_miniprot_rescue_args(self._args())
        assert a.miniprot_rescue_min_id == 0.5
        assert a.miniprot_rescue_len == (0.5, 2.0)


# --------------------------------------------------------------------------- #
# Band helper
# --------------------------------------------------------------------------- #
class TestRescueBand:
    def _args(self, band):
        return types.SimpleNamespace(miniprot_rescue_len=band)

    def test_inside_band(self):
        assert run_miniprot._miniprot_rescue_band_ok(0.7, self._args((0.5, 2.0))) is True
        assert run_miniprot._miniprot_rescue_band_ok(1.8, self._args((0.5, 2.0))) is True

    def test_outside_band(self):
        assert run_miniprot._miniprot_rescue_band_ok(0.4, self._args((0.5, 2.0))) is False
        assert run_miniprot._miniprot_rescue_band_ok(2.5, self._args((0.5, 2.0))) is False

    def test_default_band_when_attr_missing(self):
        # An args object without miniprot_rescue_len falls back to (0.5, 2.0).
        assert run_miniprot._miniprot_rescue_band_ok(0.7, types.SimpleNamespace()) is True
        assert run_miniprot._miniprot_rescue_band_ok(0.4, types.SimpleNamespace()) is False

    def test_custom_band_respected(self):
        assert run_miniprot._miniprot_rescue_band_ok(0.4, self._args((0.3, 3.0))) is True
        assert run_miniprot._miniprot_rescue_band_ok(2.5, self._args((0.3, 3.0))) is True
