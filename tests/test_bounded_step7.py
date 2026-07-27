"""Bounded-dispatch regressions for Step 7."""

from __future__ import annotations

import io
import threading
import time
from types import SimpleNamespace

import pytest

from lifton import locus_pipeline, parallel
from lifton.locus_pipeline import LocusResult, MaterialisedLocus, StepContext


class _DB:
    def __init__(self, loci):
        self.loci = loci

    def features_of_type(self, _featuretype):
        yield from self.loci


class _Gene:
    ref_gene_id = "ref"

    def __init__(self, index):
        self.index = index

    def write_entry(self, fw, _stats):
        fw.write(f"{self.index}\n")
        return True


def _context(db, *, window, output="out.gff3"):
    return StepContext(
        ref_db={}, l_feature_db=db, m_feature_db=None,
        ref_id_2_m_id_trans_dict={}, tree_dict={}, tgt_fai=None,
        ref_proteins={}, ref_trans={}, ref_features_dict={},
        fw_score=None, fw_chain=None,
        args=SimpleNamespace(
            native=False, output=output, writer_max_pending=0,
            step7_max_inflight=window,
        ),
    )


@pytest.mark.parametrize("mode", ["fused", "two_phase"])
def test_dispatch_window_is_tied_to_ordered_emission(
        monkeypatch, mode):
    """A blocked first locus cannot cause unbounded later submission."""
    count = 30
    window = 4
    loci = [SimpleNamespace(id=f"locus_{i}", index=i)
            for i in range(count)]
    db = _DB(loci)
    ctx = _context(db, window=window)
    gate = threading.Event()
    observed = {}
    submissions = {"process": 0, "prefetch": 0}
    lock = threading.Lock()

    class _Factory:
        def __init__(self, _ctx):
            self.viable = True

    def materialise(index, locus, _ctx_or_factory):
        return MaterialisedLocus(
            submission_index=index, locus=locus, locus_id=locus.id,
        )

    def process(payload, _worker_ctx):
        if payload.submission_index == 0:
            assert gate.wait(timeout=3)
        return LocusResult(
            index=payload.submission_index,
            locus_id=payload.locus_id,
            lifton_gene=_Gene(payload.submission_index),
        )

    original_executor = parallel.ThreadPoolExecutor

    class _CountingExecutor(original_executor):
        def __init__(self, *args, **kwargs):
            prefix = kwargs.get("thread_name_prefix", "")
            self._submission_kind = (
                "prefetch" if prefix == "lifton-prefetch" else "process"
            )
            super().__init__(*args, **kwargs)

        def submit(self, *args, **kwargs):
            with lock:
                submissions[self._submission_kind] += 1
            return super().submit(*args, **kwargs)

    monkeypatch.setattr(locus_pipeline, "_ThreadLocalCtxFactory", _Factory)
    monkeypatch.setattr(locus_pipeline, "materialise_locus", materialise)
    monkeypatch.setattr(
        locus_pipeline, "materialise_locus_with_factory", materialise,
    )
    monkeypatch.setattr(locus_pipeline, "process_locus_native", process)
    monkeypatch.setattr(parallel, "ThreadPoolExecutor", _CountingExecutor)
    if mode == "two_phase":
        monkeypatch.setenv("LIFTON_FUSE_STEP7", "0")

    def release_first():
        time.sleep(0.15)
        with lock:
            observed.update(submissions)
        gate.set()

    releaser = threading.Thread(target=release_first)
    releaser.start()
    output = io.StringIO()
    stats = {"coding": {}, "non-coding": {}, "other": {}}
    processed = parallel.parallel_step7(
        ["gene"], db, ctx, output, stats,
        threads=3, progress_every=0,
    )
    releaser.join(timeout=2)

    assert processed == count
    assert output.getvalue() == "".join(f"{i}\n" for i in range(count))
    assert observed["process"] == window
    if mode == "two_phase":
        # The independent prefetch stage has its own window, but it stops as
        # soon as the processing window is full instead of reading all loci.
        assert observed["prefetch"] <= 2 * window
        assert observed["prefetch"] < count
    else:
        assert observed["prefetch"] == 0


def test_inflight_limit_override_and_default(monkeypatch):
    args = SimpleNamespace()
    assert parallel._step7_inflight_limit(args, 4) == 8
    monkeypatch.setenv("LIFTON_STEP7_MAX_INFLIGHT", "5")
    assert parallel._step7_inflight_limit(args, 4) == 5
    args.step7_max_inflight = 3
    assert parallel._step7_inflight_limit(args, 4) == 3


def test_stdout_mode_progress_uses_stderr(monkeypatch, capsys):
    loci = [SimpleNamespace(id="locus_0")]
    db = _DB(loci)
    ctx = _context(db, window=1, output="stdout")

    monkeypatch.setattr(
        parallel, "process_locus",
        lambda index, _locus, *, ctx: LocusResult(
            index=index, locus_id="locus_0", lifton_gene=_Gene(index),
        ),
    )
    output = io.StringIO()
    parallel.parallel_step7(
        ["gene"], db, ctx, output,
        {"coding": {}, "non-coding": {}, "other": {}},
        threads=1, progress_every=1,
    )

    captured = capsys.readouterr()
    assert captured.out == ""
    assert "LiftOn processed: 1 features" in captured.err
    assert output.getvalue() == "0\n"
