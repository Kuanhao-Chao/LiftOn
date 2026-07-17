Development and benchmark builds
================================

Run the deterministic correctness gates before changing output, scheduling, or
annotation ingest:

.. code-block:: bash

   pytest tests/ -q
   make test-fast LIFTON_PY=python
   make benchmark-gate LIFTON_PY=python
   make test-fault LIFTON_PY=python

``make test-fault`` exercises injected write, validation, cache, and subprocess
failures. Run the longer repeated cases explicitly with
``make test-fault-stress LIFTON_PY=python``; use tmux for this and other
multi-hour work.

Guarded tmux controller
-----------------------

The controller freezes the Git/tool/input/baseline provenance in an immutable
plan and writes all new artifacts below
``benchmarks/compare/_runs/<run-id>/``. Review a plan before launching it:

.. code-block:: bash

   python -m benchmarks.compare.build_controller start --stage subset-canary --dry-run
   python -m benchmarks.compare.build_controller start --stage subset-canary
   python -m benchmarks.compare.build_controller status <run-id>
   python -m benchmarks.compare.build_controller retry <run-id>
   python -m benchmarks.compare.build_controller reconcile <run-id> --deep

Promote successful stages deliberately: ``gates``, ``subset-canary``,
``subset``, ``full-canary``, ``full``, then ``e2e-canary``/``e2e`` as needed.
The defaults allocate 8 threads per cell, admit at most 4 active cells and 2
full jobs (32 worker threads total), and launch only while one-minute load is
below 32 and available memory is at least 256 GiB. Launches are staggered by 15
seconds and resource checks repeat every 30 seconds; every value has a
corresponding ``start`` option.

A cell receives ``.success`` only after its command, result schema, LiftOn run
manifest, streaming structural check, and full GFF3 validator succeed. Subset
wall time is normalized by the paired stable executable when both measurements
are available; peak RSS and unpaired cells use absolute comparisons. A
regression above 25% receives one isolated rerun. ``retry`` only resets failed
cells, while ``reconcile --deep`` audits published markers and artifacts. The
controller writes review-only reconciled results and **never modifies the
canonical benchmark baseline**; reseeding remains the
separate, reviewed ``make benchmark-gate-update`` action.

For convenience, ``make benchmark-build-plan`` and ``make benchmark-build``
use ``BENCH_BUILD_STAGE`` (default ``subset-canary``) and optional
``BENCH_BUILD_RUN_ID`` variables.
