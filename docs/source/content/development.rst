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
``benchmarks/compare/_runs/<run-id>/``. Set ``LIFTON_PY`` to the interpreter
containing LiftOn's benchmark dependencies, then review a plan before
launching it:

.. code-block:: bash

   LIFTON_PY=/path/to/release-environment/bin/python
   "$LIFTON_PY" -m benchmarks.compare.build_controller start --stage subset-canary --dry-run
   "$LIFTON_PY" -m benchmarks.compare.build_controller start --stage subset-canary
   "$LIFTON_PY" -m benchmarks.compare.build_controller status <run-id>
   "$LIFTON_PY" -m benchmarks.compare.build_controller retry <run-id>
   "$LIFTON_PY" -m benchmarks.compare.build_controller reconcile <run-id> --deep

Promote successful stages deliberately: ``gates``, ``subset-canary``,
``subset``, ``full-canary``, ``full``, then ``e2e-canary``/``e2e`` as needed.
The defaults allocate 8 threads per cell, admit at most 4 active cells and 2
full jobs (32 worker threads total), and launch only while one-minute load is
below 32 and available memory is at least 256 GiB. Launches are staggered by 15
seconds and resource checks repeat every 30 seconds; every value has a
corresponding ``start`` option.

``status`` and ``reconcile`` are automation-safe: exit 0 means every cell and
the controller completed successfully, 1 means failed or blocked, 2 means the
plan/evidence is invalid, and 3 means work is still pending or running.

Release comparisons use the ``paired-*`` variants and require exact, clean
source trees for both arms. Create fresh, persistent clones for the candidate
and reference; development checkouts commonly contain ignored bytecode or
native modules that correctly fail the source-integrity preflight, and
``/tmp`` is unsuitable for multi-day runs. Each repetition alternates
reference/candidate order and runs exclusively so unrelated load cannot
contaminate the profile:

.. code-block:: bash

   "$LIFTON_PY" -m benchmarks.compare.build_controller start \
     --stage paired-subset-canary \
     --candidate-root /persistent/path/lifton-candidate \
     --candidate-sha <40-character-candidate-SHA> \
     --reference-root /persistent/path/lifton-v1.0.10 \
     --reference-sha 3739dfc8f73396fccab7d7be0f95e008179cea5d \
     --paired-repetitions 3 --dry-run

Canaries allow one to five diagnostic repetitions. Non-canary release stages
default to four and require an even count, giving each arm the same number of
first/second positions across the alternating AB/BA order.

For E2E diagnosis, ``--candidate-e2e-mode`` and
``--reference-e2e-mode`` select ``safe``, one of the individual or two-feature
fast-path combinations, or ``fast``. Use the same source SHA in both arms to
isolate mode effects before comparing releases.

A cell receives ``.success`` only after its command, result schema, run
evidence manifest, streaming structural check, and full GFF3 validator succeed.
Paired arms use a source-neutral manifest so releases that predate LiftOn's
native manifest remain auditable. Subset
wall time is normalized by the paired stable executable when both measurements
are available; peak RSS and unpaired cells use absolute comparisons. A
regression above 25% receives one isolated rerun. ``retry`` only resets failed
cells, while ``reconcile --deep`` audits published markers and artifacts. The
controller writes review-only reconciled results and **never modifies the
canonical benchmark baseline**; reseeding remains the
separate, reviewed ``make benchmark-gate-update`` action.

For automation, ``status`` exits 0 only when every cell and the controller are
successful, 1 for failed or blocked runs, 2 for invalid/tampered state, and 3
while work is pending, running, or otherwise incomplete.

Workers run commands in dedicated process groups. Defaults stop subset/gate
cells after 3 hours, full cells after 8 hours, E2E cells after 24 hours, and
low-CPU/no-output stalls after 2 hours; paired limits are doubled because each
cell runs two arms. TERM/KILL cleanup preserves partial artifacts and
attempt-specific logs for diagnosis.

For convenience, ``make benchmark-build-plan`` and ``make benchmark-build``
use ``BENCH_BUILD_STAGE`` (default ``subset-canary``) and optional
``BENCH_BUILD_RUN_ID`` variables.
