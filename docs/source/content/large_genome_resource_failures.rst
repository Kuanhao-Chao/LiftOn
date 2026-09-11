.. _large-genome-resource-failures:

Large-genome resource failures
==============================

This page documents LiftOn's response to `GitHub issue #71
<https://github.com/Kuanhao-Chao/LiftOn/issues/71>`_, in which minimap2 and
miniprot both stopped while mapping an approximately 20-Gb target. It separates
what the logs establish from what still requires evidence from the original
machine.

What ``-11`` means
------------------

Python represents a subprocess terminated by a POSIX signal as a negative
return code. Therefore ``-11`` means that the external process was terminated
by signal 11, ``SIGSEGV``. It is not an ordinary exit status and it does not,
by itself, prove that the operating system killed the process for excessive
memory use. An OOM conclusion requires a kernel message, scheduler reason,
cgroup event, or comparable host evidence.

In issue #71, minimap2 returned ``-11`` during the Liftoff DNA-alignment stage.
The independent miniprot process read 20,029,007,188 target bases in 657
contigs, constructed 156,477,268 blocks, printed ``collected syncmers``, and
then also returned ``-11``. Both failures occurred in native executables, not
in LiftOn's Python annotation-selection code.

The same log reports 101 input-annotation validation errors after experimental
GTF-to-GFF3 conversion. That is a separate correctness concern and should be
reviewed before accepting any eventual annotation, but the available log does
not connect those findings to either native ``SIGSEGV``. Retain
``stats/gff3_input_validation.txt`` and the converted GFF3 when diagnosing the
run.

Why v1.0.12 changes the schedule
--------------------------------

LiftOn v1.0.11 launched Liftoff/minimap2 and miniprot concurrently by default.
That improves wall time on ordinary targets, but both programs construct a
large target index. On a very large genome their high-memory intervals can
overlap. The reporter also requested 40 threads; the Liftoff wrapper created 40
Python workers even though whole-target mapping supplied only one alignment
task, while each native program independently received the full thread count.

Starting in v1.0.12:

* At or below 4,000,000,000 target bases, the default remains concurrent.
* Above 4,000,000,000 bases, Liftoff/minimap2 finishes before miniprot starts.
* ``--serial-aligners`` forces that order at every target size.
* ``--parallel-aligners`` force-enables concurrency. On a target above the
  boundary it prints a warning because the native index peaks can overlap.
* The two overrides are mutually exclusive.
* Liftoff creates at most one worker per real alignment task and divides
  ``--threads`` across those workers. One whole-target task with ``--threads
  40`` uses one worker and passes 40 threads to minimap2.

These are execution-only changes. They do not change LiftOn's mapping
thresholds, candidate scoring, rescue logic, or outcome selection.

The 4-billion-base boundary is a deterministic scheduling policy, not a claim
that every smaller target fits a particular amount of RAM. It matches the
existing point at which LiftOn selects minimap2's split-index path. Actual
memory depends on sequence composition, assembly fragmentation, native-tool
version and options, as well as target length. Users with tighter resource
limits can force serial execution below the boundary; users with measured
headroom can explicitly force concurrency above it.

Controlled reproduction
-----------------------

The following mechanism check was run on 2026-08-25 with LiftOn's bundled
50,818,468-base chromosome-22 fixture, miniprot 0.13-r248, and minimap2
2.28-r1209. Core dumps were disabled. The limits were deliberately selected
below each tool's unconstrained virtual-memory requirement; they are not
recommended production settings.

Unconstrained miniprot completed, reporting 397,020 blocks, 21,359,676
kmer-block pairs, and 0.598 GB peak RSS. With a 900-MiB address-space ceiling,
the same binary stopped immediately after ``collected syncmers`` and Python
observed ``returncode == -11``::

   python - <<'PY'
   import subprocess
   command = [
       "prlimit", "--as=943718400", "--core=0", "--",
       "miniprot", "--gff-only",
       "test/GRCh38_chr22.fa", "test/test_prot.fa",
   ]
   result = subprocess.run(
       command, stdout=subprocess.DEVNULL,
       stderr=subprocess.PIPE, text=True,
   )
   print(result.returncode)
   print(result.stderr)
   PY

Unconstrained minimap2 completed with 0.390 GB reported peak RSS. With a
400-MiB address-space ceiling, the LiftOn-equivalent alignment command returned
``-11``::

   python - <<'PY'
   import subprocess
   command = [
       "prlimit", "--as=419430400", "--core=0", "--",
       "minimap2", "-a", "--end-bonus", "5", "--eqx",
       "-N", "50", "-p", "0.5", "-t", "1",
       "-o", "/dev/null",
       "test/GRCh38_chr22.fa", "test/test.fa",
   ]
   result = subprocess.run(
       command, stdout=subprocess.DEVNULL,
       stderr=subprocess.PIPE, text=True,
   )
   print(result.returncode)
   PY

These controls show that constrained address space can reproduce both signal
signatures, including miniprot's issue-specific last completed stage. They do
not establish the resource state of the reporter's machine. miniprot's index
code allocates its final kmer-block array after reporting collected syncmers;
the relevant source is `index.c
<https://github.com/lh3/miniprot/blob/v0.18/index.c>`_.

Public 22-Gb scale surrogate
----------------------------

The private issue inputs were not available, so total-genome scaling was
tested with the public `Ptaeda2.0 loblolly-pine assembly
<https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000404065.3/>`_. The frozen
FASTA contained 22,103,635,615 bases in 1,760,464 sequences; its longest
sequence was 1,058,009 bases. The NCBI-compressed file matched MD5
``0034bdef1a866ef20b3d981bcebee6c0`` and the decompressed FASTA had SHA-256
``f81866cd46cb2d0b9517f22c18c7342233778916ee964c8d7873c85eacdb9f0f``.

Each native program was first run alone with 40 threads, its standard LiftOn
options, and a single small query. Process-group RSS was sampled every second;
core dumps were disabled, and a watchdog reserved at least 200 GiB of a 1-TiB
host. All runs completed without watchdog intervention.

.. list-table:: Native target-index scale results
   :header-rows: 1
   :widths: 24 20 18 18 20

   * - Program
     - Version
     - Exit
     - Wall seconds
     - Peak RSS (GiB)
   * - miniprot
     - 0.18-r281
     - 0
     - 830.9
     - 111.18
   * - miniprot
     - 0.13-r248
     - 0
     - 549.2
     - 111.31
   * - minimap2
     - 2.31-r1302
     - 0
     - 324.6
     - 54.95
   * - minimap2
     - 2.28-r1209
     - 0
     - 345.9
     - 52.79

The 0.13/2.28 pair was then started concurrently with identical inputs and
options. Both programs again exited 0, but combined process-group RSS reached
150.13 GiB. Running those same commands sequentially bounds their combined
peak by the larger solo value, 111.31 GiB: 38.82 GiB less, or a 25.9% reduction
relative to the concurrent peak. The concurrent trace and tool-specific logs
showed both index builders resident at the peak.

This is mechanism evidence, not an exact reproduction. The surrogate is much
more fragmented than the issue target (1,760,464 versus 657 sequences), and a
single query intentionally isolates target indexing rather than the reporter's
74,405-protein workload or later annotation processing. Its successful exits
on a 1-TiB host do not explain the reporter's machine. Wall times are
descriptive single runs on a shared host, not order-balanced performance
comparisons. The defensible conclusion is narrower: at this biological scale,
each native index alone requires substantial memory, and concurrent execution
can materially raise the peak.

Diagnosing a real run
---------------------

Every v1.0.12 run writes ``lifton_output/run_manifest.json``. For each native
execution it records the exact command, last detected stage, status, return
code, signal number/name, and at most the final 64 KiB of stderr. It also
records target sequence count, total bases, longest sequence, and why the
aligners ran concurrently or sequentially.

When a large run fails, retain that manifest and collect:

* ``lifton -V``, ``miniprot --version``, and ``minimap2 --version``;
* requested RAM and CPUs, process-tree peak RSS, scheduler accounting, OOM or
  cgroup events, and host memory/swap state;
* free space on the filesystem holding ``lifton_output/intermediate_files``;
* ``stats/gff3_input_validation.txt`` and the converted GFF3 when the input was
  GTF, with sensitive attributes redacted if necessary;
* SHA-256 hashes and FASTA sequence statistics for shareable inputs;
* a minimal reproducer, or the original files under an agreed private transfer
  mechanism.

Immediate recovery
------------------

For v1.0.11, add ``--serial-aligners`` to prevent the two native index peaks
from overlapping. In v1.0.12 this is automatic above 4 billion target bases.
If miniprot still fails when it runs alone, the miniprot index itself does not
fit the available resources; LiftOn cannot correct that by scheduling. Increase
the memory allocation, use a more suitable host, or test miniprot resource
options independently and validate any sensitivity change before adopting it.

For a sequence at least 2^31 bases long, use miniprot 0.14 or newer. LiftOn
rejects a parseable older version for that input; an unparseable version emits
a warning. The `upstream correction
<https://github.com/lh3/miniprot/commit/ec4fdba0f01485e3fce9874a9ac1c3e40d409543>`_
changed FASTA-length return values and allocation rounding from 32-bit to
64-bit forms. For minimap2 targets above 4 billion bases, LiftOn retains
minimap2's split-index path. The `minimap2 FAQ
<https://github.com/lh3/minimap2/blob/master/FAQ.md>`_ explains its temporary
storage and memory tradeoffs. Do not force a one-part index or change
miniprot's indexing parameters without measuring resource use and confirming
that the generated annotation remains suitable for the intended analysis.
