
|


.. _installation:

Installation
===============

.. _sys-reqs:

System requirements
-------------------

.. admonition:: Software dependency

   * python >= 3.10      (raised in v1.0.9; 3.9 is EOL and the networkx>=3.3 dependency requires >=3.10)
   * numpy >= 1.22.0
   * gffutils >= 0.10.1
   * biopython>=1.76
   * parasail>=1.2.4
   * intervaltree>=3.1.0
   * networkx>=3.3
   * interlap>=0.2.6
   * miniprot >= 0.10.0   (external binary — see below)
   * minimap2 >= 2.17     (external binary — see below)
   * pyfaidx>=0.5.8
   * pysam>=0.19.1
   * ujson>=3.2.0
   * duckdb>=1.0,!=1.5.3,!=1.5.4
   * pyarrow>=14
   * mappy   (optional; install ``lifton[native]`` or a prebuilt Conda package)

Pip resolves LiftOn's **Python dependencies**. It does **not** install the external
**minimap2** and **miniprot** executables. Fresh standard lifts require both on
``PATH``: minimap2 supplies Liftoff's DNA alignment and miniprot supplies protein
alignment. LiftOn checks the tools it will use before starting a run. Evaluation
and valid precomputed ``-L``/``-M`` inputs do not require the corresponding aligner.
Install the executables through Conda as shown below, or follow their upstream
`minimap2 <https://github.com/lh3/minimap2#install>`_ and
`miniprot <https://github.com/lh3/miniprot#install>`_ installation guides.

.. admonition:: Version warning
   :class: important

   LiftOn requires **Python >= 3.10** (v1.0.9 raised the floor from 3.6 to 3.10; 3.9 is EOL and the ``networkx>=3.3`` dependency requires >=3.10).

   DuckDB 1.5.3 and 1.5.4 have an upstream ``GEOMETRY`` append bug that can
   affect large ``--stream`` miniprot results. Pip will avoid those releases.
   In an existing environment, use DuckDB 1.5.2 or set
   ``LIFTON_DISABLE_RTREE=1``; LiftOn's results are unchanged because region
   queries fall back to the standard B-tree index.

   A target containing any sequence at least 2^31 bases long requires
   **miniprot >= 0.14**. LiftOn v1.0.12 fails preflight for a parseable older
   version rather than risk the known long-sequence limitation. Other
   targets retain the general miniprot >= 0.10 minimum.

   Check out the scientific python ecosystem coordination guideline `SPEC 0 <https://scientific-python.org/specs/spec-0000/>`_ — Minimum Supported Versions to configure the package version compatibility.


.. admonition:: Compiled Python dependencies and macOS
   :class: note

   Dependencies such as ``parasail``, ``pysam``, ``duckdb`` and ``pyarrow`` contain
   compiled code. Wheels are available on common Linux platforms. When a wheel
   is unavailable, particularly on macOS / Apple Silicon, use **prebuilt Conda
   packages** rather than relying on a local source build:

   .. code-block:: bash

      $ conda create -n lifton -y --override-channels -c conda-forge -c bioconda \
            --strict-channel-priority python=3.11 pip numpy biopython \
            parasail-python pysam pyfaidx gffutils intervaltree interlap \
            networkx ujson 'python-duckdb>=1.0,!=1.5.3,!=1.5.4' 'pyarrow>=14' \
            minimap2 miniprot
      $ conda activate lifton
      $ python -m pip install lifton

   The vendored ``gffbase`` backend runs **pure-Python by default** (no
   pre-built ``.so`` ships in the package), so no Rust toolchain is required to
   install or run LiftOn.

|


There are three ways that you can install LiftOn:

.. _install-through-pip:

Install through pip
-------------------------

Install the Python runtime from `PyPI <https://pypi.org/project/lifton/>`_, then
install the two external aligners. Activate the intended environment first:

.. code-block:: bash
   
   $ python -m pip install lifton
   $ conda install --override-channels -c conda-forge -c bioconda \
         --strict-channel-priority minimap2 miniprot
   $ minimap2 --version
   $ miniprot --version
   $ lifton -V

Optional experimental mappy binding
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Standard runs and ``--native`` compatibility hooks do not require mappy. Only
``--native`` **together with** ``LIFTON_NATIVE_LIFTOFF_ALIGN=1`` activates the
experimental in-process Liftoff alignment path. To supply its optional binding:

.. code-block:: bash

   $ python -m pip install 'lifton[native]'

Mappy often builds from source, requiring a **C compiler and zlib development
headers** (for example, ``build-essential`` and ``zlib1g-dev`` on Debian/Ubuntu).
These build tools are unnecessary for a standard LiftOn install on a platform
with wheels for its compiled dependencies. A prebuilt alternative is:

.. code-block:: bash

   $ conda install --override-channels -c conda-forge -c bioconda \
         --strict-channel-priority mappy

Installing mappy does not install the miniprot executable or activate the
experimental path. If the binding is absent, an explicitly requested path falls
back to subprocess minimap2 with a warning; minimap2 must then be available.

|

.. _install-through-conda: 

Install through conda
-------------------------------

A bioconda recipe for LiftOn has been **submitted and is under review**. Once it
is merged, the command below will install LiftOn together with all of its
dependencies:

.. code-block:: bash

   $ conda create -n lifton --override-channels -c conda-forge -c bioconda \
         --strict-channel-priority python=3.11 lifton
   $ conda activate lifton

The recipe includes prebuilt mappy, minimap2 and miniprot. Check the available
LiftOn version before relying on this command; a submitted recipe is not yet
an installable channel package.

|

.. _install-from-source:

Install from source
-------------------------

You can also install LiftOn from source. Check out the latest version on `GitHub <https://github.com/Kuanhao-Chao/LiftOn>`_
!

.. code-block:: bash

   $ git clone https://github.com/Kuanhao-Chao/LiftOn

   $ cd LiftOn
   $ python -m pip install .

.. _seqera-containers:

Seqera Containers / Wave
-----------------------

A pip-only container installs Python packages but lacks the aligner executables.
For a complete standard environment, select Conda Python, minimap2 and miniprot,
plus pip LiftOn. Pin Python explicitly instead of accepting Wave's latest default.
The equivalent environment specification for the packaging release is:

.. code-block:: yaml

   channels:
     - conda-forge
     - bioconda
   dependencies:
     - python=3.11
     - pip
     - minimap2
     - miniprot
     - pip:
         - lifton==1.0.13

Use the new version after its PyPI publication; older distributions still declare
mappy as mandatory and can require GCC/zlib headers. Once the updated Bioconda
recipe is published, use Conda ``lifton=1.0.13`` instead of the pip subsection.
No compiler package is needed to install its prebuilt dependencies.

Check ``lifton -V``, ``minimap2 --version`` and ``miniprot --version`` in the
container, then execute a representative lift. Help output alone does not confirm
that an annotation container is ready. Retain Wave's build report, environment
lockfile and image digest with workflow results.

|

.. _check-LiftOn-installation:

Check LiftOn installation
-------------------------------------

Run the following command to make sure LiftOn is properly installed:

.. code-block:: bash
   
   $ lifton -h


.. dropdown:: Terminal output
    :animate: fade-in-slide-down
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    .. code-block::


      ====================================================================
      An accurate homology lift-over tool between assemblies
      ====================================================================


         ██╗     ██╗███████╗████████╗ ██████╗ ███╗   ██╗
         ██║     ██║██╔════╝╚══██╔══╝██╔═══██╗████╗  ██║
         ██║     ██║█████╗     ██║   ██║   ██║██╔██╗ ██║
         ██║     ██║██╔══╝     ██║   ██║   ██║██║╚██╗██║
         ███████╗██║██║        ██║   ╚██████╔╝██║ ╚████║
         ╚══════╝╚═╝╚═╝        ╚═╝    ╚═════╝ ╚═╝  ╚═══╝

      v1.0.13

      usage: lifton [-h] [-E] [-EL] [-c] [--no-orf-search] [-o FILE] [-u FILE]
                    [-exclude_partial] [-mm2_options =STR] [-mp_options =STR] [-a A]
                    [-s S] [-min_miniprot MIN_MINIPROT] [-max_miniprot MAX_MINIPROT]
                    [-d D] [-flank F] [-V] [-D] [-t THREADS] [-m PATH] [-f TYPES]
                    [-infer-genes] [-infer_transcripts] [-chroms TXT] [-unplaced TXT]
                    [-copies] [-sc SC] [-overlap O] [-mismatch M] [-gap_open GO]
                    [-gap_extend GE] [-polish] [-cds] [-time] [--validate-output]
                    [--validate-verbose] [--allow-partial-output]
                    [--strict-completeness] [--strict-gff] [--stream]
                    [--inmemory-liftoff] [--locus-pipeline] [--no-locus-pipeline]
                    [--parallel-lift] [--no-parallel-lift] [--step7-max-inflight N]
                    [--step8-max-inflight N] [--evaluation-max-inflight N] [--native]
                    [--serial-aligners | --parallel-aligners] [--optimize]
                    [--legacy-merge] [--full-dp-align] [--fast-align] [--gene-only]
                    [--lift-gene-like] [--no-miniprot-rescue] [--miniprot-rescue]
                    [--miniprot-cross-locus-rescue] [--no-miniprot-candidate]
                    [--miniprot-candidate] [--no-adaptive-rescue-floor]
                    [--adaptive-rescue-floor] [--coverage-rescue-gate]
                    [--no-coverage-rescue-gate] [-dir PATH] [--orf-stop-completion]
                    [--no-orf-stop-completion] [--rescue-isoforms]
                    [--no-rescue-isoforms] -g GFF [-P FASTA] [-T FASTA] [-L gff]
                    [-M gff]
                    [--merge-strategy {create_unique,merge,error,warning,replace}]
                    [--id-spec ID_SPEC] [--force] [--verbose] [-ad SOURCE]
                    [--no-auto-convert-gtf]
                    target reference

      Lift features from one genome assembly to another.

      Run `lifton -h` for the complete option list. The full, current flag
      reference -- every option's default, which flags CHANGE the output vs. the
      byte-identical fast-paths, and the kept no-op aliases -- is documented in
      the User Manual / Function manual page. The most-used v1.0.12 options:

        Output-changing defaults (each ships with an opt-out flag):
          (default) lift all gene-like types ......... --gene-only
          (default) miniprot-only rescue ............. --no-miniprot-rescue
          (default) adaptive rescue floor ............ --no-adaptive-rescue-floor
          (default) miniprot merge candidate ......... --no-miniprot-candidate
          (default) best-of-outcome merge ............ --legacy-merge
          (default) banded / windowed alignment ...... --full-dp-align
          (default) protein-coverage rescue gate ..... --no-coverage-rescue-gate
          (default) isoform-aware rescue ............. --no-rescue-isoforms
          (default) terminal-stop completion ......... --no-orf-stop-completion

        Byte-identical fast-paths (output unchanged; pinned by the 24-cell matrix):
          --threads N (per-locus fan-out and the parallel Liftoff lift loop
          are automatic when N > 1; --no-locus-pipeline / --no-parallel-lift
          opt out), --stream, --inmemory-liftoff, --native
          Large-target schedule: automatic sequential execution above 4 billion
          bases; --serial-aligners / --parallel-aligners force either policy
          Memory bounds: --step7-max-inflight / --step8-max-inflight /
          --evaluation-max-inflight (default 2 x --threads)

        Run layout:
          -dir/--intermediate-dir PATH gives a run its own directory for
          intermediate files, statistics, score table and manifest (default:
          lifton_output/ beside the output file)

        Validation:
          --strict-gff (reference, input side),
          --validate-output / --validate-verbose (emitted GFF3)
          --strict-completeness (refuse to publish if any locus was skipped),
          --allow-partial-output (publish after a blocking failure)

        Core mapping thresholds (unchanged):
          -a 0.5 (coverage), -s 0.5 (sequence identity), -overlap 0.1,
          -d 2.0 (distance scaling), -flank 0.0,
          -mm2_options "-a --end-bonus 5 --eqx -N 50 -p 0.5"
|

.. _installation-complete:

Now, you are ready to go !
--------------------------
Please continue to the :ref:`Quick Start Guide`.



|
|
|
|
|


.. image:: ../_images/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image only-light
   :align: center

.. image:: ../_images/jhu-logo-white.png
   :alt: My Logo
   :class: logo, header-image only-dark
   :align: center
