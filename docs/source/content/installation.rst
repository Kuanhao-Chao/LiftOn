
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
   * cigar >= 0.1.3
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
   * mappy   (installed by pip; used only by the explicitly enabled native Liftoff path)

These dependencies are resolved automatically when you ``pip install lifton`` (a bioconda recipe has been submitted and is under review). On macOS / Apple Silicon, install the compiled dependencies via conda first (see the note below), then ``pip install lifton``. LiftOn declares **mappy** as a runtime dependency so supported pip installs are ready for the ``--native`` plus ``LIFTON_NATIVE_LIFTOFF_ALIGN=1`` path; ordinary runs do not activate that path. The two external binaries **miniprot** and **minimap2** are not on PyPI and must be installed manually and be on your ``PATH`` (LiftOn preflight-checks both at startup and exits with a clear message if either is missing). Please see the `miniprot installation guide <https://github.com/lh3/miniprot?tab=readme-ov-file#install>`_ and the `minimap2 installation guide <https://github.com/lh3/minimap2#install>`_ on GitHub. (miniprot drives the protein-to-genome alignment; minimap2 drives the Liftoff DNA lift-over.)

.. admonition:: Version warning
   :class: important

   LiftOn requires **Python >= 3.10** (v1.0.9 raised the floor from 3.6 to 3.10; 3.9 is EOL and the ``networkx>=3.3`` dependency requires >=3.10).

   DuckDB 1.5.3 and 1.5.4 have an upstream ``GEOMETRY`` append bug that can
   affect large ``--stream`` miniprot results. Pip will avoid those releases.
   In an existing environment, use DuckDB 1.5.2 or set
   ``LIFTON_DISABLE_RTREE=1``; LiftOn's results are unchanged because region
   queries fall back to the standard B-tree index.

   If your numpy version is >= 1.25.0, then it requires Python version >= 3.9.

   Check out the scientific python ecosystem coordination guideline `SPEC 0 <https://scientific-python.org/specs/spec-0000/>`_ — Minimum Supported Versions to configure the package version compatibility.


.. admonition:: Native dependencies — use conda
   :class: note

   Several runtime dependencies ship as compiled extensions (``parasail``,
   ``pysam``, ``pyfaidx``, ``gffutils``, ``duckdb``, ``pyarrow``). On
   **macOS / Apple Silicon (ARM)**, ``pip install parasail`` fails to build
   from source — install via **conda** (bioconda / conda-forge) instead, which
   ships pre-built wheels:

   .. code-block:: bash

      $ conda create -n lifton -y python=3.11
      $ conda activate lifton
      $ conda install -y -c bioconda -c conda-forge \
            parasail-python pysam pyfaidx gffutils intervaltree \
            biopython networkx ujson cigar duckdb pyarrow
      $ pip install mappy     # preinstall the declared native-path dependency
      $ pip install lifton

   ``mappy`` is declared in LiftOn's package metadata, but it is used only when
   ``--native`` and ``LIFTON_NATIVE_LIFTOFF_ALIGN=1`` explicitly select the
   experimental Liftoff alignment path. Otherwise LiftOn keeps the proven
   subprocess path. If a manually managed environment lacks mappy, the opt-in
   path falls back gracefully.

   The vendored ``gffbase`` backend runs **pure-Python by default** (no
   pre-built ``.so`` ships in the package), so no Rust toolchain is required to
   install or run LiftOn.

|


There are three ways that you can install LiftOn:

.. _install-through-pip:

Install through pip
-------------------------

LiftOn is on `PyPi <https://pypi.org/project/lifton/>`_ now. Check out all the releases `here <https://pypi.org/manage/project/lifton/releases/>`_. Pip automatically resolves and installs any dependencies required by LiftOn.

.. code-block:: bash
   
   $ pip install LiftOn

|

.. _install-through-conda: 

Install through conda
-------------------------------

A bioconda recipe for LiftOn has been **submitted and is under review**. Once it
is merged, the command below will install LiftOn together with all of its
dependencies:

.. code-block:: bash

   $ conda install -c bioconda lifton

|

.. _install-from-source:

Install from source
-------------------------

You can also install LiftOn from source. Check out the latest version on `GitHub <https://github.com/Kuanhao-Chao/LiftOn>`_
!

.. code-block:: bash

   $ git clone https://github.com/Kuanhao-Chao/LiftOn

   $ python setup.py install

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

      v1.0.10

      usage: lifton [-h] [-E] [-EL] [-c] [--no-orf-search] [-o FILE] [-u FILE]
                    [-exclude_partial] [-mm2_options =STR] [-mp_options =STR] [-a A]
                    [-s S] [-min_miniprot MIN_MINIPROT] [-max_miniprot MAX_MINIPROT]
                    [-d D] [-flank F] [-V] [-D] [-t THREADS] [-m PATH] [-f TYPES]
                    [-infer-genes] [-infer_transcripts] [-chroms TXT] [-unplaced TXT]
                    [-copies] [-sc SC] [-overlap O] [-mismatch M] [-gap_open GO]
                    [-gap_extend GE] [-polish] [-cds] [-time] [--validate-output]
                    [--validate-verbose] [--allow-partial-output] [--strict-completeness]
                    [--strict-gff] [--stream] [--inmemory-liftoff] [--locus-pipeline]
                    [--step7-max-inflight N] [--step8-max-inflight N]
                    [--evaluation-max-inflight N] [--native] [--serial-aligners]
                    [--parallel-aligners] [--optimize] [--legacy-merge] [--full-dp-align]
                    [--fast-align] [--gene-only] [--lift-gene-like] [--no-miniprot-rescue]
                    [--miniprot-rescue] [--miniprot-cross-locus-rescue]
                    [--no-miniprot-candidate] [--miniprot-candidate]
                    [--no-adaptive-rescue-floor] [--adaptive-rescue-floor] -g GFF
                    [-P FASTA] [-T FASTA] [-L gff] [-M gff]
                    [--merge-strategy {create_unique,merge,error,warning,replace}]
                    [--id-spec ID_SPEC] [--force] [--verbose] [-ad SOURCE]
                    [--no-auto-convert-gtf]
                    target reference

      Lift features from one genome assembly to another.

      Run `lifton -h` for the complete option list. The full, current flag
      reference -- every option's default, which flags CHANGE the output vs. the
      byte-identical fast-paths, and the kept no-op aliases -- is documented in
      the User Manual / Function manual page. The most-used v1.0.10 options:

        Output-changing defaults (each ships with an opt-out flag):
          (default) lift all gene-like types ......... --gene-only
          (default) miniprot-only rescue ............. --no-miniprot-rescue
          (default) adaptive rescue floor ............ --no-adaptive-rescue-floor
          (default) miniprot merge candidate ......... --no-miniprot-candidate
          (default) best-of-outcome merge ............ --legacy-merge
          (default) banded / windowed alignment ...... --full-dp-align

        Byte-identical fast-paths (output unchanged; pinned by the 24-cell matrix):
          --threads N --locus-pipeline, --stream, --inmemory-liftoff, --native,
          --serial-aligners
          Memory bounds: --step7-max-inflight / --step8-max-inflight /
          --evaluation-max-inflight (default 2 x --threads)

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
