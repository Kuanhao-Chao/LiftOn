
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
   * miniprot >= 0.10.0
   * pyfaidx>=0.5.8
   * pysam>=0.19.1
   * ujson>=3.2.0
   * duckdb>=1.0
   * pyarrow>=14
   * mappy   (*optional* — only needed for the in-process ``--native`` path)

These dependencies are resolved automatically when you ``pip install lifton`` (a bioconda recipe is planned). On macOS / Apple Silicon, install the compiled dependencies via conda first (see the note below), then ``pip install lifton``. Two exceptions: **mappy** is optional (only the ``--native`` path needs it), and **miniprot**. Since miniprot is not on PyPi, you will need to install it manually. Please check out the `miniprot installation guide <https://github.com/lh3/miniprot?tab=readme-ov-file#install>`_ on `GitHub <https://github.com/lh3/miniprot>`_.

.. admonition:: Version warning
   :class: important

   LiftOn requires **Python >= 3.10** (v1.0.9 raised the floor from 3.6 to 3.10; 3.9 is EOL and the ``networkx>=3.3`` dependency requires >=3.10).

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
      $ pip install mappy     # optional — only for the --native path
      $ pip install lifton

   ``mappy`` is **optional**: it enables the in-process ``--native`` minimap2 /
   miniprot path. If it is not installed, ``--native`` falls back gracefully to
   the subprocess path.

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

Installing LiftOn through bioconda. It will also automatically install all the dependencies required by LiftOn.

.. code-block:: bash
   
   TBC

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

      v1.0.9

      usage: lifton [-h] [-E] [-EL] [-c] [--no-orf-search] [-o FILE] [-u FILE]
                    [-exclude_partial] [-mm2_options =STR] [-mp_options =STR] [-a A]
                    [-s S] [-min_miniprot MIN_MINIPROT] [-max_miniprot MAX_MINIPROT]
                    [-d D] [-flank F] [-V] [-D] [-t THREADS] [-m PATH] [-f TYPES]
                    [-infer-genes] [-infer_transcripts] [-chroms TXT] [-unplaced TXT]
                    [-copies] [-sc SC] [-overlap O] [-mismatch M] [-gap_open GO]
                    [-gap_extend GE] [-polish] [-cds] [-time] [--validate-output]
                    [--validate-verbose] [--strict-gff] [--stream] [--inmemory-liftoff]
                    [--locus-pipeline] [--native] [--serial-aligners] [--legacy-merge]
                    [--full-dp-align] [--gene-only] [--no-miniprot-rescue]
                    [--no-auto-convert-gtf] -g GFF [-P FASTA] [-T FASTA] [-L gff]
                    [-M gff] target reference

      Lift features from one genome assembly to another.

      Run `lifton -h` for the complete option list. The full, current flag
      reference -- every option's default, which flags CHANGE the output vs. the
      byte-identical fast-paths, and the kept no-op aliases -- is documented in
      the User Manual / Function manual page. The most-used v1.0.9 options:

        Output-changing defaults (each ships with an opt-out flag):
          (default) lift all gene-like types ......... --gene-only
          (default) miniprot-only rescue ............. --no-miniprot-rescue
          (default) best-of-outcome merge ............ --legacy-merge
          (default) banded / windowed alignment ...... --full-dp-align

        Byte-identical fast-paths (output unchanged; pinned by the 24-cell matrix):
          --threads N --locus-pipeline, --stream, --inmemory-liftoff, --native,
          --serial-aligners

        Validation:
          --strict-gff (reference, input side),
          --validate-output / --validate-verbose (emitted GFF3)

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