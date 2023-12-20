.. _quick-start:

Quick Start Guide
=================

This page provides simple quick-start information for using LiftOn to lift-over annotations from chr21 GRCh38 to T2T-CHM13.

.. with :code:`BAM` and :code:`GFF` files. Please read the :ref:`alignment-detailed-section` or :ref:`annotation-detailed-section` page for more details on each step.

If you haven't already, please follow the steps in the :ref:`Installation` page to install and load LiftOn.

|

.. _super-quick-start:

Super-Quick Start (3 lines of code)
+++++++++++++++++++++++++++++++++++

There are two main use case scenarios of LiftOn. The first one is :ref:`running with an alignment file <LiftOn-bam-quick>` and second one is :ref:`running with an annotation file or assembled transcripts <LiftOn-gff-quick>`. Both of them can be done in three lines of code. 

Before you get started, make sure you have already cloned the :ref:`LiftOn GitHub repository <install-from-source>`. We provide a few examples below:

|

.. _LiftOn-bam-quick:

Example 1: clean up alignment files  (:code:`BAM`)
-----------------------------------------------------

.. code-block:: bash

    $ cd test

    # Step 1: extract splice junctions in the alignment file
    $ lifton -g GRCh38_chr22.gff3 -dir GRCh38_2_CHM13 -o GRCh38_2_CHM13_lifton.gff3 -copies chm13_chr22.fa GRCh38_chr22.fa


| 

.. _LiftOn-gff-quick:

Example 2: evaluate annotation files / assembled transcripts (:code:`GFF`)
-----------------------------------------------------------------------------

.. code-block:: bash

    $ cd test

    # Step 1: extract introns in the annotation
    $ LiftOn extract refseq_40_GRCh38.p14_chr_fixed.gff -o tmp_out_annotation

    # Step 2: score introns in the annotation
    $ LiftOn score -G chr9_subset.fa -m ../model/LiftOn_script.pt -o tmp_out_annotation tmp_out_annotation/junction.bed

    #Step 3: output statistics of each transcript
    $ LiftOn clean -o tmp_out_annotation -t 0.8

| 

LiftOn can also :ref:`run on non-human species <generalization-introduction>`. 

.. _LiftOn-generalization-example:

Example of evaluating mouse annotation files (:code:`GFF`)
----------------------------------------------------------------------

.. code-block:: bash

    $ cd test

    # Step 1: extract introns in the annotation
    LiftOn extract mouse_chr19_subset.gff -o tmp_out_generalization

    # Step 2: score introns in the annotation
    LiftOn score -A GRCm39_assembly_report.txt -G mouse_chr19.fa -m ../model/LiftOn_script.pt -o tmp_out_generalization tmp_out_generalization/junction.bed

    #Step 3: output statistics of each transcript
    LiftOn clean -o tmp_out_generalization -t 0.8

|

.. _google-colab:

Try LiftOn on Google Colab
+++++++++++++++++++++++++++++++++++

We created some reproducible and easy-to-run LiftOn examples on Google Colab. It's a good starting point, so go ahead and check them out!


.. image:: https://colab.research.google.com/assets/colab-badge.svg
    :target: https://colab.research.google.com/github/Kuanhao-Chao/LiftOn/blob/main/notebook/LiftOn_example.ipynb


|

For more detailed analysis steps, please check :

.. seealso::
    
    * :ref:`alignment-detailed-section`

    * :ref:`annotation-detailed-section`


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