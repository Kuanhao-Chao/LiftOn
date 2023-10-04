.. _quick-start:

Quick Start Guide
=================

This page provides simple quick-start information for using Splam with :code:`BAM` and :code:`GFF` files. Please read the :ref:`alignment-detailed-section` or :ref:`annotation-detailed-section` page for more details on each step.

If you haven't already, please follow the steps in the :ref:`Installation` page to install and load Splam.

|

.. _super-quick-start:

Super-Quick Start (3 lines of code)
+++++++++++++++++++++++++++++++++++

There are two main use case scenarios of Splam. The first one is :ref:`running with an alignment file <splam-bam-quick>` and second one is :ref:`running with an annotation file or assembled transcripts <splam-gff-quick>`. Both of them can be done in three lines of code. 

Before you get started, make sure you have already cloned the :ref:`Splam GitHub repository <install-from-source>`. We provide a few examples below:

|

.. _splam-bam-quick:

Example 1: clean up alignment files  (:code:`BAM`)
-----------------------------------------------------

.. code-block:: bash

    $ cd test

    # Step 1: extract splice junctions in the alignment file
    $ splam extract -P SRR1352129_chr9_sub.bam -o tmp_out_alignment

    # Step 2: score all the extracted splice junctions
    $ splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out_alignment tmp_out_alignment/junction.bed

    #Step 3: output a cleaned and sorted alignment file
    $ splam clean -P -o tmp_out_alignment -@ 5   

| 

.. _splam-gff-quick:

Example 2: evaluate annotation files / assembled transcripts (:code:`GFF`)
-----------------------------------------------------------------------------

.. code-block:: bash

    $ cd test

    # Step 1: extract introns in the annotation
    $ splam extract refseq_40_GRCh38.p14_chr_fixed.gff -o tmp_out_annotation

    # Step 2: score introns in the annotation
    $ splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out_annotation tmp_out_annotation/junction.bed

    #Step 3: output statistics of each transcript
    $ splam clean -o tmp_out_annotation -t 0.8

| 

Splam can also :ref:`run on non-human species <generalization-introduction>`. 

.. _splam-generalization-example:

Example of evaluating mouse annotation files (:code:`GFF`)
----------------------------------------------------------------------

.. code-block:: bash

    $ cd test

    # Step 1: extract introns in the annotation
    splam extract mouse_chr19_subset.gff -o tmp_out_generalization

    # Step 2: score introns in the annotation
    splam score -A GRCm39_assembly_report.txt -G mouse_chr19.fa -m ../model/splam_script.pt -o tmp_out_generalization tmp_out_generalization/junction.bed

    #Step 3: output statistics of each transcript
    splam clean -o tmp_out_generalization -t 0.8

|

.. _google-colab:

Try Splam on Google Colab
+++++++++++++++++++++++++++++++++++

We created some reproducible and easy-to-run Splam examples on Google Colab. It's a good starting point, so go ahead and check them out!


.. image:: https://colab.research.google.com/assets/colab-badge.svg
    :target: https://colab.research.google.com/github/Kuanhao-Chao/splam/blob/main/notebook/splam_example.ipynb


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