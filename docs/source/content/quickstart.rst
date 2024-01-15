.. _quick-start:

Quick Start Guide
=================

This page offers straightforward quick-start information on using LiftOn to map RefSeq v110 annotations on :code:`chr22:1-50755568` from GRCh38 version 40 patch 14 to T2T-CHM13 v2.0. If you haven't already, please follow the steps outlined on the :ref:`Installation` page to install and load LiftOn.

Before you get started, make sure you have already cloned the `LiftOn GitHub repository <https://github.com/Kuanhao-Chao/LiftOn>`_. We provide an example in `test/lifton_chr22_example.sh <https://github.com/Kuanhao-Chao/LiftOn/tree/main/test/lifton_chr22_example.sh>`_.


|

.. _super-quick-start:

Super-Quick Start (one-liner)
+++++++++++++++++++++++++++++++++++


LiftOn maps annotations from one assembly to another. To run LiftOn, all you need are three files:

1. A target assembly (**Genome** :math:`T`, FASTA Format):  `chm13_chr22.fa <https://github.com/Kuanhao-Chao/LiftOn/tree/main/test/chm13_chr22.fa>`_
2. A reference assembly (**Genome** :math:`R`, FASTA Format): `GRCh38_chr22.fa <https://github.com/Kuanhao-Chao/LiftOn/tree/main/test/GRCh38_chr22.fa>`_
3. A reference annotation (**Annotation** :math:`R_A`, GFF3 Format): `GRCh38_chr22.gff3 <https://github.com/Kuanhao-Chao/LiftOn/tree/main/test/GRCh38_chr22.gff3>`_

Run the following commands:

.. code-block:: bash

    $ cd test

    $ lifton -g GRCh38_chr22.gff3 -dir GRCh38_2_CHM13 -o GRCh38_2_CHM13_lifton.gff3 -copies chm13_chr22.fa GRCh38_chr22.fa


After this step, you will get an annotation file :code:`GRCh38_2_CHM13_lifton.gff3` in GFF3 format, and a directory :code:`lifton_output/` that contains various outputs such as sequence identity scores, unmapped features, extra copies, intermediate files, and more.

.. GRCh38_2_CHM13_lifton.gff3

.. extra_copy_features.txt
.. unmapped_features.txt
.. score.txt
.. gene.txt
.. miniprot
.. liftoff
.. intermediate_files
.. GRCh38_chr22.gff3_db
.. GRCh38_2_CHM13


We will further explain in :ref:`output_files` section.

|

.. _google-colab:

Try LiftOn on Google Colab
+++++++++++++++++++++++++++++++++++

We created a reproducible and easy-to-run LiftOn example on Google Colab. It's a good starting point, so go ahead and check it out!


.. image:: https://colab.research.google.com/assets/colab-badge.svg
    :target: https://colab.research.google.com/github/Kuanhao-Chao/LiftOn/blob/main/notebook/LiftOn_example.ipynb


|

Congratulations! You have successfully installed and run LiftOn. For more detailed analysis explaination and file format, please check:

.. seealso::
    
    * :ref:`same_species-section`

    * :ref:`close_species-section`

    * :ref:`distant_species-section`

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