
|

.. _quick-start:

Quick Start Guide
=================

This page offers straightforward quick-start information on using LiftOn to map RefSeq v220 annotations on :code:`chr22:1-50755568` from GRCh38 version 40 patch 14 to T2T-CHM13 v2.0. If you haven't already, please follow the steps outlined on the :ref:`Installation` page to install and load LiftOn.


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

    $ lifton -g GRCh38_chr22.gff3 -o GRCh38_2_CHM13/GRCh38_2_CHM13_lifton.gff3 -copies -sc 0.95 chm13_chr22.fa GRCh38_chr22.fa


After this step, you will obtain a directory, :code:`GRCh38_2_CHM13/`, which contains the output annotation file :code:`GRCh38_2_CHM13_lifton.gff3` in GFF3 format. Additionally, a directory :code:`lifton_output/` includes various outputs such as sequence identity scores, Liftoff and miniprot outputs, statistics files about unmapped features, extra copies, intermediate files, and more. Further explanations of the output file hierarchy are provided in the :ref:`output files section <output_files>`.



Interpreting outputs on the terminal
+++++++++++++++++++++++++++++++++++++

After running the command, you will see the following output on the terminal:

.. code-block:: txt

    *********************************************
    * Total features in reference	: 890
    * Lifted feature			: 879 (437 + 333 + 109)
        * Protein-coding feature	: 437
        * Non-coding feature		: 333
        * Other feature			: 109
    * Missed feature			: 11

    * Total features in target		: 906 (440 + 355 + 111)
        * Protein-coding feature	: 440 (435 + 5)
            * single copy		: 435
            * > 1 copy			: 2, 5 in total
        * Non-coding feature		: 355 (326 + 29)
            * single copy		: 326
            * > 1 copy			: 7, 29 in total
        * Other feature			: 111 (107 + 4)
            * single copy		: 107
            * > 1 copy			: 2, 4 in total
    *********************************************

There are two parts to the output:

1. It shows the total number of features in the reference assembly, the number of features lifted, and the number of features missed. The lifted features are further categorized into protein-coding, non-coding, and other features. The missed features are those that could not be transferred to the target assembly.

2. It shows the total number of features in the target assembly, including extra copies identified by LiftOn, which are also categorized into protein-coding, non-coding, and other features.

For more detailed information about our feature counting approach, check out the :ref:`gene / transcript counting section <gene_transcript_counting>`.

|

.. _google-colab:

Try LiftOn on Google Colab
+++++++++++++++++++++++++++++++++++

We created a reproducible and easy-to-run LiftOn example on Google Colab. It's a great starting point to try LiftOn, so go ahead and check it out!

.. image:: https://colab.research.google.com/assets/colab-badge.svg
    :target: https://colab.research.google.com/github/Kuanhao-Chao/LiftOn/blob/main/notebook/lifton_example.ipynb


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