
|

.. _gene_transcript_counting:

Gene / Transcript counting
============================

The counts of protein-coding and non-coding genes and transcripts reported in this study were calculated as follows:

Gene counting
---------------
Gene features were classified as "**protein-coding**", "**non-coding**", and "**others**" based on the "gene_biotype" attribute in NCBI's RefSeq :cite:p:`o2016reference` or the "gene_type" attribute in EMBL-EBI’s Ensembl/GENCODE :cite:p:`martin2023ensembl` and CHESS 3 :cite:p:`varabyou2023chess`. 

More specifically, a gene feature was categorized as a protein-coding gene if the feature type was "gene" and its type attribute was "protein_coding"; a gene was categorized as non-coding if its type attribute was either "lncRNA" or "ncRNA".

Transcript counting
-----------------------
Transcript features were also categorized into "**protein-coding**", "**non-coding**", and "**others**". The criteria for classification were based on both (1) the type of transcript and (2) the type of its parent gene. 

A transcript was counted as protein-coding if its feature type was "mRNA", and its parent gene was a protein-coding gene; a transcript was classified as non-coding if its feature type was either "lncRNA" or "ncRNA" and its parent gene was also a non-coding gene.


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