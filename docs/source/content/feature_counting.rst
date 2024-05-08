
|

.. _gene_transcript_counting:

Gene / Transcript counting
============================

The counts of protein-coding and non-coding genes and transcripts reported in this study were calculated as follows:

Gene counting
---------------
Gene features were classified as "**protein-coding**", "**non-coding**", and "**others**" based on the "gene_biotype" attribute in NCBI's RefSeq :cite:p:`o2016reference` or the "gene_type" attribute in EMBL-EBI’s Ensembl/GENCODE :cite:p:`martin2023ensembl` and CHESS 3 :cite:p:`varabyou2023chess`. Note that if "gene_biotype" or "gene_type" are missing from the attributes of a gene feature, the gene was categorized as "others".

More specifically, a gene feature was categorized as a protein-coding gene if the feature type was "gene" and its type attribute was "protein_coding"; a gene was categorized as non-coding if its type attribute was either "lncRNA" or "ncRNA".

Other gene features, including "Pseudogene" ,"miRNA", "snoRNA", "tRNA",
"V_segment", "snRNA", "J_segment", "misc_RNA", "C_region", etc., were categorized as “others”. The full list is provided in Table below.

|

.. _biotype_table:

.. list-table:: The number of genes for each biotype annotated on the main chromosomes and unplaced contigs in GRCh38 and T2T-CHM13. The T2T-CHM13 annotations were generated using LiftOn by mapping the RefSeq annotations (Release 220) from GRCh38 patch 14 to T2T-CHM13 v2.0. 
   :widths: 33 33 33
   :header-rows: 1

   * - Gene biotype
     - Number of genes in GRC38
     - Number of genes mapped onto T2T_CHM13 by LiftOn
   * - protein_coding
     - 19927
     - 20144
   * - lncRNA
     - 18008
     - 18723
   * - Pseudogene
     - 16896
     - 17855
   * - miRNA
     - 1915
     - 2383
   * - snoRNA
     - 1194
     - 1187
   * - tRNA
     - 453
     - 548
   * - V_segment
     - 239
     - 244
   * - snRNA
     - 153
     - 192
   * - J_segment
     - 98
     - 80
   * - ncRNA
     - 51
     - 49
   * - misc RNA
     - 42
     - 44
   * - C_region
     - 21
     - 23
   * - antisense RNA
     - 19
     - 19
   * - other
     - 13
     - 13
   * - Y RNA
     - 4
     - 7
   * - vault RNA
     - 4
     - 4
   * - scRNA
     - 4
     - 4
   * - telomerase RNA
     - 1
     - 1
   * - RNase P RNA
     - 1
     - 1
   * - RNase MRP RNA
     - 1
     - 1
   * - **Total**
     - **59044**
     - **61522**




Transcript counting
-----------------------
Transcript features were also categorized into "**protein-coding**", "**non-coding**", and "**others**". The criteria for classification were based on both (1) the type of transcript and (2) the type of its parent gene. 

A transcript was counted as protein-coding if its feature type was "mRNA", and its parent gene was a protein-coding gene; a transcript was classified as non-coding if its feature type was either "lncRNA" or "ncRNA" and its parent gene was also a non-coding gene. Other types of transcripts were categorized as "others".


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