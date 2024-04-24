
|

.. _gene_transcript_counting:

Gene / Transcript counting
============================

After running LiftOn, you will obtain a LiftOn GFF3 file and a :code:`lifton_output/` directory. More details are shown in the output directory hierarchy below. This page provides guidance on how to interpret your results.

.. admonition:: LiftOn output directory hierarchy 
   :class: note


   .. code-block:: 

      # LiftOn output directory hierarchy 
      .
      ├── lifton.GFF3
      └── lifton_output
         ├── intermediate_files
         │   ├── proteins.fa
         │   ├── proteins_truncated.fa
         │   ├── reference_all_genes.fa
         │   ├── reference_all_to_target_all.sam
         │   └── transcripts.fa
         │
         ├── liftoff
         │   ├── liftoff.GFF3
         │   └── liftoff.GFF3_db
         │
         ├── miniprot
         │   ├── miniprot.GFF3
         │   └── miniprot.GFF3_db
         │
         ├── extra_copy_features.txt
         ├── score.txt
         └── unmapped_features.txt


|
|


lifton.GFF3
--------------
This is the main output of LiftOn software. It is an annotation file of the target genome in GFF3 format. Following is an example of a gene locus. More details about `GFF3 file format <https://useast.ensembl.org/info/website/upload/GFF3.html>`_. 

.. dropdown:: Example
   :animate: fade-in-slide-down
   :container: + shadow
   :title: bg-light font-weight-bolder
   :body: bg-light text-left

   .. code-block:: plain-text

      NW_020825194.1  LiftOn  gene    15799   29728   .       +       .       ID=gene-Dmel_CG1483;Dbxref=FLYBASE:FBgn0002645,GeneID:43765;Name=Map205;cyt_map=100E3-100E3;description=Microtubule-associated protein 205;gbkey=Gene;gen_map=3-103 cM;gene=Map205;gene_biotype=protein_coding;gene_synonym=205-kDa MAP,205K MAP,205kD MAP,205kDa MAP,CG1483,Dmel\CG1483,map205,MAP205,MAP4;locus_tag=Dmel_CG1483
      NW_020825194.1  LiftOn  mRNA    15799   29728   .       +       .       ID=rna-NM_001276225.1;Parent=gene-Dmel_CG1483;Dbxref=FLYBASE:FBtr0334299,GeneID:43765,GenBank:NM_001276225.1,FLYBASE:FBgn0002645;Name=NM_001276225.1;Note=Map205-RC%3B Dmel\Map205-RC%3B CG1483-RC%3B Dmel\CG1483-RC;gbkey=mRNA;gene=Map205;locus_tag=Dmel_CG1483;orig_protein_id=gnl|FlyBase|CG1483-PC|gb|AGB96532;orig_transcript_id=gnl|FlyBase|CG1483-RC;product=Microtubule-associated protein 205%2C transcript variant C;transcript_id=rna-NM_001276225.1;mutation=frameshift;protein_identity=0.795;dna_identity=0.793;status=LiftOn_chaining_algorithm
      NW_020825194.1  LiftOn  exon    15799   17891   .       +       .       Parent=rna-NM_001276225.1
      NW_020825194.1  LiftOn  exon    25963   26161   .       +       .       Parent=rna-NM_001276225.1
      NW_020825194.1  LiftOn  exon    26586   29728   .       +       .       Parent=rna-NM_001276225.1
      NW_020825194.1  LiftOn  CDS     16236   17891   2375    +       0       Parent=rna-NM_001276225.1
      NW_020825194.1  LiftOn  CDS     25963   26161   221     +       0       Parent=rna-NM_001276225.1
      NW_020825194.1  LiftOn  CDS     26586   28009   1649    +       2       Parent=rna-NM_001276225.1

|
|

lifton_output/
---------------

This directory contains all LiftOn outputs, including the following:


1. score.txt
+++++++++++++++++++++++++++++++++++

It is a tsv file summarizing the LiftOn results.

.. admonition:: Column definition
   :class: note

   1. Transcript ID
   2. Liftoff transcript protein sequence identity: :math:`0-1`
   3. miniprot transcript protein sequence identity: :math:`0-1`
   4. LiftOn transcript DNA sequence identity: :math:`0-1`
   5. LiftOn transcript protein sequence identity: :math:`0-1`
   6. Status of the annotation: 
   
      * 'Liftoff_identical', 'Liftoff_synonymous', 'Liftoff_truncated', 'Liftoff_nc_transcript', 'Liftoff_no_ref_protein', 'LiftOn_chaining_algorithm', 'miniprot_identical', 'miniprot_truncated'


   7. Mutation types

      * 'synonymous', 'non-synonymous', 'in-frame insertion', 'in-frame deletion', 'frameshift', 'stop codon gain', 'stop codon loss', and 'start codon loss".

   8. Transcript locus coordinate: :code:`<chromosome>:<start>-<end>`

.. dropdown:: Example
   :animate: fade-in-slide-down
   :container: + shadow
   :title: bg-light font-weight-bolder
   :body: bg-light text-left

   .. code-block:: plain-text

      rna-NM_176598.2 0.9146141215106732      0.9507389162561576      0.8503959276018099      0.9507389162561576      LiftOn_chaining_algorithm       frameshift      NW_020825194.1:114373-268723
      rna-NM_176599.2 0.9090909090909091      0.9448051948051948      0.8489566081483935      0.9496753246753247      LiftOn_chaining_algorithm       frameshift;start_lost   NW_020825194.1:122632-268723
      rna-NM_176600.3 0.9146141215106732      0.9507389162561576      0.8896146309601568      0.9507389162561576      LiftOn_chaining_algorithm       frameshift      NW_020825194.1:112640-268723
      rna-NM_176601.3 0.9146141215106732      0.9507389162561576      0.9075364154528183      0.9507389162561576      LiftOn_chaining_algorithm       frameshift      NW_020825194.1:112640-268723

|

2. extra_copy_features.txt
+++++++++++++++++++++++++++++++++++

It is a TSV file summarizing the number of copies of a specific gene and indicating whether it is a coding or non-coding gene.

.. admonition:: Column definition
   :class: note

   1. Gene ID
   2. The number of gene copy
   3. coding or non-coding tag



.. dropdown:: Example
   :animate: fade-in-slide-down
   :container: + shadow
   :title: bg-light font-weight-bolder
   :body: bg-light text-left

   .. code-block:: plain-text

      gene-Dmel_CG32498       2       coding
      gene-Dmel_CG6998        2       coding
      gene-Dmel_CR32748       2       non-coding
      gene-Dmel_CG34417       2       coding
      gene-Dmel_CG1343        2       coding
      gene-Dmel_CR32615       2       non-coding
      gene-Dmel_CG46317       2       coding
      gene-Dmel_CG6340        2       coding
      gene-Dmel_CG46306       2       coding
      gene-Dmel_CG5004        2       coding

|

3. unmapped_features.txt
+++++++++++++++++++++++++++++++++++

It is a TSV file summarizing unmapped gene ID.

.. admonition:: Column definition
   :class: note

   1. Gene ID


.. dropdown:: Example
   :animate: fade-in-slide-down
   :container: + shadow
   :title: bg-light font-weight-bolder
   :body: bg-light text-left

   .. code-block:: plain-text

      gene-Dmel_CR40469
      gene-Dmel_CR43552
      gene-Dmel_CR45473
      gene-Dmel_CG32817
      gene-Dmel_CR43519
      gene-Dmel_CR45474
      gene-Dmel_CR45475
      gene-Dmel_CR46283
      gene-Dmel_CR44469
      gene-Dmel_CG13359
      gene-Dmel_CG14634
      gene-Dmel_CR45476

|

4. miniprot/
+++++++++++++++++++++++++++++++++++

The miniprot GFF3 file generated during the LiftOn process.

|

5. liftoff/
+++++++++++++++++++++++++++++++++++

The liftoff GFF3 annotatation generated during the LiftOn process.

|

6. intermediate_files/
+++++++++++++++++++++++++++++++++++

In this directory, it stores all intermdeiate files, including protein sequences (FASTA), truncated protein sequences (FASTA), gene seuqence to genome alignment (SAM), and transcript sequences (FASTA). 

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