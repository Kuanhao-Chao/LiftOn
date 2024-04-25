
|

.. _output_files:

Output files
=====================


After running LiftOn, you will obtain a LiftOn GFF3 file and a :code:`lifton_output/` directory. More details are shown in the output directory hierarchy below. This page provides guidance on how to interpret your results.

.. admonition:: LiftOn output directory hierarchy 
   :class: note


   .. code-block:: 

      # LiftOn output directory hierarchy 
      .
      ├── lifton.gff3
      └── lifton_output
         |
         ├── score.txt
         |
         ├── liftoff
         │   ├── liftoff.gff3
         │   ├── liftoff.gff3_db
         │   └── unmapped_features.txt
         |
         ├── miniprot
         │   ├── miniprot.gff3
         │   └── miniprot.gff3_db
         |
         ├── intermediate_files
         │   ├── proteins.fa
         │   ├── reference_all_genes.fa
         │   ├── reference_all_to_target_all.sam
         │   ├── ref_feature.txt
         │   ├── ref_transcript.txt
         │   └── transcripts.fa
         |
         └── stats
            ├── extra_copy_features.txt
            ├── mapped_feature.txt
            ├── mapped_transcript.txt
            └── unmapped_features.txt

|
|


lifton.gff3
--------------
This is the main output of LiftOn software. It is an annotation file of the target genome in GFF3 format. Following is an example of a gene locus. For more details about GFF3, check `GFF3 file format <https://useast.ensembl.org/info/website/upload/GFF3.html>`_. 

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
      6. The source of the annotation:    
      
         * 'Liftoff', 'miniprot', or 'LiftOn_chaining_algorithm'.

      7. Mutation types

         * 'synonymous', 'nonsynonymous', 'inframe_insertion', 'inframe_deletion', 'frameshift', 'stop_codon_gain', 'stop_codon_loss', and 'start_codon_loss".

      8. Transcript locus coordinate: :code:`<chromosome>:<start>-<end>`

.. dropdown:: Example
   :animate: fade-in-slide-down
   :container: + shadow
   :title: bg-light font-weight-bolder
   :body: bg-light text-left

   .. code-block:: plain-text

      rna-NR_132385.2	0	0	0.9991364421416234	0	Liftoff	non_coding	chr22:16314791-16352180
      rna-NM_014406.5	0.9982078853046595	0.996415770609319	0.9985443959243085	0.9982078853046595	LiftOn_chaining_algorithm	frameshift	chr22:17267300-17269359
      rna-NR_134584.1	0	0	0.9990338164251208	0	Liftoff	non_coding	chr22:17423569-17425137
      rna-NM_001318251.3	0.9910714285714286	0.9910714285714286	0.9981549815498155	0.9910714285714286	LiftOn_chaining_algorithm	frameshift	chr22:17460159-17498391
      rna-NM_001386955.1	0.9910714285714286	0.9910714285714286	0.9987878787878788	0.9910714285714286	LiftOn_chaining_algorithm	nonsynonymous	chr22:17460159-17502102
      rna-NM_001386956.1	0.9910714285714286	0.9910714285714286	0.997489014438167	0.9910714285714286	LiftOn_chaining_algorithm	frameshift	chr22:17460159-17502102
      rna-NM_001386957.1	0.9910714285714286	0.9910714285714286	0.998533724340176	0.9910714285714286	LiftOn_chaining_algorithm	frameshift	chr22:17460159-17502102

|


2. liftoff/
+++++++++++++++++++++++++++++++++++

The liftoff GFF3 annotatation generated during the LiftOn process.

|

3. miniprot/
+++++++++++++++++++++++++++++++++++

The miniprot GFF3 file generated during the LiftOn process.

|

4. intermediate_files/
+++++++++++++++++++++++++++++++++++

In this directory, it stores all intermdeiate files, including protein sequences (FASTA), gene seuqence to genome alignment (SAM), transcript sequences (FASTA), the type of the reference gene (coding or non-coding), and the type of the reference transcript (coding or non-coding).

         
|


5. stats/mapped_feature.txt
+++++++++++++++++++++++++++++++++++

It is a TSV file summarizing the number of features being mapped from the reference genome to the target genome.

.. admonition:: Column definition
   :class: note

      1. Feature ID
      2. The number of feature copy
      3. Feeature type: `coding`, `non-coding`, or `other`

.. dropdown:: Example
   :animate: fade-in-slide-down
   :container: + shadow
   :title: bg-light font-weight-bolder
   :body: bg-light text-left

   .. code-block:: plain-text

      gene-LOC105379428	1	non-coding
      gene-LOC124905168	1	other
      gene-OR11H1	1	coding
      gene-LOC112268291	1	non-coding
      gene-POTEH	1	coding
      gene-POTEH-AS1	1	non-coding
      gene-PSLNR	1	non-coding

|

6. stats/mapped_transcript.txt
+++++++++++++++++++++++++++++++++++

It is a TSV file summarizing the number of transcripts being mapped from the reference genome to the target genome.

.. admonition:: Column definition
   :class: note

      1. Transcript ID
      2. The number of transcript copy
      3. Transcript type: `coding`, `non-coding`, or `other`

.. dropdown:: Example
   :animate: fade-in-slide-down
   :container: + shadow
   :title: bg-light font-weight-bolder
   :body: bg-light text-left

   .. code-block:: plain-text

      rna-XM_011546114.3      1       coding
      rna-NM_001025161.3	2	coding
      rna-XR_001755413.2	3	non-coding
      rna-XR_001755415.2	3	non-coding
      rna-XR_002958735.2	3	non-coding
      rna-NR_110761.1	3	non-coding
      id-LOC124905154	2	others

|

7. stats/extra_copy_features.txt
+++++++++++++++++++++++++++++++++++

Similar to mapped_feature.txt, this TSV file summarizes the number of additional copies of genes and indicates whether they are coding or non-coding.

.. admonition:: Column definition
   :class: note

      1. Feature ID
      2. The number of feature copy
      3. Feeature type: `coding`, `non-coding`, or `other`



.. dropdown:: Example
   :animate: fade-in-slide-down
   :container: + shadow
   :title: bg-light font-weight-bolder
   :body: bg-light text-left

   .. code-block:: plain-text

      gene-LOC124905154	2	other
      gene-LOC107984037	3	non-coding
      gene-LOC102723769	3	non-coding
      gene-LOC107987323	3	non-coding
      gene-LOC105372858	4	non-coding
      gene-GGTLC3	3	coding
      gene-LOC107985584	9	non-coding
      gene-MIR650	2	other
      gene-LINC02556	5	non-coding
      gene-CYP2D6	2	coding
      gene-LOC101927372	2	non-coding

|


8. stats/unmapped_features.txt
+++++++++++++++++++++++++++++++++++

It is a TSV file summarizing unmapped gene IDs and their types.

.. admonition:: Column definition
   :class: note

      1. Gene ID
      2. Feeature type: `coding`, `non-coding`, or `other`


.. dropdown:: Example
   :animate: fade-in-slide-down
   :container: + shadow
   :title: bg-light font-weight-bolder
   :body: bg-light text-left

   .. code-block:: plain-text

      gene-LOC124905174	non-coding
      gene-LOC124900482	non-coding
      gene-FAM246B	coding
      gene-RIMBP3C	coding
      gene-IGLJ1	other
      gene-IGLJ2	other
      gene-IGLJ3	other
      gene-IGLJ4	other
      gene-IGLJ5	other
      gene-IGLJ6	other
      gene-IGLJ7	other

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