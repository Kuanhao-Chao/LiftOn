.. raw:: html

    <script type="text/javascript">

        let mutation_fuc = function(mutations) {
            var dark = document.body.dataset.theme == 'dark';

            if (document.body.dataset.theme == 'auto') {
                dark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
            }
            
            document.getElementsByClassName('sidebar_ccb')[0].src = dark ? './_static/JHU_ccb-white.png' : "./_static/JHU_ccb-dark.png";
            document.getElementsByClassName('sidebar_wse')[0].src = dark ? './_static/JHU_wse-white.png' : "./_static/JHU_wse-dark.png";



            for (let i=0; i < document.getElementsByClassName('summary-title').length; i++) {
                console.log(">> document.getElementsByClassName('summary-title')[i]: ", document.getElementsByClassName('summary-title')[i]);

                if (dark) {
                    document.getElementsByClassName('summary-title')[i].classList = "summary-title card-header bg-dark font-weight-bolder";
                    document.getElementsByClassName('summary-content')[i].classList = "summary-content card-body bg-dark text-left docutils";
                } else {
                    document.getElementsByClassName('summary-title')[i].classList = "summary-title card-header bg-light font-weight-bolder";
                    document.getElementsByClassName('summary-content')[i].classList = "summary-content card-body bg-light text-left docutils";
                }
            }

        }
        document.addEventListener("DOMContentLoaded", mutation_fuc);
        var observer = new MutationObserver(mutation_fuc)
        observer.observe(document.body, {attributes: true, attributeFilter: ['data-theme']});
        console.log(document.body);
    </script>
    <link rel="preload" href="./_images/jhu-logo-dark.png" as="image">

|

.. _main:

LiftOn's tutorial
*************************


.. raw:: html

    <embed>
        <div class="sidebar-logo-container" style="padding-bottom:-10px">
            <img class="sidebar-logo only-light" src="_static/LiftOn_color.png" alt="Light Logo">
            <img class="sidebar-logo only-dark" src="_static/LiftOn_white.png" alt="Dark Logo">
        </div>
    </embed>

.. image:: https://img.shields.io/badge/License-GPLv3-yellow.svg
    :target: https://img.shields.io/badge/License-GPLv3-yellow.svg

.. image:: https://img.shields.io/badge/version-v.0.0.1-blue
    :target: https://img.shields.io/badge/version-v.0.0.1-blue

.. .. image:: https://img.shields.io/github/downloads/Kuanhao-Chao/LiftOn/total.svg?style=social&logo=github&label=Download
..     :target: https://img.shields.io/github/downloads/Kuanhao-Chao/LiftOn/total.svg?style=social&logo=github&label=Download

.. .. image:: https://img.shields.io/badge/platform-macOS_/Linux-green.svg
..     :target: https://img.shields.io/badge/platform-macOS_/Linux-green.svg

.. .. image:: https://colab.research.google.com/assets/colab-badge.svg
..     :target: https://colab.research.google.com/github/Kuanhao-Chao/LiftOn/blob/main/notebook/LiftOn_example.ipynb

| 

.. What is LiftOn?
.. ==================


LiftOn is a homology-based lift-over tool designed to accurately map annotations in GFF or GTF between assemblies. This tool is constructed on the foundation of the fantastic `Liftoff <https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=true>`_ and integrates `miniprot <https://academic.oup.com/bioinformatics/article/39/1/btad014/6989621>`_ protein alignment for an improved homology-based lift-over process.

.. LiftOn improves the protein-coding gene annotations and improves same or closely-related species lift-over!

.. LiftOn enhances protein-coding gene annotations and facilitates lift-over for the same or closely-related species.


.. lift-over annotator that takes `Liftoff <https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=true>`_ and `miniprot <https://academic.oup.com/bioinformatics/article/39/1/btad014/6989621>`_ GFF files as input. It accurately generates gene annotations, with a particular focus on protein-coding genes. LiftOn takes consensus from both sources and generates optimal annotations that outperform both `Liftoff <https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=true>`_ and `miniprot <https://academic.oup.com/bioinformatics/article/39/1/btad014/6989621>`_!

Why LiftOn❓
==================

1. **Burgeoning number of genome assemblies**: As of December 2023, among the 15,578 distinct eukaryotic genomes, only 1,111 have been annotated (`Eukaryotic Genome Annotation at NCBI <https://www.ncbi.nlm.nih.gov/genome/annotation_euk/#graphs>`_). More and more high quality assemblies are generated. We need to accurately annotate them!

2. **Improved protein-coding gene mapping**: The popular `Liftoff <https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=true>`_ map genes only based on the DNA alignment. With the protein-to-genome alignment, LiftOn is able to further improve the lift-over annotation! LiftOn is able to improve the current released T2T-CHM13 annotation (`JHU RefSeqv110 + Liftoff v5.1 <https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz>`_). 

3. **Improved distant species lift-over**: See mapping from *Mus musculus* to *Rattus norvegicus* and *Drosophila melanogaster* to *Drosophila erecta*.

LiftOn is free, it's open source, it's easy to install , and it's in Python!

|

Who is it for❓
====================================

1. If you have sequenced and assembled a new genome and need to annotate it, LiftOn is the ideal choice for generating annotations.
2. If you wish to utilize the finest CHM13 annotation, you can run LiftOn! We have also pre-generated the `T2T_CHM13_LiftOn.gff3 <https://khchao.com>`_ file for your convenience.


|

What does LiftOn do❓
====================================

.. LiftOn takes GFF files from Liftoff and miniprot and reference protein sequences in a FASTA file, and generates a new annotation file in GFF format. 
.. LiftOn works on the same and closely-related species. 
.. <!-- We also tested LiftOn by lifting-over annotations from human to mouse, and it also does pretty good job to find the optimal protein annotations. However, there are false positives or -->

LiftOn is designed for individuals who would like to annotate a new assembly, referred to as target **Genome** :math:`T`.

The first step is to select a well-annotated genome along with its annotation, denoted as reference **Genome** :math:`R` and **Annotation** :math:`R_A`. 

The process begins by extracting protein sequences annotated by Liftoff and miniprot. These sequences are then aligned with full-length reference proteins. For each gene locus, LiftOn employs the *chaining algorithm* that compares each section of the protein alignments from Liftoff and miniprot. This algorithm corrects errors in exon and CDS boundaries, resulting in the better protein annotations that preserves the longest matching proteins.

* **Input**: 
    1. target **Genome** :math:`T` in FASTA 
    2. reference **Genome** :math:`R` in FASTA  
    3. reference **Annotation** :math:`R_A` in GFF3  
* **Output**: 
    1. LiftOn annotation file in GFF3
    2. Protein sequence identities & mutation types
    3. Features with extra copies
    4. Unmapped features

|

LiftOn's limitation
==================================
LiftOn's *chaining algorithm* currently only utilizes miniprot alignment results to fix the Liftoff annotation. However, it has the potential to expand its capabilities to chain together multiple protein-based annotation files or aasembled RNA-Seq transcripts. 

The LiftOn *chaining algorithm* now does not support multi-threading. This functionality stands as our next targeted feature on the development horizon!

|

User support
============
Please go through the :ref:`documentation <table-of-contents>` below first. If you have questions about using the package, a bug report, or a feature request, please use the GitHub issue tracker here:

https://github.com/Kuanhao-Chao/LiftOn/issues

|

Key contributors
================

LiftOn was designed and developed by `Kuan-Hao Chao <https://khchao.com/>`_.  This documentation was written by `Kuan-Hao Chao <https://khchao.com/>`_.

|

.. _table-of-contents:

Table of contents
==================

.. hide-toc:: true

.. toctree::
    :maxdepth: 2
    
    content/installation
    content/quickstart

.. toctree::
    :caption: Examples
    :maxdepth: 2
    
    content/same_species_liftover/index
    content/close_species_liftover/index
    content/distant_species_liftover/index


.. toctree::
    :caption: Info
    :maxdepth: 2
    
    content/output_explanation
    content/behind_scenes
    content/how_to_page
    content/function_manual
    content/changelog
    content/license
    content/contact

|
|
|
|
|

.. content/installation
..    content/quickstart
..    content/liftover_GRCh38_2_T2TCHM13
..    content/liftover_bee_insect
..    content/liftover_arabidopsis_plant
..    content/liftover_drosophila_erecta
..    content/liftover_mouse_2_rat
..    content/behind_scenes
..    content/how_to_page
..    content/function_manual
..    content/license
..    content/contact


.. image:: ./_images/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image only-light
   :align: center

.. image:: ./_images/jhu-logo-white.png
   :alt: My Logo
   :class: logo, header-image only-dark
   :align: center