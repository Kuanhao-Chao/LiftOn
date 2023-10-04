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

Splam's tutorial
*************************

.. _splam-logo:
.. figure::  ./_static/logo.png
   :align:   center
   :scale:   13 %

|

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
    :target: https://img.shields.io/badge/License-MIT-yellow.svg

.. image:: https://img.shields.io/badge/version-v.1.0.2-blue
    :target: https://img.shields.io/badge/version-v.1.0.2-blue

.. image:: https://img.shields.io/github/downloads/Kuanhao-Chao/splam/total.svg?style=social&logo=github&label=Download
    :target: https://img.shields.io/github/downloads/Kuanhao-Chao/splam/total.svg?style=social&logo=github&label=Download

.. image:: https://img.shields.io/badge/platform-macOS_/Linux-green.svg
    :target: https://img.shields.io/badge/platform-macOS_/Linux-green.svg

.. image:: https://colab.research.google.com/assets/colab-badge.svg
    :target: https://colab.research.google.com/github/Kuanhao-Chao/splam/blob/main/notebook/splam_example.ipynb

| 

.. What is Splam?
.. ==================

Splam is a splice junction recognition model based on a deep residual convolutional neural network that offers **fast and precise** assessment of splice junctions. 

Why Splam❓
==================

1. **We need a tool to evaluate splice junctions & spliced alignments.** Thousands of RNA-Seq datasets are generated every day, but there are no tools available for cleaning up spurious spliced alignments in these data. Splam addresses this problem!
2. **Splam-cleaned alignments lead to improved transcript assembly, which, in turn, may enhance all downstream RNA-Seq analyses**, including transcript quantification, differential gene expression analysis, and more!

|

Who is it for❓
====================================

If you are **(1) doing RNA-Seq data analysis** or **(2) seeking a trustworthy tool to evaluate splice junctions (introns)**, then Splam is the tool that you are looking for!

|

What does Splam do❓
====================================

There are two main use case scenarios:

1. Improving your **alignment file**. Splam evaluates the quality of spliced alignments and removes those containing spurious splice junctions. This significantly enhances the quality of downstream transcriptomes [:ref:`Link <alignment-detailed-section>`].

2. Evaluating the quality of introns in your **annotation file or assembled transcripts** [:ref:`Link <annotation-detailed-section>`].


Splam is free, it's open source, it's a lightweight deep learning model, and it's in Python!

.. It was trained on donor and acceptor pairs combined and focuses on a narrow window of 400 basepairs surrounding each splice site, inspired by the understanding that the splicing process primarily depends on signals within this specific region.

|

Main features
=============

* **Biologically inspired training process**: Splam was trained on combined donor and acceptor pairs, emulating the behavior of the spliceosome, with a specific emphasis on a narrow window of 400 base pairs surrounding each splice site. This approach is inspired by the understanding that the splicing process predominantly relies on signals within this specific region.
* **Generalization to non-human species**: Splam was trained exclusively using human splice junctions; however, we have demonstrated that it performs well on chimpanzee, mouse, and even the flowering plant *Arabidopsis thaliana*.
* **Python & C++ integration**: We have taken care of all the engineering work for you! Splam is easy to install and runs efficiently due to its underlying C++ implementation. You can install and run Splam with just one simple command.
* **Run Splam in three steps**: With just three lines of code, you can obtain a new alignment file that is cleaned and sorted.
* **Pytorch implementation**: Splam is implemented and trained using the popular and reliable PyTorch framework.



|

What Splam doesn't do
==================================

Splam does not have the feature to scan through the genome and score every potential splice site. Some splice site prediction tools take a DNA sequence and predict the splice sites within it, such as `SpliceAI <https://github.com/Illumina/SpliceAI>`_. However, SpliceAI was only trained on a single canonical transcript for each protein-coding gene while disregarding alternative isoforms. Splam takes a different approach, focusing on predicting at the "splice junction level" rather than the "transcript level." Splam was trained on a large collection of human splices sites taken from both "canonical" and alternative isoforms.

|

User support
============
Please go through the :ref:`documentation <table-of-contents>` below first. If you have questions about using the package, a bug report, or a feature request, please use the GitHub issue tracker here:

https://github.com/Kuanhao-Chao/splam/issues

|

Key contributors
================

Splam deep residual convolutional neural network was trained using the PyTorch framework by Kuan-Hao Chao. Kuan-Hao Chao also implemented the package that applies Splam to evaluate annotation files and clean up alignment files. This documentation was written by Kuan-Hao Chao and Alan Mao.

|

.. _table-of-contents:

Table of contents
==================

.. toctree::
   :maxdepth: 2

   content/installation
   content/quickstart
   content/alignment_evaluation
   content/annotation_evaluation
   content/generalization
   content/behind_scenes
   content/how_to_page
   content/function_manual
   content/license
   content/contact

|
|
|
|
|


.. image:: ./_images/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image only-light
   :align: center

.. image:: ./_images/jhu-logo-white.png
   :alt: My Logo
   :class: logo, header-image only-dark
   :align: center