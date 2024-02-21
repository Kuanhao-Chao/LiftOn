.. raw:: html

    <script type="text/javascript">

        let mutation_lvl_1_fuc = function(mutations) {
            var dark = document.body.dataset.theme == 'dark';

            if (document.body.dataset.theme == 'auto') {
                dark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
            }
            
            document.getElementsByClassName('sidebar_ccb')[0].src = dark ? '../../_static/JHU_ccb-white.png' : "../../_static/JHU_ccb-dark.png";
            document.getElementsByClassName('sidebar_wse')[0].src = dark ? '../../_static/JHU_wse-white.png' : "../../_static/JHU_wse-dark.png";



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
        document.addEventListener("DOMContentLoaded", mutation_lvl_1_fuc);
        var observer = new MutationObserver(mutation_lvl_1_fuc)
        observer.observe(document.body, {attributes: true, attributeFilter: ['data-theme']});
        console.log(document.body);
    </script>


|


.. _same_species_liftover_rice:

Rice
=========================================================================

*Oryza sativa*

Input files
+++++++++++++++++++++++++++++++++++

To run this example, download the following three input files.

* **Input**
    1. target **Genome** :math:`T` in FASTA : `ASM2616768v1_genomic.fna <ftp://ftp.ccb.jhu.edu/pub/LiftOn/rice/ASM2616768v1_genomic.fna>`_ 
    2. reference **Genome** :math:`R` in FASTA : `IRGSP_genomic.fna <ftp://ftp.ccb.jhu.edu/pub/LiftOn/rice/IRGSP_genomic.fna>`_
    3. reference **Annotation** :math:`R_A` in GFF3 : `IRGSP_genomic.gff <ftp://ftp.ccb.jhu.edu/pub/LiftOn/rice/IRGSP_genomic.gff>`_



.. .. important::

..     **We propose running Splam as a new step in RNA-Seq analysis pipeline to score all splice junctions.**

There is only one command you need to run LiftOn:

.. code-block:: bash

    lifton -D -g IRGSP_genomic.gff -o lifton.gff3 -copies ASM2616768v1_genomic.fna IRGSP_genomic.fna


After successfully running LiftOn, you will get the following file and output directory:

* **Output**: 
    1. LiftOn annotation file in GFF3: ftp://ftp.ccb.jhu.edu/pub/LiftOn/rice/lifton.gff3
    2. LiftOn output directory: ftp://ftp.ccb.jhu.edu/pub/LiftOn/rice/lifton_output/

       *  `score.txt <ftp://ftp.ccb.jhu.edu/pub/LiftOn/rice/lifton_output/score.txt>`_
       *  `extra_copy_features.txt <ftp://ftp.ccb.jhu.edu/pub/LiftOn/rice/lifton_output/extra_copy_features.txt>`_
       *  `unmapped_features.txt <ftp://ftp.ccb.jhu.edu/pub/LiftOn/rice/lifton_output/unmapped_features.txt>`_

|
|

Results
+++++++++++++++++++++++++++++++++++

Genome annotation evaluation
------------------------------

Here are some visualization results comparing LiftOn annotation to (1) Liftoff and (2) miniprot annotation. 


First, we calculate the protein sequence identity score for every protein-coding transcript (check :ref:`lifton_sequence_identity` section) for three annotations, LiftOn, Liftoff, and miniprot. 

:numref:`figure-rice_miniprot_vs_liftoff` compares the protein-coding gene mapping of Liftoff, based on DNA alignment, with miniprot, utilizing protein-to-DNA alignment. Dots in the lower right signify transcripts where Liftoff outperformed miniprot in protein sequence identity, while the upper left indicates transcripts where miniprot excelled. LiftOn employs the PM algorithm to enhance annotations in both, achieving improved protein-coding gene annotation, as neither approach dominates the other.

.. _figure-rice_miniprot_vs_liftoff:
.. figure::  ../../_images/rice/Liftoff_miniprot/parasail_identities.png
    :align:   center
    :scale:   25 %

    The scatter plot of protein sequence identity comparing between miniprot (y-axis) and Liftoff (x-axis). Each dot represents a protein-coding transcript.
|

Next, we individually assess LiftOn in comparison to Liftoff and miniprot. In the comparison of LiftOn versus Liftoff (:numref:`figure-rice_lifton_vs_liftoff_vs_miniprot`, left), 519 transcripts demonstrate higher protein sequence identity, with 103 achieving 100% identity. Similarly, in the LiftOn versus miniprot comparison (:numref:`figure-rice_lifton_vs_liftoff_vs_miniprot`, right), 7642 protein-coding transcripts exhibit superior matches, elevating 7412 to identical status relative to the reference.

.. _figure-rice_lifton_vs_liftoff_vs_miniprot:
.. figure::  ../../_images/rice/combined_scatter_plots.png
    :align:   center
    :scale:   21 %

    The scatter plot of protein sequence identity comparing between LiftOn (y-axis) and Liftoff (x-axis) (left) and comparing between LiftOn (y-axis) and miniprot (x-axis) (right).
|

We visualize the transcripts in a 3-D plot, incorporating LiftOn, Liftoff, and miniprot scores (see Figure :numref:`figure-rice_3D_scatter`) to provide a comprehensive comparison of the three tools. If a dot is above the :math:`x=y` plane, it indicates that the protein-coding transcript annotation of LiftOn generates a longer valid protein sequence aligning to the full-length reference protein. The 3-D plot reveals that the majority of dots are above the :math:`x=y` plane, suggesting that LiftOn annotation is better.


.. _figure-rice_3D_scatter:
.. figure::  ../../_images/rice/3d_scatter.png
    :align:   center
    :scale:   30 %

    The 3-D scatter plot of protein sequence identity comparing between LiftOn (y-axis), Liftoff (x-axis), and miniprot (z-axis).

|

Next, we check the distribution of protein sequence identities (see :numref:`figure-rice_frequency_log`). Among the three tools, LiftOn (middle) exhibits the smallest left tail, with 36 protein-coding transcripts having a protein sequence identity of :math:`< 0.4`.

.. _figure-rice_frequency_log:
.. figure::  ../../_images/rice/combined_frequency_log.png
    :align:   center
    :scale:   12 %

    Frequency plots in logarithmic scale of protein sequence identity for Liftoff (left), LiftOn (middle), and miniprot (right) for the results of rice lift-over.

|

Finding extra copies of lift-over features
-------------------------------------------------

LiftOn also has a module to find extra copies by using `intervaltree <https://github.com/chaimleib/intervaltree>`_, `Liftoff <https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=true>`_, and `miniprot <https://academic.oup.com/bioinformatics/article/39/1/btad014/6989621>`_. The Circos plot in :numref:`figure-rice_circos` shows their relative positions between the two genomes. The plot illustrates that the extra copies were predominantly located on the same chromosomes in both IRGSP and ASM2616768. The frequency plot of extra copy features are show in :numref:`figure-rice_extra_copy_fq`.

.. _figure-rice_circos:
.. figure::  ../../_images/rice/circos_plot.png
    :align:   center
    :scale:  16 %

    Circos plot illustrating the locations of extra gene copies found on ASM2616768 (left side) compared to IRGSP (right side). Each line shows the location of an extra copy, and lines are color-coded by the chromosome of the original copy.

|


.. _figure-rice_extra_copy_fq:
.. figure::  ../../_images/rice/extra_cp/frequency.png
    :align:   center
    :scale:  30 %

    Frequency plot for additional gene copy.

|

Finally, we examined the order of protein-coding genes (:numref:`figure-rice_gene_order`) between the two genomes and observed that, as expected, nearly all genes occur in the same order and orientation in both human genomes.

.. _figure-rice_gene_order:
.. figure::  ../../_images/rice/gene_order_plot.png
    :align:   center
    :scale:  30 %

    Protein-gene order plot, with the x-axis representing the reference genome (IRGSP) and the y-axis representing the target genome (ASM2616768). The protein sequence identities are color-coded on a logarithmic scale, ranging from green to red. Green represents a sequence identity score of 1, while red corresponds to a sequence identity score of 0.

|


What's next?
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

Congratulations! You have finished this tutorial.

.. seealso::
    
    * :ref:`behind-the-scenes-splam` to understand how LiftOn is designed
    * :ref:`Q&A` to check out some common questions



|
|
|
|
|

.. image:: ../../_images/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image only-light
   :align: center

.. image:: ../../_images/jhu-logo-white.png
   :alt: My Logo
   :class: logo, header-image only-dark
   :align: center