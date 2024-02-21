
|

.. _Q&A:

Q & A ...
==========

.. dropdown:: Q: What is LiftOn?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-dark font-weight-bolder howtoclass
    :body: bg-dark text-left

    LiftOn is a homology-based lift-over annotation tool designed to map annotations in GFF or GTF between assemblies.

    LiftOn is built upon the impressive tools, `Liftoff <https://github.com/agshumate/Liftoff>`_ (developed by `Dr. Alaina Shumate <https://scholar.google.com/citations?user=N3tXk7QAAAAJ&hl=en>`_) and `miniprot <https://github.com/lh3/miniprot>`_ (`Dr. Heng Li <http://liheng.org>`_). 

|


.. dropdown:: Q: Why do we need LiftOn?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Rather than repeatedly annotating novel genome assemblies, a more efficient strategy involves transferring genes from well-annotated organisms of the same or closely related species.

    `Liftoff <https://github.com/agshumate/Liftoff>`_, being entirely DNA-based, utilizes minimap2 to align gene loci DNA sequences to the genome and convert gene coordinates to the new assembly. However, when a newly assembled genome deviates significantly from the reference DNA sequence, the alignment may produce transcripts with incorrect protein-coding sequences or erroneous splice sites, posing challenges in annotation, particularly for more distantly related species.

    `miniprot <https://github.com/lh3/miniprot>`_, on the other hand, is exclusively protein-based. This approach has limitations. (1) It cannot capture untranslated regions (UTRs), (2) may miss small exons in cases of long introns, (3) is susceptible to aligning proteins to pseudogenes due to the disregard of intronic sequences, and (4) may combine coding sequences (CDSs) from distinct genes when arranged in tandem along a genome. (5) Additionally, it solely applies to protein-coding transcripts, excluding non-coding genes or other features.

    To overcome these limitations, we created LiftOn, which combines the advantages of both DNA- and protein-based approaches and applies a two-step :ref:`protein-maximization (PM) algorithm <protein-maximization_algorithm>` leading to enhanced protein-coding gene annotation.
    
|


.. dropdown:: Q: How much does LiftOn improve over Liftoff and miniprot?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Here is one example of improvement over human annotation lift-over from GRCh38 to T2T-CHM13. 

    
    .. _figure-qa-scatter-plots:

    .. figure::  ../_images/human_refseq/combined_scatter_plots.png
        :align:   center
        :scale:   21 %

    Each dot represents a protein-coding transcript. If it is above the x=y line, it indicates that the LiftOn annotation possesses a higher protein sequence identity score and corresponds to a longer protein that aligns with the proteins in the reference annotation.

    In the LiftOn versus Liftoff comparison (Figure above, left), 2,075 transcripts exhibit higher protein sequence identity, with 460 achieving 100% identity. Similarly, the LiftOn versus miniprot comparison (Figure above, right) discloses better matches for 30,276 protein-coding transcripts, improving 22,616 to identical status relative to the reference. 

    In summary, LiftOn effectively corrects quite a few protein-coding transcripts during human lift-over. The improvement is even more significant when it comes to more distant species!

    Check out the :ref:`Same species lift-over section <same_species-section>`, :ref:`Closely related species lift-over section <close_species-section>`, and :ref:`Distantly related species lift-over section  <distant_species-section>` for more details.


| 



.. dropdown:: Q: How do we interpret LiftOn output?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Check out the :ref:`Output files section <output_files>`.


| 

.. dropdown:: Q: Can you explain the new protein-maximization (PM) algorithm in LiftOn?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Check out the :ref:`Protein-maximization algorithm section <protein-maximization_algorithm>`.

|

.. dropdown:: Q: How to you evaluate the lift-over annotation?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Check out the :ref:`DNA & protein transcript sequence identity score calculation section <lifton_sequence_identity>`.


| 

.. dropdown:: Q: How does LiftOn report mutated genes?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    LiftOn compares reference and target transcripts, similar to `Liftofftools <https://github.com/agshumate/LiftoffTools>`_, generating a mutation report for mapped protein-coding transcripts. 
    
    Transcripts are considered "**identical**" if their target and reference gene DNA sequences match entirely. For mutated sequences, LiftOn categorizes changes as "**synonymous**", "**non-synonymous**", "**in-frame insertion**", "**in-frame deletion**", "**frameshift**", "**stop codon gain**", "**stop codon loss**", and "**start codon loss**".

    Check out the :ref:`Mutation report section <mutation-reporting>`.

|


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