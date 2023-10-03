# Lifton

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![version](https://img.shields.io/badge/version-v.0.0.1-blue)

LiftOn is a gene annotator that relies on Liftoff and miniprot. It enhances annotations for protein-coding genes during the lift-over process and generates optimal annotations.

##  <a name="whylifton"></a>Why LiftOn❓<a class="headerlink" href="#whylifton" title="Permalink to this heading">#</a>
1. The current approach to generate the annotation of T2T-CHM13 is to run Liftoff to lift-over annotations from GRCh38 to T2T-CHM13. However, Liftoff is not perfect. T2T-CHM13 annotation is far from perfect.
2. More and more high quality assemblies are generated. We need to annotate them.
3. The current lift-over tools mainly depend on either DNA aligners (Liftoff, minimap2) or protein aligner (miniprot). They are not flawless, as there are instances where they make errors.



## <a name="whatliftondo"></a>What does LiftOn do❓<a class="headerlink" href="#whatliftondo" title="Permalink to this heading">#</a>
LiftOn takes GFF files from Liftoff and miniprot and protein FASTA files, and generate a new annotation file in GFF format. LiftOn works on the same and closely-related species. 
<!-- We also tested LiftOn by lifting-over annotations from human to mouse, and it also does pretty good job to find the optimal protein annotations. However, there are false positives or -->


* **Input**:  Liftoff GFF file  /  miniprot GFF file  /  protein FASTA file
* **Output**: LiftOn GFF file

LiftOn utilizes gene loci coordinates obtained from Liftoff, as Liftoff employs an overlapping fixing algorithm to determine the most suitable gene locus for each gene.

First, LiftOn extracts protein sequences annotated by Liftoff and miniprot, and aligns them with the reference proteins.

Next, LiftOn employs an algorithm that compares each section of the protein alignments from Liftoff and miniprot, corrects errors in exon and CDS boundaries, and produces the optimal protein annotations.

## <a name="whosplaminterested"></a>Who is it for❓<a class="headerlink" href="#whosplaminterested" title="Permalink to this heading">#</a>
1. If you have sequenced and assembled a new human genome and need to annotate it, LiftOn is the ideal choice for generating annotations.
2. If you wish to utilize the finest CHM13 annotation, you can run LiftOn! We have also pre-generated the [T2T_CHM13_LiftOn.gff3](https://khchao.com) file for your convenience."

## <a name="installation"></a>Installation<a class="headerlink" href="#installation" title="Permalink to this heading">#</a>



