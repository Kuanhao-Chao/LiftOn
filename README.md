# LiftOn

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![version](https://img.shields.io/badge/version-v.0.0.1-blue)

LiftOn is a lift-over annotator that takes [Liftoff](https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=true) and [miniprot](https://academic.oup.com/bioinformatics/article/39/1/btad014/6989621) GFF files as input. It accurately generates gene annotations, with a particular focus on protein-coding genes. LiftOn takes consensus from both sources and generates optimal annotations that outperform both [Liftoff](https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=true) and [miniprot](https://academic.oup.com/bioinformatics/article/39/1/btad014/6989621)!




<!-- LiftOn improves the annotations of protein-coding genes provided by Liftoff. Additionally, during the lift-over process, it identifies any extra copies of protein-coding genes and generates optimal annotations.
 -->

<br>

##  <a name="whylifton"></a>Why LiftOn❓<a class="headerlink" href="#whylifton" title="Permalink to this heading">#</a>
1. The current approach to generate the annotation of [T2T-CHM13](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/) is to run Liftoff to lift-over annotations from GRCh38 to T2T-CHM13. However, Liftoff is not perfect. T2T-CHM13 annotation is far from perfect. **We need a tool to accurately generates T2T-CHM13 annotations**.
2. More and more high quality assemblies are generated. We need to annotate them.
3. The current lift-over tools mainly depend on either DNA aligners ([Liftoff](https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=true), [minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778)) or protein aligner ([miniprot](https://academic.oup.com/bioinformatics/article/39/1/btad014/6989621)). They are not perfect, as there are instances where they make mistakes.

<br>

## <a name="whatliftondo"></a>What does LiftOn do❓<a class="headerlink" href="#whatliftondo" title="Permalink to this heading">#</a>
LiftOn takes GFF files from Liftoff and miniprot and reference protein sequences in a FASTA file, and generates a new annotation file in GFF format. LiftOn works on the same and closely-related species. 
<!-- We also tested LiftOn by lifting-over annotations from human to mouse, and it also does pretty good job to find the optimal protein annotations. However, there are false positives or -->


* **Input**:  Liftoff GFF file  /  miniprot GFF file  /  protein FASTA file
* **Output**: LiftOn GFF file

LiftOn utilizes gene loci coordinates obtained from Liftoff, as Liftoff employs an overlapping fixing algorithm to determine the most suitable gene locus for each gene.

First, LiftOn extracts protein sequences annotated by Liftoff and miniprot, and aligns them to the reference proteins.

Next, LiftOn employs an algorithm that compares each section of the protein alignments from Liftoff and miniprot, corrects errors in exon and CDS boundaries, and produces the optimal protein annotations.

<br>

## <a name="whosplaminterested"></a>Who is it for❓<a class="headerlink" href="#whosplaminterested" title="Permalink to this heading">#</a>
1. If you have sequenced and assembled a new human genome and need to annotate it, LiftOn is the ideal choice for generating annotations.
2. If you wish to utilize the finest CHM13 annotation, you can run LiftOn! We have also pre-generated the [T2T_CHM13_LiftOn.gff3](https://khchao.com) file for your convenience."

<br>

<!-- ## <a name="itworks"></a>It works❗️<a class="headerlink" href="#itworks" title="Permalink to this heading">#</a> -->

<!-- <br> -->

## <a name="installation"></a>Installation<a class="headerlink" href="#installation" title="Permalink to this heading">#</a>

LiftOn is on [PyPi](https://pypi.org/). This is the easiest installation approach. Check out all the releases here.

```
$ pip install lifton
```

You can also install LiftOn from source

```
$ git clone https://github.com/Kuanhao-Chao/LiftOn --recursive

$ cd lifton

$ python setup.py install

```

<br>

## <a name="quick_start"></a>Quick Start<a class="headerlink" href="#quick-start" title="Permalink to this heading">#</a>

Running LiftOn is simple. It only requires one line of code!

<!-- See these examples on Google Colab: [![](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Kuanhao-Chao/splam/blob/main/notebook/splam_example.ipynb)

 -->

### Example 1: clean up alignment files (`BAM`)

``` bash
$ cd test

# Step 1: extract splice junctions in the alignment file
$ lifton --proteins protein.fasta --liftoffdb CHM13_MANE.sort.gff3_db --miniprotdb CHM13_MANE_miniprot.fix.sorted.gff_db -o chm13v2.0.fa
```

<br>

## <a name="citation"></a>Citation<a class="headerlink" href="#citation" title="Permalink to this heading">#</a>


Kuan-Hao Chao*, Mihaela Pertea, Steven L Salzberg*, "LiftOn: a tool to improve annotations for protein-coding genes during the lift-over process.", <i>bioRxiv</i> <b>2023.07.27.550754</b>, doi: [https://doi.org/10.1101/2023.07.27.550754](https://doi.org/10.1101/2023.07.27.550754), 2023
