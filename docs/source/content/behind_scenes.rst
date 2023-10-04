.. _behind-the-scenes-splam:

Behind the scenes
=================================================

|

.. _data-curation:

Data curation
+++++++++++++++++++++++++++++++++++

Data plays a critically important role in deep learning. The model performance largely depends on the quality, quantity, and relevance of the data used for training. 

To understand how Splam works, you first need to understand how the training data was generated. For curating input splice sites, our approach differs from the state-of-the-art tool, SpliceAI, which was trained on splice sites labeled on canonical transcripts only. Our approach focuses on a smaller region that encompasses the donor and acceptor sites. 

We extract a 400bp sequence centered around the splice sites, shifting the focus from learning the splicing pattern at the transcript level to the splice site level. This strategy enables localized prediction of splicing events, mitigating biases associated with canonical transcripts and eliminating the need for larger window sizes encompassing the entire intron and 10,000bp flanking sequences. In summary, each data point in our Splam model consists of a donor-acceptor pair with a 400bp window size for each, resulting in a total of 800bp.

|

Model input & output dimension
--------------------------------

The input to the model is a 800bp one-hot-encoded DNA sequence in the dimension of :code:`4 x 800`. The output is a :code:`3 x 800` array, with every column representing a base pair, and each row reflecting its probability of being a donor site, acceptor site, or neither (as shown in :numref:`splam-input`). The three scores sum up to 1. For the development of Splam, we used the Pytorch version 1.13.0 framework.


.. _splam-input:
.. figure::  ../_images/splam_input.png
    :align:   center
    :scale:   7 %

    Splam takes 400bp regions around a donor and acceptor site, resulting in a total of 800bp. The input sequences are one-hot encoded into vectors, where :math:`A` is :math:`[1, 0, 0, 0]`, :math:`C` is :math:`[0, 1, 0, 0]`, :math:`G` is :math:`[0, 0, 1, 0]`, :math:`T` is :math:`[0, 0, 0, 1]` and :math:`N` is :math:`[0, 0, 0, 0]`,. For the labels, a donor site is represented as :math:`[0, 0, 1]`, an acceptor site as :math:`[0, 1, 0]`, and an non-junction site as :math:`[1, 0, 0]`.

|

Positive data
---------------

We filtered the junctions by extracting only those that were supported by more than 100 alignments in the `Tiebrush <https://academic.oup.com/bioinformatics/article/37/20/3650/6272575>`_ merged BAM file. Furthermore, we used RefSeq GRCh38 patch 14 annotations to extract all protein-coding genes and intersected all 100-alignment-supported splice sites with the protein-coding splice sites observed in the RefSeq database.

The resulting splice sites were then categorized as **Positive-MANE** (blue in :numref:`splam-data-curation`) if they were in the `Ensembl GRCh38 MANE release 1.0 file <https://ftp.ncbi.nih.gov/refseq/MANE/MANE_human/release_1.0/>`_, and **Positive-ALT** (green in :numref:`splam-data-curation`) if they were only in `RefSeq GRCh38 patch 14 <https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/110/GCF_000001405.40_GRCh38.p14/>`_ but not in MANE.

This approach provides expression evidence to support splice junctions and thus eliminates transcription noise, increasing our confidence in the extracted donor and acceptor sites. Overall, we identified 180,195 splice sites in Positive MANE and 85,908 splice sites in Positive-ALT

|

Negative data
---------------

Curating high-quality negative splice junctions is a challenging task due to incomplete annotations, unknown functions, and potential transcriptional noise. Labeling unannotated junctions as negatives without caution may misclassify some spliced junctions. Various methods can be used to generate negative splice junctions, such as selecting random dinucleotide pairs or random GT-AG pairs from the genome. However, these artificial junctions may differ significantly from true splice sites, leading to the learning of non-critical patterns by the model, resulting in low sensitivity and precision.

To address this issue in training Splam, two novel approaches were adopted for generating challenging negative splice junctions. **(1)** The first approach involved selecting random GT-AG pairs on the opposite strand of protein-coding gene loci. Since overlapping genes are rare in eukaryotes, it is unlikely to have another transcribed gene on the opposite strand of a protein-coding gene. This resulted in 4,467,910 splice junctions referred to as **Negative-Random** (red in :numref:`splam-data-curation`).

To further increase the difficulty of negative junctions, **(2)** only splice junctions with 1-alignment support on the opposite strand of protein-coding gene loci were chosen. These junctions are likely to be transcriptional noise or alignment artifacts and should be eliminated, making them more suitable for Splam training. This approach generated 2,486,305 splice junctions referred to as **Negative-1** (orange in :numref:`splam-data-curation`). Accurately curating these challenging negative splice junctions allowed for effective training of Splam.


.. _splam-data-curation:
.. figure::  ../_images/splam_data_curation.png
    :align:   center
    :scale:   21 %

    Illustration of the 4 types of splice sites used for training and testing Splam: Positive-MANE, Positive-ALT, Negative-1, and Negative-Random. Positive-MANE sites (blue) are selected from the MANE database and supported by at least 100 alignments, while Positive-ALT (green) are present in the RefSeq database but missing from MANE, and also supported by at least 100 alignments. Negative-1 sites (orange) occur on the opposite strand of a known gene and are supported by only 1 alignment, and Negative Random sites (red) are random GT-AG pairs on the opposite strand that do not overlap with any known splice sites and have no alignment support. 

|


.. _model-architecture:

Model architecture
+++++++++++++++++++++++++++++++++++

Splam utilizes a deep dilated residual convolutional neural network (CNN) that incorporates grouped convolution layers within the residual units. 

|

Residual unit
---------------

Splam architecture consists of 20 residual units, each containing two convolutional layers. The model uses a grouped convolution approach with a parameter called :code:`group` set to 4. The hyperparameters of Splam include :code:`F` (number of filters), :code:`W` (window size), :code:`D` (dilation rate), and :code:`G` (groups), which are shown as (:code:`F`, :code:`W`, :code:`D`, :code:`G`) in :numref:`splam-model` (b). The concept of grouped convolution, which allows for memory saving with minimal accuracy loss, is inspired by the ResNext model.


.. For a convolutional layer, nucleotide in the sequence, it checks the region of :code:`F x (W-1)`, and therefore :code:`2F x (W-1)` neighboring positions for a residual unit. Furthermore, in a grouped convolution with :code:`G` groups, :code:`F/G` filters are applied to each :code:`F/G` of the input for a $G$× reduction in parameters used. In total, there are 651,715 parameters in Splam. 

Each convolutional layer in the residual unit follows a batch normalization and a rectified linear unit (ReLU) :numref:`splam-model` (b), and the input of the unit is residually connected to its output. He et al. :cite:p:`he2016deep` introduced residual units to address the issue of training accuracy degradation in deep learning. The inclusion of shortcut connections enables successful training of deeper models using simple stochastic gradient descent (SGD) with backpropagation, leading to improved accuracy as the depth increases.

|

Residual group
---------------

A group of four residual units forms a bigger residual group, and 20 RUs are clustered into five residual groups. Residual groups are stacked such that the output of the i \ :sup:`th`\  residual group is connected to the i+1 \ :sup:`th`\  residual group. Furthermore, the output of each residual group undergoes a convolutional layer with the parameters :code:`(64, 1, 1)`, and then being added to all the other outputs of residual groups (residual connections colored in red), which then is passed into the last convolutional layer in :code:`(3, 1, 1)` and a softmax layer. :code:`F` is set to :code:`64` for all convolutional layers, and for each residual group, :code:`W` is set to :code:`11`, :code:`11`, :code:`11`, :code:`21`, and :code:`21`, and :code:`D` is set to :code:`1`, :code:`5`, :code:`10`, :code:`15`, and :code:`20` in residual groups in sequence. :code:`G` is by default is :code:`1` for all convolutional layers, but set to :code:`4` in the residual units. We visualized the architecture of Splam in :numref:`splam-model`. For each nucleotide position, its total neighboring span of the Splam model is 

.. math::

    S=\sum_{i=1}^{20}2F_{i}\times(W_{i}-1)

.. _splam-model:
.. figure::  ../_images/splam_model_architecture.png
    :align:   center
    :scale:   18 %

    The model architecture of the Splam deep residual convolutional neural network.

|


.. _splam-train-test:

Splam training & testing
+++++++++++++++++++++++++++++++++++

After curating the gold standard dataset, we divided all splice junctions into two datasets: one for model training and the other for testing. For model training, we utilized all the splice sites on the main chromosomes, except chromosomes 1 and 9. For model testing, we used the splice sites on the held-out chromosomes 1 and 9, with the splice sites in paralogs removed.

|

Hyperparameters
-----------------

To train Splam, we used a batch size of 100 and trained it for 15 epochs. We employed the AdamW optimizer with the default learning rate of 0.03. A 1000-step warmup was utilized, with the learning rate increasing linearly from 0 to 0.03. The learning rate then decreased following the values of the cosine function between 0.03 to 0 (:numref:`train-lr`).

.. _train-lr:
.. figure::  ../_images/train_lr.png
    :align:   center
    :scale:   80 %

    The learning rate for each Splam update during training

|

Loss function
---------------

We further improved Splam's performance by changing the loss function. Instead of using the commonly used cross entropy (:ref:`Eq. 1 <equation-cel>`), we replaced it with focal loss (:ref:`Eq. 2 <equation-fl>`) :cite:p:`lin2017focal`.

.. _equation-cel:
.. math::
    :label: 1

    Loss_{CEL} = \sum_{class \in \{donor, acceptor, neither\}} I_{class}\times\log(P_{class}) 

.. _equation-fl:
.. math::
    :label: 2

    Loss_{FL} = \sum_{class \in \{donor, acceptor, neither\}} I_{class}\times (1-P_{class})^{\gamma} \times \log(P_{class}), \text{where } \gamma = 2


Focal loss puts more emphasis on the challenging data points where Splam is more likely to make incorrect predictions and penalized these data points by an additional :math:`(1-P)^{\gamma}`` scale, where :math:`\gamma = 2` and :math:`P` is the probability of each class. This scale quantifies the degree of inaccuracy in predictions, instead of simply binary misclassifications that cross entropy applies.

|

.. Training precision
.. --------------------

.. .. figure::  ../_images/train_J_threshold_precision.png
..     :align:   center
..     :scale:   80 %



.. Testing precision
.. --------------------
.. .. figure::  ../_images/test_J_threshold_precision.png
..     :align:   center
..     :scale:   80 %


|

Reference
+++++++++++++++++++++++++++++++++++

.. bibliography::


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