.. _behind-the-scenes-splam:

Behind the scenes
=================================================

|

.. _data-curation:

Matching miniprot & Liftoff genome annotation
+++++++++++++++++++++++++++++++++++++++++++++++


.. _splam-data-curation:
.. figure::  ../_images/splam_data_curation.png
    :align:   center
    :scale:   21 %

    Illustration of the 4 types of splice sites used for training and testing Splam: Positive-MANE, Positive-ALT, Negative-1, and Negative-Random. Positive-MANE sites (blue) are selected from the MANE database and supported by at least 100 alignments, while Positive-ALT (green) are present in the RefSeq database but missing from MANE, and also supported by at least 100 alignments. Negative-1 sites (orange) occur on the opposite strand of a known gene and are supported by only 1 alignment, and Negative Random sites (red) are random GT-AG pairs on the opposite strand that do not overlap with any known splice sites and have no alignment support. 

|


.. _model-architecture:

Chaining algorithm
+++++++++++++++++++++++++++++++++++

Splam utilizes a deep dilated residual convolutional neural network (CNN) that incorporates grouped convolution layers within the residual units. 

|

Mutation reporting
+++++++++++++++++++++++++++++++++++

Open-reading-frame search
+++++++++++++++++++++++++++++++++++


DNA & protein transcript sequence identity score calculation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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