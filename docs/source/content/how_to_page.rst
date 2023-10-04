.. _Q&A:

Q & A ...
==========

.. Q: What is Splam?
.. -------------------------------------------

.. <div style="padding-left:20px">

.. dropdown:: Q: What is Splam?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-dark font-weight-bolder howtoclass
    :body: bg-dark text-left

    Splam stands for two things: **(1)** Splam refers to the deep grouped residual CNN model that we designed to accurately predict splice junctions (based solely on an input DNA sequence), and **(2)** it also stands for this software which can clean up alignment files and evaluate annotation files.

|


.. Q: Why do we need Splam?
.. -------------------------------------------

.. dropdown:: Q: Why do we need Splam?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    We are concerned about the method of training splice junction predictors by relying on splice junctions in solely canonical transcripts. Designing a splice site recognition method based only on one isoform per gene may result in mislabeling alternative splice sites even when they are perfectly valid. Therefore, 

        * **we designed a biologically realistic model.** Splam was trained on combined donor and acceptor pairs, with a focus on a narrow window of 400 base pairs surrounding each splice site. This approach is inspired by the understanding that the splicing process primarily relies on signals within this specific region.


    Furthermore, there are two applications of Splam: 

    When inspecting an alignment file in IGV, it becomes apparent that some reads are spliced and aligned across different gene loci or intergenic regions. This raises the question, "Are these spliced alignments correct?" Therefore,

        * **we need a trustworthy way to evaluate all the spliced alignments in the alignment file.** Splam learns splice junction patterns, and we have demonstrated that applying Splam to remove spurious spliced alignments improves transcript assembly! :ref:`alignment evaluation section <alignment-detailed-section>`.

    Additionally, we acknowledge that annotation files are not perfect, and there are more errors in the assembled transcripts. The current approach to assessing assembled transcripts involves comparing them with the annotation.

        * **we can utilize Splam to score all introns in transcripts and provide a reference-free evalutation.**  :ref:`annotation evaluation section <annotation-detailed-section>`.



|

.. Q: What makes Splam different from SpliceAI?
.. -------------------------------------------

.. dropdown:: Q: What makes Splam different from SpliceAI?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left


    Splam and SpliceAI are both frameworks used for predicting splice junctions in DNA sequences, but they have some key differences.


    #. **Input constraints:**
 
       * **Splam**: Follows the design principle of using biologically realistic input constraints. It uses a window limited to 200 base pairs on each side of the donor and acceptor sites, totaling 800 base pairs. Furthermore, we pair each donor and acceptor as follows

       .. figure::  ../_images/splam_input.png
            :align:   center
            :scale:   7%
     
       * **SpliceAI**: The previous state-of-the-art CNN-based system, SpliceAI, relies on a window of 10,000 base pairs flanking each splice site to obtain maximal accuracy. However, this window size is much larger than what the splicing machinery in cells can recognize.


    #. **Training data**
    
       * **Splam**: Was trained using a high-quality dataset of human donor and acceptor sites. Check out the :ref:`data curation section <data-curation>`.
    
       * **SpliceAI**: Was trained with canonical transcripts only, and does not consider alternative splicing.

| 

.. dropdown:: Q: What is the difference between two released model, :code:`splam.pt` and :code:`splam_script.pt`?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    You may have noticed that we have two released Splam models: ":code:`splam.pt`" and ":code:`splam_script.pt`".

    * :code:`splam.pt` [`link <https://github.com/Kuanhao-Chao/splam/blob/main/model/splam.pt>`_] is the original model that requires the original model script to load and run.

    * :code:`splam_script.pt` [`link <https://github.com/Kuanhao-Chao/splam/blob/main/model/splam_script.pt>`_] is the Torchscripted Splam model. Torchscript serializes and optimizes PyTorch code for improved performance and deployment. Essentially, it allows you to convert PyTorch code into a more efficient intermediate representation, which can be used for Just-In-Time (JIT) compilation and deployment without the need for the Python interpreter.

    .. important::

        In sum, we strongly recommend using :code:`splam_script.pt` for all users. It provides a faster, portable, and secure way of deploying the model.



|

.. Q: Which mode should I run Splam, :code:`cpu`, :code:`cuda`, or :code:`mps`?
.. -------------------------------------------------------------------------------

.. dropdown:: Q: Which mode should I run Splam, :code:`cpu`, :code:`cuda`, or :code:`mps`?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left


    By default, Splam automatically detects your environment and runs in :code:`cuda` mode if CUDA is available. However, if your computer is running macOS, Splam will check if :code:`mps` mode is available. If neither :code:`cuda` nor :code:`mps` are available, Splam will run in :code:`cpu` mode. You can explicitly specify the mode using the :code:`-d / --device` argument.

    .. important::

        In sum, 

        1. if you are using the Apple Silicon Mac, you should run Splam with :code:`mps` mode. 


        2. If you are using Linux with CUDA installed, you should run Splam with :code:`cuda` mode.


        3. If you are none of the above cases, then you can still run Splam with :code:`cpu`` mode.


    You can check out the `Pytorch website <https://pytorch.org/docs/stable/tensor_attributes.html#torch.device>`_ for more explanation about the :code:`device` parameter.


| 

.. Q: How do I interpret Splam scores?
.. -------------------------------------

.. dropdown:: Q: How do I interpret Splam scores?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Given an input of length 800nt, Splam outputs a Tensor with dimensions (3 x 800). The first channel represents the "acceptor scores", the second channel represents the "donor scores", and the third channel represents the "non-splice site scores". Each score is between 0 and 1, representing Splam's confidence in a given site being a splice site. A score closer to one indicates a higher level of confidence in its classification.

|

.. .. Q: What is canonical transcripts? 
.. .. ------------------------------------------

.. .. dropdown:: Q: What is canonical transcripts? 
..     :animate: fade-in-slide-down
..     :container: + shadow
..     :title: bg-light font-weight-bolder
..     :body: bg-light text-left


.. |

.. .. Q: What is alternative splicing?
.. .. ------------------------------------------

.. .. dropdown:: Q: What is alternative splicing?
..     :animate: fade-in-slide-down
..     :container: + shadow
..     :title: bg-light font-weight-bolder
..     :body: bg-light text-left


.. Q: What is the model architecture of Splam?
.. -----------------------------------------


.. dropdown:: Q: What is the model architecture of Splam?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Check out the :ref:`model architecture section <model-architecture>`.

| 

.. Q: How is Splam trained?
.. --------------------------------

.. dropdown:: Q: How is Splam trained?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Check out the :ref:`splam training and testing section <splam-train-test>`.


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