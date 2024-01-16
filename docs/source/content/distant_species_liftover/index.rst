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

.. _distant_species-section:

Distant species lift-over
====================================

.. The fantastic tool, `Liftoff <https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=true>`_, does not work well mapping annotations between distant related species due to the fact that it only depends on the information of DNA aligner.  LiftOn incorporates the information from the protein aligner as well and therefore improves the annotation lift-over process, especially for mapping between distant species! Following are two examples:

The excellent tool, `Liftoff <https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=true>`_, faces limitations in effectively mapping annotations between distantly related species as it relies solely on DNA aligner information. In contrast, LiftOn integrates data from both DNA and protein aligners, enhancing the annotation lift-over process, particularly for mapping between distant species!

Following are two examples:







.. Resource Utilization: Annotating a genome is a resource-intensive process that involves experimental and computational efforts. If annotations for a closely related species already exist, lifting over those annotations can save time and resources compared to annotating the genome of a distant species from scratch.

.. Incomplete Genomic Data: In some cases, the genomic data for a particular species may be incomplete or of lower quality. Lifting over annotations from a well-annotated species can help fill in the gaps and improve the overall understanding of the genomic landscape of the less-studied species.

.. Biomedical Research: In the context of biomedical research, understanding the genomic similarities and differences between species is crucial for studying diseases and developing treatments. Annotations lift-over can contribute to identifying orthologous genes and functional elements relevant to human health.

|

.. toctree::
    :maxdepth: 1

    liftover_drosophila_erecta
    liftover_mouse_2_rat


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