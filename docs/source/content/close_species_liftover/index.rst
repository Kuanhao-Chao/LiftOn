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


.. raw:: html

    <script type="text/javascript">
        document.addEventListener("DOMContentLoaded", mutation_lvl_1_fuc);
        var observer = new MutationObserver(mutation_lvl_1_fuc)
        observer.observe(document.body, {attributes: true, attributeFilter: ['data-theme']});
        console.log(document.body);
    </script>

|


.. _close_species-section:

Closely related species lift-over
====================================

.. When you want to annotate your newly assembled genome, and you realized that there is no same species reference annotation that is previously annotated, or there is a closely related species that is way much more anotated than your assembled species. In this case, you can run LiftOn to lift-over annotations between two closely-related speceis. 

.. LiftOn excels in performing lift-over on annotations between two closely-related species. Here ase an example of running LiftOn to lift-over annotations from *Homo sapiens* GRCh38 to *Pan troglodytes* (chimpanzee).


When you intend to annotate your recently assembled genome and find that there is no pre-annotated reference for the same species or that a closely related species has significantly better annotations than your assembled species, you can use LiftOn to transfer annotations between the two species.

LiftOn is particularly effective in mapping annotation between closely-related species. Here is an example of using LiftOn to transfer annotations from *Homo sapiens* GRCh38 to *Pan troglodytes* (chimpanzee).

|


.. admonition:: LiftOn examples
    :class: note

    * :ref:`close_species_liftover_human_to_chimp`

.. toctree::
    :hidden:
    
    liftover_GRCh38_2_chimp
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