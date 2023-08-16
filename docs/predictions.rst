Retrieving predictions
======================

Oktoberfest relies on retrieving predictions from a `Koina <https://koina.proteomicsdb.org/>`_ server that hosts specific models for peptide property prediction. Users can use any publicly available community server or host their own server.

Connecting to a community server
--------------------------------

Our publicly available community server is available at `koina.proteomicsdb.org:443`.
If you want to connect to it, you need to have the following flags in your config file (default settings):

.. code-block:: json

   {
      "prediction_server": "koina.proteomicsdb.org:443",
      "ssl": true,
   }

Once more community servers become available, we will add a list here.

Currently supported models
--------------------------

This is the list of currently supported and tested models for Oktoberfest provided by our community server.

.. table::
   :class: fixed-table

   +----------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Intensity models                                                                                         |                             Description                                                                                                                                                      |
   +==========================================================================================================+==============================================================================================================================================================================================+
   | Prosit_2019_intensity                                                                                    | Developed for HCD tryptic peptides only. We recommend using the Prosit_2020_intensity_HCD model instead, since it showed slightly superior performance on tryptic peptides as well.          |
   +----------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | `Prosit_2020_intensity_HCD <https://koina.proteomicsdb.org/docs#post-/Prosit_2020_intensity_HCD/infer>`_ | Developed for HCD tryptic and non-tryptic peptides. Supported modifications are oxidation and carbamidomethylation. Latest version we recommend to use for HCD.                              |
   +----------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | `Prosit_2020_intensity_CID <https://koina.proteomicsdb.org/docs#post-/Prosit_2020_intensity_CID/infer>`_ | Developed for CID tryptic and non-tryptic peptides. Supported modifications are oxidation and carbamidomethylation. Latest version we recommend to use for CID.                              |
   +----------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | `Prosit_2020_intensity_TMT <https://koina.proteomicsdb.org/docs#post-/Prosit_2020_intensity_TMT/infer>`_ | Developed for HCD and CID, tryptic and non-tryptic peptides. Latest version we commend for TMT labeled peptides in general.                                                                  |
   +----------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. table::
   :class: fixed-table

   +-----------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
   | iRT models                                                                                    |                             Description                                                                                   |
   +===============================================================================================+===========================================================================================================================+
   | `Prosit_2019_irt <https://koina.proteomicsdb.org/docs#post-/Prosit_2019_irt/infer>`_          | While developed for tryptic peptides only, we did not observe a drop in prediction performance for non-tryptic peptides.  |
   +-----------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
   | `Prosit_2020_irt_TMT <https://koina.proteomicsdb.org/docs/#post-/Prosit_2020_irt_TMT/infer>`_ | Developed for TMT labeled peptides.                                                                                       |
   +-----------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+

Once support for additional models is implemented in Oktoberfest, they will be added here.

Hosting and adding your own models
----------------------------------

In case you are planning to host your own private or public instance of Koina or want us to host your model, please refer to the official `Koina documentation <https://koina.proteomicsdb.org/docs#overview>`_.