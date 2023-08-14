Retrieving predictions
======================

Oktoberfest relies on retrieving predictions from a `koina <https://koina.proteomicsdb.org/>`_ or any other community server that hosts specific models for peptide property prediction. server that hosts supported models for peptide property predictions. Users can use any publicly available community server or host their own server.

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

   +----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Intensity models           |                             Description                                                                                                                                                      |
   +============================+==============================================================================================================================================================================================+
   | Prosit_2019_intensity      | deprecated, please use the 2020 model                                                                                                                                                        |
   +----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Prosit_2020_intensity_HCD  | your go to model for fragment intensity prediction for HCD fragmentation, find out more about this model `here <https://koina.proteomicsdb.org/docs#post-/Prosit_2020_intensity_HCD/infer>`_ |
   +----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Prosit_2020_intensity_CID  | your go to model for fragment intensity prediction for CID fragmentation, find out more about this model `here <https://koina.proteomicsdb.org/docs#post-/Prosit_2020_intensity_CID/infer>`_ |
   +----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Prosit_2020_intensity_TMT  | your go to model for fragment intensity prediction for TMT, find out more about this model `here <https://koina.proteomicsdb.org/docs/#post-/Prosit_2020_intensity_TMT/infer>`_              |
   +----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. table::
   :class: fixed-table

   +----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | iRT models\ \              |                             Description                                                                                                                                                      |
   +============================+==============================================================================================================================================================================================+
   | Prosit_2019_irt            | all purpose model for retention time prediction, find out more about this model `here <https://koina.proteomicsdb.org/docs/#post-/Prosit_2019_irt/infer>`_                                   |
   +----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Prosit_2020_irt_TMT        | your go to model for retention time prediction for TMT, find out more about this model `here <https://koina.proteomicsdb.org/docs/#post-/Prosit_2020_irt_TMT/infer>`_                        |
   +----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Once support for additional models is implemented in Oktoberfest, they will be added here.

Hosting and adding your own models
----------------------------------

In case you are planning to host your own private or public instance of koina or want us to host your model, please refer to the official `koina documentation <https://koina.proteomicsdb.org/docs#overview>`_.


