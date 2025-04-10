Retrieving predictions
======================

Oktoberfest relies on retrieving predictions from a `Koina <https://koina.wilhelmlab.org/>`_ server that hosts specific models for peptide property prediction. Users can use any publicly available community server or host their own server.

Connecting to a community server
--------------------------------

Our publicly available community server is available at `koina.wilhelmlab.org:443`.
If you want to connect to it, you need to have the following flags in your config file (default settings):

.. code-block:: json

   {
      "prediction_server": "koina.wilhelmlab.org:443",
      "ssl": true,
   }

Once more community servers become available, we will add a list here.

Currently supported models
~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the list of currently supported and tested models for Oktoberfest provided by our community server.

.. table::
   :class: fixed-table

   +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Intensity models                                                                                                 |                             Description                                                                                                                                                      |
   +==================================================================================================================+==============================================================================================================================================================================================+
   | Prosit_2019_intensity                                                                                            | Developed for HCD tryptic peptides only. We recommend using the Prosit_2020_intensity_HCD model instead, since it showed slightly superior performance on tryptic peptides as well.          |
   +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | `Prosit_2020_intensity_HCD <https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_HCD/infer>`_           | Developed for HCD tryptic and non-tryptic peptides. Supported modifications are oxidation and carbamidomethylation. Latest version we recommend to use for HCD.                              |
   +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | `Prosit_2020_intensity_CID <https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_CID/infer>`_           | Developed for CID tryptic and non-tryptic peptides. Supported modifications are oxidation and carbamidomethylation. Latest version we recommend to use for CID.                              |
   +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | `Prosit_2020_intensity_TMT <https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_TMT/infer>`_           | Developed for HCD and CID, tryptic and non-tryptic peptides. Latest version we commend for TMT labeled peptides in general.                                                                  |
   +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | `Prosit_2023_intensity_timsTOF <https://koina.wilhelmlab.org/docs#post-/Prosit_2023_intensity_timsTOF/infer>`_   | Developed for timsTOF, tryptic and non-tryptic peptides. Latest version we commend to use for timsTOF.                                                                                       |
   +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   |  `AlphaPept_ms2_generic <https://koina.wilhelmlab.org/docs#post-/AlphaPept_ms2_generic/infer>`_                  | Developed for generic data support, including TMT, timsTOF and various instrument types.                                                                                                     |
   +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   |  `ms2pip_2021_HCD <https://koina.wilhelmlab.org/docs#post-/ms2pip_2021_HCD/infer>`_                              | Developed for HCD tryptic and non-tryptic peptides.                                                                                                                                          |
   +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. table::
   :class: fixed-table

   +-----------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
   | iRT models                                                                                    |                             Description                                                                                   |
   +===============================================================================================+===========================================================================================================================+
   | `Prosit_2019_irt <https://koina.wilhelmlab.org/docs#post-/Prosit_2019_irt/infer>`_            | While developed for tryptic peptides only, we did not observe a drop in prediction performance for non-tryptic peptides.  |
   +-----------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
   | `Prosit_2020_irt_TMT <https://koina.wilhelmlab.org/docs/#post-/Prosit_2020_irt_TMT/infer>`_   | Developed for TMT labeled peptides.                                                                                       |
   +-----------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
   | `AlphaPept_rt_generic <https://koina.wilhelmlab.org/docs#post-/AlphaPept_rt_generic/infer>`_  | Developed for for generic data support, including TMT, timsTOF and various instrument types.                              |
   +-----------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+

Once support for additional models is implemented in Oktoberfest, they will be added here.

Hosting and adding your own models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case you are planning to host your own private or public instance of Koina or want us to host your model, please refer to the official `Koina documentation <https://koina.wilhelmlab.org/docs#overview>`_.

Local intensity prediction & refinement learning
------------------------------------------------

Instead of using a pre-trained intensity predictor via Koina, you can also predict intensity locally from a `DLomix <https://github.com/wilhelm-lab/dlomix>`_ model.
If you do not have a pre-trained intensity predictor in `.keras` format, a baseline model will automatically be downloaded for you.
Importantly, this also gives you the option to refinement-learn the pre-trained predictor on your input dataset and using the refined model for intensity prediction.

For local intensity prediction and refinement learning, you need to provide either a path to a pre-trained model or the keyword `baseline`
(for the runner to automatically download a model for you) as the intensity model in the config, and specify `localPredictionOptions` as well as optional `refinementLearningOptions`.
For more details, refer to the :ref:`job <jobs:b) with refinement>` and :ref:`configuration <config:Applicable to local intensity prediction>` docs.
