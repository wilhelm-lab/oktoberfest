Retrieving predictions
======================

Oktoberfest relies on retrieving predictions from a `Koina <https://koina.wilhelmlab.org/>`_ server that hosts specific models for peptide property prediction. Users can use any publicly available community server or host their own server.

Connecting to a community server
--------------------------------

.. list-table:: Currently supported intensity models
   :class: fixed-table
   :widths: 30 70
   :header-rows: 1

   * - Intensity model
     - Description
   * - Prosit_2019_intensity
     - Developed for HCD tryptic peptides only. Use Prosit_2020_intensity_HCD instead for better performance.
   * - `Prosit_2020_intensity_HCD <https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_HCD/infer>`_
     - Developed for HCD tryptic and non-tryptic peptides. Supports oxidation and carbamidomethylation. Recommended version for HCD.
   * - `Prosit_2020_intensity_CID <https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_CID/infer>`_
     - Developed for CID tryptic and non-tryptic peptides. Supports oxidation and carbamidomethylation. Recommended for CID.
   * - `Prosit_2020_intensity_TMT <https://koina.wilhelmlab.org/docs#post-/Prosit_2020_intensity_TMT/infer>`_
     - For HCD and CID, tryptic and non-tryptic peptides. Recommended for TMT-labeled peptides.
   * - `Prosit_2023_intensity_timsTOF <https://koina.wilhelmlab.org/docs#post-/Prosit_2023_intensity_timsTOF/infer>`_
     - Developed for timsTOF. Works for both tryptic and non-tryptic peptides.
   * - `Prosit_2023_intensity_XL_CMS2 <https://koina.wilhelmlab.org/docs#post-/Prosit_2023_intensity_XL_CMS2/infer>`_
     - For HCD cleavable cross-linked peptides (DSSO, DSBU). Supports oxidation and carbamidomethylation.
   * - `Prosit_2024_intensity_XL_NMS2 <https://koina.wilhelmlab.org/docs#post-/Prosit_2024_intensity_XL_NMS2/infer>`_
     - For HCD non-cleavable cross-linked peptides (DSS, BS3). Supports oxidation and carbamidomethylation.
   * - `AlphaPept_ms2_generic <https://koina.wilhelmlab.org/docs#post-/AlphaPept_ms2_generic/infer>`_
     - For generic MS2 data, including TMT, timsTOF, and multiple instrument types.
   * - `ms2pip_2021_HCD <https://koina.wilhelmlab.org/docs#post-/ms2pip_2021_HCD/infer>`_
     - For HCD tryptic and non-tryptic peptides.


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
