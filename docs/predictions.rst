Retrieving predictions
======================

Oktoberfest relies on retrieving predictions from a `Koina <https://koina.wilhelmlab.org/>`_ server that hosts specific models for peptide property prediction. Users can use any publicly available community server or host their own server.

Connecting to a community server
--------------------------------

.. table:: Currently supported intensity models
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
    | `Prosit_2023_intensity_XL_CMS2 <https://koina.wilhelmlab.org/docs#post-/Prosit_2023_intensity_XL_CMS2/infer>`_   | Developed for HCD cleavable cross-linked peptides (DSSO, DSBU). Supports oxidation and carbamidomethylation.                                                                                 |
    +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | `Prosit_2024_intensity_XL_NMS2 <https://koina.wilhelmlab.org/docs#post-/Prosit_2024_intensity_XL_NMS2/infer>`_   | Developed for HCD non-cleavable cross-linked peptides (DSS, BS3). Supports oxidation and carbamidomethylation.                                                                               |
    +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | `AlphaPept_ms2_generic <https://koina.wilhelmlab.org/docs#post-/AlphaPept_ms2_generic/infer>`_                   | Developed for generic data support, including TMT, timsTOF and various instrument types.                                                                                                     |
    +------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | `ms2pip_2021_HCD <https://koina.wilhelmlab.org/docs#post-/ms2pip_2021_HCD/infer>`_                               | Developed for HCD tryptic and non-tryptic peptides.                                                                                                                                          |
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
