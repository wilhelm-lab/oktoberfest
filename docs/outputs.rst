Generated outputs
=================

Oktoberfest creates a variety of files and groups them within the user defined output directory. The following tables provide an overview of all subdirectories and files created depending on the job that is executed.

General directory structure
---------------------------

.. table::
    :class: fixed-table

    +-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Directory             | Description                                                                                                                                                                                                                                                                                                                                                                                                                   |
    +=======================+===============================================================================================================================================================================================================================================================================================================================================================================================================================+
    | data/                 | Contains hdf5 files that combine search results, annotated spectra, metadata in <spectra_file>.mzML.hdf5 and predictions in <spectra_file>.mzML.pred.hdf5 where <spectra_file> is replaced with the specific name of the RAW file for which information is stored. The files are updated and store the progress of the current job and enable skipping specific steps when rerunning a job. For GenerateTrainingData jobs,    |
    |                       | this folder also contains the training data for DLOmix in <spectra_file>.parquet (all annotated PSMs) and <spectra_file>.fdr1.parquet (target PSMs at 1% FDR), see `Training Data Generation <./jobs.html#training-data>`_.                                                                                                                                                                                                   |
    +-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | msms/                 | Contains combined msms.prosit and separated search results in <spectra_file>.prosit where <spectra_file> is replaced with the spectra file name for which search results are stored. The files are stored in the `internal format <./internal_format.html>`_  and are created as part of preprocessing search results from supported search engines. If a file is present, preprocessing is skipped when rerunning a job.     |
    +-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | spectra/              | Contains spectra files in mzML / hdf format (<spectra_file>.[mzML,hdf] where <spectra_file> is replaced with the specific name of the RAW file / .d folder that was converted). If spectra are provided in mzML / hdf format already, they will not show up in this folder and conversion is skipped when (re)running a job.                                                                                                  |
    +-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | proc/                 | Contains information on successful job steps in <job_step>.done where <job_step> is replaced by the specfic job step completed. If a file is present, the specific step will be skipped when rerunning a job.                                                                                                                                                                                                                 |
    +-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | results/              | Contains all results of a specific job. See below results table for an overview of all files stored in this directory.                                                                                                                                                                                                                                                                                                        |
    +-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | results/              | Subdirectory in results/ containing the collision energy calibration results and plots. See below results table for details on the files in this directory.                                                                                                                                                                                                                                                                   |
    | ce_calibration/       |                                                                                                                                                                                                                                                                                                                                                                                                                               |
    +-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | results/percolator/   | Subdirectory in results/ containing all inputs for running percolator and the corresponding results and summary plots. See below percolator table for details on the files in this directory.                                                                                                                                                                                                                                 |
    |                       |                                                                                                                                                                                                                                                                                                                                                                                                                               |
    +-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Files in the results directory
------------------------------

.. table::
    :class: fixed-table

    +---------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | File                                  | Description                                                                                                                                                                                  |
    +=======================================+==============================================================================================================================================================================================+
    | ce_calibration/<spectra_file>_ce.txt  | Contains the calibrated collision energy as JSON where <spectra_file> is replaced with the name of the spectra file for which the calibration was performed. Only HCD spectra are            |
    |                                       | calibrated: for DDA, one calibrated CE per original CE ({original CE: calibrated CE}); for DIA, one calibrated CE per precursor charge ({charge: calibrated CE}). If the file contains no    |
    |                                       | HCD spectra, the fixed CE of 35 is stored. Non-HCD spectra (e.g. CID) always use a CE of 35. The corresponding spectral angle plots (<spectra_file>_mean_spectral_angle_ce_<group>.svg,      |
    |                                       | <spectra_file>_violin_spectral_angle_ce_<group>.svg, or <spectra_file>_ce_ransac_model_<charge>.svg when using the RANSAC model) are stored in the same directory.                           |
    +---------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | predicted_library.<csv|msp|dlib>      | Contains the spectral library in either spectronaut, msp, or dlib format depending on what format was chosen in the configuration file                                                       |
    +---------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | peptides_internal.csv                 | Contains the peptide list including metadata needed for spectral library generation when providing a fasta file for in-silico digestion or a peptides list                                   |
    +---------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | spec_lib.tsv                          | Only for Rescoring jobs with acquisitionMethod "DIA": spectral library in OpenSWATH / DIA-NN compatible tsv format, built from the predicted intensities of all peptides (per precursor      |
    |                                       | charge) passing 1% FDR after rescoring with percolator. Fragments with a predicted relative intensity <= 0.05 are removed.                                                                   |
    +---------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


Files in the percolator directory
---------------------------------

This section details the files located within the percolator directory. In below file names, <method> is percolator, the tool used for target / decoy separation and FDR estimation.

.. table::
    :class: fixed-table

    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | File                                  | Description                                                                                                                                                                                 |
    +=======================================+=============================================================================================================================================================================================+
    | rescore.tab                           | Input file for percolator containing the features from peptide property prediction (see `features for target / decoy separation <./svm_features.html>`_)                                    |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | <spectra_file>.rescore.tab            | Input file for percolator containing the features from peptide property prediction for a spectra file name that replaces <spectra_file> (files are concatenated into rescore.tab)           |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | rescore.<method>.peptides.txt         | Percolator output for target peptides (features from peptide property prediction)                                                                                                           |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | rescore.<method>.psms.txt             | Percolator output for target PSMs (features from peptide property prediction)                                                                                                               |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | rescore.<method>.decoy.peptides.txt   | Percolator output for decoy peptides (features from peptide property prediction)                                                                                                            |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | rescore.<method>.decoy.psms.txt       | Percolator output for decoy PSMs (features from peptide property prediction)                                                                                                                |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | rescore.percolator.weights.csv        | Percolator output about individual feature weights (features from peptide property prediction)                                                                                              |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | rescore.log                           | Log output from percolator (features from peptide property prediction)                                                                                                                      |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | original.tab                          | Input file for percolator containing only the original features (see `features for target / decoy separation <./svm_features.html>`_)                                                       |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | <spectra_file>.original.tab           | Input file for percolator containing only the original features for a spectra file name that replaces <spectra_file> (files are concatenated into original.tab)                             |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | original.<method>.peptides.txt        | Percolator output for target peptides (original features)                                                                                                                                   |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | original.<method>.psms.txt            | Percolator output for target PSMs (original features)                                                                                                                                       |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | original.<method>.decoy.peptides.txt  | Percolator output for decoy peptides (original features)                                                                                                                                    |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | original.<method>.decoy.psms.txt      | Percolator output for decoy PSMs (original features)                                                                                                                                        |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | original.percolator.weights.csv       | Percolator output about individual feature weights (original features)                                                                                                                      |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | original.log                          | Log output from percolator (original features)                                                                                                                                              |
    +---------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

The following output files are specific to XL‑MS rescoring.

.. table::
    :class: fixed-table

    +---------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | File                                  | Description                                                                                                                                                                                                                                              |
    +=======================================+==========================================================================================================================================================================================================================================================+
    | percolator_xifdr_input.csv            | Contains Percolator scores for each CSM (cross‑linked spectrum match). Oktoberfest does **not** perform FDR estimation for XL‑MS data — this file can be uploaded directly to `xiFDR <https://github.com/Rappsilber-Laboratory/xiFDR>`_ to estimate FDR. |
    +---------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | xisearch_xifdr_input.csv              | Contains XiSearch scores for each CSM. You can feed this file to xiFDR as well, allowing a direct comparison between original XiSearch scores and rescored (Percolator) results.                                                                         |
    +---------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

PTM FDR control and localization
--------------------------------

With ``ptm_fdr_control`` or ``ptm_localization``, the percolator directory contains subdirectories for the modified
PSMs. Each subdirectory has its own percolator input (``rescore.tab``) and outputs (``rescore.percolator.*``, as in
the table above).

- **Separate FDR control** (``ptm_fdr_control``): the unmodified PSMs are rescored in place (``rescore.tab`` and
  ``rescore.percolator.*`` in the percolator directory), the modified PSMs in ``mod_fdr/``.
- **Localization** (``ptm_localization``): ``original.tab`` only contains the localizations of the search engine, so
  ``original.percolator.*`` is the search engine baseline. ``rescore.tab`` contains all site candidates, which are
  scored in ``localize_mod/`` (percolator training on the localization features) and ``localize_scores/`` (all
  candidates scored with the trained weights); with ``ptmLocalizationOptions.sdi_features``, a first pass without
  them is stored in ``localize_mod_first_pass/`` and ``localize_scores_first_pass/``. The PSMs with their most likely
  localization are rescored, with the localization probability and false localization rate in
  ``rescore.percolator_with_prob.psms.txt``:

  - phosphorylation: ``localize_phospho/`` (``localize_rt_all_<residue>/`` per residue if P or H are possible sites)
  - other PTMs and combinations of PTMs: ``mod_fdr/``, together with the PSMs without the localized PTMs

- **Citrullination** (UNIMOD:7 as the only localized PTM, without ``decoy_residues``): the PSMs are split like for the
  separate FDR control. Spectra without a citrullinated candidate are rescored in place (``rescore.tab`` and
  ``rescore.percolator.*`` in the percolator directory). Spectra with a citrullinated candidate are rescored in
  ``localize_cit_spectra/`` (best candidate per spectrum), and their citrullinated PSMs in ``localize_mod/``, the
  citrullination results. With ``decoy_residues``, citrullination is localized like the other PTMs (``mod_fdr/``).
