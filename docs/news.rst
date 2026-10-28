News
====

.. role:: date
    :class: date

Oktoberfest supports DIA data and local DLOmix models :date:`2026-10-12`
------------------------------------------------------------------------

| Oktoberfest 1.0.0 is published, adding support for data-independent acquisition (DIA). Set ``acquisitionMethod`` to ``"DIA"`` in the `configuration <config.html>`_ to calibrate the collision energy per precursor charge.
| Collision energy calibration for DDA now calibrates each original collision energy separately. Only HCD spectra are calibrated, and all results and plots are stored in ``results/ce_calibration/``.
| Fragment intensities can now be predicted locally with a trained `DLOmix <https://github.com/wilhelm-lab/dlomix>`_ model instead of Koina. Install the optional extra with ``pip install "oktoberfest[dlomix]"`` and see the DLOmix integration tutorial in the `tutorials folder on GitHub <https://github.com/wilhelm-lab/oktoberfest/tree/development/tutorials>`_.
| The new ``GenerateTrainingData`` job exports annotated spectra from your own search results as training data for DLOmix, which can be used to fine-tune any Prosit intensity model on your data (see `Training Data Generation <jobs.html#training-data>`_).
| Oktoberfest now comes with a `web UI <webui.html>`_ that can be run locally via Docker to submit and monitor jobs without writing config files.
| **Breaking changes:** Oktoberfest now requires spectrum_fundamentals >=1.0.0. The ``ExpMass`` column of the percolator / mokapot input is now the experimental precursor m/z (see `SVM features <svm_features.html>`_). The collision energy calibration results (``<spectra_file>_ce.txt``, now stored as JSON) and plots moved to ``results/ce_calibration/`` (see `outputs <outputs.html>`_).

Oktoberfest supports rescoring XL-MS datasets :date:`2025-06-14`
----------------------------------------------------------------

| Oktoberfest 0.10.0 is published, introducing support for crosslinking mass spectrometry (XL-MS) rescoring — including both cleavable (DSSO, DSBU) and non-cleavable (DSS, BS3) crosslinkers.
| The update adds compatibility with `xiSEARCH <https://www.rappsilberlab.org/software/xisearch/>`_ output and integrates the new Prosit-XL models for XL-MS data.
| Check out the publication in `Nature Communications <https://www.nature.com/articles/s41467-025-61203-4>`_ and explore the `full tutorial notebook on GitHub <https://github.com/wilhelm-lab/oktoberfest/tree/development/tutorials>`_.
| Related config files can be found here: `example_configs.py <https://github.com/wilhelm-lab/oktoberfest/blob/development/oktoberfest/utils/example_configs.py>`_.
| Oktoberfest is now also supporting python 3.12.

Oktoberfest supports OpenMS and custom modifications :date:`2025-05-21`
-----------------------------------------------------------------------

| Oktoberfest 0.9.0 is published, allowing to provide `custom modifications <custom_mods.html>`_ as long as they are supported by the respective Koina model. In addition, direct support of OpenMS searches was added.
| This release supports python 3.10, and 3.11. Support for python 3.9 was dropped!

Oktoberfest supports requantifying rescored peptides with Picked-Group-FDR and offline predictions with custom models :date:`2024-10-10`
----------------------------------------------------------------------------------------------------------------------------------------

| Oktoberfest 0.8.0 is published using Picked-Group-FDR, Oktoberfest can now requantify results from MaxQuant runs which enables downstream comparison of rescored with original results on protein level.
| Additionally, custom features can now be added and configured to be used by percolator/mokapot for rescoring.
| This release also adds preliminary support of local predictions via DLOmix. Stay tuned for a documentation update including jupyter notebooks.

Oktoberfest provides access to AlphaPept and MS2PIP predictions :date:`2024-05-29`
----------------------------------------------------------------------------------

| Oktoberfest 0.7.0 is published, providing access to predictions from AlphaPept and MS2PIP. This comes with a large API overhaul that allows easier access to underlying data for manual insights and analysis.
| Please note that you can now specify the instrument type in the configuration for AlphaPept predictions.
| Spectral libraries can now be generated using a simple peptide list using a variety of options, such as precursor charge, max. number of oxidations and collision energy, or directly in internal format. Check the new `peptides / internal format specification <./peptides_format.html>`_ and `configuration options <./config.html>`_.
| Stay tuned for a documentation update including jupyter notebooks on how to navigate and make use of your data for manual in-depth analysis.
| This release drops python 3.8 support!

Support for timsTOF added / spectral library generation overhaul :date:`2024-01-31`
-----------------------------------------------------------------------------------

| Oktoberfest 0.6.0 is out supporting the new Prosit model for timsTOF data. Check out the preprint at `BioRxiv <https://doi.org/10.1101/2023.07.17.549401>`__.
| Spectral library generation was completely overhauled, resulting in more stability, dramatic decreases in runtime and storage requirements and more customization in the config.
| Note that this release introduces changes in the configuration file that might break your previous setup! Please check the `new configuration format <./config.html>`_.

Oktoberfest now provides a low level API :date:`2023-10-03`
-----------------------------------------------------------

| With our new release 0.5.0, we refactored our code base to provide a low level API to the users. Check the new `API <./API.html>`_ section for more details.
| Please note: We dropped support for percolator version 3.05. Please update to version 3.06.1. See `installation <./installation.html>`_ for details.


Oktoberfest is published! :date:`2023-09-07`
--------------------------------------------

| Alongside the new 0.4.0 release, Oktoberfest was published in Proteomics.
| Check out the publication here: `https://doi.org/10.1002/pmic.202300112 <https://doi.org/10.1002/pmic.202300112>`_.
| Accordingly, please check `How to cite <./reference.html>`_.


The first Oktoberfest release is ready! :date:`2023-06-25`
----------------------------------------------------------

Version 0.3.0 is available on `PyPI <https://pypi.org/project/oktoberfest/>`_.

We are happy to announce that our first standalone version of Oktoberfest is finally ready.
Besides a lot of bugfixes, new features and performance enhancements over the existing online Prosit service at `ProteomicsDB <https://proteomicsdb.org/prosit>`_, the Oktoberfest codebase is now open source and available on `GitHub <https://github.com/wilhelm-lab/oktoberfest>`_.
In addition, we switched to our new online prediction service `koina <https://koina.proteomicsdb.org>`_ which allows connecting Oktoberfest to any public or self-hosted koina instance for easy prediction model integration and being independent of ProteomicsDB!
