.. include:: _key_contributors.rst

.. role:: small

.. role:: smaller

|PyPI| |Python Version| |License| |Read the Docs| |Build| |Tests| |Codecov| |pre-commit| |Black|

.. |PyPI| image:: https://img.shields.io/pypi/v/oktoberfest.svg
   :target: https://pypi.org/project/oktoberfest/
   :alt: PyPI
.. |Python Version| image:: https://img.shields.io/pypi/pyversions/oktoberfest
   :target: https://pypi.org/project/oktoberfest
   :alt: Python Version
.. |License| image:: https://img.shields.io/github/license/wilhelm-lab/oktoberfest
   :target: https://opensource.org/licenses/MIT
   :alt: License
.. |Read the Docs| image:: https://img.shields.io/readthedocs/oktoberfest/latest.svg?label=Read%20the%20Docs
   :target: https://oktoberfest.readthedocs.io/
   :alt: Read the documentation at https://oktoberfest.readthedocs.io/
.. |Build| image:: https://github.com/wilhelm-lab/oktoberfest/workflows/Build%20oktoberfest%20Package/badge.svg
   :target: https://github.com/wilhelm-lab/oktoberfest/actions?workflow=Package
   :alt: Build Package Status
.. |Tests| image:: https://github.com/wilhelm-lab/oktoberfest/workflows/Run%20oktoberfest%20Tests/badge.svg
   :target: https://github.com/wilhelm-lab/oktoberfest/actions?workflow=Tests
   :alt: Run Tests Status
.. |Codecov| image:: https://codecov.io/gh/wilhelm-lab/oktoberfest/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/wilhelm-lab/oktoberfest
   :alt: Codecov
.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit
.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Black

What you can do with oktoberfest:
=================================

CE Calibration (CollisionEnergyCalibration)
-------------------------------------------

This task estimates the optimal collision energy (CE) based on a given search result.
Prosit will:

1. Select a random subset of high-scoring PSMs.
2. Predict those for each CE from 18 to 49.
3. Calculate which CE achieves the highest correlations with the experimental spectra.

Please note: Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported. Each C is treated as Cysteine with carbamidomethylation (fixed modification in MaxQuant).

Spectral Library (SpectralLibraryGeneration)
--------------------------------------------

This task generates a spectral library either by digesting a given FASTA file or by predicting a list of peptides given in a CSV file. You need to provide a collision energy (CE) for prediction. To estimate an optimal CE for prediction, please use "CE Calibration".

When a FASTA file is provided, Oktoberfest will:

1. Digest the FASTA, considering the given parameters (i.e., protease).
2. Predict all spectra at the given collision energy.

When a CSV with peptides is provided, Prosit will directly predict all spectra.

Rescoring (Rescoring)
---------------------

This task rescores an existing search result using features generated from peptide property prediction.
Oktoberfest will:

1. Calibrate CE against the provided RAW files.
2. Predict all sequences in the search results file, e.g., msms.txt from MaxQuant.
3. Use predicted spectra to generate features for percolator.
4. Run percolator to rescore the search.

Please note: You need to provide search results that were not filtered for a given FDR (i.e., 100% FDR), otherwise valid targets may be filtered out prior to rescoring. Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported. Each C is treated as Cysteine with carbamidomethylation (fixed modification in MaxQuant).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   API
   contributing
   contributors
   reference

.. _github: https://github.com/wilhelm-lab/oktoberfest

GitHub Repository: :github:`wilhelm-lab/oktoberfest`

