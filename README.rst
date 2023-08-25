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

Oktoberfest: Rescoring and Spectral Library Generation for Proteomics
=====================================================================

Oktoberfest is a python tool for rescoring search results and generating spectral libraries for proteomics research within the Prosit ecosystem. It offers an end to end pipeline that takes search results, predicts peptide properties using koina, plots summaries and quality control figures and performs FDR estimation with either mokapot or percolator.

Installation
------------

Prerequisites
~~~~~~~~~~~~~

Oktoberfest requires python >=3.8,<=3.11. Best practise is to use a clean conda environment (`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_).

If you provide thermo raw files, make sure `ThermoRawFileParser <https://github.com/compomics/ThermoRawFileParser>`_ is installed.

If you are on linux or MacOS, make sure `mono <https://www.mono-project.com/>`_ is installed.

If you want to use percolator, make sure you install version 3.05 (`percolator <https://github.com/percolator/percolator/releases/tag/rel-3-05>`_).

Using pip (recommended)
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install oktoberfest

Docker image
~~~~~~~~~~~~

Prerequisites:
- `make <https://www.gnu.org/software/make/>`_
- `docker <https://www.docker.com/>`_

After cloning the repository of oktoberfest, checkout the branch you want to build the container from.
The latest stable version is always on the main branch.

.. code-block:: bash

   make build

Getting started
---------------

What you can do with oktoberfest:
- CE Calibration (CollisionEnergyCalibration)

  This task estimates the optimal collision energy (CE) based on a given search result.
  Prosit will:
  1. Select a random subset of high-scoring PSMs
  2. Predict those in for each CE from 18 to 49.
  3. Calculate which CE achieves highest correlations with the experimental spectra

  Please note: Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported. Each C is treated as Cysteine with carbamidomethylation (fixed modification).

- Spectral Library (SpectralLibraryGeneration)

  This task generates a spectral library either by digesting a given FASTA file, or by predicting a list of peptides given in a CSV file. You need to provide a collision energy (CE) for prediction. To estimate an optimal CE for prediction, please use "CE Calibration".
  When a FASTA file is provided, Oktoberfest will:
  1. Digest the FASTA, for the given parameters (i.e. protease).
  2. Predict all spectra at the given collision energy.

  When a CSV with peptides is provided, Prosit will directly predict all spectra.

- Rescoring (Rescoring)

  This task rescores an existing search result using features generated from peptide property prediction.
  Oktoberfest will:
  1. Calibrate CE against the provided RAW files.
  2. Predict all sequences in the search results file.
  3. Use predicted spectra and retention time to generate features for rescoring.
  4. Run percolator or mokapot to rescore the search and perform FDR estimation.

  Please note: You need to provide search results that were not filtered for a given FDR (i.e. 100% FDR), otherwise valid targets may be filtered out prior to rescoring. Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported. Each C is treated as Cysteine with carbamidomethylation (fixed modification).

Run oktoberfest
---------------

Configuration
~~~~~~~~~~~~~

Create a `config.json` file which should contain the following flags:

- `type` = "CollisionEnergyAlignment", "SpectralLibraryGeneration" or "Rescoring"
- `tag` = "tmt", "tmtpro", "itraq4" or "itraq8"; default is ""
- `fdr_estimation_method` = method used for FDR estimation on PSM and peptide level: "percolator" or "mokapot"; default = "mokapot"
- `allFeatures`` = True if all features should be used for FDR estimation; default = False
- `regressionMethod` = regression method for curve fitting (mapping from predicted iRT values to experimental retention times): "lowess", "spline" or "logistic"; default = "spline"
- `inputs`
   - `search_results` = path to the file containing the search results
   - `search_results_type` = the tool used to produce the search results, can be "Maxquant", "Msfragger", "Mascot" or "Internal"; default = "Maxquant"
   - `spectra` = path to a folder or a single file containing mass spectrometry results (raw or mzml files)
   - `spectra_type` = "raw" or "mzml"; default = "raw"
- `models`
   - `intensity` = intensity model
   - `irt` = irt model
- `prediction_server` = server for obtaining peptide property predictions
- `ssl` = Use ssl when making requests to the prediction server, can be true or false; default = true
- `numThreads` = number of raw files processed in parallel processes; default = 1
- `thermoExe` = path to ThermoRawFileParser executable; default "ThermoRawFileParser.exe"
- `massTolerance` = mass tolerance value defining the allowed tolerance between theoretical and experimentally observered fragment mass during peak filtering and annotation. Default depends on the mass analyzer: 20 (FTMS), 40 (TOF), 0.35 (ITMS)
- `unitMassTolerance` = unit for the mass tolerance, either "da" or "ppm". Default is da (mass analyzer is ITMS) and ppm (mass analyzer is FTMS or TOF)
- `output` = path to the output folder; if not provided the current working directory will be used.

For `prediction_server`, you should use the `koina <https://koina.proteomicsdb.org/>`_ instance we provide at `koina.proteomicsdb.org:443`.
For models, you should choose the models that fit your use case. You can see available models for the prediction server we offer at `https://koina.proteomicsdb.org/docs`.
For a list of currently tested models, check the "Supported Models" section below.

The following flags are relevant only for SpectralLibraryGeneration:

- `inputs`
   - `library_input` = path to the FASTA or peptides file
   - `library_input_type` = library input type: "fasta" or "peptides"
- `outputFormat` = "spectronaut" or "msp"

The following flags are relevant only if a FASTA file is provided:

- `fastaDigestOptions`
   - `fragmentation` = fragmentation method: "HCD" or "CID"
   - `digestion` = digestion mode: "full", "semi" or None; default = "full"
   - `cleavages` = number of allowed missed cleavages used in the search engine; default = 2
   - `minLength` = minimum peptide length allowed used in the search engine; default = 7
   - `maxLength` = maximum peptide length allowed used in the search engine; default = 60
   - `enzyme` = type of enzyme used in the search engine; default = "trypsin"
   - `specialAas` = special amino acids for decoy generation; default = "KR"
   - `db` = "target", "decoy" or "concat"; default = "concat"

An example of the config file can be found in `/oktoberfest/example_config.json`.

Run a job
---------

The general command for executing any job is:

.. code-block:: bash

   python oktoberfest/run_oktoberfest.py --config_path path_to_config_file

If you instead want to run oktoberfest using the docker image, run:

.. code-block:: bash

   DATA=path/to/data/dir make run_oktoberfest

Note: When using with docker, `DATA` must contain the spectra, the search results that fit the specified `search_results_type` in the config, and a `config.json` file with the configuration. The results will be written to `<DATA>/<output>/results/percolator` or `<DATA>/<output>/results/mokapot` depending on the chosen fdr estimation method.

Supported Models
----------------

This is the list of currently supported and tested models for peptide property prediction provided by `koina.proteomicsdb.org`:

- Intensity models:
   - Prosit_2019_intensity
   - Prosit_2020_intensity_HCD
   - Prosit_2020_intensity_CID
   - Prosit_2020_intensity_TMT

- iRT models:
   - Prosit_2019_irt
   - Prosit_2020_irt_TMT

Once support for additional models is added, they will be added here.

Tutorials and Documentation
---------------------------

We provide a Jupyter notebook that you can find at "tutorials/Oktoberfest Tutorial.ipynb", guiding you through the three different use cases using a public dataset.

If you want to test it inside your docker container, please refer to the README in the data/plasma subfolder.
The official Oktoberfest documentation can be found at `https://oktoberfest.readthedocs.io`.
Information about how to use koina and which models are supported by our public koina instance can be found at `https://koina.proteomicsdb.org/docs`.

License
-------

The project is licensed under the `MIT license <https://github.com/wilhelm-lab/oktoberfest/blob/main/LICENSE>`.

References
----------

[1] Gessulat S, Schmidt T, Zolg DP, Samaras P, Schnatbaum K, Zerweck J, Knaute T, Rechenberger J, Delanghe B, Huhmer A, Reimer U, Ehrlich HC, Aiche S, Kuster B, Wilhelm M: "PROSIT: Proteome-wide prediction of peptide tandem mass spectra by deep learning". Nature Methods. 2019; 16(6):509-518. doi: 10.1038/s41592-019-0426-7.

[2] Gabriel W, The M, Zolg D, Bayer FP, Shouman O, Lautenbacher L, Schnatbaum K, Zerweck J, Knaute T, Delanghe B, Huhmer A, Wenschuh H, Reimer U, Médard G, Kuster B, Wilhelm M: "Prosit-TMT: Deep Learning Boosts Identification of TMT-Labeled Peptides". Analytical Chemistry. 2022; 94(20):7181-7190. doi: 10.1021/acs.analchem.1c05435.
