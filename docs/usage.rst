Usage Principles
================

Run a job
---------

The general command for executing any job is:

.. code-block:: bash

   python oktoberfest/run_oktoberfest.py --config_path path_to_config_file

If you instead want to run oktoberfest using the docker image, run:

.. code-block:: bash

   DATA=path/to/data/dir make run_oktoberfest

Note: When using with docker, `DATA` must contain the spectra, the search results that fit the specified `search_results_type` in the config, and a `config.json` file with the configuration. The results will be written to `<DATA>/<output>/results/percolator`.

Config flags explained
----------------------

- `type` = "CollisionEnergyAlignment", "SpectralLibraryGeneration" or "Rescoring"
- `tag` = "tmt", "tmtpro", "itraq4" or "itraq8"; default is ""
- `fdr_estimation_method` = method used for FDR estimation on PSM and peptide level: "percolator" or "mokapot"; default = "mokapot"
- `allFeatures`` = True if all features should be used for FDR estimation; default = False
- `regressionMethod` = regression method for curve fitting (mapping from predicted iRT values to experimental retention times): "lowess", "spline" or "logistic"; default = "lowess"
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

Jobs
----

A. Collision Energy Calibration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This task estimates the optimal normalised collision energy (NCE) based on a given search result.
Oktoberfest will:

1. Select the 1000 highest scoring target PSMs
2. Perform peptide property prediction for NCE 18 to 49 in steps of one.
3. Calculate the spectral angle between predicted and experimentally observed fragment intensities for each NCE and report the best NCE, i.e the one that reaches the highest spectral angle.

.. note::
    Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported.

    Each C is treated as Cysteine with carbamidomethylation (fixed modification).

Example config file:

.. code-block:: python

    task_config_ce_calibration = {
        "type": "CollisionEnergyCalibration",
        "tag": "",
        "output": "./out",
        "inputs": {
            "search_results": "./msms.txt",
            "search_results_type": "Maxquant",
            "spectra": "./",
            "spectra_type": "raw"
        },
        "models": {
            "intensity": "Prosit_2020_intensity_HCD",
            "irt": "Prosit_2019_irt"
        },
        "prediction_server": "koina.proteomicsdb.org:443",
        "regressionMethod": "lowess",
        "ssl": True,
        "thermoExe": "ThermoRawFileParser.exe",
        "massTolerance": 20,
        "unitMassTolerance": "ppm"
    }

B. Spectral Library Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This task generates a spectral library either by digesting a given FASTA file, or by predicting a list of peptides given in a CSV file. You need to provide a collision energy (CE) for prediction (see above).
Oktoberfest will:
1. Digest the FASTA using a given protease and other parameters and create a peptides.csv file from that.
2. Predict all spectra at the given collision energy.

In case a CSV with peptides is provided, Oktoberfest will directly predict all spectra and skip the digestion step.

.. note::
    Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported.

    Each C is treated as Cysteine with carbamidomethylation (fixed modification).

Example config file:

.. code-block:: python

    task_config_spectral_lib = {
        "type": "SpectralLibraryGeneration",
        "tag": "",
        "output": "./out",
        "inputs": {
            "search_results": "./msms.txt",
            "search_results_type": "Maxquant",
            "library_input": "./peptides.csv",
            "library_input_type": "peptides"
        },
        "models": {
            "intensity": "Prosit_2020_intensity_HCD",
            "irt": "Prosit_2019_irt"
        },
        "outputFormat": "spectronaut",
        "prediction_server": "koina.proteomicsdb.org:443",
        "numThreads": 1,
        "ssl": True,
        "thermoExe": "ThermoRawFileParser.exe"
        "fastaDigestOptions": {
            "fragmentation": "",
            "digestion": "full",
            "missedCleavages": 2,
            "minLength": 7,
            "maxLength": 60,
            "enzyme", "trypsin",
            "specialAas": "KR",
            "db": "concat"
    }


C. Rescoring
~~~~~~~~~~~~

This task rescores an existing search result using features generated from peptide property prediction.
Oktoberfest will:

1. Calibrate CE against the provided RAW files.
2. Perform peptide property prediction for all spectra that have a match in the search results file.
3. Use predicted spectra and retention time to generate features for rescoring.
4. Run percolator or mokapot to rescore the search and perform FDR estimation.
5. Generate summary plots.

.. note::
    You need to provide search results that were not filtered for a given FDR (i.e. 100% FDR), otherwise valid targets may be filtered out prior to rescoring.

    Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported.

    Each C is treated as Cysteine with carbamidomethylation (fixed modification).

Example config file:

.. code-block:: python

    task_config_rescoring = {
        "type": "Rescoring",
        "tag": "",
        "output": "./out",
        "inputs": {
            "search_results": "./msms.txt",
            "search_results_type": "Maxquant",
            "spectra": "./",
            "spectra_type": "raw"
        },
        "models": {
            "intensity": "Prosit_2020_intensity_HCD",
            "irt": "Prosit_2019_irt"
        },
        "prediction_server": "koina.proteomicsdb.org:443",
        "numThreads": 1,
        "fdr_estimation_method": "mokapot",
        "allFeatures": False,
        "regressionMethod": "lowess",
        "ssl": True,
        "thermoExe": "ThermoRawFileParser.exe",
        "massTolerance": 20,
        "unitMassTolerance": "ppm"
    }

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

+----------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|          Model             |                             Description                                                                                                                |
+============================+========================================================================================================================================================+
| Prosit_2019_intensity      | deprecated, please use the 2020 model                                                                                                                  |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| Prosit_2020_intensity_HCD  | your go to model for fragment intensity prediction for HCD fragmentation, find out more about this model `here <https://github.com/kusterlab/prosit>`_ |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| Prosit_2020_intensity_CID  | your go to model for fragment intensity prediction for CID fragmentation, find out more about this model `here <https://github.com/kusterlab/prosit>`_ |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| Prosit_2020_intensity_TMT  | your go to model for fragment intensity prediction for TMT, find out more about this model `here <https://github.com/kusterlab/prosit>`_               |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| Prosit_2019_irt            | all purpose model for retention time prediction, find out more about this model `here <https://github.com/kusterlab/prosit>`_                          |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| Prosit_2020_irt_TMT        | your go to model for retention time prediction for TMT, find out more about this model `here <https://github.com/kusterlab/prosit>`_                   |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+

Once support for additional models is added, they will be added here.

For the `prediction_server` flag, you should use the `koina <https://koina.proteomicsdb.org/>`_ instance we provide at `koina.proteomicsdb.org:443`.
For models, you should choose the models that fit your use case. You can see available models for the prediction server we offer at `<https://koina.proteomicsdb.org/docs>`.
For a list of currently tested models, check the "Supported Models" section below.
