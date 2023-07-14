Usage Principles
================

A. CE Calibration (CollisionEnergyCalibration)
----------------------------------------------

This task estimates the optimal collision energy (CE) based on a given search result.
Prosit will:

1. Select a random subset of high-scoring PSMs
2. Predict those in for each CE from 18 to 49.
3. Calculate which CE achieves highest correlations with the experimental spectra
   Please note: Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported. Each C is treated as Cysteine with carbamidomethylation (fixed modification in MaxQuant).

Example config file:

.. code-block:: python

    task_config_ce_calibration = {
        "type": "CollisionEnergyCalibration",
        "tag": "",
        "output": "path_to_output_folder",
        "inputs": {
            "search_results": "path_to_msms",
            "search_type": "Maxquant",
            "spectra": "path_to_spectra_files",
            "spectra_type": "raw"
        },
        "models": {
            "intensity": "Prosit_2020_intensity_HCD",
            "irt": "Prosit_2019_irt"
        },
        "outputFormat": "",
        "prediction_server": "koina.proteomicsdb.org:443",
        "regressionMethod": "lowess",
        "ssl": True,
        "thermoExe": "ThermoRawFileParser.exe"
    }

B. Spectral Library (SpectralLibraryGeneration)
-----------------------------------------------

This task generates a spectral library either by digesting a given FASTA file, or by predicting a list of peptides given in a CSV file. You need to provide a collision energy (CE) for prediction. To estimate an optimal CE for prediction, please use "CE Calibration".
When a FASTA file is provided, Oktoberfest will:

1. Digest the FASTA, for the given parameters (i.e. protease).
2. Predict all spectra at the given collision energy.
   When a CSV with peptides is provided, Prosit will directly predict all spectra.

Example config file:

.. code-block:: python

    task_config_spectral_lib = {
        "type": "SpectralLibraryGeneration",
        "tag": "",
        "output": "path_to_output_folder",
        "inputs": {
            "search_results": "path_to_msms",
            "search_type": "Maxquant",
            "library_input": "path_to_peptides_csv,
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


C. Rescoring (Rescoring)
------------------------

This task rescores an existing search result using features generated from peptide property prediction.
Oktoberfest will:

1. Calibrate CE against the provided RAW files.
2. Predict all sequences in the search results file, e.g. msms.txt from MaxQuant
3. Use predicted spectra to generate features for percolator.
4. Run percolator to rescore the search.
   Please note: You need to provide search results that were not filtered for a given FDR (i.e. 100% FDR), otherwise valid targets may be filtered out prior to rescoring. Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported. Each C is treated as Cysteine with carbamidomethylation (fixed modification in MaxQuant).

Example config file:

.. code-block:: python

    task_config_rescoring = {
        "type": "Rescoring",
        "tag": "",
        "output": "path_to_output_folder",
        "inputs": {
            "search_results": "path_to_msms",
            "search_type": "Maxquant",
            "spectra": "path_to_spectra_files",
            "spectra_type": "raw"
        },
        "models": {
            "intensity": "Prosit_2020_intensity_HCD",
            "irt": "Prosit_2019_irt"
        },
        "outputFormat": "",
        "prediction_server": "koina.proteomicsdb.org:443",
        "numThreads": 4,
        "fdr_estimation_method": "mokapot",
        "allFeatures": False,
        "regressionMethod": "lowess",
        "ssl": True,
        "thermoExe": "ThermoRawFileParser.exe"
    }

