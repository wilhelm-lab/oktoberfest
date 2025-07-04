Running Oktoberfest
===================

This documentation provides an overview of the high-level job API and how to execute jobs for collision energy calibration, spectral library generation and rescoring.


Executing a job
---------------

Oktoberfest can be run in three different ways. The only input required is a configuration for setting up the job (we provide examples further down and a full documentation on the next page).

The command for executing a job from terminal:

.. code-block:: bash

   python -m oktoberfest --config_path <path/to/config.json>

The command for executing a job within python:

.. code-block:: python

   from oktoberfest.runner import run_job
   run_job("<path/to/config.json>")

If you instead want to run oktoberfest using the docker image, run:

.. code-block:: bash

   DATA=path/to/data/dir make run_oktoberfest

.. note::
    When using with docker, `DATA` must contain the spectra, the search results that fit the specified `search_results_type` in the config, and a `config.json` file with the configuration. The results will be written to `<DATA>/<output>/results/percolator`.


A. Collision Energy Calibration
-------------------------------

This task estimates the optimal normalised collision energy (NCE) based on a given search result.
Oktoberfest will:

1. Select the 1000 highest scoring target PSMs
2. Perform peptide property prediction for NCE 18 to 49 in steps of one.
3. Calculate the spectral angle between predicted and experimentally observed fragment intensities for each NCE and report the best NCE, i.e the one that reaches the highest spectral angle.

.. note::
    Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported.

    Each C is treated as Cysteine with carbamidomethylation (fixed modification).

Example config file:

.. code-block:: json

    {
        "type": "CollisionEnergyCalibration",
        "tag": "",
        "output": "./out",
        "inputs": {
            "search_results": "./msms.txt",
            "search_results_type": "Maxquant",
            "spectra": "./",
            "spectra_type": "raw",
            "instrument_type": "QE"
        },
        "models": {
            "intensity": "Prosit_2020_intensity_HCD",
            "irt": "Prosit_2019_irt"
        },
        "prediction_server": "koina.wilhelmlab.org:443",
        "numThreads": 1,
        "regressionMethod": "spline",
        "ssl": true,
        "thermoExe": "ThermoRawFileParser.exe",
        "massTolerance": 20,
        "unitMassTolerance": "ppm",
        "ce_alignment_options": {
            "ce_range": [19,50],
            "use_ransac_model": false
        }
    }

The example config can be loaded and viewed using

.. code-block:: python

    import oktoberfest as ok
    import json
    config = ok.utils.example_configs.CECALIB
    json.dumps(config, indent=4)


B. Spectral Library Generation
------------------------------

This task generates a spectral library either by digesting a given FASTA file, or by predicting a list of peptides given in a CSV file. You need to provide a collision energy (CE) for prediction (see above).
Oktoberfest will:
1. Digest the FASTA using a given protease and other parameters and create a peptides.csv file from that.
2. Predict all spectra at the given collision energy.

In case a CSV with peptides is provided, Oktoberfest will directly predict all spectra and skip the digestion step.

.. note::
    Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported.

    Each C is treated as Cysteine with carbamidomethylation (fixed modification).

Example config file:

.. code-block:: json

    {
        "type": "SpectralLibraryGeneration",
        "tag": "",
        "output": "./out",
        "inputs": {
            "library_input": "uniprot.fasta",
            "library_input_type": "fasta",
            "instrument_type": "QE"
        },
        "models": {
            "intensity": "Prosit_2020_intensity_HCD",
            "irt": "Prosit_2019_irt"
        },
        "spectralLibraryOptions": {
            "fragmentation": "HCD",
            "collisionEnergy": 30,
            "precursorCharge": [2,3],
            "minIntensity": 5e-4,
            "nrOx": 1,
            "batchsize": 10000,
            "format": "msp"
        },
        "fastaDigestOptions": {
            "digestion": "full",
            "missedCleavages": 2,
            "minLength": 7,
            "maxLength": 60,
            "enzyme": "trypsin",
            "specialAas": "KR",
            "db": "concat"
        },
        "prediction_server": "koina.wilhelmlab.org:443",
        "numThreads": 1,
        "ssl": true
    }

The example config can be loaded and viewed using

.. code-block:: python

    import oktoberfest as ok
    import json
    config = ok.utils.example_configs.LIBGEN
    json.dumps(config, indent=4)


C. Rescoring
------------

a) without refinement
~~~~~~~~~~~~~~~~~~~~~

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

.. code-block:: json

    {
        "type": "Rescoring",
        "tag": "",
        "output": "./out",
        "inputs": {
            "search_results": "./msms.txt",
            "search_results_type": "Maxquant",
            "spectra": "./",
            "spectra_type": "raw",
            "instrument_type": "QE"
        },
        "models": {
            "intensity": "Prosit_2020_intensity_HCD",
            "irt": "Prosit_2019_irt"
        },
        "prediction_server": "koina.wilhelmlab.org:443",
        "numThreads": 1,
        "fdr_estimation_method": "mokapot",
        "add_feature_cols": "none",
        "regressionMethod": "spline",
        "ssl": true,
        "thermoExe": "ThermoRawFileParser.exe",
        "massTolerance": 20,
        "unitMassTolerance": "ppm",
        "ce_alignment_options": {
            "ce_range": [19,50],
            "use_ransac_model": false
        }
    }

The example config can be loaded and viewed using

.. code-block:: python

    import oktoberfest as ok
    import json
    config = ok.utils.example_configs.RESCORING
    json.dumps(config, indent=4)


For rescoring tasks including quantification via picked-group-FDR, create a config file like this (so far only MaxQuant is supported):

.. code-block:: json

    {
        "type": "Rescoring",
        "quantification": true,
        "tag": "",
        "inputs": {
            "search_results": "mq_results/txt",
            "search_results_type": "Maxquant",
            "spectra": "./",
            "spectra_type": "raw",
            "library_input": "uniprot.fasta"
        },
        "output": "./out",
        "models": {
            "intensity": "Prosit_2020_intensity_HCD",
            "irt": "Prosit_2019_irt"
        },
        "prediction_server": "koina.proteomicsdb.org:443",
        "ssl": true,
        "thermoExe": "/opt/compomics/ThermoRawFileParser1.4.3/ThermoRawFileParser.exe",
        "numThreads": 1,
        "fdr_estimation_method": "percolator",
        "regressionMethod": "spline",
        "massTolerance": 20,
        "unitMassTolerance": "ppm",
        "fastaDigestOptions": {
            "digestion": "full",
            "missedCleavages": 2,
            "minLength": 7,
            "maxLength": 60,
            "enzyme": "trypsin",
            "specialAas": "KR",
            "db": "concat"
        }
    }


The example config can be loaded and viewed using

.. code-block:: python

    import oktoberfest as ok
    import json
    config = ok.utils.example_configs.RESCORING_WITH_QUANT
    json.dumps(config, indent=4)


In addition, XL-MS example configurations for non-cleavable and cleavable crosslinkers can be found using instead one of

.. code-block:: python

    config = ok.utils.example_configs.RESCORING_XL_NON_CLEAVABLE
    config = ok.utils.example_configs.RESCORING_XL_CLEAVABLE


b) with refinement
~~~~~~~~~~~~~~~~~~

Same as rescoring without refinement, but in addition a new intensity predictor will be trained from a baseline model using off-line reinforcement learning on the provided spectra.
The refined intensity predictor will be used along an on-line retention time predictor to generate inputs for rescoring.

.. note::
    You can either provide a baseline predictor yourself, or download a pre-trained one from GitHub automatically.

Example config file:

.. code-block:: json

    {
        "type": "Rescoring",
        "tag": "",
        "output": "./out",
        "inputs": {
            "search_results": "./msms.txt",
            "search_results_type": "Maxquant",
            "spectra": "./",
            "spectra_type": "raw",
            "instrument_type": "QE"
        },
        "models": {
            "intensity": "baseline",
            "irt": "Prosit_2019_irt"
        },
        "prediction_server": "koina.wilhelmlab.org:443",
        "numThreads": 1,
        "dlomixInferenceBatchSize": 1024,
        "refinementLearningOptions": {
            "batchSize": 1024,
            "includeOriginalSequences": false,
            "improveFurther": false,
            "datasetFilteringOptions": {
                "searchEngineScoreThreshold": 0,
                "numDuplicates": 100
            }
        },
        "fdr_estimation_method": "mokapot",
        "allFeatures": false,
        "regressionMethod": "spline",
        "ssl": true,
        "thermoExe": "ThermoRawFileParser.exe",
        "massTolerance": 20,
        "unitMassTolerance": "ppm",
        "ce_alignment_options": {
            "ce_range": [19,50],
            "use_ransac_model": false
        }
    }

The example config can be loaded and viewed using

.. code-block:: python

    import oktoberfest as ok
    import json
    config = ok.utils.example_configs.RESCORING_WITH_REFINEMENT
    json.dumps(config, indent=4)
