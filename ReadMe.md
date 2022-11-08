# Oktoberfest

## Getting started

What you can do with oktoberfest:

-   CE Calibration (CollisionEnergyAlignment)

This task estimates the optimal collision energy (CE) based on a given search result. You need to upload a RAW file as well as the MaxQuant's msms.txt for calibration.
Prosit will:

1. Select a random subset of high-scoring PSMs
2. Predict those in for each CE from 18 to 39.
3. Calculate which CE achieves highest correlations with the experimental spectra
   Please note: Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported. Each C is treated as Cysteine with carbamidomethylation (fixed modification in MaxQuant).

-   Spectral Library (SpectralLibraryGeneration)

This task generates a spectral library either by digesting a given FASTA file, or by predicting a list of peptides given in a CSV file. You need to provide a collision energy (CE) for prediction. To estimate an optimal CE for prediction, please use "CE Calibration".
When a FASTA file is provided, Prosit will:

1. Digest the FASTA, for the given parameters (i.e. protease).
2. Predict all spectra at the given collision energy.
   When a CSV with peptides is provided, Prosit will directly predict all spectra.

-   Rescoring (MaxQuantRescoring)

This task rescores an existing MaxQuant search (FDR 100%) using features generated from fragmentation prediction. You need to upload a RAW file as well as the MaxQuant's msms.txt file from a search.
Prosit will:

1. Calibrate itself against the RAW.
2. Predict all sequences in the msms.txt.
3. Use the predicted spectra to generate features for percolator.
4. Run percolator to rescore the search.
   Please note: You need a MaxQuant search at 100% FDR, otherwise targets may be filtered by MaxQuant's FDR calculation before rescoring. Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported. Each C is treated as Cysteine with carbamidomethylation (fixed modification in MaxQuant).

## Installation

After cloning oktoberfest repository, create a new conda environment:

```bash
conda create -n oktoberfest python=3.8
```

Then, go to oktoberfest repository folder and run:

```bash
pip install .
```

Create a config.json file which should contain the following flags:

-   jobType = CollisionEnergyAlignment, MaxQuantRescoring or SpectralLibraryGeneration

-   tag = tmt, tmtpro, itraq4 or itraq8; default = tmt

-   allFeatures = true if all features should be used by the percolator; default = false

-   fileUploads

    -   search_type = maxquant or internal; default = maxquant

    -   raw_type = thermo or mzml; default = thermo

    -   fasta = path to the fasta file

    -   peptides.csv =

-   models

    -   intensity = intensity model

    -   irt = irt model

    -   proteotypicity = proteotypicity model

-   prosit_server =

-   numThreads = number of threads from the config file; default = 1

-   jobId =

An example of the config file can be found in `/oktoberfest/example_config.json`.

Finally, go to `/oktoberfest/oktoberfest/` and run

```bash
python oktoberfest.py —search_dir path_to_search_dir —config_path path_to_config_file
```

Note: The search_dir should contain both the raw files and the MaxQuant's `msms.txt` from a search.

## Models
