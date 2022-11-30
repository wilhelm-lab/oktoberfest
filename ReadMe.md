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

After cloning the repository of oktoberfest, create a new conda environment:

```bash
conda create -n oktoberfest
```

After activating the newly-created conda environment, go to the main folder with the repository of oktoberfest and run:

```bash
pip install .
```

Create a `config.json` file which should contain the following flags:

-   `jobType` = CollisionEnergyAlignment, MaxQuantRescoring or SpectralLibraryGeneration

-   `tag` = tmt, tmtpro, itraq4 or itraq8; default = tmt

-   `allFeatures` = true if all features should be used by the percolator; default = false

-   `fileUploads`

    -   `search_type` = maxquant, msfragger, mascot or internal; default = maxquant

    -   `raw_type` = thermo or mzml; default = thermo

    -   `fasta` = path to the FASTA file, if FASTA file is provided

    -   `peptides.csv` = true if you like to provide the list of peptides

-   `models`

    -   `intensity` = intensity model

    -   `irt` = irt model

    -   `proteotypicity` = proteotypicity model

-   `prosit_server` = server for the Prosit prediction

-   `numThreads` = number of threads from the config file; default = 1

-   `jobId` = job ID for the Prosit prediction

-   `searchPath` = path to the search file (if the search type is msfragger, then the path to the xlsx file should be provided); default = ""

The following flags are relevant only if a FASTA file is provided:

-   `fastaDigestOptions`

    -   `fragmentation` = fragmentation method (HCD or CID)

    -   `digestion` = digestion mode (full, semi or none); default = full

    -   `cleavages` = number of allowed missed cleavages used in the search engine; default = 2

    -   `minLength` = minimum peptide length allowed used in the search engine; default = 7

    -   `maxLength` = maximum peptide length allowed used in the search engine; default = 60

    -   `enzyme` = type of enzyme used in the search engine; default = trypsin

    -   `specialAas` = special amino acids used by MaxQuant for decoy generation; default = KR

    -   `db` = Target, decoy or concat; default = concat

An example of the config file can be found in `/oktoberfest/example_config.json`.

For `prosit_server` and `jobId`: ask Wassim Gabriel (wassim.gabriel@tum.de) or Ludwig Lautenbacher (Ludwig.Lautenbacher@tum.de).

Finally, run

```bash
python oktoberfest/run_oktoberfest.py —search_dir path_to_search_dir —config_path path_to_config_file
```

Note: The search_dir should contain both the raw files and the MaxQuant's `msms.txt` from a search.

## Models

-   `intensity models`

    -   `Prosit_2019_intensity`
    -   `Prosit_2020_intensity_hcd`
    -   `Prosit_2020_intensity_cid`
    -   `Prosit_2020_intensityTMT`
    -   `Prosit_2020_intensityTMT_Phospho`

-   `irt models`
-   `Prosit_2019_irt`
-   `Prosit_2020_irt_TMT`

-   `proteotypicity models`
-   `Prosit_2020_proteotypicity`
