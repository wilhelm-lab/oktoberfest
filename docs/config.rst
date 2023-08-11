Configuration
=============

The following provides an overview of all available flags in the configuration file to use the high level API and run jobs.

Mandatory flags
---------------

- `type` = "CollisionEnergyAlignment", "SpectralLibraryGeneration" or "Rescoring"
- `tag` = "tmt", "tmtpro", "itraq4" or "itraq8"; default is ""
- `models`
   - `intensity` = intensity model
   - `irt` = irt model
- `prediction_server` = server for obtaining peptide property predictions
- `ssl` = Use ssl when making requests to the prediction server, can be true or false; default = true
- `output` = path to the output folder; if not provided the current working directory will be used.

For spectral library generation and rescoring
---------------------------------------------

- `inputs`
   - `search_results` = path to the file containing the search results
   - `search_results_type` = the tool used to produce the search results, can be "Maxquant", "Msfragger", "Mascot" or "Internal"; default = "Maxquant"
   - `spectra` = path to a folder or a single file containing mass spectrometry results (raw or mzml files)
   - `spectra_type` = "raw" or "mzml"; default = "raw"
- `numThreads` = number of raw files processed in parallel processes; default = 1
- `thermoExe` = path to ThermoRawFileParser executable; default "ThermoRawFileParser.exe"
- `massTolerance` = mass tolerance value defining the allowed tolerance between theoretical and experimentally observered fragment mass during peak filtering and annotation. Default depends on the mass analyzer: 20 (FTMS), 40 (TOF), 0.35 (ITMS)
- `unitMassTolerance` = unit for the mass tolerance, either "da" or "ppm". Default is da (mass analyzer is ITMS) and ppm (mass analyzer is FTMS or TOF)

For spectral library generation only
------------------------------------

- `inputs`
   - `library_input` = path to the FASTA or peptides file
   - `library_input_type` = library input type: "fasta" or "peptides"
- `outputFormat` = "spectronaut" or "msp"

For in-silico digestion (spectral library generation) only
----------------------------------------------------------

- `fastaDigestOptions`
   - `fragmentation` = fragmentation method: "HCD" or "CID"
   - `digestion` = digestion mode: "full", "semi" or None; default = "full"
   - `cleavages` = number of allowed missed cleavages used in the search engine; default = 2
   - `minLength` = minimum peptide length allowed used in the search engine; default = 7
   - `maxLength` = maximum peptide length allowed used in the search engine; default = 60
   - `enzyme` = type of enzyme used in the search engine; default = "trypsin"
   - `specialAas` = special amino acids for decoy generation; default = "KR"
   - `db` = "target", "decoy" or "concat"; default = "concat"

For rescoring
-------------

- `fdr_estimation_method` = method used for FDR estimation on PSM and peptide level: "percolator" or "mokapot"; default = "mokapot"
- `allFeatures`` = True if all features should be used for FDR estimation; default = False
- `regressionMethod` = regression method for curve fitting (mapping from predicted iRT values to experimental retention times): "lowess", "spline" or "logistic"; default = "spline"
