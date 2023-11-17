RESCORING = {
    "type": "Rescoring",
    "tag": "",
    "inputs": {"search_results": "msms.txt", "search_results_type": "Maxquant", "spectra": "./", "spectra_type": "raw"},
    "output": "./out",
    "models": {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"},
    "prediction_server": "koina.proteomicsdb.org:443",
    "ssl": True,
    "thermoExe": "ThermoRawFileParser.exe",
    "numThreads": 1,
    "fdr_estimation_method": "mokapot",
    "regressionMethod": "spline",
    "allFeatures": False,
    "massTolerance": 20,
    "unitMassTolerance": "ppm",
    "ce_alignment_options": {
        "ce_range": (19, 50),
        "use_ransac_model": False,
    },
}

CECALIB = {
    "type": "CollisionEnergyCalibration",
    "tag": "",
    "inputs": {"search_results": "msms.txt", "search_results_type": "Maxquant", "spectra": "./", "spectra_type": "raw"},
    "output": "./out",
    "models": {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"},
    "prediction_server": "koina.proteomicsdb.org:443",
    "ssl": True,
    "thermoExe": "ThermoRawFileParser.exe",
    "numThreads": 1,
    "massTolerance": 20,
    "unitMassTolerance": "ppm",
    "ce_alignment_options": {
        "ce_range": (19, 50),
        "use_ransac_model": False,
    },
}

LIBGEN = {
    "type": "SpectralLibraryGeneration",
    "tag": "",
    "inputs": {
        "search_results": "msms.txt",
        "search_results_type": "Maxquant",
        "library_input": "uniprot.fasta",
        "library_input_type": "fasta",
    },
    "output": "./out",
    "models": {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"},
    "outputFormat": "spectronaut",
    "prediction_server": "koina.proteomicsdb.org:443",
    "ssl": True,
    "fastaDigestOptions": {
        "fragmentation": "HCD",
        "digestion": "full",
        "missedCleavages": 2,
        "minLength": 7,
        "maxLength": 60,
        "enzyme": "trypsin",
        "specialAas": "KR",
        "db": "concat",
    },
}
