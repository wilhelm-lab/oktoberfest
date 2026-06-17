from __future__ import annotations

"""Static seed catalogs for /meta endpoints.

Model lists SHOULD be fetched from Koina when feasible; this module provides
the static fallback used locally or when the prediction server is unreachable.
"""

INTENSITY_MODELS = [
    "Prosit_2019_intensity",
    "Prosit_2020_intensity_HCD",
    "Prosit_2020_intensity_CID",
    "Prosit_2020_intensity_TMT",
    "Prosit_2023_intensity_XL_CMS2",
    "Prosit_2024_intensity_XL_NMS2",
    "AlphaPept_ms2_generic",
]

IRT_MODELS = [
    "Prosit_2019_irt",
    "Prosit_2020_irt_TMT",
    "Deeplc_hela_hf",
    "AlphaPept_rt_generic",
    "",  # none / disable
]

SEARCH_ENGINES = ["Maxquant", "Msfragger", "Mascot", "Sage", "OpenMS", "Xisearch"]

SPECTRA_TYPES = ["raw", "mzml", "d", "hdf"]

LIBRARY_FORMATS = ["msp", "spectronaut", "dlib"]

ENZYMES = ["trypsin", "trypsinp", "lys-c", "chymotrypsin", "glu-c", "arg-c", "asp-n"]

TAGS = ["", "tmt", "tmtpro", "itraq4", "itraq8"]

INSTRUMENT_TYPES = ["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"]

FRAGMENTATION_METHODS = ["HCD", "ECD", "EID", "ETciD", "UVPD"]


DEFAULTS = {
    "Rescoring": {
        "type": "Rescoring",
        "tag": "",
        "models": {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"},
        "prediction_server": "koina.wilhelmlab.org:443",
        "ssl": True,
        "numThreads": 1,
        "inputs": {
            "search_results_type": "Maxquant",
            "spectra_type": "raw",
            "instrument_type": "QE",
        },
        "fdr_estimation_method": "mokapot",
        "add_feature_cols": "none",
        "regressionMethod": "spline",
        "massTolerance": 20,
        "unitMassTolerance": "ppm",
        "thermoExe": "ThermoRawFileParser.exe",
        "ce_alignment_options": {"ce_range": [19, 50], "use_ransac_model": False},
        "quantification": False,
    },
    "CollisionEnergyCalibration": {
        "type": "CollisionEnergyCalibration",
        "tag": "",
        "models": {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"},
        "prediction_server": "koina.wilhelmlab.org:443",
        "ssl": True,
        "numThreads": 1,
        "inputs": {
            "search_results_type": "Maxquant",
            "spectra_type": "raw",
            "instrument_type": "QE",
        },
        "massTolerance": 20,
        "unitMassTolerance": "ppm",
        "thermoExe": "ThermoRawFileParser.exe",
        "ce_alignment_options": {"ce_range": [19, 50], "use_ransac_model": False},
    },
    "SpectralLibraryGeneration": {
        "type": "SpectralLibraryGeneration",
        "tag": "",
        "models": {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"},
        "prediction_server": "koina.wilhelmlab.org:443",
        "ssl": True,
        "numThreads": 1,
        "inputs": {
            "library_input_type": "fasta",
            "instrument_type": "QE",
        },
        "spectralLibraryOptions": {
            "fragmentation": "HCD",
            "collisionEnergy": 30,
            "precursorCharge": [2, 3],
            "minIntensity": 5e-4,
            "nrOx": 1,
            "batchsize": 10000,
            "format": "msp",
        },
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
    },
}
