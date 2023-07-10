import json
import logging
from pathlib import Path
from typing import Union

logger = logging.getLogger(__name__)


class Config:
    """Read config file and get information from it."""

    def __init__(self):
        """Initialize config file data."""
        self.data = {}

    def read(self, config_path: Union[str, Path]):
        """
        Read config file.

        :param config_path: path to config file as a string
        """
        logger.info(f"Reading configuration from {config_path}")
        if isinstance(config_path, str):
            config_path = Path(config_path)
        with open(config_path) as f:
            self.data = json.load(f)

    @property
    def num_threads(self) -> int:
        """Get the number of threads from the config file; if not specified return 1."""
        if "numThreads" in self.data:
            return self.data["numThreads"]
        else:
            return 1

    @property
    def fasta(self) -> str:
        """Get path to fasta file from the config file."""
        if "fileUploads" in self.data:
            return self.data["fileUploads"]["fasta"]
        elif "uploads" in self.data:
            return self.data["uploads"]["FASTA"]
        else:
            raise ValueError("No fasta file specified in config file")

    @property
    def prediction_server(self) -> str:
        """Get prosit server from the config file."""
        return self.data["prediction_server"]

    @property
    def ssl(self) -> bool:
        """Get ssl flag for prediction server."""
        return self.data.get("ssl", True)

    @property
    def models(self) -> dict:
        """Get intensity, IRT, and proteotypicity models from the config file."""
        if "models" in self.data:
            return self.data["models"]
        if "selectedIntensityModel" in self.data:
            return {
                "selectedIntensityModel": self.data["selectedIntensityModel"],
                "selectedIRTModel": self.data["selectedIRTModel"],
            }
        else:
            return self.data["models"]

    @property
    def fdr_estimation_method(self) -> str:
        """Get peptide detection method from the config file (mokapot or percolator)."""
        if "fdr_estimation_method" in self.data:
            return self.data["fdr_estimation_method"].lower()
        else:
            return "mokapot"

    @property
    def tag(self) -> str:
        """Get tag from the config file; if not specified return "tmt"."""
        if "tag" in self.data:
            return self.data["tag"].lower()
        else:
            return "tmt"

    @property
    def all_features(self) -> bool:
        """Get allFeatures flag (decides whether all features should be used by the percolator)."""
        if "allFeatures" in self.data:
            return self.data["allFeatures"]
        else:
            return False

    @property
    def curve_fitting_method(self) -> str:
        """Get regressionMethod flag (regression method for curve fitting: lowess, spline, or logistic). \
        If not specified, lowess is applied."""
        if "regressionMethod" in self.data:
            return self.data["regressionMethod"].lower()
        else:
            return "lowess"

    @property
    def job_type(self) -> str:
        """Get jobType flag (CollisionEnergyAlignment, SpectralLibraryGeneration or Rescoring) from the config file."""
        if "jobType" in self.data:
            return self.data["jobType"]
        elif "type" in self.data:
            return self.data["type"]
        else:
            raise ValueError("No job type specified in config file")

    @property
    def raw_type(self) -> str:
        """Get raw type (thermo or mzml) from the config file."""
        if "fileUploads" in self.data:
            return self.data["fileUploads"]["raw_type"].lower()
        else:
            return "thermo"

    @property
    def search_type(self) -> str:
        """Get search type (Maxquant, Msfragger, Mascot or Internal) from the config file."""
        if "fileUploads" in self.data:
            return self.data["fileUploads"]["search_type"].lower()
        else:
            return "maxquant"

    @property
    def output_format(self) -> str:
        """Get output format from the config file."""
        if "outputFormat" in self.data:
            return self.data["outputFormat"].lower()
        else:
            return ""

    @property
    def search_path(self) -> str:
        """Get search path from the config file."""
        if "searchPath" in self.data:
            return self.data["searchPath"]
        else:
            return ""

    @property
    def fragmentation(self) -> str:
        """Get fragmentation method from the config file (HCD or CID)."""
        if "fragmentation" in self.data["fastaDigestOptions"]:
            return self.data["fastaDigestOptions"]["fragmentation"].upper()
        else:
            return ""

    @property
    def digestion(self) -> str:
        """Get digestion mode (full, semi or none)."""
        if "digestion" in self.data["fastaDigestOptions"]:
            return self.data["fastaDigestOptions"]["digestion"].lower()
        else:
            return "full"

    @property
    def cleavages(self) -> int:
        """Get number of allowed missed cleavages used in the search engine."""
        if "missedCleavages" in self.data["fastaDigestOptions"]:
            return self.data["fastaDigestOptions"]["missedCleavages"]
        else:
            return 2

    @property
    def min_length(self) -> int:
        """Get minimum peptide length allowed used in the search engine."""
        if "minLength" in self.data["fastaDigestOptions"]:
            return self.data["fastaDigestOptions"]["minLength"]
        else:
            return 7

    @property
    def max_length(self) -> int:
        """Get maximum peptide length allowed used in the search engine."""
        if "maxLength" in self.data["fastaDigestOptions"]:
            return self.data["fastaDigestOptions"]["maxLength"]
        else:
            return 60

    @property
    def enzyme(self) -> str:
        """Get type of enzyme used."""
        if "enzyme" in self.data["fastaDigestOptions"]:
            return self.data["fastaDigestOptions"]["enzyme"].lower()
        else:
            return "trypsin"

    @property
    def special_aas(self) -> str:
        """Get special amino acids used by MaxQuant for decoy generation."""
        if "specialAas" in self.data["fastaDigestOptions"]:
            return self.data["fastaDigestOptions"]["specialAas"].upper()
        else:
            return "KR"

    @property
    def thermo_exe(self) -> Path:
        """Get the path to the ThermoRawFileParser executable. Returns "ThermoRawFileParser.exe" if not found."""
        return Path(self.data.get("thermoExe", "ThermoRawFileParser.exe"))

    @property
    def db(self) -> str:
        """Target, decoy or concat (relevant if fasta file provided)."""
        if "db" in self.data["fastaDigestOptions"]:
            return self.data["fastaDigestOptions"]["db"].lower()
        else:
            return "concat"
