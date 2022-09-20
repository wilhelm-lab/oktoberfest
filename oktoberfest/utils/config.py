import json
import logging
from typing import Dict

logger = logging.getLogger(__name__)


class Config:
    """Read config file and get information from it."""

    data: Dict

    def __init__(self):
        """Initialize config file data."""
        self.data = {}

    def read(self, config_path: str):
        """
        Read config file.

        :param config_path: path to config file as a string
        """
        logger.info(f"Reading configuration from {config_path}")
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
    def prosit_server(self) -> str:
        """Get prosit server from the config file."""
        return self.data["prosit_server"]

    @property
    def models(self) -> dict:
        """Get intensity, IRT, and proteotypicity models from the config file."""
        if "models" in self.data:
            return self.data["models"]
        if "selectedIntensityModel" in self.data:
            return {
                "selectedIntensityModel": self.data["selectedIntensityModel"],
                "selectedIRTModel": self.data["selectedIRTModel"],
                "selectedProteotypicityModel": self.data["selectedProteotypicityModel"],
            }
        else:
            return self.data["models"]

    @property
    def tag(self) -> str:
        """Get tag from the config file; if not specified return "tmt"."""
        if "tag" in self.data:
            return self.data["tag"]
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
    def job_type(self) -> str:
        """Get jobType flag (CollisionEnergyAlignment, MaxQuantRescoring or SpectralLibraryGeneration) from the config file."""
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
            return self.data["fileUploads"]["raw_type"]
        else:
            return "thermo"

    @property
    def search_type(self) -> str:
        """Get search type (maxquant or internal) from the config file."""
        if "fileUploads" in self.data:
            return self.data["fileUploads"]["search_type"]
        else:
            return "maxquant"

    @property
    def output_format(self) -> str:
        """Get output format from the config file."""
        if "outputFormat" in self.data:
            return self.data["outputFormat"]
        else:
            return ""
