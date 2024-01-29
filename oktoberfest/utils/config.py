import json
import logging
from pathlib import Path
from sys import platform
from typing import List, Optional, Tuple, Union

logger = logging.getLogger(__name__)


class Config:
    """Read config file and get information from it."""

    @property
    def num_threads(self) -> int:
        """Get the number of threads from the config file; if not specified return 1."""
        return self.data.get("numThreads", 1)

    @property
    def prediction_server(self) -> str:
        """Get prosit server from the config file."""
        return self.data["prediction_server"]

    @property
    def ssl(self) -> bool:
        """Get ssl flag for prediction server."""
        return self.data.get("ssl", False)

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
        """Get tag from the config file; if not specified return ""."""
        if "tag" in self.data:
            return self.data["tag"].lower()
        else:
            return ""

    @property
    def all_features(self) -> bool:
        """Get allFeatures flag (decides whether all features should be used as input for the chosen fdr estimation method)."""
        if "allFeatures" in self.data:
            return self.data["allFeatures"]
        else:
            return False

    @property
    def curve_fitting_method(self) -> str:
        """
        Get regressionMethod flag.

        Reads the regressionMethod flag that is used to determine the method for retention time alignment.
        The supported flags are "lowess", "spline", and "logistic".
        If not provided in the config file, returns "spline" by default.

        :return: a lowercase string representation of the regression method.
        """
        return self.data.get("regressionMethod", "spline").lower()

    @property
    def job_type(self) -> str:
        """Get jobType flag (CollisionEnergyAlignment, SpectralLibraryGeneration or Rescoring) from the config file."""
        if "jobType" in self.data:
            return self.data["jobType"]
        elif "type" in self.data:
            return self.data["type"]
        else:
            raise ValueError("No job type specified in config file.")

    @property
    def mass_tolerance(self) -> Optional[float]:
        """Get mass tolerance value from the config file with which to caluculate the min and max mass values."""
        return self.data.get("massTolerance", None)

    @property
    def unit_mass_tolerance(self) -> Optional[str]:
        """Get unit for the mass tolerance from the config file (da or ppm)."""
        return self.data.get("unitMassTolerance", None)

    @property
    def output(self) -> Path:
        """Get path to the output directory from the config file."""
        # check if output is absolute and if not, append with the directory
        # of the config file to resolve paths relative to the config location. This is
        # required to make sure it always works no matter from which working directory
        # oktoberfest is executed in case a relative path is provided, which would
        # otherwise be relative to the working directory.
        return self.base_path / Path(self.data.get("output", "./"))

    @property
    def thermo_exe(self) -> Path:
        """Get the path to the ThermoRawFileParser executable. Returns "ThermoRawFileParser.exe" if not found."""

        def default_thermo():
            if "linux" in platform or platform == "darwin":
                return "/opt/compomics/ThermoRawFileParser.exe"
            return "ThermoRawFileParser.exe"

        return Path(self.data.get("thermoExe", default_thermo()))

    ###########################
    # these are input options #
    ###########################

    @property
    def inputs(self) -> dict:
        """Get inputs dictionary from the config file."""
        return self.data.get("inputs", {})

    @property
    def search_results(self) -> Path:
        """Get path to the search results file from the config file."""
        search_results_path = self.inputs.get("search_results")
        if search_results_path is not None:
            # check if search_results_path is absolute and if not, append with the directory
            # of the config file to resolve paths relative to the config location. This is
            # required to make sure it always works no matter from which working directory
            # oktoberfest is executed in case a relative path is provided, which would
            # otherwise be relative to the working directory.
            return self.base_path / search_results_path
        else:
            raise ValueError("No path to a msms.txt specified in config file.")

    @property
    def search_results_type(self) -> str:
        """Get search type (Maxquant, Msfragger, Mascot or Internal) from the config file."""
        return self.inputs.get("search_results_type", "maxquant").lower()

    @property
    def spectra(self) -> Path:
        """Get path to spectra files from the config file."""
        # check if spectra is absolute and if not, append with the directory
        # of the config file to resolve paths relative to the config location. This is
        # required to make sure it always works no matter from which working directory
        # oktoberfest is executed in case a relative path is provided, which would
        # otherwise be relative to the working directory.
        return self.base_path / Path(self.inputs.get("spectra", "./"))

    @property
    def spectra_type(self) -> str:
        """Get spectra type (.raw, .mzml, .d, .hdf) from the config file."""
        return self.inputs.get("spectra_type", "raw").lower()

    @property
    def library_input(self) -> Path:
        """Get path to library input file (fasta or peptides) from the config file."""
        # check if library_input is absolute and if not, append with the directory
        # of the config file to resolve paths relative to the config location. This is
        # required to make sure it always works no matter from which working directory
        # oktoberfest is executed in case a relative path is provided, which would
        # otherwise be relative to the working directory.
        library_input_path = self.inputs.get("library_input")
        if library_input_path is not None:
            return self.base_path / Path(library_input_path)
        else:
            raise ValueError("No fasta or peptides file specified using library_input in config file.")

    @property
    def library_input_type(self) -> str:
        """Get library input file type (fasta or peptides) from the config file."""
        library_input_type = self.inputs.get("library_input_type")
        if library_input_type is not None:
            return library_input_type.lower()
        else:
            raise ValueError("No library input file type (fasta or peptides) specified in config file.")

    #####################################
    # these are fasta digestion options #
    #####################################

    @property
    def fasta_digest_options(self) -> dict:
        """Get fastaDigestOptions dictionary from the config file."""
        return self.data.get("fastaDigestOptions", {})

    @property
    def digestion(self) -> str:
        """Get digestion mode (full, semi or none)."""
        return self.fasta_digest_options.get("digestion", "full").lower()

    @property
    def missed_cleavages(self) -> int:
        """Get number of allowed missed cleavages used in the search engine."""
        return self.fasta_digest_options.get("missedCleavages", 2)

    @property
    def min_length(self) -> int:
        """Get minimum peptide length allowed used in the search engine."""
        return self.fasta_digest_options.get("minLength", 7)

    @property
    def max_length(self) -> int:
        """Get maximum peptide length allowed used in the search engine."""
        return self.fasta_digest_options.get("maxLength", 60)

    @property
    def enzyme(self) -> str:
        """Get type of enzyme used."""
        return self.fasta_digest_options.get("enzyme", "trypsin").lower()

    @property
    def special_aas(self) -> str:
        """Get special amino acids used by MaxQuant for decoy generation."""
        return self.fasta_digest_options.get("specialAas", "KR").upper()

    @property
    def db(self) -> str:
        """Target, decoy or concat (relevant if fasta file provided)."""
        return self.fasta_digest_options.get("db", "concat").lower()

    ##################################
    # these are ce alignment options #
    ##################################

    @property
    def ce_alignment_options(self) -> dict:
        """Get the ce_alignment dictionary from the config."""
        return self.data.get("ce_alignment_options", {})

    @property
    def ce_range(self) -> Tuple[int, int]:
        """Get the min and max boundaries for the CE to be used for alignment."""
        min_ce, max_ce = self.ce_alignment_options.get("ce_range", (18, 50))
        return int(min_ce), int(max_ce)

    @property
    def use_ransac_model(self) -> bool:
        """Get whether or not to perform ce calibration using a ransac model."""
        return self.ce_alignment_options.get("use_ransac_model", False)

    ######################################
    # these are spectral library options #
    ######################################

    @property
    def spec_lib_options(self) -> dict:
        """Get inputs dictionary from the config file."""
        return self.data.get("spectralLibraryOptions", {})

    @property
    def output_format(self) -> str:
        """Get output format from the config file."""
        return self.spec_lib_options.get("format", "msp").lower()

    @property
    def fragmentation(self) -> str:
        """Get output format from the config file."""
        return self.spec_lib_options.get("fragmentation", "").lower()

    @property
    def collision_energy(self) -> int:
        """Get output format from the config file."""
        return self.spec_lib_options.get("collisionEnergy", 30)

    @property
    def batch_size(self) -> int:
        """Get output format from the config file."""
        return self.spec_lib_options.get("batchsize", 10000)

    @property
    def min_intensity(self) -> float:
        """Get output format from the config file."""
        return self.spec_lib_options.get("minIntensity", 5e-4)

    @property
    def precursor_charge(self) -> List[int]:
        """Get output format from the config file."""
        return self.spec_lib_options.get("precursorCharge", [2, 3])

    ########################
    # functions start here #
    ########################

    def check(self):
        """Validate the configuration."""
        # check tmt tag and models
        if self.tag == "":
            if "tmt" in self.models["intensity"].lower():
                raise AssertionError(
                    f"You requested the intensity model {self.models['intensity']} but provided no tag. Please check."
                )
            if "tmt" in self.models["irt"].lower():
                raise AssertionError(
                    f"You requested the irt model {self.models['irt']} but provided no tag. Please check."
                )
        else:
            if "tmt" not in self.models["intensity"].lower():
                raise AssertionError(
                    f"You specified the tag {self.tag} but the chosen intensity model {self.models['intensity']} is incompatible. "
                    "Please check and use a TMT model instead."
                )
            if "tmt" not in self.models["irt"].lower():
                raise AssertionError(
                    f"You specified the tag {self.tag} but the chosen irt model {self.models['irt']} is incompatible."
                    " Please check and use a TMT model instead."
                )
        if self.job_type == "SpectralLibraryGeneration":
            self._check_for_speclib()

    def _check_for_speclib(self):
        if self.fragmentation == "hcd" and self.models["intensity"].lower().endswith("cid"):
            raise AssertionError(
                f"You requested the intensity model {self.models['intensity']} but want to create a spectral library for HCD data. "
                "Please check that the fragmentation method and the chosen model agree."
            )
        elif self.fragmentation == "cid" and self.models["intensity"].lower().endswith("hcd"):
            raise AssertionError(
                f"You requested the intensity model {self.models['intensity']} but want to create a spectral library for CID data. "
                "Please check that the fragmentation method and the chosen model agree."
            )
        elif self.fragmentation == "":
            model_end = self.models["intensity"].lower()[-3:]
            if model_end == "hcd" or model_end == "cid":
                logger.warning(
                    f"No fragmentation method was specified. Assuming {model_end} fragmentation based on chosen intensity "
                    f"model {self.models['intensity']}."
                )
            else:
                raise AssertionError(
                    f"You need to provide the fragmentation method when using the model {self.models['intensity']}."
                )

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
        self.base_path = config_path.parent
