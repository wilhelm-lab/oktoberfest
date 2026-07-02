import json
import logging
import os
from pathlib import Path
from sys import platform
from typing import Optional, Union

from koinapy.grpc import Koina

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
        self.base_path = config_path.parent

    ###########################
    # these are parameters    #
    ###########################

    @property
    def num_threads(self) -> int:
        """Get the number of threads from the config file; if not specified return 1."""
        return self.data.get("numThreads", 1)

    @property
    def prediction_server(self) -> str:
        """Get prosit server from the config file."""
        return self.data.get("prediction_server", "koina.wilhelmlab.org:443")

    @property
    def ssl(self) -> bool:
        """Get ssl flag for prediction server."""
        return self.data.get("ssl", True)

    @property
    def models(self) -> dict:
        """Get intensity, IRT, and proteotypicity models from the config file."""
        return self.data["models"]

    @property
    def tag(self) -> str:
        """Get tag from the config file; if not specified return ""."""
        return self.data.get("tag", "").lower()

    @property
    def job_type(self) -> str:
        """Get job type (CollisionEnergyAlignment, SpectralLibraryGeneration or Rescoring) from the config file."""
        job_type = self.data.get("type")
        if job_type is None:
            raise ValueError("No job type specified in config file.")
        return job_type

    @property
    def quantification(self) -> bool:
        """Get quantification flag for performing quantification using picked-group-fdr."""
        return self.data.get("quantification", False)

    @property
    def mass_tolerance(self) -> Optional[float]:
        """Get mass tolerance value from the config file with which to caluculate the min and max mass values."""
        return self.data.get("massTolerance", None)

    @property
    def p_window(self) -> Optional[float]:
        """Windows size for precursor peak removal."""
        return self.data.get("p_window", 0.0)

    @property
    def fragmentation_method(self) -> str:
        """Get fragmentation method from config file."""
        return self.data.get("fragmentation_method", "HCD")

    @property
    def ion_types(self) -> list[str]:
        """
        Return the fragment ion types used for fragment annotation and Percolator feature calculation.

        Specify fragment ion types as a concatenated string in alphabetical order,
        for example, ``"aAbcCxXyzZ"``. The default is ``"by"``.

        :returns: A list representing the fragment ion types.
        """
        raw = self.data.get("ion_types", "yb")
        # normalize to string
        seq = "".join(map(str, raw)) if isinstance(raw, list) else str(raw)

        # priority of ion in order
        allowed = "ybaAcCxXzZ"
        seen = set()
        filtered = []
        for ch in seq:
            if ch in allowed and ch not in seen:
                seen.add(ch)
                filtered.append(ch)

        filtered.sort(key=lambda ch: allowed.index(ch))
        return filtered if filtered else ["y", "b"]

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

    @property
    def mirror_plots(self) -> dict[str, list[int]]:
        """
        Get the raw files and scan numbers for which to generate mirror plots.

        This function returns a dictionary where the keys are raw file names,
        and the values are lists of scan numbers.

        :return: Dictionary mapping raw file names to a list of scan numbers.
        """
        return {str(k): list(map(int, v)) for k, v in self.data.get("mirror_plots", {}).items()}

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
    def custom_modifications(self) -> dict[str, dict[str, list[Union[int, float, str]]]]:
        """Get the custom modification dictionary from the config file."""
        return self.inputs.get("custom_modifications", {})

    @property
    def static_mods(self) -> dict[str, tuple[int, float, str]]:
        """
        Get the custom static modification information.

        This function returs a dictionary with custom mod identifiers (keys), and a tuple of
        (UNIMOD Id, modification mass delta, and optional neutral losses) (values).
        :return: dictionary mapping static mod identifiers to a tuple containing unimod id, mass,
            and optionally neutral losses
        """
        mod_items = self.custom_modifications.get("static_mods", {}).items()
        return {str(k): (int(v[0]), float(v[1]), str(v[2]) if len(v) >= 3 else "") for k, v in mod_items}

    @property
    def var_mods(self) -> dict[str, tuple[int, float, str]]:
        """
        Get the custom variable modification information.

        This function returs a dictionary with custom mod identifiers (keys), and a tuple of
        (UNIMOD Id, modification mass delta, and optional neutral losses) (values).
        :return: dictionary mapping var mod identifiers to a tuple containing unimod id, mass,
            and optionally neutral losses
        """
        mod_items = self.custom_modifications.get("var_mods", {}).items()
        return {str(k): (int(v[0]), float(v[1]), str(v[2]) if len(v) >= 3 else "") for k, v in mod_items}

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

    @property
    def instrument_type(self) -> Optional[str]:
        """Get type of mass spectrometer from the config file (superseeds value read from from mzML)."""
        _instrument_type = self.inputs.get("instrument_type")
        if _instrument_type is None:
            return None
        return _instrument_type.upper()

    @property
    def dlomix_pipeline(self):
        """Load and cache DLOmix InferencePipeline if dlomix_intensity is configured."""
        if not hasattr(self, '_dlomix_pipeline'):
            dlomix_path = self.models.get("dlomix_intensity")
            if dlomix_path:
                try:
                    # Disable GPU before importing DLOmix to avoid CUDA initialization errors in workers
                    os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
                    from dlomix.pipelines.predictor import InferencePipeline
                    pipeline_path = self.base_path / Path(dlomix_path) if not Path(dlomix_path).is_absolute() else Path(dlomix_path)
                    self._dlomix_pipeline = InferencePipeline.load(str(pipeline_path))
                except Exception as e:
                    logger.warning(f"Failed to load DLOmix pipeline from {dlomix_path}: {e}")
                    self._dlomix_pipeline = None
            else:
                self._dlomix_pipeline = None
        return self._dlomix_pipeline

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
    def ce_range(self) -> tuple[int, int]:
        """Get the min and max boundaries for the CE to be used for alignment."""
        min_ce, max_ce = self.ce_alignment_options.get("ce_range", (18, 50))
        return int(min_ce), int(max_ce)

    @property
    def use_ransac_model(self) -> bool:
        """Get whether or not to perform ce calibration using a ransac model."""
        return self.ce_alignment_options.get("use_ransac_model", False)

    ###############################
    # these are rescoring options #
    ###############################

    @property
    def use_feature_cols(self) -> Union[str, list]:
        """Get additional columns ("all" for all columns or list with column names) from the config file."""
        return self.data.get("add_feature_cols", "none")

    @property
    def all_features(self) -> bool:
        """Get allFeatures flag (decides whether all features should be used as input for the chosen fdr estimation method)."""
        return self.data.get("allFeatures", False)

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
    def fdr_estimation_method(self) -> str:
        """Get peptide detection method from the config file (mokapot or percolator). Default is percolator."""
        return self.data.get("fdr_estimation_method", "percolator").lower()

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
    def speclib_generation_batch_size(self) -> int:
        """Get output format from the config file."""
        return self.spec_lib_options.get("batchsize", 10000)

    @property
    def min_intensity(self) -> float:
        """Get output format from the config file."""
        return self.spec_lib_options.get("minIntensity", 5e-4)

    @property
    def precursor_charge(self) -> list[int]:
        """Get output format from the config file."""
        return self.spec_lib_options.get("precursorCharge", [2, 3])

    @property
    def nr_ox(self) -> int:
        """Get the maximum number of oxidations allowed on M residues in peptides during spectral library generation."""
        return self.spec_lib_options.get("nrOx", 1)

    ###########################
    # these are PTM localization options #
    ###########################

    @property
    def ptm_localization(self) -> bool:
        """Get ptm localization flag from the config file."""
        return self.data.get("ptm_localization", False)

    @property
    def ptm_localization_options(self) -> dict:
        """Get ptm localization dictionary from the config file."""
        return self.data.get("ptmLocalizationOptions", {})

    @property
    def ptm_unimod_id(self) -> int:
        """Get the unimod id required for localization."""
        unimod_id = self.ptm_localization_options.get("unimod_id", 0)
        return unimod_id

    @property
    def ptm_possible_sites(self) -> list:
        """Get which sites ptm can exist on."""
        return self.ptm_localization_options.get("possible_sites", [])

    @property
    def ptm_use_neutral_loss(self) -> bool:
        """Get neutral loss flag to indicate whether to add a score for this or not."""
        return self.ptm_localization_options.get("neutral_loss", False)

    ########################
    # functions start here #
    ########################

    def check(self):
        """Validate the configuration."""
        self._check_tmt()

        if self.job_type == "SpectralLibraryGeneration":
            self._check_for_speclib()

        # Skip alphapept and Koina checks if using DLOmix
        if not self.models.get("dlomix_intensity"):
            if "alphapept" in self.models["intensity"].lower():
                self._check_for_alphapept()

            self._check_koina_model_availability()

        if self.quantification:
            self._check_quantification()
            self._check_fasta()

    def _check_koina_model_availability(self):
        """Check if Koina model is available."""
        # This will give error automaticly in Koina if model is not available on the server.
        # Koina has function called "_is_model_ready" that checks if model is available
        _ = Koina(model_name=self.models["intensity"])

    def _check_tmt(self):
        # Skip TMT check if using DLOmix (no intensity model requirement)
        if self.models.get("dlomix_intensity"):
            return

        int_model = self.models["intensity"].lower()
        irt_model = self.models["irt"].lower()
        if self.tag == "":
            if "tmt" in int_model:
                raise AssertionError(
                    f"You requested the intensity model {self.models['intensity']} but provided no tag. Please check."
                )
            if "tmt" in irt_model:
                raise AssertionError(
                    f"You requested the irt model {self.models['irt']} but provided no tag. Please check."
                )
        else:
            if ("tmt" not in int_model) and ("ptm" not in int_model) and ("alphapept" not in int_model):
                raise AssertionError(
                    f"You specified the tag {self.tag} but the chosen intensity model {self.models['intensity']} is incompatible. "
                    "Please check and use a TMT model instead."
                )
            if ("tmt" not in irt_model) and ("ptm" not in irt_model) and ("alphapept" not in irt_model):
                raise AssertionError(
                    f"You specified the tag {self.tag} but the chosen irt model {self.models['irt']} is incompatible."
                    " Please check and use a TMT model instead."
                )

    def _check_for_alphapept(self):
        instrument_type = self.instrument_type
        valid_alphapept_instrument_types = ["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"]
        if instrument_type is not None and instrument_type not in valid_alphapept_instrument_types:
            raise ValueError(
                f"The chosen intensity model {self.models['intensity']} does not support the specified instrument type "
                f"{instrument_type}. Either let Oktoberfest read the instrument type from the mzML file, or provide one "
                f"of {valid_alphapept_instrument_types}."
            )

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
        if "alphapept" in self.models["intensity"].lower():
            self._check_for_alphapept()

    def _find_file_in_subd(self, directory: Path, filename: str) -> bool:
        """Return True if a file with the given name exists anywhere under directory."""
        return any(path.name == filename for path in directory.rglob(filename))

    def _check_quantification(self):
        if Path(self.search_results).is_file():
            path_stem = Path(self.search_results).parent
        else:
            path_stem = Path(self.search_results)

        if self.search_results_type == "maxquant" and not Path(path_stem / "evidence.txt").is_file():
            raise AssertionError(
                f"You specified the search results as {self.search_results_type} but evidence.txt is not available "
                f"at {path_stem / 'evidence.txt'}."
            )
        elif self.search_results_type == "sage":
            if not Path(path_stem / "results.sage.tsv").is_file():
                raise AssertionError(
                    f"You specified the search results as {self.search_results_type} for quantification, but "
                    f"results.sage.tsv is not available at {path_stem / 'results.sage.tsv'}."
                )
            elif not Path(path_stem / "lfq.tsv").is_file():
                raise AssertionError(
                    f"You specified the search results as {self.search_results_type} for quantification, but "
                    f"lfq.tsv is not available at {path_stem / 'lfq.tsv'}."
                )
        elif self.search_results_type == "msfragger":
            if not self._find_file_in_subd(path_stem, "psm.tsv"):
                raise AssertionError(
                    f"You specified the search results as {self.search_results_type} for quantification, but "
                    "no psm.tsv files could be found in subdirectories."
                )
            elif not Path(path_stem / "combined_ion.tsv").is_file():
                raise AssertionError(
                    f"You specified the search results as {self.search_results_type} for quantification, but "
                    f"combined_ion.tsv is not available  at {path_stem / 'combined_ion.tsv'}."
                )

    def _check_fasta(self):
        if not self.library_input_type.lower() == "fasta":
            raise AssertionError(
                f"The specified library input type is set to {self.library_input_type}. "
                "For quantification a fasta file is needed."
            )

    def check_multifrag(self) -> bool:
        """Check if rescoring will be done on multifrag options."""
        intensity_model = self.models.get("intensity", "").lower()
        return "multifrag" in intensity_model

    def custom_to_unimod(self) -> dict[str, int]:
        """
        Parse modifications to dict with custom identifier and UNIMOD integer for internal processing.

        :return: a dictionary mapping custom mod identifiers (keys) to the unimod id (values).
        """
        custom_to_unimod = {}
        for k, v in self.var_mods.items():
            custom_to_unimod[str(k)] = int(v[0])
        for k, v in self.static_mods.items():
            custom_to_unimod[str(k)] = int(v[0])
        return custom_to_unimod

    def unimod_to_mass(self) -> dict[str, float]:
        """
        Map UNIMOD Id to its mass for all static and variable modifications.

        This function maps the UNIMOD Id to its corresponding mass for each custom modifiction
        provided in the static and variable modifications.

        :return: a dictionary mapping the UNIMOD Ids (keys) to the mass(value) of a given modification
        """
        unimod_to_mass = {}
        for unimod_id, mass, _ in self.var_mods.values():
            unimod_to_mass[f"[UNIMOD:{unimod_id}]"] = mass
        for unimod_id, mass, _ in self.static_mods.values():
            unimod_to_mass[f"[UNIMOD:{unimod_id}]"] = mass
        return unimod_to_mass
