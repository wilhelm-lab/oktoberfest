import json
import logging
import os
from pathlib import Path
from sys import platform
from typing import Optional, Union

from koinapy.grpc import Koina

logger = logging.getLogger(__name__)

# residues used as decoy sites next to the possible sites (decoy_selection "adjacent"): no possible phospho sites, and
# not P, D, E, K, R, H, C and M, which change the fragmentation or are frequently modified
ADJACENT_DECOY_RESIDUES = ["A", "G", "L", "V", "I", "F", "N", "Q", "W"]

# decoy residues used to estimate the false localization rate if ptmLocalizationOptions.decoy_residues is not given
DEFAULT_DECOY_RESIDUES = ["A", "G", "L", "V"]

# options of ptmLocalizationOptions (first value is the default)
LOCALIZATION_PROBABILITIES = ["pep", "score"]
DECOY_WEIGHTINGS = ["count", "distance"]
FLR_GROUPS = ["none", "residue"]
FLR_ESTIMATORS = ["hybrid", "decoy"]

# options of quantificationOptions (first value is the default)
QUANT_PRECURSOR_METHODS = ["search_engine", "flashlfq"]
QUANT_PROTEIN_METHODS = ["directlfq", "none"]
QUANT_NORMALIZATIONS = ["drift", "mode", "median", "none"]

# DLOmix pipelines keyed by resolved model path. Module-level (not on Config) so a parent-loaded
# pipeline is inherited by fork workers via copy-on-write, while Config stays picklable.
_DLOMIX_PIPELINE_CACHE: dict = {}


def _is_number(value) -> bool:
    """Return True for int and float values (bool is not a number here)."""
    return isinstance(value, (int, float)) and not isinstance(value, bool)


class Config:
    """Read config file and get information from it."""

    def __init__(self):
        """Initialize config file data."""
        self.data = {}
        self.data["models"] = {}
        self._run_original = False

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
        """Get the number of threads from the config file; if not specified return 1.

        Forced to 1 when a DLOmix model is configured, since DLOmix's TensorFlow-based
        pipeline is not safe to run with multiple worker threads/processes (see
        ``dlomix_pipeline`` and ``check_dlomix``).
        """
        if self.models.get("dlomix_intensity"):
            return 1
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
        """Get job type (GenerateTrainingData, CollisionEnergyAlignment, SpectralLibraryGeneration or Rescoring) from the config file."""
        job_type = self.data.get("type")
        if job_type is None:
            raise ValueError("No job type specified in config file.")
        return job_type

    @property
    def acquisition_method(self) -> str:
        """Get acquisition method from the config file; default is "DDA". Possible values are "DDA" or "DIA"."""
        return self.data.get("acquisitionMethod", "DDA").lower()

    @property
    def quantification(self) -> bool:
        """Get quantification flag for performing quantification (see quantificationOptions)."""
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
        """Return the DLOmix pipeline, loading it once per model path and caching it at module level.

        Loaded once in the parent (via check_dlomix) before the worker pools fork, so workers
        inherit the single load via copy-on-write instead of each reloading it. NOTE: this relies
        on TF tolerating fork-after-init; _load_dlomix_pipeline pins TF to single-threaded execution
        as a mitigation for the fork-after-init deadlock.
        """
        dlomix_path = self.models.get("dlomix_intensity")
        if not dlomix_path:
            return None
        key = str(self.base_path / dlomix_path if not Path(dlomix_path).is_absolute() else Path(dlomix_path))
        if key not in _DLOMIX_PIPELINE_CACHE:
            _DLOMIX_PIPELINE_CACHE[key] = self._load_dlomix_pipeline()
        return _DLOMIX_PIPELINE_CACHE[key]

    def _load_dlomix_pipeline(self):
        """Load the DLOmix InferencePipeline configured via ``dlomix_intensity``, or None."""
        dlomix_path = self.models.get("dlomix_intensity")
        if not dlomix_path:
            return None
        try:
            # Disable GPU before importing DLOmix to avoid CUDA initialization errors in workers
            os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
            from dlomix.pipelines.predictor import InferencePipeline

            pipeline_path = Path(dlomix_path)
            if not pipeline_path.is_absolute():
                pipeline_path = self.base_path / pipeline_path
            logger.info(f"Loading DLOmix pipeline from {pipeline_path} (pid {os.getpid()})")
            return InferencePipeline.load(str(pipeline_path))
        except ImportError as e:
            logger.warning(
                f"DLOmix is not installed, so the 'dlomix_intensity' model could not be loaded. "
                f"Install the optional extra with `pip install oktoberfest[dlomix]` (requires python >=3.11). ({e})"
            )
            return None
        except Exception as e:
            logger.warning(f"Failed to load DLOmix pipeline from {dlomix_path}: {e}")
            return None

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

    @property
    def ce_calibration_dir(self) -> Optional[Path]:
        """Get the results/ce_calibration folder of a previous run to reuse its calibrated CE per spectra file."""
        ce_calibration_dir = self.ce_alignment_options.get("ce_calibration_dir")
        if ce_calibration_dir is None:
            return None
        return self.base_path / Path(ce_calibration_dir)

    ####################################
    # these are quantification options #
    ####################################

    @property
    def quantification_options(self) -> dict:
        """Get the quantificationOptions dictionary from the config."""
        return self.data.get("quantificationOptions", {})

    @property
    def quant_precursor(self) -> str:
        """
        Get how precursors are quantified; default is "search_engine".

        "search_engine" updates the search engine's own quantification (MaxQuant, Sage, FragPipe) with the rescored PSMs
        using picked-group-FDR. "flashlfq" measures the MS1 signal of the rescored PSMs with FlashLFQ, independent of
        the search engine.
        """
        return str(self.quantification_options.get("precursor", QUANT_PRECURSOR_METHODS[0])).lower()

    @property
    def quant_protein(self) -> str:
        """Get how protein intensities are calculated from the precursors ("directlfq" or "none"; flashlfq only)."""
        return str(self.quantification_options.get("protein", QUANT_PROTEIN_METHODS[0])).lower()

    @property
    def quant_normalization(self) -> str:
        """
        Get the normalization of the precursor intensities between runs (flashlfq only); default is "drift".

        "drift" corrects every run by a curve along the retention time, the most frequent log2 deviation from the
        precursors' median over the runs; "mode" shifts every run by the most frequent log2 ratio to the first run,
        "median" by the median log2 ratio (pulled by precursors that change between samples), "none" keeps the FlashLFQ
        intensities.
        """
        return str(self.quantification_options.get("normalization", QUANT_NORMALIZATIONS[0])).lower()

    @property
    def quant_psm_fdr(self) -> float:
        """Get the PSM q-value cutoff of the PSMs that are quantified (flashlfq only); default is 0.01."""
        return self.quantification_options.get("psmFdr", 0.01)

    @property
    def quant_peptide_fdr(self) -> float:
        """Get the Percolator peptide q-value cutoff of the peptidoform table (flashlfq only); default is 0.01."""
        return self.quantification_options.get("peptideFdr", 0.01)

    @property
    def quant_protein_fdr(self) -> float:
        """Get the protein group q-value cutoff of the quantified protein groups (flashlfq only); default is 0.01."""
        return self.quantification_options.get("proteinFdr", 0.01)

    @property
    def quant_ppm_tolerance(self) -> float:
        """Get the MS1 mass tolerance in ppm used by FlashLFQ; default is 10."""
        return self.quantification_options.get("ppmTolerance", 10)

    @property
    def quant_ppm_filter(self) -> float | None:
        """
        Get the maximum distance of a peak's apex mass error to its run's median in ppm (flashlfq only); default is 5.

        Peaks further away are not used; None keeps all peaks.
        """
        return self.quantification_options.get("ppmFilter", 5)

    @property
    def flashlfq_exe(self) -> Path:
        """Get the FlashLFQ command line executable (CMD.dll of the FlashLFQ release); default is "CMD.dll"."""
        return self.base_path / Path(self.quantification_options.get("flashlfqExe", "CMD.dll"))

    @property
    def dotnet_exe(self) -> str:
        """Get the dotnet executable used to run FlashLFQ; default is "dotnet"."""
        return self.quantification_options.get("dotnetExe", "dotnet")

    @property
    def ms1_quantification(self) -> bool:
        """Return True if precursors are quantified on MS1 by Oktoberfest, which needs MS1 spectra."""
        return self.quantification and self.quant_precursor == "flashlfq"

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

    @property
    def percolator_subset_max_train(self) -> Optional[int]:
        """
        Get the maximum number of PSMs percolator trains on (percolator --subset-max-train), None for all PSMs.

        Only used for rescoring (original and rescore percolator runs, also with separate FDR for modified peptides),
        never for the PTM localization models. Default is None.
        """
        return self.data.get("percolatorSubsetMaxTrain")

    @property
    def chimeric_peak_subtraction(self) -> bool:
        """
        Get flag whether to remove the peaks of better candidates of the same spectrum before the rescore features.

        Within every spectrum with several candidates, the peaks matched by a candidate with a higher spectral angle
        are removed before the intensity features of the next candidate are calculated (chimeric spectra, e.g. Sage
        chimera or MSFragger DDA+ searches). Not used for crosslinks. Default is False.
        """
        return self.data.get("chimericPeakSubtraction", False)

    @property
    def perform_rescoring(self) -> bool:
        """
        Get flag whether to run rescoring (percolator / mokapot), plotting and quantification. Default is True.

        If set to False, the input tab files for rescoring are still created, but rescoring is not executed.
        """
        return self.data.get("perform_rescoring", True)

    @property
    def perform_plotting(self) -> bool:
        """Get flag whether to generate the summary plots after rescoring. Default is True."""
        return self.data.get("perform_plotting", True)

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
        """
        Get ptm localization flag from the config file.

        If False, PTMs are rescored as localized by the search engine. If True, all possible localizations of the
        PTMs given in ptmLocalizationOptions are generated and scored.
        """
        return self.data.get("ptm_localization", False)

    @property
    def ptm_model(self) -> bool:
        """Whether a PTM-aware intensity model (model name containing "PTM") is used."""
        return "PTM" in self.models.get("intensity", "")

    @property
    def ptm_fdr_control(self) -> bool:
        """
        Get flag to run percolator separately on unmodified and modified PSMs (separate FDR control).

        Defaults to true if a PTM model is used without PTM localization, false otherwise. Unlike localization,
        this never generates site permutations.
        """
        return self.data.get("ptm_fdr_control", self.ptm_model and not self.ptm_localization)

    @property
    def ptm_localization_options(self) -> dict:
        """Get ptm localization dictionary from the config file."""
        return self.data.get("ptmLocalizationOptions", {})

    @property
    def ptm_sites(self) -> dict[int, list[str]]:
        """
        Get the PTMs to localize and the residues each of them can be located on.

        PTMs are given as a mapping of unimod id to possible residues, e.g. {"21": ["S", "T", "Y"], "7": ["R"]}.
        The legacy format with a single "unimod_id" and "possible_sites" is supported as well.
        Returns an empty mapping if ptm_localization is disabled.
        """
        if not self.ptm_localization:
            return {}
        options = self.ptm_localization_options
        if "ptms" in options:
            return {int(unimod_id): list(sites) for unimod_id, sites in options["ptms"].items()}
        unimod_id = options.get("unimod_id", 0)
        if unimod_id:
            return {int(unimod_id): list(options.get("possible_sites", []))}
        return {}

    @property
    def ptm_decoy_residues(self) -> list[str]:
        """
        Get the decoy residues used to estimate the false localization rate (FLR) of PTM localization.

        Localizations of the PTMs on these residues, which cannot carry the PTMs (e.g. A for phosphorylation), are
        scored as additional (passive) candidates. Spectra in which a decoy candidate would win are used to estimate
        the FLR. The residues are given in order of priority: every peptide uses the first one it contains. Default:
        DEFAULT_DECOY_RESIDUES (A, G, L, V) without the residues the localized PTMs can be located on; an empty list
        disables the decoys. Empty if ptm_localization is disabled.
        """
        if not self.ptm_localization:
            return []
        residues = self.ptm_decoy_priority_residues
        if self.ptm_decoy_selection == "adjacent":
            residues += [aa for aa in ADJACENT_DECOY_RESIDUES if aa not in residues]
        return residues

    @property
    def ptm_decoy_priority_residues(self) -> list[str]:
        """Get the decoy residues of ptmLocalizationOptions.decoy_residues (or the defaults), in order of priority."""
        if not self.ptm_localization:
            return []
        residues = self.ptm_localization_options.get("decoy_residues")
        if residues is None:
            sites = {aa for site_list in self.ptm_sites.values() for aa in site_list}
            return [aa for aa in DEFAULT_DECOY_RESIDUES if aa not in sites]
        return list(residues)

    @property
    def ptm_localization_probability(self) -> str:
        """
        Get how localization probabilities are computed (ptmLocalizationOptions.loc_probability).

        "pep" (default): from the posterior error probabilities of all candidates, (1 / PEP) / sum(1 / PEP). "score":
        from the scores with the slope of the PEP fit of the training step, robust to collapsed PEP estimates.
        """
        return self.ptm_localization_options.get("loc_probability", LOCALIZATION_PROBABILITIES[0])

    @property
    def ptm_decoy_weighting(self) -> str:
        """
        Get how many false localizations a decoy win stands for (ptmLocalizationOptions.decoy_weighting).

        "count" (default): (number of real candidates - 1) / (number of decoy candidates) of the spectrum. "distance":
        matched by the distance of the decoy and wrong real candidates to the called site.
        """
        return self.ptm_localization_options.get("decoy_weighting", DECOY_WEIGHTINGS[0])

    @property
    def ptm_flr_groups(self) -> str:
        """
        Get the groups of PSMs with separate FLR estimates (ptmLocalizationOptions.flr_groups).

        "none" (default): one estimate over all PSMs. "residue": per residue of the called sites (e.g. S, T, Y).
        """
        return self.ptm_localization_options.get("flr_groups", FLR_GROUPS[0])

    @property
    def ptm_flr_estimator(self) -> str:
        """
        Get how the FLR is estimated (ptmLocalizationOptions.flr_estimator).

        "hybrid" (default): the decoy estimate combined with model based estimates per called residue and neighbor
        class, which stay reliable above 1-2% FLR (see hybrid_localization_flr). "decoy": from the decoy candidates.
        """
        return self.ptm_localization_options.get("flr_estimator", FLR_ESTIMATORS[0])

    @property
    def ptm_sdi_features(self) -> bool:
        """
        Get whether localization candidates are scored with site-determining ion features (ptmLocalizationOptions).

        sdi_features: the stored predictions and annotation of the site-determining ions against the best competitor
        of every candidate are used as additional localization features (see localization.sdi_features).
        """
        return self.ptm_localization_options.get("sdi_features", True)

    @property
    def ptm_ambiguity_check(self) -> bool:
        """
        Get whether ambiguous localizations are ranked last in the FLR estimation (ptmLocalizationOptions).

        ambiguity_check: localizations without an observed site-determining ion against the runner-up, or whose
        observed ones favour the runner-up, are ranked after all others of their group by the FLR estimates and
        reported as localization_ambiguous with localization_runner_up (see localization.ambiguity_flags).
        """
        return self.ptm_localization_options.get("ambiguity_check", True)

    @property
    def ptm_decoy_selection(self) -> str:
        """
        Get how decoy sites are selected (ptmLocalizationOptions.decoy_selection).

        "residue" (default): one decoy residue per peptide by priority (decoy_residues). "adjacent": additionally every
        residue in ADJACENT_DECOY_RESIDUES within 2 residues of a possible site, so that wrong localizations on
        neighboring sites have decoys at the same distance (decoy residues include ADJACENT_DECOY_RESIDUES).
        """
        return self.ptm_localization_options.get("decoy_selection", "residue")

    @property
    def ptm_flr_neighbor_split(self) -> bool:
        """
        Get whether the FLR is estimated separately for localizations with and without a possible site next to them.

        Wrong localizations are mostly on a neighboring site, so localizations with a possible site next to the called
        site have a higher FLR at the same localization score (ptmLocalizationOptions.flr_neighbor_split).
        """
        return bool(self.ptm_localization_options.get("flr_neighbor_split", False))

    @property
    def ptm_permutation_sites(self) -> dict[int, list[str]]:
        """
        Get the residues to generate localization candidates on: the possible sites and the decoy residues.

        With decoy_selection "adjacent", the decoy candidates are added after the conversion of the search results
        (add_adjacent_decoy_candidates), so only the possible sites are permuted.
        """
        decoys = self.ptm_decoy_residues if self.ptm_decoy_selection == "residue" else []
        return {unimod_id: sites + [r for r in decoys if r not in sites] for unimod_id, sites in self.ptm_sites.items()}

    @property
    def ptm_decoy_proxies(self) -> dict[str, Union[str, dict[str, float]]]:
        """
        Get the replacements used to predict decoy localization candidates, e.g. {"A[UNIMOD:21]": "S[UNIMOD:21]"}.

        Prediction models do not support PTMs on decoy residues, so these candidates are predicted with the PTM on a
        possible site residue at the same position (ptmLocalizationOptions.decoy_proxy, default: the first possible
        site residue of the PTM, e.g. S for phosphorylation). The fragment ions are annotated with the decoy masses.
        With decoy_proxy "search", every decoy candidate is predicted on one of the residues in DECOY_PROXY_RESIDUES
        (S or T for phosphorylation), chosen with the frequencies of these residues in the localizations of the search
        engine (written by the preprocessing to decoy_proxy_fractions), e.g. {"A[UNIMOD:21]": {"S[UNIMOD:21]": 0.8,
        "T[UNIMOD:21]": 0.2}}.
        """
        proxy_residue = self.ptm_localization_options.get("decoy_proxy")
        fractions = self.ptm_decoy_proxy_fractions if proxy_residue == "search" else {}
        proxies: dict[str, Union[str, dict[str, float]]] = {}
        for unimod_id, sites in self.ptm_sites.items():
            if not sites:
                continue
            proxy: Union[str, dict[str, float]] = proxy_residue if proxy_residue in sites else sites[0]
            proxy = f"{proxy}[UNIMOD:{unimod_id}]"
            if fractions.get(str(unimod_id)):
                proxy = {f"{aa}[UNIMOD:{unimod_id}]": fraction for aa, fraction in fractions[str(unimod_id)].items()}
            for decoy in self.ptm_decoy_residues:
                proxies[f"{decoy}[UNIMOD:{unimod_id}]"] = proxy
        return proxies

    @property
    def decoy_proxy_fractions_file(self) -> Path:
        """Get the file with the frequencies of the proxy residues in the search engine localizations."""
        return self.output / "msms" / "decoy_proxy_fractions.json"

    @property
    def ptm_decoy_proxy_fractions(self) -> dict[str, dict[str, float]]:
        """Get the frequencies of the proxy residues per unimod id, empty if they were not determined."""
        if not self.decoy_proxy_fractions_file.is_file():
            return {}
        return json.loads(self.decoy_proxy_fractions_file.read_text())

    @property
    def ptm_unimod_id(self) -> int:
        """Get the unimod id of the PTM to localize if exactly one PTM is localized, 0 otherwise."""
        ptm_sites = self.ptm_sites
        return next(iter(ptm_sites)) if len(ptm_sites) == 1 else 0

    @property
    def ptm_possible_sites(self) -> list:
        """Get all residues any of the localized PTMs can be located on."""
        return sorted({site for sites in self.ptm_sites.values() for site in sites})

    @property
    def ptm_use_neutral_loss(self) -> bool:
        """Get neutral loss flag to indicate whether to add a score for this or not (only used for localization)."""
        return self.ptm_localization and self.ptm_localization_options.get("neutral_loss", False)

    @property
    def ptm_phospho_loss_residues(self) -> str:
        """
        Get the phosphorylated residues annotated with the neutral loss of phosphoric acid (H3PO4).

        S and T lose phosphoric acid. With decoy_neutral_loss, phosphorylations on the decoy residues are annotated
        like S and T as well, so that decoy candidates have the same neutral loss features as the S/T candidates
        they stand in for.
        """
        residues = "ST"
        if self.ptm_localization_options.get("decoy_neutral_loss", False):
            residues += "".join(self.ptm_decoy_residues)
        return residues

    @property
    def ptm_plain_loss_residues(self) -> Optional[str]:
        """
        Get the residues for counting the neutral losses of plain decoys, None if the decoys are already plain.

        With decoy_neutral_loss, the neutral losses are counted a second time without the decoy residues, so that every
        decoy candidate can also be scored as a plain decoy (used for the PSMs localized on Y, which do not lose
        phosphoric acid).
        """
        if self.ptm_phospho_loss_residues == "ST":
            return None
        return "ST"

    ########################
    # functions start here #
    ########################

    def check(self):
        """Validate the configuration."""
        self._check_tmt()

        if self.job_type == "SpectralLibraryGeneration":
            self._check_for_speclib()

        # Skip alphapept and Koina checks if using DLOmix
        if self.models.get("dlomix_intensity"):
            self.check_dlomix()
        else:
            if "alphapept" in self.models["intensity"].lower():
                self._check_for_alphapept()

            self._check_koina_model_availability()

        if self.quantification:
            self._check_quantification_options()
            if self.quant_precursor == "search_engine":
                self._check_quantification()
            self._check_fasta()

        self._check_ptm_localization()
        self._check_percolator_subset()
        self._check_chimeric_peak_subtraction()

    def _check_chimeric_peak_subtraction(self):
        """Check that chimericPeakSubtraction is a boolean."""
        if not isinstance(self.chimeric_peak_subtraction, bool):
            raise ValueError(f"chimericPeakSubtraction must be true or false, got {self.chimeric_peak_subtraction!r}.")

    def _check_percolator_subset(self):
        """Check that percolatorSubsetMaxTrain is a positive number of PSMs if given."""
        subset = self.percolator_subset_max_train
        if subset is not None and (isinstance(subset, bool) or not isinstance(subset, int) or subset <= 0):
            raise ValueError(f"percolatorSubsetMaxTrain must be a positive number of PSMs, got {subset!r}.")

    def _check_ptm_localization(self):
        """Check that the PTMs and residues to localize are specified if PTM localization is enabled."""
        options = self.ptm_localization_options
        if not self.ptm_localization:
            if "ptms" in options or options.get("unimod_id"):
                logger.warning(
                    "ptmLocalizationOptions are specified but ptm_localization is disabled. PTMs are rescored as "
                    "localized by the search engine. Set ptm_localization to true to localize them."
                )
            return
        if "ptms" in options and not isinstance(options["ptms"], dict):
            raise ValueError(
                'ptmLocalizationOptions.ptms must map unimod ids to possible residues, e.g. {"21": ["S", "T", "Y"]}.'
            )
        ptm_sites = self.ptm_sites
        if not ptm_sites:
            raise ValueError(
                "ptm_localization is enabled but no PTM to localize is specified. Specify the PTMs and the residues "
                'they can be located on in ptmLocalizationOptions, e.g. "ptms": {"21": ["S", "T", "Y"]}.'
            )
        valid_residues = set("ACDEFGHIKLMNPQRSTVWY")
        for unimod_id, sites in ptm_sites.items():
            if not sites:
                raise ValueError(
                    f"No possible residues specified for PTM UNIMOD:{unimod_id} in ptmLocalizationOptions."
                )
            invalid_sites = [site for site in sites if site not in valid_residues]
            if invalid_sites:
                raise ValueError(
                    f"Invalid residues {invalid_sites} for PTM UNIMOD:{unimod_id} in ptmLocalizationOptions. "
                    "Residues must be given as one-letter amino acid codes."
                )
        self._check_ptm_decoys(ptm_sites, valid_residues)

    def _check_ptm_decoys(self, ptm_sites: dict[int, list[str]], valid_residues: set[str]):
        """Check the decoy residues and the options of the false localization rate estimation."""
        decoys = self.ptm_decoy_residues
        invalid_decoys = [r for r in decoys if r not in valid_residues]
        if invalid_decoys:
            raise ValueError(f"Invalid decoy residues {invalid_decoys} in ptmLocalizationOptions.")
        overlapping = [r for r in decoys if any(r in sites for sites in ptm_sites.values())]
        if overlapping:
            raise ValueError(
                f"Decoy residues {overlapping} can also carry a localized PTM. Decoy residues must be residues the "
                "PTMs cannot be located on, e.g. A for phosphorylation."
            )
        for option, value, allowed in (
            ("loc_probability", self.ptm_localization_probability, LOCALIZATION_PROBABILITIES),
            ("decoy_weighting", self.ptm_decoy_weighting, DECOY_WEIGHTINGS),
            ("flr_groups", self.ptm_flr_groups, FLR_GROUPS),
            ("flr_estimator", self.ptm_flr_estimator, FLR_ESTIMATORS),
        ):
            if value not in allowed:
                raise ValueError(f"Invalid ptmLocalizationOptions.{option} {value!r}, must be one of {allowed}.")
        # read the raw options: the bool-annotated properties would trip typeguard before this check
        for option in ("sdi_features", "ambiguity_check"):
            value = self.ptm_localization_options.get(option, True)
            if not isinstance(value, bool):
                raise ValueError(f"Invalid ptmLocalizationOptions.{option} {value!r}, must be true or false.")

    def _check_koina_model_availability(self):
        """Check if Koina model is available."""
        # This will give error automaticly in Koina if model is not available on the server.
        # Koina has function called "_is_model_ready" that checks if model is available
        _ = Koina(model_name=self.models["intensity"])

    def check_dlomix(self):
        """Load the DLOmix pipeline once in the parent process (from ``check``, before the worker pools).

        Validates the path (fail fast on a bad model) and primes the module-level cache so the
        fork-based Pool workers created afterwards reuse this single load instead of each reloading.
        NOTE: loading TF in the parent and then forking can deadlock when a worker runs inference;
        _load_dlomix_pipeline pins TF to single-threaded execution to mitigate that.
        """
        if self.dlomix_pipeline is None:
            raise ValueError(
                f"Could not load DLOmix pipeline from {self.models.get('dlomix_intensity')}. "
                "Please check that 'dlomix_intensity' in your config points to a valid saved pipeline."
            )

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

    def _check_quantification_options(self):
        """Check the quantificationOptions and that MS1 quantification (flashlfq) can run with this config."""
        for option, value, allowed in (
            ("precursor", self.quant_precursor, QUANT_PRECURSOR_METHODS),
            ("protein", self.quant_protein, QUANT_PROTEIN_METHODS),
            ("normalization", self.quant_normalization, QUANT_NORMALIZATIONS),
        ):
            if value not in allowed:
                raise ValueError(f"Invalid quantificationOptions.{option} {value!r}, must be one of {allowed}.")
        fdrs = (("psmFdr", self.quant_psm_fdr), ("peptideFdr", self.quant_peptide_fdr),
                ("proteinFdr", self.quant_protein_fdr))  # fmt: skip
        for option, value in fdrs:
            if not _is_number(value) or not 0 < value <= 1:
                raise ValueError(f"quantificationOptions.{option} must be a q-value in (0, 1], got {value!r}.")
        ppm = self.quant_ppm_tolerance
        if not _is_number(ppm) or ppm <= 0:
            raise ValueError(f"quantificationOptions.ppmTolerance must be a positive number, got {ppm!r}.")
        ppm_filter = self.quant_ppm_filter
        if ppm_filter is not None and (not _is_number(ppm_filter) or ppm_filter <= 0):
            raise ValueError(f"quantificationOptions.ppmFilter must be a positive number or null, got {ppm_filter!r}.")
        if self.quant_precursor == "flashlfq":
            self._check_flashlfq()

    def _check_flashlfq(self):
        """Check that MS1 quantification with FlashLFQ can run with this config."""
        if self.fdr_estimation_method != "percolator":
            raise ValueError("Quantification with flashlfq requires fdr_estimation_method percolator.")
        if self.ptm_localization:
            raise ValueError("Quantification with flashlfq is not supported together with ptm_localization.")
        if self.spectra_type not in ("raw", "mzml"):
            raise ValueError(
                f"Quantification with flashlfq needs MS1 spectra from raw or mzml files, got spectra_type "
                f"{self.spectra_type}."
            )
        if not self.flashlfq_exe.is_file():
            raise FileNotFoundError(
                f"FlashLFQ not found at {self.flashlfq_exe}. Download the command line version (FlashLFQ.zip) from "
                "https://github.com/smith-chem-wisc/FlashLFQ/releases and set quantificationOptions.flashlfqExe to "
                "its CMD.dll."
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
