import logging
from pathlib import Path
from sys import platform
from typing import Any, List, Optional, Union

import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.annotation.annotation import annotate_spectra
from spectrum_fundamentals.fragments import compute_peptide_mass
from spectrum_fundamentals.mod_string import internal_without_mods, maxquant_to_internal
from spectrum_io.file import csv
from spectrum_io.raw import ThermoRaw
from spectrum_io.search_result import Mascot, MaxQuant, MSFragger, Sage
from spectrum_io.spectral_library import digest as digest_io

from ..data.spectra import FragmentType, Spectra

logger = logging.getLogger(__name__)


# SpectralLibrary


def gen_lib(input_file: Union[str, Path]) -> Spectra:
    """
    Generate a spectral library from a given input.

    This function reads an input file that follows the specifications provided in the usage section for
    :doc:`../../peptides_format`, creates a `Spectra` object from it and returns it.

    :param input_file: A csv file containing modified sequences, CE, precursor charge and fragmentation method.
    :return: Spectra object from the read input.
    """
    library_df = csv.read_file(input_file)
    library_df.columns = library_df.columns.str.upper()
    library = Spectra()
    library.add_columns(library_df)
    return library


def digest(
    fasta: Union[str, Path],
    output: Union[str, Path],
    fragmentation: str,
    digestion: str,
    cleavages: int,
    db: str,
    enzyme: str,
    special_aas: str,
    min_length: int,
    max_length: int,
):
    """
    Perform an in-silico digest of a given fasta file.

    This function takes sequences from a fasta file and performs digestion using a given
    protease with specific settings before writing the resulting peptides to the given output file.

    :param fasta: Path to fasta file containing sequences to digest
    :param output: Path to the output file containing the peptides
    :param fragmentation: The fragmentation method to use, can be HCD / CID
    :param digestion: TODO
    :param cleavages: The number of allowed miscleaveages
    :param db: The desired database to produce, can be target, decoy, or both
    :param enzyme: The protease to use for digestion TODO list available proteases
    :param special_aas: Special aminoacids used for decoy generation, TODO when relevant
    :param min_length: Minimal length of digested peptides
    :param max_length: Maximal length of digested peptides
    """
    if isinstance(output, str):
        output = Path(output)
    cmd = [
        "--fasta",
        f"{fasta}",
        "--prosit_input",
        f"{output / 'prosit_input.csv'}",
        "--fragmentation",
        f"{fragmentation}",
        "--digestion",
        f"{digestion}",
        "--cleavages",
        f"{cleavages}",
        "--db",
        f"{db}",
        "--enzyme",
        f"{enzyme}",
        "--special-aas",
        f"{special_aas}",
        "--min-length",
        f"{min_length}",
        "--max-length",
        f"{max_length}",
    ]
    digest_io.main(cmd)


def filter_peptides_for_model(peptides: pd.DataFrame, model: str) -> pd.DataFrame:
    """
    Filter search results to support a given peptide prediction model.

    Depending on the model used for peptide property prediction, PSMs in search results need to be filtered out.
    This function provides a shortcut using the model name for filtering instead of supplying all possible filters manually.

    :param peptides: Dataframe containing search results to be filtered
    :param model: Name of a peptide property prediction model to use as filter

    :raises ValueError: if an unsupported model is supplied

    :return: The filtered dataframe to be used with the given model.
    """
    if "prosit" in model.lower():
        filter_kwargs = {
            "min_length": 7,
            "max_length": 30,
            "max_charge": 6,
        }
    else:
        raise ValueError(f"The model {model} is not known.")

    return filter_peptides(peptides, **filter_kwargs)


def filter_peptides(peptides: pd.DataFrame, min_length: int, max_length: int, max_charge: int) -> pd.DataFrame:
    """
    Filter search results using given constraints.

    This function filters provided search results by peptide length, precursor charge,
    unsupported special aminoacids, and unsupported modifications.

    :param peptides: Dataframe containing search results to be filtered
    :param min_length: The minimal length of a peptide to be retained
    :param max_length: The maximal length of a peptide to be retained
    :param max_charge: The maximal precursor charge of a peptide to be retained

    :return: The filtered dataframe given the provided constraints.
    """
    return peptides[
        (peptides["PEPTIDE_LENGTH"] <= max_length)
        & (peptides["PEPTIDE_LENGTH"] >= min_length)
        & (peptides["PRECURSOR_CHARGE"] <= max_charge)
        & (~peptides["MODIFIED_SEQUENCE"].str.contains(r"\(ac\)"))
        & (~peptides["MODIFIED_SEQUENCE"].str.contains(r"\(Acetyl \(Protein N-term\)\)"))
        & (~peptides["MODIFIED_SEQUENCE"].str.contains(r"\[UNIMOD\:21\]"))
        & (~peptides["SEQUENCE"].str.contains("U"))
    ]


def process_and_filter_spectra_data(library: Spectra, model: str, tmt_label: Optional[str] = None) -> Spectra:
    """
    Process and filter the spectra data in the given SpectralLibrary object.

    This function applies various modifications and filters to the 'spectra_data' DataFrame
    in the provided SpectralLibrary object. It modifies the 'MODIFIED_SEQUENCE' column,
    converts the 'MODIFIED_SEQUENCE' to internal format, extracts 'SEQUENCE', and filters
    out certain entries based on specific criteria. The specification of the internal file format can be found at
    :doc:`../../internal_format`.

    :param library: A Spectra object containing the raw peptides to be proccessed and filtered.
    :param model: The peptide property prediction model to filter the spectra for
    :param tmt_label: Optional tmt-label to consider when processing peptides. If given, the corresponding
        fixed modification for the N-terminus and lysin will be added
    :return: The processed and filtered Spectra object
    """
    # add fixed mods and translate to internal format
    library.spectra_data["MODIFIED_SEQUENCE"] = library.spectra_data["MODIFIED_SEQUENCE"].apply(lambda x: "_" + x + "_")

    fixed_mods = {"C": "C[UNIMOD:4]"}
    if tmt_label is not None and tmt_label != "":
        unimod_tag = c.TMT_MODS[tmt_label]
        fixed_mods = {"C": "C[UNIMOD:4]", "^_": f"_{unimod_tag}-", "K": f"K{unimod_tag}"}

    library.spectra_data["MODIFIED_SEQUENCE"] = maxquant_to_internal(
        library.spectra_data["MODIFIED_SEQUENCE"], fixed_mods=fixed_mods
    )

    # get sequence and its length
    library.spectra_data["SEQUENCE"] = internal_without_mods(library.spectra_data["MODIFIED_SEQUENCE"])
    library.spectra_data["PEPTIDE_LENGTH"] = library.spectra_data["SEQUENCE"].apply(lambda x: len(x))

    # filter
    logger.info(f"No of sequences before Filtering is {len(library.spectra_data)}")
    library.spectra_data = filter_peptides_for_model(library.spectra_data, model)
    logger.info(f"No of sequences after Filtering is {len(library.spectra_data)}")

    library.spectra_data["MASS"] = library.spectra_data["MODIFIED_SEQUENCE"].apply(lambda x: compute_peptide_mass(x))

    return library


# CeCalibration
def load_search(input_file: Union[str, Path]) -> pd.DataFrame:
    """
    Load search results.

    Given a path to a file containing search results in Oktoberfest format, the function reads the results and returns them.
    The specification of the internal file format can be found at :doc:`../../internal_format`.

    :param input_file: Path to the file containing search results in the internal Oktoberfest format.
    :return: dataframe containing the search results.
    """
    return csv.read_file(input_file)


def convert_search(
    input_path: Union[str, Path], output_file: Union[str, Path], search_engine: str, tmt_label: str = ""
):
    """
    Convert search results to Oktoberfest format.

    Given a path to a file or directory containing search results from supported search engines,
    the function parses, converts them to the internal format used by Oktoberfest and writes them
    to a specified location. The specification of the internal file format can be found at :doc:`../../internal_format`.

    :param input_path: Path to the directory or file containing the search results.
    :param output_file: Path to the location where the converted search results should be written to.
    :param search_engine: The search engine used to produce the search results,
        currently supported are "Maxquant", "Mascot" and "MSFragger"
    :param tmt_label: Optional tmt-label to consider when processing peptides. If given, the corresponding
        fixed modification for the N-terminus and lysin will be added
    :raises ValueError: if an unsupported search engine was given
    """
    search_engine = search_engine.lower()
    search_result: Any
    if search_engine == "maxquant":
        search_result = MaxQuant
    elif search_engine == "msfragger":
        search_result = MSFragger
    elif search_engine == "mascot":
        search_result = Mascot
    elif search_engine == "sage":
        search_result = Sage
    else:
        raise ValueError(f"Unknown search engine provided: {search_engine}")

    search_result(input_path).generate_internal(tmt_labeled=tmt_label, out_path=output_file)


def list_spectra(input_dir: Union[str, Path], file_format: str) -> List[Path]:
    """
    Return a list of all spectra files of a given format.

    Given an input directory, the function searches all files containing spectra and returns a list of paths pointing to the files.
    Files are included if the extension matches the provided format (case-insensitive).
    In case the input directory is a file, the function will check if it matches the format and return it wrapped in a list.

    :param input_dir: Path to the directory to scan for spectra files
    :param file_format: Format of spectra files that match the file extension (case-insensitive), can be "mzML", "RAW" or "pkl".
    :raises NotADirectoryError: if the specified input directory does not exist
    :raises ValueError: if the specified file format is not supported
    :raises AssertionError: if no files in the provided input directory match the provided file format
    :return: A list of paths to all spectra files found in the given directory
    """
    if isinstance(input_dir, str):
        input_dir = Path(input_dir)
    raw_files = []

    if not file_format.lower() in ["mzml", "raw", "pkl"]:
        raise ValueError(f"File format {file_format} unknown. Must be one of mzML, RAW or pkl.")

    if input_dir.is_file() and input_dir.suffix.lower().endswith(file_format.lower()):
        raw_files.append(input_dir)
    elif input_dir.is_dir():
        glob_pattern = _get_glob_pattern(file_format)
        raw_files = list(input_dir.glob(glob_pattern))
    else:
        raise NotADirectoryError(f"{input_dir} does not exist.")

    if not raw_files:
        raise AssertionError(
            f"There are no spectra files with the extension {file_format.lower()} in the provided input_dir {input_dir}. "
            "Please check."
        )

    return raw_files


def _get_glob_pattern(spectra_type: str) -> str:
    """
    Get global pattern depending on spectra_type.

    :param spectra_type: The type of spectra file. Accepted values are "raw" or "mzml".
    :return: The glob pattern corresponding to the spectra_type.
    :raises ValueError: If an unsupported spectra_type is provided.
    """
    if spectra_type.lower() == "raw":
        return "*.[rR][aA][wW]"
    elif spectra_type.lower() == "mzml":
        return "*.[mM][zZ][mM][lL]"
    else:
        raise ValueError(f"{spectra_type} is not supported as rawfile-type")


def split_search(
    search_results: pd.DataFrame,
    output_dir: Union[str, Path],
    filenames: Optional[List[str]] = None,
) -> List[str]:
    """
    Split search results by spectrum file.

    Given a list of spectrum file names from which search results originate the provided search results are split
    and filename specific csv files are written to the provided output directory. The provided file names need to
    correspond to the spectrum file identifier in the "RAW_FILE" column of the provided search results. The search
    results need to be provided in internal format (see :doc:`../../internal_format`).
    If the list of file names is not provided, all spectrum file identifiers are considered, otherwise only the
    identifiers found in the list are taken into account for writing the individual csv files.
    The output file names follow the convention <filename>.rescore.
    If a file name is not found in the search results, it is ignored and a warning is printed.
    The function returns a list of file names for which search results are available, removing the ones that were
    ignored if a list of file names was provided.

    :param search_results: search results in internal format
    :param output_dir: directory in which to store individual csv files containing the search results for
        individual filenames
    :param filenames: optional list of spectrum filenames that should be considered. If not provided, all spectrum file
        identifiers in the search results are considered.

    :return: list of file names for which search results could be found
    """
    if isinstance(output_dir, str):
        output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    if filenames is None:
        filenames = search_results["RAW_FILE"].unique()

    grouped_search_results = search_results.groupby("RAW_FILE")

    filenames_found = []
    for filename in filenames:
        output_file = (output_dir / filename).with_suffix(".rescore")
        logger.info(f"Creating split msms.txt file {output_file}")
        try:
            grouped_search_results.get_group(filename).to_csv(output_file)
            filenames_found.append(filename)
        except KeyError:
            logger.warning(
                f"The search results do not contain search results for the provided file name {filename}. "
                "If this is not intended, please verify that the file names are written correctly in the "
                f"search results. {filename} is ignored."
            )
    return filenames_found


def merge_spectra_and_peptides(spectra: pd.DataFrame, search: pd.DataFrame) -> Spectra:
    """
    Merge peptides with spectra.

    This function takes spectra and and peptides from a search and merges them based on the
    RAW File identifier and the scan number using the "RAW_FILE" and "SCAN_NUMBER" column of
    the two dataframes provided.

    :param spectra: MS2 spectra of a mass spectrometry run
    :param search: Peptides from a search for the given spectra in internal format (see
        :doc:`../../internal_format`)

    :return: Dataframe containing the matched pairs of peptides and spectra (PSMs)
    """
    logger.info("Merging rawfile and search result")
    psms = search.merge(spectra, on=["RAW_FILE", "SCAN_NUMBER"])
    logger.info(f"There are {len(psms)} matched identifications")

    library = Spectra()
    library.add_columns(psms)

    return library


def annotate_spectral_library(psms: Spectra, mass_tol: Optional[float] = None, unit_mass_tol: Optional[str] = None):
    """
    Annotate spectral library with peaks and mass.

    This function annotates a given spectral library with peak intensities and mass to charge ratio,
    as well as the calculated monoisotopic mass of the precursor ion.
    The additional information is added to the provided spectral library.

    :param psms: Spectral library to be annotated.
    :param mass_tol: The mass tolerance allowed for retaining peaks
    :param unit_mass_tol: The unit in which the mass tolerance is given
    """
    logger.info("Annotating spectra...")
    df_annotated_spectra = annotate_spectra(psms.spectra_data, mass_tol, unit_mass_tol)
    logger.info("Finished annotating.")
    psms.spectra_data.drop(columns=["INTENSITIES", "MZ"], inplace=True)  # TODO check if this is needed
    psms.add_matrix(df_annotated_spectra["INTENSITIES"], FragmentType.RAW)
    psms.add_matrix(df_annotated_spectra["MZ"], FragmentType.MZ)
    psms.add_column(df_annotated_spectra["CALCULATED_MASS"].to_numpy(), "CALCULATED_MASS")


def load_spectra(filename: Union[str, Path], parser: str = "pyteomics") -> pd.DataFrame:
    """
    Read spectra from a given file.

    This function reads MS2 spectra from a given mzML or pkl file using a specified parser. The file ending
    is used to determine the correct parsing method.

    :param filename: Path to mzML / pkl file containing MS2 spectra to be loaded.
    :param parser: Name of the package to use for parsing the mzml file, can be "pyteomics" or "pymzml".
        Only used for parsing of mzML files.
    :raises ValueError: if the filename does not end in either ".pkl" or ".mzML" (case-insensitive)
    :return: measured spectra with metadata.
    """
    if isinstance(filename, str):
        filename = Path(filename)

    format_ = filename.suffix.lower()
    if format_ == ".mzml":
        return ThermoRaw.read_mzml(
            source=filename, package=parser, search_type=""
        )  # TODO in spectrum_io, remove unnecessary argument
    elif format_ == ".pkl":
        results = pd.read_pickle(filename)
        # TODO in spectrum-io in case median_RETENTION_TIME is still passed
        # TODO also change name for ion mobility
        results.rename(columns={"median_RETENTION_TIME": "RETENTION_TIME"}, inplace=True)
        return results

    else:
        raise ValueError(f"Unrecognized file format: {format_}. Only .mzML and .pkl files are supported.")


def convert_spectra_to_mzml(
    raw_file: Union[str, Path], output_file: Union[str, Path], thermo_exe: Optional[Union[str, Path]] = None
):
    """
    Convert spectra to mzML format.

    This function converts RAW files containing spectra from a mass spectrometry run into mzML format. The process
    makes use of ThermoRawFileParser and requires mono being installed if it is used on MacOS or linux.

    :param raw_file: Path to RAW file to be converted to mzML format
    :param output_file: Path to the location where the converted spectra should be written to.
    :param thermo_exe: Path to the executable of ThermoRawFileParser, default depends on the OS:
        for Windows, it is "ThermoRawFileParser.exe", while for Linux and MacOS, it is
        "/opt/compomics/ThermoRawFileParser.exe".
    """
    if thermo_exe is None:
        if "linux" in platform or platform == "darwin":
            thermo_exe = "/opt/compomics/ThermoRawFileParser.exe"
        else:
            thermo_exe = "ThermoRawFileParser.exe"

    raw = ThermoRaw()
    raw.convert_raw_mzml(input_path=raw_file, output_path=output_file, thermo_exe=thermo_exe)
