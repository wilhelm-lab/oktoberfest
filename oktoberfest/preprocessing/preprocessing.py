import logging
from itertools import chain, product, repeat
from pathlib import Path
from sys import platform
from typing import Any, Dict, List, Optional, Union

import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.annotation.annotation import annotate_spectra
from spectrum_fundamentals.fragments import compute_peptide_mass
from spectrum_fundamentals.mod_string import internal_without_mods, maxquant_to_internal
from spectrum_io.d import convert_d_hdf, read_and_aggregate_timstof
from spectrum_io.file import csv
from spectrum_io.raw import ThermoRaw
from spectrum_io.search_result import Mascot, MaxQuant, MSFragger, Sage
from spectrum_io.spectral_library.digest import get_peptide_to_protein_map

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


def generate_metadata(
    peptides: List[str],
    collision_energy: Union[int, List[int]],
    precursor_charge: Union[int, List[int]],
    fragmentation: Union[str, List[str]],
    proteins: Optional[List[List[str]]] = None,
) -> pd.DataFrame:
    """
    Create metadata about peptides for a spectral library.

    This function generates a pandas DataFrame containing metadata for peptides in a spectral library.
    Each row in the DataFrame represents a unique combination of a peptide, collision energy,
    precursor charge, and fragmentation. If multiple collision energies, precursor charges, or fragmentations
    are provided, the function creates all possible combinations for each peptide.
    An optional protein list can be provided which has to have the same length as the number of peptides.

    :param peptides: A list of peptides for which metadata is generated.
    :param collision_energy: A list of collision energies corresponding to each peptide.
    :param precursor_charge: A list of precursor charges corresponding to each peptide.
    :param fragmentation: A list of fragmentation methods corresponding to each peptide.
    :param proteins: An optional list of proteins associated with each peptide.
        If provided, it must have the same length as the number of peptides.
    :raises AssertionError: If the lengths of peptides and proteins is not the same.
    :return: A DataFrame containing metadata with the columns "modified_peptide","collision_energy",
        "precursor_charge","fragmentation" and an optional "proteins" column.
    """
    if isinstance(collision_energy, int):
        collision_energy = [collision_energy]
    if isinstance(precursor_charge, int):
        precursor_charge = [precursor_charge]
    if isinstance(fragmentation, str):
        fragmentation = [fragmentation]

    if proteins is not None and len(proteins) != len(peptides):
        raise AssertionError("Number of proteins must match the number of peptides.")

    combinations = product(peptides, collision_energy, precursor_charge, fragmentation)
    metadata = pd.DataFrame(
        combinations, columns=["modified_sequence", "collision_energy", "precursor_charge", "fragmentation"]
    )

    if proteins is not None:
        n_repeats = len(metadata) // len(proteins)
        metadata["proteins"] = list(
            chain.from_iterable([repeat(";".join(prot_list), n_repeats) for prot_list in proteins])
        )

    return metadata


def digest(
    fasta: Union[str, Path],
    digestion: str,
    missed_cleavages: int,
    db: str,
    enzyme: str,
    special_aas: str,
    min_length: int,
    max_length: int,
) -> Dict[str, List[str]]:
    """
    Digest a given fasta file with specific settings.

    This function performs an in-silico digestion of a fasta file based on the provided settings.
    It returns a dictionary that maps peptides to the list of associated protein IDs.


    :param fasta: Path to fasta file containing sequences to digest
    :param digestion: The type of digestion, one of "full, "semi", "none"
    :param missed_cleavages: The number of allowed miscleaveages
    :param db: The desired database to produce, can be target, decoy, or both
    :param enzyme: The protease to use for digestion TODO list available proteases
    :param special_aas: List of aas to be swapped with preceding aa in reverse sequences.
        This mimics the behaviour of MaxQuant when creating decoys.
    :param min_length: Minimal length of digested peptides
    :param max_length: Maximal length of digested peptides

    :return: A Dictionary that maps peptides (keys) to a list of protein IDs (values).
    """
    return get_peptide_to_protein_map(
        fasta_file=fasta,
        db=db,
        min_len=min_length,
        max_len=max_length,
        enzyme=enzyme,
        digestion=digestion,
        miscleavages=missed_cleavages,
        methionine_cleavage=True,
        use_hash_key=False,
        special_aas=special_aas,
    )


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
        & (~peptides["SEQUENCE"].str.contains("U|X"))
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
    logger.info(f"No of sequences before filtering is {len(library.spectra_data)}")
    library.spectra_data = filter_peptides_for_model(library.spectra_data, model)
    logger.info(f"No of sequences after filtering is {len(library.spectra_data)}")

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
    input_path: Union[str, Path],
    search_engine: str,
    tmt_label: str = "",
    output_file: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Convert search results to Oktoberfest format.

    Given a path to a file or directory containing search results from supported search engines,
    the function parses, converts them to the internal format used by Oktoberfest and returns it as a dataframe.
    If a path to an output file is provided, the converted results are also stored to the specified location.
    The specification of the internal file format can be found at :doc:`../../internal_format`.

    :param input_path: Path to the directory or file containing the search results.
    :param search_engine: The search engine used to produce the search results,
        currently supported are "Maxquant", "Mascot" and "MSFragger"
    :param tmt_label: Optional tmt-label to consider when processing peptides. If given, the corresponding
        fixed modification for the N-terminus and lysin will be added
    :param output_file: Optional path to the location where the converted search results should be written to.
        If this is omitted, the results are not stored.

    :raises ValueError: if an unsupported search engine was given
    :return: A dataframe containing the converted results.
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

    return search_result(input_path).generate_internal(tmt_labeled=tmt_label, out_path=output_file)


def convert_timstof_metadata(
    input_path: Union[str, Path], search_engine: str, output_file: Optional[Union[str, Path]] = None
) -> pd.DataFrame:
    """
    Convert metadata for timstof to Oktoberfest format.

    Given a path to a directory containing search results from supported search engines,
    the function parses, converts metadata relevant for timstof to the internal format used by Oktoberfest and
    returns it as a dataframe.
    If a path to an output file is provided, the converted results are also stored to the specified location.

    :param input_path: Path to the directory or file containing the metadata.
    :param search_engine: The search engine used to produce the search results,
        currently supported is "Maxquant"
    :param output_file: Optional path to the location where the converted metadata should be written to.
        If this is omitted, the metadata are not stored.
    :raises ValueError: if an unsupported search engine was given

    :return: dataframe containing metadata that maps scan numbers to precursors
    """
    search_engine = search_engine.lower()
    if search_engine == "maxquant":
        metadata_df = MaxQuant(input_path).generate_internal_timstof_metadata()
    else:
        raise ValueError(f"Unsupported search engine provided for reading timstof metadata: {search_engine}")

    return metadata_df


def list_spectra(input_dir: Union[str, Path], input_format: str) -> List[Path]:
    """
    Return a list of all spectra files of a given format.

    Given an input directory, the function searches all files containing spectra and returns a list of paths pointing to the files.
    Files are included if the extension matches the provided format (case-insensitive).
    In case the input directory is a file, the function will check if it matches the format and return it wrapped in a list.
    If the format is "d" and the input directory ends with ".d", the function will return the input directory wrapped in a list.

    :param input_dir: Path to the directory to scan for spectra files
    :param input_format: Format of the input for the provided directory. This must match the file extension (mzml, raw, hdf) or
        directory extension (d). Matching is case-insensitive.
    :raises NotADirectoryError: if the specified input directory does not exist
    :raises ValueError: if the specified file format is not supported
    :raises AssertionError: if the provided input directory (d) does not match the provided format or if none of the
        files within the provided input directory (mzml, raw, hdf) match the provided format
    :return: A list of paths to all spectra files found in the given directory
    """
    if isinstance(input_dir, str):
        input_dir = Path(input_dir)
    raw_files = []

    input_format = input_format.lower()

    if input_format not in ["mzml", "raw", "hdf", "d"]:
        raise ValueError(f"Input format {input_format} unknown. Must be one of mzml, raw, d, hdf.")

    if input_dir.is_file() and input_dir.suffix.lower().endswith(input_format):
        raw_files.append(input_dir)
    elif input_dir.is_dir():
        if input_dir.suffix == ".d":
            if input_format == "d":
                raw_files = [input_dir]
            else:
                raise AssertionError(
                    f"Provided a '.d' input directory but the provided input format is {input_format}. Please check."
                )
        else:
            glob_pattern = _get_glob_pattern(input_format)
            raw_files = list(input_dir.glob(glob_pattern))
    else:
        raise NotADirectoryError(f"{input_dir} does not exist.")

    if not raw_files:
        raise AssertionError(
            f"There are no files / directories with the extension {input_format} in the provided input directory {input_dir}. "
            "Please check."
        )

    if len(raw_files) > 1:
        input_type_str = "directories" if input_format == "d" else "files"
    else:
        input_type_str = "directory" if input_format == "d" else "file"
    logger.info(f"Found {len(raw_files)} {input_format} {input_type_str} in the spectra input directory.")

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
    elif spectra_type.lower() == "hdf":
        return "*.[hH][dD][fF]"
    elif spectra_type.lower() == "d":
        return "*.[dD]"
    else:
        raise ValueError(f"{spectra_type} is not supported as rawfile-type")


def split_search(
    search_results: pd.DataFrame, output_dir: Union[str, Path], filenames: Optional[List[str]] = None
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
        logger.info(f"Creating split search results file {output_file}")
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


def split_timstof_metadata(
    timstof_metadata: pd.DataFrame, output_dir: Union[str, Path], filenames: Optional[List[str]] = None
) -> List[str]:
    """
    Split timstof metadata by spectrum file.

    Given a list of spectrum file names from which timstof metadata originate the provided timstof metadata are split
    and filename specific csv files are written to the provided output directory. The provided file names need to
    correspond to the spectrum file identifier in the "RAW_FILE" column of the provided timstof_metadata. The timstof
    metadata need to be provided in internal format  #TODO provided documentation.
    If the list of file names is not provided, all spectrum file identifiers are considered, otherwise only the
    identifiers found in the list are taken into account for writing the individual csv files.
    The output file names follow the convention <filename>.timsmeta.
    If a file name is not found in the timstof metadata, it is ignored and a warning is printed.
    The function returns a list of file names for which timstof metadata are available, removing the ones that were
    ignored if a list of file names was provided.

    :param timstof_metadata: timstof metadata in internal format
    :param output_dir: directory in which to store individual csv files containing the timstof metadata for
        individual filenames
    :param filenames: optional list of spectrum filenames that should be considered. If not provided, all spectrum file
        identifiers in the timstof metadata are considered.

    :return: list of file names for which timstof metadata could be found
    """
    if isinstance(output_dir, str):
        output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    if filenames is None:
        filenames = timstof_metadata["RAW_FILE"].unique()

    grouped_timstof_metadata = timstof_metadata.groupby("RAW_FILE")

    filenames_found = []
    for filename in filenames:
        output_file = (output_dir / filename).with_suffix(".timsmeta")
        logger.info(f"Creating split timstof metadata file {output_file}")
        try:
            grouped_timstof_metadata.get_group(filename).to_csv(output_file)
            filenames_found.append(filename)
        except KeyError:
            logger.warning(
                f"No timstof metadata could be found for the provided file name {filename}. "
                "If this is not intended, please verify that the file names are written correctly in the "
                f"timstof. {filename} is ignored."
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


def load_spectra(
    filename: Union[str, Path], parser: str = "pyteomics", tims_meta_file: Optional[Union[str, Path]] = None
) -> pd.DataFrame:
    """
    Read spectra from a given file.

    This function reads MS2 spectra from a given mzML or hdf file using a specified parser. The file ending
    is used to determine the correct parsing method.

    :param filename: Path to mzML / hdf file containing MS2 spectra to be loaded.
    :param parser: Name of the package to use for parsing the mzml file, can be "pyteomics" or "pymzml".
        Only used for parsing of mzML files.
    :param tims_meta_file: Optional path to timstof metadata file in internal format. This is only required
        when loading timstof spectra and used for summation of spectra.
    :raises ValueError: if the filename does not end in either ".hdf" or ".mzML" (case-insensitive)
    :raises AssertionError: if no tims_meta_file was provided when loading timsTOF hdf data
    :return: measured spectra with metadata.
    """
    if isinstance(filename, str):
        filename = Path(filename)

    format_ = filename.suffix.lower()
    if format_ == ".mzml":
        return ThermoRaw.read_mzml(
            source=filename, package=parser, search_type=""
        )  # TODO in spectrum_io, remove unnecessary argument
    elif format_ == ".hdf":
        if tims_meta_file is None:
            raise AssertionError(
                "Loading spectra from a timsTOF hdf file requires metadata provided by tims_meta_file."
            )
        results = read_and_aggregate_timstof(source=filename, tims_meta_file=Path(tims_meta_file))
        return results

    else:
        raise ValueError(f"Unrecognized file format: {format_}. Only .mzML and .hdf files are supported.")


def convert_raw_to_mzml(
    raw_file: Union[str, Path], output_file: Union[str, Path], thermo_exe: Optional[Union[str, Path]] = None
):
    """
    Convert raw to mzML format.

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


def convert_d_to_hdf(d_dir: Union[str, Path], output_file: Union[str, Path]):
    """
    Convert d to hdf format.

    This function converts spectra within a d folder from a mass spectrometry run into hdf format.

    :param d_dir: Path to d folder with spectra to be converted to hdf format
    :param output_file: Path to the location where the converted spectra should be written to
    """
    d_dir = Path(d_dir)
    output_file = Path(output_file)
    convert_d_hdf(input_path=d_dir, output_path=output_file)
