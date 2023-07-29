import logging
from pathlib import Path

import numpy as np
import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.annotation.annotation import annotate_spectra
from spectrum_fundamentals.fragments import compute_peptide_mass
from spectrum_fundamentals.mod_string import internal_without_mods, maxquant_to_internal
from spectrum_io.file import csv
from spectrum_io.raw import ThermoRaw
from spectrum_io.search_result import Mascot, MaxQuant, MSFragger
from spectrum_io.spectral_library import digest

from ..data.spectra import FragmentType, Spectra
from ..utils.config import Config
from ..utils.process_step import ProcessStep

logger = logging.getLogger(__name__)


# SpectralLibrary
def gen_lib(config: Config, library: Spectra) -> Spectra:
    """Read input csv file and add it to library."""
    library_input_type = config.library_input_type
    if library_input_type == "fasta":
        read_fasta(config, config.output)
        library_file = config.output / "prosit_input.csv"
    elif library_input_type == "peptides":
        library_file = config.library_input
    library_df = csv.read_file(library_file)
    library_df.columns = library_df.columns.str.upper()
    library.add_columns(library_df)
    return library


def read_fasta(config: Config, out_path: Path):
    """Read fasta file."""
    cmd = [
        "--fasta",
        f"{config.library_input}",
        "--prosit_input",
        f"{config.output / 'prosit_input.csv'}",
        "--fragmentation",
        f"{config.fragmentation}",
        "--digestion",
        f"{config.digestion}",
        "--cleavages",
        f"{config.cleavages}",
        "--db",
        f"{config.db}",
        "--enzyme",
        f"{config.enzyme}",
        "--special-aas",
        f"{config.special_aas}",
        "--min-length",
        f"{config.min_length}",
        "--max-length",
        f"{config.max_length}",
    ]
    digest.main(cmd)


def process_and_filter_spectra_data(config: Config, library: Spectra) -> Spectra:
    """
    Process and filter the spectra data in the given SpectralLibrary object.

    This function applies various modifications and filters to the 'spectra_data' DataFrame
    in the provided SpectralLibrary object. It modifies the 'MODIFIED_SEQUENCE' column,
    converts the 'MODIFIED_SEQUENCE' to internal format, extracts 'SEQUENCE', and filters
    out certain entries based on specific criteria.
    """
    library.spectra_data["MODIFIED_SEQUENCE"] = library.spectra_data["MODIFIED_SEQUENCE"].apply(lambda x: "_" + x + "_")
    library.spectra_data["MODIFIED_SEQUENCE"] = maxquant_to_internal(
        library.spectra_data["MODIFIED_SEQUENCE"], fixed_mods={}
    )
    library.spectra_data["SEQUENCE"] = internal_without_mods(library.spectra_data["MODIFIED_SEQUENCE"])
    library.spectra_data["PEPTIDE_LENGTH"] = library.spectra_data["SEQUENCE"].apply(lambda x: len(x))

    logger.info(f"No of sequences before Filtering is {len(library.spectra_data['PEPTIDE_LENGTH'])}")
    library.spectra_data = library.spectra_data[(library.spectra_data["PEPTIDE_LENGTH"] <= 30)]
    library.spectra_data = library.spectra_data[(~library.spectra_data["MODIFIED_SEQUENCE"].str.contains(r"\(ac\)"))]
    library.spectra_data = library.spectra_data[
        (~library.spectra_data["MODIFIED_SEQUENCE"].str.contains(r"\(Acetyl \(Protein N-term\)\)"))
    ]
    library.spectra_data = library.spectra_data[(~library.spectra_data["SEQUENCE"].str.contains("U"))]
    library.spectra_data = library.spectra_data[library.spectra_data["PRECURSOR_CHARGE"] <= 6]
    library.spectra_data = library.spectra_data[library.spectra_data["PEPTIDE_LENGTH"] >= 7]
    logger.info(f"No of sequences after Filtering is {len(library.spectra_data['PEPTIDE_LENGTH'])}")

    tmt_model = False
    for _, value in config.models.items():
        if value and "TMT" in value:
            tmt_model = True

    if tmt_model and config.tag != "":
        unimod_tag = c.TMT_MODS[config.tag]
        library.spectra_data["MODIFIED_SEQUENCE"] = maxquant_to_internal(
            library.spectra_data["MODIFIED_SEQUENCE"],
            fixed_mods={"C": "C[UNIMOD:4]", "^_": f"_{unimod_tag}", "K": f"K{unimod_tag}"},
        )
    else:
        library.spectra_data["MODIFIED_SEQUENCE"] = maxquant_to_internal(library.spectra_data["MODIFIED_SEQUENCE"])

    library.spectra_data["MASS"] = library.spectra_data["MODIFIED_SEQUENCE"].apply(lambda x: compute_peptide_mass(x))

    return library


# CeCalibration
def load_search(config: Config) -> pd.DataFrame:
    """Load search type."""
    switch = config.search_results_type
    logger.info(f"search_type is {switch}")
    if switch in ["maxquant", "msfragger", "mascot"]:
        search_path = _gen_internal_search_result_from_msms(config)
        switch = "internal"
    if switch == "internal":
        return csv.read_file(search_path)
    else:
        raise ValueError(f"{switch} is not a supported search type. Convert to internal format manually.")


def _gen_internal_search_result_from_msms(config: Config) -> Path:
    """Generate internal search result from msms.txt."""
    logger.info(f"Converting msms data at {config.search_results} to internal search result.")

    search_type = config.search_results_type
    if search_type == "maxquant":
        search_result = MaxQuant(config.search_results)
    elif search_type == "msfragger":
        search_result = MSFragger(config.search_results)
    elif search_type == "mascot":
        search_result = Mascot(config.search_results)
    else:
        raise ValueError(f"Unknown search_type provided in config: {search_type}")

    tmt_labeled = config.tag if any("TMT" in value for value in config.models.values()) else ""
    search_path = search_result.generate_internal(tmt_labeled=tmt_labeled)
    return search_path


# ReScore
def get_raw_files(config: Config):
    """
    Obtains raw files by scanning through the raw_path directory.

    If raw_path is a file, only process this one.
    :raises ValueError: raw_type is not supported as rawfile-type
    :raises FileNotFoundError: if raw file could not be found
    """
    raw_files = []
    if config.spectra.is_file():
        raw_files = [config.spectra]
    elif config.spectra.is_dir():
        spectra_type = config.spectra_type
        glob_pattern = _get_glob_pattern(spectra_type)
        raw_files = list(config.spectra.glob(glob_pattern))
        logger.info(f"Found {len(raw_files)} raw files in the search directory")
    else:
        raise FileNotFoundError(f"{config.spectra} does not exist.")


def _get_glob_pattern(spectra_type: str) -> str:
    """
    Get global pattern depending on spectra_type.

    :param spectra_type: The type of spectra file. Accepted values are "raw" or "mzml".
    :return: The glob pattern corresponding to the spectra_type.
    :raises ValueError: If an unsupported spectra_type is provided.
    """
    if spectra_type == "raw":
        return "*.[rR][aA][wW]"
    elif spectra_type == "mzml":
        return "*.[mM][zZ][mM][lL]"
    else:
        raise ValueError(f"{spectra_type} is not supported as rawfile-type")


def split_msms(config: Config, split_msms_step: ProcessStep):
    """Splits msms.txt file per raw file such that we can process each raw file in parallel \
    without reading the entire msms.txt."""
    out_path = config.output
    out_path.mkdir(exist_ok=True)

    if split_msms_step.is_done():
        return
    get_msms_folder_path(config).mkdir(exist_ok=True)

    df_search = load_search(config)
    logger.info(f"Read {len(df_search.index)} PSMs from {config.search_results}")
    for raw_file, df_search_split in df_search.groupby("RAW_FILE"):
        spectra_file_path = config.spectra / raw_file
        spectra_type = config.spectra_type
        glob_pattern = _get_glob_pattern(spectra_type)
        matched_files = [file for file in spectra_file_path.parent.glob(glob_pattern)]
        if not any(path.is_file() for path in matched_files):
            logger.info(f"Did not find {raw_file} in search directory, skipping this file")
            continue
        split_msms = _get_split_msms_path(config, raw_file + ".rescore")
        logger.info(f"Creating split msms.txt file {split_msms}")
        df_search_split = df_search_split[(df_search_split["PEPTIDE_LENGTH"] <= 30)]
        df_search_split = df_search_split[(~df_search_split["MODIFIED_SEQUENCE"].str.contains(r"\(ac\)"))]
        df_search_split = df_search_split[
            (~df_search_split["MODIFIED_SEQUENCE"].str.contains(r"\(Acetyl \(Protein N-term\)\)"))
        ]
        df_search_split = df_search_split[(~df_search_split["SEQUENCE"].str.contains("U"))]
        df_search_split = df_search_split[df_search_split["PRECURSOR_CHARGE"] <= 6]
        df_search_split = df_search_split[df_search_split["PEPTIDE_LENGTH"] >= 7]
        df_search_split.to_csv(split_msms, sep="\t", index=False)

    split_msms_step.mark_done()


def get_msms_folder_path(config: Config) -> Path:
    """Get folder path to msms."""
    return config.output / "msms"


def _get_split_msms_path(config: Config, raw_file: str) -> Path:
    """
    Get path to split msms.

    :param raw_file: path to raw file as a string
    :return: path to split msms file
    """
    return get_msms_folder_path(config) / raw_file


def merge_mzml_and_msms(config: Config, library: Spectra, df_search: pd.DataFrame) -> Spectra:
    """
    Read input search and mzml and add it to library.

    :param df_search: search result as pd.DataFrame
    """
    df_raw = _load_rawfile(config)
    # return df_search
    logger.info("Merging rawfile and search result")
    df_join = df_search.merge(df_raw, on=["RAW_FILE", "SCAN_NUMBER"])
    logger.info(f"There are {len(df_join)} matched identifications")
    logger.info("Annotating raw spectra")
    df_annotated_spectra = annotate_spectra(df_join)
    df_join.drop(columns=["INTENSITIES", "MZ"], inplace=True)
    # return df_annotated_spectra["INTENSITIES"]
    logger.info("Preparing library")
    library.add_columns(df_join)
    library.add_matrix(df_annotated_spectra["INTENSITIES"], FragmentType.RAW)
    library.add_matrix(df_annotated_spectra["MZ"], FragmentType.MZ)
    library.add_column(df_annotated_spectra["CALCULATED_MASS"], "CALCULATED_MASS")
    return library


def _load_rawfile(config: Config):
    """Load raw file."""
    switch = config.spectra_type
    search_engine = config.search_results_type
    logger.info(f"raw_type is {switch}")
    if switch == "raw":
        raw_path = _gen_mzml_from_thermo(config)
        switch = "mzml"
    if switch == "mzml":
        # TODO self.mzml_reader_package
        return ThermoRaw.read_mzml(source=raw_path, package="pyteomics", search_type=search_engine)
    else:
        raise ValueError(f"{switch} is not supported as rawfile-type")


def _gen_mzml_from_thermo(config: Config):
    """Generate mzml from thermo raw file."""
    logger.info("Converting thermo rawfile to mzml.")
    raw = ThermoRaw()
    raw_path = raw.convert_raw_mzml(
        input_path=config.spectra, output_path=get_mzml_path(config), thermo_exe=config.thermo_exe
    )
    return raw_path


def get_mzml_path(config: Config) -> Path:
    """Get path to mzml file."""
    if config.spectra_type == "mzml":
        spectra_path = config.spectra
    else:
        spectra_path = config.output / "mzML"
    spectra_path.mkdir(exist_ok=True)
    return spectra_path / config.spectra.with_suffix(".mzML").name


def prepare_alignment_df(library: Spectra) -> Spectra:
    """
    Prepare an alignment DataFrame from the given Spectra library.

    This function creates an alignment DataFrame by removing decoy and HCD fragmented spectra
    from the input library, selecting the top 1000 highest-scoring spectra, and repeating the
    DataFrame for each collision energy (CE) in the range [18, 50].
    """
    alignment_library = Spectra()
    alignment_library.spectra_data = library.spectra_data.copy()

    # Remove decoy and HCD fragmented spectra
    alignment_library.spectra_data = alignment_library.spectra_data[
        (alignment_library.spectra_data["FRAGMENTATION"] == "HCD") & (~alignment_library.spectra_data["REVERSE"])
    ]
    # Select the 1000 highest scoring or all if there are less than 1000
    alignment_library.spectra_data = alignment_library.spectra_data.sort_values(by="SCORE", ascending=False).iloc[:1000]

    # Repeat dataframe for each CE
    ce_range = range(18, 50)
    nrow = len(alignment_library.spectra_data)
    alignment_library.spectra_data = pd.concat([alignment_library.spectra_data for _ in ce_range], axis=0)
    alignment_library.spectra_data["COLLISION_ENERGY"] = np.repeat(ce_range, nrow)
    alignment_library.spectra_data.reset_index(inplace=True)
    return alignment_library
