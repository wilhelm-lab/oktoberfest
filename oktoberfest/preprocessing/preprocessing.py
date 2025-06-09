import logging
import re
from itertools import chain, combinations, product, repeat
from pathlib import Path
from sys import platform
from typing import Any, Optional, Union

import numpy as np
import pandas as pd
import spectrum_fundamentals.constants as c
from anndata import AnnData
from spectrum_fundamentals.annotation.annotation import annotate_spectra
from spectrum_fundamentals.fragments import compute_peptide_mass, retrieve_ion_types
from spectrum_fundamentals.mod_string import internal_without_mods, maxquant_to_internal
from spectrum_io.d import convert_d_hdf, read_and_aggregate_timstof
from spectrum_io.file import csv
from spectrum_io.raw import ThermoRaw
from spectrum_io.search_result import Mascot, MaxQuant, MSAmanda, MSFragger, OpenMS, Sage, Scout, Xisearch
from spectrum_io.spectral_library.digest import get_peptide_to_protein_map

from ..data.spectra import FragmentType, Spectra

logger = logging.getLogger(__name__)


# SpectralLibrary


def gen_lib(input_file: Union[str, Path]) -> Spectra:
    r"""
    Generate a spectral library from a given input.

    This function reads an input file that follows the specifications provided in the usage section for
    :doc:`../../peptides_format`, creates a `Spectra` object from it and returns it.

    :param input_file: A csv file containing modified sequences, CE, precursor charge and fragmentation method.
    :return: Spectra object from the read input.

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> import pandas as pd
        >>> peptides_df = pd.DataFrame({'modified_sequence': ['KIIDRAITSL', 'KIEKLKVEL', 'KINQQKLKL', 'TYDDATKTFTVTE', 'ASPTQPIQL'],
        >>>                             'collision_energy': [34, 35, 34, 34, 35],
        >>>                             'precursor_charge': [2, 2, 3, 2, 1]})
        >>> peptides_df.to_csv("./tests/doctests/input/peptides.csv", sep='\t',index=False)
        >>> library = pp.gen_lib("./tests/doctests/input/peptides.csv")
        >>> print(library)
    """
    library_df = csv.read_file(input_file)
    library_df.columns = library_df.columns.str.upper()
    if "PROTEINS" not in library_df.columns:
        library_df["PROTEINS"] = "unknown"
    var_df = Spectra._gen_vars_df()
    spec = Spectra(obs=library_df, var=var_df)

    spec.var_names = var_df.index
    return spec


def generate_metadata(
    peptides: list[str],
    collision_energy: Union[int, list[int]],
    precursor_charge: Union[int, list[int]],
    fragmentation: Union[str, list[str]],
    nr_ox: int,
    instrument_type: Optional[str] = None,
    proteins: Optional[list[list[str]]] = None,
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
    :param nr_ox: Maximal number of allowed oxidations.
    :param instrument_type: The type of mass spectrometeter. Only required when predicting intensities
        with AlphaPept. Choose one of ["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"].
    :param proteins: An optional list of proteins associated with each peptide.
        If provided, it must have the same length as the number of peptides.
    :raises AssertionError: If the lengths of peptides and proteins is not the same.
    :return: A DataFrame containing metadata with the columns "modified_peptide","collision_energy",
        "precursor_charge","fragmentation" and an optional "proteins" column.

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> metadata = pp.generate_metadata(peptides=["AAACRFVQ","RMPCHKPYL"],
        >>>                                             collision_energy=[30,35],
        >>>                                             precursor_charge=[1,2],
        >>>                                             fragmentation=["HCD","HCD"],
        >>>                                             nr_ox=1)
        >>> print(metadata)
    """
    if isinstance(collision_energy, int):
        collision_energy = [collision_energy]
    if isinstance(precursor_charge, int):
        precursor_charge = [precursor_charge]
    if isinstance(fragmentation, str):
        fragmentation = [fragmentation]

    if proteins is not None and len(proteins) != len(peptides):
        raise AssertionError("Number of proteins must match the number of peptides.")

    combinations_product = product(peptides, collision_energy, precursor_charge, fragmentation)

    metadata = pd.DataFrame(
        combinations_product, columns=["modified_sequence", "collision_energy", "precursor_charge", "fragmentation"]
    )
    metadata["peptide_length"] = metadata["modified_sequence"].str.len()
    metadata["instrument_types"] = instrument_type

    if proteins is not None:
        n_repeats = len(metadata) // len(proteins)
        metadata["proteins"] = list(
            chain.from_iterable([repeat(";".join(prot_list), n_repeats) for prot_list in proteins])
        )
    else:
        metadata["proteins"] = "unknown"

    modified_peptides = []
    for _, row in metadata.iterrows():
        peptide = row["modified_sequence"]
        res = [i.start() for i in re.finditer("M", peptide)]
        res.reverse()
        for i in range(1, min(len(res), nr_ox) + 1):
            possible_indices = list(combinations(res, i))
            for index in possible_indices:
                string_mod = peptide
                for j in index:
                    string_mod = string_mod[: j + 1] + "[UNIMOD:35]" + string_mod[j + 1 :]
                new_row = row.copy()
                new_row["modified_sequence"] = string_mod
                modified_peptides.append(new_row)
    metadata = pd.concat([metadata, pd.DataFrame(modified_peptides)], ignore_index=True)

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
) -> dict[str, list[str]]:
    r"""
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

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> peptides = [
        >>>     (">Peptide1 Example peptide 1", "MKTIIALSYIFCLVFAD"),
        >>>     (">Peptide2 Example peptide 2", "GILGFVFRTLTVPS"),
        >>>     (">Peptide3 Example peptide 3", "LLGATCMFV")
        >>> ]
        >>> with open("./tests/doctests/input/peptides.fasta", "w") as file:
        >>>     for header, sequence in peptides:
        >>>         file.write(f"{header}\n")
        >>>         file.write(f"{sequence}\n")
        >>> digest_dict = pp.digest(fasta="./tests/doctests/input/peptides.fasta",
        >>>                         digestion="full",
        >>>                         missed_cleavages=2,
        >>>                         db="concat",
        >>>                         enzyme="trypsin",
        >>>                         special_aas="KR",
        >>>                         min_length=7,
        >>>                         max_length=60)
        >>> print(digest_dict)
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


def filter_peptides_for_model(peptides: Union[pd.DataFrame, AnnData], model: str) -> Union[pd.DataFrame, AnnData]:
    """
    Filter search results to support a given peptide prediction model.

    Depending on the model used for peptide property prediction, PSMs in search results need to be filtered out.
    This function provides a shortcut using the model name for filtering instead of supplying all possible filters manually.

    :param peptides: Dataframe containing search results to be filtered
    :param model: Name of a peptide property prediction model to use as filter

    :raises ValueError: if an unsupported model is supplied

    :return: The filtered dataframe or AnnData object to be used with the given model.

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> import pandas as pd
        >>> search_results = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL","TAIASPEK"],
        >>>                     "SEQUENCE": ["AAACRFVQ","RMPCHKPYL","TAIASPEK"],
        >>>                     "PEPTIDE_LENGTH": [8,9,8],
        >>>                     "PRECURSOR_CHARGE": [1,2,7]})
        >>> filtered_peptides = pp.filter_peptides_for_model(peptides=search_results, model="prosit")
        >>> print(filtered_peptides)
    """
    if "prosit" in model.lower() and "xl" in model.lower():
        filter_kwargs = {
            "min_length": 6,
            "max_length": 30,
            "max_charge": 6,
        }
        return filter_xl_peptides(peptides, **filter_kwargs)
    if "prosit" in model.lower():
        filter_kwargs = {
            "min_length": 7,
            "max_length": 30,
            "max_charge": 6,
        }
    elif "ms2pip" in model.lower():
        filter_kwargs = {
            "min_length": 2,
            "max_length": 100,
            "max_charge": 6,
        }
    elif "alphapept" in model.lower():
        filter_kwargs = {
            "min_length": 7,
            "max_length": 35,
            "max_charge": 4,
        }
    else:
        raise ValueError(f"The model {model} is not known.")

    return filter_peptides(peptides, **filter_kwargs)


def filter_peptides(
    peptides: Union[pd.DataFrame, AnnData], min_length: int, max_length: int, max_charge: int
) -> Union[pd.DataFrame, AnnData]:
    """
    Filter search results using given constraints.

    This function filters provided search results by peptide length, precursor charge,
    unsupported special aminoacids, and unsupported modifications.

    :param peptides: Dataframe containing search results to be filtered
    :param min_length: The minimal length of a peptide to be retained
    :param max_length: The maximal length of a peptide to be retained
    :param max_charge: The maximal precursor charge of a peptide to be retained

    :return: The filtered dataframe or AnnData object given the provided constraints.

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> import pandas as pd
        >>> search_results = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL","TAIASPEK"],
        >>>                     "SEQUENCE": ["AAACRFVQ","RMPCHKPYL","TAIASPEK"],
        >>>                     "PEPTIDE_LENGTH": [8,9,8],
        >>>                     "PRECURSOR_CHARGE": [1,2,7]})
        >>> filtered_peptides = pp.filter_peptides(peptides=search_results, min_length=7, max_length=30, max_charge=6)
        >>> print(filtered_peptides)
    """
    if isinstance(peptides, AnnData):
        df = peptides.obs
    else:
        df = peptides
    peptide_filter = (
        (df["PEPTIDE_LENGTH"] <= max_length)
        & (df["PEPTIDE_LENGTH"] >= min_length)
        & (df["PRECURSOR_CHARGE"] <= max_charge)
        & (~df["SEQUENCE"].str.contains(r"B|\*|\.|O|U|X|Z"))
    )
    return peptides[peptide_filter.values]


def filter_xl_peptides(peptides: pd.DataFrame, min_length: int, max_length: int, max_charge: int) -> pd.DataFrame:
    """
    Filter xl search results using given constraints.

    This function filters provided search results by peptide length, precursor charge,
    unsupported special aminoacids, and unsupported modifications.

    :param peptides: Dataframe containing search results to be filtered
    :param min_length: The minimal length of a peptide to be retained
    :param max_length: The maximal length of a peptide to be retained
    :param max_charge: The maximal precursor charge of a peptide to be retained

    :return: The filtered dataframe or AnnData object given the provided constraints.
    """
    if isinstance(peptides, AnnData):
        df = peptides.obs
    else:
        df = peptides

    peptide_filter = (
        (df["PEPTIDE_LENGTH_A"] <= max_length)
        & (df["PEPTIDE_LENGTH_B"] <= max_length)
        & (df["PEPTIDE_LENGTH_A"] >= min_length)
        & (df["PEPTIDE_LENGTH_B"] >= min_length)
        & (df["PRECURSOR_CHARGE"] <= max_charge)
        & (df["MODIFIED_SEQUENCE_A"].str.contains(r"\[UNIMOD\:1896\]|\[UNIMOD\:1884\]"))
        & (df["MODIFIED_SEQUENCE_B"].str.contains(r"\[UNIMOD\:1896\]|\[UNIMOD\:1884\]"))
    )

    return peptides[peptide_filter]


def process_and_filter_spectra_data(library: Spectra, model: str, tmt_label: Optional[str] = None) -> Spectra:
    """
    Process and filter the spectra data in the given SpectralLibrary object.

    This function applies various modifications and filters to the obs DataFrame
    in the provided SpectralLibrary object. It modifies the 'MODIFIED_SEQUENCE' column,
    converts the 'MODIFIED_SEQUENCE' to internal format, extracts 'SEQUENCE', and filters
    out certain entries based on specific criteria. The specification of the internal file format can be found at
    :doc:`../../internal_format`.

    :param library: A Spectra object containing the raw peptides to be proccessed and filtered.
    :param model: The peptide property prediction model to filter the spectra for
    :param tmt_label: Optional tmt-label to consider when processing peptides. If given, the corresponding
        fixed modification for the N-terminus and lysin will be added
    :return: The processed and filtered Spectra object

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> import pandas as pd
        >>> meta_df = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL","TAIASPEK"],
        >>>                         "PRECURSOR_CHARGE": [1,2,4]})
        >>> var = Spectra._gen_vars_df()
        >>> peptide_library = Spectra(obs=meta_df, var=var)
        >>> processed_peptides = pp.process_and_filter_spectra_data(library=peptide_library, model="prosit")
        >>> print(processed_peptides)
    """
    # add fixed mods and translate to internal format
    library.obs["MODIFIED_SEQUENCE"] = library.obs["MODIFIED_SEQUENCE"].apply(lambda x: "_" + x + "_")

    fixed_mods = {"C": "C[UNIMOD:4]"}
    if tmt_label is not None and tmt_label != "":
        unimod_tag = c.TMT_MODS[tmt_label]
        fixed_mods = {"C": "C[UNIMOD:4]", "^_": f"_{unimod_tag}-", "K": f"K{unimod_tag}"}

    # we use this method since we expect the input to be similar to MQ in that fixed modifications are
    # not written. This needs to be changed once we allow arbitrary modifications for the spectral library
    # generation, not just a number of oxidations and fixed carbamidomethylation / + TMT.
    library.obs["MODIFIED_SEQUENCE"] = maxquant_to_internal(library.obs["MODIFIED_SEQUENCE"], mods=fixed_mods)

    # get sequence and its length
    library.obs["SEQUENCE"] = internal_without_mods(library.obs["MODIFIED_SEQUENCE"])
    library.obs["PEPTIDE_LENGTH"] = library.obs["SEQUENCE"].apply(lambda x: len(x))

    # filter
    library = filter_peptides_for_model(library, model)
    library.obs["MASS"] = library.obs["MODIFIED_SEQUENCE"].apply(lambda x: compute_peptide_mass(x))

    return library


def load_search(
    input_file: Union[str, Path],
) -> pd.DataFrame:
    """
    Load search results.

    Given a path to a file containing search results in Oktoberfest format, the function reads the results and returns them.
    The specification of the internal file format can be found at :doc:`../../internal_format`.

    :param input_file: Path to the file containing search results in the internal Oktoberfest format.
    :return: dataframe containing the search results.

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> search_results = pd.DataFrame({"MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL","TAIASPEK"],
        >>>                     "SEQUENCE": ["AAACRFVQ","RMPCHKPYL","TAIASPEK"],
        >>>                     "PEPTIDE_LENGTH": [8,9,8],
        >>>                     "PRECURSOR_CHARGE": [1,2,7]})
        >>> search_results.to_csv("./tests/doctests/input/search_results.csv",index=False)
        >>> library = pp.load_search("./tests/doctests/input/search_results.csv")
        >>> print(library)
    """
    search_results = csv.read_file(input_file)
    return search_results


def convert_search(
    input_path: Union[str, Path],
    search_engine: str,
    tmt_label: str = "",
    custom_mods: Optional[dict[str, int]] = None,
    output_file: Optional[Union[str, Path]] = None,
    ptm_unimod_id: Optional[int] = 0,
    ptm_sites: Optional[list] = None,
) -> pd.DataFrame:
    r"""
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
    :param custom_mods: Optional dictionary parameter given when input_file is not in internal Oktoberfest format with
        static and variable mods as keys. The values are the integer values of the respective unimod identifier
    :param output_file: Optional path to the location where the converted search results should be written to.
        If this is omitted, the results are not stored.
    :param ptm_unimod_id: unimod id used for site localization
    :param ptm_sites: possible sites that the ptm can exist on
    :raises ValueError: if an unsupported search engine was given
    :return: A dataframe containing the converted results.

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> import pandas as pd
        >>> msms = pd.DataFrame({'Raw file': ['GN20170722_SK_HLA_G0103_R1_01', 'GN20170722_SK_HLA_G0103_R2_02'],
        >>> 'Scan number': [21329, 20501],
        >>> 'Scan index': [18847, 17998],
        >>> 'Sequence': ['AAAAVVSGPKRGRKKP', 'AAAAVVSGPKRGRKKP'],
        >>> 'Length': [16, 16],
        >>> 'Missed cleavages': ['', ''],
        >>> 'Modifications': ['Unmodified', 'Unmodified'],
        >>> 'Modified sequence': ['_AAAAVVSGPKRGRKKP_', '_AAAAVVSGPKRGRKKP_'],
        >>> 'Oxidation (M) Probabilities': ['', ''],
        >>> 'Oxidation (M) Score Diffs': ['', ''],
        >>> 'Oxidation (M)': [0, 0],
        >>> 'Proteins': ['', ''],
        >>> 'Charge': [3, 3],
        >>> 'Fragmentation': ['HCD', 'HCD'],
        >>> 'Mass analyzer': ['FTMS', 'FTMS'],
        >>> 'Type': ['MULTI-SECPEP', 'MULTI-SECPEP'],
        >>> 'Scan event number': [9, 5],
        >>> 'Isotope index': [2, 2],
        >>> 'm/z': [531.66176, 531.66176],
        >>> 'Mass': [1591.9634, 1591.9634],
        >>> 'Mass Error [ppm]': [-2.1109999999999998, -1.1018],
        >>> 'Simple Mass Error [ppm]': [1259.2803, 1259.2803],
        >>> 'Retention time': [46.272, 46.388000000000005],
        >>> 'PEP': [0.57389, 0.57389],
        >>> 'Score': [7.9138, 4.7582],
        >>> 'Delta score': [3.5652, 1.4401],
        >>> 'Score diff': ['', ''],
        >>> 'Localization prob': [1, 1],
        >>> 'Combinatorics': [0, 0],
        >>> 'PIF': [0, 0],
        >>> 'Fraction of total spectrum': [0, 0],
        >>> 'Base peak fraction': [0, 0],
        >>> 'Precursor Full ScanNumber': [-1, -1],
        >>> 'Precursor Intensity': [0, 0],
        >>> 'Precursor Apex Fraction': [0, 0],
        >>> 'Precursor Apex Offset': [0, 0],
        >>> 'Precursor Apex Offset Time': [0, 0],
        >>> 'Matches Intensities': ['y5;y10;y5-NH3;a2;b12(2+)', 'y5;y5-NH3;b12(2+)'],
        >>> 'Mass Deviations [Da]': ['34666.4;2191.7;88570.6;2148.7;89073.6', '10544.1;36224.8;73327.7'],
        >>> 'Mass Deviations [ppm]': ['0.008335659;-0.01799215;-0.002397317;-0.0004952438;-0.004926575',
        >>>                             '0.009286639;-0.004650567;-0.002822918'],
        >>> 'Masses': ['14.23987;-16.19888;-4.217963;-4.303209;-9.237617', '15.86446;-8.182414;-5.293158'],
        >>> 'Number of Matches': [5, 3],
        >>> 'Intensity coverage': [0.1016966, 0.1564349],
        >>> 'Peak coverage': [0.04166667, 0.04477612],
        >>> 'Neutral loss level': ['None', 'None'],
        >>> 'ETD identification type': ['Unknown', 'Unknown'],
        >>> 'Reverse': ['Unknown +', 'Unknown +'],
        >>> 'All scores': ['7.913836;4.348669;4.097387', '4.758178;3.318045;2.968256'],
        >>> 'All sequences': ['AAAAVVSGPKRGRKKP;GVVAKGALTPKLSPVVG;GVVPSLKPTLAGKAVVG',
        >>>                     'AAAAVVSGPKRGRKKP;VMKLLRHDKLVQL;QEILRKILPLGELA'],
        >>> 'All modified sequences': ['_AAAAVVSGPKRGRKKP_;_GVVAKGALTPKLSPVVG_;_GVVPSLKPTLAGKAVVG_',
        >>>                             '_AAAAVVSGPKRGRKKP_;_VMKLLRHDKLVQL_;_QEILRKILPLGELA_'],
        >>> 'id': [1378, 1379],
        >>> 'Protein group IDs': ['42625', '42625'],
        >>> 'Peptide ID': [533, 533],
        >>> 'Mod. peptide ID': [537, 537],
        >>> 'Evidence ID': [1075, 1076],
        >>> 'Oxidation (M) site IDs': ['', '']})
        >>> msms.to_csv("./tests/doctests/input/msms.txt",sep='\t',index=False)
        >>> converted_results = pp.convert_search(input_path="./tests/doctests/input/", search_engine="maxquant")
        >>> print(converted_results)
    """
    search_engine = search_engine.lower()
    search_result: Any
    xl = False
    if search_engine == "maxquant":
        search_result = MaxQuant
    elif search_engine == "msfragger":
        search_result = MSFragger
    elif search_engine == "mascot":
        search_result = Mascot
    elif search_engine == "sage":
        search_result = Sage
    elif search_engine == "openms":
        search_result = OpenMS
    elif search_engine == "xisearch":
        search_result = Xisearch
        xl = True
    elif search_engine == "scout":
        search_result = Scout
        xl = True
    elif search_engine == "msamanda":
        search_result = MSAmanda
    else:
        raise ValueError(f"Unknown search engine provided: {search_engine}")

    return search_result(input_path).generate_internal(
        tmt_label=tmt_label,
        out_path=output_file,
        custom_mods=custom_mods,
        ptm_unimod_id=ptm_unimod_id,
        ptm_sites=ptm_sites,
        xl=xl,
    )


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


def list_spectra(input_dir: Union[str, Path], input_format: str) -> list[Path]:
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

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> import os
        >>> # creating minimum viable example .mzml file
        >>> filecontent = '''<?xml version="1.0" encoding="UTF-8"?>
        >>> <mzML xmlns="http://example" version="1.1.0">
        >>>   <cvList count="2">
        >>>     <cv id="MS" fullName="Mass Spectrometry Ontology" version="4.1.0" URI="https://example"/>
        >>>     <cv id="UO" fullName="Unit Ontology" version="1.23" URI="http://example"/>
        >>>   </cvList>
        >>>   <fileDescription>
        >>>     <fileContent>
        >>>       <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum"/>
        >>>     </fileContent>
        >>>   </fileDescription>
        >>>   <referenceableParamGroupList count="1">
        >>>     <referenceableParamGroup id="commonInstrumentParams">
        >>>       <cvParam cvRef="MS" accession="MS:1000031" name="instrument model" value="Example Instrument"/>
        >>>     </referenceableParamGroup>
        >>>   </referenceableParamGroupList>
        >>>   <run id="run1" defaultInstrumentConfigurationRef="IC1">
        >>>     <spectrumList count="1">
        >>>       <spectrum index="0" id="scan=1" defaultArrayLength="5">
        >>>         <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
        >>>         <binaryDataArrayList count="2">
        >>>           <binaryDataArray encodedLength="20">
        >>>             <cvParam cvRef="MS" accession="4" name="m/z" unitCvRef="MS" unitAccession="0" unitName="m/z"/>
        >>>             <binary>...</binary>
        >>>           </binaryDataArray>
        >>>           <binaryDataArray encodedLength="20">
        >>>             <cvParam cvRef="MS" accession="5" name="i" unitCvRef="MS" unitAccession="1" unitName="c"/>
        >>>             <binary>...</binary>
        >>>           </binaryDataArray>
        >>>         </binaryDataArrayList>
        >>>       </spectrum>
        >>>     </spectrumList>
        >>>   </run>
        >>> </mzML>'''
        >>> os.makedirs("./tests/doctests/input/spectra", exist_ok=True)
        >>> with open("./tests/doctests/input/spectra/File1.mzml","w+") as f:
        >>>     f.writelines(filecontent)
        >>> paths = pp.list_spectra(input_dir="./tests/doctests/input/spectra/", input_format="mzml")
        >>> print(paths)
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
    search_results: pd.DataFrame, output_dir: Union[str, Path], filenames: Optional[list[str]] = None
) -> list[str]:
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

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> import pandas as pd
        >>> search_results = pd.DataFrame({"RAW_FILE": ["File1","File2"],
        >>>                             "SCAN_NUMBER": [5123,4012],
        >>>                             "MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
        >>>                             "PRECURSOR_CHARGE": [1,2],
        >>>                             "SCAN_EVENT_NUMBER": [4,10],
        >>>                             "MASS": [1000.41,1589.1],
        >>>                             "SCORE": [3.64,5.45],
        >>>                             "REVERSE": [False,False],
        >>>                             "SEQUENCE": ["AAACRFVQ","RMPCHKPYL"],
        >>>                             "PEPTIDE_LENGTH": [8,9]})
        >>> pp.split_search(search_results=search_results, output_dir="./tests/doctests/output/")
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
    timstof_metadata: pd.DataFrame, output_dir: Union[str, Path], filenames: Optional[list[str]] = None
) -> list[str]:
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

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> import pandas as pd
        >>> timstod_meta = pd.DataFrame({"RAW_FILE": ["220331_NHG_malignant_CLL_02_Tue39L243_17%_DDA_Rep3",
        >>>                                 "220331_NHG_malignant_CLL_02_Tue39L243_17%_DDA_Rep4"],
        >>>                                 "FRAME": [2733,2824],
        >>>                                 "PRECURSOR": [2195,2299],
        >>>                                 "SCAN_NUM_BEGIN": [1416,1488],
        >>>                                 "SCAN_NUM_END": [1439,1511],
        >>>                                 "COLLISION_ENERGY": [26.25,24.57],
        >>>                                 "SCAN_NUMBER": [8646,4879]})
        >>> pp.split_timstof_metadata(timstof_metadata=timstod_meta, output_dir="./tests/doctests/output/")
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


def merge_spectra_and_peptides(spectra: pd.DataFrame, search: pd.DataFrame) -> pd.DataFrame:
    """
    Merge peptides with spectra.

    This function takes spectra and and peptides from a search and merges them based on the
    RAW File identifier and the scan number using the "RAW_FILE" and "SCAN_NUMBER" column of
    the two dataframes provided.

    :param spectra: MS2 spectra of a mass spectrometry run
    :param search: Peptides from a search for the given spectra in internal format (see
        :doc:`../../internal_format`)

    :return: Dataframe containing the matched pairs of peptides and spectra (PSMs)

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> import numpy as np
        >>> import pandas as pd
        >>> search_results = pd.DataFrame({"RAW_FILE": ["File1","File2"],
        >>>                                 "SCAN_NUMBER": [5123,4012],
        >>>                                 "MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
        >>>                                 "PRECURSOR_CHARGE": [1,2],
        >>>                                 "SCAN_EVENT_NUMBER": [4,10],
        >>>                                 "MASS": [1000.41,1589.1],
        >>>                                 "SCORE": [3.64,5.45],
        >>>                                 "REVERSE": [False,False],
        >>>                                 "SEQUENCE": ["AAACRFVQ","RMPCHKPYL"],
        >>>                                 "PEPTIDE_LENGTH": [8,9]})
        >>> spectra = pd.DataFrame({"RAW_FILE": ["File1","File2"],
        >>>                         "SCAN_NUMBER": [5123,4012],
        >>>                         "INTENSITIES": [np.random.rand(174),np.random.rand(174)],
        >>>                         "MZ": [np.random.rand(174),np.random.rand(174)],
        >>>                         "MZ_RANGE": ["100.0-385.0","100.0-402.0"],
        >>>                         "RETENTION_TIME": [59.1, 110.42],
        >>>                         "MASS_ANALYZER": ["FTMS","FTMS"],
        >>>                         "FRAGMENTATION": ["HCD","HCD"],
        >>>                         "COLLISION_ENERGY": [27,27],
        >>>                         "INSTRUMENT_TYPES": ["Q Exactive Plus","Q Exactive Plus"]})
        >>> psms = pp.merge_spectra_and_peptides(spectra=spectra, search=search_results)
        >>> print(psms)
    """
    logger.info("Merging rawfile and search result")
    psms = search.merge(spectra, on=["RAW_FILE", "SCAN_NUMBER"])
    return psms


def annotate_spectral_library(
    psms: pd.DataFrame,
    fragmentation_method: str = "HCD",
    mass_tol: Optional[float] = None,
    unit_mass_tol: Optional[str] = None,
    custom_mods: Optional[dict[str, float]] = None,
    annotate_neutral_loss: Optional[bool] = False,
) -> Spectra:
    """
    Annotate all specified ion peaks of given PSMs (Default b and y ions).

    This function annotates the b any ion peaks of given psms by matching the mzs
    of all peaks to the theoretical mzs and discards all other peaks. It also calculates
    the theoretical monoisotopic mass of each b and y ion fragment.
    The function thenr returns a Spectra object containing the mzs and intensities of
    all b and y ions in charge states 1-3 and the additional metadata.

    :param psms: Spectral library to be annotated.
    :param mass_tol: The mass tolerance allowed for retaining peaks
    :param unit_mass_tol: The unit in which the mass tolerance is given
    :param fragmentation_method: fragmentation method that was used
    :param custom_mods: mapping of custom UNIMOD string identifiers ('[UNIMOD:xyz]') to their mass
    :param annotate_neutral_loss: flag to indicate whether to annotate neutral loss peaks or not

    :return: Spectra object containing the annotated b and y ion peaks including metadata

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> import pandas as pd
        >>> psms = pd.DataFrame({"RAW_FILE": ["File1","File2"],
        >>>                     "SCAN_NUMBER": [5123,4012],
        >>>                     "MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
        >>>                     "PRECURSOR_CHARGE": [1,2],
        >>>                     "PEPTIDE_LENGTH": [8,9],
        >>>                     "MASS_ANALYZER": ["FTMS","FTMS"],
        >>>                     "INTENSITIES": [np.random.rand(174),np.random.rand(174)],
        >>>                     "MZ": [np.random.rand(174),np.random.rand(174)]})
        >>> library = pp.annotate_spectral_library(psms=psms, mass_tol=15, unit_mass_tol="ppm")
        >>> print(library)
    """
    logger.info("Annotating spectra...")
    df_annotated_spectra = annotate_spectra(
        un_annot_spectra=psms,
        mass_tolerance=mass_tol,
        unit_mass_tolerance=unit_mass_tol,
        fragmentation_method=fragmentation_method,
        custom_mods=custom_mods,
        annotate_neutral_loss=annotate_neutral_loss,
    )

    ion_types = retrieve_ion_types(fragmentation_method)
    var_df = Spectra._gen_vars_df(ion_types)
    aspec = Spectra(obs=psms.drop(columns=["INTENSITIES", "MZ"]), var=var_df)
    aspec.uns["ion_types"] = ion_types
    aspec.add_intensities(
        np.stack(df_annotated_spectra["INTENSITIES"]), aspec.var_names.values[None, ...], FragmentType.RAW
    )
    aspec.add_mzs(np.stack(df_annotated_spectra["MZ"]), FragmentType.MZ)
    aspec.add_column(df_annotated_spectra["CALCULATED_MASS"].values, "CALCULATED_MASS")
    aspec.add_column(df_annotated_spectra["EXPECTED_NL_COUNT"].values, "EXPECTED_NL_COUNT")
    aspec.add_column(df_annotated_spectra["ANNOTATED_NL_COUNT"].values, "ANNOTATED_NL_COUNT")
    aspec.strings_to_categoricals()

    logger.info("Finished annotating.")

    return aspec


def annotate_spectral_library_xl(
    psms: pd.DataFrame, mass_tol: Optional[float] = None, unit_mass_tol: Optional[str] = None
):
    """
    Annotate spectral library with peaks and mass for cross-linked peptides.

    This function annotates a given spectral library with peak intensities and mass to charge ratio,
    as well as the calculated monoisotopic mass of the precursor ion.
    The additional information is added to the provided spectral library.

    :param psms: Spectral library to be annotated.
    :param mass_tol: The mass tolerance allowed for retaining peaks
    :param unit_mass_tol: The unit in which the mass tolerance is given
    :return: Spectra object containing the annotated b and y ion peaks including metadata
    """
    logger.info("Annotating spectra...")
    df_annotated_spectra = annotate_spectra(psms, mass_tol, unit_mass_tol)
    aspec = Spectra(obs=psms.drop(columns=["INTENSITIES", "MZ"]), var=Spectra._gen_vars_df(xl=True))
    aspec.add_intensities(
        np.stack(df_annotated_spectra["INTENSITIES_A"]), aspec.var_names.values[None, ...], FragmentType.RAW_A
    )
    aspec.add_intensities(
        np.stack(df_annotated_spectra["INTENSITIES_B"]), aspec.var_names.values[None, ...], FragmentType.RAW_B
    )
    aspec.add_mzs(np.stack(df_annotated_spectra["MZ_A"]), FragmentType.MZ_A)
    aspec.add_mzs(np.stack(df_annotated_spectra["MZ_B"]), FragmentType.MZ_B)
    aspec.add_column(df_annotated_spectra["CALCULATED_MASS_A"].values, "CALCULATED_MASS_A")
    aspec.add_column(df_annotated_spectra["CALCULATED_MASS_B"].values, "CALCULATED_MASS_B")

    logger.info("Finished annotating.")
    return aspec


def load_spectra(
    filenames: Union[str, Path, list[Union[str, Path]]],
    parser: str = "pyteomics",
    tims_meta_file: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Read spectra from a given file.

    This function reads MS2 spectra from a given mzML or hdf file using a specified parser. The file ending
    is used to determine the correct parsing method.

    :param filenames: Path(s) to files containing MS2 spectra. Filenames need to end in ".mzML" (case-insensitive).
        For timstof data, a single hdf5 path ending in ".hdf" (case-insensitive) needs to be provided.
        Multiple paths are not yet supported for timstof.
    :param parser: Name of the package to use for parsing the mzml file, can be "pyteomics" or "pymzml".
        Only used for parsing of mzML files.
    :param tims_meta_file: Optional path to timstof metadata file in internal format. This is only required
        when loading timstof spectra and used for summation of spectra.
    :raises TypeError: if not all filenames are provided as str or Path objects.
    :raises ValueError: if the filename does not end in either ".hdf" or ".mzML" (case-insensitive)
    :raises AssertionError: if no tims_meta_file was provided when loading timsTOF hdf data
    :return: measured spectra with metadata.

    :Example:

    .. code-block:: python

        >>> from oktoberfest import preprocessing as pp
        >>> filecontent = '''<?xml version="1.0" encoding="UTF-8"?>
        >>> <mzML xmlns="http://example" version="1.1.0">
        >>>   <cvList count="2">
        >>>     <cv id="MS" fullName="Mass Spectrometry Ontology" version="4.1.0" URI="https://example"/>
        >>>     <cv id="UO" fullName="Unit Ontology" version="1.23" URI="http://example"/>
        >>>   </cvList>
        >>>   <fileDescription>
        >>>     <fileContent>
        >>>       <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum"/>
        >>>     </fileContent>
        >>>   </fileDescription>
        >>>   <referenceableParamGroupList count="1">
        >>>     <referenceableParamGroup id="commonInstrumentParams">
        >>>       <cvParam cvRef="MS" accession="MS:1000031" name="instrument model" value="Example Instrument"/>
        >>>     </referenceableParamGroup>
        >>>   </referenceableParamGroupList>
        >>>   <run id="run1" defaultInstrumentConfigurationRef="IC1">
        >>>     <spectrumList count="1">
        >>>       <spectrum index="0" id="scan=1" defaultArrayLength="5">
        >>>         <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
        >>>         <binaryDataArrayList count="2">
        >>>           <binaryDataArray encodedLength="20">
        >>>             <cvParam cvRef="MS" accession="4" name="m/z" unitCvRef="MS" unitAccession="0" unitName="m/z"/>
        >>>             <binary>...</binary>
        >>>           </binaryDataArray>
        >>>           <binaryDataArray encodedLength="20">
        >>>             <cvParam cvRef="MS" accession="5" name="i" unitCvRef="MS" unitAccession="1" unitName="c"/>
        >>>             <binary>...</binary>
        >>>           </binaryDataArray>
        >>>         </binaryDataArrayList>
        >>>       </spectrum>
        >>>     </spectrumList>
        >>>   </run>
        >>> </mzML>'''
        >>> with open("./tests/doctests/input/File1.mzml","w+") as f:
        >>>     f.writelines(filecontent)
        >>> spectra = pp.load_spectra(filenames=["./tests/doctests/input/File1.mzml"], parser="pyteomics")
        >>> print(spectra)
    """
    if isinstance(filenames, (str, Path)):
        internal_filenames = [Path(filenames)]
    elif isinstance(filenames, list):
        internal_filenames = [Path(filename) for filename in filenames]
    else:
        raise TypeError("Type of filenames not understood.")

    format_ = internal_filenames[0].suffix.lower()
    if format_ == ".mzml":
        return ThermoRaw.read_mzml(source=filenames, package=parser)
    elif format_ == ".hdf":
        if tims_meta_file is None:
            raise AssertionError(
                "Loading spectra from a timsTOF hdf file requires metadata provided by tims_meta_file."
            )
        results = read_and_aggregate_timstof(source=internal_filenames[0], tims_meta_file=Path(tims_meta_file))
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
