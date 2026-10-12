from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import mokapot
import numpy as np
import pandas as pd
import scipy.sparse as sp
from spectrum_fundamentals.metrics.percolator import Percolator

from ..data import FragmentType
from ..utils import JobPool

if TYPE_CHECKING:
    from ..data import Spectra

logger = logging.getLogger(__name__)


def generate_features(
    library: Spectra,
    search_type: str,
    output_file: str | Path,
    additional_columns: str | list | None = None,
    all_features: bool = False,
    xl: bool = False,
    cms2: bool = False,
    regression_method: str = "spline",
    add_neutral_loss_features: bool = False,
    remove_miss_cleavage_features: bool = False,
    task: str = "default",
    featured_ions: Optional[list] = None,
):
    """
    Generate features to be used for percolator or mokapot target decoy separation.

    The function calculates a range of metrics and features on the provided library for the chosen
    fdr estimation method, then writes the input tab file to the chosen output file.

    :param library: the library to perform feature generation on
    :param search_type: One of "original" and "rescore", which determines the generated features
    :param output_file: the location to the generated tab file to be used for percolator / mokapot
    :param additional_columns: additional columns supplied in the search results to be used as features (either a list or "all")
    :param all_features: whether to use all features or only the standard set TODO
    :param xl: crosslinked or linear peptide
    :param cms2: cleavable or non-cleavable crosslinker
    :param regression_method: The regression method to use for iRT alignment
    :param add_neutral_loss_features: Flag to indicate whether to add neutral loss features to percolator or not
    :param remove_miss_cleavage_features: Flag to indicate whether to remove miss cleavage features from percolator or not
    :param task: Flag to indicate whether to use multifrag features or not
    :param featured_ions: The ion series to use for calculating percolator features

    :Example:

    .. code-block:: python

        >>> from oktoberfest import rescore as re
        >>> from oktoberfest import predict as pr
        >>> from oktoberfest.data import Spectra, FragmentType
        >>> import pandas as pd
        >>> import numpy as np
        >>> # Required columns: RAW_FILE, MODIFIED_SEQUENCE, SEQUENCE, CALCULATED_MASS, SCAN_NUMBER,
        >>> # COLLISION_ENERGY, PRECURSOR_CHARGE, REVERSE and SCORE
        >>> meta_df = pd.DataFrame({"RAW_FILE": ["File1","File1"],
        >>>                         "MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
        >>>                         "SEQUENCE": ["AAACRFVQ","RMPCHKPYL"],
        >>>                         "CALCULATED_MASS": [1000,4000],
        >>>                         "SCAN_NUMBER": [1,2],
        >>>                         "COLLISION_ENERGY": [30,35],
        >>>                         "PRECURSOR_CHARGE": [1,2],
        >>>                         "FRAGMENTATION": ["HCD","HCD"],
        >>>                         "REVERSE": [False,False],
        >>>                         "SCORE": [0,0]})
        >>> var = Spectra._gen_vars_df()
        >>> library = Spectra(obs=meta_df, var=var)
        >>> raw_intensities = np.random.rand(2,174)
        >>> mzs = np.random.rand(2,174)*1000
        >>> annotation = np.array([var.index,var.index])
        >>> library.add_intensities(raw_intensities, annotation, FragmentType.RAW)
        >>> library.add_mzs(mzs, FragmentType.MZ)
        >>> library.strings_to_categoricals()
        >>> intensity_predictor = pr.Predictor.from_koina(
        >>>                         model_name="Prosit_2020_intensity_HCD",
        >>>                         server_url="koina.wilhelmlab.org:443",
        >>>                         ssl=True,
        >>> intensity_predictor.predict_intensities(data=library)
        >>> re.generate_features(library=library,
        >>>                         search_type="original",
        >>>                         regression_method="spline",
        >>>                         output_file="./tests/doctests/output/original.tab")
    """
    meta = library.get_meta_data().reset_index(drop=True)

    if search_type == "rescore" and xl:
        pred = sp.hstack([library.get_matrix(FragmentType.PRED_A), library.get_matrix(FragmentType.PRED_B)])
        raw = sp.hstack([library.get_matrix(FragmentType.RAW_A), library.get_matrix(FragmentType.RAW_B)])
        mz = sp.hstack([library.get_matrix(FragmentType.MZ_A), library.get_matrix(FragmentType.MZ_B)])
    elif search_type == "rescore":
        pred = library.get_matrix(FragmentType.PRED)
        raw = library.get_matrix(FragmentType.RAW)
        mz = library.get_matrix(FragmentType.MZ)
    else:
        pred = None
        raw = None
        mz = None

    perc_features = Percolator(
        metadata=meta,
        pred_intensities=pred,
        true_intensities=raw,
        mz=mz,
        input_type=search_type,
        all_features_flag=all_features,
        regression_method=regression_method,
        neutral_loss_flag=add_neutral_loss_features,
        drop_miss_cleavage_flag=remove_miss_cleavage_features,
        cms2=cms2,
        additional_columns=additional_columns,
        task=task,
        featured_ions=featured_ions,
    )

    perc_features.calc()
    perc_features.write_to_file(str(output_file))


def _fillna_tab_file(tab_file: Path, include_header: bool) -> bytes:
    """Read a tab file, fill NaN values with 0, and return it serialized as tab-separated bytes."""
    df = pd.read_csv(tab_file, sep="\t")
    df = df.fillna(0)
    content = df.to_csv(sep="\t", index=False, header=include_header).encode()
    tab_file.unlink()
    return content


def merge_input(
    tab_files: list[Path],
    output_file: str | Path,
    num_threads: int = 1,
):
    r"""
    Merge spectra file identifier specific tab files into one large file for combined percolation.

    The function reads each tab file once, fills NaN values with 0 (optionally in parallel across files),
    and writes the concatenated result directly to the chosen output file location.

    Fastest solution according to:
    https://stackoverflow.com/questions/44211461/what-is-the-fastest-way-to-combine-100-csv-files-with-headers-into-one

    :param tab_files: list of paths pointing to the individual tab files to be concatenated
    :param output_file: path to the generated output tab file
    :param num_threads: number of threads used in parallel for filling NaN values per tab file

    :Example:

    .. code-block:: python

        >>> from oktoberfest import rescore as re
        >>> from pathlib import Path
        >>> import pandas as pd
        >>> rescore_df1 = pd.DataFrame({'SpecId': ["F1-81-AAAAAAALQAK-2-5","F1-15-VGVFQHGK-3-2"],
        >>>                           'Label': [1,0],
        >>>                           'ScanNr': [81,15],
        >>>                           'filename': ["F1","F1"],
        >>>                           'CID': [0,0],
        >>>                           'Charge1': [0,0],
        >>>                           'Charge2': [1,0],
        >>>                           'Charge3': [0,1],
        >>>                           'Charge4': [0,0],
        >>>                           'Charge5': [0,0],
        >>>                           'Charge6': [0,0],
        >>>                           'HCD': [1,1],
        >>>                           'KR': [1,1],
        >>>                           'Mass': [1402.18,1103.54],
        >>>                           'spectral_angle': [0.71,0.23],
        >>>                           'sequence_length': [11,8],
        >>>                           'Peptide': ["_.AAAAAAALQAK._","_.VGVFQHGK._"],
        >>>                           'Proteins': ["AAAAAAALQAK","VGVFQHGK"],
        >>>                           'RT': [64.79,57.84],
        >>>                           'iRT': [65.99,56.22],
        >>>                           'pred_RT': [58.86,55.34]})
        >>> rescore_df2 = pd.DataFrame({'SpecId': ["F2-13-AEAEQEKDQLR-1-11","F2-27-TGFLEQLK-2-7"],
        >>>                           'Label': [1,0],
        >>>                           'ScanNr': [13,27],
        >>>                           'filename': ["F2","F2"],
        >>>                           'CID': [0,0],
        >>>                           'Charge1': [1,0],
        >>>                           'Charge2': [0,1],
        >>>                           'Charge3': [0,0],
        >>>                           'Charge4': [0,0],
        >>>                           'Charge5': [0,0],
        >>>                           'Charge6': [0,0],
        >>>                           'HCD': [1,1],
        >>>                           'KR': [2,1],
        >>>                           'Mass': [1202.43,1009.14],
        >>>                           'spectral_angle': [0.55,0.12],
        >>>                           'sequence_length': [11,8],
        >>>                           'Peptide': ["_.AEAEQEKDQLR._","_.TGFLEQLK._"],
        >>>                           'Proteins': ["AEAEQEKDQLR","TGFLEQLK"],
        >>>                           'RT': [62.33,51.23],
        >>>                           'iRT': [63.98,53.24],
        >>>                           'pred_RT': [59.16,50.76]})
        >>> rescore_df1.to_csv("./tests/doctests/input/rescore1.tab",sep='\t',index=False)
        >>> rescore_df2.to_csv("./tests/doctests/input/rescore2.tab",sep='\t',index=False)
        >>> tabfile1 = Path("./tests/doctests/input/rescore1.tab")
        >>> tabfile2 = Path("./tests/doctests/input/rescore2.tab")
        >>> filelist = [tabfile1,tabfile2]
        >>> re.merge_input(tab_files=filelist, output_file="./tests/doctests/output/merged_rescore.tab")
    """
    if num_threads > 1:
        fillna_pool = JobPool(processes=num_threads)
        for i, tab_file in enumerate(tab_files):
            fillna_pool.apply_async(_fillna_tab_file, [tab_file, i == 0])
        contents = fillna_pool.check_pool()
    else:
        contents = [_fillna_tab_file(tab_file, i == 0) for i, tab_file in enumerate(tab_files)]

    with open(output_file, "wb") as fout:
        for content in contents:
            fout.write(content)


def rescore_with_percolator(
    input_file: str | Path,
    output_folder: str | Path | None = None,
    num_threads: int = 3,
    test_fdr: float = 0.01,
    train_fdr: float = 0.01,
    xl: bool = False,
):
    """
    Rescore using percolator.

    The function takes an input file location, as well as an output folder location for percolator outputs and
    executes percolator with additional optional parameters.

    :param input_file: Path to percolator tab file
    :param output_folder: An optional output folder for all percolator files, default is the parent directory of the input_file
    :param num_threads: The number of threads used in parallel for percolator
    :param test_fdr: the fdr cutoff for the test set
    :param train_fdr: the fdr cutoff for the train set
    :param xl: crosslinked or linear peptide
    :raises FileNotFoundError: if the input file does not exist
    """
    if isinstance(input_file, str):
        input_file = Path(input_file)

    if not input_file.is_file():
        raise FileNotFoundError(f"{input_file} does not exist.")

    if output_folder is None:
        output_folder = input_file.parent
    if isinstance(output_folder, str):
        output_folder = Path(output_folder)

    file_prefix = input_file.stem
    weights_file = output_folder / f"{file_prefix}.percolator.weights.csv"
    target_psms = output_folder / f"{file_prefix}.percolator.psms.txt"
    decoy_psms = output_folder / f"{file_prefix}.percolator.decoy.psms.txt"
    target_peptides = output_folder / f"{file_prefix}.percolator.peptides.txt"
    decoy_peptides = output_folder / f"{file_prefix}.percolator.decoy.peptides.txt"
    log_file = output_folder / f"{file_prefix}.log"
    if xl:
        cmd = f"percolator --weights {weights_file} \
                    --results-psms {target_psms} \
                    --decoy-results-psms {decoy_psms} \
                    --only-psms \
                    {input_file} 2> {log_file}"

    else:
        cmd = f"percolator --weights {weights_file} \
                        --num-threads {num_threads} \
                        --subset-max-train 500000 \
                        --post-processing-tdc \
                        --testFDR {test_fdr} \
                        --trainFDR {train_fdr} \
                        --results-psms {target_psms} \
                        --decoy-results-psms {decoy_psms} \
                        --results-peptides {target_peptides} \
                        --decoy-results-peptides {decoy_peptides} \
                        {input_file} 2> {log_file}"

    logger.info(f"Starting percolator with command {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    logger.info("Finished rescoring using percolator.")


def rescore_with_mokapot(
    input_file: str | Path,
    output_folder: str | Path | None = None,
    test_fdr: float = 0.01,
    xl: bool = False,
):
    """
    Rescore using mokapot.

    The function takes an input file location, as well as an output folder location for mokapot outputs and
    executes mokapot with additional optional parameters.

    :param input_file: Path to percolator tab file
    :param output_folder: An optional output folder for all percolator files, default is the parent directory of the input_file
    :param test_fdr: the fdr cutoff for the test set
    :param xl: crosslinked or linear peptide (currently unused)
    :raises FileNotFoundError: if the input file does not exist
    """
    if isinstance(input_file, str):
        input_file = Path(input_file)

    if not input_file.is_file():
        raise FileNotFoundError(f"{input_file} does not exist.")

    if output_folder is None:
        output_folder = input_file.parent
    if isinstance(output_folder, str):
        output_folder = Path(output_folder)

    mokapot_logger = logging.getLogger("mokapot")
    mokapot_logger.setLevel(logging.INFO)
    log_file = output_folder / f"{input_file.stem}.log"

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    log_formatter = logging.Formatter("%(levelname)s: %(message)s")
    file_handler.setFormatter(log_formatter)
    mokapot_logger.addHandler(file_handler)

    np.random.seed(123)

    psms = mokapot.read_pin(input_file)
    logger.info(f"Number of PSMs used for training: {len(psms)}")

    logger.info("Running mokapot...")
    results, models = mokapot.brew(psms, test_fdr=test_fdr)
    results.to_txt(dest_dir=output_folder, file_root=f"{input_file.stem}", decoys=True)
    logger.info("Finished rescoring using mokapot.")
