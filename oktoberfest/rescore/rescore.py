import logging
import subprocess
from pathlib import Path
from typing import List, Optional, Union

import mokapot
import numpy as np
import pandas as pd
from spectrum_fundamentals.metrics.percolator import Percolator

from ..data import FragmentType, Spectra

logger = logging.getLogger(__name__)


def generate_features(
    library: Spectra,
    search_type: str,
    output_file: Union[str, Path],
    all_features: bool = False,
    regression_method: str = "spline",
):
    """
    Generate features to be used for percolator or mokapot target decoy separation.

    The function calculates a range of metrics and features on the provided library for the chosen
    fdr estimation method, then writes the input tab file to the chosen output file.

    :param library: the library to perform feature generation on
    :param search_type: One of "original" and "rescore", which determines the generated features
    :param output_file: the location to the generated tab file to be used for percolator / mokapot
    :param all_features: whether to use all features or only the standard set TODO
    :param regression_method: The regression method to use for iRT alignment
    """
    perc_features = Percolator(
        metadata=library.get_meta_data(),
        pred_intensities=library.get_matrix(FragmentType.PRED)[0],
        true_intensities=library.get_matrix(FragmentType.RAW)[0],
        mz=library.get_matrix(FragmentType.MZ)[0],
        input_type=search_type,
        all_features_flag=all_features,
        regression_method=regression_method,
    )
    perc_features.calc()
    perc_features.write_to_file(str(output_file))


def merge_input(
    tab_files: List[Path],
    output_file=Union[str, Path],
):
    """
    Merge spectra file identifier specific tab files into one large file for combined percolation.

    The function takes a list of tab files and concatenates them before writing a combined tab file back to
    the chosen output file location.

    Fastest solution according to:
    https://stackoverflow.com/questions/44211461/what-is-the-fastest-way-to-combine-100-csv-files-with-headers-into-one

    :param tab_files: list of paths pointing to the individual tab files to be concatenated
    :param output_file: path to the generated output tab file
    """
    with open(output_file, "wb") as fout:
        first = True
        for tab_file in tab_files:
            with open(tab_file, "rb") as f:
                if not first:
                    next(f)  # skip the header
                else:
                    first = False
                fout.write(f.read())
            tab_file.unlink()

    # TODO make this more efficient
    df_prosit = pd.read_csv(output_file, sep="\t")
    df_prosit = df_prosit.fillna(0)

    # We exploit the expmass column here by assigning a unique id per filename+scannr group.
    # This ensures percolator will not deduplicate scannrs between filenames, as it uses only
    # filename+ExpMass for TDC.
    df_prosit.insert(loc=4, column="ExpMass", value=df_prosit.groupby(["filename", "ScanNr"]).ngroup())

    df_prosit.to_csv(output_file, sep="\t", index=False)


def rescore_with_percolator(
    input_file: Union[str, Path],
    output_folder: Optional[Union[str, Path]] = None,
    num_threads: int = 3,
    test_fdr: float = 0.01,
    train_fdr: float = 0.01,
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
    input_file: Union[str, Path],
    output_folder: Optional[Union[str, Path]] = None,
    test_fdr: float = 0.01,
):
    """
    Rescore using mokapot.

    The function takes an input file location, as well as an output folder location for mokapot outputs and
    executes mokapot with additional optional parameters.

    :param input_file: Path to percolator tab file
    :param output_folder: An optional output folder for all percolator files, default is the parent directory of the input_file
    :param test_fdr: the fdr cutoff for the test set
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

    df = pd.read_csv(input_file, sep="\t")

    # TODO remove this if not necessary
    df = df.rename(columns={"Protein": "Proteins"})
    df.to_csv(input_file, sep="\t")

    psms = mokapot.read_pin(input_file)
    logger.info(f"Number of PSMs used for training: {len(psms)}")

    logger.info("Running mokapot...")
    results, models = mokapot.brew(psms, test_fdr=test_fdr)
    results.to_txt(dest_dir=output_folder, file_root=f"{input_file.stem}", decoys=True)
    logger.info("Finished rescoring using mokapot.")
