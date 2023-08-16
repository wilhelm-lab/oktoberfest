import logging
import subprocess
from pathlib import Path
from typing import List

import mokapot
import numpy as np
import pandas as pd

from ..calculate_features import CalculateFeatures
from ..utils.config import Config
from ..utils.multiprocessing_pool import JobPool
from ..utils.process_step import ProcessStep

logger = logging.getLogger(__name__)


# This function cannot be a function inside ReScore since the multiprocessing pool does not work with class member functions
def calculate_features_single(
    raw_file_path: Path,
    split_msms_path: Path,
    percolator_rescore_path: Path,
    percolator_orig_path: Path,
    mzml_path: Path,
    config_path: Path,
    calc_feature_step: ProcessStep,
):
    """Create CalculateFeatures object and calculate features for a given raw file."""
    logger.info(f"Calculating features for {raw_file_path}")
    features = CalculateFeatures(search_path="", raw_path=raw_file_path, out_path=mzml_path, config_path=config_path)

    df_search = pd.read_csv(split_msms_path, delimiter="\t")
    features.predict_with_aligned_ce(df_search)
    features.gen_perc_metrics("rescore", percolator_rescore_path)
    features.gen_perc_metrics("original", percolator_orig_path)

    calc_feature_step.mark_done()


def calculate_features(config: Config, config_path: Path, raw_files: List[Path]):
    """Calculates percolator input features per raw file using multiprocessing."""
    num_threads = config.num_threads
    if num_threads > 1:
        processing_pool = JobPool(processes=num_threads)

    mzml_path = get_mzml_folder_path(config)
    mzml_path.mkdir(exist_ok=True)

    data_path = config.output / "data"
    data_path.mkdir(exist_ok=True)

    results_path = config.output / "results"
    results_path.mkdir(exist_ok=True)

    perc_path = get_percolator_folder_path(config)
    perc_path.mkdir(exist_ok=True)

    for raw_file in raw_files:
        calc_feature_step = ProcessStep(config.output, "calculate_features." + raw_file.stem)
        if calc_feature_step.is_done():
            continue

        percolator_rescore_path = _get_split_perc_input_path(config, raw_file.stem, "rescore")
        percolator_orig_path = _get_split_perc_input_path(config, raw_file.stem, "original")

        split_msms_path = _get_split_msms_path(config, raw_file.with_suffix(".rescore").name)

        args = [
            raw_file,
            split_msms_path,
            percolator_rescore_path,
            percolator_orig_path,
            config.output,
            config_path,
            calc_feature_step,
        ]
        if num_threads > 1:
            processing_pool.apply_async(calculate_features_single, args)
        else:
            calculate_features_single(*args)

    if num_threads > 1:
        processing_pool.check_pool(print_progress_every=1)


def get_mzml_folder_path(config: Config) -> Path:
    """Get folder path to mzml."""
    if config.spectra_type == "mzml":
        return config.spectra
    else:
        return config.output / "mzML"


def get_percolator_folder_path(config: Config) -> Path:
    """Get folder path to percolator."""
    return config.output / "results" / "percolator"


def _get_split_perc_input_path(config: Config, raw_file: str, search_type: str) -> Path:
    """
    Specify search_type to differentiate between percolator and andromeda output.

    :param raw_file: path to raw file as a string
    :param search_type: model (rescore or original) as a string
    :return: path to split percolator input file
    """
    return get_percolator_folder_path(config) / f"{raw_file}.{search_type}.tab"


def _get_split_msms_path(config: Config, raw_file: str) -> Path:
    """
    Get path to split msms.

    :param raw_file: path to raw file as a string
    :return: path to split msms file
    """
    return get_msms_folder_path(config) / raw_file


def get_msms_folder_path(config: Config) -> Path:
    """Get folder path to msms."""
    return config.output / "msms"


def _get_merged_perc_input_path(config: Config, search_type: str) -> Path:
    """
    Get merged percolator input path.

    :param search_type: model (rescore or original) as a string
    :return: path to merged percolator input folder
    """
    return get_percolator_folder_path(config) / f"{search_type}.tab"


def merge_input(
    config: Config,
    merge_input_step_prosit: ProcessStep,
    merge_input_step_andromeda: ProcessStep,
    raw_files: List[Path],
    search_type: str = "rescore",
):
    """
    Merge percolator input files into one large file for combined percolation.

    Fastest solution according to:
    https://stackoverflow.com/questions/44211461/what-is-the-fastest-way-to-combine-100-csv-files-with-headers-into-one

    :param search_type: choose either rescore or original to merge percolator files for this.
    """
    if search_type == "rescore":
        if merge_input_step_prosit.is_done():
            return
    else:
        if merge_input_step_andromeda.is_done():
            return

    merged_perc_input_file_prosit = _get_merged_perc_input_path(config, search_type)

    logger.info(f"Merging percolator input files for {search_type}")
    with open(merged_perc_input_file_prosit, "wb") as fout:
        first = True
        for raw_file in raw_files:
            percolator_input_path = _get_split_perc_input_path(config, raw_file.stem, search_type)
            with open(percolator_input_path, "rb") as f:
                if not first:
                    next(f)  # skip the header
                else:
                    first = False
                fout.write(f.read())
            percolator_input_path.unlink()

    df_prosit = pd.read_csv(merged_perc_input_file_prosit, sep="\t")
    df_prosit = df_prosit.fillna(0)
    df_prosit.to_csv(merged_perc_input_file_prosit, sep="\t", index=False)

    if search_type == "rescore":
        merge_input_step_prosit.mark_done()
    else:
        merge_input_step_andromeda.mark_done()
    return merge_input_step_prosit, merge_input_step_andromeda


def rescore(
    config: Config,
    rescore_step_prosit: ProcessStep,
    rescore_step_andromeda: ProcessStep,
    search_type: str = "rescore",
    test_fdr: float = 0.01,
    train_fdr: float = 0.01,
):
    """Use percolator to re-score library."""
    if search_type == "rescore" and rescore_step_prosit.is_done():
        return
    elif rescore_step_andromeda.is_done():
        return

    perc_path = get_percolator_folder_path(config)
    weights_file = perc_path / f"{search_type}.weights.csv"
    target_psms = perc_path / f"{search_type}.target.psms"
    decoy_psms = perc_path / f"{search_type}.decoy.psms"
    target_peptides = perc_path / f"{search_type}.target.peptides"
    decoy_peptides = perc_path / f"{search_type}.decoy.peptides"
    log_file = perc_path / f"{search_type}.log"

    fdr_estimation_method = config.fdr_estimation_method
    if fdr_estimation_method == "percolator":
        cmd = f"percolator --weights {weights_file} \
                        --num-threads {config.num_threads} \
                        --subset-max-train 500000 \
                        --post-processing-tdc \
                        --testFDR {test_fdr} \
                        --trainFDR {train_fdr} \
                        --results-psms {target_psms} \
                        --decoy-results-psms {decoy_psms} \
                        --results-peptides {target_peptides} \
                        --decoy-results-peptides {decoy_peptides} \
                        {_get_merged_perc_input_path(config, search_type)} 2> {log_file}"
        logger.info(f"Starting percolator with command {cmd}")
        subprocess.run(cmd, shell=True, check=True)
    elif fdr_estimation_method == "mokapot":
        logger.info("Starting mokapot rescoring")
        np.random.seed(123)
        file_path = perc_path / f"{search_type}.tab"
        df = pd.read_csv(file_path, sep="\t")
        df = df.rename(columns={"Protein": "Proteins"})
        df.to_csv(file_path, sep="\t")
        psms = mokapot.read_pin(file_path)
        results, models = mokapot.brew(psms, test_fdr=test_fdr)
        results.to_txt(dest_dir=perc_path, file_root=f"{search_type}", decoys=True)
    else:
        raise ValueError(
            f"Unknown fdr estimation method: {fdr_estimation_method}. Choose between mokapot and percolator."
        )

    if search_type == "rescore":
        rescore_step_prosit.mark_done()
    else:
        rescore_step_andromeda.mark_done()
    logger.info(f"Finished rescoring using {fdr_estimation_method}.")
    return rescore_step_prosit, rescore_step_andromeda
