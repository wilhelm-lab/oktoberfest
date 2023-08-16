import logging
import multiprocessing
import subprocess
from pathlib import Path
from typing import List, Union

import mokapot
import numpy as np
import pandas as pd

from .calculate_features import CalculateFeatures
from .utils.multiprocessing_pool import JobPool
from .utils.process_step import ProcessStep

logger = logging.getLogger(__name__)


# This function cannot be a function inside ReScore since the multiprocessing pool does not work with class member functions
def calculate_features_single(
    raw_file_path: Path,
    split_msms_path: Path,
    fdr_estimation_method_rescore_path: Path,
    fdr_estimation_method_orig_path: Path,
    mzml_path: Path,
    config_path: Path,
    calc_feature_step: ProcessStep,
):
    """Create CalculateFeatures object and calculate features for a given raw file."""
    logger.info(f"Calculating features for {raw_file_path}")
    features = CalculateFeatures(search_path="", raw_path=raw_file_path, out_path=mzml_path, config_path=config_path)

    df_search = pd.read_csv(split_msms_path, delimiter="\t")
    features.predict_with_aligned_ce(df_search)
    features.gen_perc_metrics("rescore", fdr_estimation_method_rescore_path)
    features.gen_perc_metrics("original", fdr_estimation_method_orig_path)

    calc_feature_step.mark_done()


def get_glob_pattern(spectra_type: str) -> str:
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


class ReScore(CalculateFeatures):
    """
    Main to init a re-score obj and go through the steps.

    1- get_raw_files
    2- split_msms
    3- calculate_features
    4- merge_input
    5- rescore_with_percolator or mokaput
    """

    raw_files: List[Path]

    split_msms_step: ProcessStep
    merge_input_step_prosit: ProcessStep
    merge_input_step_andromeda: ProcessStep
    rescore_step_prosit: ProcessStep
    rescore_step_andromeda: ProcessStep

    def __init__(
        self,
        search_path: Union[str, Path],
        raw_path: Union[str, Path],
        out_path: Union[str, Path],
        config_path: Union[str, Path],
        mzml_reader_package: str = "pymzml",
    ):
        """
        Initialize a ReScore object and go through the steps.

        1- get_raw_files
        2- split_msms
        3- calculate_features
        4- merge_input
        5- rescore_with_percolator or mokapot

        :param search_path: path to search directory
        :param raw_path: path to raw file as a string
        :param out_path: path to output folder
        :param config_path: path to config file
        :param mzml_reader_package: mzml reader (pymzml or pyteomics)
        """
        super().__init__(
            search_path, raw_path, out_path, config_path=config_path, mzml_reader_package=mzml_reader_package
        )
        self.split_msms_step = ProcessStep(out_path, "split_msms")
        self.merge_input_step_prosit = ProcessStep(out_path, "merge_input_prosit")
        self.merge_input_step_andromeda = ProcessStep(out_path, "merge_input_andromeda")
        self.rescore_step_prosit = ProcessStep(out_path, f"{self.config.fdr_estimation_method}_prosit")
        self.rescore_step_andromeda = ProcessStep(out_path, f"{self.config.fdr_estimation_method}_andromeda")

    def get_raw_files(self):
        """
        Obtains raw files by scanning through the raw_path directory.

        If raw_path is a file, only process this one.
        :raises FileNotFoundError: if raw file could not be found
        """
        self.raw_files = []
        if self.raw_path.is_file():
            self.raw_files = [self.raw_path]
            self.raw_path = self.raw_path.parent
        elif self.raw_path.is_dir():
            spectra_type = self.config.spectra_type
            glob_pattern = get_glob_pattern(spectra_type)
            self.raw_files = list(self.raw_path.glob(glob_pattern))
            logger.info(f"Found {len(self.raw_files)} raw files in the search directory")
        else:
            raise FileNotFoundError(f"{self.raw_path} does not exist.")

    def split_msms(self):
        """Splits search results per raw file identifier such that we can process each raw file in parallel \
        without reading the entire search results multiple times."""
        out_path = self.out_path
        out_path.mkdir(exist_ok=True)

        if self.split_msms_step.is_done():
            return
        self.get_msms_folder_path().mkdir(exist_ok=True)

        df_search = self._load_search()
        logger.info(f"Read {len(df_search.index)} PSMs from {self.search_path}")
        for raw_file, df_search_split in df_search.groupby("RAW_FILE"):
            spectra_file_path = self.raw_path / raw_file
            spectra_type = self.config.spectra_type
            glob_pattern = get_glob_pattern(spectra_type)
            matched_files = [file for file in spectra_file_path.parent.glob(glob_pattern)]
            if not any(path.is_file() for path in matched_files):
                logger.info(f"Did not find {raw_file} in search directory, skipping this file")
                continue
            split_msms = self._get_split_msms_path(raw_file + ".rescore")
            logger.info(f"Creating split search results file {split_msms}")
            df_search_split = df_search_split[(df_search_split["PEPTIDE_LENGTH"] <= 30)]
            df_search_split = df_search_split[(~df_search_split["MODIFIED_SEQUENCE"].str.contains(r"\(ac\)"))]
            df_search_split = df_search_split[
                (~df_search_split["MODIFIED_SEQUENCE"].str.contains(r"\(Acetyl \(Protein N-term\)\)"))
            ]
            df_search_split = df_search_split[(~df_search_split["SEQUENCE"].str.contains("U"))]
            df_search_split = df_search_split[df_search_split["PRECURSOR_CHARGE"] <= 6]
            df_search_split = df_search_split[df_search_split["PEPTIDE_LENGTH"] >= 7]
            df_search_split.to_csv(split_msms, sep="\t", index=False)

        self.split_msms_step.mark_done()

    def calculate_features(self):
        """Calculates percolator/mokapot input features per raw file using multiprocessing."""
        num_threads = self.config.num_threads
        self.config
        if num_threads > 1:
            processing_pool = JobPool(processes=num_threads)

        mzml_path = self.get_mzml_folder_path()
        mzml_path.mkdir(exist_ok=True)

        data_path = self.out_path / "data"
        data_path.mkdir(exist_ok=True)

        fdr_estimation_method_path = self.get_fdr_estimation_method_folder_path()
        fdr_estimation_method_path.mkdir(exist_ok=True)

        for raw_file in self.raw_files:
            calc_feature_step = ProcessStep(self.out_path, "calculate_features." + raw_file.stem)
            if calc_feature_step.is_done():
                continue

            fdr_estimation_method_rescore_path = self._get_split_perc_input_path(raw_file.stem, "rescore")
            fdr_estimation_method_orig_path = self._get_split_perc_input_path(raw_file.stem, "original")

            split_msms_path = self._get_split_msms_path(raw_file.with_suffix(".rescore").name)

            args = [
                raw_file,
                split_msms_path,
                fdr_estimation_method_rescore_path,
                fdr_estimation_method_orig_path,
                self.out_path,
                self.config_path,
                calc_feature_step,
            ]
            if num_threads > 1:
                processing_pool.apply_async(calculate_features_single, args)
            else:
                calculate_features_single(*args)

        if num_threads > 1:
            processing_pool.check_pool()

    def merge_input(self, search_type: str = "rescore"):
        """
        Merge percolator/mokapot input files into one large file for combined percolation.

        Fastest solution according to:
        https://stackoverflow.com/questions/44211461/what-is-the-fastest-way-to-combine-100-csv-files-with-headers-into-one

        :param search_type: choose either rescore or original to merge percolator/mokapot files for this.
        """
        if search_type == "rescore":
            if self.merge_input_step_prosit.is_done():
                return
        else:
            if self.merge_input_step_andromeda.is_done():
                return

        merged_perc_input_file_prosit = self._get_merged_perc_input_path(search_type)

        logger.info(f"Merging {self.config.fdr_estimation_method} input files for {search_type}")
        with open(merged_perc_input_file_prosit, "wb") as fout:
            first = True
            for raw_file in self.raw_files:
                fdr_estimation_method_input_path = self._get_split_perc_input_path(raw_file.stem, search_type)
                with open(fdr_estimation_method_input_path, "rb") as f:
                    if not first:
                        next(f)  # skip the header
                    else:
                        first = False
                    fout.write(f.read())

        df_prosit = pd.read_csv(merged_perc_input_file_prosit, sep="\t")
        df_prosit = df_prosit.fillna(0)
        df_prosit.to_csv(merged_perc_input_file_prosit, sep="\t", index=False)

        if search_type == "rescore":
            self.merge_input_step_prosit.mark_done()
        else:
            self.merge_input_step_andromeda.mark_done()

    def rescore(self, search_type: str = "rescore", test_fdr: float = 0.01, train_fdr: float = 0.01):
        """Use percolator/mokapot to re-score library."""
        if search_type == "rescore" and self.rescore_step_prosit.is_done():
            return
        elif self.rescore_step_andromeda.is_done():
            return

        fdr_estimation_method_path = self.get_fdr_estimation_method_folder_path()
        weights_file = fdr_estimation_method_path / f"{search_type}.percolator.weights.csv"
        target_psms = fdr_estimation_method_path / f"{search_type}.percolator.psms.txt"
        decoy_psms = fdr_estimation_method_path / f"{search_type}.percolator.decoy.psms.txt"
        target_peptides = fdr_estimation_method_path / f"{search_type}.percolator.peptides.txt"
        decoy_peptides = fdr_estimation_method_path / f"{search_type}.percolator.decoy.peptides.txt"
        log_file = fdr_estimation_method_path / f"{search_type}.log"

        fdr_estimation_method = self.config.fdr_estimation_method
        if fdr_estimation_method == "percolator":
            cmd = f"percolator --weights {weights_file} \
                            --num-threads {self.config.num_threads} \
                            --subset-max-train 500000 \
                            --post-processing-tdc \
                            --testFDR {test_fdr} \
                            --trainFDR {train_fdr} \
                            --results-psms {target_psms} \
                            --decoy-results-psms {decoy_psms} \
                            --results-peptides {target_peptides} \
                            --decoy-results-peptides {decoy_peptides} \
                            {self._get_merged_perc_input_path(search_type)} 2> {log_file}"
            logger.info(f"Starting percolator with command {cmd}")
            subprocess.run(cmd, shell=True, check=True)
        elif fdr_estimation_method == "mokapot":
            logger.info("Starting mokapot rescoring")
            mokapot_logger = logging.getLogger("mokapot")
            mokapot_logger.setLevel(logging.INFO)
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(logging.INFO)
            log_formatter = logging.Formatter("%(levelname)s: %(message)s")
            file_handler.setFormatter(log_formatter)
            mokapot_logger.addHandler(file_handler)
            np.random.seed(123)
            file_path = fdr_estimation_method_path / f"{search_type}.tab"
            df = pd.read_csv(file_path, sep="\t")
            df = df.rename(columns={"Protein": "Proteins"})
            df.to_csv(file_path, sep="\t")
            psms = mokapot.read_pin(file_path)
            num_cores = multiprocessing.cpu_count()
            logger.info(f"Running mokapot on {num_cores} cores")
            results, models = mokapot.brew(psms, test_fdr=test_fdr)
            results.to_txt(dest_dir=fdr_estimation_method_path, file_root=f"{search_type}", decoys=True)
            logger.info(f"Number of PSMs used for training: {len(psms)}")

        else:
            raise ValueError(
                f"Unknown fdr estimation method: {fdr_estimation_method}. Choose between mokapot and percolator."
            )

        if search_type == "rescore":
            self.rescore_step_prosit.mark_done()
        else:
            self.rescore_step_andromeda.mark_done()
        logger.info(f"Finished rescoring using {fdr_estimation_method}.")

    def _get_split_msms_path(self, raw_file: str) -> Path:
        """
        Get path to split search results file.

        :param raw_file: path to raw file as a string
        :return: path to split search results file
        """
        return self.get_msms_folder_path() / raw_file

    def get_mzml_folder_path(self) -> Path:
        """Get folder path to mzml."""
        if self.config.spectra_type == "mzml":
            return self.config.spectra
        else:
            return self.out_path / "mzML"

    def get_fdr_estimation_method_folder_path(self) -> Path:
        """Get folder path to percolator/mokapot."""
        return self.results_path / f"{self.config.fdr_estimation_method}"

    def _get_split_perc_input_path(self, raw_file: str, search_type: str) -> Path:
        """
        Specify search_type to differentiate between percolator/mokapot and andromeda output.

        :param raw_file: path to raw file as a string
        :param search_type: model (rescore or original) as a string
        :return: path to split percolator/mokapot input file
        """
        return self.get_fdr_estimation_method_folder_path() / f"{raw_file}.{search_type}.tab"

    def _get_merged_perc_input_path(self, search_type: str) -> Path:
        """
        Get merged percolator/mokapot input path.

        :param search_type: model (rescore or original) as a string
        :return: path to merged percolator/mokapot input folder
        """
        return self.get_fdr_estimation_method_folder_path() / f"{search_type}.tab"
