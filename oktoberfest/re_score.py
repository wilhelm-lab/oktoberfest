import logging
import os
import subprocess
import tracemalloc
from typing import List

import numpy as np
import pandas as pd

from .calculate_features import CalculateFeatures
from .utils.process_step import ProcessStep

logger = logging.getLogger(__name__)


# This function cannot be a function inside ReScore since the multiprocessing pool does not work with class member functions
def calculate_features_single(
    raw_file_path, split_msms_path, percolator_input_path, mzml_path, config_path, calc_feature_step
):
    logger.info(f"Calculating features for {raw_file_path}")
    print(raw_file_path, split_msms_path, percolator_input_path, mzml_path)
    features = CalculateFeatures(search_path="", raw_path=raw_file_path, out_path=mzml_path, config_path=config_path)

    df_search = pd.read_csv(split_msms_path, delimiter="\t")
    features.predict_with_aligned_ce(df_search)
    features.gen_perc_metrics("prosit", percolator_input_path)
    features.gen_perc_metrics("andromeda", percolator_input_path.replace("prosit", "andromeda"))

    calc_feature_step.mark_done()


class ReScore(CalculateFeatures):
    """
    main to init a re-score obj and go through the steps:
    1- get_raw_files
    2- split_msms
    3- calculate_features
    4- merge_input
    5- rescore_with_perc
    """

    raw_files: List[str]

    split_msms_step: ProcessStep
    merge_input_step_prosit: ProcessStep
    merge_input_step_andromeda: ProcessStep
    percolator_step_prosit: ProcessStep
    percolator_step_andromeda: ProcessStep

    def __init__(self, search_path, raw_path, out_path, config_path=None, mzml_reader_package="pymzml"):
        super().__init__(
            search_path, raw_path, out_path, config_path=config_path, mzml_reader_package=mzml_reader_package
        )
        self.split_msms_step = ProcessStep(out_path, "split_msms")
        self.merge_input_step_prosit = ProcessStep(out_path, "merge_input_prosit")
        self.merge_input_step_andromeda = ProcessStep(out_path, "merge_input_andromeda")
        self.percolator_step_prosit = ProcessStep(out_path, "percolator_prosit")
        self.percolator_step_andromeda = ProcessStep(out_path, "percolator_andromeda")

    def get_raw_files(self):
        """
        Obtains raw files by scanning through the raw_path directory. If raw_path is a file, only process this one.
        """
        self.raw_files = []
        if os.path.isfile(self.raw_path):
            self.raw_files = [self.raw_path]
            self.raw_path = os.path.dirname(self.raw_path)
        elif os.path.isdir(self.raw_path):
            raw_type = self.config.get_raw_type()
            if raw_type == "thermo":
                extension = ".raw"
            elif raw_type == "mzml":
                extension = ".mzml"
            else:
                raise ValueError(f"{raw_type} is not supported as rawfile-type")

            self.raw_files = [os.path.basename(f) for f in os.listdir(self.raw_path) if f.lower().endswith(extension)]
            logger.info(f"Found {len(self.raw_files)} raw files in the search directory")

    def split_msms(self):
        """
        Splits msms.txt file per raw file such that we can process each raw file in parallel without reading the entire msms.txt
        """
        if self.split_msms_step.is_done():
            return

        msms_path = self.get_msms_folder_path()
        if not os.path.isdir(msms_path):
            os.makedirs(msms_path)

        df_search = self._load_search()
        logger.info(f"Read {len(df_search.index)} PSMs from {self.search_path}")
        for raw_file, df_search_split in df_search.groupby("RAW_FILE"):

            if not os.path.isfile(os.path.join(self.raw_path, raw_file + ".raw")):
                logger.info(f"Did not find {raw_file} in search directory, skipping this file")
                continue

            split_msms = self._get_split_msms_path(raw_file)
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

        self.split_msms_step.mark_done()

    def calculate_features(self):
        """
        Calculates percolator input features per raw file using multiprocessing
        """
        num_threads = self.config.get_num_threads()
        if num_threads > 1:
            from .utils.multiprocessing_pool import JobPool

            processingPool = JobPool(processes=num_threads)

        mzml_path = self.get_mzml_folder_path()
        if not os.path.isdir(mzml_path):
            os.makedirs(mzml_path)

        perc_path = self.get_percolator_folder_path()
        if not os.path.isdir(perc_path):
            os.makedirs(perc_path)

        for raw_file in self.raw_files:
            calc_feature_step = ProcessStep(self.out_path, "calculate_features." + raw_file)
            if calc_feature_step.is_done():
                continue

            raw_file_path = os.path.join(self.raw_path, raw_file)

            mzml_file_path = os.path.join(mzml_path, raw_file.replace(".raw", ".mzML"))

            percolator_input_path = self._get_split_perc_input_path(raw_file, "prosit")
            split_msms_path = self._get_split_msms_path(raw_file)

            if num_threads > 1:
                processingPool.applyAsync(
                    calculate_features_single,
                    (
                        raw_file_path,
                        split_msms_path,
                        percolator_input_path,
                        mzml_file_path,
                        self.config_path,
                        calc_feature_step,
                    ),
                )
            else:
                calculate_features_single(
                    raw_file_path,
                    split_msms_path,
                    percolator_input_path,
                    mzml_file_path,
                    self.config_path,
                    calc_feature_step,
                )

        if num_threads > 1:
            processingPool.checkPool(printProgressEvery=1)

    def merge_input(self, search_type: str = "prosit"):
        """
        Merges percolator input files into one large file for combined percolation.
        Fastest solution according to: https://stackoverflow.com/questions/44211461/what-is-the-fastest-way-to-combine-100-csv-files-with-headers-into-one
        search_type: choose either prosit or andromeda to merge percolator files for this.
        """
        if search_type == "prosit":
            if self.merge_input_step_prosit.is_done():
                return
        else:
            if self.merge_input_step_andromeda.is_done():
                return

        merged_perc_input_file_prosit = self._get_merged_perc_input_path(search_type)

        logger.info(f"Merging percolator input files for " + search_type)
        with open(merged_perc_input_file_prosit, "wb") as fout:
            first = True
            for raw_file in self.raw_files:
                percolator_input_path = self._get_split_perc_input_path(raw_file, search_type)
                with open(percolator_input_path, "rb") as f:
                    if not first:
                        next(f)  # skip the header
                    else:
                        first = False
                    fout.write(f.read())
                os.remove(percolator_input_path)

        if search_type == "prosit":
            self.merge_input_step_prosit.mark_done()
        else:
            self.merge_input_step_andromeda.mark_done()

    def rescore_with_perc(self, search_type: str = "prosit", test_fdr: float = 0.01, train_fdr: float = 0.01):
        """
        Use percolator to re-score library.
        """
        if search_type == "prosit":
            if self.percolator_step_prosit.is_done():
                return
        else:
            if self.percolator_step_andromeda.is_done():
                return

        perc_path = self.get_percolator_folder_path()
        weights_file = os.path.join(perc_path, f"{search_type}_weights.csv")
        target_psms = os.path.join(perc_path, f"{search_type}_target.psms")
        decoy_psms = os.path.join(perc_path, f"{search_type}_decoy.psms")
        target_peptides = os.path.join(perc_path, f"{search_type}_target.peptides")
        decoy_peptides = os.path.join(perc_path, f"{search_type}_decoy.peptides")
        log_file = os.path.join(perc_path, f"{search_type}.log")

        cmd = f"percolator --weights {weights_file} \
                          --num-threads {self.config.get_num_threads()} \
                          --subset-max-train 500000 \
                          --post-processing-tdc \
                          --search-input concatenated \
                          --testFDR {test_fdr} \
                          --trainFDR {train_fdr} \
                          --results-psms {target_psms} \
                          --decoy-results-psms {decoy_psms} \
                          --results-peptides {target_peptides} \
                          --decoy-results-peptides {decoy_peptides} \
                          {self._get_merged_perc_input_path(search_type)} 2> {log_file}"
        logger.info(f"Starting percolator with command {cmd}")
        subprocess.run(cmd, shell=True, check=True)

        if search_type == "prosit":
            self.percolator_step_prosit.mark_done()
        else:
            self.percolator_step_andromeda.mark_done()

    def get_msms_folder_path(self):
        return os.path.join(self.out_path, "msms")

    def _get_split_msms_path(self, raw_file: str):
        return os.path.join(self.get_msms_folder_path(), os.path.splitext(raw_file)[0] + ".prosit")

    def get_mzml_folder_path(self):
        return os.path.join(self.out_path, "mzML")

    def get_percolator_folder_path(self):
        return os.path.join(self.results_path, "percolator")

    def _get_split_perc_input_path(self, raw_file: str, search_type: str):
        """
        Specify search_type to differentiate between percolator and andromeda output.
        """
        return os.path.join(
            self.get_percolator_folder_path(), os.path.splitext(raw_file)[0] + "_" + search_type + ".tab"
        )

    def _get_merged_perc_input_path(self, search_type: str):
        return os.path.join(self.get_percolator_folder_path(), search_type + ".tab")
