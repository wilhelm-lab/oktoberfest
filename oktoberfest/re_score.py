import os
import logging
from typing import List
import subprocess

import numpy as np
import pandas as pd

from .calculate_features import CalculateFeatures
from .utils.process_step import ProcessStep

logger = logging.getLogger(__name__)


# This function cannot be a function inside ReScore since the multiprocessing pool does not work with class member functions
def calculate_features_single(raw_file_path, split_msms_path, percolator_input_path, calc_feature_step):
    logger.info(f"Calculating features for {raw_file_path}")
    features = CalculateFeatures(search_path = "",
                      raw_path = raw_file_path)
    df_search = pd.read_csv(split_msms_path, delimiter = '\t')
    features.predict_with_aligned_ce(df_search)
    features.gen_perc_metrics(percolator_input_path)
    
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
    merge_input_step: ProcessStep
    percolator_step: ProcessStep
    
    def __init__(self, search_path, raw_path, config_path=None):
        super().__init__(search_path, raw_path, config_path=config_path)
        self.split_msms_step = ProcessStep(raw_path, "split_msms")
        self.merge_input_step = ProcessStep(raw_path, "merge_input")
        self.percolator_step = ProcessStep(raw_path, "percolator")

    def get_raw_files(self):
        self.raw_files = []
        if os.path.isfile(self.raw_path):
            self.raw_files = [self.raw_path]
            self.raw_path = os.path.dirname(self.raw_path)
        elif os.path.isdir(self.raw_path):
            switch = self.config["fileUploads"]["raw_type"]
            if switch == "thermo":
                extension = ".raw"
            elif switch == "mzml":
                extension = ".mzml"
            else:
                raise ValueError(f"{switch} is not supported as rawfile-type")

            self.raw_files = [os.path.basename(f) for f in os.listdir(self.raw_path) if f.lower().endswith(extension)]
            logger.info(f"Found {len(self.raw_files)} raw files in the search directory")
        
    def get_msms_folder_path(self):
        return os.path.join(self.raw_path, "msms")
    
    def _get_split_msms_path(self, raw_file: str):
        return os.path.join(self.get_msms_folder_path(), os.path.splitext(raw_file)[0] + ".txt")

    def get_percolator_folder_path(self):
        return os.path.join(self.raw_path, "percolator")
    
    def _get_split_perc_input_path(self, raw_file: str):
        return os.path.join(self.get_percolator_folder_path(), os.path.splitext(raw_file)[0] + '.tab')
    
    def _get_merged_perc_input_path(self):
        return os.path.join(self.get_percolator_folder_path(), 'prosit.tab')
    
    def split_msms(self):
        if self.split_msms_step.is_done():
            return
        
        df_search = self._load_search()
        logger.info(f"Read {len(df_search.index)} PSMs from {self.search_path}")
        for raw_file, df_search_split in df_search.groupby('RAW_FILE'):
            logger.info(f"Found raw file {raw_file} in msms.txt")
            
            if not os.path.isfile(os.path.join(self.raw_path, raw_file + ".raw")):
                logger.info(f"Did not find {raw_file} in search directory, skipping this file")
                continue
            
            split_msms = self._get_split_msms_path(raw_file)
            logger.info(f"Creating split msms.txt file {split_msms}")
            df_search_split.to_csv(split_msms, sep='\t',index=False)    
        
        self.split_msms_step.mark_done()
    
    def calculate_features(self):
        num_threads = self.config['numThreads']
        if num_threads > 1:
            from .utils.multiprocessing_pool import JobPool
            processingPool = JobPool(processes = num_threads)

        for raw_file in self.raw_files:
            calc_feature_step = ProcessStep(self.raw_path, "calculate_features." + raw_file)
            if calc_feature_step.is_done():
                continue

            raw_file_path = os.path.join(self.raw_path, raw_file)
            percolator_input_path = self._get_split_perc_input_path(raw_file)
            split_msms_path = self._get_split_msms_path(raw_file)
            
            if num_threads > 1:
                processingPool.applyAsync(calculate_features_single, (raw_file_path, split_msms_path, percolator_input_path, calc_feature_step))
            else:
                calculate_features_single(raw_file_path, split_msms_path, percolator_input_path, calc_feature_step)
            
        if num_threads > 1:
            processingPool.checkPool(printProgressEvery = 1)

    def merge_input(self):
        if self.merge_input_step.is_done():
            return
        
        merged_perc_input_file = self._get_merged_perc_input_path()

        logger.info(f"Merging percolator input files")
        with open(merged_perc_input_file, "wb") as fout:
            first = True
            for percolator_input_path in map(self._get_split_perc_input_path, self.raw_files):
                with open(percolator_input_path, "rb") as f:
                    if not first:
                        next(f) # skip the header
                    else:
                        first = False
                    fout.write(f.read())
        
        self.merge_input_step.mark_done()
    
    def rescore_with_perc(self, search_type: str = "prosit", 
                          test_fdr: float = 0.01, train_fdr: float = 0.01):
        """
        Use percolator to re-score library.
        """
        if self.percolator_step.is_done():
            return

        perc_path = self.get_percolator_folder_path()
        if not os.path.isdir(perc_path):
            os.makedirs(perc_path)
        
        weights_file = os.path.join(perc_path, f"{search_type}_weights.csv")
        target_psms = os.path.join(perc_path, f"{search_type}_target.psms")
        decoy_psms = os.path.join(perc_path, f"{search_type}_decoy.psms")
        target_peptides = os.path.join(perc_path, f"{search_type}_target.peptides")
        decoy_peptides = os.path.join(perc_path, f"{search_type}_decoy.peptides")
        log_file = os.path.join(perc_path, f"{search_type}.log")

        cmd = f"percolator --weights {weights_file} \
                          --num-threads {self.config['numThreads']} \
                          --post-processing-tdc \
                          --search-input concatenated \
                          --testFDR {test_fdr} \
                          --trainFDR {train_fdr} \
                          --results-psms {target_psms} \
                          --decoy-results-psms {decoy_psms} \
                          --results-peptides {target_peptides} \
                          --decoy-results-peptides {decoy_peptides} \
                          {self._get_merged_perc_input_path()} 2> {log_file}"
        logger.info(f"Starting percolator with command {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        
        self.percolator_step.mark_done()

