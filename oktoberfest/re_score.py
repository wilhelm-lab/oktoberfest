import os
import logging
from typing import List
import subprocess

import numpy as np
import pandas as pd

from .calculate_features import CalculateFeatures

logger = logging.getLogger(__name__)


class ReScore(CalculateFeatures):
    """
        main to init a re-score obj and go through the steps:
        1- predict_with_aligned_ce
        2- gen_perc_metrics
        3- rescore_with_perc
        4- write output
    """
    raw_path: str
    raw_files: List[str]
    
    def get_raw_files(self):
        self.raw_files = []
        if os.path.isfile(self.raw_path):
            self.raw_files = [self.raw_path]
            self.raw_path = os.path.dirname(self.raw_path)
        elif os.path.isdir(self.raw_path):
            # TODO: read extension from config file
            self.raw_files = [os.path.basename(f) for f in os.listdir(self.raw_path) if f.endswith('.raw')]
            logger.info(f"Found {len(self.raw_files)} raw files in the search directory")
    
    def _get_split_msms_path(self, raw_file: str):
        return os.path.join(self.raw_path, "msms_" + os.path.splitext(raw_file)[0] + ".txt")
    
    def _get_split_prosit_path(self, raw_file: str):
        return os.path.join(self.raw_path, 'prosit_' + os.path.splitext(raw_file)[0] + '.tab')
    
    def split_msms(self):
        split_msms_done_file = os.path.join(self.raw_path, 'split_msms.done')
        if os.path.isfile(split_msms_done_file):
            logger.info(f"Found {split_msms_done_file}, skipping split_msms")
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
        
        open(split_msms_done_file, 'w').close()
    
    def calculate_features(self):
        for raw_file in self.raw_files:
            raw_file_path = os.path.join(self.raw_path, raw_file)
            percolator_input_path = self._get_split_prosit_path(raw_file)
            
            if not os.path.isfile(percolator_input_path):
                logger.info(f"Calculating features for {raw_file}")
                features = CalculateFeatures(search_path = "",
                                  raw_path = raw_file_path)
                df_search = pd.read_csv(self._get_split_msms_path(raw_file), delimiter = '\t')
                features.predict_with_aligned_ce(df_search)
                features.gen_perc_metrics(percolator_input_path)
            else:
                logger.info(f"Found percolator input file {percolator_input_path} for raw file {raw_file}, skipping feature calculation.")
        
    def merge_input(self, merged_input_path: str):
        with open(merged_input_path, "wb") as fout:
            first = True
            for percolator_input_path in map(self._get_split_prosit_path, self.raw_files):
                with open(percolator_input_path, "rb") as f:
                    if not first:
                        next(f) # skip the header
                    else:
                        first = False
                    fout.write(f.read())
    
    def rescore_with_perc(self, percolator_input_file: str, output_path: str, search_type: str = "prosit", 
                          test_fdr: float = 0.01, train_fdr: float = 0.01):
        """
        Use percolator to re-score library.
        """
        
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        
        weights_file = os.path.join(output_path, f"{search_type}_weights.csv")
        target_psms = os.path.join(output_path, f"{search_type}_target.psms")
        decoy_psms = os.path.join(output_path, f"{search_type}_decoy.psms")
        target_peptides = os.path.join(output_path, f"{search_type}_target.peptides")
        decoy_peptides = os.path.join(output_path, f"{search_type}_decoy.peptides")
        log_file = os.path.join(output_path, f"{search_type}.log")
        
        # TODO: add number of threads as a parameter, otherwise it will use the default maximum of 3
        cmd = f"percolator --weights {weights_file} \
                          --post-processing-tdc \
                          --search-input concatenated \
                          --testFDR {test_fdr} \
                          --trainFDR {train_fdr} \
                          --results-psms {target_psms} \
                          --decoy-results-psms {decoy_psms} \
                          --results-peptides {target_peptides} \
                          --decoy-results-peptides {decoy_peptides} \
                          {percolator_input_file} 2> {log_file}"
        logger.info(f"Starting percolator with command {cmd}")
        subprocess.run(cmd, shell=True, check=True)

