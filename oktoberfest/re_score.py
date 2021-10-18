import os
import logging
from typing import List
import subprocess

import numpy as np
import pandas as pd
import tracemalloc


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
    features: CalculateFeatures

    def calculate_features_single(self, raw_file_path, split_msms_path, percolator_input_path):
        logger.info(f"Calculating features for {raw_file_path}")
        self.features = CalculateFeatures(search_path="",
                                     raw_path=raw_file_path)
        df_search = pd.read_csv(split_msms_path, delimiter='\t')
        self.features.predict_with_aligned_ce(df_search)
        self.features.gen_perc_metrics(percolator_input_path)

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
        return os.path.join(self.raw_path, "msms_" + os.path.splitext(raw_file)[0] + ".prosit")
    
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
            print(df_search_split['PEPTIDE_LENGTH'].dtype)
            df_search_split = df_search_split[(df_search_split['PEPTIDE_LENGTH'] <= 30)]
            df_search_split = df_search_split[(~df_search_split['MODIFIED_SEQUENCE'].str.contains('\(ac\)'))]
            df_search_split = df_search_split[(~df_search_split['SEQUENCE'].str.contains('U'))]
            df_search_split = df_search_split[df_search_split['PRECURSOR_CHARGE'] <= 6]
            df_search_split = df_search_split[df_search_split['PEPTIDE_LENGTH'] >= 7]
            df_search_split.to_csv(split_msms, sep='\t',index=False)    
        
        open(split_msms_done_file, 'w').close()
    
    def calculate_features(self, num_threads = 5):
        if num_threads > 1:
            tracemalloc.start()
            from . import multiprocessing_pool as pool
            processingPool = pool.MyPool(processes = num_threads)
        for raw_file in self.raw_files:
            raw_file_path = os.path.join(self.raw_path, raw_file)
            percolator_input_path = self._get_split_prosit_path(raw_file)
            split_msms_path = self._get_split_msms_path(raw_file)
            
            if not os.path.isfile(percolator_input_path):
                if num_threads > 1:
                    processingPool.applyAsync(self.calculate_features_single, (raw_file_path, split_msms_path, percolator_input_path))
                else:
                    self.calculate_features_single(raw_file_path, split_msms_path, percolator_input_path)
            else:
                logger.info(f"Found percolator input file {percolator_input_path} for raw file {raw_file}, skipping feature calculation.")
            
        if num_threads > 1:
            processingPool.checkPool(printProgressEvery = 1)

    def merge_input(self, merged_input_path: str):
        logger.info(f"Merging percolator input files")
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
        
        # TODO: read number of threads from config file, otherwise percolator will use the default maximum of 3
        cmd = f"percolator --weights {weights_file} \
                          --num-threads 6 \
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
        pass



if __name__ == "main":
    ce_cal = CeCalibration(search_path = "D:/Compmass/workDir/HCD_OT/msms.txt",
                          raw_path = "D:/Compmass/workDir/HCD_OT/190416_FPTMT_MS3_HCDOT_R1.mzml")
    df_search = ce_cal._load_search()
    grouped_search = df_search.groupby('RAW_FILE')
    raw_files = grouped_search.groups.keys()
    re_score_raw = {}
    for raw_file in raw_files:
        re_score_raw[raw_file] = ReScore(search_path="D:/Compmass/workDir/HCD_OT/msms.txt",
                                             raw_path="D:/Compmass/workDir/HCD_OT/" + raw_file + ".mzml")
        msms_raw = grouped_search.get_group(raw_file)
        msms_raw = msms_raw[~msms_raw['SEQUENCE'].str.contains('U')]
        msms_raw = msms_raw[msms_raw['PRECURSOR_CHARGE']<=6]
        re_score_raw[raw_file].predict_with_aligned_ce(msms_raw)
        re_score_raw[raw_file].percolator = re_score_raw[raw_file].gen_perc_metrics()
        re_score_raw[raw_file].percolator = re_score_raw.percolator[ReScore.align_percolator_cols()]
        re_score_raw[raw_file].percolator.to_csv('prosit_' + raw_file +'.tab', sep='\t',index=False)

