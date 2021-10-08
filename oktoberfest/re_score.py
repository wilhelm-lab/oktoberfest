import numpy as np
import pandas as pd
import glob


from .ce_calibration import CeCalibration
from fundamentals.metrics.percolator import Percolator
from .data.spectra import FragmentType


class ReScore(CeCalibration):
    """
        main to init a re-score obj and go through the steps:
        1- predict_aligned_ce
        2- gen_perc_metrics
        3- rescore_with_perc
        4- write output
    """
    percolator: pd.DataFrame
    def predict_with_aligned_ce(self, df_search):
        """
        Get best ce with ce_calibration then use it for prediction.
        """
        self.perform_alignment(df_search)
        self.library.spectra_data['COLLISION_ENERGY'] = self.best_ce
        self.grpc_predict(self.library)

    def gen_perc_metrics(self):
        """
        get all percolator metrics and add it to library
        """

        perc_features = Percolator(self.library.get_meta_data(),
                                   self.library.get_matrix(FragmentType.PRED),
                                   self.library.get_matrix(FragmentType.RAW),
                                   'Prosit')
        perc_features.calc()
        return perc_features.metrics_val

    def rescore_with_perc(self):
        """
        Use percolator to re-score library.
        """
        pass

    def align_percolator_cols(self):
        all_columns = self.percolator.columns
        first_columns = ['SpecId', 'Label', 'ScanNr']
        last_columns = ['Peptide', 'Protein']
        mid_columns = list(set(all_columns) - set(first_columns) - set(last_columns))
        new_columns = first_columns + mid_columns + last_columns
        return new_columns

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
        re_score_raw[raw_file].percolator = ce_cal_raw[raw_file].gen_perc_metrics()
        re_score_raw[raw_file].percolator = ce_cal_raw.percolator[ReScore.align_percolator_cols()]
        re_score_raw[raw_file].percolator.to_csv('prosit_' + raw_file +'.tab', sep='\t',index=False)

    raw_files_list = list(raw_files)
    with open("merged_prosit.tab", "wb") as fout:
        # first file:
        with open("prosit_" + raw_files_list[0] + ".tab", "rb") as f:
            fout.write(f.read())
        # now the rest:
        for raw_file in list(raw_files)[1:]:
            with open("prosit_" + rawfile + ".tab", "rb") as f:
                next(f)  # skip the header
                fout.write(f.read())