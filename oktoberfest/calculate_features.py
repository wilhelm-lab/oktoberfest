import numpy as np
import pandas as pd
import glob
import logging

from .ce_calibration import CeCalibration
from fundamentals.metrics.percolator import Percolator
from .data.spectra import FragmentType

logger = logging.getLogger(__name__)


class CalculateFeatures(CeCalibration):
    """
        main to init a re-score obj and go through the steps:
        1- predict_with_aligned_ce
        2- gen_perc_metrics
    """
    def predict_with_aligned_ce(self, df_search):
        """
        Get best ce with ce_calibration then use it for prediction.
        """
        self.perform_alignment(df_search)
        self.library.spectra_data['COLLISION_ENERGY'] = self.best_ce
        self.grpc_predict(self.library)

    def gen_perc_metrics(self, file_path = None):
        """
        get all percolator metrics and add it to library
        """

        perc_features = Percolator(self.library.get_meta_data(),
                                   self.library.get_matrix(FragmentType.PRED),
                                   self.library.get_matrix(FragmentType.RAW),
                                   'Prosit')
        perc_features.calc()
        if file_path:
            perc_features.write_to_file(file_path)
