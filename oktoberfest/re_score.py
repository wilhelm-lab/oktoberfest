import numpy as np
import pandas as pd

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

    def predict_with_aligned_ce(self):
        """
        Get best ce with ce_calibration then use it for prediction.
        """
        self.perform_alignment()
        self.library.spectra_data['COLLISION_ENERGY'] = self.best_ce
        self.grpc_predict(self.library)

    def gen_perc_metrics(self):
        """
        get all percolator metrics and add it to library
        """

        perc_features = Percolator(self.library.get_meta_data, self.library.get_matrix(FragmentType.PRED), self.library.get_matrix(FragmentType.RAW),
                                    'Prosit')
        perc_features.calc()
        return perc_features.metrics_val

    def rescore_with_perc(self):
        """
        Use percolator to re-score library.
        """
        pass
