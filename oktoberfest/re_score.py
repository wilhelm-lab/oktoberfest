import numpy as np
import pandas as pd

from .ce_calibration import CeCalibrarion

class ReScore(CeCalibrarion):
    """
        main to init a re-score obj and go through the steps:
        1- predict_aligned_ce
        2- gen_perc_metrics
        3- rescore_with_perc
        4- write output
    """
    def predict_aligned_ce(self):
        """
        Get best ce with ce_calibration then use it for prediction.
        """
        pass

    def gen_perc_metrics(self):
        """
        get all percolator metrics and add it to library
        """
        pass

    def rescore_with_perc(self):
        """
        Use percolator to re-score library.
        """
        pass
