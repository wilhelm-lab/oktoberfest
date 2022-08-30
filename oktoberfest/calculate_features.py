import glob
import logging

import numpy as np
import pandas as pd
from fundamentals.metrics.percolator import Percolator

from .ce_calibration import CeCalibration
from .data.spectra import FragmentType
from .utils.config import Config

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
        self.library.spectra_data["COLLISION_ENERGY"] = self.best_ce
        self.grpc_predict(self.library)
        self.library.write_pred_as_hdf5(self.get_pred_path())

    def gen_perc_metrics(self, search_type, file_path=None):
        """
        get all percolator metrics and add it to library
        """

        perc_features = Percolator(
            self.library.get_meta_data(),
            self.library.get_matrix(FragmentType.PRED),
            self.library.get_matrix(FragmentType.RAW),
            search_type,
            self.config.get_all_features(),
        )
        perc_features.calc()
        if file_path:
            perc_features.write_to_file(file_path)
