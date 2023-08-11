import logging
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from spectrum_fundamentals.metrics.percolator import Percolator

from .ce_calibration import CeCalibration
from .data.spectra import FragmentType

logger = logging.getLogger(__name__)


class CalculateFeatures(CeCalibration):
    """
    Main to init a re-score obj and go through the steps.

    1- predict_with_aligned_ce
    2- gen_perc_metrics
    """

    def predict_with_aligned_ce(self, df_search: pd.DataFrame):
        """
        Get best collision energy with ce_calibration then use it for prediction.

        :param df_search: a msms matrix as a pd.DataFrame
        """
        self.perform_alignment(df_search)
        self.library.spectra_data["COLLISION_ENERGY"] = self.best_ce
        self.grpc_predict(self.library)
        self.library.write_pred_as_hdf5(self.get_pred_path())

    def gen_perc_metrics(self, search_type: str, file_path: Optional[Union[str, Path]] = None):
        """
        Get all percolator / mokapot metrics and add them to library.

        :param search_type: model (rescore or original)
        :param file_path: Optional path to percolator / mokapot input file
        """
        perc_features = Percolator(
            metadata=self.library.get_meta_data(),
            pred_intensities=self.library.get_matrix(FragmentType.PRED),
            true_intensities=self.library.get_matrix(FragmentType.RAW),
            mz=self.library.get_matrix(FragmentType.MZ),
            input_type=search_type,
            all_features_flag=self.config.all_features,
            regression_method=self.config.curve_fitting_method,
        )
        perc_features.calc()
        if file_path:
            perc_features.write_to_file(str(file_path))
