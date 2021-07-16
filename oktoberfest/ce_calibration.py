import os

import numpy as np
import pandas as pd
import scipy.sparse

from prosit_io.search_result import MaxQuant
from prosit_io.raw import ThermoRaw
from prosit_io.file import csv
from fundamentals.annotation.annotation import annotate_spectra
from fundamentals.metrics.metric import Metric
from fundamentals.metrics.similarity import SimilarityMetrics

from .spectral_library import SpectralLibrary
from .data.spectra import Spectra
from .data.spectra import FragmentType

class CeCalibration(SpectralLibrary):
    """
        main to init a CeCalibrarion obj and go through the steps:
        1- gen_lib
        2- allign_ce
        3- get_best_ce
        4- write output
    """
    search_path: str
    raw_path: str
    best_ce: float

    def __init__(self, search_path, raw_path):
        super().__init__(search_path)
        self.search_path = search_path
        self.raw_path = raw_path
        self.library = Spectra()

    def _gen_internal_search_result_from_msms(self):
        mxq = MaxQuant(self.search_path)
        self.search_path = mxq.generate_internal()

    def _gen_mzml_from_thermo(self):
        raw = ThermoRaw()
        self.raw_path = raw.raw_mzml(self.raw_path)

    def _load_search(self):
        switch = self.config["fileUploads"]["search_type"]
        if switch == "maxquant":
            self._gen_internal_search_result_from_msms()
        elif switch == "internal":
            pass
        else:
            raise ValueError(f"{switch} is not supported as search-type")

        return MaxQuant.read_internal(self.search_path)


    def _load_rawfile(self):
        switch = self.config["fileUploads"]["raw_type"]
        if switch == "thermo":
            ThermoRaw.raw_mzml(self.raw_path)
        elif switch == "mzml":
            pass
        else:
            raise ValueError(f"{switch} is not supported as rawfile-type")

        return ThermoRaw.read_mzml(self.raw_path)


    def gen_lib(self):
        """
        Read input search and raw and add it to library
        """

        df_search = self._load_search()
        df_raw = self._load_rawfile()

        df_join = df_search.merge(df_raw, on=["RAW_FILE", "SCAN_NUMBER"])
        df_annotated_spectra = annotate_spectra(df_join)
        df_join.drop(columns=["INTENSITIES", "MZ"], inplace=True)
        self.library.add_columns(df_join)
        #return df_annotated_spectra
        self.library.add_matrix(df_annotated_spectra["INTENSITIES"],FragmentType.RAW)
        self.library.add_matrix(df_annotated_spectra["MZ"],FragmentType.MZ)


    def _prepare_alignment_df(self):
        self.alignment_library = Spectra()
        self.alignment_library.spectra_data = self.library.spectra_data.copy()

        # Remove decoy and HCD fragmented spectra
        self.alignment_library.spectra_data = self.alignment_library.spectra_data[(self.alignment_library.spectra_data["FRAGMENTATION"] == "HCD") & (~self.alignment_library.spectra_data["REVERSE"])]
        # Select the 1000 highest scoring or all if there are less than 1000
        self.alignment_library.spectra_data = self.alignment_library.spectra_data.sort_values(by="SCORE", ascending=False).iloc[:1000]

        # Old code to repeat dataframe for each CE
        # self.alignment_library.spectra_data["COLLISION_ENERGY"] = np.array([x for x in range(18,50)] for _ in range(len(self.alignment_library.spectra_data)))
        # self.alignment_library.spectra_data = self.alignment_library.spectra_data.explode("COLLISION_ENERGY")
        # self.alignment_library.spectra_data.reset_index(inplace=True)
        # self.alignment_library.spectra_data = self.alignment_library.spectra_data.copy()

        # Repeat dataframe for each CE
        CE_RANGE = range(18,50)
        nrow = len(self.alignment_library.spectra_data)
        self.alignment_library.spectra_data = pd.concat([self.alignment_library.spectra_data for _ in CE_RANGE])
        self.alignment_library.spectra_data["COLLISION_ENERGY"] = np.repeat(CE_RANGE, nrow)
        self.alignment_library.spectra_data.reset_index(inplace=True)

    def _predict_alignment(self):
        self.grpc_predict(self.alignment_library)

    def _alignment(self):
        """
        Edit library to try different ranges of ce 15-50.
        then predict with the new library.
        Check https://gitlab.lrz.de/proteomics/prosit_tools/oktoberfest/-/blob/develop/oktoberfest/ce_calibration/grpc_alignment.py
        """
        pred_intensity = self.alignment_library.get_matrix(FragmentType.PRED)
        raw_intensity = self.alignment_library.get_matrix(FragmentType.RAW)

        sm = SimilarityMetrics(pred_intensity,raw_intensity)
        self.alignment_library.spectra_data["SPECTRAL_ANGLE"] = sm.spectral_angle(raw_intensity,pred_intensity)

        self.ce_alignment = self.alignment_library.spectra_data.groupby(by=["COLLISION_ENERGY"])["SPECTRAL_ANGLE"].mean()

    def _get_best_ce(self):
        """
        Get aligned ce for this lib.
        """
        return self.ce_alignment.idxmax()


    def perform_alignment(self):
        self.gen_lib()
        self._prepare_alignment_df()
        self._predict_alignment()
        self._alignment()
        return self._get_best_ce()


if __name__ == "main":
    pass
