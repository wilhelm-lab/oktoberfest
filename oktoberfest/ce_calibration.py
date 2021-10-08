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
from prosit_io.file import hdf5

import logging
logger = logging.getLogger(__name__)


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
        self.best_ce = 0


    def _gen_internal_search_result_from_msms(self):
        logger.info(f"Converting msms.txt at location {self.search_path} to internal search result.")
        mxq = MaxQuant(self.search_path)
        self.search_path = mxq.generate_internal()


    def _gen_mzml_from_thermo(self):
        logger.info("Converting thermo rawfile to mzml.")
        raw = ThermoRaw()
        self.raw_path = raw.convert_raw_mzml(self.raw_path)


    def _load_search(self):
        switch = self.config["fileUploads"]["search_type"]
        logger.info(f"search_type is {switch}")
        if switch == "maxquant":
            self._gen_internal_search_result_from_msms()
        elif switch == "internal":
            pass
        else:
            raise ValueError(f"{switch} is not supported as search-type")

        return MaxQuant.read_internal(self.search_path)


    def _load_rawfile(self):
        switch = self.config["fileUploads"]["raw_type"]
        logger.info(f"raw_type is {switch}")
        if switch == "thermo":
            ThermoRaw.convert_raw_mzml(self.raw_path)
        elif switch == "mzml":
            pass
        else:
            raise ValueError(f"{switch} is not supported as rawfile-type")

        return ThermoRaw.read_mzml(self.raw_path)


    def gen_lib(self, df_search):
        """
        Read input search and raw and add it to library
        """
        #df_search = self._load_search()
        df_raw = self._load_rawfile()
        #return df_search
        logger.info("Merging rawfile and search result")
        df_join = df_search.merge(df_raw, on=["RAW_FILE", "SCAN_NUMBER"])
        logger.info(f"There are {len(df_join)} matched identifications")

        logger.info("Annotating raw spectra")
        df_annotated_spectra = annotate_spectra(df_join)
        df_join.drop(columns=["INTENSITIES", "MZ"], inplace=True)
        #return df_annotated_spectra["INTENSITIES"]
        logger.info("Preparing library")
        self.library.add_columns(df_join)
        self.library.add_matrix(df_annotated_spectra["INTENSITIES"],FragmentType.RAW)
        self.library.add_matrix(df_annotated_spectra["MZ"],FragmentType.MZ)
        self.library.add_column(df_annotated_spectra['CALCULATED_MASS'],'CALCULATED_MASS')

    def write_metadata_annotation(self):
        #hdf5_path = os.path.join(self.out_path, raw_file_name + '.hdf5')
        hdf5_path = self.raw_path+'.hdf5'
        data_set_names = [hdf5.META_DATA_KEY, hdf5.INTENSITY_RAW_KEY, hdf5.MZ_RAW_KEY]

        sparse_matrix_intensity_raw, columns_intensity = self.library.get_matrix(FragmentType.RAW, True)
        sparse_matrix_mz, columns_mz = self.library.get_matrix(FragmentType.MZ, True)
        data_sets = [self.library.get_meta_data(), sparse_matrix_intensity_raw, sparse_matrix_mz]
        column_names = [columns_intensity, columns_mz]

        hdf5.write_file(data_sets, hdf5_path, data_set_names, column_names)


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
        #print(pred_intensity.toarray())
        #return pred_intensity.toarray(), raw_intensity.toarray()
        sm = SimilarityMetrics(pred_intensity,raw_intensity)
        self.alignment_library.spectra_data["SPECTRAL_ANGLE"] = sm.spectral_angle(raw_intensity,pred_intensity)

        self.ce_alignment = self.alignment_library.spectra_data.groupby(by=["COLLISION_ENERGY"])["SPECTRAL_ANGLE"].mean()


    def _get_best_ce(self):
        """
        Get aligned ce for this lib.
        """
        self.best_ce = self.ce_alignment.idxmax()

    def perform_alignment(self, df_search):        
        self.gen_lib(df_search)
        self.write_metadata_annotation()
        self._prepare_alignment_df()
        self._predict_alignment()
        self._alignment()
        self._get_best_ce()

if __name__ == "main":
    ce_cal = CeCalibration(search_path = "D:/Compmass/workDir/HCD_OT/msms.txt",
                      raw_path = "D:/Compmass/workDir/HCD_OT/190416_FPTMT_MS3_HCDOT_R1.mzml")
    df_search = ce_cal._load_search()
    grouped_search = df_search.groupby('RAW_FILE')
    raw_files = grouped_search.groups.keys()
    ce_cal_raw = {}
    for raw_file in raw_files:
        ce_cal_raw[raw_file] = CeCalibration(search_path="D:/Compmass/workDir/HCD_OT/msms.txt",
                                             raw_path="D:/Compmass/workDir/HCD_OT/" + raw_file + ".mzml")
        ce_cal_raw[raw_file].perform_alignment(grouped_search.get_group(raw_file))
