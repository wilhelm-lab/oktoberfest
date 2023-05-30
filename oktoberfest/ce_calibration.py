import logging
import os
from typing import Optional

import numpy as np
import pandas as pd
from spectrum_fundamentals.annotation.annotation import annotate_spectra
from spectrum_fundamentals.metrics.similarity import SimilarityMetrics
from spectrum_io.file import csv
from spectrum_io.raw import ThermoRaw
from spectrum_io.search_result import Mascot, MaxQuant, MSFragger

from .data.spectra import FragmentType, Spectra
from .spectral_library import SpectralLibrary
from .utils.plotting import plot_mean_sa_ce

logger = logging.getLogger(__name__)


class CeCalibration(SpectralLibrary):
    """
    Main to init a CeCalibrarion obj and go through the steps.

    1- merge_mzml_and_msms
    2- allign_ce
    3- get_best_ce
    4- write output
    """

    raw_path: str
    out_path: str
    best_ce: float

    def __init__(
        self,
        search_path: str,
        raw_path: str,
        out_path: str,
        config_path: Optional[str],
        mzml_reader_package: str = "pyteomics",
    ):
        """
        Initialize a CeCalibration object.

        :param search_path: path to search directory
        :param raw_path: path to directory containing the msms.txt and raw files
        :param out_path: path to output folder
        :param config_path: path to configuration file
        :param mzml_reader_package: mzml reader (pymzml or pyteomics)
        """
        super().__init__(search_path, out_path, config_path=config_path)
        self.search_path = search_path
        self.raw_path = raw_path
        self.out_path = out_path
        self.mzml_reader_package = mzml_reader_package
        self.best_ce = 0

    def _gen_internal_search_result_from_msms(self):
        """Generate internal search result from msms.txt."""
        logger.info(f"Converting msms data at {self.search_path} to internal search result.")

        search_type = self.config.search_type
        if search_type == "maxquant":
            search_result = MaxQuant(self.search_path)
        elif search_type == "msfragger":
            search_result = MSFragger(self.search_path)
        elif search_type == "mascot":
            search_result = Mascot(self.search_path)
        else:
            raise ValueError(f"Unknown search_type provided in config: {search_type}")

        tmt_labeled = self.config.tag if any("TMT" in value for value in self.config.models.values()) else ""
        self.search_path = search_result.generate_internal(tmt_labeled=tmt_labeled)

    def _gen_mzml_from_thermo(self):
        """Generate mzml from thermo raw file."""
        logger.info("Converting thermo rawfile to mzml.")
        raw = ThermoRaw()
        if not (self.out_path.endswith(".mzML")) and (not (self.out_path.endswith(".raw"))):
            self.out_path = os.path.join(self.out_path, self.raw_path.split("/")[-1].split(".")[0] + ".mzml")
        self.raw_path = raw.convert_raw_mzml(input_path=self.raw_path, output_path=self.out_path)

    def _load_search(self):
        """Load search type."""
        switch = self.config.search_type
        logger.info(f"search_type is {switch}")
        if switch in ["maxquant", "msfragger", "mascot"]:
            self._gen_internal_search_result_from_msms()
            switch = "internal"
        if switch == "internal":
            return csv.read_file(self.search_path)
        else:
            raise ValueError(f"{switch} is not a supported search type. Convert to internal format manually.")

    def _load_rawfile(self):
        """Load raw file."""
        switch = self.config.raw_type
        search_engine = self.config.search_type
        logger.info(f"raw_type is {switch}")
        if switch == "thermo":
            self._gen_mzml_from_thermo()
        elif switch == "mzml":
            pass
        else:
            raise ValueError(f"{switch} is not supported as rawfile-type")
        print(self.raw_path)
        self.raw_path = self.raw_path.as_posix().replace(".raw", ".mzml")
        return ThermoRaw.read_mzml(source=self.out_path, package=self.mzml_reader_package, search_type=search_engine)

    def merge_mzml_and_msms(self, df_search: pd.DataFrame):
        """
        Read input search and mzml and add it to library.

        :param df_search: search result as pd.DataFrame
        """
        df_raw = self._load_rawfile()
        # return df_search
        logger.info("Merging rawfile and search result")
        df_join = df_search.merge(df_raw, on=["RAW_FILE", "SCAN_NUMBER"])
        logger.info(f"There are {len(df_join)} matched identifications")

        logger.info("Annotating raw spectra")
        df_annotated_spectra = annotate_spectra(df_join)
        df_join.drop(columns=["INTENSITIES", "MZ"], inplace=True)
        # return df_annotated_spectra["INTENSITIES"]
        logger.info("Preparing library")
        self.library.add_columns(df_join)
        self.library.add_matrix(df_annotated_spectra["INTENSITIES"], FragmentType.RAW)
        self.library.add_matrix(df_annotated_spectra["MZ"], FragmentType.MZ)
        self.library.add_column(df_annotated_spectra["CALCULATED_MASS"], "CALCULATED_MASS")

    def get_hdf5_path(self) -> str:
        """Get path to hdf5 file."""
        return self.out_path + ".hdf5"

    def get_pred_path(self) -> str:
        """Get path to prediction hdf5 file."""
        return self.out_path + "_pred.hdf5"

    def _prepare_alignment_df(self):
        self.alignment_library = Spectra()
        self.alignment_library.spectra_data = self.library.spectra_data.copy()

        # Remove decoy and HCD fragmented spectra
        self.alignment_library.spectra_data = self.alignment_library.spectra_data[
            (self.alignment_library.spectra_data["FRAGMENTATION"] == "HCD")
            & (~self.alignment_library.spectra_data["REVERSE"])
        ]
        # Select the 1000 highest scoring or all if there are less than 1000
        self.alignment_library.spectra_data = self.alignment_library.spectra_data.sort_values(
            by="SCORE", ascending=False
        ).iloc[:1000]

        # Repeat dataframe for each CE
        ce_range = range(18, 50)
        nrow = len(self.alignment_library.spectra_data)
        self.alignment_library.spectra_data = pd.concat([self.alignment_library.spectra_data for _ in ce_range], axis=0)
        self.alignment_library.spectra_data["COLLISION_ENERGY"] = np.repeat(ce_range, nrow)
        self.alignment_library.spectra_data.reset_index(inplace=True)

    def _alignment(self):
        """
        Edit library to try different ranges of ce 15-50. then predict with the new library.

        Check https://gitlab.lrz.de/proteomics/prosit_tools/oktoberfest/-/blob/develop/oktoberfest/ce_calibration/grpc_alignment.py
        """
        pred_intensity = self.alignment_library.get_matrix(FragmentType.PRED)
        raw_intensity = self.alignment_library.get_matrix(FragmentType.RAW)
        # return pred_intensity.toarray(), raw_intensity.toarray()
        sm = SimilarityMetrics(pred_intensity, raw_intensity)
        self.alignment_library.spectra_data["SPECTRAL_ANGLE"] = sm.spectral_angle(raw_intensity, pred_intensity, 0)

        self.ce_alignment = self.alignment_library.spectra_data.groupby(by=["COLLISION_ENERGY"])[
            "SPECTRAL_ANGLE"
        ].mean()

        plot_mean_sa_ce(
            sa_ce_df=self.ce_alignment,
            directory=self.raw_path.parent / "results",
            raw_file_name=self.raw_path.name,
        )

    def perform_alignment(self, df_search: pd.DataFrame):
        """
        Perform alignment and get the best CE.

        :param df_search: search result as pd.DataFrame
        """
        hdf5_path = self.get_hdf5_path()
        logger.info(f"Path to hdf5 file with annotations for {self.out_path}: {hdf5_path}")
        if os.path.isfile(hdf5_path):
            self.library.read_from_hdf5(hdf5_path)
        else:
            self.merge_mzml_and_msms(df_search)
            self.library.write_as_hdf5(hdf5_path)  # write_metadata_annotation
        # Check if all data is HCD no need to align and return the best ce as 35
        hcd_df = self.library.spectra_data[(self.library.spectra_data["FRAGMENTATION"] == "HCD")]
        if len(hcd_df.index) == 0:
            self.best_ce = 35.0
            return
        self._prepare_alignment_df()
        self.grpc_predict(self.alignment_library, alignment=True)  # predict alignment
        self._alignment()
        self.best_ce = self.ce_alignment.idxmax()
