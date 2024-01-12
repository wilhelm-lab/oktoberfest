import logging
from enum import Enum
from pathlib import Path
from threading import Thread
from typing import List, Tuple, Type, TypeVar, Union

import numpy as np
import pandas as pd
import scipy
import spectrum_fundamentals.constants as c
from scipy.sparse import coo_matrix, csr_matrix
from spectrum_io.file import hdf5

logger = logging.getLogger(__name__)


SpectraT = TypeVar("SpectraT", bound="Spectra")


class FragmentType(Enum):
    """FragmentType class to enumerate pred, raw, and mz."""

    PRED = 1
    RAW = 2
    MZ = 3


class Spectra:
    """Main to init spectra data."""

    INTENSITY_COLUMN_PREFIX = "INTENSITY_RAW"
    INTENSITY_PRED_PREFIX = "INTENSITY_PRED"
    MZ_COLUMN_PREFIX = "MZ_RAW"
    COLUMNS_FRAGMENT_ION = ["Y1+", "Y1++", "Y1+++", "B1+", "B1++", "B1+++"]

    spectra_data: pd.DataFrame

    def __init__(self):
        """Initialize spectra data as a pd.DataFrame."""
        self.spectra_data = pd.DataFrame()

    @staticmethod
    def _gen_column_names(fragment_type: FragmentType) -> List[str]:
        """
        Get column names of the spectra data.

        :param fragment_type: choose predicted, raw, or mz
        :return: A list of column names
        """
        prefix = Spectra._resolve_prefix(fragment_type)
        columns = []
        for i in range(1, 30):
            for column in Spectra.COLUMNS_FRAGMENT_ION:
                columns.append(prefix + "_" + column.replace("1", str(i)))
        return columns

    @staticmethod
    def _resolve_prefix(fragment_type: FragmentType) -> str:
        """
        Resolve prefix given fragment type (1 for pred, 2 for raw, 3 for mz).

        :param fragment_type: choose predicted, raw, or mz
        :return: prefix as string
        """
        logger.debug(fragment_type)
        if fragment_type.value == 1:
            prefix = Spectra.INTENSITY_PRED_PREFIX
        elif fragment_type.value == 2:
            prefix = Spectra.INTENSITY_COLUMN_PREFIX
        else:
            prefix = Spectra.MZ_COLUMN_PREFIX
        return prefix

    def add_column(self, column_data: np.ndarray, name: str) -> None:
        """
        Add column to spectra data.

        :param column_data: data for a column as np
        :param name: name of the column
        """
        self.spectra_data[name] = column_data

    def add_columns(self, columns_data: pd.DataFrame) -> None:
        """
        Add columns to spectra data.

        :param columns_data: a pandas data frame to add can be metrics or metadata
        """
        # Check if columns already exist
        self.spectra_data = pd.concat([self.spectra_data, columns_data], axis=1)

    def get_meta_data(self) -> pd.DataFrame:
        """Get meta data with intensity, mz and intensity predictions as pd.DataFrame."""
        columns = list(self.spectra_data.columns)
        meta_data_columns = list(
            filter(
                lambda x: not x.startswith(
                    tuple([Spectra.INTENSITY_COLUMN_PREFIX, Spectra.MZ_COLUMN_PREFIX, Spectra.INTENSITY_PRED_PREFIX])
                ),
                columns,
            )
        )

        return self.spectra_data[meta_data_columns]

    def add_matrix_from_hdf5(self, intensity_data: pd.DataFrame, fragment_type: FragmentType) -> None:
        """
        Concatenate intensity df as a sparse matrix to our data.

        :param intensity_data: intensity sparse matrix
        :param fragment_type: choose predicted, raw, or mz
        """
        # generate column names and build dataframe from sparse matrix
        columns = self._gen_column_names(fragment_type)
        intensity_data.columns = columns
        self.add_columns(intensity_data)

    def add_matrix(self, intensity_data: pd.Series, fragment_type: FragmentType) -> None:
        """
        Concatenate intensity df as a sparse matrix to our data.

        :param intensity_data: intensity numpy array to add
        :param fragment_type: choose predicted, raw, or mz
        """
        intensity_df = intensity_data.explode()

        # reshape based on the number of fragments
        intensity_array = intensity_df.values.astype(np.float32).reshape(-1, c.VEC_LENGTH)

        # Change zeros to epislon to keep the info of invalid values
        # change the -1 values to 0 (for better performance when converted to sparse representation)
        intensity_array[intensity_array == 0] = c.EPSILON
        intensity_array[intensity_array == -1] = 0.0

        # generate column names and build dataframe from sparse matrix
        intensity_df = pd.DataFrame.sparse.from_spmatrix(coo_matrix(intensity_array, dtype=np.float32))
        columns = self._gen_column_names(fragment_type)
        intensity_df.columns = columns
        self.add_columns(intensity_df)

    def get_columns(self, fragment_type: FragmentType, return_column_names: bool = False) -> coo_matrix:
        """
        Get intensities sparse matrix from dataframe.

        :param fragment_type: choose predicted, raw, or mz
        :param return_column_names: whether column names should be returned
        :return: sparse matrix with the required data
        """
        prefix = Spectra._resolve_prefix(fragment_type)
        logger.debug(prefix)
        columns_to_select = list(filter(lambda c: c.startswith(prefix), self.spectra_data.columns))
        if return_column_names:
            return scipy.sparse.csr_matrix(self.spectra_data[columns_to_select].values), columns_to_select
        # Check if conversion is low change to coo then csr from coo
        return self.spectra_data[columns_to_select]

    def get_matrix(self, fragment_type: FragmentType) -> Tuple[csr_matrix, List[str]]:
        """
        Get intensities sparse matrix from dataframe.

        :param fragment_type: choose predicted, raw, or mz
        :return: sparse matrix with the required data
        """
        prefix = Spectra._resolve_prefix(fragment_type)
        logger.debug(prefix)
        columns_to_select = list(filter(lambda c: c.startswith(prefix), self.spectra_data.columns))
        return scipy.sparse.csr_matrix(self.spectra_data[columns_to_select].values), columns_to_select

    def write_as_hdf5(self, output_file: Union[str, Path]) -> Thread:
        """
        Write intensity and mz data as hdf5.

        :param output_file: path to output file
        :return: the thread object from the hdf5 writer for later joining
        """
        data_set_names = [hdf5.META_DATA_KEY, hdf5.INTENSITY_RAW_KEY, hdf5.MZ_RAW_KEY]

        sparse_matrix_intensity_raw, columns_intensity = self.get_matrix(FragmentType.RAW)
        sparse_matrix_mz, columns_mz = self.get_matrix(FragmentType.MZ)
        data_sets = [self.get_meta_data(), sparse_matrix_intensity_raw, sparse_matrix_mz]
        column_names = [columns_intensity, columns_mz]

        return hdf5.write_file(data_sets, output_file, data_set_names, column_names)

    def write_pred_as_hdf5(self, output_file: Union[str, Path]) -> Thread:
        """
        Write intensity, mz, and pred data as hdf5.

        :param output_file: path to output file
        :return: the thread object from the hdf5 writer for later joining

        """
        data_set_names = [hdf5.META_DATA_KEY, hdf5.INTENSITY_RAW_KEY, hdf5.MZ_RAW_KEY, hdf5.INTENSITY_PRED_KEY]

        sparse_matrix_intensity_raw, columns_intensity = self.get_matrix(FragmentType.RAW)
        sparse_matrix_mz, columns_mz = self.get_matrix(FragmentType.MZ)
        sparse_matrix_pred, columns_pred = self.get_matrix(FragmentType.PRED)
        data_sets = [self.get_meta_data(), sparse_matrix_intensity_raw, sparse_matrix_mz, sparse_matrix_pred]
        column_names = [columns_intensity, columns_mz, columns_pred]

        return hdf5.write_file(data_sets, output_file, data_set_names, column_names)

    @classmethod
    def from_hdf5(cls: Type[SpectraT], input_file: Union[str, Path]) -> SpectraT:
        """
        Read from hdf5 file.

        :param input_file: path to input file
        :return: a spectra instance
        """
        input_file = str(input_file)
        spectra = cls()
        spectra.add_columns(hdf5.read_file(input_file, hdf5.META_DATA_KEY))
        sparse_raw_intensities = hdf5.read_file(input_file, f"sparse_{hdf5.INTENSITY_RAW_KEY}")
        if not sparse_raw_intensities.empty:
            spectra.add_matrix_from_hdf5(sparse_raw_intensities, FragmentType.RAW)
        sparse_raw_mzs = hdf5.read_file(input_file, f"sparse_{hdf5.MZ_RAW_KEY}")
        if not sparse_raw_mzs.empty:
            spectra.add_matrix_from_hdf5(sparse_raw_mzs, FragmentType.MZ)

        return spectra

    @classmethod
    def from_csv(cls: Type[SpectraT], input_file: Union[str, Path]) -> SpectraT:
        """
        Read from hdf5 file.

        :param input_file: path to input file
        :return: a spectra instance
        """
        input_file = str(input_file)
        spectra = cls()
        spectra.spectra_data = pd.read_csv(input_file)

        return spectra
