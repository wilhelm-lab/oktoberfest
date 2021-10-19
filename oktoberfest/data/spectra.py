import logging
from enum import Enum

import pandas as pd
import scipy
from scipy.sparse import coo_matrix, spmatrix
import numpy as np
import fundamentals.constants as C
from prosit_io.file import hdf5

logger = logging.getLogger(__name__)


class FragmentType(Enum):
    PRED = 1
    RAW = 2
    MZ = 3


class Spectra:
    NO_OF_FRAGMENTS = 174
    INTENSITY_COLUMN_PREFIX = 'INTENSITY_RAW'
    INTENSITY_PRED_PREFIX = 'INTENSITY_PRED'
    MZ_COLUMN_PREFIX = 'MZ_RAW'
    EPSILON = 1e-7
    COLUMNS_FRAGMENT_ION = ['Y1+', 'Y1++', 'Y1+++', 'B1+', 'B1++', 'B1+++']

    spectra_data: pd.DataFrame

    def __init__(self):
        self.spectra_data = pd.DataFrame()

    @staticmethod
    def _gen_column_names(fragment_type: FragmentType):
        prefix = Spectra._resolve_prefix(fragment_type)
        columns = []
        for i in range(1, 30):
            for column in Spectra.COLUMNS_FRAGMENT_ION:
                columns.append(prefix + '_' + column.replace('1', str(i)))
        return columns

    @staticmethod
    def _resolve_prefix(fragment_type):
        logger.debug(fragment_type)
        if fragment_type.value == 1:
            prefix = Spectra.INTENSITY_PRED_PREFIX
        elif fragment_type.value == 2:
            prefix = Spectra.INTENSITY_COLUMN_PREFIX
        else:
            prefix = Spectra.MZ_COLUMN_PREFIX
        return prefix

    def add_column(self, column_data: np, name: str):
        self.spectra_data[name] = column_data

    def add_columns(self, columns_data: pd.DataFrame):
        """
        Add columns to spectra data.
        :param columns_data: a pandas data frame to add can be metrics or metadata.
        """
        # Check if columns already exist
        self.spectra_data = pd.concat([self.spectra_data, columns_data], axis=1)

    def get_meta_data(self):
        columns = self.spectra_data.columns
        meta_data_columns = list(filter(lambda x: not x.startswith(
            tuple([Spectra.INTENSITY_COLUMN_PREFIX, Spectra.MZ_COLUMN_PREFIX, Spectra.INTENSITY_PRED_PREFIX])), columns))

        return self.spectra_data[meta_data_columns]

    def add_matrix(self, intensity_data, fragment_type):
        """
        concat intensity df as a sparse matrix to our data
        :param intensity_data: Intensity numpy array to add
        :param fragment_type: Choose type of fragments predicted or raw
        """
        intensity_df = intensity_data.explode()

        # reshape based on the number of fragments
        intensity_array = intensity_df.values.astype(np.float32).reshape(-1, Spectra.NO_OF_FRAGMENTS)

        # Change zeros to epislon to keep the info of invalid values
        # change the -1 values to 0 (for better performance when converted to sparse representation)
        intensity_array[intensity_array == 0] = Spectra.EPSILON
        intensity_array[intensity_array == -1] = 0

        # generate column names and build dataframe from sparse matrix
        intensity_df = pd.DataFrame.sparse.from_spmatrix(coo_matrix(intensity_array)).astype(np.float32)
        columns = self._gen_column_names(fragment_type)
        intensity_df.columns = columns
        self.add_columns(intensity_df)

    def get_columns(self, fragment_type, return_column_names=False) -> spmatrix:
        """
        Get intensities sparse matrix from dataframe.
        :param fragment_type: choose predicted, raw, or mz
        :return: sparse matrix with the required data
        """

        prefix = Spectra._resolve_prefix(fragment_type)
        logger.debug(prefix)
        columns_to_select = list(filter(lambda c: c.startswith(prefix), self.spectra_data.columns))
        if return_column_names:
            return scipy.sparse.csr_matrix(self.spectra_data[columns_to_select].values), columns_to_select
        # Check if conversion is low change to coo then csr from coo
        return self.spectra_data[columns_to_select]

    def get_matrix(self, fragment_type, return_column_names=False) -> spmatrix:
        """
        Get intensities sparse matrix from dataframe.
        :param fragment_type: choose predicted, raw, or mz
        :return: sparse matrix with the required data
        """

        prefix = Spectra._resolve_prefix(fragment_type)
        logger.debug(prefix)
        columns_to_select = list(filter(lambda c: c.startswith(prefix), self.spectra_data.columns))
        if return_column_names:
            return scipy.sparse.csr_matrix(self.spectra_data[columns_to_select].values), columns_to_select
        # Check if conversion is low change to coo then csr from coo
        return scipy.sparse.csr_matrix(self.spectra_data[columns_to_select].values)
    
    def write_as_hdf5(self, output_file):
        data_set_names = [hdf5.META_DATA_KEY, hdf5.INTENSITY_RAW_KEY, hdf5.MZ_RAW_KEY]

        sparse_matrix_intensity_raw, columns_intensity = self.get_matrix(FragmentType.RAW, True)
        sparse_matrix_mz, columns_mz = self.get_matrix(FragmentType.MZ, True)
        data_sets = [self.get_meta_data(), sparse_matrix_intensity_raw, sparse_matrix_mz]
        column_names = [columns_intensity, columns_mz]

        hdf5.write_file(data_sets, output_file, data_set_names, column_names)
    
    def read_from_hdf5(self, input_file):
        self.add_columns(hdf5.read_file(input_file, hdf5.META_DATA_KEY))
        self.add_columns(hdf5.read_file(input_file, f"sparse_{hdf5.INTENSITY_RAW_KEY}"))
        self.add_columns(hdf5.read_file(input_file, f"sparse_{hdf5.MZ_RAW_KEY}"))

