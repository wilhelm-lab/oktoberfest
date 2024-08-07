import logging
from enum import Enum
from pathlib import Path
from typing import List, Optional, Tuple, Type, TypeVar, Union

import anndata
import numpy as np
import pandas as pd
import scipy
import spectrum_fundamentals.constants as c
from scipy.sparse import csr_matrix, dok_matrix

logger = logging.getLogger(__name__)


SpectraT = TypeVar("SpectraT", bound="Spectra")


class FragmentType(Enum):
    """FragmentType class to enumerate pred, raw, and mz."""

    PRED = 1
    RAW = 2
    MZ = 3


class Spectra(anndata.AnnData):
    """Main to init spectra data."""

    INTENSITY_COLUMN_PREFIX = "INTENSITY_RAW"
    INTENSITY_PRED_PREFIX = "INTENSITY_PRED"
    MZ_COLUMN_PREFIX = "MZ_RAW"
    INTENSITY_PRED_LAYER_NAME = "pred_int"
    INTENSITY_LAYER_NAME = "raw_int"
    MZ_LAYER_NAME = "mz"
    COLUMNS_FRAGMENT_ION = ["Y1+", "Y1++", "Y1+++", "B1+", "B1++", "B1+++"]
    MAX_CHARGE = 3

    @staticmethod
    def _gen_vars_df(specified_ion_types: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Creates Annotation dataframe for vars in AnnData object.

        :param specified_ion_types: ion types that are expected to be in the spectra. If None defaults to y,b
        :return: pd.Dataframe of Frgment Annotation
        """
        if not specified_ion_types:
            specified_ion_types = ["y", "b"]

        number_of_ion_types = len(specified_ion_types)
        ion_nums = np.repeat(np.arange(1, 30, dtype=np.int32), 3 * number_of_ion_types)
        ion_charge = np.tile(np.arange(1, 4, dtype=np.int32), 29 * number_of_ion_types)
        temp_cols = []
        for size in range(1, 30):
            for type in specified_ion_types:
                for charge in ["+1", "+2", "+3"]:
                    temp_cols.append(f"{type}{size}{charge}")
        ion_types = [frag[0] for frag in temp_cols]
        var_df = pd.DataFrame({"ion": temp_cols, "num": ion_nums, "type": ion_types, "charge": ion_charge})
        var_df = var_df.set_index("ion")
        return var_df

    @staticmethod
    def _gen_column_names(fragment_type: FragmentType):  # , fragmentation_methods: Set[str]) -> List[str]:
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
        elif fragment_type.value == 3:
            prefix = Spectra.MZ_COLUMN_PREFIX
        return prefix

    @staticmethod
    def _resolve_layer_name(fragment_type: FragmentType) -> str:
        """
        Resolve layer name given fragment type (1 for pred, 2 for raw, 3 for mz).

        :param fragment_type: choose predicted, raw, or mz
        :return: layer name as string
        """
        if fragment_type.value == 1:
            layer = Spectra.INTENSITY_PRED_LAYER_NAME
        elif fragment_type.value == 2:
            layer = Spectra.INTENSITY_LAYER_NAME
        elif fragment_type.value == 3:
            layer = Spectra.MZ_LAYER_NAME
        return layer

    def __getitem__(self, index: anndata._core.index.Index):
        """Returns a sliced view of the object with this type to avoid returning AnnData instances when slicing."""
        oidx, vidx = self._normalize_indices(index)
        return Spectra(self, oidx=oidx, vidx=vidx, asview=True)

    def add_column(self, data: Union[np.ndarray, pd.Series], name: Optional[str] = None) -> None:
        """
        Add column to spectra data.

        :param data: data for a column as np
        :param name: Optional name of the column, required if providing data as a numpy array.
            In case of a pd.Series, providing a name replaces the series' name.

        :raises AssertionError: if data is not 1-dimensional or a column name is missing when
            providing a numpy array.
        :raises TypeError: if the data type is not understood
        """
        if isinstance(data, np.ndarray):
            if name is None:
                raise AssertionError("Missing column name.")
            if data.ndim != 1:
                raise AssertionError("Column data must be supplied as a 1D numpy array.")
            self.obs[name] = data
        elif isinstance(data, pd.Series):
            self.obs[name or data.name] = data
        else:
            raise TypeError(f"Unsupported data type provided: {type(data)}")

    def get_meta_data(self) -> pd.DataFrame:
        """Get meta data with intensity, mz and intensity predictions as pd.DataFrame."""
        return self.obs

    def add_matrix_from_hdf5(self, intensity_data: pd.DataFrame, fragment_type: FragmentType) -> None:
        """
        Concatenate intensity df as a sparse matrix to our data.

        :param intensity_data: intensity sparse matrix
        :param fragment_type: choose predicted, raw, or mz
        """
        # add sparse matrix of intensities to corresponding layer
        layer = self._resolve_layer_name(fragment_type)
        self.layers[layer] = scipy.sparse.csr_matrix(intensity_data)

    def add_mzs(self, mzs: np.ndarray, fragment_type: FragmentType):
        """
        Add mass to charge ratios.

        This function adds a matrix of mass to charge ratios of shape (PSMs x fragment ions)
        to this data object.

        :param mzs: the mass to charge ratio array
        :param fragment_type: the type of mzs to add. Currently, only FragmentType.MZ is supported.
        """
        layer = self._resolve_layer_name(fragment_type)
        self.layers[layer] = csr_matrix(mzs)

    def add_intensities(self, intensities: np.ndarray, annotation: np.ndarray, fragment_type: FragmentType):
        """
        Add predicted intensities and convert to sparse matrix.

        This function takes two numpy arrays, containing intensities and the fragment ion annotations
        in ProForma notation representing the column index.
        Each intensity array is reordered using the annotation to match the order of the
        fragment ion annotations in self.var_names and stored as a csr_matrix.

        :param intensities: intensity numpy array to add with shapes (n x m)
        :param annotation: fragment ion annotation numpy array in ProForma notation with shape (... x m). Only
            the first row of the annotation array is used, i.e. intensities for all PSMs must be provided in
            the same order.
        :param fragment_type: the type of intensities to add. Can be FragmentType.RAW or FragmentType.PRED.
        """
        intensities[intensities == 0] = c.EPSILON
        intensities[intensities == -1] = 0.0

        annotation_to_index = {annot: index for index, annot in enumerate(self.var_names)}
        col_index = np.vectorize(annotation_to_index.get)(annotation[0].astype(str))
        sparse_intensity_matrix = dok_matrix(self.shape)
        sparse_intensity_matrix[:, col_index] = intensities

        layer = self._resolve_layer_name(fragment_type)
        self.layers[layer] = csr_matrix(sparse_intensity_matrix)

    def add_list_of_predicted_intensities(
        self,
        intensities: List[np.ndarray],
        annotations: List[np.ndarray],
        chunk_indices: List[np.ndarray],
    ):
        """
        Add chunks of predicted intensities and convert to sparse matrix.

        This function takes three lists of numpy arrays, containing intensities, the fragment ion annotations
        in ProForma notation representing the column index, and a numeric index representing the row index.
        Each intensity array is reordered using the corresponding annotation element to match the order of the
        fragment ion annotations in self.var_names and stored to the appropriate rows of a dok_matrix,
        incrementally creating the full, sparse intensity matrix ordered by fragment types. The function then
        converts the matrix to csr format.

        :param intensities: List of intensity numpy arrays to add with shapes (n_1 x m_1), ..., (n_N x m_N)
        :param annotations: List of fragment ion annotations in ProForma notation with shapes (m_1), ..., (m_N)
        :param chunk_indices: List of row numbers with shapes (n_1), ..., (n_N)
        """
        sparse_intensity_matrix = dok_matrix(self.shape)

        for ints, annots, chunk in zip(intensities, annotations, chunk_indices):
            self._add_predicted_intensites(
                mat=sparse_intensity_matrix,
                intensity_data=ints,
                annotation=annots.astype(str),
                index=chunk,
            )

        layer = self._resolve_layer_name(FragmentType.PRED)
        self.layers[layer] = csr_matrix(sparse_intensity_matrix)

    def _add_predicted_intensites(
        self,
        mat: dok_matrix,
        intensity_data: np.ndarray,
        annotation: np.ndarray,
        index: np.ndarray,
    ) -> None:
        """
        Concatenate intensity df as a sparse matrix to our data.

        :param mat: The dok_matrix into which to store the data
        :param intensity_data: Intensity numpy array to add with shape (n x m)
        :param annotation: Fragment ion annotations in ProForma notation with shape (m)
        :param index: Row numbers with shape (n)
        """
        # ensure intensities are properly masked where required (alphapept does not do that)
        annotation_to_index = {annotation: index for index, annotation in enumerate(self.var_names)}
        col_index = np.vectorize(annotation_to_index.get)(annotation[0])
        fragment_charges = self.var.loc[annotation[0], "charge"].values
        precursor_charges = self.obs.iloc[index][["PRECURSOR_CHARGE"]].values
        intensity_data = np.where(fragment_charges <= precursor_charges, intensity_data, -1)
        row_index = self.obs.index.get_indexer(index)[..., None]
        # Change zeros to epsilon to keep the info of invalid values
        # change the -1 values to 0 (for better performance when converted to sparse representation)
        intensity_data[intensity_data == 0] = c.EPSILON
        intensity_data[intensity_data == -1] = 0.0
        mat[row_index, col_index] = intensity_data

        # self.obs.iloc[index]["done"] = True

    def get_matrix(self, fragment_type: FragmentType) -> Tuple[csr_matrix, List[str]]:
        """
        Get intensities sparse matrix from AnnData object.

        :param fragment_type: choose predicted, raw, or mz
        :return: sparse matrix with the required data
        """
        prefix = self._resolve_prefix(fragment_type)
        logger.debug(prefix)

        layer = self._resolve_layer_name(fragment_type)
        matrix = self.layers[layer]

        return matrix, self._gen_column_names(fragment_type)  # , set(self.obs["FRAGMENTATION"]))

    def write_as_hdf5(self, output_file: Union[str, Path]):
        """
        Write spectra_data to hdf5 file.

        :param output_file: path to output file
        """
        self.write(output_file, compression="gzip")

    @classmethod
    def from_hdf5(cls: Type[SpectraT], input_file: Union[str, Path]) -> SpectraT:
        """
        Read from hdf5 file.

        :param input_file: path to input file
        :return: a spectra instance
        """
        return cls(anndata.read_h5ad(str(input_file)))

    def remove_decoys(self):
        """Remove decoys in-place."""
        self.__dict__ = Spectra(self[~self.obs.REVERSE].copy()).__dict__

    def filter_by_score(self, threshold: float):
        """Filter out peptides with Andromeda score below threshold in-place."""
        self.__dict__ = Spectra(self[self.obs.SCORE >= threshold].copy()).__dict__

    def convert_to_df(self) -> pd.DataFrame:
        """
        Gives back spectra_data instance as a pandas Dataframe.

        :return: a pandas DataFrame
        """
        df_merged = self.obs
        logger.debug(self.obs.columns)

        if "mz" in list(self.layers):
            mz_cols = pd.DataFrame(self.get_matrix(FragmentType.MZ)[0].toarray())
            mz_cols.columns = self._gen_column_names(FragmentType.MZ)  # , set(self.obs["FRAGMENTATION"]))
            df_merged = pd.concat([df_merged, mz_cols], axis=1)
        if "raw_int" in list(self.layers):
            raw_cols = pd.DataFrame(self.get_matrix(FragmentType.RAW)[0].toarray())
            raw_cols.columns = self._gen_column_names(FragmentType.RAW)  # , set(self.obs["FRAGMENTATION"]))
            df_merged = pd.concat([df_merged, raw_cols], axis=1)
        if "pred_int" in list(self.layers):
            pred_cols = pd.DataFrame(self.get_matrix(FragmentType.PRED)[0].toarray())
            pred_cols.columns = self._gen_column_names(FragmentType.PRED)  # , set(self.obs["FRAGMENTATION"]))
            df_merged = pd.concat([df_merged, pred_cols], axis=1)
        return df_merged

    def preprocess_for_machine_learning(
        self,
        filter_peptides: bool = False,
        include_intensities: bool = True,
        include_additional_columns: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """Filter and preprocess for machine learning applications and transform into a Parquet-serializable dataframe.

        :param filter_peptides: Whether to filter peptides based on Andromeda score
        :param include_intensities: Whether to include intensity (label) column
        :param include_additional_columns: Additional column names that are not required by DLomix to include in output.
            Capitalization does not matter - internal column names are all uppercase, whereas returned column names are
            all lowercase.

        :return: Pandas DataFrame with column names and dtypes corresponding to those required by DLomix
            - modified_sequence (str)
            - precursor_charge_onehot (list[int])
            - collision_energy_aligned_normal (int)
            - method_nbr (int)
            [- intensities_raw (list[float]) (if `include_intensities == True`)]
            [additional columns (if specified via `include_additional_columns`)]
        """
        self.remove_decoys()

        if filter_peptides:
            logger.warning("Peptide filtering by score not implemented, including all peptides")
            threshold = 0.0  # TODO actually set threshold
            self.filter_by_score(threshold)

        df = pd.DataFrame()
        df["modified_sequence"] = self.obs["MODIFIED_SEQUENCE"]
        df["precursor_charge_onehot"] = list(np.eye(6, dtype=int)[self.obs["PRECURSOR_CHARGE"].to_numpy() - 1])
        df["collision_energy_aligned_normed"] = self.obs["COLLISION_ENERGY"]
        df["method_nbr"] = self.obs["FRAGMENTATION"].apply(lambda x: c.FRAGMENTATION_ENCODING[x])

        if include_intensities:
            raw_int = self.layers["raw_int"].toarray()
            raw_int[raw_int == 0] = -1
            raw_int[raw_int == c.EPSILON] = 0
            df["intensities_raw"] = list(raw_int)

        if include_additional_columns:
            for column_name in include_additional_columns:
                if column_name.upper() in self.obs:
                    df[column_name.lower()] = self.obs[column_name.upper()]
                else:
                    logger.warning(f"Column '{column_name.upper()}' not present in spectrum, excluded from output")

        return df
