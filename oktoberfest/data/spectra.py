from __future__ import annotations

import logging
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Optional, TypeVar

import anndata
import numpy as np
import pandas as pd
import scipy
import spectrum_fundamentals.constants as c
from scipy.sparse import csr_matrix, dok_matrix
from spectrum_fundamentals.fragments import format_fragment_ion_annotation, generate_fragment_ion_annotations

if TYPE_CHECKING:
    from anndata.compat import Index

logger = logging.getLogger(__name__)


SpectraT = TypeVar("SpectraT", bound="Spectra")


class FragmentType(Enum):
    """FragmentType class to enumerate pred, raw, and mz."""

    PRED = 1
    PRED_A = 2
    PRED_B = 3
    RAW = 4
    RAW_A = 5
    RAW_B = 6
    MZ = 7
    MZ_A = 8
    MZ_B = 9


class Spectra(anndata.AnnData):
    """Main to init spectra data."""

    INTENSITY_COLUMN_PREFIX = "INTENSITY_RAW"
    INTENSITY_COLUMN_PREFIX_A = "INTENSITY_RAW_A"
    INTENSITY_COLUMN_PREFIX_B = "INTENSITY_RAW_B"
    INTENSITY_PRED_PREFIX = "INTENSITY_PRED"
    INTENSITY_PRED_PREFIX_A = "INTENSITY_PRED_A"
    INTENSITY_PRED_PREFIX_B = "INTENSITY_PRED_B"
    MZ_COLUMN_PREFIX = "MZ_RAW"
    MZ_COLUMN_PREFIX_A = "MZ_RAW_A"
    MZ_COLUMN_PREFIX_B = "MZ_RAW_B"
    INTENSITY_PRED_LAYER_NAME = "pred_int"
    INTENSITY_PRED_LAYER_NAME_A = "pred_int_A"
    INTENSITY_PRED_LAYER_NAME_B = "pred_int_B"
    INTENSITY_LAYER_NAME = "raw_int"
    INTENSITY_LAYER_NAME_A = "raw_int_A"
    INTENSITY_LAYER_NAME_B = "raw_int_B"
    MZ_LAYER_NAME = "mz"
    MZ_LAYER_NAME_A = "mz_A"
    MZ_LAYER_NAME_B = "mz_B"
    COLUMNS_FRAGMENT_ION = ["Y1+", "Y1++", "Y1+++", "B1+", "B1++", "B1+++"]
    MAX_CHARGE = 3

    @staticmethod
    def _gen_vars_df(ion_types: list[str] = c.FRAGMENTATION_TO_IONS_BY_PAIRS["HCD"], cms2: bool = False) -> pd.DataFrame:
        """
        Create annotation dataframe for vars in AnnData object.

        :param ion_types: ion types that are expected to be in the spectra
        :param cms2: cleavable crosslinked or linear peptide
        :return: pd.Dataframe of fragment annotations
        """
        df = pd.DataFrame(
            [
                {"ion": f"{ion_type}{pos}+{charge}", "num": pos, "type": ion_type, "charge": charge}
                for pos in (c.POSITIONS_XL if cms2 else c.POSITIONS)
                for ion_type in ion_types
                for charge in c.CHARGES
            ]
        )
        df.set_index("ion", inplace=True)
        return df

    @staticmethod
    def _gen_column_names(fragment_type: FragmentType, cms2: bool = False) -> list[str]:
        """
        Get column names of the spectra data.

        :param fragment_type: choose predicted, raw, or mz
        :param cms2: cleavable crosslinked or linear peptide
        :return: A list of column names
        """
        prefix = Spectra._resolve_prefix(fragment_type)
        columns = []
        if cms2:
            max_range = 59
        else:
            max_range = 30
        for i in range(1, max_range):
            for column in Spectra.COLUMNS_FRAGMENT_ION:
                columns.append(prefix + "_" + column.replace("1", str(i)))
        return columns

    @staticmethod
    def _resolve_prefix(fragment_type: FragmentType) -> str:
        """
        Resolve prefix given fragment type.

        (1 for pred, 2 for xl_pred_a, 3 for xl_pred_a, 4 for raw, 5 for xl_raw_a,
        6 for xl_raw_b, 7 for mz, 8 for xl_mz_a, 9 for xl_mz_b).

        :param fragment_type: choose predicted, raw, or mz
        :return: prefix as string
        """
        logger.debug(fragment_type)
        if fragment_type.value == 1:
            prefix = Spectra.INTENSITY_PRED_PREFIX
        elif fragment_type.value == 2:
            prefix = Spectra.INTENSITY_PRED_PREFIX_A
        elif fragment_type.value == 3:
            prefix = Spectra.INTENSITY_PRED_PREFIX_B
        elif fragment_type.value == 4:
            prefix = Spectra.INTENSITY_COLUMN_PREFIX
        elif fragment_type.value == 5:
            prefix = Spectra.INTENSITY_COLUMN_PREFIX_A
        elif fragment_type.value == 6:
            prefix = Spectra.INTENSITY_COLUMN_PREFIX_B
        elif fragment_type.value == 7:
            prefix = Spectra.MZ_COLUMN_PREFIX
        elif fragment_type.value == 8:
            prefix = Spectra.MZ_COLUMN_PREFIX_A
        else:
            prefix = Spectra.MZ_COLUMN_PREFIX_B
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
            layer = Spectra.INTENSITY_PRED_LAYER_NAME_A
        elif fragment_type.value == 3:
            layer = Spectra.INTENSITY_PRED_LAYER_NAME_B
        elif fragment_type.value == 4:
            layer = Spectra.INTENSITY_LAYER_NAME
        elif fragment_type.value == 5:
            layer = Spectra.INTENSITY_LAYER_NAME_A
        elif fragment_type.value == 6:
            layer = Spectra.INTENSITY_LAYER_NAME_B
        elif fragment_type.value == 7:
            layer = Spectra.MZ_LAYER_NAME
        elif fragment_type.value == 8:
            layer = Spectra.MZ_LAYER_NAME_A
        elif fragment_type.value == 9:
            layer = Spectra.MZ_LAYER_NAME_B
        return layer

    def __getitem__(self, index: Index):
        """Returns a sliced view of the object with this type to avoid returning AnnData instances when slicing."""
        oidx, vidx = self._normalize_indices(index)
        return Spectra(self, oidx=oidx, vidx=vidx, asview=True)

    def add_column(self, data: np.ndarray | pd.Series, name: str | None = None) -> None:
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

    def add_intensities_without_mapping(self, intensities: np.ndarray, fragment_type: FragmentType):
        """
        Add predicted intensities and convert to sparse matrix.

        This function takes a numpy array, containing intensities.
        The intensity array is expected to have the same shape as this object and will be added to
        the respective lazer without checking the order of fragment annotations.

        :param intensities: intensity numpy array to add with shapes (n x m)
        :param fragment_type: the type of intensities to add. Can be FragmentType.RAW or FragmentType.PRED.
        """
        intensities[intensities == 0] = c.EPSILON
        intensities[intensities == -1] = 0.0

        layer = self._resolve_layer_name(fragment_type)
        self.layers[layer] = csr_matrix(intensities)

    def add_list_of_predicted_intensities(
        self,
        intensities: list[np.ndarray],
        annotations: list[np.ndarray],
        chunk_indices: list[np.ndarray],
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

    def get_matrix(self, fragment_type: FragmentType) -> csr_matrix:
        """
        Get intensities sparse matrix from AnnData object.

        :param fragment_type: choose predicted, raw, or mz
        :return: sparse matrix with the required data
        """
        prefix = self._resolve_prefix(fragment_type)
        logger.debug(prefix)

        layer = self._resolve_layer_name(fragment_type)
        matrix = self.layers[layer]

        return matrix

    def write_as_hdf5(self, output_file: str | Path):
        """
        Write spectra_data to hdf5 file.

        :param output_file: path to output file
        """
        self.write(output_file, compression="gzip")

    @classmethod
    def from_hdf5(cls: type[SpectraT], input_file: str | Path) -> SpectraT:
        """
        Read from hdf5 file.

        :param input_file: path to input file
        :return: a spectra instance
        """
        return cls(anndata.read_h5ad(str(input_file)))

    def remove_decoys(self) -> None:
        """Remove decoys in-place."""
        self.__dict__ = Spectra(self[~self.obs.REVERSE].copy()).__dict__

    def filter_by_score(self, threshold: float) -> None:
        """
        Filter out peptides with search engine score below threshold in-place.

        :param threshold: The threshold to use below which peptides are filtered out.
        """
        self.__dict__ = Spectra(self[self.obs.SCORE >= threshold].copy()).__dict__

    def remove_duplicates(self, num_duplicates: int) -> None:
        """Filter out (peptide, charge, collision energy) duplicates if there's more than n_duplicates."""
        self.obs["duplicate_count"] = self.obs.groupby(["SEQUENCE", "PRECURSOR_CHARGE", "COLLISION_ENERGY"]).cumcount()
        self.__dict__ = Spectra(self[self.obs["duplicate_count"] < num_duplicates].copy()).__dict__
        self.obs.drop(columns="duplicate_count", inplace=True)

    def convert_to_df(self) -> pd.DataFrame:
        """
        Gives back spectra_data instance as a pandas Dataframe.

        :return: a pandas DataFrame
        """
        df_merged = self.obs
        logger.debug(self.obs.columns)

        if "mz" in list(self.layers):
            mz_cols = pd.DataFrame(self.get_matrix(FragmentType.MZ).toarray())
            mz_cols.columns = self._gen_column_names(FragmentType.MZ)
            df_merged = pd.concat([df_merged, mz_cols], axis=1)
        if "raw_int" in list(self.layers):
            raw_cols = pd.DataFrame(self.get_matrix(FragmentType.RAW).toarray())
            raw_cols.columns = self._gen_column_names(FragmentType.RAW)
            df_merged = pd.concat([df_merged, raw_cols], axis=1)
        if "pred_int" in list(self.layers):
            pred_cols = pd.DataFrame(self.get_matrix(FragmentType.PRED).toarray())
            pred_cols.columns = self._gen_column_names(FragmentType.PRED)
            df_merged = pd.concat([df_merged, pred_cols], axis=1)
        return df_merged

    def preprocess_for_machine_learning(
        self,
        include_intensities: bool = True,
        include_additional_columns: Optional[list[str]] = None,
        ion_type_order: Optional[list[str]] = None,
        remove_decoys: bool = False,
        search_engine_score_threshold: Optional[float] = None,
        num_duplicates: Optional[int] = None,
    ) -> pd.DataFrame:
        """Filter and preprocess for machine learning applications and transform into a Parquet-serializable dataframe.

        :param include_intensities: Whether to include intensity (label) column
        :param include_additional_columns: Additional column names that are not required by DLomix to include in output.
            Capitalization does not matter - internal column names are all uppercase, whereas returned column names are
            all lowercase.
        :param ion_type_order: Ion type order in which to save output intensity values.
        :param remove_decoys: Whether to remove decoys
        :param search_engine_score_threshold: Search engine score cutoff for peptides included in output
        :param num_duplicates: Number of (sequence, charge, collision energy) duplicates to keep in output

        :return: Pandas DataFrame with column names and dtypes corresponding to those required by DLomix
            - modified_sequence (str)
            - precursor_charge_onehot (list[int])
            - collision_energy_aligned_normal (int)
            - method_nbr (int)
            [- intensities_raw (list[float]) (if `include_intensities == True`)]
            [additional columns (if specified via `include_additional_columns`)]
        """
        if remove_decoys:
            self.remove_decoys()

        if search_engine_score_threshold:
            self.filter_by_score(search_engine_score_threshold)

        if num_duplicates:
            self.remove_duplicates(num_duplicates)

        df = pd.DataFrame()
        df["modified_sequence"] = self.obs["MODIFIED_SEQUENCE"]
        df["precursor_charge_onehot"] = list(
            np.eye(c.NUM_CHARGES_ONEHOT, dtype=int)[self.obs["PRECURSOR_CHARGE"].to_numpy() - 1]
        )
        df["collision_energy_aligned_normed"] = self.obs["COLLISION_ENERGY"]
        df["method_nbr"] = self.obs["FRAGMENTATION"].apply(lambda x: c.FRAGMENTATION_ENCODING[x])

        if include_intensities:
            intensities = self.to_df(layer=self._resolve_layer_name(FragmentType.RAW))
            intensities[intensities == 0] = -1
            intensities[intensities == c.EPSILON] = 0

            if ion_type_order:
                fragment_ion_order = [
                    format_fragment_ion_annotation(ann)
                    for ann in generate_fragment_ion_annotations(
                        ion_type_order, order=("position", "ion_type", "charge")
                    )
                ]
                intensities = intensities[fragment_ion_order]

            df["intensities_raw"] = list(intensities.to_numpy())

        if include_additional_columns:
            for column_name in include_additional_columns:
                if column_name.upper() in self.obs:
                    df[column_name.lower()] = self.obs[column_name.upper()]
                else:
                    logger.warning(f"Column {column_name.upper()!r} not present in spectrum, excluded from output")

        return df
