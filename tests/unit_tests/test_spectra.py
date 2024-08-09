import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import spectrum_fundamentals.constants as c
from anndata.tests.helpers import assert_equal
from spectrum_fundamentals.fragments import format_fragment_ion_annotation, generate_fragment_ion_annotations

from oktoberfest.data.spectra import Spectra


class TestSpectra(unittest.TestCase):
    """Test class for Spectra."""

    @classmethod
    def setUpClass(cls):  # noqa: D102
        cls.mini_spectra = Spectra.from_hdf5(Path(__file__).parent / "data/spectra/test_spectra.hdf5")
        cls.temp_dir = Path(tempfile.mkdtemp())

        df = pd.read_csv(Path(__file__).parent / "data/spectra/df_for_parquet.csv", sep="\t", index_col="Unnamed: 0")
        df = df.astype({"method_nbr": "category", "modified_sequence": "category"})
        df["intensities_raw"] = df["intensities_raw"].map(
            lambda intens: np.fromstring(
                intens.replace("\n", "").replace("[", "").replace("]", "").replace("  ", " "), sep=" "
            )
        )
        df["precursor_charge_onehot"] = df["precursor_charge_onehot"].map(
            lambda intens: np.fromstring(
                intens.replace("\n", "").replace("[", "").replace("]", "").replace("  ", " "), dtype="int", sep=" "
            )
        )
        cls.df_for_parquet = df

    @classmethod
    def tearDownClass(cls):  # noqa: D102
        shutil.rmtree(cls.temp_dir)

    def test_read_and_write_hdf5(self):
        """Test writing hdf5 files and reading them method."""
        spectra_path = Path(self.temp_dir / "spectra.hdf5")
        spec1 = self.mini_spectra
        spec1.write_as_hdf5(spectra_path)
        spec2 = Spectra.from_hdf5(spectra_path)
        assert_equal(spec1, spec2)

    def test_gen_vardf(self):
        """Test gen_vardf method."""
        for ion_types in [c.HCD_IONS, c.ETD_IONS, c.ETCID_IONS, c.UVPD_IONS]:
            with self.subTest(ion_types=ion_types):
                var_df = Spectra._gen_vars_df(ion_types=ion_types)
                ion_type_annotations = generate_fragment_ion_annotations(
                    ion_types, order=("position", "ion_type", "charge")
                )
                df = pd.DataFrame.from_records(ion_type_annotations, columns=["type", "num", "charge"])[
                    ["num", "type", "charge"]
                ]
                df["ion"] = [format_fragment_ion_annotation(ann) for ann in ion_type_annotations]
                df.set_index("ion", inplace=True)
                pd.testing.assert_frame_equal(var_df, df)

    def test_preprocess_for_machine_learning(self):
        """Test preprocess_for_machine_learning method."""
        df = self.mini_spectra.preprocess_for_machine_learning()
        df = df.astype({"modified_sequence": "object"})
        df = df.astype({"modified_sequence": "category"})
        pd.testing.assert_frame_equal(
            df.reset_index(drop=True), self.df_for_parquet.reset_index(drop=True), check_column_type=False
        )
