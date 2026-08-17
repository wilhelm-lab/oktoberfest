import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import spectrum_fundamentals.constants as c

from oktoberfest import pp
from oktoberfest.data.spectra import FragmentType, Spectra


class TestProcessing(unittest.TestCase):
    """Test class for preprocessing functions."""

    def test_list_spectra(self):
        """Test listing of spectra with expected user input."""
        spectra_path = Path(__file__).parent
        spectra_file = spectra_path / "test.mzml"
        spectra_file.open("w").close()
        self.assertEqual([spectra_path / "test.mzml"], pp.list_spectra(spectra_path, input_format="mzml"))
        spectra_file.unlink()

    def test_list_spectra_with_empty_string_folder(self):
        """Test listing spectra in a string folder without matching files."""
        self.assertRaises(AssertionError, pp.list_spectra, str(Path(__file__).parent), "raw")

    def test_list_spectra_with_wrong_folder(self):
        """Test listing spectra in a folder that does not exist."""
        self.assertRaises(NotADirectoryError, pp.list_spectra, Path(__file__).parent / "noexist", "raw")

    def test_list_spectra_with_wrong_format(self):
        """Test listing spectra with a format that isn't allowed."""
        self.assertRaises(ValueError, pp.list_spectra, Path(__file__).parent, "mzm")

    def test_convert_anndata_to_parquet_with_optional_columns(self):
        """Test converting AnnData hdf5 to parquet when optional columns are available."""
        var_df = Spectra._gen_vars_df()
        obs = pd.DataFrame(
            {
                "RAW_FILE": ["file_a", "file_b"],
                "SCAN_NUMBER": [1, 2],
                "MODIFIED_SEQUENCE": [
                    "LYEADFVLFGYPKPENLLRD",
                    "[UNIMOD:1]-FC[UNIMOD:4]SETIDIIQPALVLATGDLTDAK",
                ],
                "PRECURSOR_CHARGE": [2, 3],
                "COLLISION_ENERGY": [30.0, 32.5],
                "FRAGMENTATION": ["HCD", "CID"],
            }
        )
        spectra = Spectra(obs=obs, var=var_df)
        spectra.var_names = var_df.index
        intensities = np.ones((2, len(var_df)))
        intensities[0, 0] = 0
        intensities[1, 0] = -1
        spectra.add_intensities_without_mapping(intensities, FragmentType.RAW)
        mzs = np.arange(2 * len(var_df), dtype=float).reshape(2, len(var_df))
        spectra.add_mzs(mzs, FragmentType.MZ)

        with tempfile.TemporaryDirectory() as temp_dir:
            hdf5_path = Path(temp_dir) / "spectra.hdf5"
            parquet_path = Path(temp_dir) / "spectra.parquet"
            spectra.write_as_hdf5(hdf5_path)

            pp.convert_anndata_to_parquet(hdf5_path, parquet_path)

            df = pd.read_parquet(parquet_path)

        self.assertEqual(
            list(df.columns),
            [
                "raw_file",
                "scan_number",
                "method_nbr",
                "precursor_charge_onehot",
                "collision_energy_aligned_normed",
                "intensities_raw",
                "mz_raw",
                "package",
                "modified_sequence",
            ],
        )
        self.assertEqual(df["raw_file"].tolist(), ["file_a", "file_b"])
        self.assertEqual(df["scan_number"].dtype, "int64")
        self.assertEqual(df["method_nbr"].tolist(), [c.FRAGMENTATION_ENCODING["HCD"], c.FRAGMENTATION_ENCODING["CID"]])
        self.assertEqual(df["collision_energy_aligned_normed"].dtype, "float64")
        self.assertEqual(df["collision_energy_aligned_normed"].tolist(), [0.3, 0.325])
        self.assertEqual(df["package"].tolist(), ["HCD", "CID"])
        self.assertEqual(
            df["modified_sequence"].tolist(),
            ["[]-LYEADFVLFGYPKPENLLRD-[]", "[UNIMOD:1]-FC[UNIMOD:4]SETIDIIQPALVLATGDLTDAK-[]"],
        )
        np.testing.assert_array_equal(np.asarray(df.loc[0, "precursor_charge_onehot"]), np.array([0, 1, 0, 0, 0, 0]))
        np.testing.assert_array_equal(np.asarray(df.loc[1, "precursor_charge_onehot"]), np.array([0, 0, 1, 0, 0, 0]))
        self.assertEqual(np.asarray(df.loc[0, "intensities_raw"]).shape, (len(var_df),))
        self.assertEqual(np.asarray(df.loc[0, "intensities_raw"])[0], 0)
        self.assertEqual(np.asarray(df.loc[1, "intensities_raw"])[0], -1)
        np.testing.assert_array_equal(np.asarray(df.loc[0, "mz_raw"]), mzs[0])
        np.testing.assert_array_equal(np.asarray(df.loc[1, "mz_raw"]), mzs[1])

    def test_convert_anndata_to_parquet_omits_unavailable_optional_columns(self):
        """Test converting AnnData hdf5 to parquet with only required columns."""
        var_df = Spectra._gen_vars_df()
        spectra = Spectra(
            obs=pd.DataFrame({"MODIFIED_SEQUENCE": ["_FC[UNIMOD:4]SETIDIIQPALVLATGDLTDAK_"]}),
            var=var_df,
        )
        spectra.var_names = var_df.index
        spectra.add_intensities_without_mapping(np.ones((1, len(var_df))), FragmentType.RAW)

        with tempfile.TemporaryDirectory() as temp_dir:
            hdf5_path = Path(temp_dir) / "spectra.hdf5"
            parquet_path = Path(temp_dir) / "spectra.parquet"
            spectra.write_as_hdf5(hdf5_path)

            pp.convert_anndata_to_parquet(hdf5_path, parquet_path)

            df = pd.read_parquet(parquet_path)

        self.assertEqual(list(df.columns), ["intensities_raw", "modified_sequence"])
        self.assertEqual(df["modified_sequence"].tolist(), ["[]-FC[UNIMOD:4]SETIDIIQPALVLATGDLTDAK-[]"])
        np.testing.assert_array_equal(np.asarray(df.loc[0, "intensities_raw"]), np.ones(len(var_df)))

    def test_convert_anndata_to_parquet_keeps_terminal_modifications(self):
        """Test converting modified sequences with explicit terminal modifications."""
        var_df = Spectra._gen_vars_df()
        spectra = Spectra(
            obs=pd.DataFrame({"MODIFIED_SEQUENCE": ["_PEPTIDE-[UNIMOD:2]_", "_[UNIMOD:1]-PEPTIDE-[UNIMOD:2]_"]}),
            var=var_df,
        )
        spectra.var_names = var_df.index
        spectra.add_intensities_without_mapping(np.ones((2, len(var_df))), FragmentType.RAW)

        with tempfile.TemporaryDirectory() as temp_dir:
            hdf5_path = Path(temp_dir) / "spectra.hdf5"
            parquet_path = Path(temp_dir) / "spectra.parquet"
            spectra.write_as_hdf5(hdf5_path)

            pp.convert_anndata_to_parquet(hdf5_path, parquet_path)

            df = pd.read_parquet(parquet_path)

        self.assertEqual(
            df["modified_sequence"].tolist(),
            ["[]-PEPTIDE-[UNIMOD:2]", "[UNIMOD:1]-PEPTIDE-[UNIMOD:2]"],
        )
