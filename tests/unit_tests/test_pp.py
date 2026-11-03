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

    def test_calibrate_fragment_masses(self):
        """Test the fragment mass offset is removed and the tolerance follows the spread of the mass errors."""
        from spectrum_fundamentals.fragments import fragment_arrays

        rng = np.random.default_rng(0)
        sequences = ["PEPTIDEK", "LGEYGFQNALIVR", "VEEAQAIASK"] * 20
        mzs = []
        for sequence in sequences:
            fragments = fragment_arrays(sequence, "FTMS", 2, featured_ions=["y", "b"])[0]
            theoretical = fragments["mass"][(fragments["charge"] == 1) & ~fragments["neutral_loss"]]
            ppm = 4 + rng.normal(0, 1, len(theoretical))  # offset 4 ppm, SD 1 ppm
            mzs.append(np.sort(theoretical * (1 + ppm * 1e-6)))
        psms = pd.DataFrame(
            {
                "MZ": mzs,
                "MODIFIED_SEQUENCE": sequences,
                "PRECURSOR_CHARGE": 2,
                "MASS_ANALYZER": "FTMS",
                "REVERSE": False,
            }
        )

        calibrated, tolerance = pp.calibrate_fragment_masses(psms, featured_ions=["y", "b"])

        self.assertAlmostEqual(tolerance, 1.645, delta=0.3)  # 90th percentile of |N(0, 1)|
        np.testing.assert_allclose(calibrated["MZ"][0], mzs[0] * (1 - 4e-6), rtol=1e-6)
        self.assertIs(psms["MZ"][0], mzs[0])  # the input is not changed

        itms = psms.assign(MASS_ANALYZER="ITMS")
        self.assertEqual(pp.calibrate_fragment_masses(itms, featured_ions=["y", "b"]), (itms, None))

    def test_convert_search_custom_input(self):
        """Test that search results in internal format are loaded as they are when using custom_input."""
        search_results = pd.DataFrame(
            {
                "MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ", "TAIASPEK"],
                "SEQUENCE": ["AAACRFVQ", "TAIASPEK"],
                "PEPTIDE_LENGTH": [8, 8],
                "PRECURSOR_CHARGE": [1, 2],
            }
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            search_file = Path(tmp_dir) / "search_results.csv"
            search_results.to_csv(search_file, index=False)
            converted = pp.convert_search(search_file, "custom_input")
            pd.testing.assert_frame_equal(converted, pp.load_search(search_file))

    def test_list_spectra(self):
        """Test listing of spectra with expected user input."""
        spectra_path = Path(__file__).parent
        spectra_file = spectra_path / "test.mzml"
        spectra_file.open("w").close()
        self.assertEqual([spectra_path / "test.mzml"], pp.list_spectra(spectra_path, input_format="mzml"))
        spectra_file.unlink()

    def test_merge_spectra_and_peptides_numeric_raw_files(self):
        """Test raw files named by numbers only match whether they are read as int or str."""
        search = pd.DataFrame({"RAW_FILE": [1, 2], "SCAN_NUMBER": [10, 20], "MODIFIED_SEQUENCE": ["PEPK", "PEPR"]})
        spectra = pd.DataFrame({"RAW_FILE": ["1", "2"], "SCAN_NUMBER": [10, 20], "RETENTION_TIME": [1.0, 2.0]})
        psms = pp.merge_spectra_and_peptides(spectra=spectra, search=search)
        self.assertEqual(psms["MODIFIED_SEQUENCE"].tolist(), ["PEPK", "PEPR"])
        self.assertEqual(psms["RAW_FILE"].tolist(), ["1", "2"])

    def test_merge_spectra_and_peptides_named_raw_files(self):
        """Test raw files with regular names are merged unchanged."""
        search = pd.DataFrame({"RAW_FILE": ["a", "b"], "SCAN_NUMBER": [10, 20], "MODIFIED_SEQUENCE": ["PEPK", "PEPR"]})
        spectra = pd.DataFrame({"RAW_FILE": ["b", "a"], "SCAN_NUMBER": [20, 10], "RETENTION_TIME": [2.0, 1.0]})
        psms = pp.merge_spectra_and_peptides(spectra=spectra, search=search)
        self.assertEqual(psms["RETENTION_TIME"].tolist(), [1.0, 2.0])
        self.assertEqual(psms["RAW_FILE"].tolist(), ["a", "b"])

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

    def test_add_intensities_with_arrow_backed_var_names(self):
        """Test that add_intensities works when var_names is backed by ArrowStringArray (pandas 3 default)."""
        var_df = Spectra._gen_vars_df()
        obs = pd.DataFrame(
            {
                "RAW_FILE": ["file_a"],
                "SCAN_NUMBER": [1],
                "MODIFIED_SEQUENCE": ["PEPTIDEK"],
                "PRECURSOR_CHARGE": [2],
                "COLLISION_ENERGY": [30.0],
                "FRAGMENTATION": ["HCD"],
            }
        )
        spectra = Spectra(obs=obs, var=var_df)
        # Force Arrow-backed string index to simulate pandas 3 default behavior
        arrow_index = pd.Index(pd.array(var_df.index.tolist(), dtype="string[pyarrow]"))
        spectra.var_names = arrow_index

        intensities = np.ones((1, len(var_df)))
        annotation = np.array(arrow_index)[None, ...]
        # Verify the annotation array has the correct shape and the mapping works
        self.assertEqual(annotation.shape, (1, len(var_df)))
        # add_intensities must not raise IndexError with Arrow-backed var_names
        spectra.add_intensities(intensities, annotation, FragmentType.RAW)
        self.assertIn("raw_int", spectra.layers)


class TestSplitSearch(unittest.TestCase):
    """Test splitting search results and timsTOF metadata by spectra file."""

    def test_file_names_with_dots(self):
        """Test spectra file names with dots keep their full name (the runner reads <spectra file stem>.rescore)."""
        search_results = pd.DataFrame({"RAW_FILE": ["run.trimmed", "run.trimmed", "other"], "SCAN_NUMBER": [1, 2, 3]})
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir)
            found = pp.split_search(search_results, output_dir, filenames=["run.trimmed", "missing"])
            found_meta = pp.split_timstof_metadata(search_results, output_dir, filenames=["run.trimmed"])

            self.assertListEqual(found, ["run.trimmed"])
            self.assertListEqual(found_meta, ["run.trimmed"])
            spectra_file = Path("run.trimmed.mzML")
            rescore = pd.read_csv(output_dir / spectra_file.with_suffix(".rescore").name, index_col=0)
            self.assertListEqual(rescore["SCAN_NUMBER"].tolist(), [1, 2])
            self.assertTrue((output_dir / spectra_file.with_suffix(".timsmeta").name).exists())
            self.assertListEqual(
                sorted(p.name for p in output_dir.iterdir()), ["run.trimmed.rescore", "run.trimmed.timsmeta"]
            )
