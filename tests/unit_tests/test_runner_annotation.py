import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from oktoberfest.data import Spectra
from oktoberfest.runner import _annotate_and_get_library
from oktoberfest.utils import Config


def _make_annotated_library(n: int = 2, instrument_type: str = "QE") -> Spectra:
    obs = pd.DataFrame(
        {
            "SCAN_NUMBER": list(range(1, n + 1)),
            "INSTRUMENT_TYPES": [instrument_type] * n,
            "start_pos_p1": [1] * n,
            "start_pos_p2": [2] * n,
        }
    )
    obs.index = obs.index.astype(str)
    return Spectra(obs=obs, var=Spectra._gen_vars_df())


@patch("oktoberfest.runner.pp.annotate_spectral_library_xl")
@patch("oktoberfest.runner.pp.annotate_spectral_library")
@patch("oktoberfest.runner.pp.merge_spectra_and_peptides")
@patch("oktoberfest.runner.pp.load_search")
@patch("oktoberfest.runner.pp.load_spectra")
@patch("oktoberfest.runner.pp.convert_d_to_hdf")
@patch("oktoberfest.runner.pp.convert_raw_to_mzml")
class TestAnnotateAndGetLibrary(unittest.TestCase):
    """Test loading, merging and annotating spectra with the search results."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.config = Config()
        self.config.base_path = Path(self.tmp_dir.name)
        self.config.data["output"] = "out"
        self.config.data["inputs"] = {"search_results_type": "maxquant"}
        self.config.data["models"] = {"intensity": "Prosit_2020_intensity_HCD"}
        self.config.output.mkdir()
        self.hdf5_path = self.config.output / "data" / "run_a.mzml.hdf5"

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    def _set_up_mocks(self, mock_load_spectra, mock_annotate, fragmentation=("HCD", "CID")):
        mock_load_spectra.return_value = pd.DataFrame(
            {"SCAN_NUMBER": [1, 2], "FRAGMENTATION": list(fragmentation), "INSTRUMENT_TYPES": ["QE", "QE"]}
        )
        mock_annotate.return_value = _make_annotated_library()

    def _loaded_spectra(self, mock_merge) -> pd.DataFrame:
        """Return the spectra that were merged with the search results."""
        return mock_merge.call_args.args[0]

    def test_keeps_fragmentation_of_spectra_if_not_configured(
        self, mock_raw, mock_d, mock_load_spectra, mock_load_search, mock_merge, mock_annotate, mock_annotate_xl
    ):
        """Test CID spectra are not relabeled as HCD without fragmentation_method in the config."""
        self._set_up_mocks(mock_load_spectra, mock_annotate)

        _annotate_and_get_library(Path("run_a.mzML"), self.config)

        self.assertEqual(self._loaded_spectra(mock_merge)["FRAGMENTATION"].tolist(), ["HCD", "CID"])

    def test_annotates_with_fragmentation_of_spectra_if_not_configured(
        self, mock_raw, mock_d, mock_load_spectra, mock_load_search, mock_merge, mock_annotate, mock_annotate_xl
    ):
        """Test the fragmentation shared by all spectra is used for annotation without fragmentation_method."""
        for fragmentation, expected in [("CID", "CID"), ("HCD", "HCD"), ("uvpd", "UVPD")]:
            with self.subTest(fragmentation=fragmentation):
                self.hdf5_path.unlink(missing_ok=True)
                self._set_up_mocks(mock_load_spectra, mock_annotate, fragmentation=(fragmentation, fragmentation))
                mock_merge.side_effect = lambda spectra, search: spectra

                with self.assertNoLogs("oktoberfest.runner", level="WARNING"):
                    _annotate_and_get_library(Path("run_a.mzML"), self.config)

                self.assertEqual(mock_annotate.call_args.kwargs["fragmentation_method"], expected)

    def test_unknown_unsupported_or_mixed_fragmentation_falls_back_to_hcd(
        self, mock_raw, mock_d, mock_load_spectra, mock_load_search, mock_merge, mock_annotate, mock_annotate_xl
    ):
        """
        Test fragment ions are annotated as HCD with a warning if the spectra do not share one supported method.

        The fragmentation of the spectra is kept as is. Downstream, non HCD spectra are excluded from CE calibration,
        unknown ones get the UnknownFragmentationMethod percolator feature and are rejected by DLOmix models using
        method_onehot.
        """
        for fragmentation in [("HCD", "unknown"), ("unknown", "unknown"), ("HCD", "CID"), ("ETD", "ETD")]:
            with self.subTest(fragmentation=fragmentation):
                self.hdf5_path.unlink(missing_ok=True)
                self._set_up_mocks(mock_load_spectra, mock_annotate, fragmentation=fragmentation)
                mock_merge.side_effect = lambda spectra, search: spectra

                with self.assertLogs("oktoberfest.runner", level="WARNING") as logs:
                    _annotate_and_get_library(Path("run_a.mzML"), self.config)

                self.assertIn("annotated as HCD", logs.output[0])
                self.assertEqual(self._loaded_spectra(mock_merge)["FRAGMENTATION"].tolist(), list(fragmentation))
                self.assertEqual(mock_annotate.call_args.kwargs["fragmentation_method"], "HCD")

    def test_configured_fragmentation_overrides_spectra(
        self, mock_raw, mock_d, mock_load_spectra, mock_load_search, mock_merge, mock_annotate, mock_annotate_xl
    ):
        """Test fragmentation_method in the config overrides the fragmentation of all spectra."""
        self._set_up_mocks(mock_load_spectra, mock_annotate, fragmentation=("HCD", "unknown"))
        self.config.data["fragmentation_method"] = "CID"

        _annotate_and_get_library(Path("run_a.mzML"), self.config)

        self.assertEqual(self._loaded_spectra(mock_merge)["FRAGMENTATION"].tolist(), ["CID", "CID"])
        self.assertEqual(mock_annotate.call_args.kwargs["fragmentation_method"], "CID")

    def test_configured_instrument_type_overrides_spectra(
        self, mock_raw, mock_d, mock_load_spectra, mock_load_search, mock_merge, mock_annotate, mock_annotate_xl
    ):
        """Test instrument_type in the config overrides the instrument type read from the spectra."""
        self._set_up_mocks(mock_load_spectra, mock_annotate)
        self.config.data["inputs"]["instrument_type"] = "lumos"

        _annotate_and_get_library(Path("run_a.mzML"), self.config)

        self.assertEqual(self._loaded_spectra(mock_merge)["INSTRUMENT_TYPES"].tolist(), ["LUMOS", "LUMOS"])

    def test_spectra_formats(
        self, mock_raw, mock_d, mock_load_spectra, mock_load_search, mock_merge, mock_annotate, mock_annotate_xl
    ):
        """Test RAW files are converted to mzML, .d folders to hdf and mzML files are loaded directly."""
        spectra_dir = self.config.output / "spectra"
        for spectra_file, converter, expected_file in [
            (Path("run_a.raw"), mock_raw, spectra_dir / "run_a.mzML"),
            (Path("run_a.d"), mock_d, spectra_dir / "run_a.hdf"),
            (Path("run_a.mzML"), None, Path("run_a.mzML")),
        ]:
            with self.subTest(spectra_file=spectra_file):
                for mock in (mock_raw, mock_d, mock_load_spectra):
                    mock.reset_mock()
                self._set_up_mocks(mock_load_spectra, mock_annotate)
                self.hdf5_path.unlink(missing_ok=True)

                _annotate_and_get_library(spectra_file, self.config)

                for mock in (mock_raw, mock_d):
                    if mock is converter:
                        self.assertEqual(mock.call_args.args[:2], (spectra_file, expected_file))
                    else:
                        mock.assert_not_called()
                self.assertEqual(mock_load_spectra.call_args.args[0], expected_file)
                mock_load_search.assert_called_with(self.config.output / "msms" / "run_a.rescore")
                self.assertTrue(self.hdf5_path.is_file())

    def test_raw_conversion_keeps_ms1_for_ms1_quantification(
        self, mock_raw, mock_d, mock_load_spectra, mock_load_search, mock_merge, mock_annotate, mock_annotate_xl
    ):
        """Test RAW files keep their MS1 scans only if precursors are quantified on MS1."""
        for data, ms_level in [
            ({}, 2),
            ({"quantification": True}, 2),
            ({"quantification": True, "quantificationOptions": {"precursor": "flashlfq"}}, [1, 2]),
        ]:
            with self.subTest(data=data):
                self.config.data.update(data)
                self._set_up_mocks(mock_load_spectra, mock_annotate)
                self.hdf5_path.unlink(missing_ok=True)

                _annotate_and_get_library(Path("run_a.raw"), self.config)

                self.assertEqual(mock_raw.call_args.kwargs["ms_level"], ms_level)

    def test_stored_library_is_reused(
        self, mock_raw, mock_d, mock_load_spectra, mock_load_search, mock_merge, mock_annotate, mock_annotate_xl
    ):
        """Test an existing annotated library is loaded instead of annotating again."""
        self.hdf5_path.parent.mkdir()
        _make_annotated_library(3).write_as_hdf5(self.hdf5_path)

        library = _annotate_and_get_library(Path("run_a.mzML"), self.config)

        self.assertEqual(library.n_obs, 3)
        for mock in (mock_raw, mock_load_spectra, mock_merge, mock_annotate, mock_annotate_xl):
            mock.assert_not_called()

    def test_stored_library_gets_configured_instrument_type(
        self, mock_raw, mock_d, mock_load_spectra, mock_load_search, mock_merge, mock_annotate, mock_annotate_xl
    ):
        """Test a changed instrument_type in the config is applied to and stored in an existing library."""
        self.hdf5_path.parent.mkdir()
        _make_annotated_library(instrument_type="QE").write_as_hdf5(self.hdf5_path)
        self.config.data["inputs"]["instrument_type"] = "LUMOS"

        library = _annotate_and_get_library(Path("run_a.mzML"), self.config)

        self.assertTrue((library.obs["INSTRUMENT_TYPES"] == "LUMOS").all())
        self.assertTrue((Spectra.from_hdf5(self.hdf5_path).obs["INSTRUMENT_TYPES"] == "LUMOS").all())
        mock_annotate.assert_not_called()

    def test_crosslinking_models_use_xl_annotation(
        self, mock_raw, mock_d, mock_load_spectra, mock_load_search, mock_merge, mock_annotate, mock_annotate_xl
    ):
        """Test XL models use the crosslinking annotation and xiSEARCH start positions are stored as strings."""
        self._set_up_mocks(mock_load_spectra, mock_annotate_xl)
        self.config.data["inputs"]["search_results_type"] = "xisearch"
        for model, cms2 in [("Prosit_2024_intensity_XL_NMS2", False), ("Prosit_2023_intensity_XL_CMS2", True)]:
            with self.subTest(model=model):
                self.config.data["models"]["intensity"] = model
                self.hdf5_path.unlink(missing_ok=True)
                mock_annotate_xl.return_value = _make_annotated_library()

                library = _annotate_and_get_library(Path("run_a.mzML"), self.config)

                self.assertEqual(mock_annotate_xl.call_args.kwargs["cms2"], cms2)
                self.assertEqual(library.obs["start_pos_p1"].tolist(), ["1", "1"])
                self.assertEqual(library.obs["start_pos_p2"].tolist(), ["2", "2"])
        mock_annotate.assert_not_called()
