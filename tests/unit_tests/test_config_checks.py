import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, PropertyMock, patch

from oktoberfest.utils import Config


def _config(**data) -> Config:
    config = Config()
    config.data.update(data)
    return config


class TestCheckTMT(unittest.TestCase):
    """Test the compatibility check of mass tags and models."""

    def test_tmt_model_without_tag(self):
        """Test TMT models require a tag."""
        for models in [
            {"intensity": "Prosit_2020_intensity_TMT", "irt": "Prosit_2019_irt"},
            {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2020_irt_TMT"},
        ]:
            with self.subTest(models=models):
                with self.assertRaisesRegex(AssertionError, "provided no tag"):
                    _config(models=models)._check_tmt()

    def test_tag_with_incompatible_model(self):
        """Test a tag requires TMT compatible models."""
        for models in [
            {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2020_irt_TMT"},
            {"intensity": "Prosit_2020_intensity_TMT", "irt": "Prosit_2019_irt"},
        ]:
            with self.subTest(models=models):
                with self.assertRaisesRegex(AssertionError, "tmt"):
                    _config(tag="tmt", models=models)._check_tmt()

    def test_compatible_combinations(self):
        """Test valid combinations of tags and models pass."""
        for tag, models in [
            ("", {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"}),
            ("tmt", {"intensity": "Prosit_2020_intensity_TMT", "irt": "Prosit_2020_irt_TMT"}),
            ("tmt", {"intensity": "Prosit_2023_intensity_PTMs_gl", "irt": "Prosit_2023_irt_PTMs_gl"}),
            ("tmt", {"intensity": "AlphaPept_ms2_generic", "irt": "AlphaPept_rt_generic"}),
            ("tmt", {"dlomix_intensity": "some/path", "irt": "Prosit_2019_irt"}),
        ]:
            with self.subTest(tag=tag, models=models):
                _config(tag=tag, models=models)._check_tmt()


class TestCheckForSpeclib(unittest.TestCase):
    """Test the compatibility check of fragmentation and model for spectral library generation."""

    def test_fragmentation_and_model_mismatch(self):
        """Test the configured fragmentation must match the model."""
        for fragmentation, model in [("HCD", "Prosit_2020_intensity_CID"), ("CID", "Prosit_2020_intensity_HCD")]:
            with self.subTest(fragmentation=fragmentation):
                config = _config(models={"intensity": model}, spectralLibraryOptions={"fragmentation": fragmentation})
                with self.assertRaisesRegex(AssertionError, "fragmentation method and the chosen model agree"):
                    config._check_for_speclib()

    def test_missing_fragmentation(self):
        """Test the fragmentation is derived from HCD / CID models and required for all other models."""
        for model in ["Prosit_2020_intensity_HCD", "Prosit_2020_intensity_CID"]:
            with self.subTest(model=model):
                with self.assertLogs("oktoberfest.utils.config", level="WARNING") as logs:
                    _config(models={"intensity": model})._check_for_speclib()
                self.assertIn("No fragmentation method was specified", logs.output[0])

        with self.assertRaisesRegex(AssertionError, "need to provide the fragmentation method"):
            _config(models={"intensity": "Prosit_2023_intensity_timsTOF"})._check_for_speclib()

    def test_matching_fragmentation(self):
        """Test matching fragmentation and model pass."""
        _config(
            models={"intensity": "Prosit_2020_intensity_HCD"}, spectralLibraryOptions={"fragmentation": "HCD"}
        )._check_for_speclib()


class TestCheckQuantification(unittest.TestCase):
    """Test the input file checks for quantification with picked group FDR."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.search_dir = Path(self.tmp_dir.name)

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    def _config(self, search_results_type: str) -> Config:
        config = Config()
        config.base_path = self.search_dir
        config.data["inputs"] = {"search_results": str(self.search_dir), "search_results_type": search_results_type}
        return config

    def test_maxquant(self):
        """Test MaxQuant requires evidence.txt."""
        with self.assertRaisesRegex(AssertionError, "evidence.txt"):
            self._config("maxquant")._check_quantification()
        (self.search_dir / "evidence.txt").touch()
        self._config("maxquant")._check_quantification()

    def test_sage(self):
        """Test Sage requires results.sage.tsv and lfq.tsv."""
        with self.assertRaisesRegex(AssertionError, "results.sage.tsv"):
            self._config("sage")._check_quantification()
        (self.search_dir / "results.sage.tsv").touch()
        with self.assertRaisesRegex(AssertionError, "lfq.tsv"):
            self._config("sage")._check_quantification()
        (self.search_dir / "lfq.tsv").touch()
        self._config("sage")._check_quantification()

    def test_msfragger(self):
        """Test MSFragger requires psm.tsv files in subdirectories and combined_ion.tsv."""
        with self.assertRaisesRegex(AssertionError, "psm.tsv"):
            self._config("msfragger")._check_quantification()
        (self.search_dir / "experiment1").mkdir()
        (self.search_dir / "experiment1" / "psm.tsv").touch()
        with self.assertRaisesRegex(AssertionError, "combined_ion.tsv"):
            self._config("msfragger")._check_quantification()
        (self.search_dir / "combined_ion.tsv").touch()
        self._config("msfragger")._check_quantification()

    def test_search_results_file(self):
        """Test the files are searched next to the search results if a file is given."""
        (self.search_dir / "msms.txt").touch()
        (self.search_dir / "evidence.txt").touch()
        config = self._config("maxquant")
        config.data["inputs"]["search_results"] = str(self.search_dir / "msms.txt")
        config._check_quantification()

    def test_fasta_required(self):
        """Test quantification requires a fasta file as library input."""
        config = self._config("maxquant")
        config.data["inputs"]["library_input_type"] = "peptides"
        with self.assertRaisesRegex(AssertionError, "fasta"):
            config._check_fasta()
        config.data["inputs"]["library_input_type"] = "fasta"
        config._check_fasta()


@patch.object(Config, "_check_fasta")
@patch.object(Config, "_check_quantification")
@patch.object(Config, "_check_for_speclib")
@patch.object(Config, "_check_koina_model_availability")
class TestCheck(unittest.TestCase):
    """Test which validations are run by Config.check."""

    def test_koina_model(self, mock_koina, mock_speclib, mock_quant, mock_fasta):
        """Test Koina models are checked for availability and job specific checks only run when needed."""
        _config(type="Rescoring", models={"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"}).check()

        mock_koina.assert_called_once()
        mock_speclib.assert_not_called()
        mock_quant.assert_not_called()
        mock_fasta.assert_not_called()

    def test_speclib_and_quantification_checks(self, mock_koina, mock_speclib, mock_quant, mock_fasta):
        """Test spectral library and quantification checks run for the respective jobs."""
        models = {"intensity": "Prosit_2020_intensity_HCD", "irt": "Prosit_2019_irt"}
        _config(type="SpectralLibraryGeneration", models=models).check()
        mock_speclib.assert_called_once()

        _config(type="Rescoring", quantification=True, models=models).check()
        mock_quant.assert_called_once()
        mock_fasta.assert_called_once()

    def test_dlomix_model(self, mock_koina, mock_speclib, mock_quant, mock_fasta):
        """Test DLOmix models are loaded instead of checking Koina and fail clearly if they cannot be loaded."""
        config = _config(type="Rescoring", models={"dlomix_intensity": "some/path", "irt": "Prosit_2019_irt"})

        with patch.object(Config, "dlomix_pipeline", new_callable=PropertyMock, return_value=MagicMock()):
            config.check()
        mock_koina.assert_not_called()

        with patch.object(Config, "dlomix_pipeline", new_callable=PropertyMock, return_value=None):
            with self.assertRaisesRegex(ValueError, "Could not load DLOmix pipeline from some/path"):
                config.check()
