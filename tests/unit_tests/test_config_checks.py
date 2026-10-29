import json
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


class TestQuantificationOptions(unittest.TestCase):
    """Test the quantificationOptions and the checks for MS1 quantification with FlashLFQ."""

    def setUp(self):  # noqa: D102
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.flashlfq = Path(self.tmp_dir.name) / "CMD.dll"
        self.flashlfq.touch()

    def tearDown(self):  # noqa: D102
        self.tmp_dir.cleanup()

    def _flashlfq_config(self, **data) -> Config:
        options = {"precursor": "flashlfq", "flashlfqExe": str(self.flashlfq)}
        config = _config(quantification=True, quantificationOptions=options, **data)
        config.base_path = Path(self.tmp_dir.name)
        config.data["inputs"] = {"spectra_type": "mzml"}
        return config

    def test_defaults(self):
        """Test the search engine's quantification is used by default and MS1 options default to sensible values."""
        config = _config(quantification=True)
        self.assertEqual(config.quant_precursor, "search_engine")
        self.assertEqual(config.quant_protein, "directlfq")
        self.assertEqual(config.quant_normalization, "drift")
        self.assertEqual((config.quant_psm_fdr, config.quant_peptide_fdr, config.quant_protein_fdr), (0.01, 0.01, 0.01))
        self.assertEqual((config.quant_ppm_tolerance, config.quant_ppm_filter), (10, 5))
        self.assertFalse(config.ms1_quantification)
        config._check_quantification_options()

    def test_flashlfq(self):
        """Test a valid FlashLFQ configuration passes and enables MS1 quantification."""
        config = self._flashlfq_config()
        config._check_quantification_options()
        self.assertTrue(config.ms1_quantification)
        self.assertFalse(
            _config(quantification=False, quantificationOptions={"precursor": "flashlfq"}).ms1_quantification
        )

    def test_invalid_options(self):
        """Test unknown methods and invalid q-values or tolerances are rejected."""
        for option, value in [
            ("precursor", "ionquant"),
            ("protein", "maxlfq"),
            ("normalization", "quantile"),
            ("psmFdr", 0),
            ("psmFdr", 1.5),
            ("proteinFdr", "0.01"),
            ("proteinFdr", True),
            ("ppmTolerance", -10),
            ("peptideFdr", 2),
            ("ppmFilter", 0),
            ("ppmFilter", "5"),
        ]:
            with self.subTest(option=option, value=value):
                with self.assertRaisesRegex(ValueError, option):
                    _config(quantificationOptions={option: value})._check_quantification_options()

    def test_ppm_filter_off(self):
        """Test the mass error filter can be switched off with null."""
        config = _config(quantificationOptions={"ppmFilter": None})
        self.assertIsNone(config.quant_ppm_filter)
        config._check_quantification_options()

    def test_flashlfq_requirements(self):
        """Test FlashLFQ needs Percolator, MS1 spectra, no PTM localization and the FlashLFQ executable."""
        with self.assertRaisesRegex(ValueError, "percolator"):
            self._flashlfq_config(fdr_estimation_method="mokapot")._check_quantification_options()
        with self.assertRaisesRegex(ValueError, "ptm_localization"):
            self._flashlfq_config(ptm_localization=True)._check_quantification_options()
        config = self._flashlfq_config()
        config.data["inputs"]["spectra_type"] = "d"
        with self.assertRaisesRegex(ValueError, "spectra_type"):
            config._check_quantification_options()
        self.flashlfq.unlink()
        with self.assertRaisesRegex(FileNotFoundError, "FlashLFQ"):
            self._flashlfq_config()._check_quantification_options()


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

        # MS1 quantification does not need the search engine's quantification files
        with patch.object(Config, "_check_quantification_options") as mock_options:
            _config(
                type="Rescoring", quantification=True, quantificationOptions={"precursor": "flashlfq"}, models=models
            ).check()
        mock_options.assert_called_once()
        mock_quant.assert_called_once()
        self.assertEqual(mock_fasta.call_count, 2)

    def test_dlomix_model(self, mock_koina, mock_speclib, mock_quant, mock_fasta):
        """Test DLOmix models are loaded instead of checking Koina and fail clearly if they cannot be loaded."""
        config = _config(type="Rescoring", models={"dlomix_intensity": "some/path", "irt": "Prosit_2019_irt"})

        with patch.object(Config, "dlomix_pipeline", new_callable=PropertyMock, return_value=MagicMock()):
            config.check()
        mock_koina.assert_not_called()

        with patch.object(Config, "dlomix_pipeline", new_callable=PropertyMock, return_value=None):
            with self.assertRaisesRegex(ValueError, "Could not load DLOmix pipeline from some/path"):
                config.check()


class TestPTMLocalization(unittest.TestCase):
    """Test the PTM settings with and without localization."""

    def test_without_localization(self):
        """Test PTMs are not localized if localization is disabled, even if PTMs are specified."""
        config = _config(ptm_localization=False, ptmLocalizationOptions={"ptms": {"21": ["S", "T", "Y"]}})
        self.assertEqual(config.ptm_sites, {})
        self.assertEqual(config.ptm_unimod_id, 0)
        self.assertEqual(config.ptm_possible_sites, [])
        self.assertFalse(config.ptm_use_neutral_loss)
        with self.assertLogs("oktoberfest.utils.config", level="WARNING"):
            config._check_ptm_localization()

    def test_multiple_ptms(self):
        """Test several PTMs with their own residues can be localized."""
        config = _config(
            ptm_localization=True,
            ptmLocalizationOptions={"ptms": {"21": ["S", "T", "Y"], "7": ["R"]}, "neutral_loss": True},
        )
        config._check_ptm_localization()
        self.assertEqual(config.ptm_sites, {21: ["S", "T", "Y"], 7: ["R"]})
        self.assertEqual(config.ptm_unimod_id, 0)
        self.assertEqual(config.ptm_possible_sites, ["R", "S", "T", "Y"])
        self.assertTrue(config.ptm_use_neutral_loss)

    def test_legacy_single_ptm(self):
        """Test the legacy format with a single unimod id and possible sites."""
        config = _config(
            ptm_localization=True, ptmLocalizationOptions={"unimod_id": 7, "possible_sites": ["R", "N", "Q"]}
        )
        config._check_ptm_localization()
        self.assertEqual(config.ptm_sites, {7: ["R", "N", "Q"]})
        self.assertEqual(config.ptm_unimod_id, 7)

    def test_invalid_localization_settings(self):
        """Test localization requires the user to specify PTMs and valid residues."""
        for options, message in [
            ({}, "no PTM to localize"),
            ({"ptms": {}}, "no PTM to localize"),
            ({"ptms": ["21"]}, "must map unimod ids"),
            ({"ptms": {"21": []}}, "No possible residues"),
            ({"ptms": {"21": ["S", "X"]}}, "Invalid residues"),
        ]:
            with self.subTest(options=options):
                with self.assertRaisesRegex(ValueError, message):
                    _config(ptm_localization=True, ptmLocalizationOptions=options)._check_ptm_localization()

    def test_ptm_decoy4_defaults(self):
        """Test the decoy4 defaults: A/G/L/V decoys (without PTM sites), PEP probabilities, count weights, one group."""
        config = _config(ptm_localization=True, ptmLocalizationOptions={"ptms": {"21": ["S", "T", "Y"]}})
        self.assertEqual(config.ptm_decoy_residues, ["A", "G", "L", "V"])
        self.assertEqual(config.ptm_localization_probability, "pep")
        self.assertEqual(config.ptm_decoy_weighting, "count")
        self.assertEqual(config.ptm_flr_groups, "none")
        self.assertEqual(config.ptm_flr_estimator, "hybrid")
        config._check_ptm_localization()
        # residues the PTM can be located on are no decoys
        config = _config(ptm_localization=True, ptmLocalizationOptions={"ptms": {"1": ["A", "K"]}})
        self.assertEqual(config.ptm_decoy_residues, ["G", "L", "V"])
        # an empty list disables the decoys
        options = {"ptms": {"21": ["S", "T", "Y"]}, "decoy_residues": []}
        self.assertEqual(_config(ptm_localization=True, ptmLocalizationOptions=options).ptm_decoy_residues, [])
        self.assertEqual(_config(ptmLocalizationOptions={"ptms": {"21": ["S"]}}).ptm_decoy_residues, [])

    def test_ptm_flr_options_invalid(self):
        """Test invalid FLR estimation options are rejected."""
        for option, value in (
            ("loc_probability", "max"),
            ("decoy_weighting", "equal"),
            ("flr_groups", "peptide"),
            ("flr_estimator", "model"),
        ):
            with self.subTest(option=option):
                options = {"ptms": {"21": ["S", "T", "Y"]}, option: value}
                with self.assertRaisesRegex(ValueError, option):
                    _config(ptm_localization=True, ptmLocalizationOptions=options)._check_ptm_localization()

    def test_ptm_ambiguity_options(self):
        """Test sdi_features and ambiguity_check are on by default and must be booleans."""
        options = {"ptms": {"21": ["S", "T", "Y"]}}
        config = _config(ptm_localization=True, ptmLocalizationOptions=options)
        self.assertTrue(config.ptm_sdi_features)
        self.assertTrue(config.ptm_ambiguity_check)
        config = _config(
            ptm_localization=True, ptmLocalizationOptions={**options, "sdi_features": False, "ambiguity_check": False}
        )
        config._check_ptm_localization()
        self.assertFalse(config.ptm_sdi_features)
        self.assertFalse(config.ptm_ambiguity_check)
        for option in ("sdi_features", "ambiguity_check"):
            with self.subTest(option=option):
                config = _config(ptm_localization=True, ptmLocalizationOptions={**options, option: "yes"})
                with self.assertRaisesRegex(ValueError, option):
                    config._check_ptm_localization()

    def test_ptm_phospho_loss_residues(self):
        """Test the phosphoric acid loss is annotated on S and T, and on the decoy residues with decoy_neutral_loss."""
        options = {"ptms": {"21": ["S", "T", "Y"]}, "decoy_residues": ["A", "G"]}
        self.assertEqual(_config(ptm_localization=True, ptmLocalizationOptions=options).ptm_phospho_loss_residues, "ST")
        options["decoy_neutral_loss"] = True
        config = _config(ptm_localization=True, ptmLocalizationOptions=options)
        self.assertEqual(config.ptm_phospho_loss_residues, "STAG")

    def test_ptm_decoy_proxies(self):
        """Test decoys are predicted on the first possible site residue, or on the residue given in decoy_proxy."""
        options = {"ptms": {"21": ["S", "T", "Y"]}, "decoy_residues": ["A", "G"]}
        config = _config(ptm_localization=True, ptmLocalizationOptions=options)
        self.assertEqual(config.ptm_decoy_proxies, {"A[UNIMOD:21]": "S[UNIMOD:21]", "G[UNIMOD:21]": "S[UNIMOD:21]"})
        options["decoy_proxy"] = "T"
        config = _config(ptm_localization=True, ptmLocalizationOptions=options)
        self.assertEqual(config.ptm_decoy_proxies, {"A[UNIMOD:21]": "T[UNIMOD:21]", "G[UNIMOD:21]": "T[UNIMOD:21]"})

    def test_ptm_decoy_proxy_search(self):
        """Test decoys are predicted as a mixture of pS and pT with the frequencies of the search engine results."""
        import pandas as pd

        from oktoberfest.predict.predictor import _choose_proxy
        from oktoberfest.rescore import decoy_proxy_fractions

        search = pd.DataFrame(
            {
                "RAW_FILE": ["a", "a", "a", "a", "a"],
                "SCAN_NUMBER": [1, 1, 2, 3, 4],
                "REVERSE": [False, False, False, False, True],
                "SEARCH_MODIFIED_SEQUENCE": ["PES[UNIMOD:21]K"] * 2
                + ["PET[UNIMOD:21]K", "PEY[UNIMOD:21]K", "S[UNIMOD:21]K"],
            }
        )
        # scan 1 counted once, Y is no proxy residue, the decoy PSM (scan 4) is not counted
        fractions = decoy_proxy_fractions(search, {21: ["S", "T", "Y"]})
        self.assertEqual(fractions, {"21": {"S": 0.5, "T": 0.5}})
        with tempfile.TemporaryDirectory() as tmp:
            options = {"ptms": {"21": ["S", "T", "Y"]}, "decoy_residues": ["A"], "decoy_proxy": "search"}
            config = _config(ptm_localization=True, ptmLocalizationOptions=options, output=tmp)
            config.base_path = Path(tmp)
            config.decoy_proxy_fractions_file.parent.mkdir(parents=True)
            config.decoy_proxy_fractions_file.write_text(json.dumps(fractions))
            proxies = config.ptm_decoy_proxies["A[UNIMOD:21]"]
        self.assertEqual(proxies, {"S[UNIMOD:21]": 0.5, "T[UNIMOD:21]": 0.5})
        chosen = [_choose_proxy(f"PEPA[UNIMOD:21]{i}K", proxies) for i in range(2000)]
        self.assertAlmostEqual(chosen.count("T[UNIMOD:21]") / len(chosen), 0.5, delta=0.05)
        self.assertEqual(_choose_proxy("PEA[UNIMOD:21]K", proxies), _choose_proxy("PEA[UNIMOD:21]K", proxies))

    def test_ptm_fdr_control_default(self):
        """Test PTM models split unmodified and modified PSMs by default, but never permute without localization."""
        ptm_models = {"intensity": "Prosit_2024_intensity_PTMs_gl", "irt": "Prosit_2024_irt_PTMs_gl"}
        self.assertTrue(_config(models=ptm_models).ptm_fdr_control)
        self.assertEqual(_config(models=ptm_models).ptm_sites, {})
        self.assertFalse(_config(models=ptm_models, ptm_fdr_control=False).ptm_fdr_control)
        self.assertFalse(_config(models={"intensity": "Prosit_2020_intensity_HCD"}).ptm_fdr_control)
        localization = {"ptms": {"21": ["S", "T", "Y"]}}
        config = _config(models=ptm_models, ptm_localization=True, ptmLocalizationOptions=localization)
        self.assertFalse(config.ptm_fdr_control)
        self.assertEqual(config.ptm_sites, {21: ["S", "T", "Y"]})


class TestPercolatorSubset(unittest.TestCase):
    """Test the maximum number of PSMs percolator trains on."""

    def test_default(self):
        """Test percolator trains on all PSMs by default."""
        config = _config()
        self.assertIsNone(config.percolator_subset_max_train)
        config._check_percolator_subset()

    def test_subset(self):
        """Test a positive number of PSMs is accepted and anything else is rejected."""
        config = _config(percolatorSubsetMaxTrain=500000)
        self.assertEqual(config.percolator_subset_max_train, 500000)
        config._check_percolator_subset()
        for value in [0, -1, 1.5, "500000", True]:
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "percolatorSubsetMaxTrain"):
                    _config(percolatorSubsetMaxTrain=value)._check_percolator_subset()


class TestChimericPeakSubtraction(unittest.TestCase):
    """Test the flag for removing the peaks of better candidates in chimeric spectra."""

    def test_default(self):
        """Test the peaks are not subtracted by default."""
        config = _config()
        self.assertFalse(config.chimeric_peak_subtraction)
        config._check_chimeric_peak_subtraction()

    def test_flag(self):
        """Test a boolean is accepted and anything else is rejected."""
        self.assertTrue(_config(chimericPeakSubtraction=True).chimeric_peak_subtraction)
        for value in [1, "true", None]:
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "chimericPeakSubtraction"):
                    _config(chimericPeakSubtraction=value)._check_chimeric_peak_subtraction()
