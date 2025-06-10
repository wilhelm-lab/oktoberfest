import shutil
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from pandas.testing import assert_frame_equal

import oktoberfest as ok
from oktoberfest.__main__ import main
from oktoberfest.runner import _calculate_features, _ce_calib, _preprocess, prepare_rescore_xl_psm_level
from oktoberfest.utils import Config


class TestRunner(unittest.TestCase):
    """Test class for use cases of runner."""

    def test_speclib_digest(self):
        """Test the runner for a spectral library generation with a fasta digest."""
        config_path = Path(__file__).parent / "configs" / "spectral_library_with_digest.json"
        with patch("sys.argv", ["oktoberfest", f"--config_path={config_path}"]):
            main()

        config = Config()
        config.read(config_path)
        shutil.rmtree(config.output)

    def test_rescoring_cms2_xl(self):
        """Test the runner for a rescoring run with cleavable crosslinking."""
        config_path = Path(__file__).parent / "configs" / "rescoring_cleavable_xl.json"
        # with patch("sys.argv", ["oktoberfest", f"--config_path={config_path}"]):
        config = Config()
        config.read(config_path)
        config.check()

        # load spectra file names
        spectra_files = ok.pp.list_spectra(input_dir=config.spectra, input_format=config.spectra_type)

        proc_dir = config.output / "proc"
        proc_dir.mkdir(parents=True, exist_ok=True)

        spectra_file = _preprocess(spectra_files, config)[0]
        _ = _ce_calib(spectra_file, config)

        _calculate_features(spectra_file, config, xl=True)

        # prepare rescoring

        fdr_dir = config.output / "results" / config.fdr_estimation_method
        original_tab_files = [
            fdr_dir / spectra_file.with_suffix(".original.tab").name for spectra_file in spectra_files
        ]
        rescore_tab_files = [fdr_dir / spectra_file.with_suffix(".rescore.tab").name for spectra_file in spectra_files]

        ok.re.merge_input(tab_files=original_tab_files, output_file=fdr_dir / "original.tab")
        ok.re.merge_input(tab_files=rescore_tab_files, output_file=fdr_dir / "rescore.tab")

        # prepare rescoring file

        rescore_features_path = fdr_dir / "rescore_features_csm.tab"
        if not rescore_features_path.exists():
            shutil.copy(fdr_dir / "rescore.tab", rescore_features_path)
            input_psm_rescore = prepare_rescore_xl_psm_level(str(fdr_dir), "rescore")
            input_psm_rescore.to_csv(str(fdr_dir) + "/rescore.tab", sep="\t", index=None)

        original_features_path = fdr_dir / "original_features_csm.tab"
        if not original_features_path.exists():
            shutil.copy(fdr_dir / "original.tab", original_features_path)
            input_psm_original = prepare_rescore_xl_psm_level(str(fdr_dir), "original")
            input_psm_original.to_csv(str(fdr_dir) + "/original.tab", sep="\t", index=None)

        expected_perc_tab_file = pd.read_csv(
            Path(__file__).parent / "data" / "xl" / "cleavable" / "expected_outputs" / "expected_rescore.tab", sep="\t"
        )

        created_perc_tab_file = pd.read_csv(
            fdr_dir / "rescore.tab",
            sep="\t",
        )

        try:
            assert_frame_equal(
                expected_perc_tab_file, created_perc_tab_file, check_dtype=True, check_exact=False, rtol=1e-1
            )
        except AssertionError as e:
            print("DataFrames are not equal:", e)
            raise  # Re-raise the assertion error for the test framework to catch

        config = Config()
        config.read(config_path)
        shutil.rmtree(Path(__file__).parent / "data" / "xl" / "cleavable" / "out")

    def test_rescoring_nms2_xl(self):
        """Test the runner for a rescoring run with non-cleavable crosslinking."""
        config_path = Path(__file__).parent / "configs" / "rescoring_non_cleavable_xl.json"
        # with patch("sys.argv", ["oktoberfest", f"--config_path={config_path}"]):
        config = Config()
        config.read(config_path)
        config.check()

        # load spectra file names
        spectra_files = ok.pp.list_spectra(input_dir=config.spectra, input_format=config.spectra_type)

        proc_dir = config.output / "proc"
        proc_dir.mkdir(parents=True, exist_ok=True)

        spectra_file = _preprocess(spectra_files, config)[0]
        _ = _ce_calib(spectra_file, config)

        _calculate_features(spectra_file, config, xl=True)

        # prepare rescoring

        fdr_dir = config.output / "results" / config.fdr_estimation_method
        original_tab_files = [
            fdr_dir / spectra_file.with_suffix(".original.tab").name for spectra_file in spectra_files
        ]
        rescore_tab_files = [fdr_dir / spectra_file.with_suffix(".rescore.tab").name for spectra_file in spectra_files]

        ok.re.merge_input(tab_files=original_tab_files, output_file=fdr_dir / "original.tab")
        ok.re.merge_input(tab_files=rescore_tab_files, output_file=fdr_dir / "rescore.tab")

        # prepare rescoring file

        rescore_features_path = fdr_dir / "rescore_features_csm.tab"
        if not rescore_features_path.exists():
            shutil.copy(fdr_dir / "rescore.tab", rescore_features_path)
            input_psm_rescore = prepare_rescore_xl_psm_level(str(fdr_dir), "rescore")
            input_psm_rescore.to_csv(str(fdr_dir) + "/rescore.tab", sep="\t", index=None)

        original_features_path = fdr_dir / "original_features_csm.tab"
        if not original_features_path.exists():
            shutil.copy(fdr_dir / "original.tab", original_features_path)
            input_psm_original = prepare_rescore_xl_psm_level(str(fdr_dir), "original")
            input_psm_original.to_csv(str(fdr_dir) + "/original.tab", sep="\t", index=None)

        expected_perc_tab_file = pd.read_csv(
            Path(__file__).parent / "data" / "xl" / "non-cleavable" / "expected_outputs" / "expected_rescore.tab", sep="\t"
        )

        created_perc_tab_file = pd.read_csv(
            fdr_dir / "rescore.tab",
            sep="\t",
        )

        try:
            assert_frame_equal(
                expected_perc_tab_file, created_perc_tab_file, check_dtype=True, check_exact=False, rtol=1e-1
            )
        except AssertionError as e:
            print("DataFrames are not equal:", e)
            raise  # Re-raise the assertion error for the test framework to catch

        config = Config()
        config.read(config_path)
        shutil.rmtree(Path(__file__).parent / "data" / "xl" / "non-cleavable" / "out")
