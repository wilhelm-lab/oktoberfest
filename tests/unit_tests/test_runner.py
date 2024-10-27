import shutil
import unittest
from pathlib import Path
from unittest.mock import patch

from oktoberfest.__main__ import main
from oktoberfest.utils import Config
import numpy as np
import pandas as pd


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

    def test_rescoring_xl(self):
        """Test the runner for a rescoring run with crosslinking."""
        config_path = Path(__file__).parent / "configs" / "rescoring_cleavable_xl.json"
        with patch("sys.argv", ["oktoberfest", f"--config_path={config_path}"]):
            main()

        #expected_perc_tab_file = pd.read_csv(Path(__file__).parent / "data" / "xl" / "cleavable" / "expected_outputs" / "expected_rescore.tab", sep="\t")
        #created_perc_tab_file = pd.read_csv(Path(__file__).parent / "data" / "xl" / "cleavable" / "out" / "results" / "percolator" / "rescore.tab", sep="\t")
        #np.testing.assert_almost_equal(expected_perc_tab_file.values, created_perc_tab_file.values)
        #pd.testing.assert_frame_equal(expected_perc_tab_file, created_perc_tab_file)

        config = Config()
        config.read(config_path)
        shutil.rmtree(Path(__file__).parent / "data" / "xl" / "cleavable" / "out")
        