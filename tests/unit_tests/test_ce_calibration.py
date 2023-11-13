import shutil
import unittest
from pathlib import Path
from unittest.mock import patch

from oktoberfest.__main__ import main
from oktoberfest.utils import Config


class TestCECalibration(unittest.TestCase):
    """Test class for testing ce calibration in the case of HCD, CID and ransac regression."""

    def test_ce_calibration_cid(self):
        """Test ce calibration bypassing for CID data."""
        pass

    def test_ce_calibration_hcd(self):
        """Test ce calibration by highest scoring target PSMS without ransac for HCD data."""
        pass

    def test_ce_calibration_ransac(self):
        """Test ce calibration with a ransac regressor for timsTOF data."""
        config_path = Path(__file__).parent / "configs" / "ce_calib_ransac.json"
        with patch("sys.argv", ["oktoberfest", f"--config_path={config_path}"]):
            main()
        config = Config()
        config.read(config_path)
        shutil.rmtree(config.output)
