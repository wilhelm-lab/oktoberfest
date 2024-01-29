import unittest


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
        pass
