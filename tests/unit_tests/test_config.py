import unittest

from oktoberfest.utils.config import Config


class TestConfig(unittest.TestCase):
    """Class to test spectra."""

    def test_num_threads(self):
        """Test num_threads."""
        config = Config()
        self.assertEqual(config.num_threads, 1)
        config.data = {"numThreads": 4}
        self.assertEqual(config.num_threads, 4)

    def test_prosit_server(self):
        """Test prosit_server."""
        config = Config()
        config.data = {"prosit_server": "10.152.135.57:8500"}
        self.assertEqual(config.prosit_server, "10.152.135.57:8500")

    def test_tag(self):
        """Test tag."""
        config = Config()
        self.assertEqual(config.tag, "tmt")
        config.data = {"tag": "itraq4"}
        self.assertEqual(config.tag, "itraq4")

    def test_all_features(self):
        """Test all_features."""
        config = Config()
        self.assertEqual(config.all_features, False)
        config.data = {"allFeatures": True}
        self.assertEqual(config.all_features, True)

    def test_curve_fitting_method(self):
        """Test curve_fitting_method."""
        config = Config()
        self.assertEqual(config.curve_fitting_method, "lowess")
        config.data = {"regressionMethod": "spline"}
        self.assertEqual(config.curve_fitting_method, "spline")

    def test_job_type(self):
        """Test job_type."""
        config = Config()
        config.data = {"jobType": "CollisionEnergyAlignment"}
        self.assertEqual(config.job_type, "CollisionEnergyAlignment")
        config.data = {"jobType": "Rescoring"}
        self.assertEqual(config.job_type, "Rescoring")

    def test_raw_type(self):
        """Test raw_type."""
        config = Config()
        self.assertEqual(config.raw_type, "thermo")
        config.data = {"fileUploads": {"raw_type": "mzml"}}
        self.assertEqual(config.raw_type, "mzml")

    def test_models(self):
        """Test models."""
        config = Config()
        config.data = {
            "models": {
                "intensity": "Prosit_2020_intensityTMT",
                "irt": "Prosit_2020_irt_TMT",
                "proteotypicity": "Prosit_2020_proteotypicity",
            }
        }
        self.assertEqual(
            config.models,
            {
                "intensity": "Prosit_2020_intensityTMT",
                "irt": "Prosit_2020_irt_TMT",
                "proteotypicity": "Prosit_2020_proteotypicity",
            },
        )

    def test_search_type(self):
        """Test search_type."""
        config = Config()
        self.assertEqual(config.search_type, "maxquant")
        config.data = {"fileUploads": {"search_type": "Msfragger"}}
        self.assertEqual(config.search_type, "msfragger")
