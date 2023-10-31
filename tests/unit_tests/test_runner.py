import unittest
from pathlib import Path
from unittest.mock import patch

from oktoberfest.__main__ import main


class TestRunner(unittest.TestCase):
    """Test class for use cases of runner."""

    def test_speclib_digest(self):
        """Test the runner for a spectral library generation with a fasta digest."""
        config_path = Path(__file__).parent / "configs" / "spectral_library_with_digest.json"
        with patch("sys.argv", ["oktoberfest", f"--config_path={config_path}"]):
            main()
