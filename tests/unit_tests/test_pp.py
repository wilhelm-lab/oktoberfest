import unittest
from pathlib import Path

from oktoberfest import pp


class TestProcessing(unittest.TestCase):
    """Test class for preprocessing functions."""

    def test_list_spectra(self):
        """Test listing of spectra with expected user input."""
        spectra_path = Path(__file__).parent
        spectra_file = spectra_path / "test.mzml"
        spectra_file.open("w").close()
        self.assertEqual([spectra_path / "test.mzml"], pp.list_spectra(spectra_path, file_format="mzml"))
        spectra_file.unlink()

    def test_list_spectra_with_empty_string_folder(self):
        """Test listing spectra in a string folder without matching files."""
        self.assertRaises(AssertionError, pp.list_spectra, str(Path(__file__).parent), "raw")

    def test_list_spectra_with_wrong_folder(self):
        """Test listing spectra in a folder that does not exist."""
        self.assertRaises(NotADirectoryError, pp.list_spectra, Path(__file__).parent / "noexist", "raw")

    def test_list_spectra_with_wrong_format(self):
        """Test listing spectra with a format that isn't allowed."""
        self.assertRaises(ValueError, pp.list_spectra, Path(__file__).parent, "mzm")
