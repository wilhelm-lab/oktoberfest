import unittest
from pathlib import Path

from oktoberfest.data import Spectra
from oktoberfest.predict.alignment import _alignment, _prepare_alignment_df


class TestAlignment(unittest.TestCase):
    """Test alignment utils."""

    @classmethod
    def setUpClass(cls):  # noqa: D102
        cls.spectra = Spectra.from_hdf5(Path(__file__).parent / "data/spectra/test_spectra.hdf5")
        cls.ce_range = (18, 50)

    def test_alignment(self):
        """Test alignment of predicted vs. raw intensities."""
        library = None  # TODO
        alignment_library = _prepare_alignment_df(library, ce_range=self.ce_range, group_by_charge=group_by_charge)
        self.predict_intensities(data=alignment_library, chunk_idx=chunk_idx, keep_dataset=False, **kwargs)
        _alignment(self.spectra.copy())
