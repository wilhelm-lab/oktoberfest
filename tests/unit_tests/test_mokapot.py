import unittest
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import mokapot
import numpy as np

import oktoberfest.plotting as pl
from oktoberfest.utils import Config


def _run_mokapot(path: Path, search_type: str):
    file_path = path / f"{search_type}.tab"
    psms = mokapot.read_pin(file_path)
    results, models = mokapot.brew(psms, test_fdr=0.01)
    results.to_txt(dest_dir=path, file_root=f"{search_type}", decoys=True)


class TestMokapot(unittest.TestCase):
    """Test various plotting functions."""

    @unittest.skipIf(np.__version__ >= "2.0.0", "Skip mokapot test on Python 3.12")
    def test_mokapot_and_plot_all(self):
        """Test the mokapot execution and subsequent plotting."""
        path = Path(__file__).parent / "data/mokapot/"

        _run_mokapot(path, "original")
        _run_mokapot(path, "rescore")

        # switch to non-Gui, preventing plots being displayed
        plt.switch_backend("Agg")
        # suppress UserWarning that agg cannot show plots
        warnings.filterwarnings("ignore", "Matplotlib is currently using agg")

        pl.plot_all(path, Config())

        for file in path.glob("*"):
            if file.suffix != ".tab":
                file.unlink()
