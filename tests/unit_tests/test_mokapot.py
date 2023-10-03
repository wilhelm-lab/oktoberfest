import unittest
from pathlib import Path

import mokapot
import pandas as pd

import oktoberfest.plotting as pl


def _run_mokapot(path: Path, search_type: str):
    file_path = path / f"{search_type}.tab"
    df = pd.read_csv(file_path, sep="\t")
    df = df.rename(columns={"Protein": "Proteins"})
    df.to_csv(file_path, sep="\t")
    psms = mokapot.read_pin(file_path)
    results, models = mokapot.brew(psms, test_fdr=0.01)
    results.to_txt(dest_dir=path, file_root=f"{search_type}", decoys=True)


class TestMokapot(unittest.TestCase):
    """Test various plotting functions."""

    def test_mokapot_and_plot_all(self):
        """Test the mokapot execution and subsequent plotting."""
        path = Path(__file__).parent / "data/mokapot/"

        _run_mokapot(path, "original")
        _run_mokapot(path, "rescore")

        pl.plot_all(path)
