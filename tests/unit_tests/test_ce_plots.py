"""Tests for the CE calibration plots in oktoberfest.plotting."""

from unittest.mock import patch

import matplotlib
import pandas as pd
import pytest

from oktoberfest.plotting import plot_violin_sa_ce

matplotlib.use("Agg")


def _violin_line_x(sa_ce_df, tmp_path):
    """Return the x position of the red line drawn by plot_violin_sa_ce and the CE tick labels."""
    with patch("oktoberfest.plotting.plotting.plt.savefig"), patch("oktoberfest.plotting.plotting.plt.close"):
        plot_violin_sa_ce(sa_ce_df, tmp_path / "violin.svg")
        ax = matplotlib.pyplot.gcf().axes[0]
        line = [line for line in ax.get_lines() if line.get_color() == "red"][0]
        labels = [float(t.get_text()) for t in ax.get_xticklabels()]
    matplotlib.pyplot.close("all")
    return line.get_xdata()[0], labels


@pytest.mark.parametrize(
    "ces, sas, best",
    [
        # best mean at 30, but the single highest spectral angle is at 31
        ([30, 30, 31, 31], [0.8, 0.8, 0.95, 0.2], 30),
        # CEs not consecutive: 32 is the second violin, not 32 - 30 = 2
        ([30, 30, 32, 32, 34, 34], [0.5, 0.5, 0.9, 0.9, 0.4, 0.4], 32),
    ],
)
def test_violin_marks_best_mean_ce(ces, sas, best, tmp_path):
    """The red line sits on the violin of the CE with the highest mean spectral angle."""
    x, labels = _violin_line_x(pd.DataFrame({"COLLISION_ENERGY": ces, "SPECTRAL_ANGLE": sas}), tmp_path)
    assert labels[int(x)] == best
