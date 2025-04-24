"""Tests for plot_mirror_spectrum in oktoberfest.plotting."""

from io import BytesIO
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
import pytest
from matplotlib.backends.backend_pdf import PdfPages

from oktoberfest.data.spectra import Spectra
from oktoberfest.plotting import plot_mirror_spectrum


@pytest.fixture
def sample_data():
    """Sample Spectra, mzML, prosit, target, decoy, and config."""
    mock_spectra = Mock(spec=Spectra)
    mock_spectra.obs = pd.DataFrame(
        [
            {
                "SCAN_NUMBER": 123,
                "MODIFIED_SEQUENCE": "PEPTIDE",
                "PRECURSOR_CHARGE": 2,
                "MASS": 500.0,
                "MASS_ANALYZER": "FTMS",
                "FRAGMENTATION": "HCD",
                "RAW_FILE": "file.mzML",
                "RETENTION_TIME": 35.4,
                "COLLISION_ENERGY": 30.0,
            }
        ]
    )
    mock_spectra.layers = {
        "mz": Mock(toarray=lambda: np.array([[100, 200, 300]])),
        "pred_int": Mock(toarray=lambda: np.array([[0.1, 0.5, 1.0]])),
    }

    mzml = pd.DataFrame(
        {
            "SCAN_NUMBER": [123],
            "MZ": [np.array([100.0, 200.0, 300.0])],
            "INTENSITIES": [np.array([10.0, 20.0, 30.0])],
        }
    )
    prosit_df = pd.DataFrame(
        {
            "ScanNr": [123],
            "filename": ["file.mzML"],
            "spectral_angle": [0.85],
            "abs_rt_diff": [0.5],
            "collision_energy_aligned": [30.0],
        }
    )
    target_df = pd.DataFrame(
        {
            "ScanNr": [123],
            "SpecId": [123],
            "filename": ["file.mzML"],
            "mokapot score": [0.9],
            "mokapot q-value": [0.01],
            "Peptide": ["PEPTIDE"],
        }
    )
    decoy_df = pd.DataFrame(
        {
            "ScanNr": [123],
            "SpecId": [123],
            "filename": ["file.mzML"],
            "mokapot score": [-0.4],
            "mokapot q-value": [0.05],
            "Peptide": ["DECOY"],
        }
    )

    config = Mock()
    config.models = {"intensity": "default_model"}
    config.ion_types = "by"

    return mock_spectra, mzml, prosit_df, target_df, decoy_df, config


@patch(
    "oktoberfest.plotting._check_columns",
    return_value=("mokapot score", "mokapot q-value", "Peptide", "SpecId"),
)
@patch("spectrum_utils.spectrum.MsmsSpectrum")
@patch("spectrum_utils.plot.mirror")
@patch("seaborn.kdeplot")
def test_plot_mirror_spectrum(
    mock_kdeplot: Mock,
    mock_mirror: Mock,
    mock_msms_spectrum: Mock,
    mock_check_columns: Mock,
    sample_data: tuple,
):
    """
    Verify plot_mirror_spectrum calls the expected plotting routines.

    :param mock_kdeplot: mocked seaborn.kdeplot
    :param mock_mirror: mocked spectrum_utils.plot.mirror
    :param mock_msms_spectrum: mocked MsmsSpectrum constructor
    :param mock_check_columns: mocked column checker
    :param sample_data: fixture returning test inputs
    """
    spec_pred, mzml, prosit_df, target_df, decoy_df, config = sample_data
    pdf_buffer = BytesIO()
    pdf = PdfPages(pdf_buffer)

    plot_mirror_spectrum(
        spec_pred=spec_pred,
        mzml=mzml,
        raw_file="file.mzML",
        scan_number=123,
        config=config,
        prosit_df=prosit_df,
        target_df=target_df,
        decoy_df=decoy_df,
        pdf=pdf,
    )

    assert mock_mirror.called
    assert mock_kdeplot.call_count >= 2
    assert mock_msms_spectrum.call_count == 2

    pdf.close()
