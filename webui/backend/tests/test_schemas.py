"""Phase 2 tests: Pydantic config model validation and serialization."""

from __future__ import annotations

import pytest

from app.schemas.rescoring import RescoringConfig
from app.schemas.ce_calibration import CeCalibrationConfig
from app.schemas.speclib import SpeclibConfig
from app.schemas.common import CeAlignmentOptions, FastaDigestOptions


# ── Rescoring ─────────────────────────────────────────────────────────────────


def test_rescoring_accepts_valid_config():
    cfg = RescoringConfig(
        inputs={"search_results_type": "Maxquant", "spectra_type": "raw"},
        fdr_estimation_method="mokapot",
    )
    assert cfg.type == "Rescoring"
    assert cfg.fdr_estimation_method == "mokapot"


def test_rescoring_rejects_bad_fdr_method():
    with pytest.raises(Exception, match="fdr_estimation_method"):
        RescoringConfig(fdr_estimation_method="invalid")


def test_rescoring_rejects_bad_search_engine():
    with pytest.raises(Exception):
        RescoringConfig(inputs={"search_results_type": "NotAnEngine"})


def test_rescoring_quant_requires_fasta_digest():
    with pytest.raises(Exception, match="fastaDigestOptions"):
        RescoringConfig(quantification=True)


def test_rescoring_serializes_to_oktoberfest_shape():
    cfg = RescoringConfig(
        inputs={
            "search_results_type": "Maxquant",
            "spectra_type": "mzml",
            "search_results": "/data/job1/inputs/msms.txt",
            "spectra": "/data/job1/inputs/run.mzML",
        },
    )
    d = cfg.to_oktoberfest_dict("job1", "/data/job1/output")
    assert d["type"] == "Rescoring"
    assert d["output"] == "/data/job1/output"
    assert "inputs" in d
    assert d["inputs"]["search_results"] == "/data/job1/inputs/msms.txt"
    # ce_range must be a list (not tuple)
    assert isinstance(d["ce_alignment_options"]["ce_range"], list)


# ── CE Calibration ────────────────────────────────────────────────────────────


def test_ce_calibration_accepts_valid():
    cfg = CeCalibrationConfig()
    assert cfg.type == "CollisionEnergyCalibration"


def test_ce_calibration_serializes_correctly():
    cfg = CeCalibrationConfig()
    d = cfg.to_oktoberfest_dict("j2", "/out")
    assert d["type"] == "CollisionEnergyCalibration"
    assert isinstance(d["ce_alignment_options"]["ce_range"], list)


# ── Spectral Library Generation ───────────────────────────────────────────────


def test_speclib_accepts_valid():
    cfg = SpeclibConfig(
        inputs={"library_input_type": "fasta"},
        spectralLibraryOptions={"format": "msp"},
    )
    assert cfg.type == "SpectralLibraryGeneration"
    assert cfg.fastaDigestOptions is not None  # auto-populated


def test_speclib_rejects_bad_format():
    with pytest.raises(Exception):
        SpeclibConfig(
            inputs={"library_input_type": "fasta"},
            spectralLibraryOptions={"format": "invalid_format"},
        )


def test_speclib_precursor_charge_must_be_positive():
    with pytest.raises(Exception, match=">="):
        SpeclibConfig(
            inputs={"library_input_type": "fasta"},
            spectralLibraryOptions={"precursorCharge": [0, 1]},
        )


def test_speclib_serializes_precursor_charge_as_list():
    cfg = SpeclibConfig(inputs={"library_input_type": "fasta"})
    d = cfg.to_oktoberfest_dict("j3", "/out")
    assert isinstance(d["spectralLibraryOptions"]["precursorCharge"], list)
    assert d["fastaDigestOptions"]["enzyme"] == "trypsin"


# ── CeAlignmentOptions ────────────────────────────────────────────────────────


def test_ce_range_min_must_be_less_than_max():
    with pytest.raises(Exception, match="min must be less than max"):
        CeAlignmentOptions(ce_range=[50, 19])


def test_ce_range_must_have_two_elements():
    with pytest.raises(Exception, match="exactly 2"):
        CeAlignmentOptions(ce_range=[19])


def test_fasta_digest_maxlength_must_ge_minlength():
    with pytest.raises(Exception, match="maxLength"):
        FastaDigestOptions(minLength=10, maxLength=5)
