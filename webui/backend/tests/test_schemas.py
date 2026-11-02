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
        fdr_estimation_method="percolator",
    )
    assert cfg.type == "Rescoring"
    assert cfg.fdr_estimation_method == "percolator"


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
    assert d["percolatorSubsetMaxTrain"] == 500000


def test_rescoring_percolator_subset_can_be_disabled():
    d = RescoringConfig(percolatorSubsetMaxTrain=None).to_oktoberfest_dict("job1", "/out")
    assert "percolatorSubsetMaxTrain" not in d
    with pytest.raises(Exception, match="percolatorSubsetMaxTrain"):
        RescoringConfig(percolatorSubsetMaxTrain=0)


def test_rescoring_quant_search_engine_serializes_fasta_and_options():
    cfg = RescoringConfig(
        inputs={"search_results_type": "Maxquant", "library_input": "/data/job1/inputs/db.fasta"},
        quantification=True,
        fastaDigestOptions={},
    )
    d = cfg.to_oktoberfest_dict("job1", "/out")
    assert d["inputs"]["library_input"] == "/data/job1/inputs/db.fasta"
    assert d["inputs"]["library_input_type"] == "fasta"
    assert d["quantificationOptions"]["precursor"] == "search_engine"
    assert "flashlfqExe" not in d["quantificationOptions"]


def test_rescoring_quant_search_engine_rejects_unsupported_engine():
    with pytest.raises(Exception, match="search engine quantification"):
        RescoringConfig(inputs={"search_results_type": "Msfragger"}, quantification=True, fastaDigestOptions={})


def test_rescoring_quant_flashlfq(monkeypatch):
    from app.config import settings

    monkeypatch.setattr(settings, "flashlfq_exe", "/tools/CMD.dll")
    cfg = RescoringConfig(
        inputs={"search_results_type": "Msfragger", "spectra_type": "mzml"},
        quantification=True,
        quantificationOptions={"precursor": "flashlfq", "ppmFilter": "none", "mbrFdr": 0.05},
    )
    quant = cfg.to_oktoberfest_dict("job1", "/out")["quantificationOptions"]
    assert quant["ppmFilter"] is None
    assert quant["mbrFdr"] == 0.05
    assert quant["flashlfqExe"] == "/tools/CMD.dll"


def test_rescoring_quant_flashlfq_checks(monkeypatch):
    from app.config import settings

    monkeypatch.setattr(settings, "flashlfq_exe", "")
    with pytest.raises(Exception, match="FLASHLFQ_EXE"):
        RescoringConfig(inputs={"spectra_type": "mzml"}, quantification=True, quantificationOptions={"precursor": "flashlfq"})
    monkeypatch.setattr(settings, "flashlfq_exe", "/tools/CMD.dll")
    with pytest.raises(Exception, match="raw or mzml"):
        RescoringConfig(inputs={"spectra_type": "d"}, quantification=True, quantificationOptions={"precursor": "flashlfq"})
    with pytest.raises(Exception, match="ptm_localization"):
        RescoringConfig(
            inputs={"spectra_type": "mzml"},
            models={"intensity": "Prosit_2024_intensity_PTMs_gl", "irt": "Prosit_2024_irt_PTMs_gl"},
            quantification=True,
            quantificationOptions={"precursor": "flashlfq"},
            ptm_localization=True,
        )


def test_rescoring_ptm_localization_serializes():
    cfg = RescoringConfig(
        models={"intensity": "Prosit_2024_intensity_PTMs_gl", "irt": "Prosit_2024_irt_PTMs_gl"},
        ptm_localization=True,
        ptmLocalizationOptions={"ptms": {"21": ["S", "T", "Y"]}, "neutral_loss": True, "flr_estimator": "decoy"},
    )
    d = cfg.to_oktoberfest_dict("job1", "/out")
    assert d["ptm_localization"] is True
    assert d["ptmLocalizationOptions"]["ptms"] == {"21": ["S", "T", "Y"]}
    assert d["ptmLocalizationOptions"]["flr_estimator"] == "decoy"
    assert d["ptmLocalizationOptions"]["decoy_residues"] == ["A", "G", "L", "V"]


def test_rescoring_ptm_localization_defaults_match_published_runs():
    cfg = RescoringConfig(
        models={"intensity": "Prosit_2024_intensity_PTMs_gl", "irt": "Prosit_2024_irt_PTMs_gl"},
        ptm_localization=True,
    )
    loc = cfg.to_oktoberfest_dict("job1", "/out")["ptmLocalizationOptions"]
    assert loc["ptms"] == {"21": ["S", "T", "Y"]}
    assert loc["neutral_loss"] is True
    assert loc["decoy_residues"] == ["A", "G", "L", "V"]
    assert (loc["decoy_selection"], loc["loc_probability"], loc["decoy_weighting"]) == ("residue", "pep", "count")
    assert (loc["flr_groups"], loc["flr_estimator"]) == ("none", "hybrid")
    assert loc["sdi_features"] and loc["ambiguity_check"]
    assert not loc["flr_neighbor_split"] and not loc["decoy_neutral_loss"]
    assert "decoy_proxy" not in loc


def test_rescoring_ptm_localization_checks():
    ptm_models = {"intensity": "Prosit_2024_intensity_PTMs_gl", "irt": "Prosit_2024_irt_PTMs_gl"}
    with pytest.raises(Exception, match="PTM intensity model"):
        RescoringConfig(ptm_localization=True)
    with pytest.raises(Exception, match="can also carry"):
        RescoringConfig(models=ptm_models, ptm_localization=True,
                        ptmLocalizationOptions={"ptms": {"21": ["S", "T"]}, "decoy_residues": ["S"]})
    with pytest.raises(Exception, match="invalid residues"):
        RescoringConfig(models=ptm_models, ptm_localization=True, ptmLocalizationOptions={"ptms": {"21": ["B"]}})
    with pytest.raises(Exception, match="decoy_neutral_loss"):
        RescoringConfig(models=ptm_models, ptm_localization=True,
                        ptmLocalizationOptions={"neutral_loss": False, "decoy_neutral_loss": True})


def test_rescoring_new_rescoring_options():
    d = RescoringConfig(massTolerance="auto", chimericPeakSubtraction=True, ptm_fdr_control=False).to_oktoberfest_dict(
        "job1", "/out"
    )
    assert d["massTolerance"] == "auto"
    assert d["chimericPeakSubtraction"] is True
    assert d["ptm_fdr_control"] is False
    assert d["ptm_localization"] is False
    with pytest.raises(Exception, match="massTolerance"):
        RescoringConfig(massTolerance=0)
    with pytest.raises(Exception, match="ppm"):
        RescoringConfig(massTolerance="auto", unitMassTolerance="da")


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


def test_thermo_exe_comes_from_server_settings(monkeypatch):
    from app.config import settings

    monkeypatch.setattr(settings, "thermo_exe", "/opt/ThermoRawFileParser/ThermoRawFileParser.exe")
    # a path in the job config is ignored: the server decides which executable runs
    cfg = RescoringConfig(thermoExe="/tmp/evil.exe")
    assert cfg.to_oktoberfest_dict("job1", "/out")["thermoExe"] == "/opt/ThermoRawFileParser/ThermoRawFileParser.exe"


def test_prediction_server_default_follows_settings(monkeypatch):
    from app.config import settings

    monkeypatch.setattr(settings, "prediction_server", "koina.example.org:443")
    assert RescoringConfig().prediction_server == "koina.example.org:443"
    assert RescoringConfig(prediction_server="other:443").prediction_server == "other:443"
