"""
The configs the web UI writes must pass oktoberfest's own Config.check() (needs oktoberfest installed).

Koina is not contacted: the model availability check is patched out.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from app.schemas.rescoring import RescoringConfig

oktoberfest_config = pytest.importorskip("oktoberfest.utils.config")

PTM_MODELS = {"intensity": "Prosit_2024_intensity_PTMs_gl", "irt": "Prosit_2024_irt_PTMs_gl"}


@pytest.fixture(autouse=True)
def no_koina(monkeypatch):
    monkeypatch.setattr(oktoberfest_config.Config, "_check_koina_model_availability", lambda self: None)


def _inputs(tmp_path: Path, search_results_type: str, *files: str) -> dict:
    for name in files:
        (tmp_path / name).write_text("")
    return {
        "search_results_type": search_results_type,
        "search_results": str(tmp_path / files[0]),
        "spectra_type": "mzml",
        "spectra": str(tmp_path),
    }


def _check(tmp_path: Path, cfg: RescoringConfig) -> oktoberfest_config.Config:
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(cfg.to_oktoberfest_dict("job1", str(tmp_path / "out"))))
    config = oktoberfest_config.Config()
    config.read(config_path)
    config.check()
    return config


def test_default_rescoring_passes_oktoberfest_check(tmp_path):
    _check(tmp_path, RescoringConfig(inputs=_inputs(tmp_path, "Maxquant", "msms.txt")))


def test_search_engine_quantification_passes_oktoberfest_check(tmp_path):
    inputs = _inputs(tmp_path, "Maxquant", "msms.txt", "evidence.txt", "db.fasta")
    inputs["library_input"] = str(tmp_path / "db.fasta")
    config = _check(tmp_path, RescoringConfig(inputs=inputs, quantification=True, fastaDigestOptions={}))
    assert config.quant_precursor == "search_engine"


def test_flashlfq_quantification_passes_oktoberfest_check(tmp_path, monkeypatch):
    from app.config import settings

    flashlfq = tmp_path / "CMD.dll"
    flashlfq.write_text("")
    monkeypatch.setattr(settings, "flashlfq_exe", str(flashlfq))
    inputs = _inputs(tmp_path, "Msfragger", "psm.tsv", "db.fasta")
    inputs["library_input"] = str(tmp_path / "db.fasta")
    cfg = RescoringConfig(
        inputs=inputs,
        quantification=True,
        fastaDigestOptions={},
        quantificationOptions={"precursor": "flashlfq", "ppmFilter": "none", "mbrFdr": 0.05},
    )
    config = _check(tmp_path, cfg)
    assert config.quant_precursor == "flashlfq"
    assert config.quant_ppm_filter is None


def test_ptm_localization_passes_oktoberfest_check(tmp_path):
    cfg = RescoringConfig(
        inputs=_inputs(tmp_path, "Maxquant", "msms.txt"),
        models=PTM_MODELS,
        ptm_localization=True,
        massTolerance="auto",
        chimericPeakSubtraction=True,
    )
    config = _check(tmp_path, cfg)
    assert config.ptm_localization
    assert config.mass_tolerance == "auto"
