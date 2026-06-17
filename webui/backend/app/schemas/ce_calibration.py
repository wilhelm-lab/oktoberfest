from __future__ import annotations

from typing import Any, Optional

from pydantic import BaseModel, Field, field_validator

from app.schemas.common import CeAlignmentOptions, CommonConfigMixin
from app.schemas.rescoring import RescoringInputs


class CeCalibrationConfig(CommonConfigMixin):
    type: str = "CollisionEnergyCalibration"
    inputs: RescoringInputs = Field(default_factory=RescoringInputs)
    massTolerance: float = Field(default=20.0, gt=0)
    unitMassTolerance: str = "ppm"
    thermoExe: str = "ThermoRawFileParser.exe"
    ce_alignment_options: CeAlignmentOptions = Field(default_factory=CeAlignmentOptions)

    @field_validator("unitMassTolerance")
    @classmethod
    def validate_unit(cls, v):
        if v not in ("ppm", "da"):
            raise ValueError("unitMassTolerance must be ppm or da")
        return v

    def to_oktoberfest_dict(self, job_id: str, output_dir: str) -> dict[str, Any]:
        d: dict[str, Any] = {
            "type": self.type,
            "tag": self.tag.value if hasattr(self.tag, "value") else self.tag,
            "inputs": {
                "search_results": self.inputs.search_results or "",
                "search_results_type": self.inputs.search_results_type,
                "spectra": self.inputs.spectra or "",
                "spectra_type": self.inputs.spectra_type,
            },
            "output": output_dir,
            "models": {"intensity": self.models.intensity, "irt": self.models.irt},
            "prediction_server": self.prediction_server,
            "ssl": self.ssl,
            "thermoExe": self.thermoExe,
            "numThreads": self.numThreads,
            "massTolerance": self.massTolerance,
            "unitMassTolerance": self.unitMassTolerance,
            "ce_alignment_options": {
                "ce_range": self.ce_alignment_options.ce_range,
                "use_ransac_model": self.ce_alignment_options.use_ransac_model,
            },
        }
        if self.inputs.instrument_type:
            d["inputs"]["instrument_type"] = self.inputs.instrument_type
        if self.instrument_type:
            d["instrument_type"] = self.instrument_type.value
        if self.fragmentation_method:
            d["fragmentation_method"] = self.fragmentation_method.value
        if self.p_window is not None:
            d["p_window"] = self.p_window
        if self.static_mods:
            d["static_mods"] = self.static_mods
        if self.var_mods:
            d["var_mods"] = self.var_mods
        if self.ion_types:
            d["ion_types"] = self.ion_types
        return d
