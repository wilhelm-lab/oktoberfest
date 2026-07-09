from __future__ import annotations

from typing import Any, Optional, Union

from pydantic import BaseModel, Field, field_validator, model_validator

from app.schemas.common import (
    CeAlignmentOptions,
    CommonConfigMixin,
    FastaDigestOptions,
)


class RescoringInputs(BaseModel):
    search_results_type: str = "Maxquant"
    spectra_type: str = "raw"
    instrument_type: Optional[str] = "QE"
    library_input: Optional[str] = None  # FASTA, populated server-side

    # populated server-side after upload
    search_results: Optional[str] = None
    spectra: Optional[str] = None

    @field_validator("search_results_type")
    @classmethod
    def validate_search_engine(cls, v):
        allowed = {"Maxquant", "Msfragger", "Mascot", "Sage", "OpenMS", "Xisearch", ""}
        if v not in allowed:
            raise ValueError(f"search_results_type must be one of {sorted(allowed)}")
        return v

    @field_validator("spectra_type")
    @classmethod
    def validate_spectra_type(cls, v):
        if v not in ("raw", "mzml", "d", "hdf"):
            raise ValueError("spectra_type must be raw, mzml, d, or hdf")
        return v


class RescoringConfig(CommonConfigMixin):
    type: str = "Rescoring"
    inputs: RescoringInputs = Field(default_factory=RescoringInputs)
    fdr_estimation_method: str = "percolator"
    regressionMethod: str = "spline"
    add_feature_cols: Union[str, list[str]] = "none"
    quantification: bool = False
    massTolerance: float = Field(default=20.0, gt=0)
    unitMassTolerance: str = "ppm"
    thermoExe: str = "/cmnfs/scratch/oktoberfest_webserver/ThermoRawFileParser1.4.4/ThermoRawFileParser.exe"
    ce_alignment_options: CeAlignmentOptions = Field(default_factory=CeAlignmentOptions)
    fastaDigestOptions: Optional[FastaDigestOptions] = None

    @field_validator("fdr_estimation_method")
    @classmethod
    def validate_fdr(cls, v):
        if v not in ("percolator",):
            raise ValueError("fdr_estimation_method must be percolator")
        return v

    @field_validator("regressionMethod")
    @classmethod
    def validate_regression(cls, v):
        if v not in ("spline", "lowess", "logistic"):
            raise ValueError("regressionMethod must be spline, lowess, or logistic")
        return v

    @field_validator("unitMassTolerance")
    @classmethod
    def validate_unit(cls, v):
        if v not in ("ppm", "da"):
            raise ValueError("unitMassTolerance must be ppm or da")
        return v

    @model_validator(mode="after")
    def validate_quantification_deps(self) -> "RescoringConfig":
        if self.quantification:
            if not self.fastaDigestOptions:
                raise ValueError("fastaDigestOptions required when quantification=true")
        return self

    def to_oktoberfest_dict(self, job_id: str, output_dir: str) -> dict[str, Any]:
        """Serialize to exact Oktoberfest JSON shape."""
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
            "fdr_estimation_method": self.fdr_estimation_method,
            "add_feature_cols": self.add_feature_cols,
            "regressionMethod": self.regressionMethod,
            "massTolerance": self.massTolerance,
            "unitMassTolerance": self.unitMassTolerance,
            "ce_alignment_options": {
                "ce_range": self.ce_alignment_options.ce_range,
                "use_ransac_model": self.ce_alignment_options.use_ransac_model,
            },
        }
        if self.inputs.instrument_type:
            d["inputs"]["instrument_type"] = self.inputs.instrument_type
        if self.quantification:
            d["quantification"] = True
            if self.inputs.library_input:
                d["inputs"]["library_input"] = self.inputs.library_input
            if self.fastaDigestOptions:
                d["fastaDigestOptions"] = self.fastaDigestOptions.model_dump(exclude_none=True)
        else:
            d["quantification"] = False
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
