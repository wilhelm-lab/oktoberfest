from __future__ import annotations

from typing import Any, Optional, Union

from pydantic import BaseModel, Field, field_validator, model_validator

from app.schemas.common import CommonConfigMixin, FastaDigestOptions


class SpectralLibraryOptions(BaseModel):
    fragmentation: str = "HCD"
    collisionEnergy: float = Field(default=30.0, gt=0, le=100)
    precursorCharge: list[int] = Field(default=[2, 3])
    minIntensity: float = Field(default=5e-4, gt=0)
    nrOx: int = Field(default=1, ge=0)
    batchsize: int = Field(default=10000, ge=1)
    format: str = "msp"

    @field_validator("fragmentation")
    @classmethod
    def validate_fragmentation(cls, v):
        if v not in ("HCD", "CID", ""):
            raise ValueError("fragmentation must be HCD or CID")
        return v

    @field_validator("format")
    @classmethod
    def validate_format(cls, v):
        if v not in ("msp", "spectronaut", "dlib"):
            raise ValueError("format must be msp, spectronaut, or dlib")
        return v

    @field_validator("precursorCharge")
    @classmethod
    def validate_charges(cls, v):
        if not v:
            raise ValueError("precursorCharge must not be empty")
        for c in v:
            if c < 1:
                raise ValueError("precursorCharge values must be >= 1")
        return v


class SpeclibInputs(BaseModel):
    library_input_type: str = "fasta"  # "fasta", "peptides", ""
    instrument_type: Optional[str] = "QE"
    library_input: Optional[str] = None  # populated server-side

    @field_validator("library_input_type")
    @classmethod
    def validate_lib_type(cls, v):
        if v not in ("fasta", "peptides", ""):
            raise ValueError("library_input_type must be 'fasta', 'peptides', or ''")
        return v


class SpeclibConfig(CommonConfigMixin):
    type: str = "SpectralLibraryGeneration"
    inputs: SpeclibInputs = Field(default_factory=SpeclibInputs)
    spectralLibraryOptions: SpectralLibraryOptions = Field(default_factory=SpectralLibraryOptions)
    fastaDigestOptions: Optional[FastaDigestOptions] = None

    @model_validator(mode="after")
    def validate_fasta_deps(self) -> "SpeclibConfig":
        if self.inputs.library_input_type == "fasta" and not self.fastaDigestOptions:
            # Auto-populate with defaults instead of raising
            self.fastaDigestOptions = FastaDigestOptions(fragmentation=self.spectralLibraryOptions.fragmentation)
        return self

    def to_oktoberfest_dict(self, job_id: str, output_dir: str) -> dict[str, Any]:
        d: dict[str, Any] = {
            "type": self.type,
            "tag": self.tag.value if hasattr(self.tag, "value") else self.tag,
            "inputs": {
                "library_input": self.inputs.library_input or "",
                "library_input_type": self.inputs.library_input_type,
            },
            "output": output_dir,
            "models": {"intensity": self.models.intensity, "irt": self.models.irt},
            "prediction_server": self.prediction_server,
            "ssl": self.ssl,
            "numThreads": self.numThreads,
            "spectralLibraryOptions": {
                "fragmentation": self.spectralLibraryOptions.fragmentation,
                "collisionEnergy": int(self.spectralLibraryOptions.collisionEnergy),
                "precursorCharge": self.spectralLibraryOptions.precursorCharge,
                "minIntensity": self.spectralLibraryOptions.minIntensity,
                "nrOx": self.spectralLibraryOptions.nrOx,
                "batchsize": self.spectralLibraryOptions.batchsize,
                "format": self.spectralLibraryOptions.format,
            },
        }
        if self.inputs.instrument_type:
            d["inputs"]["instrument_type"] = self.inputs.instrument_type
        if self.fastaDigestOptions:
            fd = self.fastaDigestOptions
            d["fastaDigestOptions"] = {
                "fragmentation": fd.fragmentation,
                "digestion": fd.digestion,
                "missedCleavages": fd.missedCleavages,
                "minLength": fd.minLength,
                "maxLength": fd.maxLength,
                "enzyme": fd.enzyme,
                "specialAas": fd.specialAas,
                "db": fd.db,
            }
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
        return d
