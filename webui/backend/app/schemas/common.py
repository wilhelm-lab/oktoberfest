from __future__ import annotations

from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, Field, field_validator, model_validator


class TagEnum(str, Enum):
    none = ""
    tmt = "tmt"
    tmtpro = "tmtpro"
    itraq4 = "itraq4"
    itraq8 = "itraq8"


class InstrumentTypeEnum(str, Enum):
    QE = "QE"
    LUMOS = "LUMOS"
    TIMSTOF = "TIMSTOF"
    SCIEXTOF = "SCIEXTOF"


class FragmentationMethodEnum(str, Enum):
    HCD = "HCD"
    ECD = "ECD"
    EID = "EID"
    ETciD = "ETciD"
    UVPD = "UVPD"


class ModelsBlock(BaseModel):
    intensity: str = "Prosit_2020_intensity_HCD"
    irt: str = "Prosit_2019_irt"


class CeAlignmentOptions(BaseModel):
    ce_range: list[int] = Field(default=[19, 50])
    use_ransac_model: bool = False

    @field_validator("ce_range")
    @classmethod
    def validate_ce_range(cls, v: list[int]) -> list[int]:
        if len(v) != 2:
            raise ValueError("ce_range must be a list of exactly 2 integers [min, max]")
        if v[0] >= v[1]:
            raise ValueError("ce_range min must be less than max")
        if not (1 <= v[0] <= 100 and 1 <= v[1] <= 100):
            raise ValueError("ce_range values must be between 1 and 100")
        return v


class FastaDigestOptions(BaseModel):
    digestion: Optional[str] = "full"  # "full", "semi", null
    missedCleavages: int = 2
    minLength: int = 7
    maxLength: int = 60
    enzyme: str = "trypsin"
    specialAas: str = "KR"
    db: str = "concat"  # "target", "decoy", "concat"
    fragmentation: str = "HCD"

    @field_validator("digestion")
    @classmethod
    def validate_digestion(cls, v):
        if v is not None and v not in ("full", "semi", "none", "null"):
            raise ValueError("digestion must be 'full', 'semi', 'none', or null")
        return v

    @field_validator("db")
    @classmethod
    def validate_db(cls, v):
        if v not in ("target", "decoy", "concat"):
            raise ValueError("db must be 'target', 'decoy', or 'concat'")
        return v

    @model_validator(mode="after")
    def validate_lengths(self) -> "FastaDigestOptions":
        if self.maxLength < self.minLength:
            raise ValueError("maxLength must be >= minLength")
        return self


class CommonConfigMixin(BaseModel):
    """Fields shared by all three job types (excluding type, output, file paths)."""

    tag: TagEnum = TagEnum.none
    models: ModelsBlock = Field(default_factory=ModelsBlock)
    prediction_server: str = "koina.wilhelmlab.org:443"
    ssl: bool = True
    numThreads: int = Field(default=1, ge=1)
    instrument_type: Optional[InstrumentTypeEnum] = None
    fragmentation_method: Optional[FragmentationMethodEnum] = None
    p_window: Optional[float] = Field(default=None, gt=0)
    static_mods: Optional[dict[str, Any]] = None
    var_mods: Optional[dict[str, Any]] = None
    ion_types: Optional[str] = None


# Job status enum
class JobStatus(str, Enum):
    CREATED = "CREATED"
    QUEUED = "QUEUED"
    RUNNING = "RUNNING"
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"


TERMINAL_STATUSES = {JobStatus.SUCCEEDED, JobStatus.FAILED, JobStatus.CANCELLED}
LEGAL_TRANSITIONS: dict[JobStatus, set[JobStatus]] = {
    JobStatus.CREATED: {JobStatus.QUEUED},
    JobStatus.QUEUED: {JobStatus.RUNNING, JobStatus.CANCELLED},
    JobStatus.RUNNING: {JobStatus.SUCCEEDED, JobStatus.FAILED, JobStatus.CANCELLED},
    JobStatus.SUCCEEDED: set(),
    JobStatus.FAILED: set(),
    JobStatus.CANCELLED: set(),
}
