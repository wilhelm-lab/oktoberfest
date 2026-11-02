from __future__ import annotations

from typing import Any, Literal, Optional, Union

from pydantic import BaseModel, Field, field_validator, model_validator

from app.config import settings
from app.schemas.common import (
    CeAlignmentOptions,
    CommonConfigMixin,
    FastaDigestOptions,
)

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# search engines whose own quantification can be updated with the rescored PSMs from flat uploads (FragPipe needs
# psm.tsv files in per-experiment subfolders, which uploads cannot reproduce)
SEARCH_ENGINE_QUANT_ENGINES = {"Maxquant", "Sage"}


class QuantificationOptions(BaseModel):
    """quantificationOptions of Oktoberfest; flashlfqExe and dotnetExe are set by the server."""

    precursor: Literal["search_engine", "flashlfq"] = "search_engine"
    protein: Literal["directlfq", "none"] = "directlfq"
    normalization: Literal["stable_drift", "drift", "mode", "median", "none"] = "stable_drift"
    psmFdr: float = Field(default=0.01, gt=0, le=1)
    peptideFdr: float = Field(default=0.01, gt=0, le=1)
    proteinFdr: float = Field(default=0.01, gt=0, le=1)
    ppmTolerance: float = Field(default=10.0, gt=0)
    # "auto", a positive number of ppm, or "none" to keep all peaks (null in the Oktoberfest config)
    ppmFilter: Union[Literal["auto", "none"], float] = "auto"
    # PIP q-value cutoff of match-between-runs peaks, None disables match-between-runs
    mbrFdr: Optional[float] = Field(default=None, gt=0, le=1)

    @field_validator("ppmFilter")
    @classmethod
    def validate_ppm_filter(cls, v):
        if not isinstance(v, str) and v <= 0:
            raise ValueError('ppmFilter must be a positive number, "auto" or "none"')
        return v

    def to_oktoberfest_dict(self) -> dict[str, Any]:
        d = self.model_dump(exclude_none=True)
        if self.ppmFilter == "none":
            d["ppmFilter"] = None
        if self.precursor == "flashlfq":
            d["flashlfqExe"] = settings.flashlfq_exe
            d["dotnetExe"] = settings.dotnet_exe
        return d


class PtmLocalizationOptions(BaseModel):
    """ptmLocalizationOptions of Oktoberfest, defaults as in the published phospho benchmark runs."""

    # unimod id -> residues the PTM can be located on, e.g. {"21": ["S", "T", "Y"]}
    ptms: dict[str, list[str]] = Field(default_factory=lambda: {"21": ["S", "T", "Y"]})
    neutral_loss: bool = True
    # decoy residues by priority; None: Oktoberfest default (A, G, L, V without the possible sites); [] disables them
    decoy_residues: Optional[list[str]] = Field(default_factory=lambda: ["A", "G", "L", "V"])
    decoy_neutral_loss: bool = False
    # residue decoy candidates are predicted on, or "search"; None: first possible site
    decoy_proxy: Optional[str] = None
    decoy_selection: Literal["residue", "adjacent"] = "residue"
    flr_neighbor_split: bool = False
    loc_probability: Literal["pep", "score"] = "pep"
    decoy_weighting: Literal["count", "distance"] = "count"
    flr_groups: Literal["none", "residue"] = "none"
    flr_estimator: Literal["hybrid", "decoy"] = "hybrid"
    sdi_features: bool = True
    ambiguity_check: bool = True

    @field_validator("ptms")
    @classmethod
    def validate_ptms(cls, v: dict[str, list[str]]) -> dict[str, list[str]]:
        if not v:
            raise ValueError("ptms must contain at least one PTM to localize")
        for unimod_id, sites in v.items():
            if not str(unimod_id).isdigit() or int(unimod_id) <= 0:
                raise ValueError(f"ptms keys must be unimod ids, got {unimod_id!r}")
            if not sites:
                raise ValueError(f"no possible residues given for UNIMOD:{unimod_id}")
            invalid = [aa for aa in sites if aa not in AMINO_ACIDS]
            if invalid:
                raise ValueError(f"invalid residues {invalid} for UNIMOD:{unimod_id}")
        return v

    @model_validator(mode="after")
    def validate_decoys(self) -> "PtmLocalizationOptions":
        sites = {aa for site_list in self.ptms.values() for aa in site_list}
        decoys = self.decoy_residues or []
        invalid = [aa for aa in decoys if aa not in AMINO_ACIDS]
        if invalid:
            raise ValueError(f"invalid decoy residues {invalid}")
        overlapping = [aa for aa in decoys if aa in sites]
        if overlapping:
            raise ValueError(f"decoy residues {overlapping} can also carry a localized PTM")
        if self.decoy_proxy is not None and self.decoy_proxy != "search" and self.decoy_proxy not in sites:
            raise ValueError('decoy_proxy must be one of the possible sites or "search"')
        if self.decoy_neutral_loss and not self.neutral_loss:
            raise ValueError("decoy_neutral_loss requires neutral_loss")
        return self


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
    # positive number or "auto" (fragment tolerance estimated per spectra file, ppm only)
    massTolerance: Union[Literal["auto"], float] = 20.0
    unitMassTolerance: str = "ppm"
    thermoExe: str = "/cmnfs/scratch/oktoberfest_webserver/ThermoRawFileParser1.4.4/ThermoRawFileParser.exe"
    ce_alignment_options: CeAlignmentOptions = Field(default_factory=CeAlignmentOptions)
    fastaDigestOptions: Optional[FastaDigestOptions] = None
    # train percolator on at most this many PSMs (rescoring only, never localization): ~6x faster and ~3x less RAM
    # on large inputs at the same identifications and entrapment FDP; None trains on all PSMs
    percolatorSubsetMaxTrain: Optional[int] = Field(default=500000, ge=1)
    chimericPeakSubtraction: bool = False
    quantificationOptions: QuantificationOptions = Field(default_factory=QuantificationOptions)
    ptm_localization: bool = False
    # None: Oktoberfest default (true for PTM models without localization)
    ptm_fdr_control: Optional[bool] = None
    ptmLocalizationOptions: Optional[PtmLocalizationOptions] = None

    @field_validator("massTolerance")
    @classmethod
    def validate_mass_tolerance(cls, v):
        if not isinstance(v, str) and v <= 0:
            raise ValueError('massTolerance must be a positive number or "auto"')
        return v

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
        if not self.quantification:
            return self
        if self.quantificationOptions.precursor == "search_engine":
            if not self.fastaDigestOptions:
                raise ValueError("fastaDigestOptions required when quantification=true")
            if self.inputs.search_results_type not in SEARCH_ENGINE_QUANT_ENGINES:
                raise ValueError(
                    f"search engine quantification supports {sorted(SEARCH_ENGINE_QUANT_ENGINES)}; "
                    "use quantificationOptions.precursor flashlfq for other search engines"
                )
        else:
            if self.ptm_localization:
                raise ValueError("flashlfq quantification is not supported together with ptm_localization")
            if self.inputs.spectra_type not in ("raw", "mzml"):
                raise ValueError("flashlfq quantification needs raw or mzml spectra (MS1 scans)")
            if not settings.flashlfq_exe:
                raise ValueError("flashlfq quantification is not available: FLASHLFQ_EXE is not set on this server")
        return self

    @model_validator(mode="after")
    def validate_ptm_localization(self) -> "RescoringConfig":
        if self.ptm_localization:
            if self.ptmLocalizationOptions is None:
                self.ptmLocalizationOptions = PtmLocalizationOptions()
            if "PTM" not in self.models.intensity:
                raise ValueError("ptm_localization requires a PTM intensity model, e.g. Prosit_2024_intensity_PTMs_gl")
        if self.massTolerance == "auto" and self.unitMassTolerance != "ppm":
            raise ValueError('massTolerance "auto" requires unitMassTolerance ppm')
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
            d["quantificationOptions"] = self.quantificationOptions.to_oktoberfest_dict()
            if self.inputs.library_input:
                d["inputs"]["library_input"] = self.inputs.library_input
                d["inputs"]["library_input_type"] = "fasta"
            if self.fastaDigestOptions:
                d["fastaDigestOptions"] = self.fastaDigestOptions.model_dump(exclude_none=True)
        else:
            d["quantification"] = False
        if self.percolatorSubsetMaxTrain is not None:
            d["percolatorSubsetMaxTrain"] = self.percolatorSubsetMaxTrain
        if self.chimericPeakSubtraction:
            d["chimericPeakSubtraction"] = True
        d["ptm_localization"] = self.ptm_localization
        if self.ptm_localization and self.ptmLocalizationOptions:
            d["ptmLocalizationOptions"] = self.ptmLocalizationOptions.model_dump(exclude_none=True)
        if self.ptm_fdr_control is not None:
            d["ptm_fdr_control"] = self.ptm_fdr_control
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
