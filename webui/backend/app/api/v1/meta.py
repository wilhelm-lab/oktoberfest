from __future__ import annotations

from fastapi import APIRouter

from app.services.catalog import (
    DEFAULTS,
    ENZYMES,
    FRAGMENTATION_METHODS,
    INSTRUMENT_TYPES,
    INTENSITY_MODELS,
    IRT_MODELS,
    LIBRARY_FORMATS,
    SEARCH_ENGINES,
    SPECTRA_TYPES,
    TAGS,
)

router = APIRouter(tags=["meta"])


@router.get("/models")
def get_models():
    return {"intensity": INTENSITY_MODELS, "irt": IRT_MODELS}


@router.get("/search-engines")
def get_search_engines():
    return SEARCH_ENGINES


@router.get("/spectra-types")
def get_spectra_types():
    return SPECTRA_TYPES


@router.get("/library-formats")
def get_library_formats():
    return LIBRARY_FORMATS


@router.get("/enzymes")
def get_enzymes():
    return ENZYMES


@router.get("/tags")
def get_tags():
    return TAGS


@router.get("/instrument-types")
def get_instrument_types():
    return INSTRUMENT_TYPES


@router.get("/fragmentation-methods")
def get_fragmentation_methods():
    return FRAGMENTATION_METHODS


@router.get("/defaults/{job_type}")
def get_defaults(job_type: str):
    if job_type not in DEFAULTS:
        from fastapi import HTTPException

        raise HTTPException(status_code=404, detail=f"Unknown job type: {job_type}")
    return DEFAULTS[job_type]
