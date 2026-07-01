from __future__ import annotations

import json
from datetime import datetime
from typing import Any, Optional

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel
from sqlalchemy.orm import Session

from app.api.deps import get_current_user, check_quota, AppUser, check_job_access
from app.config import settings
from app.db import get_db
from app.models import Job
from app.schemas.common import JobStatus, TERMINAL_STATUSES
from app.schemas.rescoring import RescoringConfig
from app.schemas.ce_calibration import CeCalibrationConfig
from app.schemas.speclib import SpeclibConfig
from app.services import job_service
from app.services.storage import storage
from app.worker.tasks import run_oktoberfest_job

router = APIRouter(tags=["jobs"])

# ── Request / Response models ──────────────────────────────────────────────────


class CreateJobRequest(BaseModel):
    job_type: str
    config: dict[str, Any]


class JobResponse(BaseModel):
    job_id: str
    job_type: str
    status: str
    created_at: Optional[datetime]
    started_at: Optional[datetime]
    finished_at: Optional[datetime]
    has_results: bool
    error: Optional[str]
    config: Optional[dict[str, Any]]
    required_files: Optional[list[dict[str, Any]]] = None
    progress_phase: Optional[str] = None

    model_config = {"from_attributes": True}


def _job_to_response(job: Job, include_required: bool = False) -> dict[str, Any]:
    d = {
        "job_id": job.id,
        "job_type": job.job_type,
        "status": job.status,
        "created_at": job.created_at,
        "started_at": job.started_at,
        "finished_at": job.finished_at,
        "has_results": job.has_results,
        "error": job.error,
        "config": json.loads(job.config_json) if job.config_json else None,
        "progress_phase": job_service.get_job_progress_phase(job.id) if job.status in ("RUNNING", "SUCCEEDED", "FAILED") else None,
    }
    if include_required:
        cfg = json.loads(job.config_json) if job.config_json else {}
        d["required_files"] = job_service.get_required_files(job.job_type, cfg)
    return d


def _validate_config(job_type: str, config: dict[str, Any]):
    """Parse and validate config against the matching Pydantic model."""
    try:
        if job_type == "Rescoring":
            return RescoringConfig(**config)
        elif job_type == "CollisionEnergyCalibration":
            return CeCalibrationConfig(**config)
        elif job_type == "SpectralLibraryGeneration":
            return SpeclibConfig(**config)
        else:
            raise HTTPException(status_code=422, detail=f"Unknown job_type: {job_type}")
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc))


# ── Endpoints ──────────────────────────────────────────────────────────────────


@router.post("", status_code=201)
async def create_job(
    body: CreateJobRequest,
    db: Session = Depends(get_db),
    user: AppUser = Depends(get_current_user),
):
    """
    Create a new job with the specified type and configuration.
    The job starts in the CREATED status and requires files to be uploaded before submission.
    """
    allowed = {"Rescoring", "CollisionEnergyCalibration", "SpectralLibraryGeneration"}
    if body.job_type not in allowed:
        raise HTTPException(status_code=422, detail=f"job_type must be one of {sorted(allowed)}")

    await check_quota(user)

    # Validate config (file paths not yet resolved — that's ok)
    validated = _validate_config(body.job_type, body.config)

    job = job_service.create_job(db, body.job_type, validated.model_dump(exclude_none=True))
    # Store owner_id for hosted-mode scoping
    job.owner_id = user.id if settings.app_mode == "hosted" else None
    db.commit()
    return _job_to_response(job, include_required=True)


@router.post("/{job_id}/submit", status_code=202)
def submit_job(
    job_id: str,
    db: Session = Depends(get_db),
    user: AppUser = Depends(get_current_user),
):
    """
    Submit a previously created job for execution.
    Validates that all required files have been uploaded and transitions the job to QUEUED.
    """
    job = job_service.get_job(db, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    check_job_access(job, user)
    if job.status != JobStatus.CREATED.value:
        raise HTTPException(status_code=409, detail=f"Job is in status {job.status}, cannot submit")

    cfg_raw = json.loads(job.config_json) if job.config_json else {}

    # Check required files exist
    missing = job_service.check_uploaded_files(job_id, job.job_type, cfg_raw)
    if missing:
        raise HTTPException(status_code=422, detail=f"Missing required uploads: {missing}")

    # Resolve uploaded file paths into config
    inputs_dir = storage.inputs_dir(job_id)
    uploaded_files = list(inputs_dir.iterdir()) if inputs_dir.exists() else []
    output_dir = str(storage.output_dir(job_id))

    if "inputs" not in cfg_raw:
        cfg_raw["inputs"] = {}

    if job.job_type in ("Rescoring", "CollisionEnergyCalibration"):
        # search_results: first uploaded file matching typical search result extensions
        search_exts = {".txt", ".tsv", ".csv", ".pepxml", ".mzid", ".xml"}
        spectra_exts = {".raw", ".mzml", ".hdf", ".hdf5", ".d"}
        search_files = [f for f in uploaded_files if f.suffix.lower() in search_exts]
        spectra_files = [f for f in uploaded_files if f.suffix.lower() in spectra_exts]
        # If a file doesn't match either, assign by position
        if not search_files and not spectra_files and uploaded_files:
            search_files = [uploaded_files[0]]
            spectra_files = uploaded_files[1:] or [uploaded_files[0]]
        if search_files:
            cfg_raw["inputs"]["search_results"] = str(search_files[0])
        if spectra_files:
            cfg_raw["inputs"]["spectra"] = str(spectra_files[0]) if len(spectra_files) == 1 else str(inputs_dir)
    elif job.job_type == "SpectralLibraryGeneration":
        if uploaded_files:
            cfg_raw["inputs"]["library_input"] = str(uploaded_files[0])

    # Final validation with resolved paths
    validated = _validate_config(job.job_type, cfg_raw)

    # Serialize to Oktoberfest JSON shape
    config_dict = validated.to_oktoberfest_dict(job_id, output_dir)

    # Write config.json
    config_path = storage.config_path(job_id)
    with open(config_path, "w") as f:
        json.dump(config_dict, f, indent=2)

    # Transition to QUEUED
    job_service.transition_status(
        db,
        job,
        JobStatus.QUEUED,
        config_json=json.dumps(config_dict),
    )

    # Trigger scheduler to start jobs
    from app.services.scheduler import trigger_scheduler
    trigger_scheduler()

    return {"job_id": job_id, "status": "QUEUED"}


@router.get("/{job_id}")
def get_job(
    job_id: str,
    db: Session = Depends(get_db),
    user: AppUser = Depends(get_current_user),
):
    """Retrieve details and status of a specific job by its ID."""
    job = job_service.get_job(db, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    check_job_access(job, user)
    return _job_to_response(job)


@router.get("/{job_id}/log")
def get_job_log(
    job_id: str,
    db: Session = Depends(get_db),
    user: AppUser = Depends(get_current_user),
):
    """Retrieve the captured stdout/stderr log of the job execution."""
    job = job_service.get_job(db, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    check_job_access(job, user)
    log_path = storage.captured_log_path(job_id)
    if not log_path.exists():
        return {"log": ""}
    return {"log": log_path.read_text(errors="replace")}


@router.get("/{job_id}/results")
def download_results(
    job_id: str,
    db: Session = Depends(get_db),
    user: AppUser = Depends(get_current_user),
):
    """Download the final results ZIP archive if the job has successfully finished."""
    job = job_service.get_job(db, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    check_job_access(job, user)
    if not storage.results_exists(job_id):
        raise HTTPException(status_code=404, detail="Results not ready")
    from fastapi.responses import FileResponse

    return FileResponse(
        path=str(storage.results_path(job_id)),
        media_type="application/zip",
        filename=f"results_{job_id}.zip",
    )


@router.get("")
def list_jobs(
    limit: int = 50,
    offset: int = 0,
    db: Session = Depends(get_db),
    user: AppUser = Depends(get_current_user),
):
    """
    List jobs with pagination.
    In hosted mode, non-admin users only see their own jobs.
    """
    if settings.app_mode == "hosted" and not user.is_global:
        jobs = job_service.list_jobs_for_user(db, user_id=user.id, limit=limit, offset=offset)
    else:
        jobs = job_service.list_jobs(db, limit=limit, offset=offset)
    return [_job_to_response(j) for j in jobs]


@router.post("/{job_id}/cancel", status_code=200)
def cancel_job(
    job_id: str,
    db: Session = Depends(get_db),
    user: AppUser = Depends(get_current_user),
):
    """
    Cancel an active job.
    Attempts to revoke the Celery task and marks the job status as CANCELLED.
    """
    job = job_service.get_job(db, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    check_job_access(job, user)

    if job.status in (s.value for s in TERMINAL_STATUSES):
        raise HTTPException(status_code=409, detail=f"Job already in terminal status {job.status}")

    # Revoke Celery task (best-effort)
    if job.celery_task_id:
        try:
            from app.worker.celery_app import celery_app

            celery_app.control.revoke(job.celery_task_id, terminate=True)
        except Exception:
            pass

    job.status = JobStatus.CANCELLED.value
    db.commit()

    # Trigger scheduler to start waiting jobs
    from app.services.scheduler import trigger_scheduler
    trigger_scheduler()

    return {"job_id": job_id, "status": "CANCELLED"}
