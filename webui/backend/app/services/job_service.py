from __future__ import annotations

import json
import uuid
from datetime import datetime
from typing import Any, Optional

from sqlalchemy.orm import Session

from app.config import settings
from app.models import Job
from app.schemas.common import JobStatus, LEGAL_TRANSITIONS, TERMINAL_STATUSES
from app.services.storage import storage


# Allowed file extensions per upload role
ROLE_EXTENSIONS: dict[str, list[str]] = {
    "search_results": [".txt", ".tsv", ".csv", ".pepxml", ".mzid", ".xml"],
    "spectra": [".raw", ".mzml", ".mzML", ".hdf", ".hdf5", ".d", ".zip"],
    "fasta": [".fasta", ".fa", ".faa"],
    "peptides": [".csv", ".tsv", ".txt"],
}

# Required upload roles per job type
REQUIRED_ROLES: dict[str, list[str]] = {
    "Rescoring": ["search_results", "spectra"],
    "CollisionEnergyCalibration": ["search_results", "spectra"],
    "SpectralLibraryGeneration": ["fasta"],  # or peptides — checked dynamically
}


def create_job(db: Session, job_type: str, config_data: dict[str, Any]) -> Job:
    """Create a new Job row and its filesystem directories."""
    job_id = str(uuid.uuid4())
    job = Job(
        id=job_id,
        job_type=job_type,
        status=JobStatus.CREATED.value,
        config_json=json.dumps(config_data),
        created_at=datetime.utcnow(),
    )
    db.add(job)
    db.commit()
    db.refresh(job)
    storage.create_job_dirs(job_id)
    return job


def transition_status(
    db: Session,
    job: Job,
    new_status: JobStatus,
    **kwargs: Any,
) -> None:
    """Update job status, enforcing legal transitions."""
    current = JobStatus(job.status)
    if new_status not in LEGAL_TRANSITIONS[current]:
        raise ValueError(f"Cannot transition from {current} to {new_status}")
    job.status = new_status.value
    for k, v in kwargs.items():
        setattr(job, k, v)
    db.commit()


def get_job(db: Session, job_id: str) -> Optional[Job]:
    return db.get(Job, job_id)


def list_jobs(db: Session, limit: int = 50, offset: int = 0) -> list[Job]:
    return db.query(Job).order_by(Job.created_at.desc()).offset(offset).limit(limit).all()


def get_required_files(job_type: str, config_data: dict[str, Any]) -> list[dict[str, Any]]:
    """Return expected upload roles for the frontend."""
    result = []
    if job_type in ("Rescoring", "CollisionEnergyCalibration"):
        engine = config_data.get("inputs", {}).get("search_results_type", "Maxquant")
        result.append({"role": "search_results", "accept": ROLE_EXTENSIONS["search_results"], "multiple": True})
        result.append({"role": "spectra", "accept": ROLE_EXTENSIONS["spectra"], "multiple": True})
    elif job_type == "SpectralLibraryGeneration":
        lib_type = config_data.get("inputs", {}).get("library_input_type", "fasta")
        if lib_type in ("fasta", ""):
            result.append({"role": "fasta", "accept": ROLE_EXTENSIONS["fasta"], "multiple": False})
        else:
            result.append({"role": "peptides", "accept": ROLE_EXTENSIONS["peptides"], "multiple": False})
    return result


def check_uploaded_files(job_id: str, job_type: str, config_data: dict[str, Any]) -> list[str]:
    """Return list of missing required roles."""
    inputs_dir = storage.inputs_dir(job_id)
    if not inputs_dir.exists():
        return [r["role"] for r in get_required_files(job_type, config_data)]

    uploaded = {f.name for f in inputs_dir.iterdir() if f.is_file()}
    missing = []

    if job_type in ("Rescoring", "CollisionEnergyCalibration"):
        # Check at least one search_results file and one spectra file
        # We store them with original name; at least one file must exist in either group
        # (simplified: check inputs dir is non-empty)
        if not uploaded:
            missing = ["search_results", "spectra"]
    elif job_type == "SpectralLibraryGeneration":
        if not uploaded:
            missing = ["fasta"]

    return missing


def get_job_progress_phase(job_id: str) -> Optional[str]:
    """Scan the captured.log file to determine the specific sub-phase of the execution."""
    log_path = storage.captured_log_path(job_id)
    if not log_path.exists():
        return None

    try:
        content = log_path.read_text(errors="replace")
    except Exception:
        return None

    # Milestones from last to first
    phases = [
        # Rescoring phases
        ("Finished quantification", "Finishing quantification"),
        ("Starting quantification", "Performing quantification"),
        ("Finished rescoring.", "Generating summary plots"),
        ("Generating summary plots...", "Generating summary plots"),
        ("Finished Generating xiFDR input.", "Plotting / Finishing rescoring"),
        ("Generating xiFDR input.", "Generating xiFDR input"),
        ("Finished rescoring.", "Finishing rescoring"),
        ("Starting rescoring", "Rescoring / FDR Estimation"),
        ("Merging input tab files for rescoring with peptide property prediction", "Preparing files / predicting peptide properties"),
        ("Merging input tab files for rescoring without peptide property prediction", "Preparing files / merging search results"),
        
        # Spectral Library Generation phases
        ("Finished writing the library to disk", "Library generated successfully"),
        ("Writing library", "Writing library to disk"),
        ("Getting predictions", "Getting predictions from model"),
        ("speclib_digested", "Generating internal metadata"),
        ("speclib_created_internal", "Generating internal metadata"),
        
        # CE Calibration / Common phases
        ("Performing RANSAC regression", "Collision energy alignment (RANSAC)"),
        ("Converting search results from", "Converting search results"),
        ("Job executed with the following config:", "Initializing Job"),
    ]

    for marker, label in phases:
        if marker in content:
            return label

    return "Initializing"

