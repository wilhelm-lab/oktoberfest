from __future__ import annotations

import os
from pathlib import Path

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile
from sqlalchemy.orm import Session

from app.api.deps import get_current_user, AppUser, check_job_access
from app.config import settings
from app.db import get_db
from app.schemas.common import JobStatus
from app.services import job_service
from app.services.storage import storage, _safe_filename
from app.services.job_service import ROLE_EXTENSIONS as _ROLE_EXTENSIONS

router = APIRouter(tags=["files"])


@router.post("/{job_id}/files")
async def upload_file(
    job_id: str,
    role: str = Form(...),
    file: UploadFile = File(...),
    db: Session = Depends(get_db),
    user: AppUser = Depends(get_current_user),
):
    job = job_service.get_job(db, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    check_job_access(job, user)
    if job.status != JobStatus.CREATED.value:
        raise HTTPException(status_code=409, detail=f"Job is {job.status}; cannot upload files")

    allowed_roles = {"search_results", "spectra", "fasta", "peptides"}
    if role not in allowed_roles:
        raise HTTPException(status_code=422, detail=f"role must be one of {sorted(allowed_roles)}")

    # Extension check
    ext = Path(file.filename or "").suffix.lower()
    allowed_exts = _ROLE_EXTENSIONS.get(role, [])
    if allowed_exts and ext not in [e.lower() for e in allowed_exts]:
        raise HTTPException(
            status_code=422,
            detail=f"Extension {ext!r} not allowed for role {role!r}. Allowed: {allowed_exts}",
        )

    # Size cap — stream to disk, counting bytes
    safe_name = _safe_filename(file.filename or "upload")
    dest = storage.inputs_dir(job_id) / safe_name
    dest.parent.mkdir(parents=True, exist_ok=True)

    # Path safety
    try:
        dest.resolve().relative_to(storage.inputs_dir(job_id).resolve())
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid filename (path traversal)")

    total = 0
    with open(dest, "wb") as f:
        while chunk := await file.read(1 << 20):  # 1 MB chunks
            total += len(chunk)
            if total > settings.max_upload_bytes:
                dest.unlink(missing_ok=True)
                raise HTTPException(status_code=413, detail="File exceeds size limit")
            f.write(chunk)

    return {
        "job_id": job_id,
        "role": role,
        "filename": safe_name,
        "stored_path": str(dest),
        "size": total,
    }
