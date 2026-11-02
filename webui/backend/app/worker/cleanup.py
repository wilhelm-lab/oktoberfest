"""
Periodic cleanup of old jobs; the worker runs it every hour (beat_schedule in celery_app).

- RESULT_RETENTION_DAYS > 0: delete jobs that finished (or drafts created) more than that many days ago.
- MAX_DATA_BYTES > 0: when the data folder reaches this size, delete jobs that finished (or drafts created) more than
  CLEANUP_MIN_AGE_HOURS ago, oldest first, until the folder is below the limit again.

Queued and running jobs are never deleted. Both rules are off by default (keep everything).
"""

from __future__ import annotations

import logging
import os
import shutil
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional

from sqlalchemy import func

from app.worker.celery_app import celery_app
from app.config import settings

logger = logging.getLogger(__name__)


def _dir_size(path: Path) -> int:
    """Total size of the files below path in bytes (symlinks are not followed)."""
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += os.lstat(os.path.join(root, name)).st_size
            except OSError:
                pass
    return total


def cleanup_old_jobs(now: Optional[datetime] = None) -> list[str]:
    """Apply the retention and size rules; return the ids of the deleted jobs."""
    from app.db import SessionLocal
    from app.models import Job
    from app.services.storage import storage
    from app.schemas.common import JobStatus, TERMINAL_STATUSES

    now = now or datetime.utcnow()
    deletable = [s.value for s in TERMINAL_STATUSES] + [JobStatus.CREATED.value]  # finished jobs and drafts
    age_reference = func.coalesce(Job.finished_at, Job.created_at)
    deleted: list[str] = []

    def delete(job) -> int:
        job_dir = storage.job_dir(job.id)
        size = _dir_size(job_dir) if job_dir.exists() else 0
        shutil.rmtree(job_dir, ignore_errors=True)
        db.delete(job)
        deleted.append(job.id)
        return size

    db = SessionLocal()
    try:
        candidates = db.query(Job).filter(Job.status.in_(deletable)).order_by(age_reference.asc())

        if settings.result_retention_days > 0:
            cutoff = now - timedelta(days=settings.result_retention_days)
            for job in candidates.filter(age_reference < cutoff).all():
                delete(job)
            db.commit()

        if settings.max_data_bytes > 0:
            used = _dir_size(Path(settings.data_dir))
            if used >= settings.max_data_bytes:
                cutoff = now - timedelta(hours=settings.cleanup_min_age_hours)
                for job in candidates.filter(age_reference < cutoff).all():
                    if used < settings.max_data_bytes:
                        break
                    used -= delete(job)
                db.commit()
                if used >= settings.max_data_bytes:
                    logger.warning(
                        f"Data folder still uses {used / 1e9:.1f} GB (limit {settings.max_data_bytes / 1e9:.1f} GB): "
                        f"the remaining jobs are younger than {settings.cleanup_min_age_hours} h or still running"
                    )
    finally:
        db.close()

    if deleted:
        logger.info(f"Cleanup deleted {len(deleted)} old jobs")
    return deleted


@celery_app.task(name="cleanup_expired_jobs")
def cleanup_expired_jobs() -> int:
    """Hourly cleanup task, see the module docstring."""
    return len(cleanup_old_jobs())
