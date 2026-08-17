"""
Periodic cleanup task for hosted mode (result retention).

Disabled locally (result_retention_days = 0 → keep forever).
To enable: set result_retention_days > 0 AND configure Celery beat schedule.
"""

from __future__ import annotations

import shutil
from datetime import datetime, timedelta

from app.worker.celery_app import celery_app
from app.config import settings


@celery_app.task(name="cleanup_expired_jobs")
def cleanup_expired_jobs():
    """Delete job dirs and DB rows older than result_retention_days.

    No-op if result_retention_days == 0 (default).
    Hosted mode: configure Celery beat to run this task periodically.
    """
    if settings.result_retention_days == 0:
        return

    from app.db import SessionLocal
    from app.models import Job
    from app.services.storage import storage
    from app.schemas.common import TERMINAL_STATUSES

    cutoff = datetime.utcnow() - timedelta(days=settings.result_retention_days)
    db = SessionLocal()
    try:
        expired = (
            db.query(Job)
            .filter(Job.created_at < cutoff)
            .filter(Job.status.in_([s.value for s in TERMINAL_STATUSES]))
            .all()
        )
        for job in expired:
            job_dir = storage.job_dir(job.id)
            if job_dir.exists():
                shutil.rmtree(job_dir, ignore_errors=True)
            db.delete(job)
        db.commit()
    finally:
        db.close()
