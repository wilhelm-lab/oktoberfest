from __future__ import annotations

import json
import logging
from datetime import datetime
import redis

from app.config import settings
from app.db import SessionLocal
from app.models import Job
from app.schemas.common import JobStatus
from app.worker.tasks import run_oktoberfest_job

logger = logging.getLogger(__name__)


def trigger_scheduler() -> None:
    """Acquires a lock and runs the scheduler to start queued jobs if resources are available."""
    r = redis.from_url(settings.redis_url)
    lock = r.lock("oktoberfest_scheduler_lock", timeout=10)

    # Try to acquire lock
    if not lock.acquire(blocking=True, blocking_timeout=2):
        logger.warning("Could not acquire scheduler lock, skipping scheduling run.")
        return

    db = SessionLocal()
    try:
        # 1. Fetch running jobs to determine current usage
        running_jobs = db.query(Job).filter(Job.status == JobStatus.RUNNING.value).all()

        # Calculate currently running threads and running user IDs
        running_user_ids = set()
        running_threads = 0
        for job in running_jobs:
            if job.owner_id:
                running_user_ids.add(job.owner_id)

            # Extract numThreads from the config_json
            num_threads = 1
            if job.config_json:
                try:
                    cfg = json.loads(job.config_json)
                    num_threads = int(cfg.get("numThreads", 1))
                except Exception:
                    pass
            running_threads += num_threads

        logger.info(
            f"Scheduler state: running_threads={running_threads}/{settings.max_total_threads}, active_users={running_user_ids}"
        )

        # 2. Fetch all queued jobs, ordered by created_at (FIFO)
        queued_jobs = db.query(Job).filter(Job.status == JobStatus.QUEUED.value).order_by(Job.created_at.asc()).all()

        for job in queued_jobs:
            # Check user limit: one job per user at a time
            owner_id = job.owner_id
            if settings.app_mode == "hosted" and owner_id:
                if owner_id in running_user_ids:
                    logger.debug(f"Skipping job {job.id} because user {owner_id} already has a running job.")
                    continue

            # Determine how many threads this job wants
            job_threads = 1
            if job.config_json:
                try:
                    cfg = json.loads(job.config_json)
                    job_threads = int(cfg.get("numThreads", 1))
                except Exception:
                    pass

            # Force max threads per user (4 threads)
            if settings.app_mode == "hosted" and job_threads > settings.max_threads_per_user:
                job_threads = settings.max_threads_per_user

            # Check total thread limit
            if running_threads + job_threads > settings.max_total_threads:
                logger.info(
                    f"Skipping job {job.id} (needs {job_threads} threads) because it would exceed max total threads ({running_threads} + {job_threads} > {settings.max_total_threads})."
                )
                continue

            # If we get here, the job can run!
            logger.info(f"Starting job {job.id} for user {owner_id} using {job_threads} threads.")

            # Start Celery task
            task = run_oktoberfest_job.delay(job.id)

            # Update job status to RUNNING and set Celery task ID
            job.status = JobStatus.RUNNING.value
            job.started_at = datetime.utcnow()
            job.celery_task_id = task.id

            # Add to running sets
            if owner_id:
                running_user_ids.add(owner_id)
            running_threads += job_threads

        db.commit()
    except Exception as exc:
        logger.exception(f"Error during scheduling: {exc}")
        db.rollback()
    finally:
        db.close()
        try:
            lock.release()
        except Exception:
            pass
