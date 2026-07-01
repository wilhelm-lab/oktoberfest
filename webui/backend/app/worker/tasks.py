from __future__ import annotations

import zipfile
from datetime import datetime
from pathlib import Path

from app.worker.celery_app import celery_app
from app.config import settings
from app.db import SessionLocal
from app.schemas.common import JobStatus
from app.services.storage import storage
from app.worker.backends.local import get_execution_backend


def _read_tail(path: Path, n_lines: int = 200) -> str:
    """Read the last n_lines from a file to capture error logs."""
    if not path.exists():
        return ""
    lines = path.read_text(errors="replace").splitlines()
    return "\n".join(lines[-n_lines:])


def _package_results(job_dir: Path, output_dir: Path, config_path: Path, captured_log: Path) -> Path:
    """Zip the output directory, configuration file, and captured logs into a results.zip archive."""
    zip_path = job_dir / "results.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        # Everything in output/
        if output_dir.exists():
            for f in output_dir.rglob("*"):
                if f.is_file():
                    zf.write(f, f.relative_to(job_dir))
        # config.json for reproducibility
        if config_path.exists():
            zf.write(config_path, config_path.name)
        # captured.log
        if captured_log.exists():
            zf.write(captured_log, captured_log.name)
    return zip_path


@celery_app.task(bind=True, name="run_oktoberfest_job", acks_late=True)
def run_oktoberfest_job(self, job_id: str):
    """
    Execute an Oktoberfest job as a Celery background task.
    Updates the job status in the database, runs the execution backend, 
    packages the results, and triggers the scheduler upon completion.
    """
    from app.models import Job

    db = SessionLocal()
    try:
        job = db.get(Job, job_id)
        if job is None:
            raise RuntimeError(f"Job {job_id} not found in DB")

        # Idempotency: if already succeeded, skip
        if job.status == JobStatus.SUCCEEDED.value and storage.results_exists(job_id):
            return

        # Transition to RUNNING
        job.status = JobStatus.RUNNING.value
        job.started_at = datetime.utcnow()
        db.commit()

        job_dir = storage.job_dir(job_id)
        config_path = storage.config_path(job_id)
        output_dir = storage.output_dir(job_id)
        captured_log = storage.captured_log_path(job_id)

        backend = get_execution_backend(settings.execution_backend)
        result = backend.run(job_dir=job_dir, config_path=config_path)

        if result.returncode != 0:
            log_tail = _read_tail(captured_log)
            has_results = False
            try:
                _package_results(job_dir, output_dir, config_path, captured_log)
                has_results = True
            except Exception:
                pass
            job.status = JobStatus.FAILED.value
            job.finished_at = datetime.utcnow()
            job.error = f"Process exited with code {result.returncode}\n\n{log_tail}"
            job.has_results = has_results
            db.commit()
            raise RuntimeError(f"Oktoberfest exited with code {result.returncode}")

        # Package results
        zip_path = _package_results(job_dir, output_dir, config_path, captured_log)

        job.status = JobStatus.SUCCEEDED.value
        job.finished_at = datetime.utcnow()
        job.has_results = True
        db.commit()

    except Exception as exc:
        captured_log = storage.captured_log_path(job_id)
        log_tail = _read_tail(captured_log)
        job = db.get(Job, job_id)
        if job and job.status not in (JobStatus.SUCCEEDED.value, JobStatus.FAILED.value):
            has_results = False
            try:
                _package_results(storage.job_dir(job_id), storage.output_dir(job_id), storage.config_path(job_id), captured_log)
                has_results = True
            except Exception:
                pass
            job.status = JobStatus.FAILED.value
            job.finished_at = datetime.utcnow()
            job.error = f"{type(exc).__name__}: {exc}\n\n{log_tail}"
            job.has_results = has_results
            db.commit()
        raise
    finally:
        db.close()
        try:
            from app.services.scheduler import trigger_scheduler
            trigger_scheduler()
        except Exception:
            pass
