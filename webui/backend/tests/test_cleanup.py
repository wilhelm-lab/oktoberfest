"""Hourly cleanup of old jobs: retention days and the data folder size limit."""

from __future__ import annotations

from datetime import datetime, timedelta
from unittest.mock import patch

from app.config import settings
from app.models import Job
from app.schemas.common import JobStatus
from app.worker.cleanup import cleanup_old_jobs
from tests.test_scheduler import test_db, tmp_store  # noqa: F401

NOW = datetime(2026, 11, 3, 12, 0)


def _job(session, store, job_id, status, hours_ago, size=0):
    finished = NOW - timedelta(hours=hours_ago)
    job = Job(id=job_id, job_type="Rescoring", status=status.value, created_at=finished - timedelta(hours=1))
    if status != JobStatus.CREATED:
        job.finished_at = finished if status not in (JobStatus.QUEUED, JobStatus.RUNNING) else None
    session.add(job)
    store.create_job_dirs(job_id)
    (store.inputs_dir(job_id) / "data.bin").write_bytes(b"x" * size)


def _setup(test_db, tmp_store, jobs):
    tmp_path, store = tmp_store
    session = test_db()
    for job in jobs:
        _job(session, store, *job)
    session.commit()
    session.close()
    return tmp_path


def _remaining(test_db):
    session = test_db()
    ids = sorted(j.id for j in session.query(Job).all())
    session.close()
    return ids


def test_nothing_is_deleted_by_default(test_db, tmp_store):
    _setup(test_db, tmp_store, [("old", JobStatus.SUCCEEDED, 24 * 400, 10)])
    assert cleanup_old_jobs(NOW) == []


def test_retention_days(test_db, tmp_store):
    tmp_path = _setup(test_db, tmp_store, [
        ("old", JobStatus.SUCCEEDED, 24 * 31),
        ("old_draft", JobStatus.CREATED, 24 * 31),
        ("old_running", JobStatus.RUNNING, 24 * 31),
        ("recent", JobStatus.FAILED, 24 * 29),
    ])
    with patch.object(settings, "result_retention_days", 30):
        assert sorted(cleanup_old_jobs(NOW)) == ["old", "old_draft"]
    assert _remaining(test_db) == ["old_running", "recent"]
    assert not (tmp_path / "old").exists()


def test_size_limit_deletes_oldest_jobs_older_than_a_day(test_db, tmp_store):
    mb = 1_000_000  # job sizes well above the few KB of the test database in the same folder
    _setup(test_db, tmp_store, [
        ("oldest", JobStatus.SUCCEEDED, 72, mb),
        ("older", JobStatus.CANCELLED, 48, mb),
        ("today", JobStatus.SUCCEEDED, 2, mb),
        ("queued", JobStatus.QUEUED, 100, mb),
    ])
    with patch.object(settings, "data_dir", str(tmp_store[0])):
        # 4 MB: deleting the oldest finished job is enough to get below 3.5 MB
        with patch.object(settings, "max_data_bytes", 3.5 * mb):
            assert cleanup_old_jobs(NOW) == ["oldest"]
        # still above: only jobs older than a day are deleted, never queued or running ones
        with patch.object(settings, "max_data_bytes", 0.5 * mb):
            assert cleanup_old_jobs(NOW) == ["older"]
    assert _remaining(test_db) == ["queued", "today"]


def test_size_limit_not_reached(test_db, tmp_store):
    _setup(test_db, tmp_store, [("old", JobStatus.SUCCEEDED, 72, 1_000_000)])
    with patch.object(settings, "max_data_bytes", 10_000_000), patch.object(settings, "data_dir", str(tmp_store[0])):
        assert cleanup_old_jobs(NOW) == []
