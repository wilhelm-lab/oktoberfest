"""Phase 4 tests: Celery worker task with execution backend mocked."""

from __future__ import annotations

import json
import os
import zipfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

os.environ.setdefault("DATABASE_URL", "sqlite:///:memory:")
os.environ.setdefault("REDIS_URL", "redis://localhost:6379/0")
os.environ.setdefault("DATA_DIR", "/tmp/okt_worker_test")

from app.db import Base
from app.models import Job
from app.schemas.common import JobStatus
from app.worker.backends.base import ExecutionResult


@pytest.fixture
def tmp_db(tmp_path):
    db_url = f"sqlite:///{tmp_path}/test.db"
    engine = create_engine(db_url, connect_args={"check_same_thread": False})
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session, tmp_path


def _make_job(session_factory, job_id: str, tmp_path: Path):
    db = session_factory()
    job = Job(
        id=job_id,
        job_type="SpectralLibraryGeneration",
        status=JobStatus.QUEUED.value,
        config_json=json.dumps({"type": "SpectralLibraryGeneration"}),
    )
    db.add(job)
    db.commit()
    db.close()
    (tmp_path / job_id / "inputs").mkdir(parents=True)
    (tmp_path / job_id / "output").mkdir(parents=True)
    (tmp_path / job_id / "config.json").write_text('{"type":"SpectralLibraryGeneration"}')


def _run_task(job_id: str, Session, fake_storage, mock_backend):
    """Call the Celery task body directly using .run() (bypasses broker)."""
    import app.worker.tasks as task_module

    with (
        patch.object(task_module, "SessionLocal", Session),
        patch.object(task_module, "storage", fake_storage),
        patch.object(task_module, "get_execution_backend", return_value=mock_backend),
    ):
        # For bind=True tasks, .run is already bound — just pass job_id
        task_module.run_oktoberfest_job.run(job_id)


def test_task_success_path(tmp_db):
    Session, tmp_path = tmp_db
    job_id = "success-job"
    _make_job(Session, job_id, tmp_path)

    from app.services.storage import LocalStorage

    fake_storage = LocalStorage(str(tmp_path))
    (tmp_path / job_id / "output" / "result.msp").write_text("library")

    mock_backend = MagicMock()
    mock_backend.run.return_value = ExecutionResult(returncode=0, stdout="done", stderr="")

    _run_task(job_id, Session, fake_storage, mock_backend)

    db = Session()
    job = db.get(Job, job_id)
    assert job.status == JobStatus.SUCCEEDED.value
    assert job.has_results is True
    assert job.finished_at is not None
    db.close()

    zip_path = tmp_path / job_id / "results.zip"
    assert zip_path.exists()
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
    assert any("result.msp" in n for n in names)
    assert "config.json" in names


def test_task_failure_path(tmp_db):
    Session, tmp_path = tmp_db
    job_id = "fail-job"
    _make_job(Session, job_id, tmp_path)
    (tmp_path / job_id / "captured.log").write_text("Error: something went wrong\n")

    from app.services.storage import LocalStorage

    fake_storage = LocalStorage(str(tmp_path))

    mock_backend = MagicMock()
    mock_backend.run.return_value = ExecutionResult(returncode=1, stdout="", stderr="error")

    with pytest.raises(RuntimeError):
        _run_task(job_id, Session, fake_storage, mock_backend)

    db = Session()
    job = db.get(Job, job_id)
    assert job.status == JobStatus.FAILED.value
    assert "exited with code 1" in (job.error or "")
    assert "something went wrong" in (job.error or "")
    db.close()


def test_task_idempotency(tmp_db):
    """Re-running a SUCCEEDED job is a no-op."""
    Session, tmp_path = tmp_db
    job_id = "idempotent-job"
    _make_job(Session, job_id, tmp_path)

    db = Session()
    job = db.get(Job, job_id)
    job.status = JobStatus.SUCCEEDED.value
    job.has_results = True
    db.commit()
    db.close()

    (tmp_path / job_id / "results.zip").write_bytes(b"PK existing")

    from app.services.storage import LocalStorage

    fake_storage = LocalStorage(str(tmp_path))
    mock_backend = MagicMock()

    _run_task(job_id, Session, fake_storage, mock_backend)

    mock_backend.run.assert_not_called()

    db = Session()
    job = db.get(Job, job_id)
    assert job.status == JobStatus.SUCCEEDED.value
    db.close()
