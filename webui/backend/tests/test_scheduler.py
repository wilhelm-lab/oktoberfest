from __future__ import annotations

import os
import json
import pytest
from unittest.mock import MagicMock, patch

from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from app.main import app
from app.db import Base, get_db
from app.config import settings
from app.api.deps import get_current_user, AppUser, check_job_access
from app.models import Job
from app.schemas.common import JobStatus
from app.services.scheduler import trigger_scheduler
from app.services import job_service


@pytest.fixture(scope="function")
def tmp_store(tmp_path):
    from app.services.storage import LocalStorage
    ls = LocalStorage(str(tmp_path))
    with (
        patch("app.services.storage.storage", ls),
        patch("app.api.v1.jobs.storage", ls),
        patch("app.api.v1.files.storage", ls),
        patch("app.services.job_service.storage", ls),
        patch("app.worker.tasks.storage", ls),
    ):
        yield tmp_path, ls


@pytest.fixture(scope="function")
def test_db(tmp_path):
    db_url = f"sqlite:///{tmp_path}/test.db"
    engine = create_engine(db_url, connect_args={"check_same_thread": False})
    Base.metadata.create_all(engine)
    TestSession = sessionmaker(autocommit=False, autoflush=False, bind=engine)

    def override_get_db():
        db = TestSession()
        try:
            yield db
        finally:
            db.close()

    app.dependency_overrides[get_db] = override_get_db
    with (
        patch("app.db.SessionLocal", TestSession),
        patch("app.services.scheduler.SessionLocal", TestSession),
        patch("app.worker.tasks.SessionLocal", TestSession),
    ):
        yield TestSession
    app.dependency_overrides.clear()
    engine.dispose()


@pytest.fixture
def client(test_db, tmp_store):
    return TestClient(app, raise_server_exceptions=True)


# Mock Redis for trigger_scheduler
@pytest.fixture(autouse=True)
def mock_redis():
    with patch("app.services.scheduler.redis") as mock_redis_module:
        mock_r = MagicMock()
        mock_lock = MagicMock()
        mock_lock.acquire.return_value = True
        mock_r.lock.return_value = mock_lock
        mock_redis_module.from_url.return_value = mock_r
        yield mock_r, mock_lock


def test_auth_headers_and_role_scoping(client, test_db):
    # Enable hosted mode for this test
    with patch.object(settings, "app_mode", "hosted"):
        with patch.object(settings, "global_users", "admin_user"):
            # 1. Unauthenticated request should fail
            resp = client.post("/api/v1/jobs", json={
                "job_type": "Rescoring",
                "config": {"inputs": {"search_results_type": "Maxquant", "spectra_type": "mzml"}},
            })
            assert resp.status_code == 401

            # 2. Authenticated request as regular user1
            resp = client.post(
                "/api/v1/jobs",
                json={
                    "job_type": "Rescoring",
                    "config": {"inputs": {"search_results_type": "Maxquant", "spectra_type": "mzml"}},
                },
                headers={"X-User-Id": "user1"}
            )
            assert resp.status_code == 201
            job1_id = resp.json()["job_id"]

            # 3. Authenticated request as regular user2
            resp = client.post(
                "/api/v1/jobs",
                json={
                    "job_type": "Rescoring",
                    "config": {"inputs": {"search_results_type": "Maxquant", "spectra_type": "mzml"}},
                },
                headers={"X-User-Id": "user2"}
            )
            assert resp.status_code == 201
            job2_id = resp.json()["job_id"]

            # 4. User1 gets user1's job -> 200
            resp = client.get(f"/api/v1/jobs/{job1_id}", headers={"X-User-Id": "user1"})
            assert resp.status_code == 200

            # 5. User1 gets user2's job -> 404 (isolation)
            resp = client.get(f"/api/v1/jobs/{job2_id}", headers={"X-User-Id": "user1"})
            assert resp.status_code == 404

            # 6. Admin user (global) gets user2's job -> 200 (access to all)
            resp = client.get(f"/api/v1/jobs/{job2_id}", headers={"X-User-Id": "admin_user"})
            assert resp.status_code == 200

            # 7. Listing jobs for user1 should only return user1's job
            resp = client.get("/api/v1/jobs", headers={"X-User-Id": "user1"})
            assert resp.status_code == 200
            jobs = resp.json()
            assert len(jobs) == 1
            assert jobs[0]["job_id"] == job1_id

            # 8. Listing jobs for admin should return all jobs
            resp = client.get("/api/v1/jobs", headers={"X-User-Id": "admin_user"})
            assert resp.status_code == 200
            jobs = resp.json()
            assert len(jobs) == 2

            # 9. User3 in group 1491800192 gets user2's job -> 200 (belongs to global group)
            resp = client.get(
                f"/api/v1/jobs/{job2_id}",
                headers={"X-User-Id": "user3", "X-User-Groups": "some_group, 1491800192, other"}
            )
            assert resp.status_code == 200

            # 10. Listing jobs for user3 in group 1491800192 should return all jobs
            resp = client.get(
                "/api/v1/jobs",
                headers={"X-User-Id": "user3", "X-User-Groups": "1491800192"}
            )
            assert resp.status_code == 200
            jobs = resp.json()
            assert len(jobs) == 2


def test_thread_capping_in_hosted_mode(client):
    with patch.object(settings, "app_mode", "hosted"):
        with patch.object(settings, "max_threads_per_user", 4):
            # Create a job asking for 8 threads
            resp = client.post(
                "/api/v1/jobs",
                json={
                    "job_type": "Rescoring",
                    "config": {
                        "inputs": {"search_results_type": "Maxquant", "spectra_type": "mzml"},
                        "numThreads": 8,
                    },
                },
                headers={"X-User-Id": "user1"}
            )
            assert resp.status_code == 201
            job = resp.json()
            assert job["config"]["numThreads"] == 4  # capped at 4!


def test_scheduler_limits_and_queuing(test_db):
    Session = test_db
    db = Session()

    # Clear any leftover jobs
    db.query(Job).delete()
    db.commit()

    with patch.object(settings, "app_mode", "hosted"):
        with patch.object(settings, "max_threads_per_user", 4):
            with patch.object(settings, "max_total_threads", 6):
                with patch("app.services.scheduler.run_oktoberfest_job") as mock_run_task:
                    mock_delay = MagicMock()
                    mock_delay.id = "fake-task-id"
                    mock_run_task.delay.return_value = mock_delay

                    # 1. Create a job for user1 using 4 threads
                    job1 = Job(
                        id="job1",
                        job_type="Rescoring",
                        status=JobStatus.QUEUED.value,
                        owner_id="user1",
                        config_json=json.dumps({"numThreads": 4})
                    )
                    db.add(job1)

                    # 2. Create another job for user1 using 2 threads (should be skipped because user1 already has a running job)
                    job2 = Job(
                        id="job2",
                        job_type="Rescoring",
                        status=JobStatus.QUEUED.value,
                        owner_id="user1",
                        config_json=json.dumps({"numThreads": 2})
                    )
                    db.add(job2)

                    # 3. Create a job for user2 using 4 threads (would make total threads = 8 > 6, so should be skipped)
                    job3 = Job(
                        id="job3",
                        job_type="Rescoring",
                        status=JobStatus.QUEUED.value,
                        owner_id="user2",
                        config_json=json.dumps({"numThreads": 4})
                    )
                    db.add(job3)

                    # 4. Create a job for user3 using 2 threads (can run! total threads = 4 + 2 = 6 <= 6)
                    job4 = Job(
                        id="job4",
                        job_type="Rescoring",
                        status=JobStatus.QUEUED.value,
                        owner_id="user3",
                        config_json=json.dumps({"numThreads": 2})
                    )
                    db.add(job4)

                    db.commit()

                    # Trigger scheduler
                    trigger_scheduler()

                    # Reload jobs from DB
                    db.refresh(job1)
                    db.refresh(job2)
                    db.refresh(job3)
                    db.refresh(job4)

                    # Assert correct execution states
                    assert job1.status == JobStatus.RUNNING.value  # Ran (FIFO)
                    assert job2.status == JobStatus.QUEUED.value   # Skipped (User1 already running)
                    assert job3.status == JobStatus.QUEUED.value   # Skipped (Would exceed max 6 threads)
                    assert job4.status == JobStatus.RUNNING.value  # Ran (Fits total thread budget)

                    # Total running threads is exactly 6 (4 + 2)
                    assert mock_run_task.delay.call_count == 2
                    mock_run_task.delay.assert_any_call("job1")
                    mock_run_task.delay.assert_any_call("job4")

    db.close()
