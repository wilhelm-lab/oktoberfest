"""Phase 3 API tests using FastAPI TestClient with Celery mocked."""

from __future__ import annotations

import json
import os
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

os.environ.setdefault("DATABASE_URL", "sqlite:///:memory:")
os.environ.setdefault("REDIS_URL", "redis://localhost:6379/0")

from app.main import app
from app.db import Base, get_db
from app.config import settings


@pytest.fixture(scope="function")
def tmp_store(tmp_path):
    """Patch the storage singleton globally for this test."""
    from app.services.storage import LocalStorage

    ls = LocalStorage(str(tmp_path))
    with (
        patch("app.services.storage.storage", ls),
        patch("app.api.v1.jobs.storage", ls),
        patch("app.api.v1.files.storage", ls),
        patch("app.services.job_service.storage", ls),
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
    yield tmp_path
    app.dependency_overrides.clear()
    engine.dispose()


@pytest.fixture
def client(test_db, tmp_store):
    return TestClient(app, raise_server_exceptions=True)


@pytest.fixture
def tmp_path_combined(test_db, tmp_store):
    """Return the storage tmp_path (same for both fixtures when we control it)."""
    store_path, ls = tmp_store
    return store_path


def test_health(client):
    resp = client.get("/api/v1/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_meta_models(client):
    resp = client.get("/api/v1/meta/models")
    assert resp.status_code == 200
    data = resp.json()
    assert "intensity" in data
    assert len(data["intensity"]) > 0


def test_meta_defaults_rescoring(client):
    resp = client.get("/api/v1/meta/defaults/Rescoring")
    assert resp.status_code == 200
    assert resp.json()["type"] == "Rescoring"


def test_meta_defaults_unknown(client):
    resp = client.get("/api/v1/meta/defaults/UnknownType")
    assert resp.status_code == 404


def test_create_rescoring_job(client):
    body = {
        "job_type": "Rescoring",
        "config": {
            "inputs": {"search_results_type": "Maxquant", "spectra_type": "mzml"},
            "fdr_estimation_method": "mokapot",
        },
    }
    resp = client.post("/api/v1/jobs", json=body)
    assert resp.status_code == 201
    data = resp.json()
    assert "job_id" in data
    assert data["status"] == "CREATED"
    assert data["required_files"] is not None


def test_create_job_invalid_type(client):
    resp = client.post("/api/v1/jobs", json={"job_type": "Fake", "config": {}})
    assert resp.status_code == 422


def test_create_job_invalid_fdr(client):
    resp = client.post(
        "/api/v1/jobs",
        json={"job_type": "Rescoring", "config": {"fdr_estimation_method": "bad_value"}},
    )
    assert resp.status_code == 422


def test_submit_without_files_returns_422(client):
    body = {"job_type": "SpectralLibraryGeneration", "config": {"inputs": {"library_input_type": "fasta"}}}
    resp = client.post("/api/v1/jobs", json=body)
    job_id = resp.json()["job_id"]
    resp = client.post(f"/api/v1/jobs/{job_id}/submit")
    assert resp.status_code == 422


def test_full_speclib_happy_path(client, tmp_store):
    """POST job → POST file → POST submit (mocked) → GET status."""
    store_path, ls = tmp_store

    # 1. Create job
    body = {
        "job_type": "SpectralLibraryGeneration",
        "config": {
            "inputs": {"library_input_type": "fasta"},
        },
    }
    resp = client.post("/api/v1/jobs", json=body)
    assert resp.status_code == 201
    job_id = resp.json()["job_id"]

    # 2. Upload a fake FASTA file
    fake_fasta = b">sp|P12345|PROT Human\nMKSLIVALLVGAVGA\n"
    resp = client.post(
        f"/api/v1/jobs/{job_id}/files",
        data={"role": "fasta"},
        files={"file": ("test.fasta", fake_fasta, "text/plain")},
    )
    assert resp.status_code == 200, resp.text
    assert resp.json()["role"] == "fasta"

    # 3. Submit (mock Celery task's delay)
    with patch("app.api.v1.jobs.run_oktoberfest_job") as mock_task:
        mock_delay = MagicMock()
        mock_delay.id = "fake-celery-task-id"
        mock_task.delay.return_value = mock_delay
        resp = client.post(f"/api/v1/jobs/{job_id}/submit")
    assert resp.status_code == 202, resp.text
    assert resp.json()["status"] == "QUEUED"
    mock_task.delay.assert_called_once_with(job_id)

    # 4. Check config.json was written in the storage tmp dir
    config_path = store_path / job_id / "config.json"
    assert config_path.exists(), f"config.json missing at {config_path}"
    config = json.loads(config_path.read_text())
    assert config["type"] == "SpectralLibraryGeneration"
    assert "library_input" in config["inputs"]
    assert config["output"] == str(store_path / job_id / "output")

    # 5. GET status
    resp = client.get(f"/api/v1/jobs/{job_id}")
    assert resp.status_code == 200
    assert resp.json()["status"] == "QUEUED"


def test_results_404_before_completion(client):
    body = {"job_type": "SpectralLibraryGeneration", "config": {"inputs": {"library_input_type": "fasta"}}}
    resp = client.post("/api/v1/jobs", json=body)
    job_id = resp.json()["job_id"]
    resp = client.get(f"/api/v1/jobs/{job_id}/results")
    assert resp.status_code == 404


def test_results_200_when_zip_present(client, tmp_store):
    store_path, ls = tmp_store
    body = {"job_type": "SpectralLibraryGeneration", "config": {"inputs": {"library_input_type": "fasta"}}}
    resp = client.post("/api/v1/jobs", json=body)
    job_id = resp.json()["job_id"]
    # Simulate completed results by writing results.zip via the storage instance
    ls.results_path(job_id).write_bytes(b"PK fake zip")
    resp = client.get(f"/api/v1/jobs/{job_id}/results")
    assert resp.status_code == 200


def test_list_jobs(client):
    resp = client.get("/api/v1/jobs")
    assert resp.status_code == 200
    assert isinstance(resp.json(), list)


def test_double_submit_returns_409(client, tmp_store):
    store_path, ls = tmp_store
    body = {"job_type": "SpectralLibraryGeneration", "config": {"inputs": {"library_input_type": "fasta"}}}
    resp = client.post("/api/v1/jobs", json=body)
    job_id = resp.json()["job_id"]

    # Upload file
    client.post(
        f"/api/v1/jobs/{job_id}/files",
        data={"role": "fasta"},
        files={"file": ("x.fasta", b">P\nM\n", "text/plain")},
    )

    with patch("app.api.v1.jobs.run_oktoberfest_job") as mock_task:
        mock_delay = MagicMock()
        mock_delay.id = "t1"
        mock_task.delay.return_value = mock_delay
        resp1 = client.post(f"/api/v1/jobs/{job_id}/submit")
    assert resp1.status_code == 202

    resp2 = client.post(f"/api/v1/jobs/{job_id}/submit")
    assert resp2.status_code == 409
