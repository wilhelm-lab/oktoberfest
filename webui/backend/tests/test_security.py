"""Protections against users reaching files, hosts or resources they should not."""

from __future__ import annotations

from unittest.mock import patch

from app.config import settings
from tests.test_scheduler import RESCORING_JOB, _as, client, hosted, mock_redis, test_db, tmp_store  # noqa: F401


def _create(client, session="user1"):
    resp = _as(client, session).post("/api/v1/jobs", json=RESCORING_JOB)
    assert resp.status_code == 201, resp.text
    return resp.json()["job_id"]


def _upload(client, job_id, role, name, data=b"x"):
    return client.post(f"/api/v1/jobs/{job_id}/files", data={"role": role}, files={"file": (name, data)})


def test_submit_ignores_input_paths_from_the_config(client):
    job_id = _create(client)
    assert _upload(client, job_id, "spectra", "run.mzML").status_code == 200
    # only spectra uploaded: a search_results path in the config must not be used (e.g. another job's files)
    config = {"inputs": {"search_results_type": "Maxquant", "search_results": "/etc/passwd"}}
    resp = client.post(f"/api/v1/jobs/{job_id}/submit", json={"config": config})
    assert resp.status_code == 422
    assert "search_results" in resp.json()["detail"]


def test_prediction_server_is_restricted_in_hosted_mode(client, hosted):
    job = {**RESCORING_JOB, "config": {**RESCORING_JOB["config"], "prediction_server": "10.0.0.1:22"}}
    assert _as(client, "user1").post("/api/v1/jobs", json=job).status_code == 422
    with patch.object(settings, "allowed_prediction_servers", "koina.example.org:443"):
        job["config"]["prediction_server"] = "koina.example.org:443"
        assert _as(client, "user1").post("/api/v1/jobs", json=job).status_code == 201


def test_open_jobs_per_client_are_limited(client, hosted):
    with patch.object(settings, "max_open_jobs_per_client", 2):
        _create(client, "user1")
        _create(client, "user1")
        # a new session from the same address counts too
        assert _as(client, "user2").post("/api/v1/jobs", json=RESCORING_JOB).status_code == 429


def test_upload_size_per_job_is_limited(client, hosted):
    job_id = _create(client)
    with patch.object(settings, "max_job_upload_bytes", 10):
        assert _upload(client, job_id, "search_results", "msms.txt", b"x" * 6).status_code == 200
        assert _upload(client, job_id, "spectra", "run.mzML", b"x" * 6).status_code == 413
        # replacing a file does not count its old size
        assert _upload(client, job_id, "search_results", "msms.txt", b"x" * 8).status_code == 200


def test_cors_only_allows_local_pages(client):
    def allowed_origin(origin):
        resp = client.get("/api/v1/health", headers={"Origin": origin})
        return resp.headers.get("access-control-allow-origin")

    assert allowed_origin("http://localhost:5173") == "http://localhost:5173"
    assert allowed_origin("https://evil.example.org") is None


def test_client_address_is_the_one_added_by_the_proxy(client, hosted):
    with patch.object(settings, "max_open_jobs_per_client", 1):
        headers = {"X-Forwarded-For": "1.1.1.1, 10.0.0.5"}
        assert _as(client, "user1").post("/api/v1/jobs", json=RESCORING_JOB, headers=headers).status_code == 201
        # a forged first entry and a new session do not give a new client
        headers = {"X-Forwarded-For": "2.2.2.2, 10.0.0.5"}
        assert _as(client, "user2").post("/api/v1/jobs", json=RESCORING_JOB, headers=headers).status_code == 429


def test_admin_sign_in(client, hosted):
    key = "k" * 48
    with patch.object(settings, "global_users", key):
        _create(client, "user1")
        assert _as(client).post("/api/v1/session/admin", json={"key": "wrong"}).status_code == 403
        assert _as(client).post("/api/v1/session/admin", json={"key": ""}).status_code == 403

        resp = _as(client, "user2").post("/api/v1/session/admin", json={"key": key})
        assert resp.status_code == 200
        assert f"oktoberfest_session={key}" in resp.headers["set-cookie"]
        assert _as(client, key).get("/api/v1/session").json()["is_admin"] is True
        assert len(_as(client, key).get("/api/v1/jobs").json()) == 1  # sees user1's job

        resp = _as(client, key).post("/api/v1/session/signout")
        assert key not in resp.headers["set-cookie"]
        assert _as(client, "user2").get("/api/v1/session").json()["is_admin"] is False


def test_admin_sign_in_without_admins_or_in_local_mode(client, hosted):
    assert _as(client).post("/api/v1/session/admin", json={"key": ""}).status_code == 403
    with patch.object(settings, "app_mode", "local"):
        assert _as(client).post("/api/v1/session/admin", json={"key": "x"}).status_code == 400
