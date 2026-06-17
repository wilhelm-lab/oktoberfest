"""
End-to-end integration test: runs a real SpectralLibraryGeneration job using
a tiny FASTA fixture and the LocalExecutionBackend (no Redis needed — we call
the task body directly).

Mark: slow / requires-oktoberfest. Skip in CI unless ENABLE_E2E=1.
"""

from __future__ import annotations

import json
import os
import zipfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("ENABLE_E2E") != "1", reason="E2E test requires ENABLE_E2E=1 and oktoberfest installed"
)

TINY_FASTA = """>sp|P12345|TINY_HUMAN Tiny test protein OS=Homo sapiens
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQNTIVKQVKQGIGEFYISGQYKGFEFMLKPTPKIVHKNKVLKIDAYISGKIPARIVDSTTNHSQGLKIMNETLQEAGTKDIVDNRFGKALNQYLDSDPEYIQKLMEHYQLNRDQKDVLGGIFKFLDQIFQSQEKKLNGQELDFEPSVTTDEKAAQDMWKQVTLNHQFMKDLESQ
"""


@pytest.fixture
def tiny_fasta(tmp_path: Path) -> Path:
    f = tmp_path / "tiny.fasta"
    f.write_text(TINY_FASTA)
    return f


def test_speclib_e2e(tiny_fasta: Path, tmp_path: Path):
    """Full pipeline: config.json → run_job → results.zip (non-empty)."""
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker
    from app.db import Base
    from app.models import Job
    from app.schemas.common import JobStatus
    from app.services.storage import LocalStorage
    from app.schemas.speclib import SpeclibConfig

    # Setup test DB
    engine = create_engine(f"sqlite:///{tmp_path}/test_e2e.db", connect_args={"check_same_thread": False})
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)

    storage = LocalStorage(str(tmp_path))
    job_id = "e2e-speclib-job"
    storage.create_job_dirs(job_id)

    # Copy FASTA into inputs
    import shutil

    shutil.copy(tiny_fasta, storage.inputs_dir(job_id) / "tiny.fasta")

    # Build and write config
    cfg = SpeclibConfig(
        inputs={
            "library_input_type": "fasta",
            "library_input": str(storage.inputs_dir(job_id) / "tiny.fasta"),
        },
        spectralLibraryOptions={
            "fragmentation": "HCD",
            "collisionEnergy": 30,
            "precursorCharge": [2, 3],
            "minIntensity": 5e-4,
            "nrOx": 1,
            "batchsize": 100,
            "format": "msp",
        },
    )
    config_dict = cfg.to_oktoberfest_dict(job_id, str(storage.output_dir(job_id)))
    config_path = storage.config_path(job_id)
    with open(config_path, "w") as f:
        json.dump(config_dict, f, indent=2)

    # Insert job into DB as QUEUED
    db = Session()
    job = Job(
        id=job_id,
        job_type="SpectralLibraryGeneration",
        status=JobStatus.QUEUED.value,
        config_json=json.dumps(config_dict),
    )
    db.add(job)
    db.commit()
    db.close()

    # Run via task body
    import app.worker.tasks as task_module

    with (
        patch.object(task_module, "SessionLocal", Session),
        patch.object(task_module, "storage", storage),
    ):
        task_module.run_oktoberfest_job.run(job_id)

    # Verify results
    db = Session()
    job = db.get(Job, job_id)
    assert job.status == JobStatus.SUCCEEDED.value, f"Job failed: {job.error}"
    assert job.has_results is True
    db.close()

    zip_path = storage.results_path(job_id)
    assert zip_path.exists()
    assert zip_path.stat().st_size > 100, "results.zip should be non-trivial"

    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
    print(f"results.zip contents: {names}")
    assert "config.json" in names
    assert any(n.endswith(".msp") or "output" in n for n in names), f"No .msp file in {names}"
