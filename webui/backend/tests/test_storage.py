"""Phase 1 tests: LocalStorage and path safety."""

from __future__ import annotations

import io
import pytest
from pathlib import Path

from app.services.storage import LocalStorage, _safe_filename


@pytest.fixture
def tmp_storage(tmp_path):
    return LocalStorage(data_dir=str(tmp_path))


def test_create_job_dirs(tmp_storage, tmp_path):
    tmp_storage.create_job_dirs("abc123")
    assert (tmp_path / "abc123" / "inputs").is_dir()
    assert (tmp_path / "abc123" / "output").is_dir()


def test_save_upload_roundtrip(tmp_storage, tmp_path):
    tmp_storage.create_job_dirs("job1")
    data = io.BytesIO(b"hello world")
    dest = tmp_storage.save_upload("job1", "test.txt", data)
    assert dest.exists()
    assert dest.read_bytes() == b"hello world"


def test_results_not_exists_before_creation(tmp_storage):
    tmp_storage.create_job_dirs("job2")
    assert not tmp_storage.results_exists("job2")


def test_results_exists_after_creation(tmp_storage, tmp_path):
    tmp_storage.create_job_dirs("job3")
    results = tmp_storage.results_path("job3")
    results.write_bytes(b"fake zip")
    assert tmp_storage.results_exists("job3")


def test_safe_filename_strips_path():
    assert _safe_filename("../../../etc/passwd") == "passwd"
    assert _safe_filename("foo bar.txt") == "foo_bar.txt"
    assert _safe_filename("normal_file.mzML") == "normal_file.mzML"


def test_path_safety_rejects_traversal(tmp_storage, tmp_path):
    tmp_storage.create_job_dirs("job4")
    data = io.BytesIO(b"bad")
    # A traversal attempt via filename is neutralized by _safe_filename
    # so the save should succeed (safe name is used)
    dest = tmp_storage.save_upload("job4", "../../evil.txt", data)
    assert str(tmp_path / "job4") in str(dest)
