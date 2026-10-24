from __future__ import annotations

import os
import re
from pathlib import Path
from typing import BinaryIO, Protocol

from app.config import settings


def _safe_filename(name: str) -> str:
    """Sanitize a filename: keep only alphanumerics, dots, dashes, underscores."""
    name = os.path.basename(name)
    name = re.sub(r"[^\w.\-]", "_", name)
    return name or "upload"


def _assert_within_dir(path: Path, root: Path, job_id: str) -> None:
    """Raise ValueError if path escapes root/job_id."""
    job_dir = root / job_id
    try:
        path.resolve().relative_to(job_dir.resolve())
    except ValueError:
        raise ValueError(f"Path {path} escapes job directory for {job_id}")


class StorageBackend(Protocol):
    def job_dir(self, job_id: str) -> Path: ...
    def inputs_dir(self, job_id: str) -> Path: ...
    def output_dir(self, job_id: str) -> Path: ...
    def config_path(self, job_id: str) -> Path: ...
    def results_path(self, job_id: str) -> Path: ...
    def captured_log_path(self, job_id: str) -> Path: ...
    def save_upload(self, job_id: str, filename: str, data: BinaryIO) -> Path: ...
    def results_exists(self, job_id: str) -> bool: ...
    def open_results(self, job_id: str) -> BinaryIO: ...
    def create_job_dirs(self, job_id: str) -> None: ...


class LocalStorage:
    def __init__(self, data_dir: str | None = None):
        self._root = Path(data_dir or settings.data_dir).resolve()

    def job_dir(self, job_id: str) -> Path:
        return self._root / job_id

    def inputs_dir(self, job_id: str) -> Path:
        return self.job_dir(job_id) / "inputs"

    def output_dir(self, job_id: str) -> Path:
        return self.job_dir(job_id) / "output"

    def config_path(self, job_id: str) -> Path:
        return self.job_dir(job_id) / "config.json"

    def results_path(self, job_id: str) -> Path:
        return self.job_dir(job_id) / "results.zip"

    def captured_log_path(self, job_id: str) -> Path:
        return self.job_dir(job_id) / "captured.log"

    def create_job_dirs(self, job_id: str) -> None:
        self.inputs_dir(job_id).mkdir(parents=True, exist_ok=True)
        self.output_dir(job_id).mkdir(parents=True, exist_ok=True)

    def save_upload(self, job_id: str, filename: str, data: BinaryIO) -> Path:
        safe_name = _safe_filename(filename)
        dest = self.inputs_dir(job_id) / safe_name
        _assert_within_dir(dest, self._root, job_id)
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "wb") as f:
            while chunk := data.read(1 << 20):  # 1 MB chunks
                f.write(chunk)
        return dest

    def results_exists(self, job_id: str) -> bool:
        return self.results_path(job_id).exists()

    def open_results(self, job_id: str) -> BinaryIO:
        return open(self.results_path(job_id), "rb")


# Module-level singleton
storage = LocalStorage()
