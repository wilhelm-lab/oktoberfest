from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from app.worker.backends.base import ExecutionBackend, ExecutionResult


class LocalExecutionBackend:
    """Runs oktoberfest as a child subprocess, capturing stdout+stderr to a log file."""

    def run(self, job_dir: Path, config_path: Path) -> ExecutionResult:
        log_path = job_dir / "captured.log"
        cmd = [sys.executable, "-m", "oktoberfest", "--config_path", str(config_path)]
        with open(log_path, "w") as log_file:
            proc = subprocess.run(
                cmd,
                cwd=str(job_dir),
                stdout=log_file,
                stderr=subprocess.STDOUT,
                text=True,
            )
        stdout = log_path.read_text(errors="replace") if log_path.exists() else ""
        return ExecutionResult(returncode=proc.returncode, stdout=stdout, stderr="")

    def cancel(self, handle: str) -> None:
        # Subprocess already terminated; no-op for completed/queued tasks
        pass


# Registry — add new backends here; no other file needs to change
_BACKENDS: dict[str, type] = {
    "local": LocalExecutionBackend,
}


def get_execution_backend(backend_name: str = "local") -> ExecutionBackend:
    cls = _BACKENDS.get(backend_name)
    if cls is None:
        raise ValueError(f"Unknown execution backend: {backend_name!r}. Available: {list(_BACKENDS)}")
    return cls()
