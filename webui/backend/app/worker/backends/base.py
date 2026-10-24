from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable


@dataclass
class ExecutionResult:
    returncode: int
    stdout: str
    stderr: str


@runtime_checkable
class ExecutionBackend(Protocol):
    def run(self, job_dir: Path, config_path: Path) -> ExecutionResult: ...
    def cancel(self, handle: str) -> None: ...
