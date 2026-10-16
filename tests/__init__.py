"""Test suite for the oktoberfest package."""

from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / "data"
CONFIGS_DIR = TESTS_DIR / "configs"

__all__ = ["TESTS_DIR", "DATA_DIR", "CONFIGS_DIR"]
