"""Pytest configuration shared by all gates.

Ensures the repository root is importable without an editable install, and
provides the fiducial config as a fixture.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from config import Config, load_config  # noqa: E402


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return _REPO_ROOT


@pytest.fixture(scope="session")
def cfg() -> Config:
    """Fiducial run config (SPEC §8). Frozen, so safe to share across tests."""
    return load_config(_REPO_ROOT / "config" / "fiducial.yaml")
