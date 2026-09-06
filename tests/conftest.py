"""Shared fixtures. Tests run against the deposited data in this repository."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from helpers import NETWORKS

# Tests cover the analysis entry points in scripts/, which is not a package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))


@pytest.fixture(scope="session", params=sorted(NETWORKS))
def network(request: pytest.FixtureRequest) -> str:
    return request.param
