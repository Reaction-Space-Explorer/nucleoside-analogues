"""Shared fixtures. Tests run against the deposited data in this repository."""

from __future__ import annotations

import pytest
from helpers import NETWORKS


@pytest.fixture(scope="session", params=sorted(NETWORKS))
def network(request: pytest.FixtureRequest) -> str:
    return request.param
