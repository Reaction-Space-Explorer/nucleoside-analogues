"""Spontaneity classification logic.

These need no eQuilibrator install: they exercise the decision rule that turns
an estimate plus its uncertainty into a call, which is where the scientific
claim actually lives.
"""

from __future__ import annotations

import pytest

from nucleoside_analogues.thermo import DEFAULT_Z, ReactionEnergy, classify, summarise


def test_confidently_negative_is_spontaneous() -> None:
    assert classify(-100.0, 5.0) == "spontaneous"


def test_confidently_positive_is_not_spontaneous() -> None:
    assert classify(100.0, 5.0) == "non_spontaneous"


def test_interval_spanning_zero_is_undetermined() -> None:
    """The case the original binary filter silently called spontaneous."""
    assert classify(-10.0, 20.0) == "undetermined"


@pytest.mark.parametrize("dg", [-0.4, -5.0, -19.0])
def test_small_negative_values_are_undetermined_at_realistic_error(dg: float) -> None:
    """Component contribution routinely carries ~10-40 kJ/mol on exotic compounds.

    Roughly half the deposited generation-3 reactions classified spontaneous sit
    in this regime.
    """
    assert classify(dg, 20.0) == "undetermined"


def test_without_uncertainty_it_degrades_to_the_sign() -> None:
    """Reproduces the original filter, for comparison against deposited results."""
    assert classify(-0.4, None) == "spontaneous"
    assert classify(0.4, None) == "non_spontaneous"


def test_boundary_is_exactly_z_sigma() -> None:
    sigma = 10.0
    just_inside = -(DEFAULT_Z * sigma) - 1e-6
    just_outside = -(DEFAULT_Z * sigma) + 1e-6
    assert classify(just_inside, sigma) == "spontaneous"
    assert classify(just_outside, sigma) == "undetermined"


def test_failures_stay_countable_rather_than_vanishing() -> None:
    energies = [
        ReactionEnergy("r1", -100.0, 5.0, "ok"),
        ReactionEnergy("r2", -10.0, 20.0, "ok"),
        ReactionEnergy("r3", 50.0, 5.0, "ok"),
        ReactionEnergy("r4", None, None, "compound_missing", "C#N"),
        ReactionEnergy("r5", None, None, "estimation_failed", "solver"),
    ]
    counts = summarise(energies)
    assert counts["total"] == 5
    assert counts["estimated"] == 3
    assert counts["spontaneous"] == 1
    assert counts["undetermined"] == 1
    assert counts["non_spontaneous"] == 1
    assert counts["compound_missing"] == 1
    assert counts["estimation_failed"] == 1
    accounted = (
        counts["estimated"]
        + counts["compound_missing"]
        + counts["decomposition_failed"]
        + counts["estimation_failed"]
    )
    assert accounted == counts["total"], "every reaction must be accounted for"


def test_reaction_energy_reports_its_own_call() -> None:
    assert ReactionEnergy("r", -100.0, 5.0, "ok").spontaneity() == "spontaneous"
    assert ReactionEnergy("r", None, None, "compound_missing").spontaneity() is None
