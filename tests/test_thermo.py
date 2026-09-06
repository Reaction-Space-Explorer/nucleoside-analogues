"""Spontaneity classification logic.

These need no eQuilibrator install: they exercise the decision rule that turns
an estimate plus its uncertainty into a call, which is where the scientific
claim actually lives.
"""

from __future__ import annotations

import pytest

from nucleoside_analogues.thermo import DEFAULT_Z, ReactionEnergy, classify, summarise


@pytest.mark.parametrize(
    ("dg", "sigma", "expected"),
    [
        (-100.0, 5.0, "spontaneous"),
        (100.0, 5.0, "non_spontaneous"),
        # the case the original binary filter silently called spontaneous;
        # component contribution routinely carries 10-40 kJ/mol on exotic
        # compounds, and roughly half the deposited generation-3 reactions it
        # called spontaneous sit in this regime
        (-10.0, 20.0, "undetermined"),
        (-0.4, 20.0, "undetermined"),
        (-19.0, 20.0, "undetermined"),
    ],
)
def test_classify_uses_the_whole_interval(dg: float, sigma: float, expected: str) -> None:
    assert classify(dg, sigma) == expected


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


def test_null_estimates_are_not_spontaneous():
    """Component contribution returns 0 +/- 0 where reagents and products share
    a decomposition. A point value of -1e-05 must not pass the 95% test on
    rounding, so such reactions count as unestimable."""
    from make_si_tables import is_null

    assert is_null({"dG_prime_kJ_mol": "-1.018336507740969e-05", "sigma_kJ_mol": "0.0"})
    assert is_null({"dG_prime_kJ_mol": "0.0", "sigma_kJ_mol": "0.0"})
    # a real estimate that happens to be small keeps its uncertainty
    assert not is_null({"dG_prime_kJ_mol": "-0.0001", "sigma_kJ_mol": "1.4"})
    # a real estimate that is large and certain is not null
    assert not is_null({"dG_prime_kJ_mol": "-13.4", "sigma_kJ_mol": "0.0"})
