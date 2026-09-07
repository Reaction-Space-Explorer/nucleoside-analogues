"""The max-min driving force LP, checked against cases with known answers."""

import math

import pytest
from driving_force import C_MAX, C_MIN, RT, mdf


def test_single_reaction_matches_the_closed_form() -> None:
    """One reaction, so the optimum puts every reactant at the top of the range
    and the product at the bottom, and the answer can be written down."""
    value, bottleneck, conc = mdf({"r": {"A": -2.0, "B": 1.0}}, {"r": -16.4})
    assert value is not None and conc is not None
    expected = -(-16.4 + RT * math.log(C_MIN / C_MAX**2))
    assert bottleneck == "r"
    assert value == pytest.approx(expected, abs=1e-9)
    assert math.exp(conc["A"]) == pytest.approx(C_MAX)
    assert math.exp(conc["B"]) == pytest.approx(C_MIN)


def test_the_bottleneck_is_the_least_favourable_step() -> None:
    steps = {"a": {"X": -1.0, "Y": 1.0}, "b": {"Y": -1.0, "Z": 1.0}}
    _, bottleneck, _ = mdf(steps, {"a": -40.0, "b": 5.0})
    assert bottleneck == "b"


def test_a_route_that_cannot_be_driven_gives_a_negative_value() -> None:
    steps = {"a": {"X": -1.0, "Y": 1.0}, "b": {"Y": -1.0, "Z": 1.0}}
    value, _, _ = mdf(steps, {"a": -1.0, "b": 60.0})
    assert value is not None and value < 0


def test_water_is_held_at_unit_activity() -> None:
    """Water is not a variable, so adding it to a reaction cannot change the
    answer; the transformed energies already assume it."""
    dry, _, _ = mdf({"r": {"A": -1.0, "B": 1.0}}, {"r": -20.0})
    wet, _, _ = mdf({"r": {"A": -1.0, "B": 1.0, "O": -1.0}}, {"r": -20.0})
    assert dry is not None and wet is not None
    assert dry == pytest.approx(wet)
