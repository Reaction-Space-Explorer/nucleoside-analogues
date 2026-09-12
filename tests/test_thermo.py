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


@pytest.mark.slow
def test_a_reaction_and_its_reverse_have_opposite_free_energies(network: str) -> None:
    """Thermodynamic consistency over the whole deposited set.

    Where the network contains both A -> B and B -> A, component contribution
    must give exactly opposite free energies. It does, to floating point: over
    40,007 such pairs in Formose the largest deviation is 2e-13 kJ/mol.
    """
    import csv

    import pandas as pd
    from helpers import ORIGINAL
    from make_si_tables import FULL, deepest, is_null

    from nucleoside_analogues.rels import pivot_rels

    generation = deepest(network)
    energies = FULL / f"{network}_G{generation}_energies_pH7.4.csv"
    if not energies.exists():
        pytest.skip("full-depth energies not present")
    rels = pivot_rels(
        pd.read_csv(ORIGINAL / "Rels" / network / f"{network}Rels_{generation}.tsv", sep="\t")
    )
    rels["Index"] = rels["Index"].astype(str)
    table = {r["Index"]: r for r in csv.DictReader(energies.open())}
    forward = {
        (tuple(sorted(a)), tuple(sorted(b))): i
        for i, a, b in zip(rels["Index"], rels["Reagents"], rels["Products"], strict=True)
    }
    checked = 0
    for (reagents, products), i in forward.items():
        j = forward.get((products, reagents))
        if j is None or i >= j:
            continue
        a, b = table.get(i), table.get(j)
        if not a or not b or a["estimable"] != "True" or b["estimable"] != "True":
            continue
        if is_null(a) or is_null(b):
            continue
        checked += 1
        assert abs(float(a["dG_prime_kJ_mol"]) + float(b["dG_prime_kJ_mol"])) < 1e-6, (i, j)
    assert checked, f"no reversible pairs found in {network}"


def test_unusable_reason_separates_the_three_causes():
    """Three different things are reported as unestimable and are not interchangeable.

    A sentinel-variance reaction still carries a dG, and a null estimate carries
    the most confident uncertainty in the file, so neither is caught by testing
    the value or the error alone.
    """
    from nucleoside_analogues.thermo import unusable_reason

    assert unusable_reason(373.0, 1e5) == "unbounded_variance"
    assert unusable_reason(-373.0, 1e5) == "unbounded_variance"
    assert unusable_reason(0.0, 0.0) == "null_estimate"
    assert unusable_reason(-1e-5, 0.0) == "null_estimate"
    assert unusable_reason(None, None, "compound_missing") == "no_estimate"
    # a real estimate, however small, is usable
    assert unusable_reason(-13.4, 0.0) is None
    assert unusable_reason(-40.0, 1.9) is None


@pytest.mark.slow
def test_the_three_causes_account_for_every_unestimable_formose_reaction():
    """Pins the counts the notebooks README quotes, from the deposited energies."""
    import csv
    from collections import Counter

    from helpers import REPO

    from nucleoside_analogues.thermo import unusable_reason

    path = REPO / "ProcessedData" / "SI" / "full" / "Formose_G6_energies_pH7.4.csv"
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    causes = Counter()
    for row in rows:
        dg = float(row["dG_prime_kJ_mol"]) if row["dG_prime_kJ_mol"] else None
        sigma = float(row["sigma_kJ_mol"]) if row["sigma_kJ_mol"] else None
        causes[unusable_reason(dg, sigma, row["status"])] += 1

    assert len(rows) == 306244
    assert causes["unbounded_variance"] == 55721
    assert causes["null_estimate"] == 5466
    assert causes["no_estimate"] == 553
    assert sum(v for k, v in causes.items() if k) == 61740
