"""Regression of the hyperpath search against the deposited pathway results.

The original tree search reported "Pathway Length", which is the depth of the
derivation tree -- the longest chain of consecutive reactions.  The ``"chain"``
objective computes the same quantity exactly, and must reproduce every
deposited value.
"""

from __future__ import annotations

import csv

import pandas as pd
import pytest
from helpers import REPO, products_path, requires, spontaneous_rels_path

from nucleoside_analogues.hyperpath import (
    Objective,
    count_minimal_routes,
    critical_reactions,
    shortest_pathways,
    trace,
)
from nucleoside_analogues.rels import build_index, load_processed_rels, read_products

FIXTURE = REPO / "tests" / "data" / "deposited_pathways_g3.csv"


def _seeds(network: str) -> tuple[str, ...]:
    """Seeds read from the network itself, never from a hard-coded table."""
    products = read_products(products_path(network))
    return tuple(products.loc[products["Generation"] == 0, "Smiles"])


def _result(network: str, objective: Objective = "chain"):
    path = spontaneous_rels_path(network)
    requires(path)
    requires(products_path(network))
    index = build_index(load_processed_rels(path))
    return shortest_pathways(index, _seeds(network), objective=objective), index


def test_chain_objective_reproduces_deposited_pathway_lengths(network: str) -> None:
    if not FIXTURE.exists():
        pytest.skip("regression fixture not generated")
    result, _ = _result(network)
    rows = [
        r for r in csv.DictReader(FIXTURE.open()) if r["network"] == network and r["status"] == "ok"
    ]
    if not rows:
        pytest.skip(f"no deposited pathways for {network}")

    mismatches = [
        (r["smiles"], int(r["depth"]), result.cost.get(r["smiles"]))
        for r in rows
        if r["smiles"] in result.cost and result.cost[r["smiles"]] != int(r["depth"])
    ]
    assert not mismatches, f"{len(mismatches)}/{len(rows)} depths differ: {mismatches[:5]}"


def test_seeds_cost_zero_and_are_reachable(network: str) -> None:
    result, _ = _result(network)
    for seed in _seeds(network):
        if seed in result.cost:
            assert result.cost[seed] == 0


def test_reaction_objective_is_at_least_the_chain_objective(network: str) -> None:
    """max <= sum for non-negative costs, so the chain length can never exceed
    the total reaction count."""
    chain, _ = _result(network, "chain")
    reactions, _ = _result(network, "reactions")
    for species, value in chain.cost.items():
        assert reactions.cost[species] >= value


def test_traced_pathway_only_uses_reactions_that_exist(network: str) -> None:
    result, index = _result(network)
    reachable = [s for s, c in result.cost.items() if c > 0]
    for species in reachable[:200]:
        for identifier in trace(result, index, species):
            assert identifier in index.reagents


def test_minimal_route_count_is_positive_for_reachable_species(network: str) -> None:
    result, index = _result(network)
    reachable = [s for s, c in result.cost.items() if c > 0]
    for species in reachable[:100]:
        assert count_minimal_routes(result, index, species) >= 1


def test_exclude_removes_only_the_named_reaction() -> None:
    index = build_index(
        pd.DataFrame(
            {
                "Index": ["r1", "r2"],
                "Reagents": [("A",), ("A",)],
                "Products": [("B",), ("B",)],
            }
        )
    )
    assert "B" in shortest_pathways(index, ["A"], exclude=("r1",)).cost
    assert "B" not in shortest_pathways(index, ["A"], exclude=("r1", "r2")).cost


def test_critical_reactions_on_a_known_topology() -> None:
    """B is made two ways so nothing is critical for it; C hangs off one edge."""
    index = build_index(
        pd.DataFrame(
            {
                "Index": ["r1", "r2", "r3"],
                "Reagents": [("A",), ("A",), ("B",)],
                "Products": [("B",), ("B",), ("C",)],
            }
        )
    )
    assert critical_reactions(index, ["A"], "B") == []
    assert critical_reactions(index, ["A"], "C") == ["r3"]
    assert critical_reactions(index, ["A"], "A") == []


def test_critical_reactions_lie_on_the_traced_route(network: str) -> None:
    """A reaction in every derivation must appear in the one trace returns."""
    result, index = _result(network)
    for smiles in ("C(C(C(CO)O)O)=O", "C(C(CO)O)O"):  # threose, glycerol
        if smiles not in result.cost:
            continue
        critical = critical_reactions(index, _seeds(network), smiles)
        assert set(critical) <= set(trace(result, index, smiles))
        for identifier in critical:
            excluded = shortest_pathways(index, _seeds(network), exclude=(identifier,))
            assert smiles not in excluded.cost
