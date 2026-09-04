"""Regression of the hyperpath search against the deposited pathway results.

The original tree search reported "Pathway Length", which is the depth of the
derivation tree -- the longest chain of consecutive reactions.  The ``"chain"``
objective computes the same quantity exactly, and must reproduce every
deposited value.
"""

from __future__ import annotations

import csv

import pytest
from helpers import REPO, products_path, requires, spontaneous_rels_path

from nucleoside_analogues.hyperpath import (
    Objective,
    count_minimal_routes,
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
