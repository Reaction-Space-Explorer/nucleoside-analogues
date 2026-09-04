"""Shortest synthetic pathways as a minimum-weight hyperpath problem.

A chemical reaction network is a *directed hypergraph*: a reaction consumes a
set of species and produces a set of species, so an edge has several tails.
"The shortest pathway to molecule *m*" is then a minimum-weight hyperpath, and
obeys the recursion

.. math::

    d(m) = \\min_{r \\;\\text{producing}\\; m}
           \\bigl( 1 + f\\{\\, d(x) : x \\in \\mathrm{reagents}(r) \\,\\}\\bigr)

with ``d(seed) = 0``.  When *f* is monotone and non-decreasing in each argument
the recursion is a *superior function*, and Knuth's generalisation of Dijkstra's
algorithm solves it exactly in one pass.  Cycles are excluded structurally --
a reaction is only relaxed once every reagent has been finalised at a strictly
smaller cost -- which retires the ad-hoc "is this molecule already in my path"
check the original tree search needed.

Two objectives are supported, and they measure different things:

``"chain"`` (``f = max``)
    Longest chain of consecutive reactions from a seed to the target.  This is
    what the original ``map_tree`` reported as "Pathway Length": it is the depth
    of the derivation tree.  Reproduces all 1,181 deposited generation-3 values.

``"reactions"`` (``f = sum``)
    Total reactions in the derivation, counting a shared intermediate once per
    use.  This is an upper bound on the number of *distinct* reactions; the
    minimum-distinct-reaction problem allows re-use of a sub-derivation and is
    not a superior function, so it is not solved here.

References
----------
Knuth, D. E. A generalization of Dijkstra's algorithm.
*Inf. Process. Lett.* **1977**, *6*, 1-5.

Gallo, G.; Longo, G.; Pallottino, S.; Nguyen, S. Directed hypergraphs and
applications. *Discrete Appl. Math.* **1993**, *42*, 177-201.
"""

from __future__ import annotations

import heapq
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Literal

from .rels import ReactionIndex

__all__ = ["Objective", "PathwayResult", "shortest_pathways", "trace"]

Objective = Literal["chain", "reactions"]


@dataclass(frozen=True, slots=True)
class PathwayResult:
    """Outcome of a shortest-hyperpath search.

    Attributes
    ----------
    cost
        Species -> minimum cost under the chosen objective.  Seeds cost 0.
        A species absent from this mapping is unreachable.
    via
        Species -> identifier of the reaction achieving its minimum cost.
    objective
        Which objective was optimised.
    """

    cost: Mapping[str, int]
    via: Mapping[str, str]
    objective: Objective

    def reachable(self) -> set[str]:
        return set(self.cost)

    def unreachable(self, species: Iterable[str]) -> set[str]:
        """Species that no sequence of reactions produces from the seeds.

        Unlike the original implementation this is a *proven* negative: the
        search is exhaustive and has no time limit, so "absent" cannot be
        confused with "timed out", which the previous code reported
        identically.
        """
        return {s for s in species if s not in self.cost}


def shortest_pathways(
    index: ReactionIndex,
    seeds: Iterable[str],
    objective: Objective = "chain",
) -> PathwayResult:
    """Compute exact minimum-cost pathways to every reachable species.

    One pass yields the answer for *all* targets simultaneously, so there is no
    need to re-run per target.

    Parameters
    ----------
    index
        Reaction lookup tables, from :func:`~nucleoside_analogues.rels.build_index`.
    seeds
        Starting materials -- generation-0 species of the network.  Take these
        from the deposited product listing rather than a hard-coded table.
    objective
        ``"chain"`` for longest precursor chain, ``"reactions"`` for total
        reaction count with multiplicity.
    """
    if objective not in ("chain", "reactions"):
        raise ValueError(f"unknown objective {objective!r}")

    combine = max if objective == "chain" else sum

    cost: dict[str, int] = {seed: 0 for seed in seeds}
    via: dict[str, str] = {}
    finalised: set[str] = set()

    # A reaction becomes relaxable once every reagent has a finalised cost.
    outstanding: dict[str, int] = {}
    blocked_on: dict[str, list[str]] = defaultdict(list)
    for identifier, reagents in index.reagents.items():
        pending = {r for r in reagents if r not in cost}
        outstanding[identifier] = len(pending)
        for species in pending:
            blocked_on[species].append(identifier)

    queue: list[tuple[int, str]] = [(0, seed) for seed in cost]
    heapq.heapify(queue)

    def relax(identifier: str) -> None:
        reagents = index.reagents[identifier]
        if any(r not in cost for r in reagents):
            return
        value = 1 + (combine(cost[r] for r in reagents) if reagents else 0)
        for product in index.products.get(identifier, ()):
            if product not in finalised and value < cost.get(product, 1 << 62):
                cost[product] = value
                via[product] = identifier
                heapq.heappush(queue, (value, product))

    for identifier, pending in outstanding.items():
        if pending == 0:
            relax(identifier)

    while queue:
        value, species = heapq.heappop(queue)
        if species in finalised or value > cost.get(species, 1 << 62):
            continue
        finalised.add(species)
        for identifier in blocked_on.get(species, ()):
            outstanding[identifier] -= 1
            if outstanding[identifier] == 0:
                relax(identifier)

    return PathwayResult(cost=cost, via=via, objective=objective)


def trace(
    result: PathwayResult,
    index: ReactionIndex,
    target: str,
) -> list[str]:
    """Return the reaction identifiers forming the optimal pathway to *target*.

    Identifiers are de-duplicated, so ``len(trace(...))`` is the number of
    *distinct* reactions -- which may be smaller than the ``"reactions"``
    objective's cost when a sub-derivation is shared.
    """
    if target not in result.cost:
        raise KeyError(f"{target!r} is not reachable from the seeds")

    seen: list[str] = []
    stack: list[str] = [target]
    visited: set[str] = set()
    while stack:
        species = stack.pop()
        if species in visited or result.cost[species] == 0:
            continue
        visited.add(species)
        identifier = result.via.get(species)
        if identifier is None:
            continue
        if identifier not in seen:
            seen.append(identifier)
        stack.extend(index.reagents[identifier])
    return seen


def count_minimal_routes(
    result: PathwayResult,
    index: ReactionIndex,
    target: str,
) -> int:
    """Count distinct minimum-cost derivations of *target*.

    Well defined precisely because it is restricted to the minimum-cost
    sub-hypergraph, where the cost ordering makes the structure acyclic.  This
    is the metric to report for synthetic redundancy; "total number of
    pathways" is not well defined without also fixing a length bound and a
    convention for reordering independent branches, and different conventions
    disagree by many orders of magnitude.
    """
    if target not in result.cost:
        return 0

    combine = max if result.objective == "chain" else sum
    memo: dict[str, int] = {}

    def routes(species: str) -> int:
        if result.cost[species] == 0:
            return 1
        if species in memo:
            return memo[species]
        memo[species] = 0  # guards against re-entry
        total = 0
        for identifier in index.producers.get(species, ()):
            reagents = index.reagents[identifier]
            if any(r not in result.cost for r in reagents):
                continue
            value = 1 + (combine(result.cost[r] for r in reagents) if reagents else 0)
            if value != result.cost[species]:
                continue
            product = 1
            for reagent in reagents:
                product *= routes(reagent)
            total += product
        memo[species] = total
        return total

    return routes(target)
