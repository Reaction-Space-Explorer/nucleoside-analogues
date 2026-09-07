"""Max-min driving force for each traced route.

    uv run --extra thermo python scripts/driving_force.py

A route being spontaneous under standard conditions says nothing about whether
one set of concentrations can drive every step of it at once. The max-min
driving force asks exactly that: choose concentrations within physiological
bounds so as to maximise the smallest driving force along the route, and report
that smallest value. A positive MDF means some single concentration profile
makes every step favourable; the step attaining it is the bottleneck.

This is the linear program of Noor et al., written out directly rather than
taken from a package, since the routes here are short and the formulation is
worth being able to read:

    maximise  B
    over      ln c_i,  B
    subject to  dGr'° _j + RT sum_i S_ij ln c_i <= -B   for every reaction j
                ln c_min <= ln c_i <= ln c_max

Water is held at unit activity and left out of the variables, matching the
convention the transformed energies already use.

Writes ProcessedData/SI/driving_force.csv.
"""

import csv
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import linprog

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import RELS, REPO, SI, deepest

from nucleoside_analogues.rels import pivot_rels

RT = 8.314462618e-3 * 298.15  # kJ/mol
C_MIN, C_MAX = 1e-6, 1e-2  # molar, the usual physiological window
WATER = "O"
ROUTES = REPO / "figures" / "routes"


def mdf(stoich: dict[str, dict[str, float]], dg: dict[str, float]):
    """Max-min driving force, and the reaction attaining it."""
    reactions = sorted(stoich)
    species = sorted({s for r in reactions for s in stoich[r]} - {WATER})
    n = len(species)
    index = {s: k for k, s in enumerate(species)}

    # variables: ln c for each species, then B
    a_ub = np.zeros((len(reactions), n + 1))
    b_ub = np.zeros(len(reactions))
    for j, r in enumerate(reactions):
        for s, coeff in stoich[r].items():
            if s != WATER:
                a_ub[j, index[s]] = RT * coeff
        a_ub[j, n] = 1.0
        b_ub[j] = -dg[r]

    bounds = [(math.log(C_MIN), math.log(C_MAX))] * n + [(None, None)]
    result = linprog(
        c=np.r_[np.zeros(n), -1.0], A_ub=a_ub, b_ub=b_ub, bounds=bounds, method="highs"
    )
    if not result.success:
        return None, None, None
    value = result.x[n]
    slack = b_ub - a_ub @ result.x
    return value, reactions[int(np.argmin(slack))], dict(zip(species, result.x[:n], strict=True))


def route_steps(spec: Path) -> list[str]:
    doc = yaml.safe_load(spec.read_text())
    out: list[str] = []

    def walk(node):
        reaction = node.get("reaction")
        if reaction and reaction["id"] not in out:
            out.append(reaction["id"])
        for child in node.get("from", []):
            walk(child)

    walk(doc["target"])
    return out


def main() -> None:
    rows = []
    cache: dict[str, tuple] = {}
    for spec in sorted(ROUTES.glob("*.yaml")):
        network, target = spec.stem.split("_", 1)
        if network not in cache:
            generation = deepest(network)
            rels = pivot_rels(
                pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t")
            )
            rels["Index"] = rels["Index"].astype(str)
            cache[network] = (
                {
                    i: (tuple(a), tuple(b))
                    for i, a, b in zip(
                        rels["Index"], rels["Reagents"], rels["Products"], strict=True
                    )
                },
                {
                    r["Index"]: float(r["dG_prime_kJ_mol"])
                    for r in csv.DictReader(
                        (SI / "full" / f"{network}_G{generation}_energies_pH7.4.csv").open()
                    )
                    if r["estimable"] == "True" and r["dG_prime_kJ_mol"]
                },
            )
        lut, energies = cache[network]
        steps = route_steps(spec)
        stoich, dg = {}, {}
        for rid in steps:
            reagents, products = lut[rid]
            counts: dict[str, float] = {}
            for s in reagents:
                counts[s] = counts.get(s, 0.0) - 1.0
            for s in products:
                counts[s] = counts.get(s, 0.0) + 1.0
            stoich[rid] = counts
            dg[rid] = energies[rid]
        value, bottleneck, _ = mdf(stoich, dg)
        rows.append(
            {
                "network": network,
                "target": target,
                "reactions": len(steps),
                "sum_dG_kJ_mol": round(sum(dg.values()), 1),
                "mdf_kJ_mol": round(value, 2) if value is not None else "",
                "bottleneck_reaction": bottleneck or "",
                "bottleneck_dG_kJ_mol": round(dg[bottleneck], 1) if bottleneck else "",
            }
        )
        print(
            f"  {network:12s} {target:12s} {len(steps):2d} reactions  "
            f"MDF {value:7.2f} kJ/mol  bottleneck {bottleneck} "
            f"(dGr'o {dg[bottleneck]:+.1f})",
            flush=True,
        )
    out = SI / "driving_force.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
