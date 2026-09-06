"""Length and free energy of the spontaneous route to every matched analogue.

    uv run python scripts/make_pathway_energetics.py

For each network at its deepest generation, over spontaneous reactions only,
every matched nucleoside analogue that is reachable contributes one route: its
chain depth, the sum of dGr'° over the distinct reactions of that route, and
that sum divided by the number of reactions.
Writes ProcessedData/SI/figure_pathway_energetics.csv, one row per network and
length. Energies are kJ/mol, as eQuilibrator returns them.
"""

import csv
import statistics
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import PRODUCTS, RELS, SI, admitted, deepest, energy_file

from nucleoside_analogues.hyperpath import shortest_pathways, trace
from nucleoside_analogues.rels import build_index, pivot_rels, read_products

REPO = Path(__file__).resolve().parent.parent


def main() -> None:
    rows = []
    for network, products_file in PRODUCTS.items():
        generation = deepest(network)
        rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
        rels["Index"] = rels["Index"].astype(str)
        products = read_products(
            REPO / "OriginalData" / "OriginalNetworkData" / "Products" / products_file
        )
        seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
        species = set(products["Smiles"].astype(str))
        matched = set(
            pd.read_csv(
                REPO / "ProcessedData" / "MatchesFiles" / f"{network}Matches.tsv", sep="\t"
            )["NetworkSmiles"].astype(str)
        )
        energies = {
            r["Index"]: float(r["dG_prime_kJ_mol"])
            for r in csv.DictReader(energy_file(network, generation).open())
            if r["estimable"] == "True" and r["dG_prime_kJ_mol"]
        }
        index = build_index(
            rels[rels["Index"].isin(admitted(network, rels, "estimable_only", generation))]
        )
        result = shortest_pathways(index, seeds)

        by_length: dict[int, list[float]] = {}
        for smiles in matched & species:
            cost = result.cost.get(smiles)
            if cost is None or cost == 0:
                continue
            reactions = trace(result, index, smiles)
            if any(r not in energies for r in reactions):
                continue  # every step of a spontaneous route has a value; skip if not
            total = sum(energies[r] for r in reactions)
            by_length.setdefault(int(cost), []).append((total, total / len(reactions)))

        for length in sorted(by_length):
            totals = [t for t, _ in by_length[length]]
            per_step = [p for _, p in by_length[length]]
            rows.append(
                {
                    "network": network,
                    "generation": generation,
                    "length": length,
                    "pathways": len(totals),
                    "mean_dg_kJ_mol": round(statistics.fmean(totals), 2),
                    "sd_dg_kJ_mol": round(statistics.stdev(totals), 2) if len(totals) > 1 else 0.0,
                    "median_dg_kJ_mol": round(statistics.median(totals), 2),
                    "mean_dg_per_step_kJ_mol": round(statistics.fmean(per_step), 2),
                    "sd_dg_per_step_kJ_mol": (
                        round(statistics.stdev(per_step), 2) if len(per_step) > 1 else 0.0
                    ),
                }
            )
        total = sum(len(v) for v in by_length.values())
        spread = f"{min(by_length)}-{max(by_length)}" if by_length else "none"
        print(f"  {network:12s} G{generation}  routes {total:6,d}  lengths {spread}", flush=True)

    out = SI / "figure_pathway_energetics.csv"
    with out.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
