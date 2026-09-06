"""Numbers behind the reachability funnel and the sink/hill classification.

    uv run python scripts/make_figure_data.py

Writes CSVs to ProcessedData/SI/; the figures are drawn separately so that
plotting and computation can run on different machines.
"""

import ast
import collections
import csv

import pandas as pd

import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from make_si_tables import PRODUCTS, RELS, SI, admitted, deepest, energy_file

from nucleoside_analogues.hyperpath import shortest_pathways
from nucleoside_analogues.rels import build_index, pivot_rels, read_products


def literal(v):
    return tuple(ast.literal_eval(v)) if isinstance(v, str) else tuple(v)


def main() -> None:
    funnel, thermo = [], []
    for network, products_file in PRODUCTS.items():
        generation = deepest(network)
        rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
        rels["Index"] = rels["Index"].astype(str)
        products = read_products(
            SI.parent.parent / "OriginalData" / "OriginalNetworkData" / "Products" / products_file
        )
        seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
        species = set(products["Smiles"].astype(str))
        matched = set(
            pd.read_csv(
                SI.parent / "MatchesFiles" / f"{network}Matches.tsv", sep="\t"
            )["NetworkSmiles"].astype(str)
        )

        row = {"network": network, "generation": generation,
               "products": len(species), "matched": len(matched & species)}
        reach = {}
        for basis in ("with_unestimable", "estimable_only"):
            index = build_index(rels[rels["Index"].isin(admitted(network, rels, basis, generation))])
            result = shortest_pathways(index, seeds)
            reach[basis] = result
            row[f"reachable_{basis}"] = len(matched & set(result.cost))
        funnel.append(row)
        print(f"  {network:12s} G{generation} products {row['products']:7,d} | matched {row['matched']:6,d} "
              f"| spontaneous {row['reachable_with_unestimable']:6,d} | estimable-only "
              f"{row['reachable_estimable_only']:6,d}", flush=True)

        # sinks and hills, after Wolos et al. every incoming exothermic and every
        # outgoing endothermic is a sink; the reverse is a hill
        energies = {r["Index"]: r for r in csv.DictReader(energy_file(network, generation).open())}
        incoming, outgoing = collections.defaultdict(list), collections.defaultdict(list)
        for identifier, reagents, products_ in zip(
            rels["Index"], rels["Reagents"], rels["Products"], strict=True
        ):
            e = energies.get(str(identifier))
            if e is None or e["estimable"] != "True" or not e["dG_prime_kJ_mol"]:
                continue
            dg = float(e["dG_prime_kJ_mol"])
            for s in literal(products_):
                incoming[s].append(dg)
            for s in literal(reagents):
                outgoing[s].append(dg)
        sink = hill = neither = 0
        for s in matched & species:
            i, o = incoming.get(s), outgoing.get(s)
            if not i or not o:
                continue
            if all(x < 0 for x in i) and all(x > 0 for x in o):
                sink += 1
            elif all(x > 0 for x in i) and all(x < 0 for x in o):
                hill += 1
            else:
                neither += 1
        thermo.append({"network": network, "generation": generation, "sinks": sink,
                       "hills": hill, "neither": neither, "classified": sink + hill + neither})
        print(f"               sinks {sink:5,d} | hills {hill:5,d} | neither {neither:6,d}", flush=True)

    for rows, name in ((funnel, "figure_funnel.csv"), (thermo, "figure_sinks_hills.csv")):
        with (SI / name).open("w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        print("wrote", name)


if __name__ == "__main__":
    main()
