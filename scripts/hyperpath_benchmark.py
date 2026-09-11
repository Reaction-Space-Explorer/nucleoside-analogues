"""Scale and running time of Algorithm 1 on the full network of each CRNR.

Times the search itself, not the reshape that precedes it. The largest network
is Formose at generation six; the reshape cost quoted in ``rels`` is for the
generation-four FormoseAmm network and is a different operation.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import pandas as pd

from nucleoside_analogues.hyperpath import shortest_pathways
from nucleoside_analogues.rels import build_index, pivot_rels, read_products

REPO = Path(__file__).resolve().parent.parent
RELS = REPO / "OriginalData" / "OriginalNetworkData" / "Rels"
OUT = REPO / "ProcessedData" / "SI" / "hyperpath_benchmark.csv"
PRODUCTS = {
    "Formose": "formose_output.tsv",
    "FormoseAmm": "formose_amm_output.tsv",
    "Glucose": "glucose_degradation_output.tsv",
    "GlucoseAmm": "glucose_amm_output.tsv",
    "PyruvicAcid": "pyruvic_output.tsv",
}
REPEATS = 5


def deepest(network: str) -> int:
    return max(int(f.stem.split("_")[-1]) for f in (RELS / network).glob("*Rels_*.tsv"))


def main() -> None:
    rows = []
    for network, products_file in PRODUCTS.items():
        generation = deepest(network)
        rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
        rels["Index"] = rels["Index"].astype(str)
        index = build_index(rels)
        products = read_products(
            REPO / "OriginalData" / "OriginalNetworkData" / "Products" / products_file
        )
        seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
        species = {s for r in index.reagents.values() for s in r}
        species |= {s for p in index.products.values() for s in p}

        row = {
            "network": network,
            "generation": generation,
            "reactions": len(index.reagents),
            "species": len(species),
            "seeds": len(seeds),
        }
        for objective in ("chain", "reactions"):
            best = min(_time(index, seeds, objective) for _ in range(REPEATS))
            row[f"seconds_{objective}"] = round(best, 3)
            row[f"reached_{objective}"] = len(
                shortest_pathways(index, seeds, objective=objective).cost
            )
        rows.append(row)
        print(
            f"  {network:<12} G{generation}  {row['reactions']:>7,} reactions"
            f"  {row['species']:>7,} species"
            f"  chain {row['seconds_chain']:>6.3f}s"
            f"  reactions {row['seconds_reactions']:>6.3f}s"
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    largest = max(rows, key=lambda r: r["reactions"])
    slowest = max(rows, key=lambda r: max(r["seconds_chain"], r["seconds_reactions"]))
    print(f"\n  largest: {largest['network']} G{largest['generation']}")
    print(f"  slowest: {slowest['network']} G{slowest['generation']}")
    print(f"  wrote {OUT.relative_to(REPO)}")


def _time(index, seeds, objective: str) -> float:
    start = time.perf_counter()
    shortest_pathways(index, seeds, objective=objective)
    return time.perf_counter() - start


if __name__ == "__main__":
    main()
