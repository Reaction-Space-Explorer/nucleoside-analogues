"""SI File 1: the matched analogues each CRNR reaches by spontaneous reactions.

    uv run python scripts/spontaneous_smiles.py

The deposited version of this listing was built in 2023 from the generation-3
networks under a hard dGr'o < 0 test on the point estimate. The manuscript now
reports each network at the deepest generation it was generated to, and calls a
reaction spontaneous only when the whole 95% interval lies below zero. The two
differ by an order of magnitude, so the listing is regenerated here from the
same functions that produce the reachability funnel.

Both admission bases are written, because the paper reports both: `estimable_only`
uses only reactions component contribution can evaluate, and `with_unestimable`
additionally allows reactions carrying no usable estimate to pass.

Writes ProcessedData/SpontaneousSMILES/.
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import PRODUCTS, RELS, REPO, admitted, deepest, species_generations

from nucleoside_analogues.hyperpath import shortest_pathways
from nucleoside_analogues.rels import build_index, pivot_rels, read_products

OUT = REPO / "ProcessedData" / "SpontaneousSMILES"
BASES = ("estimable_only", "with_unestimable")


def write(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reached: dict[str, dict[str, set[str]]] = {b: {} for b in BASES}
    matched_by: dict[str, set[str]] = {}
    gens: dict[str, dict[str, int]] = {}

    for network, products_file in PRODUCTS.items():
        generation = deepest(network)
        rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
        rels["Index"] = rels["Index"].astype(str)
        products = read_products(
            REPO / "OriginalData" / "OriginalNetworkData" / "Products" / products_file
        )
        seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
        gens[network] = species_generations(network)
        matched = set(
            pd.read_csv(
                REPO / "ProcessedData" / "MatchesFiles" / f"{network}Matches.tsv", sep="\t"
            )["NetworkSmiles"].astype(str)
        )
        matched_by[network] = matched & set(gens[network])

        for basis in BASES:
            index = build_index(
                rels[rels["Index"].isin(admitted(network, rels, basis, generation))]
            )
            reached[basis][network] = matched_by[network] & set(
                shortest_pathways(index, seeds).cost
            )
        print(
            f"  {network:12s} G{generation}  matched {len(matched_by[network]):6,d}"
            f"  spontaneous {len(reached['estimable_only'][network]):6,d}"
            f"  (with unestimable {len(reached['with_unestimable'][network]):6,d})",
            flush=True,
        )

    # per-network and combined listings, on the basis the paper leads with
    combined = []
    for network in PRODUCTS:
        rows = sorted(
            (network, f"G{gens[network][s]}", s) for s in reached["estimable_only"][network]
        )
        write(OUT / f"{network}.tsv", ["Network", "Generation", "Smiles"], [list(r) for r in rows])
        combined += rows
    write(
        OUT / "AllSpontaneousSmiles.tsv",
        ["Network", "Generation", "Smiles"],
        [list(r) for r in combined],
    )

    permissive = sorted(
        (network, f"G{gens[network][s]}", s)
        for network in PRODUCTS
        for s in reached["with_unestimable"][network]
    )
    write(
        OUT / "AllSpontaneousSmiles_with_unestimable.tsv",
        ["Network", "Generation", "Smiles"],
        [list(r) for r in permissive],
    )

    # spontaneous against non-spontaneous, per generation
    header = ["Generation"]
    for network in PRODUCTS:
        header += [f"{network}_S", f"{network}_NS"]
    depth = max(max(g.values()) for g in gens.values())
    counts = []
    for generation in range(1, depth + 1):
        row = [f"G{generation}"]
        for network in PRODUCTS:
            at = {s for s in matched_by[network] if gens[network][s] == generation}
            spont = at & reached["estimable_only"][network]
            row += [len(spont), len(at - spont)]
        counts.append(row)
    write(OUT / "SpontVsNonSpont.tsv", header, counts)

    # species reached spontaneously in exactly one network
    seen = Counter(s for network in PRODUCTS for s in reached["estimable_only"][network])
    unique = sorted(
        (s, f"G{gens[network][s]}", network)
        for network in PRODUCTS
        for s in reached["estimable_only"][network]
        if seen[s] == 1
    )
    write(OUT / "UniqueSmiles.tsv", ["Smiles", "Generation", "Network"], [list(r) for r in unique])

    total = sum(len(v) for v in reached["estimable_only"].values())
    print(f"\n  {total:,} spontaneous matches in all, {len(unique):,} reached by one network only")
    print(f"  wrote {OUT.relative_to(REPO)}/")


if __name__ == "__main__":
    main()
