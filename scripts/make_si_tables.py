"""Regenerate SI Tables 1 and 2 from the recomputed generation-3 energies.

uv run python scripts/make_si_tables.py
"""

import csv
from pathlib import Path

import pandas as pd

from nucleoside_analogues.hyperpath import (
    count_minimal_routes,
    critical_reactions,
    shortest_pathways,
)
from nucleoside_analogues.rels import build_index, read_products

REPO = Path(__file__).resolve().parent.parent
SI = REPO / "ProcessedData" / "SI"
PRODUCTS = {
    "Formose": "formose_output.tsv",
    "FormoseAmm": "formose_amm_output.tsv",
    "Glucose": "glucose_degradation_output.tsv",
    "GlucoseAmm": "glucose_amm_output.tsv",
    "PyruvicAcid": "pyruvic_output.tsv",
}
TARGETS = {
    "Deoxyribose": "C(CC(C(CO)O)O)=O",
    "Ribose": "C(C(C(C(CO)O)O)O)=O",
    "Threose": "C(C(C(CO)O)O)=O",
    "Glycerol": "C(C(CO)O)O",
}


#: A reaction counts as spontaneous only where the whole 95% interval is
#: below zero, matching the Methods. The tables are identical on the point
#: estimate, so nothing here depends on the choice.
SPONTANEITY = "spontaneous_95"


def admitted(network: str, rels: pd.DataFrame, basis: str) -> set[str]:
    """Reaction ids allowed under a basis. Reactions with no energy row could
    not be built at all, and count as unestimable rather than vanishing."""
    energies = {
        row["Index"]: row
        for row in csv.DictReader((SI / f"{network}_G3_energies_pH7.4.csv").open())
    }
    keep = set()
    for identifier in rels["Index"].astype(str):
        row = energies.get(identifier)
        estimable = row is not None and row["estimable"] == "True"
        spontaneous = estimable and row[SPONTANEITY] == "True"
        if spontaneous or (not estimable and basis == "with_unestimable"):
            keep.add(identifier)
    return keep


def main() -> None:
    t1, t2 = [], []
    for network, products_file in PRODUCTS.items():
        rels = pd.read_csv(
            REPO / "ProcessedData" / "RelsFiles" / network / f"{network}G3ProcessedRels.tsv",
            sep="\t",
        )
        rels["Index"] = rels["Index"].astype(str)
        products = read_products(
            REPO / "OriginalData" / "OriginalNetworkData" / "Products" / products_file
        )
        seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
        for basis in ("estimable_only", "with_unestimable"):
            keep = admitted(network, rels, basis)
            index = build_index(rels[rels["Index"].isin(keep)])
            depth = shortest_pathways(index, seeds)
            count = shortest_pathways(index, seeds, objective="reactions")
            for name, smiles in TARGETS.items():
                ok = smiles in depth.cost
                t2.append(
                    {
                        "network": network,
                        "basis": basis,
                        "target": name,
                        "reachable": "yes" if ok else "no",
                        "chain_depth": depth.cost.get(smiles, ""),
                        "reactions": count.cost.get(smiles, ""),
                    }
                )
                t1.append(
                    {
                        "network": network,
                        "basis": basis,
                        "target": name,
                        "reachable": "yes" if ok else "no",
                        "minimal_routes": count_minimal_routes(depth, index, smiles) if ok else "",
                        "critical_reactions": (
                            len(critical_reactions(index, seeds, smiles)) if ok else ""
                        ),
                    }
                )
            print(f"  {network:12s} {basis:16s} reactions admitted {len(keep):5,d}")

    for rows, name in ((t1, "SI_Table1_routes.csv"), (t2, "SI_Table2_steps.csv")):
        with (SI / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print("wrote", (SI / name).relative_to(REPO))


if __name__ == "__main__":
    main()
