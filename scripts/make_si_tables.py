"""Regenerate SI Tables 1 and 2, at each network's deepest generation.

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
from nucleoside_analogues.rels import build_index, pivot_rels, read_products

REPO = Path(__file__).resolve().parent.parent
SI = REPO / "ProcessedData" / "SI"
#: Full-depth energies live in a subdirectory; the generation-3 files remain alongside.
FULL = SI / "full"
RELS = REPO / "OriginalData" / "OriginalNetworkData" / "Rels"
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


#: A reaction counts as spontaneous only where the whole 95% interval lies
#: below zero, matching the Methods. The tables are identical on the point
#: estimate, so nothing here depends on the choice.
Z = 1.96

#: Component contribution returns exactly 0 +/- 0 where reagents and products
#: share a decomposition, as keto-enol migrations do. That is an unresolved
#: estimate, not a spontaneous one, and a point value of -1e-05 must not pass
#: the test on rounding. Such reactions are counted unestimable.
NULL_DG = 1e-3


def is_null(row: dict) -> bool:
    try:
        return float(row["sigma_kJ_mol"]) == 0.0 and abs(float(row["dG_prime_kJ_mol"])) < NULL_DG
    except (ValueError, KeyError):
        return False


def is_spontaneous(row: dict) -> bool:
    """dG + 1.96 sigma < 0, computed rather than read.

    The generation-3 files carry precomputed spontaneity columns and the
    full-depth files do not, so deriving it keeps one definition for both.
    """
    if row.get("spontaneous_95") is not None:
        return row["spontaneous_95"] == "True"
    try:
        return float(row["dG_prime_kJ_mol"]) + Z * float(row["sigma_kJ_mol"]) < 0
    except (ValueError, KeyError):
        return False


def species_generations(network: str) -> dict[str, int]:
    """Every species in the network, and the generation it first appears in.

    Taken from the rels rather than the product listing: Formose's listing
    stops at generation five while its reactions run to six, so the 94,415
    species of that last tranche are absent from it. They take the deepest
    generation, which is where they appear.
    """
    generation = deepest(network)
    products = read_products(
        REPO / "OriginalData" / "OriginalNetworkData" / "Products" / PRODUCTS[network]
    )
    known = {
        str(s): int(g) for s, g in zip(products["Smiles"], products["Generation"], strict=True)
    }
    rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
    for column in ("Reagents", "Products"):
        for row in rels[column]:
            for smiles in row:
                known.setdefault(str(smiles), generation)
    return known


def deepest(network: str) -> int:
    return max(int(f.stem.split("_")[-1]) for f in (RELS / network).glob("*Rels_*.tsv"))


def energy_file(network: str, generation: int) -> Path:
    """Full-depth energies where present, else the generation-3 set."""
    full = FULL / f"{network}_G{generation}_energies_pH7.4.csv"
    return full if full.exists() else SI / f"{network}_G3_energies_pH7.4.csv"


def admitted(network: str, rels: pd.DataFrame, basis: str, generation: int) -> set[str]:
    """Reaction ids allowed under a basis. Reactions with no energy row could
    not be built at all, and count as unestimable rather than vanishing."""
    energies = {
        row["Index"]: row for row in csv.DictReader(energy_file(network, generation).open())
    }
    keep = set()
    for identifier in rels["Index"].astype(str):
        row = energies.get(identifier)
        estimable = row is not None and row["estimable"] == "True" and not is_null(row)
        spontaneous = estimable and is_spontaneous(row)
        if spontaneous or (not estimable and basis == "with_unestimable"):
            keep.add(identifier)
    return keep


def main() -> None:
    t1, t2 = [], []
    for network, products_file in PRODUCTS.items():
        generation = deepest(network)
        rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
        rels["Index"] = rels["Index"].astype(str)
        products = read_products(
            REPO / "OriginalData" / "OriginalNetworkData" / "Products" / products_file
        )
        seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
        for basis in ("estimable_only", "with_unestimable"):
            keep = admitted(network, rels, basis, generation)
            index = build_index(rels[rels["Index"].isin(keep)])
            depth = shortest_pathways(index, seeds)
            count = shortest_pathways(index, seeds, objective="reactions")
            for name, smiles in TARGETS.items():
                ok = smiles in depth.cost
                t2.append(
                    {
                        "network": network,
                        "generation": generation,
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
                        "generation": generation,
                        "basis": basis,
                        "target": name,
                        "reachable": "yes" if ok else "no",
                        "minimal_routes": count_minimal_routes(depth, index, smiles) if ok else "",
                        "critical_reactions": (
                            len(critical_reactions(index, seeds, smiles)) if ok else ""
                        ),
                    }
                )
            print(f"  {network:12s} G{generation} {basis:16s} admitted {len(keep):7,d}")

    for rows, name in ((t1, "SI_Table1_routes.csv"), (t2, "SI_Table2_steps.csv")):
        with (SI / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print("wrote", (SI / name).relative_to(REPO))


if __name__ == "__main__":
    main()
