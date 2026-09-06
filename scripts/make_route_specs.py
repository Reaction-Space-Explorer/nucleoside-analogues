"""Emit autocycle route specs for the traced pathway to each target.

    uv run python scripts/make_route_specs.py

One YAML per network and target, at the network's deepest generation, over the
spontaneous-only network -- the same basis as SI Table 1. Draw them with

    autocycle route figures/routes/Formose_Ribose.yaml --style rich --drop O -o ribose.pdf

The rich style scales each arrow by the magnitude of its free energy and
colours it by sign, and prints the route total. autocycle is not a dependency
here; it reads the YAML this writes.
"""

import csv
import sys
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import PRODUCTS, RELS, REPO, TARGETS, admitted, deepest, energy_file

from nucleoside_analogues.hyperpath import shortest_pathways
from nucleoside_analogues.rels import build_index, pivot_rels, read_products

RDLogger.DisableLog("rdApp.*")
OUT = REPO / "figures" / "routes"

#: Names for the molecules a reader is expected to recognise. The networks
#: carry no stereocentres, so a node with more than one diastereomer is named
#: for its constitution -- "tetrose", not "threose" -- and the target keeps the
#: manuscript's name. Anything absent is drawn unlabelled rather than guessed.
NAMES = {
    "C=O": "formaldehyde",
    "OCC=O": "glycolaldehyde",
    "OCC(O)C=O": "glyceraldehyde",
    "OCC(=O)CO": "dihydroxyacetone",
    "OCCO": "ethylene glycol",
    "CC=O": "acetaldehyde",
    "OCC(O)C(O)C=O": "aldotetrose",
    "OCC(O)C(O)C(O)C=O": "aldopentose",
    "OCC(O)C(O)CC=O": "2-deoxyaldopentose",
    "OCC(O)C(O)CO": "tetritol",
    "OCC(O)CO": "glycerol",
    "O": "water",
    "OC=O": "formic acid",
    "CC(=O)C(=O)O": "pyruvic acid",
    "OCC(O)C(O)C(O)C(O)C=O": "aldohexose",
}


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else smiles


NAMES = {canonical(k): v for k, v in NAMES.items()}


def node(species, result, index, rules, energies, path, label=None):
    """One route node, recursing through the optimal back-pointers."""
    mol = {"smiles": species}
    name = label or NAMES.get(canonical(species))
    if name:
        mol["label"] = name

    if result.cost.get(species) == 0:
        return {"mol": mol, "seed": True}
    identifier = result.via.get(species)
    if identifier is None or species in path:
        return {"mol": mol, "terminal": "unknown"}

    reaction = {"id": identifier, "rule": rules.get(identifier, "")}
    row = energies.get(identifier)
    if row and row["estimable"] == "True":
        reaction["dg"] = round(float(row["dG_prime_kJ_mol"]), 1)
    coproducts = [p for p in index.products.get(identifier, ()) if p != species]
    if coproducts:
        reaction["produces"] = coproducts

    return {
        "mol": mol,
        "reaction": reaction,
        "from": [
            node(r, result, index, rules, energies, path | {species})
            for r in index.reagents[identifier]
        ],
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    written = 0
    for network, products_file in PRODUCTS.items():
        generation = deepest(network)
        rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
        rels["Index"] = rels["Index"].astype(str)
        rules = dict(zip(rels["Index"], rels["Rule"], strict=True))
        energies = {
            row["Index"]: row for row in csv.DictReader(energy_file(network, generation).open())
        }
        products = read_products(
            REPO / "OriginalData" / "OriginalNetworkData" / "Products" / products_file
        )
        seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
        keep = admitted(network, rels, "estimable_only", generation)
        index = build_index(rels[rels["Index"].isin(keep)])
        result = shortest_pathways(index, seeds)

        for target, smiles in TARGETS.items():
            if smiles not in result.cost:
                continue
            spec = {
                "title": (
                    f"{target} from {network}, G{generation}, spontaneous reactions only. "
                    "Structures are constitutional; the network carries no stereochemistry."
                ),
                "target": node(smiles, result, index, rules, energies, frozenset(), target),
            }
            path = OUT / f"{network}_{target}.yaml"
            path.write_text(yaml.safe_dump(spec, sort_keys=False, width=100))
            written += 1
            print(f"  {path.relative_to(REPO)}  depth {result.cost[smiles]}")
    print(f"wrote {written} specs to {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
