"""The hexose branch point, against the experiment of Yi et al.

    uv run python scripts/ketohexose_case.py

Yi et al. condensed glycolaldehyde alone at pH 8.5 and found 2-ketohexoses as
the major six-carbon products, with no aldohexoses detected, and assigned the
route as aldol, carbonyl migration, aldol, carbonyl migration: two glycolaldehyde
to a tetrose, the tetrose carried to a tetrulose, aldol with a third
glycolaldehyde to a 3-ketohexose, and a final migration to the 2-ketohexose.

That is a test of minimum-step accessibility against a measured product
distribution, on the one feedstock both share, and it cuts both ways. This
script traces it and writes what it finds:

  * The route to the 3-ketohexose is theirs exactly -- aldol, migration, aldol --
    found by shortest-path search over the spontaneous network, from glycolaldehyde
    alone. Every elementary step they assign is present and admitted.

  * The route to the 2-ketohexose is not. Ours is a step shorter and runs through
    the aldohexose, the species they report absent. Their final step exists in the
    network but returns a null estimate, zero free energy with zero uncertainty,
    so the spontaneous network cannot express their mechanism at all.

Writes ProcessedData/SI/ketohexose_case.csv and the two autocycle specs.
"""

import csv
import sys
from pathlib import Path

import pandas as pd
import yaml
from rdkit import RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_route_specs import ROUTES, canonical, node
from make_si_tables import (
    PRODUCTS,
    RELS,
    REPO,
    SI,
    admitted,
    deepest,
    energy_file,
    is_null,
    species_generations,
)

from nucleoside_analogues.hyperpath import shortest_pathways
from nucleoside_analogues.rels import build_index, pivot_rels, read_products

RDLogger.DisableLog("rdApp.*")
NETWORK = "Formose"
CASES = {
    "3ketohexose": ("OCC(O)C(=O)C(O)C(O)CO", "3-ketohexose"),
    "2ketohexose": ("OCC(=O)C(O)C(O)C(O)CO", "2-ketohexose"),
    "aldohexose": ("OCC(O)C(O)C(O)C(O)C=O", "aldohexose"),
}


def main() -> None:
    generation = deepest(NETWORK)
    rels = pivot_rels(pd.read_csv(RELS / NETWORK / f"{NETWORK}Rels_{generation}.tsv", sep="\t"))
    rels["Index"] = rels["Index"].astype(str)
    rules = dict(zip(rels["Index"], rels["Rule"], strict=True))
    energies = {r["Index"]: r for r in csv.DictReader(energy_file(NETWORK, generation).open())}
    products = read_products(
        REPO / "OriginalData" / "OriginalNetworkData" / "Products" / PRODUCTS[NETWORK]
    )
    seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
    admitted_ids = set(admitted(NETWORK, rels, "estimable_only", generation))
    index = build_index(rels[rels["Index"].isin(admitted_ids)])
    result = shortest_pathways(index, seeds)
    # cost is keyed by the network's own SMILES spelling, not a canonical one
    spelling = {canonical(str(s)): str(s) for s in species_generations(NETWORK)}

    rows = []
    for key, (smiles, label) in CASES.items():
        target = spelling.get(canonical(smiles), canonical(smiles))
        depth = result.cost.get(target)
        rows.append(
            {
                "species": label,
                "smiles": target,
                "spontaneous_depth": depth if depth is not None else "unreached",
                "observed_by_Yi": {
                    "2-ketohexose": "major product",
                    "3-ketohexose": "assigned intermediate",
                    "aldohexose": "not detected",
                }[label],
            }
        )
        if depth:
            spec = {
                "title": (
                    f"{label} from {NETWORK}, G{generation}, spontaneous reactions only. "
                    "Structures are constitutional; the network carries no stereochemistry."
                ),
                "target": node(target, result, index, rules, energies, frozenset(), label),
            }
            (ROUTES / f"case_{key}.yaml").write_text(
                yaml.safe_dump(spec, sort_keys=False, width=100)
            )

    # their final step: 3-ketohexose -> 2-ketohexose
    three = canonical(CASES["3ketohexose"][0])
    two = canonical(CASES["2ketohexose"][0])
    for i, a, b in zip(rels["Index"], rels["Reagents"], rels["Products"], strict=True):
        if three in {canonical(str(x)) for x in a} and two in {canonical(str(x)) for x in b}:
            row = energies.get(i)
            rows.append(
                {
                    "species": "their final migration (3-keto to 2-keto)",
                    "smiles": i,
                    "spontaneous_depth": "admitted" if i in admitted_ids else "excluded",
                    "observed_by_Yi": (
                        f"rule {rules[i]}; dG {row['dG_prime_kJ_mol']}, sigma "
                        f"{row['sigma_kJ_mol']}, null estimate {is_null(row)}"
                    ),
                }
            )

    nulls = [r for r in energies.values() if r["estimable"] == "True" and is_null(r)]
    migration = sum(1 for r in nulls if rules[r["Index"]] == "Keto-enol migration twice")
    rows.append(
        {
            "species": "null estimates in this network",
            "smiles": "",
            "spontaneous_depth": len(nulls),
            "observed_by_Yi": f"{migration} of them are Keto-enol migration twice",
        }
    )

    out = SI / "ketohexose_case.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for r in rows:
        print(f"  {r['species']:44s} {str(r['spontaneous_depth']):>10s}  {r['observed_by_Yi']}")
    print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
