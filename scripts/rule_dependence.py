"""How far the traced routes depend on carbonyl migration.

    uv run python scripts/rule_dependence.py

Sutton et al. (Chem 2025, 11, 102553) find by 13C NMR that formaldehyde aldol
addition dominates the formose reaction, that tetroses and pentoses are not
observed, and that carbonyl migration is therefore absent; they question the
Breslow autocatalytic pathway on the same evidence. The rule set used here
encodes carbonyl migration freely, so this asks which results survive without
it: every target is re-searched over the spontaneous network with the
migration rules removed.

Writes ProcessedData/SI/rule_dependence.csv.
"""

import csv
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import PRODUCTS, RELS, REPO, SI, TARGETS, admitted, deepest

from nucleoside_analogues.hyperpath import shortest_pathways
from nucleoside_analogues.rels import build_index, pivot_rels, read_products

#: The two rules that move a carbonyl along the chain.
MIGRATION = {"Keto-enol migration twice", "Elimination + enol to keto"}


def main() -> None:
    rows = []
    for network in PRODUCTS:
        generation = deepest(network)
        rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
        rels["Index"] = rels["Index"].astype(str)
        products = read_products(
            REPO / "OriginalData" / "OriginalNetworkData" / "Products" / PRODUCTS[network]
        )
        seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
        keep = rels[rels["Index"].isin(admitted(network, rels, "estimable_only", generation))]
        without = keep[~keep["Rule"].isin(MIGRATION)]
        full = shortest_pathways(build_index(keep), seeds).cost
        cut = shortest_pathways(build_index(without), seeds).cost
        share = len(keep) - len(without)
        for name, smiles in TARGETS.items():
            a, b = full.get(smiles), cut.get(smiles)
            rows.append(
                {
                    "network": network,
                    "generation": generation,
                    "target": name,
                    "depth_all_rules": a if a is not None else "",
                    "depth_without_migration": b if b is not None else "",
                    "survives": "yes" if b is not None else ("no" if a is not None else "n/a"),
                    "migration_reactions": share,
                    "spontaneous_reactions": len(keep),
                }
            )
        print(
            f"  {network:12s} {share:7,d} of {len(keep):7,d} spontaneous reactions are migrations; "
            f"targets surviving: "
            f"{sum(1 for r in rows[-len(TARGETS) :] if r['survives'] == 'yes')}/{len(TARGETS)}",
            flush=True,
        )
    out = SI / "rule_dependence.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
