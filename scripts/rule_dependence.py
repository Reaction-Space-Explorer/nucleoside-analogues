"""How far the traced routes depend on mechanisms currently in dispute.

    uv run python scripts/rule_dependence.py

Sutton et al. (Chem 2025, 11, 102553) find by 13C NMR that formaldehyde aldol
addition dominates the formose reaction and that carbonyl migration is absent,
and question the Breslow autocatalytic pathway on the same evidence. Tabata et
al. (Chem. Sci. 2023, 14, 13475) build an autocatalytic sugar cycle under
neutral conditions that requires the aldose-ketose transformation and that
suppresses the Cannizzaro reaction, which consumes both formaldehyde and
sugars under base.

The rule set used here encodes both freely. This removes each family from the
spontaneous network and re-searches every target, so that the results can be
read against either position.

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

MIGRATION = {"Keto-enol migration twice", "Elimination + enol to keto"}
FAMILIES = {
    "carbonyl migration": lambda rule: rule in MIGRATION,
    "Cannizzaro": lambda rule: "Cannizarro" in rule,
    "both": lambda rule: rule in MIGRATION or "Cannizarro" in rule,
}


def main() -> None:
    rows = []
    for network in PRODUCTS:
        generation = deepest(network)
        rels = pivot_rels(
            pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t")
        )
        rels["Index"] = rels["Index"].astype(str)
        products = read_products(
            REPO / "OriginalData" / "OriginalNetworkData" / "Products" / PRODUCTS[network]
        )
        seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
        keep = rels[rels["Index"].isin(admitted(network, rels, "estimable_only", generation))]
        full = shortest_pathways(build_index(keep), seeds).cost
        reachable = sum(1 for s in TARGETS.values() if s in full)
        for family, matches in FAMILIES.items():
            removed = keep["Rule"].map(matches)
            cut = shortest_pathways(build_index(keep[~removed]), seeds).cost
            survived = 0
            for name, smiles in TARGETS.items():
                a, b = full.get(smiles), cut.get(smiles)
                survived += b is not None
                rows.append(
                    {
                        "network": network,
                        "generation": generation,
                        "removed": family,
                        "target": name,
                        "depth_all_rules": a if a is not None else "",
                        "depth_after_removal": b if b is not None else "",
                        "survives": "yes" if b is not None else ("no" if a is not None else "n/a"),
                        "reactions_removed": int(removed.sum()),
                        "spontaneous_reactions": len(keep),
                    }
                )
            print(
                f"  {network:12s} without {family:18s} {int(removed.sum()):7,d} removed of "
                f"{len(keep):7,d};  {survived}/{reachable} reachable targets survive",
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
