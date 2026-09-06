"""Are the canonical nucleosides distinguished among comparable analogues?

    uv run python scripts/matched_controls.py

Comparing four targets against several thousand matched analogues would confound
the answer with size and composition. Following Wolos et al., each target is
instead compared against controls drawn from the same network: the analogues
most similar to it by Morgan-fingerprint Tanimoto, and those closest to it in
exact mass. The question is then whether the target is reached in fewer steps
than structures the network finds comparably easy to make.
"""

import csv
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Descriptors import ExactMolWt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import PRODUCTS, RELS, SI, TARGETS, admitted, deepest

from nucleoside_analogues.hyperpath import shortest_pathways
from nucleoside_analogues.rels import build_index, pivot_rels, read_products

RDLogger.DisableLog("rdApp.*")
REPO = Path(__file__).resolve().parent.parent
#: Controls per target, from each of the two matching criteria.
N_CONTROLS = 25
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


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
        index = build_index(
            rels[rels["Index"].isin(admitted(network, rels, "with_unestimable", generation))]
        )
        depth = shortest_pathways(index, seeds).cost

        matches = pd.read_csv(REPO / "ProcessedData" / "MatchesFiles" / f"{network}Matches.tsv",
                              sep="\t")
        pool = []
        for smiles in dict.fromkeys(matches["NetworkSmiles"].astype(str)):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None or smiles not in depth:
                continue
            pool.append((smiles, GEN.GetFingerprint(mol), ExactMolWt(mol), depth[smiles]))

        for name, target in TARGETS.items():
            mol = Chem.MolFromSmiles(target)
            if target not in depth:
                continue
            fp, mass = GEN.GetFingerprint(mol), ExactMolWt(mol)
            others = [p for p in pool if p[0] != target]
            by_tan = sorted(others, key=lambda p: -DataStructs.TanimotoSimilarity(fp, p[1]))
            by_mass = sorted(others, key=lambda p: abs(p[2] - mass))
            controls = {p[0]: p for p in by_tan[:N_CONTROLS]}
            controls.update({p[0]: p for p in by_mass[:N_CONTROLS]})
            depths = sorted(p[3] for p in controls.values())
            target_depth = depth[target]
            rank = sum(1 for d in depths if d < target_depth)
            rows.append({
                "network": network, "target": name, "target_steps": target_depth,
                "controls": len(controls),
                "control_median_steps": depths[len(depths) // 2],
                "control_min_steps": depths[0], "control_max_steps": depths[-1],
                "controls_strictly_faster": rank,
                "percentile": round(100 * rank / len(depths), 1),
                "mean_tanimoto_of_controls": round(
                    sum(DataStructs.TanimotoSimilarity(fp, p[1]) for p in controls.values())
                    / len(controls), 3),
            })
            print(f"  {network:12s} {name:12s} target {target_depth} steps | controls "
                  f"n={len(controls)} median {depths[len(depths)//2]} "
                  f"[{depths[0]}-{depths[-1]}] | {rank} faster ({rows[-1]['percentile']}%)",
                  flush=True)

    with (SI / "matched_controls.csv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("wrote ProcessedData/SI/matched_controls.csv")


if __name__ == "__main__":
    main()
