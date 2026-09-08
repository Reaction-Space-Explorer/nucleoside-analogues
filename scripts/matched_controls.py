"""Are the canonical nucleosides distinguished among comparable analogues?

    uv run python scripts/matched_controls.py

Comparing four targets against several thousand matched analogues would confound
the answer with size and composition. Following Wolos et al., each target is
instead compared against controls drawn from the same network: the analogues
most similar to it by Morgan-fingerprint Tanimoto, and those closest to it in
exact mass. The question is then whether the target is reached in fewer steps
than structures the network finds comparably easy to make.

Both admission bases are reported. The spontaneous-only basis is the one used
everywhere else in this work, but on it the PA CRNR reaches none of the five
targets, so its comparison would be empty; admitting reactions of unestimable
free energy brings PA in. The conclusion does not depend on the choice, which is
why both are written out rather than one being picked silently.

Significance is the conservative rank estimator of Phipson and Smyth, p =
(r + 1) / (n + 1) for r controls strictly faster than the target out of n, which
cannot return zero and is bounded below by 1/(n + 1); with fifty controls the
smallest attainable value is 0.02. The per-pair values are combined by Fisher's
method, and the formose CRNRs are compared against the rest by enumerating all
C(25, 10) = 3,268,760 splits, so the permutation p is exact rather than sampled.
"""

import csv
import itertools
import math
import sys
from pathlib import Path

import numpy as np
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


BASES = ("estimable_only", "with_unestimable")


def main() -> None:
    rows = []
    for basis in BASES:
        for network, products_file in PRODUCTS.items():
            generation = deepest(network)
            rels = pivot_rels(
                pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t")
            )
            rels["Index"] = rels["Index"].astype(str)
            products = read_products(
                REPO / "OriginalData" / "OriginalNetworkData" / "Products" / products_file
            )
            seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
            index = build_index(
                rels[rels["Index"].isin(admitted(network, rels, basis, generation))]
            )
            depth = shortest_pathways(index, seeds).cost

            matches = pd.read_csv(
                REPO / "ProcessedData" / "MatchesFiles" / f"{network}Matches.tsv", sep="\t"
            )
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
                rows.append(
                    {
                        "basis": basis,
                        "network": network,
                        "target": name,
                        "target_steps": target_depth,
                        "controls": len(controls),
                        "control_median_steps": depths[len(depths) // 2],
                        "control_min_steps": depths[0],
                        "control_max_steps": depths[-1],
                        "controls_strictly_faster": rank,
                        "p_value": round((rank + 1) / (len(controls) + 1), 4),
                        "percentile": round(100 * rank / len(depths), 1),
                        "mean_tanimoto_of_controls": round(
                            sum(DataStructs.TanimotoSimilarity(fp, p[1]) for p in controls.values())
                            / len(controls),
                            3,
                        ),
                    }
                )
                print(
                    f"  {network:12s} {name:12s} target {target_depth} steps | controls "
                    f"n={len(controls)} median {depths[len(depths) // 2]} "
                    f"[{depths[0]}-{depths[-1]}] | {rank} faster ({rows[-1]['percentile']}%)",
                    flush=True,
                )

    with (SI / "matched_controls.csv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote ProcessedData/SI/matched_controls.csv")
    for basis in BASES:
        write_statistics([r for r in rows if r["basis"] == basis], basis)


FORMOSE = ("Formose", "FormoseAmm")


def write_statistics(rows: list[dict], basis: str) -> None:
    """Fisher's combined p, and the exact permutation test separating the formose CRNRs."""
    values = np.array([r["p_value"] for r in rows])
    chi2 = -2 * float(np.log(values).sum())
    degrees = 2 * len(values)
    # survival function of chi2 with even df, in closed form, so that summarising
    # the result costs no dependency the rest of the analysis does not already have
    half = chi2 / 2
    term = 1.0
    total = 1.0
    for k in range(1, degrees // 2):
        term *= half / k
        total += term
    combined = math.exp(-half) * total

    block = [i for i, r in enumerate(rows) if r["network"] in FORMOSE]
    k = len(block)
    observed = float(values[block].mean())
    if math.comb(len(values), k) > 5_000_000:
        raise SystemExit(f"{math.comb(len(values), k):,} splits is too many to enumerate")
    splits = np.array(list(itertools.combinations(range(len(values)), k)), dtype=np.int16)
    means = values[splits].sum(axis=1) / k
    at_least = int((means <= observed + 1e-12).sum())
    exact = at_least / len(means)

    summary = [
        {"statistic": "basis", "value": basis},
        {"statistic": "pairs", "value": len(values)},
        {"statistic": "fisher_chi2", "value": round(chi2, 2)},
        {"statistic": "fisher_df", "value": degrees},
        {"statistic": "fisher_combined_p", "value": f"{combined:.2e}"},
        {"statistic": "mean_p_formose", "value": round(observed, 4)},
        {"statistic": "mean_p_other", "value": round(float(np.delete(values, block).mean()), 4)},
        {"statistic": "permutation_splits", "value": len(means)},
        {"statistic": "permutation_p_exact", "value": round(exact, 5)},
    ]
    path = SI / f"control_statistics_{basis}.csv"
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["statistic", "value"])
        w.writeheader()
        w.writerows(summary)
    print(f"  --- {basis} ---")
    for row in summary:
        print(f"  {row['statistic']:22s} {row['value']}")
    print(f"wrote {path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
