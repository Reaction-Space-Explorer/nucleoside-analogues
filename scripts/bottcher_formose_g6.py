"""Böttcher complexity for the Formose G6 matches.

    uv run python scripts/bottcher_formose_g6.py [workers]

ComplexityData is MatchesFiles with a Complexity column, so extending the
matching to generation six left the complexity file short by the same 5,901
species. This is the notebook implementation
(notebooks/pipeline/BottcherComplexity.ipynb) applied to those rows.

It is verified before it is used: recomputing the deposited G1-G5 rows must
reproduce their complexities, or the run aborts.
"""

import math
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import RDConfig

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.append(os.path.join(RDConfig.RDContribDir, "ChiralPairs"))

import ChiralDescriptors  # noqa: E402
from make_si_tables import REPO  # noqa: E402

RDLogger.DisableLog("rdApp.*")
MATCHES = REPO / "ProcessedData" / "MatchesFiles" / "FormoseMatches.tsv"
COMPLEXITY = REPO / "ProcessedData" / "ComplexityData" / "FormoseComplexityData.tsv"

VALENCE = {
    1: "H Li Na K Rb Cs Fr",
    2: "Be Mg Ca Sr Ba Ra",
    3: "B Al Ga In Tl Nh",
    4: "C Si Ge Sn Pb Fl",
    5: "N P As Sb Bi Mc",
    6: "O S Se Te Po Lv",
    7: "F Cl Br I At Ts",
    8: "He Ne Ar Kr Xe Rn Og",
}
VALENCE = {s: k for k, v in VALENCE.items() for s in v.split()}


def nonequivs(atom, mol) -> int:
    subs = [[] for _ in range(8)]
    found = ChiralDescriptors.determineAtomSubstituents(
        atom.GetIdx(), mol, Chem.GetDistanceMatrix(mol)
    )[0]
    for item, key in enumerate(found):
        for sub in found[key]:
            subs[item].append(mol.GetAtomWithIdx(sub).GetSymbol())
    return len({tuple(s) for s in subs if s})


def local_diversity(atom) -> int:
    neighbours = {n.GetSymbol() for n in atom.GetNeighbors()}
    return len(neighbours) if atom.GetSymbol() in neighbours else len(neighbours) + 1


def isomeric(atom) -> int:
    try:
        return 2 if atom.GetProp("_CIPCode") else 1
    except KeyError:
        return 1


def bond_index(atom) -> int:
    bonds = [str(b.GetBondType()) for b in atom.GetBonds()]
    rank = sum({"SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3}.get(b, 0) for b in bonds)
    if "AROMATIC" in bonds:
        rank += 3 if atom.GetSymbol() == "C" else 2 if atom.GetSymbol() in ("N", "S") else 0
    return rank


def split_atoms(mol, ranks: list[str]):
    """Symmetry-equivalent atoms count once each per pair, the odd one at half."""
    full, half, remaining = [], [], list(ranks)
    for atom in mol.GetAtoms():
        n = remaining.count(atom.GetProp("_CIPRank"))
        if n == 1:
            full.append(atom)
            continue
        full.extend([atom] * (n // 2))
        half.extend([atom] * (n % 2))
        remaining = [r for r in remaining if r != atom.GetProp("_CIPRank")]
    return full, half


def complexity(smiles: str) -> float | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    ranks = [a.GetProp("_CIPRank") for a in mol.GetAtoms()]
    full, half = split_atoms(mol, ranks)
    total = 0.0
    for atoms, weight in ((full, 1.0), (half, 0.5)):
        for atom in atoms:
            d = nonequivs(atom, mol)
            e = local_diversity(atom)
            s = isomeric(atom)
            v = VALENCE.get(atom.GetSymbol(), 0)
            b = bond_index(atom)
            if v * b <= 0:
                continue
            total += weight * d * e * s * math.log(v * b, 2)
    return total


def main() -> None:
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    deposited = pd.read_csv(COMPLEXITY, sep="\t", index_col=0)
    matches = pd.read_csv(MATCHES, sep="\t")

    check = deposited.head(400)
    with Pool(workers) as pool:
        got = pool.map(complexity, check["NetworkSmiles"].astype(str).tolist(), chunksize=20)
    bad = [
        (s, a, b)
        for s, a, b in zip(check["NetworkSmiles"], check["Complexity"], got, strict=True)
        if b is None or abs(float(a) - b) > 1e-6
    ]
    if bad:
        print(f"  ABORT: {len(bad)} of {len(check)} deposited complexities not reproduced")
        for s, a, b in bad[:3]:
            print(f"    {s}  deposited {a}  recomputed {b}")
        raise SystemExit(1)
    print(f"  method verified on {len(check)} deposited rows", flush=True)

    new = matches[matches["Generation"] == "G6"].reset_index(drop=True)
    print(f"  computing {len(new):,} G6 complexities on {workers} workers", flush=True)
    with Pool(workers) as pool:
        values = pool.map(complexity, new["NetworkSmiles"].astype(str).tolist(), chunksize=20)
    new["Complexity"] = values
    missing = new["Complexity"].isna().sum()
    if missing:
        print(f"  {missing} could not be computed and are left out")
        new = new.dropna(subset=["Complexity"])

    out = pd.concat([deposited, new], ignore_index=True)
    out.to_csv(COMPLEXITY, sep="\t")
    print(f"  wrote {len(out):,} rows to ProcessedData/ComplexityData/FormoseComplexityData.tsv")


if __name__ == "__main__":
    main()
