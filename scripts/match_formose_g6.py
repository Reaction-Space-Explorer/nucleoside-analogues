"""Match the Formose G6 species against the nucleoside analogue libraries.

    uv run python scripts/match_formose_g6.py

The Formose network's reactions run to generation six, but its deposited
product list stops at five, so the analogue matching stopped there too. This
matches the 94,415 species that appear only in FormoseRels_6.

The method is verified before it is used: it must reproduce the deposited
G1-G5 matches exactly, or the run aborts. Matching is on the first 14
characters of the InChIKey, as everywhere else here.
"""

import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import RELS, REPO

from nucleoside_analogues.rels import pivot_rels, read_products

RDLogger.DisableLog("rdApp.*")
LIB = REPO / "OriginalData" / "OriginalNucleosideAnalogueData"
OUT = REPO / "ProcessedData" / "MatchesFiles"


def key14(smiles: str) -> tuple[str, str] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    try:
        return smiles, Chem.MolToInchiKey(mol)[:14]
    except Exception:  # noqa: BLE001 - one bad structure must not stop the library
        return None


def keys(smiles, workers: int, chunk: int = 2000):
    with Pool(workers) as pool:
        for r in pool.imap_unordered(key14, smiles, chunksize=chunk):
            if r is not None:
                yield r


def main() -> None:
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    analogues: dict[str, set[str]] = defaultdict(set)
    for name in ("CHO_Smiles.tsv", "CHNO_Smiles.tsv"):
        smiles = pd.read_csv(LIB / name, sep="\t")["SMILES"].astype(str).tolist()
        print(f"  {name}: {len(smiles):,} structures", flush=True)
        for smi, k in keys(smiles, workers):
            analogues[k].add(smi)
    print(f"  analogue library: {len(analogues):,} distinct skeletons", flush=True)

    products = read_products(
        REPO / "OriginalData" / "OriginalNetworkData" / "Products" / "formose_output.tsv"
    )
    known = {
        str(s): int(g) for s, g in zip(products["Smiles"], products["Generation"], strict=True)
    }

    rels = pivot_rels(pd.read_csv(RELS / "Formose" / "FormoseRels_6.tsv", sep="\t"))
    in_rels = set()
    for column in ("Reagents", "Products"):
        for row in rels[column]:
            in_rels.update(row)
    g6 = sorted(in_rels - set(known))
    print(f"  species only in Rels_6: {len(g6):,}", flush=True)

    def match(smiles_list):
        out = []
        for smi, k in keys(smiles_list, workers):
            if k in analogues:
                out.append((k, smi, analogues[k]))
        return out

    # verify the method against the deposited matches before extending anything
    deposited = pd.read_csv(OUT / "FormoseMatches.tsv", sep="\t")
    want = set(
        zip(deposited["INCHIKEY"].astype(str), deposited["NetworkSmiles"].astype(str), strict=True)
    )
    got = {(k, s) for k, s, _ in match(sorted(known))}
    if got != want:
        print(
            f"  ABORT: reproduced {len(got):,} of the deposited {len(want):,} matches; "
            f"missing {len(want - got):,}, extra {len(got - want):,}"
        )
        raise SystemExit(1)
    print(f"  method verified: reproduces all {len(want):,} deposited G1-G5 matches", flush=True)

    rows = [
        {"Generation": "G6", "INCHIKEY": k, "NetworkSmiles": s, "AnalogueSmiles": str(a)}
        for k, s, a in match(g6)
    ]
    print(f"  new G6 matches: {len(rows):,}", flush=True)
    frame = pd.concat([deposited, pd.DataFrame(rows)], ignore_index=True)
    frame.to_csv(OUT / "FormoseMatches.tsv", sep="\t", index=False)
    print(f"  wrote {len(frame):,} rows to ProcessedData/MatchesFiles/FormoseMatches.tsv")


if __name__ == "__main__":
    main()
