"""Match the Formose G6 species against the nucleoside analogue library.

    uv run python scripts/match_formose_g6.py [workers]

The Formose network's reactions run to generation six, but its deposited
product list stops at five, so the analogue matching stopped there too. This
matches the 94,415 species that appear only in FormoseRels_6.

The library is ProcessedData/Nucleoside_Stereoisomers.tsv, which already
carries the InChIKey first block against the enumerated stereoisomers; that
first block is what a match is. The method is verified before it is used: it
must reproduce the deposited G1-G5 matches exactly, or the run aborts.
"""

import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import RELS, REPO

from nucleoside_analogues.rels import pivot_rels, read_products

RDLogger.DisableLog("rdApp.*")
OUT = REPO / "ProcessedData" / "MatchesFiles"
LIBRARY = REPO / "ProcessedData" / "Nucleoside_Stereoisomers.tsv"


def key14(smiles: str) -> tuple[str, str] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    try:
        return smiles, Chem.MolToInchiKey(mol)[:14]
    except Exception:  # noqa: BLE001 - one bad structure must not stop the run
        return None


def main() -> None:
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    library = pd.read_csv(LIBRARY, sep="\t")
    analogues = dict(
        zip(library["INCHIKEY"].astype(str), library["SMILES"].astype(str), strict=True)
    )
    print(f"  analogue library: {len(analogues):,} skeletons", flush=True)

    products = read_products(
        REPO / "OriginalData" / "OriginalNetworkData" / "Products" / "formose_output.tsv"
    )
    known = set(products["Smiles"].astype(str))

    rels = pivot_rels(pd.read_csv(RELS / "Formose" / "FormoseRels_6.tsv", sep="\t"))
    in_rels = set()
    for column in ("Reagents", "Products"):
        for row in rels[column]:
            in_rels.update(row)
    g6 = sorted(in_rels - known)
    print(f"  species only in Rels_6: {len(g6):,}", flush=True)

    def match(species):
        with Pool(workers) as pool:
            done = pool.map(key14, sorted(species), chunksize=500)
        return [(k, s) for r in done if r for s, k in [r] if k in analogues]

    deposited = pd.read_csv(OUT / "FormoseMatches.tsv", sep="\t")
    want = set(
        zip(deposited["INCHIKEY"].astype(str), deposited["NetworkSmiles"].astype(str), strict=True)
    )
    got = set(match(known))
    if got != want:
        print(
            f"  ABORT: reproduced {len(got):,} of the deposited {len(want):,}; "
            f"missing {len(want - got):,}, extra {len(got - want):,}"
        )
        raise SystemExit(1)
    print(f"  method verified: reproduces all {len(want):,} deposited G1-G5 matches", flush=True)

    rows = [
        {"Generation": "G6", "INCHIKEY": k, "NetworkSmiles": s, "AnalogueSmiles": analogues[k]}
        for k, s in sorted(match(g6))
    ]
    print(f"  new G6 matches: {len(rows):,}", flush=True)
    frame = pd.concat([deposited, pd.DataFrame(rows)], ignore_index=True)
    frame.to_csv(OUT / "FormoseMatches.tsv", sep="\t", index=False)
    print(f"  wrote {len(frame):,} rows to ProcessedData/MatchesFiles/FormoseMatches.tsv")


if __name__ == "__main__":
    main()
