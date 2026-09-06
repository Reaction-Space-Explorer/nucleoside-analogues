"""Regenerate the analogue/network matches and diff against the deposited files.

    uv run python scripts/verify_matches.py

The deposited MatchesFiles were used throughout without ever being re-derived.
This rebuilds them from the raw library and the network product listings and
reports any disagreement, rather than assuming they are right.

Each library entry carries a placeholder Cl that stands for either an OH or an
NH2 substituent, so every entry yields two molecules. Matching is on the first
14 InChIKey characters, which encode constitution and charge but not
stereochemistry.
"""

import csv
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import species_generations  # noqa: E402

RDLogger.DisableLog("rdApp.*")
REPO = Path(__file__).resolve().parent.parent
LIBRARY = REPO / "OriginalData" / "OriginalNucleosideAnalogueData"
PRODUCTS = REPO / "OriginalData" / "OriginalNetworkData" / "Products"
MATCHES = REPO / "ProcessedData" / "MatchesFiles"
OUT = REPO / "ProcessedData" / "SI" / "matches_verification.csv"
NETWORKS = {
    "Formose": "formose_output.tsv",
    "FormoseAmm": "formose_amm_output.tsv",
    "Glucose": "glucose_degradation_output.tsv",
    "GlucoseAmm": "glucose_amm_output.tsv",
    "PyruvicAcid": "pyruvic_output.tsv",
}
SUBSTITUTIONS = ("O[H]", "N([H])[H]")


def key(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    return None if mol is None else Chem.MolToInchiKey(mol)[:14]


def library_keys() -> set[str]:
    keys: set[str] = set()
    for name in ("CHNO_Smiles.tsv", "CHO_Smiles.tsv"):
        frame = pd.read_csv(LIBRARY / name, sep="\t")
        for n, entry in enumerate(frame["SMILES"], 1):
            for replacement in SUBSTITUTIONS:
                k = key(entry.replace("Cl", replacement))
                if k is not None:
                    keys.add(k)
            if n % 50_000 == 0:
                print(f"  {name}: {n:,} entries, {len(keys):,} unique keys", flush=True)
        print(f"  {name}: done, {len(keys):,} unique keys", flush=True)
    return keys


def main() -> None:
    print("building library keys (both substitutions, deduplicated)", flush=True)
    keys = library_keys()

    rows = []
    for network in NETWORKS:
        first: dict[str, int] = {}
        for smiles, generation in species_generations(network).items():
            k = key(str(smiles))
            if k is not None and (k not in first or generation < first[k]):
                first[k] = int(generation)
        mine = {k: g for k, g in first.items() if k in keys}

        deposited_frame = pd.read_csv(MATCHES / f"{network}Matches.tsv", sep="\t")
        deposited = {
            str(k): int(str(g)[1:])
            for k, g in zip(deposited_frame["INCHIKEY"], deposited_frame["Generation"], strict=True)
        }
        only_mine = set(mine) - set(deposited)
        only_dep = set(deposited) - set(mine)
        gen_diff = {k for k in set(mine) & set(deposited) if mine[k] != deposited[k]}
        rows.append(
            {
                "network": network,
                "regenerated": len(mine),
                "deposited": len(deposited),
                "only_regenerated": len(only_mine),
                "only_deposited": len(only_dep),
                "generation_differs": len(gen_diff),
                "agree": len(only_mine) == 0 and len(only_dep) == 0 and len(gen_diff) == 0,
            }
        )
        print(
            f"{network:12s} regenerated {len(mine):5,d} | deposited {len(deposited):5,d} | "
            f"only-new {len(only_mine):4,d} | only-dep {len(only_dep):4,d} | "
            f"gen differs {len(gen_diff):4,d}",
            flush=True,
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("wrote", OUT.relative_to(REPO))


if __name__ == "__main__":
    main()
