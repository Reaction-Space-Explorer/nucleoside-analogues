"""Match network products against curated metabolite databases, by generation.

    uv run python scripts/database_matches.py

ChEBI is downloaded on first use. HMDB must be placed by hand at
OriginalData/reference_databases/hmdb_structures.sdf, since it refuses scripted
download; it is skipped with a warning when absent, rather than silently
dropped from the output.

Matching is on the first 14 InChIKey characters, as everywhere else here, so a
hit means the database contains some stereoisomer of the same constitution.
"""

import csv
import gzip
import sys
import urllib.request
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import PRODUCTS, SI

RDLogger.DisableLog("rdApp.*")
REPO = Path(__file__).resolve().parent.parent
REF = REPO / "OriginalData" / "reference_databases"
CHEBI_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi_lite_3_stars.sdf.gz"


def skeleton_keys_from_sdf(path: Path) -> set[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    keys: set[str] = set()
    with opener(path, "rb") as handle:
        for mol in Chem.ForwardSDMolSupplier(handle):
            if mol is None or mol.GetNumAtoms() == 0:
                continue
            try:
                keys.add(Chem.MolToInchiKey(mol)[:14])
            except Exception:  # noqa: BLE001 - a bad record must not abort the file
                continue
    return keys


def main() -> None:
    REF.mkdir(parents=True, exist_ok=True)
    chebi = REF / "chebi_lite_3_stars.sdf.gz"
    if not chebi.exists():
        print(f"downloading ChEBI to {chebi.name}", flush=True)
        urllib.request.urlretrieve(CHEBI_URL, chebi)  # noqa: S310 - fixed EBI URL

    databases = {"ChEBI": chebi}
    hmdb = REF / "hmdb_structures.sdf"
    if hmdb.exists():
        databases["HMDB"] = hmdb
    else:
        print("HMDB not present; download 'Structures' from https://www.hmdb.ca/downloads "
              f"and unzip to {hmdb}. Continuing without it.", flush=True)

    reference = {}
    for name, path in databases.items():
        reference[name] = skeleton_keys_from_sdf(path)
        print(f"  {name}: {len(reference[name]):,} distinct skeletons", flush=True)

    rows = []
    for network, products_file in PRODUCTS.items():
        products = pd.read_csv(
            REPO / "OriginalData" / "OriginalNetworkData" / "Products" / products_file, sep="\t"
        )
        seen: dict[int, set[str]] = {}
        for smiles, generation in zip(products["Smiles"], products["Generation"], strict=True):
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is None:
                continue
            g = int(str(generation).lstrip("Gg"))
            seen.setdefault(g, set()).add(Chem.MolToInchiKey(mol)[:14])
        cumulative: set[str] = set()
        for g in sorted(seen):
            cumulative |= seen[g]
            row = {"network": network, "generation": g, "products_cumulative": len(cumulative)}
            for name, keys in reference.items():
                row[name] = len(cumulative & keys)
            rows.append(row)
            print(f"  {network:12s} G{g}  products {len(cumulative):7,d}  "
                  + "  ".join(f"{n} {row[n]:4,d}" for n in reference), flush=True)

    with (SI / "database_matches.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("wrote ProcessedData/SI/database_matches.csv")


if __name__ == "__main__":
    main()
