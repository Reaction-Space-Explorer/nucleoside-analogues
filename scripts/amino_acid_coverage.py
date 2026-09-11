"""What the five CRNRs make in place of nucleobases.

None of these networks builds a nucleobase, which the manuscript reports. The
nitrogen chemistry they do carry ends elsewhere: the two ammonia-seeded CRNRs
produce alpha-amino acids, and the three seeded without nitrogen cannot.

Counts are substructure matches, not analogue-library matches -- no amino acid
library was enumerated for this work -- so they bound what is present rather
than characterising it.

Writes ProcessedData/SI/amino_acid_coverage.csv.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import PRODUCTS, REPO, SI

from nucleoside_analogues.rels import read_products

RDLogger.DisableLog("rdApp.*")

#: Primary amine or secondary amine on the carbon alpha to a carboxylate,
#: excluding amides. Matches glycine and sarcosine, not beta-alanine.
ALPHA = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)][CX4][CX3](=O)[OX2H1,OX1-]")
CANONICAL = {
    "glycine": "NCC(O)=O",
    "alanine": "CC(N)C(O)=O",
    "serine": "OCC(N)C(O)=O",
    "threonine": "CC(O)C(N)C(O)=O",
    "aspartate": "OC(=O)CC(N)C(O)=O",
    "glutamate": "OC(=O)CCC(N)C(O)=O",
    "valine": "CC(C)C(N)C(O)=O",
    "leucine": "CC(C)CC(N)C(O)=O",
    "isoleucine": "CCC(C)C(N)C(O)=O",
    "proline": "OC(=O)C1CCCN1",
    "asparagine": "NC(=O)CC(N)C(O)=O",
    "glutamine": "NC(=O)CCC(N)C(O)=O",
    "cysteine": "SCC(N)C(O)=O",
    "lysine": "NCCCCC(N)C(O)=O",
    "arginine": "NC(N)=NCCCC(N)C(O)=O",
    "histidine": "OC(=O)C(N)Cc1cnc[nH]1",
    "phenylalanine": "OC(=O)C(N)Cc1ccccc1",
    "tyrosine": "OC(=O)C(N)Cc1ccc(O)cc1",
    "methionine": "CSCCC(N)C(O)=O",
    "tryptophan": "OC(=O)C(N)Cc1c[nH]c2ccccc12",
}


def main() -> None:
    keys = {n: Chem.MolToInchiKey(Chem.MolFromSmiles(s))[:14] for n, s in CANONICAL.items()}
    rows = []
    for network, listing in PRODUCTS.items():
        products = read_products(
            REPO / "OriginalData" / "OriginalNetworkData" / "Products" / listing
        )
        seen: dict[str, int] = {}
        alpha = parsed = 0
        for smiles, generation in zip(products["Smiles"], products["Generation"], strict=True):
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is None:
                continue
            parsed += 1
            if mol.HasSubstructMatch(ALPHA):
                alpha += 1
            seen.setdefault(Chem.MolToInchiKey(mol)[:14], int(generation))
        present = sorted(n for n, k in keys.items() if k in seen)
        rows.append(
            {
                "network": network,
                "species": parsed,
                "alpha_amino_acids": alpha,
                "percent": round(100 * alpha / parsed, 2),
                "canonical_present": len(present),
                "which": ";".join(present),
            }
        )
        print(f"  {network:<12} {alpha:>6,} of {parsed:>7,} ({rows[-1]['percent']:>5.2f}%)  {';'.join(present) or '-'}")

    out = SI / "amino_acid_coverage.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
