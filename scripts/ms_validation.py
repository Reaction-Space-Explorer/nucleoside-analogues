"""Match FT-ICR MS formulas against the network products.

    uv run python scripts/ms_validation.py

The spectra are negative-ESI, and MIDAS lists the composition of the observed
[M-H]- ion: the measured mass exceeds the listed formula mass by one electron
mass (0.00055 Da), which is how the ion assignment is confirmed here rather
than assumed. The neutral molecule is therefore the listed formula plus one
hydrogen, and that is what is compared against the network.

Only CHO and CHNO formulas are considered; sulfur- and 13C-containing
assignments are counted and excluded, not silently dropped.
"""

import csv
import json
from pathlib import Path

from make_si_tables import species_generations
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")
REPO = Path(__file__).resolve().parent.parent
MS = Path(
    "/private/tmp/claude-501/-Users-sid-Downloads-popvax/c2e5a6c4-f67b-4512-aa86-d6aa1aa364aa/scratchpad/ms/NHMFLMSData.11.16.20"
)
PRODUCTS = REPO / "OriginalData" / "OriginalNetworkData" / "Products"
OUT = REPO / "ProcessedData" / "SI"

#: Sample number -> (network, product listing). From SampleExplanations.docx.
SAMPLES = {
    "39": ("Formose", "formose_output.tsv", "Formose 150 C"),
    "50": ("Formose", "formose_output.tsv", "Formose 85 C"),
    "40": ("FormoseAmm", "formose_amm_output.tsv", "Formose + NH3 150 C"),
    "34": ("Glucose", "glucose_degradation_output.tsv", "Glucose dry 150 C"),
    "38": ("Glucose", "glucose_degradation_output.tsv", "Dextrose wet 150 C"),
    "37": ("GlucoseAmm", "glucose_amm_output.tsv", "Glucose + NH4OH wet"),
    "46": ("PyruvicAcid", "pyruvic_output.tsv", "Pyruvic acid, neat"),
    "47": ("PyruvicAcid", "pyruvic_output.tsv", "Pyruvic acid in water"),
}
ORGANIC = {"C", "H", "N", "O"}
#: The networks only build small molecules, so comparison is restricted to the
#: mass range they can reach. Reported at several cutoffs rather than one.
CUTOFFS = (150.0, 200.0, 250.0, 300.0)


def read_midas(path: Path) -> list[dict]:
    """Peaks with a CHNO assignment, converted from ion to neutral formula."""
    peaks = []
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("Peak Number"):
                break
        for fields in csv.reader(handle):
            if len(fields) < 14:
                continue
            try:
                mass, abundance = float(fields[1]), float(fields[2])
            except ValueError:
                continue
            tail = [x.strip() for x in fields[13:] if x.strip()]
            counts: dict[str, int] = {}
            ok = True
            for i in range(0, len(tail) - 1, 2):
                element, number = tail[i], tail[i + 1]
                if not number.lstrip("-").isdigit():
                    ok = False
                    break
                counts[element] = counts.get(element, 0) + int(number)
            if not ok or not counts:
                continue
            peaks.append(
                {
                    "mass": mass,
                    "abundance": abundance,
                    "ion": counts,
                    "organic": set(counts) <= ORGANIC,
                }
            )
    return peaks


def neutral_formula(ion: dict[str, int]) -> str:
    """Listed ion formula plus one hydrogen, in Hill order."""
    counts = dict(ion)
    counts["H"] = counts.get("H", 0) + 1
    order = ["C", "H"] + sorted(k for k in counts if k not in ("C", "H"))
    return "".join(f"{k}{counts[k]}" for k in order if counts.get(k))


def network_formulas(network: str) -> dict[str, int]:
    """Molecular formula -> earliest generation it appears in."""
    first: dict[str, int] = {}
    for smiles, generation in species_generations(network).items():
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            continue
        formula = CalcMolFormula(mol).replace("+", "").replace("-", "")
        if formula not in first or generation < first[formula]:
            first[formula] = int(generation)
    return first


def main() -> None:
    rows = []
    detail = {}
    for number, (network, _products_file, label) in SAMPLES.items():
        matches = list(MS.glob(f"*_{number}_*"))
        if not matches:
            print(f"sample {number}: FILE MISSING")
            continue
        peaks = read_midas(matches[0])
        organic = [p for p in peaks if p["organic"]]
        net = network_formulas(network)
        row = {
            "sample": number,
            "label": label,
            "network": network,
            "peaks_assigned": len(peaks),
            "peaks_CHNO": len(organic),
        }
        parts = []
        for cutoff in CUTOFFS:
            formulas = {neutral_formula(p["ion"]) for p in organic if p["mass"] < cutoff}
            hit = formulas & set(net)
            pct = round(100 * len(hit) / len(formulas), 1) if formulas else 0.0
            row[f"formulas_lt{int(cutoff)}"] = len(formulas)
            row[f"matched_lt{int(cutoff)}"] = len(hit)
            row[f"percent_lt{int(cutoff)}"] = pct
            parts.append(f"<{int(cutoff)}: {len(hit):3d}/{len(formulas):4d} ({pct:4.1f}%)")
            if cutoff == 250.0:
                detail[number] = {"matched": sorted(hit), "unmatched": sorted(formulas - set(net))}
        rows.append(row)
        print(f"{number:2s} {label:22s} {network:12s} " + "  ".join(parts), flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "ms_validation.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (OUT / "ms_validation_formulas.json").write_text(json.dumps(detail, indent=1))
    print("wrote", (OUT / "ms_validation.csv").relative_to(REPO))


if __name__ == "__main__":
    main()
