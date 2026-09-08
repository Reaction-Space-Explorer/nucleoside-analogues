"""The uronate route to pentoses, and what the estimator can see of it.

    uv run python scripts/uronate_case.py

Yi et al. (JACS Au 2023, 3, 2522) make pentoses without the formose reaction:
nonselective oxidation of a C6 aldonate gives a mixture of oxo-uronates,
carbonyl migration carries the carbonyl to C3, and beta-decarboxylation then
yields ribose, arabinose, xylose and lyxose. It is the transformation the
pentose phosphate pathway performs on 6-phosphogluconate.

The networks here contain the species and the decisive step, but not one
instance of that step is admitted: each carries the infinite-variance sentinel,
so each is unestimable. The route is excluded for want of an estimate rather
than for want of a favourable one, which is the same outcome as the hexose
branch point reached by the other mechanism -- a null estimate there, unbounded
uncertainty here.

This also breaks the unestimable count into its causes, which SI Table 3 reports
only as a total.

Writes ProcessedData/SI/uronate_case.csv and estimator_coverage.csv.
"""

import csv
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import (
    PRODUCTS,
    RELS,
    REPO,
    SI,
    admitted,
    deepest,
    energy_file,
    is_null,
    species_generations,
)

from nucleoside_analogues.rels import pivot_rels

RDLogger.DisableLog("rdApp.*")
ALDONATE, URONATE, PENTOSE = "C6H12O7", "C6H10O7", "C5H10O5"
RULE = "Beta Decarboxylation"


def formula(smiles) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    return CalcMolFormula(mol).replace("+", "").replace("-", "") if mol else "?"


def main() -> None:
    route, coverage = [], []
    for network in PRODUCTS:
        generation = deepest(network)
        rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
        rels["Index"] = rels["Index"].astype(str)
        energies = {r["Index"]: r for r in csv.DictReader(energy_file(network, generation).open())}
        keep = set(admitted(network, rels, "estimable_only", generation))

        counts = {ALDONATE: 0, URONATE: 0, PENTOSE: 0}
        for smiles in species_generations(network):
            f = formula(smiles)
            if f in counts:
                counts[f] += 1

        steps = admit = estimable = 0
        worst = ""
        for _, row in rels[rels["Rule"] == RULE].iterrows():
            if URONATE not in {formula(x) for x in row["Reagents"]}:
                continue
            if PENTOSE not in {formula(x) for x in row["Products"]}:
                continue
            steps += 1
            record = energies.get(row["Index"])
            if record and record["estimable"] == "True" and not is_null(record):
                estimable += 1
            if row["Index"] in keep:
                admit += 1
            if record and not worst:
                worst = f"dG {float(record['dG_prime_kJ_mol']):.0f}, sigma {float(record['sigma_kJ_mol']):.0f}"
        route.append(
            {
                "network": network,
                "aldonates": counts[ALDONATE],
                "oxo_uronates": counts[URONATE],
                "pentoses": counts[PENTOSE],
                "uronate_to_pentose_steps": steps,
                "estimable": estimable,
                "admitted": admit,
                "example_estimate": worst,
            }
        )

        total = len(energies)
        null = inf = missing = ok = 0
        for record in energies.values():
            if record["status"] != "ok":
                missing += 1
            elif record["estimable"] != "True":
                inf += 1
            elif is_null(record):
                null += 1
            else:
                ok += 1
        coverage.append(
            {
                "network": network,
                "reactions": total,
                "estimable": ok,
                "null_estimate": null,
                "infinite_variance": inf,
                "compound_missing": missing,
                "percent_unestimable": round(100 * (total - ok) / total, 1),
            }
        )
        print(
            f"  {network:12s} aldonates {counts[ALDONATE]:3d}  uronates {counts[URONATE]:4d}  "
            f"pentoses {counts[PENTOSE]:3d} | uronate->pentose {steps:2d}, estimable {estimable}, "
            f"admitted {admit} | unestimable {coverage[-1]['percent_unestimable']:4.1f}%",
            flush=True,
        )

    for table, name in ((route, "uronate_case.csv"), (coverage, "estimator_coverage.csv")):
        path = SI / name
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(table[0]))
            writer.writeheader()
            writer.writerows(table)
        print(f"wrote {path.relative_to(REPO)}")
    total_steps = sum(r["uronate_to_pentose_steps"] for r in route)
    print(
        f"\n  {total_steps} uronate-to-pentose decarboxylations in all, "
        f"{sum(r['admitted'] for r in route)} admitted"
    )


if __name__ == "__main__":
    main()
