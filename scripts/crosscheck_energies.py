"""Does another machine reproduce the deposited free energies exactly?

    uv run --extra thermo python scripts/crosscheck_energies.py

Component contribution is deterministic, but it depends on a downloaded
compound cache and on the local eQuilibrator build, so agreement across
machines is worth demonstrating rather than assuming. This recomputes the
generation-three energies and compares them against the deposited file,
reporting the largest difference and any reaction whose estimability changed.

Run on a 64-core Linux host against a set computed on macOS, the two agreed to
2.3e-13 kJ/mol over 653 estimable reactions with no estimability change; that
transcript is in ProcessedData/SI/logs/xcheck.txt.
"""

import ast
import csv
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import REPO, SI

from nucleoside_analogues.thermo import reaction_energies

NETWORKS = ("Formose", "PyruvicAcid")
TOLERANCE = 1e-9


def literal(value):
    return ast.literal_eval(value) if isinstance(value, str) else tuple(value)


def main() -> None:
    from equilibrator_api import Q_, ComponentContribution

    cc = ComponentContribution()
    cc.p_h = Q_(7.4)
    print(f"pH {cc.p_h}, I {cc.ionic_strength}, pMg {cc.p_mg}, T {cc.temperature}", flush=True)

    failures = 0
    for network in NETWORKS:
        rels = pd.read_csv(
            REPO / "ProcessedData" / "RelsFiles" / network / f"{network}G3ProcessedRels.tsv",
            sep="\t",
        )
        rows = [
            (str(i), literal(a), literal(b))
            for i, a, b in zip(rels["Index"], rels["Reagents"], rels["Products"], strict=True)
        ]
        recomputed = {e.index: e for e in reaction_energies(rows, cc=cc)}
        with (SI / f"{network}_G3_energies_pH7.4.csv").open() as handle:
            deposited = {r["Index"]: r for r in csv.DictReader(handle)}

        worst = 0.0
        changed = compared = differing = 0
        for index, old in deposited.items():
            new = recomputed.get(index)
            if new is None:
                changed += 1
                continue
            was = old["estimable"] == "True"
            now = new.dg_prime is not None and new.uncertainty is not None and new.uncertainty < 1e4
            if was != now:
                changed += 1
            elif was:
                compared += 1
                gap = abs(float(old["dG_prime_kJ_mol"]) - new.dg_prime)
                worst = max(worst, gap)
                differing += gap > TOLERANCE
        ok = changed == 0 and worst < TOLERANCE
        failures += not ok
        print(
            f"{network:12s} {len(deposited):5d} rows, {compared:5d} estimable compared | "
            f"estimability changed {changed} | differing {differing} | "
            f"max diff {worst:.3e} kJ/mol | {'IDENTICAL' if ok else 'DIFFERS'}",
            flush=True,
        )
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
