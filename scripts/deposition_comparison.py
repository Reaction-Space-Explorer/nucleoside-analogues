"""How the recomputed free energies differ from the earlier deposition, at G3.

    uv run python scripts/deposition_comparison.py

All free energies here were recomputed rather than carried over, so the two
sets can disagree. Where a reaction the earlier deposition called spontaneous
is no longer called so, the question is whether the estimate changed sign or
whether component contribution cannot evaluate the reaction at all; the latter
dominates, and is reported rather than hidden, because a reaction dropped for
want of an estimate is a different statement from one dropped on its energy.

The earlier deposition used the sign of the point estimate; the criterion here
is that the 95% interval lie below zero, which is stricter, so some of the
difference is the criterion rather than the numbers.

Writes ProcessedData/SI/deposition_comparison.csv.
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import PRODUCTS, REPO, SI, is_null

OLD = REPO / "ProcessedData" / "SpontaneousRelsWithThermoFiles"


def main() -> None:
    rows = []
    for network in PRODUCTS:
        with (OLD / f"Spontaneous{network}G3RelsWithThermo.tsv").open() as handle:
            old = {r["Index"] for r in csv.DictReader(handle, delimiter="\t")}
        with (SI / f"{network}_G3_energies_pH7.4.csv").open() as handle:
            new = {r["Index"]: r for r in csv.DictReader(handle)}

        lost = unestimable = missing = 0
        for index in old:
            row = new.get(index)
            if row is None:
                missing += 1
                continue
            estimable = row["estimable"] == "True" and not is_null(row)
            spontaneous = (
                estimable and float(row["dG_prime_kJ_mol"]) + 1.96 * float(row["sigma_kJ_mol"]) < 0
            )
            if not spontaneous:
                lost += 1
                if not estimable:
                    unestimable += 1
        rows.append(
            {
                "network": network,
                "spontaneous_in_deposition": len(old),
                "not_spontaneous_now": lost,
                "of_which_unestimable": unestimable,
                "share_unestimable": round(100 * unestimable / lost, 1) if lost else "",
                "absent_from_recomputation": missing,
            }
        )
        print(
            f"  {network:12s} deposition {len(old):6,d}  lost {lost:5,d}  "
            f"unestimable {unestimable:5,d} ({rows[-1]['share_unestimable']}%)"
            + (f"  absent {missing}" if missing else ""),
            flush=True,
        )

    out = SI / "deposition_comparison.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    shares = [r["share_unestimable"] for r in rows if r["share_unestimable"] != ""]
    print(f"\n  share unestimable ranges {min(shares):.0f}% to {max(shares):.0f}% across the CRNRs")
    print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
