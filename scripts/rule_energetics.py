"""Reaction free energy by rule and by pH, pooled over the five CRNRs.

    uv run python scripts/rule_energetics.py

The rule is the unit here, not the network: a rule is a general transformation,
so pooling gives its energetic signature and how far that signature moves with
pH. Reactions the estimator cannot resolve are excluded, as everywhere.

`distinct_dG` counts how many different free energies a rule's reactions take.
It is small for rules whose transformation is the same wherever it applies:
group contribution assigns the same increments regardless of the rest of the
molecule, so those reactions are indistinguishable to the estimator however
different their substrates.

Writes ProcessedData/SI/rule_energetics.csv.
"""

import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import FULL, PH_VALUES, PRODUCTS, RELS, REPO, SI, deepest, is_null

from nucleoside_analogues.rels import pivot_rels


def main() -> None:
    values: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for network in PRODUCTS:
        generation = deepest(network)
        rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
        rels["Index"] = rels["Index"].astype(str)
        rule = dict(zip(rels["Index"], rels["Rule"], strict=True))
        for ph in PH_VALUES:
            path = FULL / f"{network}_G{generation}_energies_pH{ph}.csv"
            if not path.exists():
                continue
            for row in csv.DictReader(path.open()):
                if row["estimable"] == "True" and not is_null(row):
                    values[rule[row["Index"]]][ph].append(float(row["dG_prime_kJ_mol"]))
        print(f"  {network} read", flush=True)

    rows = []
    for name in sorted(values, key=lambda r: -len(values[r].get("7.4", []))):
        at = values[name]
        if not at.get("7.4"):
            continue
        base = statistics.median(at["7.4"])
        low, high = at.get("7.0"), at.get("11.0")
        rows.append(
            {
                "rule": name,
                "reactions": len(at["7.4"]),
                "median_dG_pH7.4": round(base, 2),
                "q1_dG_pH7.4": round(statistics.quantiles(at["7.4"], n=4)[0], 2)
                if len(at["7.4"]) > 3
                else "",
                "q3_dG_pH7.4": round(statistics.quantiles(at["7.4"], n=4)[2], 2)
                if len(at["7.4"]) > 3
                else "",
                "distinct_dG": len({round(v, 3) for v in at["7.4"]}),
                "share_exergonic": round(sum(1 for v in at["7.4"] if v < 0) / len(at["7.4"]), 3),
                "median_shift_pH7_to_11": round(statistics.median(high) - statistics.median(low), 2)
                if low and high
                else "",
            }
        )
    out = SI / "rule_energetics.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n  {len(rows)} rules, {sum(r['reactions'] for r in rows):,} reactions at pH 7.4")
    flat = sorted(
        (r for r in rows if r["reactions"] >= 1000),
        key=lambda r: r["distinct_dG"] / r["reactions"],
    )[:4]
    print("  rules the estimator barely distinguishes within:")
    for r in flat:
        print(
            f"    {r['rule'][:40]:42s} {r['reactions']:7,d} reactions, "
            f"{r['distinct_dG']:5,d} distinct values"
        )
    shifted = [r for r in rows if r["median_shift_pH7_to_11"] not in ("", 0.0)]
    big = sorted(shifted, key=lambda r: -abs(r["median_shift_pH7_to_11"]))[:5]
    print("  largest pH sensitivity:")
    for r in big:
        print(f"    {r['rule'][:44]:46s} {r['median_shift_pH7_to_11']:+7.1f} kJ/mol")
    print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
