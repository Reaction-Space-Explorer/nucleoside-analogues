"""How committed each traced route is, step by step.

    uv run python scripts/route_commitment.py

Lauber et al. (ALIFE 2025, 10.1162/isal.a.847) assign kinetic barriers in a
way consistent with detailed balance: an exergonic step carries an intrinsic
barrier set by its reaction rule, and the reverse of that step carries the
same barrier plus the magnitude of the reaction free energy.

Two things follow for the routes traced here. Every step is exergonic, so
forward barriers are fixed by rule identity and no ordering of rates is
possible without barrier parameters this work does not supply. The reverse
direction is bounded, though: the ratio of reverse to forward rate is
exp(dGr/RT), which needs only the free energies already computed.

Writes ProcessedData/SI/route_commitment.csv.
"""

import csv
import math
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import REPO, SI

RT = 8.314462618e-3 * 298.15  # kJ/mol
ROUTES = REPO / "figures" / "routes"


def steps(node) -> list[tuple[str, float]]:
    out = []
    reaction = node.get("reaction")
    if reaction and "dg" in reaction:
        out.append((reaction["id"], float(reaction["dg"])))
    for child in node.get("from", []):
        out += steps(child)
    return out


def main() -> None:
    rows, every = [], []
    for spec in sorted(ROUTES.glob("*.yaml")):
        found = steps(yaml.safe_load(spec.read_text())["target"])
        if not found:
            continue
        network, target = spec.stem.split("_", 1)
        weakest_id, weakest = max(found, key=lambda x: x[1])
        every += [dg for _, dg in found]
        rows.append(
            {
                "network": network,
                "target": target,
                "steps": len(found),
                "median_dG_kJ_mol": round(sorted(dg for _, dg in found)[len(found) // 2], 1),
                "least_committed_step": weakest_id,
                "least_committed_dG_kJ_mol": round(weakest, 1),
                "reverse_over_forward": f"{math.exp(weakest / RT):.2g}",
            }
        )
        print(
            f"  {network:12s} {target:12s} {len(found):2d} steps, "
            f"least committed {weakest:6.1f} kJ/mol "
            f"(reverse/forward {math.exp(weakest / RT):.1e})",
            flush=True,
        )
    committed = sum(1 for dg in every if dg < -25)
    print(
        f"\n  {len(every)} steps in all; median {sorted(every)[len(every) // 2]:.1f} kJ/mol; "
        f"{committed} ({100 * committed / len(every):.0f}%) below -25 kJ/mol, "
        f"reverse flux under {math.exp(-25 / RT):.0e} of forward"
    )
    out = SI / "route_commitment.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
