"""Does the reachability ordering of the routes survive as a kinetic ordering?

    uv run python scripts/kinetic_ordering.py

The paper orders routes by how many spontaneous reactions they need. That is a
count, not a rate, and the obvious objection is that a short route built from
sluggish steps could be slower than a long route built from fast ones.

Rates cannot be computed here: barriers are not available. But they need not be.
Rule-based stochastic treatments of this chemistry assign barriers by a
Bell-Evans-Polanyi relation, Ea = alpha * dGr'° + beta, with alpha and beta
global rather than fitted per rule. Two consequences follow without ever
choosing values for them:

  1. Within a route, Ea is monotone increasing in dGr'°, so the highest barrier
     belongs to the least exergonic step. The rate-limiting step is therefore
     the least committed step, which is already computed.
  2. Across routes, ordering by rate-limiting barrier is ordering by that step's
     dGr'°. alpha and beta are global and the map is monotone, so they cancel:
     the ordering is the same for every parameterisation of the relation.

Consequence 1 is used in the paper: it gives the least committed step a kinetic
reading at no cost.

Consequence 2 is NOT used, and this script is the reason. Ordering by
rate-limiting barrier does correlate with route length, but the correlation is
mostly an artefact: the least exergonic of n steps drifts upward with n, so a
longer route has a weaker weakest step by construction. Against a null that
draws each route's steps at random from the pooled distribution and keeps the
lengths, the observed correlation is not significant on the step-firing basis
the paper uses (rho 0.748, p 0.24) and only marginally so counting distinct
reactions (rho 0.841, p 0.04). A result that changes verdict with a defensible
bookkeeping choice is not a result, so the paper claims no kinetic ordering.

Writes ProcessedData/SI/kinetic_ordering.csv, both bases, as a deposited
negative.
"""

import csv
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import REPO, SI

ROUTES = REPO / "figures" / "routes"
DRAWS = 20000
SEED = 0


def steps(node, seen: set[str] | None = None) -> list[float]:
    """Step free energies; with `seen`, each reaction counted once per route."""
    out = []
    reaction = node.get("reaction")
    if reaction and "dg" in reaction and (seen is None or reaction["id"] not in seen):
        if seen is not None:
            seen.add(reaction["id"])
        out.append(float(reaction["dg"]))
    for child in node.get("from", []):
        out += steps(child, seen)
    return out


def analyse(distinct: bool) -> dict:
    routes = {}
    for spec in sorted(ROUTES.glob("*.yaml")):
        found = steps(yaml.safe_load(spec.read_text())["target"], set() if distinct else None)
        if found:
            routes[spec.stem] = found

    lengths = np.array([len(v) for v in routes.values()])
    # Ea is monotone increasing in dGr'°, so the largest dGr'° is the highest barrier
    limiting = np.array([max(v) for v in routes.values()])
    rho, p_rho = stats.spearmanr(lengths, limiting)

    pool = np.array([x for v in routes.values() for x in v])
    rng = np.random.default_rng(SEED)
    null = np.empty(DRAWS)
    for i in range(DRAWS):
        null[i] = stats.spearmanr(
            lengths, [pool[rng.integers(0, len(pool), k)].max() for k in lengths]
        )[0]
    p_null = (np.sum(null >= rho) + 1) / (DRAWS + 1)

    return {
        "basis": "distinct_reactions" if distinct else "step_firings",
        "routes": len(routes),
        "steps_total": int(lengths.sum()),
        "spearman_rho": round(float(rho), 3),
        "spearman_p": f"{p_rho:.2e}",
        "null_mean_rho": round(float(null.mean()), 3),
        "null_95th_percentile": round(float(np.percentile(null, 95)), 3),
        "p_against_length_artefact": round(float(p_null), 4),
        "survives_null": "yes" if p_null < 0.05 else "no",
    }


def main() -> None:
    rows = [analyse(False), analyse(True)]
    with (SI / "kinetic_ordering.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(
            f"  {row['basis']:20s} n={row['steps_total']:4d}  rho {row['spearman_rho']:+.3f}  "
            f"null mean {row['null_mean_rho']:+.3f}  p {row['p_against_length_artefact']:.3f}  "
            f"survives null: {row['survives_null']}"
        )
    print("\n  the verdict changes with the bookkeeping, so no kinetic ordering is claimed")
    print(f"wrote {(SI / 'kinetic_ordering.csv').relative_to(REPO)}")


if __name__ == "__main__":
    main()
