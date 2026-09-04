"""Reaction free energies for whole networks, to their maximum generation.

    uv run --extra thermo python scripts/compute_energies.py [network ...]

Rels files are cumulative, so one pass at the deepest generation covers every
shallower one; slice by generation afterwards rather than recomputing.

Resumable: reactions already present in the output CSV are skipped, so the run
can be interrupted and restarted without losing work.
"""

import ast
import csv
import sys
import time
from pathlib import Path

import pandas as pd

from nucleoside_analogues.rels import pivot_rels
from nucleoside_analogues.thermo import compound_cache, reaction_energies

REPO = Path(__file__).resolve().parent.parent
RELS = REPO / "OriginalData" / "OriginalNetworkData" / "Rels"
OUT = REPO / "ProcessedData" / "SI" / "full"
NETWORKS = ("Formose", "FormoseAmm", "Glucose", "GlucoseAmm", "PyruvicAcid")
FIELDS = ["Index", "dG_prime_kJ_mol", "sigma_kJ_mol", "estimable", "status"]
CHUNK = 2000
INFINITE_VARIANCE = 1e4


def literal(value):
    return tuple(ast.literal_eval(value)) if isinstance(value, str) else tuple(value)


def max_generation(network: str) -> int:
    found = [int(f.stem.split("_")[-1]) for f in (RELS / network).glob("*Rels_*.tsv")]
    if not found:
        raise SystemExit(f"no rels files for {network} under {RELS / network}")
    return max(found)


def run(network: str, cc) -> None:
    generation = max_generation(network)
    frame = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
    rows = [
        (str(i), literal(a), literal(b))
        for i, a, b in zip(frame["Index"], frame["Reagents"], frame["Products"], strict=True)
    ]
    target = OUT / f"{network}_G{generation}_energies_pH7.4.csv"
    done: set[str] = set()
    if target.exists():
        with target.open() as handle:
            done = {r["Index"] for r in csv.DictReader(handle)}
    todo = [r for r in rows if r[0] not in done]
    print(f"[{network}] G{generation}: {len(rows):,} reactions, {len(done):,} done, "
          f"{len(todo):,} to do", flush=True)
    if not todo:
        return

    species = list(dict.fromkeys(s for _, a, b in rows for s in (*a, *b)))
    start = time.perf_counter()
    print(f"[{network}] resolving {len(species):,} species", flush=True)
    compounds = compound_cache(species, cc=cc)
    print(f"[{network}] resolved in {time.perf_counter() - start:.0f}s "
          f"({len(compounds[1]):,} unresolved)", flush=True)

    target.parent.mkdir(parents=True, exist_ok=True)
    new = not target.exists()
    with target.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if new:
            writer.writeheader()
        for offset in range(0, len(todo), CHUNK):
            batch = todo[offset : offset + CHUNK]
            for energy in reaction_energies(batch, cc=cc, compounds=compounds):
                estimable = (
                    energy.dg_prime is not None
                    and energy.uncertainty is not None
                    and energy.uncertainty < INFINITE_VARIANCE
                )
                writer.writerow({
                    "Index": energy.index,
                    "dG_prime_kJ_mol": "" if energy.dg_prime is None else energy.dg_prime,
                    "sigma_kJ_mol": "" if energy.uncertainty is None else energy.uncertainty,
                    "estimable": estimable,
                    "status": energy.status,
                })
            handle.flush()
            elapsed = time.perf_counter() - start
            seen = offset + len(batch)
            print(f"[{network}] {seen:,}/{len(todo):,} reactions, {elapsed / 60:.1f} min, "
                  f"eta {elapsed / seen * (len(todo) - seen) / 60:.0f} min", flush=True)
    print(f"[{network}] wrote {target.relative_to(REPO)}", flush=True)


def main() -> None:
    from equilibrator_api import Q_, ComponentContribution

    cc = ComponentContribution()
    cc.p_h = Q_(7.4)
    print(f"conditions: pH {cc.p_h}, I {cc.ionic_strength}, pMg {cc.p_mg}", flush=True)
    for network in sys.argv[1:] or NETWORKS:
        run(network, cc)


if __name__ == "__main__":
    main()
