"""Reaction free energies for whole networks, to their deepest generation.

    uv run --extra thermo python scripts/compute_energies.py [--workers N] [network ...]

Rels files are cumulative, so one pass at the deepest generation covers every
shallower one; slice by generation afterwards rather than recomputing.

Two things make this tractable. Energies come from ``standard_dg_prime_multi``,
which is about twelve times faster than scoring reactions one at a time and
agrees with it to 7e-15 kJ/mol on every estimable reaction. Compound
resolution, which dominates, is spread over worker processes: each takes a
disjoint block of reactions and resolves only the species it needs, returning
plain numbers rather than ORM objects.

Resumable: reactions already in the output are skipped, so an interrupted run
restarts where it stopped.
"""

import argparse
import ast
import csv
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from nucleoside_analogues.rels import pivot_rels

REPO = Path(__file__).resolve().parent.parent
RELS = REPO / "OriginalData" / "OriginalNetworkData" / "Rels"
OUT = REPO / "ProcessedData" / "SI" / "full"
NETWORKS = ("Formose", "FormoseAmm", "Glucose", "GlucoseAmm", "PyruvicAcid")
FIELDS = ["Index", "dG_prime_kJ_mol", "sigma_kJ_mol", "estimable", "status"]
#: Component contribution flags an unusable estimate with a huge variance.
INFINITE_VARIANCE = 1e4
_CC = None


def _init() -> None:
    global _CC
    from equilibrator_api import Q_, ComponentContribution

    _CC = ComponentContribution()
    _CC.p_h = Q_(7.4)


def _score(block: list[tuple[str, tuple[str, ...], tuple[str, ...]]]) -> list[dict]:
    """Resolve this block's species and score its reactions. Returns plain rows."""
    from equilibrator_api import Reaction

    from nucleoside_analogues.thermo import compound_cache

    species = list(dict.fromkeys(s for _, a, b in block for s in (*a, *b)))
    resolved, missing = compound_cache(species, cc=_CC)
    unresolved = set(missing)

    rows: list[dict] = []
    scorable, reactions = [], []
    for identifier, reagents, products in block:
        blocked = [s for s in (*reagents, *products) if s in unresolved]
        if blocked:
            rows.append(
                {
                    "Index": identifier,
                    "dG_prime_kJ_mol": "",
                    "sigma_kJ_mol": "",
                    "estimable": False,
                    "status": "compound_missing",
                }
            )
            continue
        stoichiometry: dict = {}
        for s in reagents:
            stoichiometry[resolved[s]] = stoichiometry.get(resolved[s], 0.0) - 1.0
        for s in products:
            stoichiometry[resolved[s]] = stoichiometry.get(resolved[s], 0.0) + 1.0
        scorable.append(identifier)
        reactions.append(Reaction(stoichiometry))

    if reactions:
        try:
            values, sqrt_factor = _CC.standard_dg_prime_multi(
                reactions, uncertainty_representation="sqrt"
            )
            dg = np.asarray(values.m_as("kJ/mol"), dtype=float)
            sigma = np.sqrt((np.asarray(sqrt_factor.m_as("kJ/mol"), dtype=float) ** 2).sum(axis=1))
            for identifier, value, error in zip(scorable, dg, sigma, strict=True):
                rows.append(
                    {
                        "Index": identifier,
                        "dG_prime_kJ_mol": float(value),
                        "sigma_kJ_mol": float(error),
                        "estimable": bool(error < INFINITE_VARIANCE),
                        "status": "ok",
                    }
                )
        except Exception as error:  # noqa: BLE001 - reported, never swallowed
            for identifier in scorable:
                rows.append(
                    {
                        "Index": identifier,
                        "dG_prime_kJ_mol": "",
                        "sigma_kJ_mol": "",
                        "estimable": False,
                        "status": f"estimation_failed: {str(error)[:80]}",
                    }
                )
    return rows


def literal(value):
    return tuple(ast.literal_eval(value)) if isinstance(value, str) else tuple(value)


def run(network: str, workers: int) -> None:
    generation = max(int(f.stem.split("_")[-1]) for f in (RELS / network).glob("*Rels_*.tsv"))
    frame = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
    rows = [
        (str(i), literal(a), literal(b))
        for i, a, b in zip(frame["Index"], frame["Reagents"], frame["Products"], strict=True)
    ]

    target = OUT / f"{network}_G{generation}_energies_pH7.4.csv"
    target.parent.mkdir(parents=True, exist_ok=True)
    done: set[str] = set()
    if target.exists():
        with target.open() as handle:
            done = {r["Index"] for r in csv.DictReader(handle)}
    todo = [r for r in rows if r[0] not in done]
    print(
        f"[{network}] G{generation}: {len(rows):,} reactions, {len(done):,} done, "
        f"{len(todo):,} to do",
        flush=True,
    )
    if not todo:
        return

    blocks = [todo[i::workers] for i in range(workers)]
    start = time.perf_counter()
    new = not target.exists()
    with target.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if new:
            writer.writeheader()
        with ProcessPoolExecutor(workers, initializer=_init) as pool:
            for n, block_rows in enumerate(pool.map(_score, blocks), 1):
                writer.writerows(block_rows)
                handle.flush()
                elapsed = (time.perf_counter() - start) / 60
                print(f"[{network}] block {n}/{workers} written, {elapsed:.1f} min", flush=True)
    print(
        f"[{network}] wrote {target.relative_to(REPO)} "
        f"in {(time.perf_counter() - start) / 60:.1f} min",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("networks", nargs="*", default=list(NETWORKS))
    args = parser.parse_args()
    for network in args.networks or NETWORKS:
        run(network, args.workers)


if __name__ == "__main__":
    main()
