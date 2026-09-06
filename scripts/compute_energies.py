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
#: pH values reported. 7.4 is the reference; the others show how far the
#: spontaneous set moves with proton activity. Compound resolution dominates
#: the cost and does not depend on pH, so extra values are nearly free.
PH_VALUES = (7.0, 7.4, 9.0, 11.0)
_CC = None


def _init() -> None:
    global _CC
    from equilibrator_api import ComponentContribution

    _CC = ComponentContribution()


def _score(block: list[tuple[str, tuple[str, ...], tuple[str, ...]]]) -> dict[float, list[dict]]:
    """Resolve this block once, then score it at every pH. Returns pH -> rows."""
    from equilibrator_api import Q_, Reaction

    from nucleoside_analogues.thermo import compound_cache

    species = list(dict.fromkeys(s for _, a, b in block for s in (*a, *b)))
    resolved, missing = compound_cache(species, cc=_CC)
    unresolved = set(missing)

    blocked_rows: list[dict] = []
    scorable, reactions = [], []
    for identifier, reagents, products in block:
        blocked = [s for s in (*reagents, *products) if s in unresolved]
        if blocked:
            blocked_rows.append(
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

    out: dict[float, list[dict]] = {}
    for p_h in PH_VALUES:
        rows = list(blocked_rows)
        if reactions:
            _CC.p_h = Q_(p_h)
            try:
                values, sqrt_factor = _CC.standard_dg_prime_multi(
                    reactions, uncertainty_representation="sqrt")
                dg = np.asarray(values.m_as("kJ/mol"), dtype=float)
                sigma = np.sqrt(
                    (np.asarray(sqrt_factor.m_as("kJ/mol"), dtype=float) ** 2).sum(axis=1))
                for identifier, value, error in zip(scorable, dg, sigma, strict=True):
                    rows.append({"Index": identifier, "dG_prime_kJ_mol": float(value),
                                 "sigma_kJ_mol": float(error),
                                 "estimable": bool(error < INFINITE_VARIANCE), "status": "ok"})
            except Exception as error:  # noqa: BLE001 - reported, never swallowed
                for identifier in scorable:
                    rows.append({"Index": identifier, "dG_prime_kJ_mol": "", "sigma_kJ_mol": "",
                                 "estimable": False,
                                 "status": f"estimation_failed: {str(error)[:80]}"})
        out[p_h] = rows
    return out


def literal(value):
    return tuple(ast.literal_eval(value)) if isinstance(value, str) else tuple(value)


def run(network: str, workers: int, generation: int | None = None, out: Path = OUT) -> None:
    if generation is None:
        generation = max(int(f.stem.split("_")[-1]) for f in (RELS / network).glob("*Rels_*.tsv"))
    frame = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
    rows = [
        (str(i), literal(a), literal(b))
        for i, a, b in zip(frame["Index"], frame["Reagents"], frame["Products"], strict=True)
    ]

    out.mkdir(parents=True, exist_ok=True)
    targets = {p: out / f"{network}_G{generation}_energies_pH{p}.csv" for p in PH_VALUES}
    reference = targets[7.4]
    done: set[str] = set()
    if reference.exists():
        with reference.open() as handle:
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
    handles, writers = {}, {}
    for p, path in targets.items():
        fresh = not path.exists()
        handles[p] = path.open("a", newline="")
        writers[p] = csv.DictWriter(handles[p], fieldnames=FIELDS)
        if fresh:
            writers[p].writeheader()
    try:
        with ProcessPoolExecutor(workers, initializer=_init) as pool:
            for n, per_ph in enumerate(pool.map(_score, blocks), 1):
                for p, block_rows in per_ph.items():
                    writers[p].writerows(block_rows)
                    handles[p].flush()
                elapsed = (time.perf_counter() - start) / 60
                print(f"[{network}] block {n}/{workers} written, {elapsed:.1f} min", flush=True)
    finally:
        for h in handles.values():
            h.close()
    print(
        f"[{network}] wrote {len(targets)} pH files in "
        f"{(time.perf_counter() - start) / 60:.1f} min",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--generation",
        type=int,
        default=None,
        help="fix the generation instead of using the deepest available",
    )
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("networks", nargs="*", default=list(NETWORKS))
    args = parser.parse_args()
    for network in args.networks or NETWORKS:
        run(network, args.workers, args.generation, args.out)


if __name__ == "__main__":
    main()
