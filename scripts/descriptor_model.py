"""Do molecular descriptors predict which matched analogues are reachable?

    uv run --extra ml python scripts/descriptor_model.py

The complexity analysis found that one descriptor, Bottcher complexity, does
not separate the analogues a CRNR reaches spontaneously from those it does not.
This asks the same of a whole panel, and guards against the obvious confound:
descriptors track molecular size, size tracks the generation a species first
appears in, and generation tracks reachability. Generation alone already gives
an area of 0.54 to 0.61, so a descriptor model has to beat that to be saying
anything about chemistry rather than about depth.

Models are compared per CRNR under the same folds: descriptors alone,
generation alone, and both. PyruvicAcid is excluded, having one reachable
matched species in 4,971 and so no second class to predict.

The result does not mean what it first appears to. Reachability on the
estimable-only basis is entangled with what component contribution can
evaluate at all, since a species nothing estimable produces is unreachable by
construction. So the same descriptors are also asked to predict estimator
coverage. They do it better than they predict reachability, which is the
finding: the panel is largely reading which structures fall inside the
group-contribution basis, not which are prebiotically accessible.

Writes ProcessedData/SI/descriptor_model.csv and, with --shap, the mean
absolute SHAP value per descriptor.
"""

import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_si_tables import (
    PRODUCTS,
    RELS,
    REPO,
    SI,
    admitted,
    deepest,
    is_null,
    species_generations,
)

from nucleoside_analogues.descriptors import calc_descriptors
from nucleoside_analogues.hyperpath import shortest_pathways
from nucleoside_analogues.rels import build_index, pivot_rels, read_products

SEED = 0
FOLDS = 5
EXCLUDE = {"PyruvicAcid"}


def dataset(network: str) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    generation = deepest(network)
    rels = pivot_rels(pd.read_csv(RELS / network / f"{network}Rels_{generation}.tsv", sep="\t"))
    rels["Index"] = rels["Index"].astype(str)
    products = read_products(
        REPO / "OriginalData" / "OriginalNetworkData" / "Products" / PRODUCTS[network]
    )
    seeds = tuple(products.loc[products["Generation"] == 0, "Smiles"])
    cost = shortest_pathways(
        build_index(
            rels[rels["Index"].isin(admitted(network, rels, "estimable_only", generation))]
        ),
        seeds,
    ).cost
    energies = {
        r["Index"]: (r["estimable"] == "True" and not is_null(r))
        for r in csv.DictReader(
            (SI / "full" / f"{network}_G{generation}_energies_pH7.4.csv").open()
        )
    }
    made_by_estimable: set[str] = set()
    for i, made in zip(rels["Index"], rels["Products"], strict=True):
        if energies.get(i, False):
            made_by_estimable.update(made)
    matched = pd.read_csv(
        REPO / "ProcessedData" / "MatchesFiles" / f"{network}Matches.tsv", sep="\t"
    )["NetworkSmiles"].astype(str)
    gens = species_generations(network)
    smiles = [s for s in dict.fromkeys(matched) if s in gens]
    frame = calc_descriptors(smiles)
    keep = [c for c in frame.columns if frame[c].nunique(dropna=True) > 1]
    dropped = sorted(set(frame.columns) - set(keep))
    if dropped:
        print(f"    dropped {len(dropped)} constant descriptors: {', '.join(dropped)}", flush=True)
    return (
        frame[keep],
        np.array([int(s in cost) for s in smiles]),
        np.array([gens[s] for s in smiles], dtype=float),
        np.array([int(s in made_by_estimable) for s in smiles]),
    )


def cross_val_auc(x: np.ndarray, y: np.ndarray) -> float:
    folds = StratifiedKFold(n_splits=FOLDS, shuffle=True, random_state=SEED)
    scores = []
    for train, test in folds.split(x, y):
        model = XGBClassifier(
            n_estimators=200,
            max_depth=4,
            learning_rate=0.1,
            random_state=SEED,
            n_jobs=4,
            eval_metric="logloss",
        )
        model.fit(x[train], y[train])
        scores.append(roc_auc_score(y[test], model.predict_proba(x[test])[:, 1]))
    return float(np.mean(scores))


def main() -> None:
    want_shap = "--shap" in sys.argv
    rows, shap_rows = [], []
    for network in PRODUCTS:
        if network in EXCLUDE:
            print(f"  {network:12s} excluded: one reachable species, no second class", flush=True)
            continue
        print(f"  {network}", flush=True)
        frame, y, gen, estimable = dataset(network)
        d = frame.to_numpy(dtype=float)
        both = np.column_stack([d, gen])
        result = {
            "network": network,
            "matched": len(y),
            "reachable": int(y.sum()),
            "descriptors": frame.shape[1],
            "auc_descriptors": round(cross_val_auc(d, y), 3),
            "auc_generation": round(cross_val_auc(gen.reshape(-1, 1), y), 3),
            "auc_both": round(cross_val_auc(both, y), 3),
            "auc_estimator_coverage": round(cross_val_auc(d, estimable), 3),
        }
        rows.append(result)
        print(
            f"    n={result['matched']:,}  reachable {100 * result['reachable'] / result['matched']:.0f}%"
            f"   descriptors {result['auc_descriptors']:.3f}"
            f"   generation {result['auc_generation']:.3f}"
            f"   both {result['auc_both']:.3f}"
            f"   estimator coverage {result['auc_estimator_coverage']:.3f}",
            flush=True,
        )
        if want_shap:
            import shap

            model = XGBClassifier(
                n_estimators=200,
                max_depth=4,
                learning_rate=0.1,
                random_state=SEED,
                n_jobs=4,
                eval_metric="logloss",
            ).fit(d, y)
            sample = d[np.random.default_rng(SEED).choice(len(d), min(2000, len(d)), replace=False)]
            # additivity is checked, not disabled: attributions that do not sum
            # to the prediction are not attributions
            values = shap.TreeExplainer(model).shap_values(sample, check_additivity=True)
            for name, mean_abs in zip(frame.columns, np.abs(values).mean(axis=0), strict=True):
                shap_rows.append(
                    {"network": network, "descriptor": name, "mean_abs_shap": round(mean_abs, 5)}
                )

    for table, name in ((rows, "descriptor_model.csv"), (shap_rows, "descriptor_shap.csv")):
        if not table:
            continue
        out = SI / name
        with out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(table[0]))
            writer.writeheader()
            writer.writerows(table)
        print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
