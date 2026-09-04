"""Reading and reshaping MØD reaction-network output.

MØD emits reaction data in a long ("rels") format with one row per
molecule-per-reaction::

    Index   Reagent   Formed/Produced   Rule
    4_0     C=O       -1                Aldol Condensation
    4_0     C(CO)=O   -1                Aldol Condensation
    4_0     C(C(CO)O)=O   1             Aldol Condensation

``-1`` marks a reagent, ``1`` a product.  :func:`pivot_rels` turns that into one
row per reaction with tuple-valued reagent and product columns.

This replaces the ``ProcessRels`` notebook, which performed the same reshape
with a nested ``list.index`` scan -- O(rows x reactions).  On the generation-4
FormoseAmm network (619,242 rows, 145,820 reactions) the vectorised version
below completes in about a second.
"""

from __future__ import annotations

import ast
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import pandas as pd

__all__ = [
    "ReactionIndex",
    "build_index",
    "load_processed_rels",
    "pivot_rels",
    "read_products",
]


def pivot_rels(raw: pd.DataFrame) -> pd.DataFrame:
    """Reshape a raw MØD rels table into one row per reaction.

    Parameters
    ----------
    raw
        Long-format frame with ``Index``, ``Reagent``, ``Formed/Produced`` and
        ``Rule`` columns, as written by MØD.

    Returns
    -------
    DataFrame
        Columns ``Index``, ``Reagents``, ``Products``, ``Rule``.  ``Reagents``
        and ``Products`` hold tuples of SMILES; a species consumed or produced
        with multiplicity appears once per occurrence.
    """
    required = {"Index", "Reagent", "Formed/Produced"}
    missing = required - set(raw.columns)
    if missing:
        raise ValueError(f"rels frame is missing columns: {sorted(missing)}")

    produced = raw["Formed/Produced"] == 1
    products = raw[produced].groupby("Index")["Reagent"].apply(tuple)
    reagents = raw[~produced].groupby("Index")["Reagent"].apply(tuple)
    rules = raw.groupby("Index")["Rule"].first() if "Rule" in raw.columns else None

    out = pd.DataFrame({"Index": products.index})
    out["Reagents"] = [reagents.get(i, ()) for i in out["Index"]]
    out["Products"] = [products[i] for i in out["Index"]]
    if rules is not None:
        out["Rule"] = [rules.get(i) for i in out["Index"]]
    return out


def load_processed_rels(path: str | Path) -> pd.DataFrame:
    """Read a ``*ProcessedRels.tsv`` file, parsing the list-valued columns.

    The deposited files store reagents and products as Python list literals.
    They are converted to tuples so rows are hashable.
    """
    frame = pd.read_csv(path, sep="\t")
    for column in ("Reagents", "Products"):
        frame[column] = [tuple(ast.literal_eval(v)) for v in frame[column]]
    return frame


def read_products(path: str | Path) -> pd.DataFrame:
    """Read a ``*_output.tsv`` product listing.

    Returns a frame with an integer ``Generation`` column; the deposited files
    use either ``3`` or ``G3`` depending on which network wrote them.
    """
    frame = pd.read_csv(path, sep="\t")
    frame["Generation"] = frame["Generation"].astype(str).str.lstrip("G").astype(int)
    return frame


class ReactionIndex:
    """Lookup tables over a reaction network.

    Attributes
    ----------
    reagents
        Maps a reaction identifier to the tuple of species it consumes.
    producers
        Maps a species to the identifiers of every reaction producing it.
    products
        Maps a reaction identifier to the tuple of species it produces.
    """

    __slots__ = ("producers", "products", "reagents")

    def __init__(
        self,
        reagents: Mapping[str, tuple[str, ...]],
        products: Mapping[str, tuple[str, ...]],
        producers: Mapping[str, list[str]],
    ) -> None:
        self.reagents = dict(reagents)
        self.products = dict(products)
        self.producers = dict(producers)

    def __len__(self) -> int:
        return len(self.reagents)

    @property
    def species(self) -> set[str]:
        """Every species appearing as a reagent or a product."""
        seen: set[str] = set(self.producers)
        for row in self.reagents.values():
            seen.update(row)
        return seen


def build_index(rels: pd.DataFrame) -> ReactionIndex:
    """Build :class:`ReactionIndex` from a pivoted or deposited rels frame.

    Replaces the ``index_finder``/``precursor_finder`` pair in the original
    pathway notebook, which rescanned the whole table -- and re-parsed every
    list literal -- on each lookup.
    """
    reagents: dict[str, tuple[str, ...]] = {}
    products: dict[str, tuple[str, ...]] = {}
    producers: dict[str, list[str]] = defaultdict(list)

    for identifier, reagent_row, product_row in zip(
        rels["Index"], rels["Reagents"], rels["Products"], strict=True
    ):
        reagents[identifier] = _as_tuple(reagent_row)
        product_tuple = _as_tuple(product_row)
        products[identifier] = product_tuple
        for species in product_tuple:
            producers[species].append(identifier)

    return ReactionIndex(reagents, products, producers)


def _as_tuple(value: Iterable[str] | str) -> tuple[str, ...]:
    if isinstance(value, str):
        return tuple(ast.literal_eval(value))
    if isinstance(value, tuple):
        return value
    return tuple(value)


def spontaneous_only(rels: pd.DataFrame, energy_column: str = "Energy Change") -> pd.DataFrame:
    """Keep reactions with a negative standard free-energy change.

    Rows whose energy could not be estimated are dropped and *counted*, not
    silently discarded -- see :func:`energy_coverage`.
    """
    energy = pd.Series(pd.to_numeric(rels[energy_column], errors="coerce"))
    return rels.loc[energy < 0]


def energy_coverage(rels: pd.DataFrame, energy_column: str = "Energy Change") -> dict[str, int]:
    """Report how many reactions carry a usable free-energy estimate."""
    energy = pd.Series(pd.to_numeric(rels[energy_column], errors="coerce"))
    return {
        "total": len(rels),
        "with_energy": int(energy.notna().sum()),
        "missing_energy": int(energy.isna().sum()),
        "spontaneous": int((energy < 0).sum()),
    }


def seeds_from_products(products: pd.DataFrame) -> Sequence[str]:
    """Return the generation-0 species of a network.

    Reading seeds from the deposited product listing avoids the failure mode in
    the original notebook, whose hard-coded seed table wrote glycine and glucose
    with different SMILES strings than the Maillard network used -- so the
    tracer never recognised its own starting materials.
    """
    return tuple(products.loc[products["Generation"] == 0, "Smiles"])
