"""Matching nucleoside-analogue libraries against network products.

Both sides are reduced to the first 14 characters of the InChIKey, which encode
molecular constitution and charge but not stereochemistry.  The network
expansions do not track stereochemistry, so comparing full keys would miss
every match; flattening both sides to the skeleton layer is what makes the
comparison meaningful.

The flattening is worth stating plainly in any downstream analysis: a match
means the network produced *some* stereoisomer of the target scaffold, not that
it produced the specific one.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

from rdkit import Chem

from ._rdkit import silence_rdkit

silence_rdkit()

__all__ = ["find_matches", "skeleton_key", "skeleton_keys"]


def skeleton_key(smiles: str) -> str | None:
    """First 14 InChIKey characters -- constitution and charge, no stereochemistry."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol)[:14]


def skeleton_keys(smiles: Iterable[str]) -> dict[str, str]:
    """Map each parseable SMILES to its skeleton key."""
    out: dict[str, str] = {}
    for entry in smiles:
        key = skeleton_key(entry)
        if key is not None:
            out[entry] = key
    return out


def find_matches(
    library: Iterable[str],
    products: Mapping[str, int],
) -> dict[str, tuple[str, int]]:
    """Intersect an analogue library with network products.

    Parameters
    ----------
    library
        SMILES of the target analogue library.
    products
        Network product SMILES mapped to generation of first appearance.

    Returns
    -------
    dict
        Skeleton key -> ``(network SMILES, generation)``, keeping the earliest
        generation when several products share a skeleton.
    """
    by_key: dict[str, tuple[str, int]] = {}
    for smiles, generation in products.items():
        key = skeleton_key(smiles)
        if key is None:
            continue
        current = by_key.get(key)
        if current is None or generation < current[1]:
            by_key[key] = (smiles, generation)

    matches: dict[str, tuple[str, int]] = {}
    for entry in library:
        key = skeleton_key(entry)
        if key is not None and key in by_key:
            matches[key] = by_key[key]
    return matches
