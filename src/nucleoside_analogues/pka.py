"""pKa assignment for the transformed free-energy calculation."""

from __future__ import annotations

from collections.abc import Iterable

from rdkit import Chem

from ._rdkit import silence_rdkit

silence_rdkit()

__all__ = ["ASSIGNMENTS", "assign", "titratable", "titrates_in_range"]

# SMARTS -> (pKa, source). Aliphatic values from Perrin, Dissociation Constants
# of Organic Bases in Aqueous Solution (1965); acid values from Lide, CRC Handbook.
# (site, SMARTS, pKa, label). Within a site the first match wins, so specific
# patterns precede general ones.
ASSIGNMENTS: tuple[tuple[str, str, float, str], ...] = (
    ("amine", "[NX3;H2][CX4][CX3](=O)[OX2H1]", 9.60, "alpha-amino acid amine"),
    ("amine", "[NX3;H1;!$(NC=O)]([CX4])[CX4]", 10.80, "secondary aliphatic amine"),
    ("amine", "[NX3;H2;!$(NC=O)][CX4]", 10.60, "primary aliphatic amine"),
    ("acid", "[CX3](=O)[OX2H1]", 3.80, "carboxylic acid"),
    # CO2 hydrates to carbonic acid; both constants lie inside pH 7-11, which is
    # why decarboxylations are pH-dependent.
    ("carbonate", "[CX2H0](=[OX1])=[OX1]", 10.33, "carbonic acid, second"),
    ("carbonate2", "[CX2H0](=[OX1])=[OX1]", 6.35, "carbonic acid, first"),
)

_PATTERNS = tuple(
    (site, Chem.MolFromSmarts(smarts), pka, label) for site, smarts, pka, label in ASSIGNMENTS
)


def assign(smiles: str) -> list[float]:
    """pKa values for a molecule, highest first. Empty if none apply."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    values: list[float] = []
    seen: set[str] = set()
    for site, pattern, pka, _ in _PATTERNS:
        if site in seen or pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            values.append(pka)
            seen.add(site)
    return sorted(values, reverse=True)


def titrates_in_range(smiles: str, low: float = 7.0, high: float = 11.0) -> bool:
    return any(low <= pka <= high for pka in assign(smiles))


def titratable(smiles: Iterable[str], low: float = 7.0, high: float = 11.0) -> set[str]:
    """Subset whose protonation state changes between *low* and *high*."""
    return {s for s in smiles if titrates_in_range(s, low, high)}
