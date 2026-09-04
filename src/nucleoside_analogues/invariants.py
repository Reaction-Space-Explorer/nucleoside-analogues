"""Chemical soundness checks for a generated reaction network.

Three tiers, in increasing order of what they can actually tell you:

1. **Structural invariants** -- every SMILES parses, every reaction conserves
   atoms and charge, every molecule referenced by a reaction exists in the
   product listing, and no reaction consumes a reagent from a later generation.
   Cheap, and they catch a broken rule set immediately.
2. **Implausibility screens** -- motifs a prebiotic aqueous network should not
   be emitting (peroxides, N-N and N-O bonds, orthoesters, strained rings).
   A blocklist cannot prove a molecule is real, only that it is not obviously
   wrong.
3. **Experimental recall** -- does the network contain species actually
   observed in the laboratory reaction it models?  This is the only tier that
   validates rather than self-checks; see :mod:`nucleoside_analogues.recall`.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

from rdkit import Chem

from ._rdkit import silence_rdkit

silence_rdkit()

__all__ = [
    "IMPLAUSIBLE_MOTIFS",
    "BalanceReport",
    "atom_balance",
    "composition",
    "screen_motifs",
]


#: SMARTS for substructures that should be absent, or vanishingly rare, in an
#: aqueous prebiotic network.  Values are human-readable names.
IMPLAUSIBLE_MOTIFS: dict[str, str] = {
    "peroxide (O-O)": "[OX2][OX2]",
    "N-N single bond": "[NX3][NX3]",
    "N-O single bond": "[NX3][OX2]",
    "geminal diol": "[CX4]([OX2H1])([OX2H1])",
    "orthoester": "[CX4]([OX2])([OX2])([OX2])",
    "orthocarbonate": "[CX4]([OX2])([OX2])([OX2])[OX2]",
    "hemiaminal": "[CX4]([OX2H1])([NX3])",
    "geminal diamine": "[CX4]([NX3])([NX3])",
    "sp2 carbon in 3-ring": "[CX3;R3]",
    "sp2 carbon in 4-ring": "[CX3;R4]",
    "allene": "[CX2](=[CX3])=[CX3]",
}


def composition(smiles: str) -> Counter[str] | None:
    """Element counts including implicit hydrogens, plus formal charge.

    Charge is stored under the key ``"charge"`` so a single :class:`Counter`
    comparison checks both atom and charge balance.  Returns ``None`` if the
    SMILES does not parse.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    counts: Counter[str] = Counter()
    for atom in mol.GetAtoms():
        counts[atom.GetSymbol()] += 1
        counts["H"] += atom.GetTotalNumHs()
    counts["charge"] = Chem.GetFormalCharge(mol)
    return counts


@dataclass(slots=True)
class BalanceReport:
    """Result of an atom- and charge-balance sweep over a network."""

    checked: int = 0
    balanced: int = 0
    unparseable: list[str] = field(default_factory=list)
    unbalanced: list[tuple[str, dict[str, int]]] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.unparseable and not self.unbalanced

    def __str__(self) -> str:
        return (
            f"{self.balanced}/{self.checked} reactions balanced; "
            f"{len(self.unbalanced)} unbalanced, "
            f"{len(self.unparseable)} unparseable SMILES"
        )


def atom_balance(
    reactions: Iterable[tuple[str, Sequence[str], Sequence[str]]],
) -> BalanceReport:
    """Check that every reaction conserves atoms and charge.

    Parameters
    ----------
    reactions
        Triples of ``(identifier, reagents, products)`` where reagents and
        products are sequences of SMILES.

    Returns
    -------
    BalanceReport
        Offending reactions are listed with their per-element discrepancy
        (reagents minus products) rather than being counted and dropped.
    """
    report = BalanceReport()
    cache: dict[str, Counter[str] | None] = {}

    def lookup(smiles: str) -> Counter[str] | None:
        if smiles not in cache:
            cache[smiles] = composition(smiles)
        return cache[smiles]

    for identifier, reagents, products in reactions:
        report.checked += 1
        left: Counter[str] = Counter()
        right: Counter[str] = Counter()
        broken = False
        for smiles in reagents:
            counts = lookup(smiles)
            if counts is None:
                report.unparseable.append(smiles)
                broken = True
                break
            left.update(counts)
        if broken:
            continue
        for smiles in products:
            counts = lookup(smiles)
            if counts is None:
                report.unparseable.append(smiles)
                broken = True
                break
            right.update(counts)
        if broken:
            continue

        if left == right:
            report.balanced += 1
        else:
            keys = set(left) | set(right)
            delta = {k: left[k] - right[k] for k in keys if left[k] != right[k]}
            report.unbalanced.append((identifier, delta))

    return report


def screen_motifs(
    smiles: Iterable[str],
    motifs: dict[str, str] | None = None,
) -> dict[str, int]:
    """Count how many molecules match each implausibility SMARTS.

    A non-zero count is not automatically an error -- geminal diols are real
    species in water, and formaldehyde is almost entirely hydrated -- but a
    sharp change between rule-set versions is worth investigating.
    """
    patterns = {
        name: Chem.MolFromSmarts(smarts) for name, smarts in (motifs or IMPLAUSIBLE_MOTIFS).items()
    }
    hits = dict.fromkeys(patterns, 0)
    for entry in smiles:
        mol = Chem.MolFromSmiles(entry)
        if mol is None:
            continue
        for name, query in patterns.items():
            if query is not None and mol.HasSubstructMatch(query):
                hits[name] += 1
    return hits


def generation_monotonic(
    rels: Iterable[tuple[str, Sequence[str], Sequence[str]]],
    generation: dict[str, int],
    step: int,
) -> list[str]:
    """Reactions in generation *step* that consume a too-recent reagent.

    Every reagent of a generation-*n* reaction must already exist at generation
    ``n - 1`` or earlier.  Products may legitimately be older than *n*: a later
    reaction re-forming an existing compound is ordinary chemistry, and is what
    makes the reaction graph cyclic.
    """
    offenders = []
    for identifier, reagents, _products in rels:
        if any(generation.get(r, step) > step - 1 for r in reagents):
            offenders.append(identifier)
    return offenders
