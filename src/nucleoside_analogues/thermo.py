"""Reaction free energies, with uncertainty carried through.

Estimates standard transformed reaction Gibbs energies (ΔrG′°) with
eQuilibrator's component-contribution method, and — unlike the original
pipeline — keeps the uncertainty attached to every value.

Why that matters here
---------------------
Downstream analysis filters reactions on ``ΔrG′° < 0``. Applied to a point
estimate that is a hard threshold on a noisy quantity. Component contribution
carries a per-reaction standard error that is frequently large for compounds
outside its training set — which is most of a prebiotic network. On the
deposited generation-3 data roughly half of all reactions classified
spontaneous have ``|ΔrG′°| < 20 kJ/mol``, comparable to that error.

:func:`classify` therefore returns three outcomes rather than two:
``spontaneous``, ``non_spontaneous`` and ``undetermined``. A reaction whose
confidence interval straddles zero is not evidence either way, and saying so is
a stronger result than picking a side.

Failure handling
----------------
Every reaction gets a :class:`ReactionEnergy` with an explicit ``status``.
Nothing is silently dropped or written as a bare ``NaN``: a compound the cache
cannot decompose is reported as ``compound_missing`` and stays countable. The
original code caught every exception and appended the string ``"NaN"``, which
made an unusable estimate indistinguishable from a genuine result of zero.

Units are kJ/mol throughout, and asserted rather than assumed.

Requires the ``thermo`` extra::

    uv sync --extra thermo
"""

# equilibrator-api ships no type information and is an optional extra, so it is
# not resolvable in the default environment.
# pyright: reportMissingImports=false, reportMissingTypeStubs=false

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:  # pragma: no cover - import cost only matters at runtime
    from equilibrator_api import ComponentContribution

__all__ = [
    "DEFAULT_Z",
    "ReactionEnergy",
    "Spontaneity",
    "classify",
    "compound_cache",
    "reaction_energies",
    "summarise",
]

Spontaneity = Literal["spontaneous", "non_spontaneous", "undetermined"]
Status = Literal["ok", "compound_missing", "decomposition_failed", "estimation_failed"]

#: Two-sided 95% coverage. A reaction counts as spontaneous only when the whole
#: interval lies below zero.
DEFAULT_Z = 1.96


@dataclass(frozen=True, slots=True)
class ReactionEnergy:
    """A single reaction's free-energy estimate.

    Attributes
    ----------
    index
        Reaction identifier, matching the rels table.
    dg_prime
        ΔrG′° in kJ/mol, or ``None`` when it could not be estimated.
    uncertainty
        One standard error in kJ/mol, or ``None``.
    status
        Why a value is missing, when it is.
    detail
        Free text — typically the offending SMILES.
    """

    index: str
    dg_prime: float | None
    uncertainty: float | None
    status: Status
    detail: str = ""

    def spontaneity(self, z: float = DEFAULT_Z) -> Spontaneity | None:
        if self.dg_prime is None:
            return None
        return classify(self.dg_prime, self.uncertainty, z=z)


def classify(
    dg_prime: float,
    uncertainty: float | None,
    z: float = DEFAULT_Z,
) -> Spontaneity:
    """Three-way spontaneity call.

    ``spontaneous`` when the upper end of the interval is below zero,
    ``non_spontaneous`` when the lower end is above it, and ``undetermined``
    when the interval spans zero.

    With ``uncertainty=None`` this degrades to the sign of the point estimate,
    reproducing the original binary filter — useful for comparing against the
    deposited results, but not what should be reported.
    """
    if uncertainty is None:
        return "spontaneous" if dg_prime < 0 else "non_spontaneous"
    margin = z * uncertainty
    if dg_prime + margin < 0:
        return "spontaneous"
    if dg_prime - margin > 0:
        return "non_spontaneous"
    return "undetermined"


def compound_cache(
    smiles: Sequence[str],
    cc: ComponentContribution | None = None,
    bypass_chemaxon: bool = True,
) -> tuple[dict[str, object], list[str]]:
    """Resolve SMILES to eQuilibrator compounds.

    Returns the mapping and, separately, the SMILES that could not be resolved,
    so the caller can report coverage instead of discovering gaps as silent
    ``NaN`` values later.

    ``bypass_chemaxon`` skips the licensed ChemAxon pKa determination, as the
    original pipeline did. Protonation states are then assumed rather than
    computed, which is a real caveat for the pH 7.4 transform and should be
    stated wherever these numbers appear.
    """
    from equilibrator_api import ComponentContribution
    from equilibrator_assets.generate_compound import get_or_create_compound

    engine = ComponentContribution() if cc is None else cc
    unique = list(dict.fromkeys(smiles))
    compounds = get_or_create_compound(
        engine.ccache, unique, mol_format="smiles", bypass_chemaxon=bypass_chemaxon
    )

    resolved: dict[str, object] = {}
    missing: list[str] = []
    for entry, compound in zip(unique, compounds, strict=True):
        if compound is None or getattr(compound, "inchi_key", None) is None:
            missing.append(entry)
        else:
            resolved[entry] = compound
    return resolved, missing


def reaction_energies(
    reactions: Iterable[tuple[str, Sequence[str], Sequence[str]]],
    cc: ComponentContribution | None = None,
    bypass_chemaxon: bool = True,
) -> list[ReactionEnergy]:
    """Estimate ΔrG′° with uncertainty for each reaction.

    Parameters
    ----------
    reactions
        Triples of ``(identifier, reagent SMILES, product SMILES)``. Species
        appearing more than once are counted with their multiplicity.
    cc
        An existing :class:`ComponentContribution`; one is built if omitted.
        Constructing it downloads a ~1.3 GB compound cache on first use.

    Returns
    -------
    list[ReactionEnergy]
        One entry per input reaction, in order, always. Failures carry a status
        rather than being omitted.
    """
    from equilibrator_api import ComponentContribution, Reaction

    engine = ComponentContribution() if cc is None else cc
    rows = list(reactions)

    every_smiles = [s for _, reagents, products in rows for s in (*reagents, *products)]
    resolved, missing = compound_cache(every_smiles, cc=engine, bypass_chemaxon=bypass_chemaxon)
    unresolved = set(missing)

    out: list[ReactionEnergy] = []
    for identifier, reagents, products in rows:
        blocked = [s for s in (*reagents, *products) if s in unresolved]
        if blocked:
            out.append(
                ReactionEnergy(identifier, None, None, "compound_missing", "; ".join(blocked[:3]))
            )
            continue

        stoichiometry: dict[object, int] = {}
        for smiles in reagents:
            compound = resolved[smiles]
            stoichiometry[compound] = stoichiometry.get(compound, 0) - 1
        for smiles in products:
            compound = resolved[smiles]
            stoichiometry[compound] = stoichiometry.get(compound, 0) + 1

        try:
            estimate = engine.standard_dg_prime(Reaction(stoichiometry))
        except Exception as error:  # noqa: BLE001 - reported, never swallowed
            out.append(
                ReactionEnergy(identifier, None, None, "estimation_failed", str(error)[:200])
            )
            continue

        value = estimate.value.m_as("kJ/mol")
        error_bar = estimate.error.m_as("kJ/mol")
        out.append(ReactionEnergy(identifier, float(value), float(error_bar), "ok"))

    return out


def summarise(energies: Iterable[ReactionEnergy], z: float = DEFAULT_Z) -> dict[str, int]:
    """Count outcomes, so coverage can be reported rather than inferred."""
    counts = {
        "total": 0,
        "estimated": 0,
        "spontaneous": 0,
        "non_spontaneous": 0,
        "undetermined": 0,
        "compound_missing": 0,
        "decomposition_failed": 0,
        "estimation_failed": 0,
    }
    for energy in energies:
        counts["total"] += 1
        if energy.dg_prime is None:
            counts[energy.status] += 1
            continue
        counts["estimated"] += 1
        counts[classify(energy.dg_prime, energy.uncertainty, z=z)] += 1
    return counts
