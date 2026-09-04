"""Prebiotic nucleoside-analogue reachability analysis.

Public API
----------
:mod:`~nucleoside_analogues.rels`
    Read and reshape MØD reaction-network output.
:mod:`~nucleoside_analogues.hyperpath`
    Exact shortest synthetic pathways via minimum-weight hyperpaths.
:mod:`~nucleoside_analogues.matching`
    Stereochemistry-flattened InChIKey matching.
:mod:`~nucleoside_analogues.invariants`
    Structural and chemical soundness checks.
:mod:`~nucleoside_analogues.descriptors`
    Physicochemical descriptor panel.
:mod:`~nucleoside_analogues.thermo`
    Reaction free energies with uncertainty (requires the ``thermo`` extra).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from . import descriptors, hyperpath, invariants, matching, pka, rels

if TYPE_CHECKING:
    from . import thermo

__version__ = "0.2.0"

__all__ = ["descriptors", "hyperpath", "invariants", "matching", "pka", "rels", "thermo"]


def __getattr__(name: str):
    """Import :mod:`thermo` lazily; it needs the optional ``thermo`` extra."""
    if name == "thermo":
        from . import thermo as module

        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
