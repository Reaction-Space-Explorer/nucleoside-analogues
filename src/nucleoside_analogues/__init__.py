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
"""

from __future__ import annotations

from . import descriptors, hyperpath, invariants, matching, rels

__version__ = "0.2.0"

__all__ = ["descriptors", "hyperpath", "invariants", "matching", "rels"]
