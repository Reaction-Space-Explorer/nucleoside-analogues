"""RDKit interop helpers.

RDKit ships incomplete type stubs; ``RDLogger.DisableLog`` exists at runtime but
is not declared.  Confining the suppression here keeps a single ``type: ignore``
in the codebase instead of one per module.
"""

from __future__ import annotations

from rdkit import RDLogger

__all__ = ["silence_rdkit"]

_SILENCED = False


def silence_rdkit() -> None:
    """Suppress RDKit's parser warnings. Idempotent."""
    global _SILENCED
    if _SILENCED:
        return
    RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]
    _SILENCED = True
