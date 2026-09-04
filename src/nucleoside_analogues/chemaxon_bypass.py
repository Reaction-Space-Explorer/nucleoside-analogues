"""Allow eQuilibrator compound creation without ChemAxon's cxcalc.

This module is what makes the reported energies reproducible without a
licence, so it belongs in the package rather than beside it.

equilibrator-assets 0.6 computes pKa mappings unconditionally in
``create_compounds`` and only *uses* them on the ``chemaxon`` branch, so a
missing ``cxcalc`` binary aborts compound creation even when
``bypass_chemaxon=True`` is requested.

This makes the bypass path work as documented: compounds are built from the
supplied SMILES with no dissociation constants, exactly as
``_populate_compound_information`` does on its ``bypass`` branch. It changes no
chemistry relative to what ``bypass_chemaxon=True`` already means -- protonation
states are assumed rather than computed, which must be stated wherever the
resulting energies are reported.
"""

# equilibrator-assets ships no type information and is an optional extra.
# pyright: reportMissingImports=false, reportMissingTypeStubs=false
from equilibrator_assets import chemaxon, generate_compound, thermodynamics

_ORIGINAL = generate_compound._populate_compound_information


def _bypass_only(row):
    if getattr(row, "compound_dict", None) is None and row.method == "chemaxon":
        return None  # fall through to the bypass method
    return _ORIGINAL(row)


def enable() -> None:
    if chemaxon.get_chemaxon_status() == 0:
        return  # a licensed install is present; use it
    # generate_compound reaches this through the thermodynamics module, so one
    # patch covers both call sites.
    thermodynamics.get_compound_mappings = lambda molecules, *a, **k: [None] * len(molecules)
    generate_compound._populate_compound_information = _bypass_only
