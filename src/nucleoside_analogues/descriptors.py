"""Physicochemical descriptor panel for analogue scaffolds.

Adapted from the MycoPermeNet project (Nevbarunegbe/Mycomembrane-permeability-project,
``revised/shap_analyses.py``), MIT licensed.  The descriptor selection follows
that work; the 3D descriptors are made optional here for two reasons specific
to this dataset.

**Cost.** ``globularity`` and ``plane of best fit`` need an embedded conformer
ensemble.  Generating 20 ETKDG conformers per molecule costs about 400 ms,
against 0.36 ms for the other twenty-one descriptors combined -- a factor of
roughly 1100.  On the full 161k-scaffold library that is 18 core-hours versus
one minute.

**Validity.** They are also the descriptors least defensible on this data.  The
analogue library and the network products are both flattened to constitutional
isomers, so any conformer is generated for an arbitrary stereoisomer of a
structure whose stereochemistry was deliberately discarded.  ``n_chiral_centers``
and ``fcsp3_bm`` carry the same caveat.

Two properties from the original 25-descriptor panel are unavailable here:
``logD`` and ``logS`` came from a commercial tool.  Crippen ``logP`` is included
as the closest open substitute.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import (
    QED,
    AllChem,
    FindMolChiralCenters,
    rdDistGeom,
    rdForceFieldHelpers,
    rdFreeSASA,
    rdMolDescriptors,
)
from rdkit.Chem.rdmolops import GetFormalCharge
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol

from ._rdkit import silence_rdkit

silence_rdkit()

__all__ = ["DESCRIPTOR_COLUMNS", "SHAPE_COLUMNS", "calc_descriptors"]

#: Descriptors computable from the 2D graph alone.
DESCRIPTOR_COLUMNS: tuple[str, ...] = (
    "HBA",
    "HBD",
    "HBA+HBD",
    "NumRings",
    "RTB",
    "NumAmideBonds",
    "TPSA",
    "logP",
    "MR",
    "exact_mass",
    "Csp3",
    "fmf",
    "QED",
    "HAC",
    "NumRingsFused",
    "unique_HBAD",
    "max_ring_size",
    "n_chiral_centers",
    "fcsp3_bm",
    "formal_charge",
    "abs_charge",
)

#: Descriptors that require an embedded 3D conformer ensemble.
SHAPE_COLUMNS: tuple[str, ...] = ("Globularity", "PBF")


def _fused_ring_count(mol: Chem.Mol) -> int:
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    return sum(1 for i, a in enumerate(rings) for b in rings[i + 1 :] if len(a & b) >= 2)


def _unique_hbad_atoms(mol: Chem.Mol) -> int:
    donors = {m[0] for m in mol.GetSubstructMatches(Chem.MolFromSmarts("[$([N;!H0]),$([O,S;H1])]"))}
    acceptors = {
        m[0] for m in mol.GetSubstructMatches(Chem.MolFromSmarts("[$([O,S;H0;v2]),$([N;v3])]"))
    }
    return len(donors | acceptors)


def _max_ring_size(mol: Chem.Mol) -> int:
    rings = mol.GetRingInfo().AtomRings()
    return max((len(r) for r in rings), default=0)


def _shape(mol: Chem.Mol, n_conformers: int = 20, prune_rms: float = 0.1) -> tuple[float, float]:
    """Mean globularity and plane-of-best-fit over an ETKDG conformer ensemble."""
    mol_h = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv2()
    params.pruneRmsThresh = prune_rms
    ids = rdDistGeom.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)
    if not len(ids):
        return float("nan"), float("nan")
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol_h)

    globularities: list[float] = []
    pbfs: list[float] = []
    for conf_id in ids:
        try:
            radii = rdFreeSASA.classifyAtoms(mol_h)
            sasa = rdFreeSASA.CalcSASA(mol_h, radii, confIdx=conf_id)
            volume = AllChem.ComputeMolVolume(mol_h, confId=conf_id)
            if sasa > 0:
                globularities.append(((volume * 3 / (4 * np.pi)) ** (2 / 3)) * 4 * np.pi / sasa)
            pbfs.append(rdMolDescriptors.CalcPBF(mol_h, confId=conf_id))
        except Exception:  # noqa: BLE001 - a failed conformer should not abort the molecule
            continue
    return (
        float(np.nanmean(globularities)) if globularities else float("nan"),
        float(np.nanmean(pbfs)) if pbfs else float("nan"),
    )


def calc_descriptors(
    smiles: Iterable[str],
    include_shape: bool = False,
    n_conformers: int = 20,
) -> pd.DataFrame:
    """Compute the descriptor panel for a collection of SMILES.

    Parameters
    ----------
    smiles
        Structures to describe.  Rows align with the input order.
    include_shape
        Also compute ``Globularity`` and ``PBF``.  These need conformer
        generation and cost roughly 1100x the rest of the panel combined; see
        the module docstring before enabling them.
    n_conformers
        Conformers per molecule when ``include_shape`` is set.

    Returns
    -------
    DataFrame
        One row per input.  Molecules whose SMILES do not parse yield a row of
        ``NaN`` rather than being dropped, so the frame stays aligned with the
        input and failures remain visible.
    """
    entries: Sequence[str] = list(smiles)
    columns = list(DESCRIPTOR_COLUMNS)
    if include_shape:
        columns += list(SHAPE_COLUMNS)

    rows: list[list[float]] = []
    for entry in entries:
        mol = Chem.MolFromSmiles(entry)
        if mol is None:
            rows.append([float("nan")] * len(columns))
            continue

        hba = rdMolDescriptors.CalcNumHBA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        logp, mr = rdMolDescriptors.CalcCrippenDescriptors(mol)
        scaffold = GetScaffoldForMol(mol)
        heavy = mol.GetNumHeavyAtoms()

        row: list[float] = [
            hba,
            hbd,
            hba + hbd,
            rdMolDescriptors.CalcNumRings(mol),
            rdMolDescriptors.CalcNumRotatableBonds(mol),
            rdMolDescriptors.CalcNumAmideBonds(mol),
            rdMolDescriptors.CalcTPSA(mol),
            logp,
            mr,
            rdMolDescriptors.CalcExactMolWt(mol),  # monoisotopic, not average
            rdMolDescriptors.CalcFractionCSP3(mol),
            (scaffold.GetNumHeavyAtoms() / heavy) if scaffold is not None and heavy else 0.0,
            QED.qed(mol),
            heavy,
            _fused_ring_count(mol),
            _unique_hbad_atoms(mol),
            _max_ring_size(mol),
            len(FindMolChiralCenters(mol, includeUnassigned=True)),
            rdMolDescriptors.CalcFractionCSP3(scaffold) if scaffold is not None else float("nan"),
            GetFormalCharge(mol),
            abs(GetFormalCharge(mol)),
        ]
        if include_shape:
            row.extend(_shape(mol, n_conformers=n_conformers))
        rows.append(row)

    return pd.DataFrame(rows, columns=columns)
