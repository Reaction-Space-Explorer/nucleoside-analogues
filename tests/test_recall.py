"""Recall against an experimentally derived structure set.

Self-consistency checks cannot show that a network models real chemistry.  This
one can: it asks whether the formose network contains the species reported from
laboratory formose reactions.

The reference set is Omran/Decker, as used in Arya et al., Chem. Sci. 2022 and
distributed with the reac-space-exp repository.
"""

from __future__ import annotations

import pytest
from helpers import REPO, products_path, requires
from rdkit import Chem

from nucleoside_analogues._rdkit import silence_rdkit
from nucleoside_analogues.matching import skeleton_key
from nucleoside_analogues.rels import read_products

silence_rdkit()

TEST_SET = REPO / "tests" / "data" / "OmranDeckerFormoseTestSet.sdf"

#: Recall observed on the deposited network.  One reference structure needs
#: methanol chemistry the rule set does not cover, and one SDF record is empty.
MINIMUM_RECALL = 0.85


@pytest.mark.parametrize("network", ["Formose", "FormoseAmm"])
def test_formose_network_recovers_the_literature_set(network: str) -> None:
    if not TEST_SET.exists():
        pytest.skip("literature test set not vendored")
    requires(products_path(network))

    reference: set[str] = set()
    for mol in Chem.SDMolSupplier(str(TEST_SET)):
        if mol is None or mol.GetNumAtoms() == 0:
            continue
        reference.add(Chem.MolToInchiKey(mol)[:14])

    produced = {
        key
        for smiles in read_products(products_path(network))["Smiles"]
        if (key := skeleton_key(smiles)) is not None
    }
    recovered = reference & produced
    recall = len(recovered) / len(reference)
    assert recall >= MINIMUM_RECALL, (
        f"{network} recovers {len(recovered)}/{len(reference)} "
        f"literature structures ({recall:.0%}), below {MINIMUM_RECALL:.0%}"
    )
