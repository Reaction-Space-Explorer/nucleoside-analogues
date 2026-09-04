from __future__ import annotations

import pytest
from helpers import products_path, requires

from nucleoside_analogues.pka import assign, titratable, titrates_in_range
from nucleoside_analogues.rels import read_products

CLEAN = ["Formose", "Glucose", "PyruvicAcid"]
AMINE = ["FormoseAmm", "GlucoseAmm"]


def test_glycine_gets_both_constants() -> None:
    assert assign("NCC(=O)O") == [9.60, 3.80]


def test_sugars_do_not_titrate_in_range() -> None:
    for s in ("C(C(C(CO)O)O)=O", "C(C(C(C(C(CO)O)O)O)O)=O", "OCC=O"):
        assert not titrates_in_range(s)


def test_carboxylic_acid_is_flat_above_ph_5() -> None:
    assert assign("C(CO)(O)=O") == [3.80]
    assert not titrates_in_range("C(CO)(O)=O")


def test_amines_titrate_in_range() -> None:
    assert titrates_in_range("NCC(=O)O")
    assert titrates_in_range("NCCO")


@pytest.mark.parametrize("network", CLEAN)
def test_clean_networks_have_no_titratable_compounds(network: str) -> None:
    requires(products_path(network))
    smiles = read_products(products_path(network))["Smiles"]
    assert titratable(smiles) == set()


@pytest.mark.parametrize("network", AMINE)
def test_amine_networks_do_have_them(network: str) -> None:
    requires(products_path(network))
    smiles = read_products(products_path(network))["Smiles"]
    assert titratable(smiles), f"{network} was expected to contain amines"
