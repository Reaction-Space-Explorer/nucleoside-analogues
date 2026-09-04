from __future__ import annotations

import pytest
from helpers import REPO, products_path, requires

from nucleoside_analogues.pka import assign, titratable, titrates_in_range
from nucleoside_analogues.rels import read_products

# Every network contains CO2, whose carbonic acid constants (6.35, 10.33) fall in
# range, so no network is free of titratable species. What separates them is how
# many: the ammonia-seeded ones carry aliphatic amines throughout.
LOW = ["Formose", "Glucose", "PyruvicAcid"]
AMINE = ["FormoseAmm", "GlucoseAmm"]


def test_glycine_gets_both_constants() -> None:
    assert assign("NCC(=O)O") == [9.60, 3.80]


def test_carbon_dioxide_is_the_only_titratable_species_in_sugar_chemistry() -> None:
    assert assign("O=C=O") == [10.33, 6.35]
    assert titrates_in_range("O=C=O")


def test_sugars_do_not_titrate_in_range() -> None:
    for s in ("C(C(C(CO)O)O)=O", "C(C(C(C(C(CO)O)O)O)O)=O", "OCC=O"):
        assert not titrates_in_range(s)


def test_carboxylic_acid_is_flat_above_ph_5() -> None:
    assert assign("C(CO)(O)=O") == [3.80]
    assert not titrates_in_range("C(CO)(O)=O")


def test_amines_titrate_in_range() -> None:
    assert titrates_in_range("NCC(=O)O")
    assert titrates_in_range("NCCO")


@pytest.mark.parametrize("network", LOW)
def test_carbon_networks_titrate_only_through_co2(network: str) -> None:
    requires(products_path(network))
    smiles = read_products(products_path(network))["Smiles"]
    found = titratable(smiles)
    assert {s for s in found if assign(s) != [10.33, 6.35]} == set()


@pytest.mark.parametrize("network", AMINE)
def test_amine_networks_titrate_broadly(network: str) -> None:
    requires(products_path(network))
    smiles = read_products(products_path(network))["Smiles"]
    found = titratable(smiles)
    assert len(found) > 0.5 * len(smiles), f"{network} was expected to be amine-rich"


def test_si_table3_titratable_column_matches_pka(network: str) -> None:
    """The deposited count is recomputed here, so it cannot go stale again.

    It did once: the column was written before the carbonic acid constants
    were added, and so missed CO2 in four of the five networks.
    """
    import ast
    import csv

    import pandas as pd

    table = REPO / "ProcessedData" / "SI" / "SI_Table3_pH_robustness.csv"
    rels = REPO / "ProcessedData" / "RelsFiles" / network / f"{network}G3ProcessedRels.tsv"
    requires(table)
    requires(rels)

    frame = pd.read_csv(rels, sep="\t")
    literal = lambda value: ast.literal_eval(value) if isinstance(value, str) else tuple(value)  # noqa: E731
    species: set[str] = set()
    for reagents, products in zip(frame["Reagents"], frame["Products"], strict=True):
        species |= set(literal(reagents)) | set(literal(products))
    expected = sum(1 for s in species if titrates_in_range(s))

    rows = [r for r in csv.DictReader(table.open()) if r["network"] == network]
    assert rows, f"no rows for {network}"
    for row in rows:
        assert int(row["compounds"]) == len(species)
        assert int(row["titratable_7_11"]) == expected


#: Species where an independent predictor should agree with ASSIGNMENTS.
CROSS_CHECK = (
    ("NCC(=O)O", True, "glycine"),
    ("CCN", True, "ethylamine"),
    ("OCC(O)CO", False, "glycerol"),
    ("CC(=O)O", False, "acetic acid"),
    ("OCC(O)C(O)C=O", False, "threose"),
)


def test_assignments_agree_with_dimorphite_dl() -> None:
    """Cross-check the SMARTS table against an independent pKa predictor.

    dimorphite-dl enumerates protonation states directly; a species whose set
    of states differs between pH 7 and pH 11 is titratable in range.
    """
    dimorphite = pytest.importorskip("dimorphite_dl")
    for smiles, expected, name in CROSS_CHECK:
        low = set(dimorphite.protonate_smiles(smiles, ph_min=7.0, ph_max=7.0))
        high = set(dimorphite.protonate_smiles(smiles, ph_min=11.0, ph_max=11.0))
        assert (low != high) is expected, f"dimorphite-dl disagrees on {name}"
        assert titrates_in_range(smiles) is expected, f"ASSIGNMENTS disagree on {name}"


def test_carbon_dioxide_is_the_documented_exception() -> None:
    """dimorphite-dl treats O=C=O as neutral CO2 and sees no titration.

    eQuilibrator's CO2 is total dissolved inorganic carbon, whose second
    constant (10.33) does fall in range, so ASSIGNMENTS and dimorphite-dl
    differ here by design. The energies are unaffected either way: CO2 is
    taken from the compound cache, so specified_pkas never applies to it.
    """
    dimorphite = pytest.importorskip("dimorphite_dl")
    low = set(dimorphite.protonate_smiles("O=C=O", ph_min=7.0, ph_max=7.0))
    high = set(dimorphite.protonate_smiles("O=C=O", ph_min=11.0, ph_max=11.0))
    assert low == high
    assert titrates_in_range("O=C=O") is True
