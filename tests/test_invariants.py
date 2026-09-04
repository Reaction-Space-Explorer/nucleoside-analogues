"""Structural invariants over the deposited networks.

These all pass on the data as deposited.  They exist to fail loudly if a rule
set change, a re-run, or a reformatting silently alters the chemistry.
"""

from __future__ import annotations

import pytest
from helpers import processed_rels_path, products_path, requires

from nucleoside_analogues.invariants import atom_balance, generation_monotonic, screen_motifs
from nucleoside_analogues.rels import build_index, load_processed_rels, read_products


@pytest.mark.parametrize("generation", [1, 2, 3])
def test_reactions_conserve_atoms_and_charge(network: str, generation: int) -> None:
    path = processed_rels_path(network, generation)
    requires(path)
    rels = load_processed_rels(path)
    report = atom_balance(zip(rels["Index"], rels["Reagents"], rels["Products"], strict=True))
    assert not report.unparseable, f"unparseable SMILES: {report.unparseable[:5]}"
    assert not report.unbalanced, f"unbalanced reactions: {report.unbalanced[:5]}"
    assert report.balanced == report.checked


def test_every_species_appears_in_the_product_listing(network: str) -> None:
    rels_path = processed_rels_path(network, 3)
    requires(rels_path)
    requires(products_path(network))
    rels = load_processed_rels(rels_path)
    known = set(read_products(products_path(network))["Smiles"])
    index = build_index(rels)
    orphans = index.species - known
    assert not orphans, f"{len(orphans)} species absent from the product listing"


@pytest.mark.parametrize("generation", [1, 2, 3])
def test_reagents_never_come_from_a_later_generation(network: str, generation: int) -> None:
    rels_path = processed_rels_path(network, generation)
    requires(rels_path)
    requires(products_path(network))
    rels = load_processed_rels(rels_path)
    products = read_products(products_path(network))
    gen = dict(zip(products["Smiles"], products["Generation"], strict=True))
    offenders = generation_monotonic(
        zip(rels["Index"], rels["Reagents"], rels["Products"], strict=True), gen, generation
    )
    assert not offenders, f"{len(offenders)} reactions consume a too-recent reagent"


def test_no_implausible_motifs_beyond_known_trace_levels(network: str) -> None:
    """Screen for substructures an aqueous prebiotic network should not emit.

    Geminal diols and enols are real aqueous species, so a small count is
    expected; the threshold guards against a rule set that starts emitting them
    wholesale, or emitting peroxides and N-N bonds at all.
    """
    requires(products_path(network))
    smiles = read_products(products_path(network))["Smiles"]
    hits = screen_motifs(smiles)
    total = len(smiles)
    for motif in ("peroxide (O-O)", "N-N single bond", "N-O single bond", "orthoester"):
        assert hits[motif] == 0, f"{network}: {hits[motif]} molecules contain {motif}"
    for motif, count in hits.items():
        assert count / total < 0.01, f"{network}: {motif} in {count}/{total} molecules"


#: Independent RDKit calls, deliberately not routed through descriptors.py.
def test_descriptor_panel_matches_rdkit() -> None:
    """Every descriptor checked against a separate RDKit call.

    ``exact_mass`` is monoisotopic, not average molecular weight: it was
    called ``MW``, which reads as the latter and differs by 0.12 Da on
    aspirin. The name now says which it is.
    """
    from rdkit import Chem
    from rdkit.Chem import QED, Crippen, Descriptors, rdMolDescriptors  # noqa: F401

    from nucleoside_analogues.descriptors import calc_descriptors

    molecules = {
        "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "ribose": "C(C(C(C(CO)O)O)O)=O",
        "glycine": "NCC(=O)O",
    }
    frame = calc_descriptors(list(molecules.values()))
    frame.index = list(molecules)
    reference = {
        "HBA": rdMolDescriptors.CalcNumHBA,
        "HBD": rdMolDescriptors.CalcNumHBD,
        "TPSA": rdMolDescriptors.CalcTPSA,
        "logP": Crippen.MolLogP,  # pyright: ignore[reportAttributeAccessIssue]
        "MR": Crippen.MolMR,  # pyright: ignore[reportAttributeAccessIssue]
        "exact_mass": Descriptors.ExactMolWt,  # pyright: ignore[reportAttributeAccessIssue]
        "RTB": rdMolDescriptors.CalcNumRotatableBonds,
        "NumRings": rdMolDescriptors.CalcNumRings,
        "NumAmideBonds": rdMolDescriptors.CalcNumAmideBonds,
        "Csp3": rdMolDescriptors.CalcFractionCSP3,
        "QED": QED.qed,
        "HAC": lambda m: m.GetNumHeavyAtoms(),
        "formal_charge": Chem.GetFormalCharge,
    }
    for column, function in reference.items():
        for name, smiles in molecules.items():
            expected = function(Chem.MolFromSmiles(smiles))
            assert abs(float(frame.loc[name, column]) - expected) < 1e-6, f"{column}/{name}"

    assert [a + b for a, b in zip(frame["HBA"], frame["HBD"], strict=True)] == list(
        frame["HBA+HBD"]
    )
    assert [abs(v) for v in frame["formal_charge"]] == list(frame["abs_charge"])


def test_ring_and_stereo_descriptors() -> None:
    """Fused rings, ring size and stereocentres on molecules with known answers."""
    from nucleoside_analogues.descriptors import calc_descriptors

    cases = {
        "naphthalene": "c1ccc2ccccc2c1",
        "biphenyl": "c1ccc(-c2ccccc2)cc1",
        "cyclopropane": "C1CC1",
        "L-alanine": "C[C@@H](N)C(=O)O",
        "glycine": "NCC(=O)O",
    }
    frame = calc_descriptors(list(cases.values()))
    frame.index = list(cases)
    assert frame.loc["naphthalene", "NumRingsFused"] == 1
    assert frame.loc["biphenyl", "NumRingsFused"] == 0  # linked, not fused
    assert frame.loc["cyclopropane", "max_ring_size"] == 3
    assert frame.loc["L-alanine", "n_chiral_centers"] == 1
    assert frame.loc["glycine", "n_chiral_centers"] == 0
