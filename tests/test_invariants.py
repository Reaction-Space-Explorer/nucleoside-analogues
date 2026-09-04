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
