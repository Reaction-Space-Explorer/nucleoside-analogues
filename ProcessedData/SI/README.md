# Supporting information tables

Regenerate with `uv run python scripts/make_si_tables.py`. Nothing here is
edited by hand.

| File | Contents |
|---|---|
| `SI_Table1_routes.csv` | Minimum-cost route count and critical-reaction count per target |
| `SI_Table2_steps.csv` | Chain depth and reaction count per target |
| `SI_Table3_pH_robustness.csv` | Spontaneity of all five G3 networks at pH 7.0, 7.4, 9.0 and 11.0 |
| `<network>_G3_energies_pH7.4.csv` | Per-reaction dGr'o, sigma, estimability and both spontaneity calls |

## SI Table 3

Free energies from eQuilibrator component contribution at I = 0.25 M, pMg 3.0,
298.15 K, with pKa values assigned by `nucleoside_analogues.pka` for compounds
created here and taken from the compound cache otherwise.

`unestimable` counts reactions whose stoichiometric vector leaves a residual
outside the span of both the reactant- and group-contribution spaces; component
contribution assigns these infinite variance and they are excluded rather than
given a value.

`identical_to_pH7` compares the *membership* of the spontaneous set, not its
size.

`titratable_7_11` counts species whose protonation changes between pH 7 and 11
according to `nucleoside_analogues.pka`. CO2 counts: eQuilibrator's CO2 is
total dissolved inorganic carbon and its second carbonate constant (10.33) is
in range, which is why Glucose and PyruvicAcid show one titratable species
rather than none. `tests/test_pka.py` recomputes this column, so it cannot
drift from `pka.py` again -- it did once, having been written before the
carbonate constants were added.

It is identical at every pH for Formose and Glucose. It changes for the two
ammonia-seeded networks, whose products carry aliphatic amines (pKa ~9-11), and
for PyruvicAcid, where decarboxylation is common and the carbonic acid constants
(6.35, 10.33) both fall in range.

A reaction only changes classification if it is titratable *and* already within
the pH-induced shift of zero. For CO2-releasing reactions that shift is about
23 kJ/mol; in PyruvicAcid exactly 11 reactions meet both conditions, and exactly
11 change.

## SI Tables 1 and 2

Both are reported on two bases: `estimable_only` admits a reaction only when
component contribution returns a free energy and that energy is negative;
`with_unestimable` additionally admits reactions whose energy cannot be
estimated. Reactions absent from the energy file could not be constructed at
all and count as unestimable — one reaction in PyruvicAcid.

`chain_depth` is the longest chain of consecutive reactions from a seed
molecule; `reactions` is the total in the derivation, counting a shared
intermediate once per use. Both are exact minima of their objectives, and they
differ whenever a route is convergent.

`critical_reactions` counts reactions whose individual removal leaves the
target unreachable by any route of any length — robustness in the sense of
network expansion, not merely of the shortest route. A critical reaction lies
on every derivation and so on the traced route; `tests/test_hyperpath.py`
asserts that, and re-verifies each removal.
