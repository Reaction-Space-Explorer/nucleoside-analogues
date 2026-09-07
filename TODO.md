# TODO

## Running on bizon

- **Four-pH energies** (`compute_energies.py --workers 8`, pH 7.0/7.4/9.0/11.0). Four of five
  networks written: GlucoseAmm, FormoseAmm, Glucose, PyruvicAcid. **Formose G6 is last in the
  queue and still going**, about an hour and a quarter into its own pass; its four files exist but
  are empty, because each is written only when its network finishes. The local pH 7.4 copies are
  intact and committed, so nothing already computed is at risk.

## Blocked on that run

- **SI Table 3** still holds the generation-3 set (582 Formose reactions). Its caption says so and
  must change with it.
- dGr distribution by pH and by reaction rule.
- Max-min driving force on the traced routes. The spike passed at 42-106 kJ/mol on seven routes.
- Descriptors and SHAP, dropping the constant features.

## Open, needing a decision rather than work

- **Do the route figures replace Figure 6?** Sixteen are current in `figures/routes/`, one per
  spontaneously reachable network and target, drawn from the data by autocycle. Figure 6 is still
  the hand-drawn version. Replacing it changes that figure's panel structure, which is why it has
  not been done unasked.

## Open, needing something from outside

- **HMDB.** Download "Structures" (`structures.zip`) from <https://www.hmdb.ca/downloads> and unzip
  to `OriginalData/reference_databases/hmdb_structures.sdf`. HMDB returns 403 to scripted requests.
  `scripts/database_matches.py` picks it up automatically and is otherwise finished.
- **KEGG** cannot be re-derived; its bulk data is licensed. The deposited counts stay, flagged.

## Author-side

- ORCIDs (8 comments outstanding).
- TOC graphic.
- Conflict of interest and funding statements.
- Zenodo DOI for the archived repository.
- Delete the reviewer-suggestion paragraph (p4).
- Figure 5 content.
- Put `Jim_NA/` under version control. It holds the manuscript, the red-marked copy and INDEX.md,
  none of it versioned.

## Offshoots this work could support

- **Stochastic kinetics on these networks.** Lauber et al. need two global parameters, not one per
  rule, and take dfG from eQuilibrator exactly as we do. The obstacle is that this repository holds
  MØD's output and not its grammar; the rules live in reac-space-exp, which is ours.
- **Kinetic ordering of the traced routes** once barriers exist, against the reachability ordering.
- **Integer hyperflow**, following Abel et al. (ref 73), asking which traced routes are
  stoichiometrically realizable, closing the caveat rather than stating it.
- **The rule-dependence method itself** applied more widely: removing a disputed mechanism and
  re-searching is general, and the disagreement it settles here is not specific to sugars.

## Settled this round, for the record

98 references, first-citation order verified. Glyoxylate added as a fifth target (Bean 81,
Krishnamurthy & Liotta 83). Sutton 86, Tabata 87, Lauber 88, Cruz 91, Abel 73. Figure 2 rebuilt
from the reaction listings after its Formose row proved stale; the reachability funnel became
Figure 4; Figure 7 uncrowded to two panels; Figure S2 (Böttcher) cut, its content now in the
discussion as percentiles and an AUROC, and the remaining SI figures renumbered S1-S4. Formose G6
matching closed, 3,404 matched species to 9,305. Apiose, ethylene glycol, threitol, glyceric and
glycolic acid were checked for backbone literature and have none, so they are not targets.
