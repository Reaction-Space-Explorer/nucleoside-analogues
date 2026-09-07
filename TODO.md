# TODO

## Running

- **Four-pH energies on bizon** (`compute_energies.py --workers 8`, pH 7.0/7.4/9.0/11.0).
  GlucoseAmm G4 done; FormoseAmm, Formose G6, Glucose G5, PyruvicAcid G6 to go.
  Compound resolution is pH-independent, so the four pH values cost one resolution pass.
  Blocks: SI Table 3 at full depth, dGr by pH and rule, MDF.

## Blocked on that run

- SI Table 3 refresh. The embedded table is still the generation-3 set (582 Formose
  reactions); the caption says so, and must be updated with it.
- dGr distribution by pH and by reaction rule.
- Max-min driving force on the traced routes. Spike passed, 42-106 kJ/mol on 7 routes.
- Descriptors + SHAP, dropping the constant features.

## Open, not blocked

- **HMDB.** Download **"Structures"** (`structures.zip`) from <https://www.hmdb.ca/downloads>
  and unzip to `OriginalData/reference_databases/hmdb_structures.sdf`. HMDB returns 403 to
  scripted requests, so this is manual. `scripts/database_matches.py` picks it up
  automatically and is otherwise finished. ChEBI is already re-derived at full depth.
- KEGG cannot be re-derived; its bulk data is licensed. The deposited counts in
  `ProcessedData/DatabaseMatches/` stay as they are, flagged in the SI.

## Done since the last update

- SI Table 4, sinks and hills, with the share of species classified alongside. That share
  runs from 9% to 96%, so the counts cannot be read without it.
- Figure S4, chain depth against reaction count, worked through the formose route to
  deoxyribose: longest chain three, whole derivation four.
- Two miscitations corrected. The SI methods cited Robertson & Miller for component
  contribution; it now cites Noor and Beber, 69 and 70, as the main text already did.
  Figure S2 cited Rogers & Hahn for Bottcher complexity; it now cites 86. Both were
  introduced by earlier insertion shifts. All 92 references remain cited.

- Oro citation checked and corrected. Ref 64 is Oro's 1961 adenine-from-HCN paper and
  said nothing about deoxyribose; the claim now cites ref 79, Teichert, Kruse & Trapp,
  Angew. Chem. Int. Ed. 2019, 58, 9944-9947, already in the reference list and verified
  against PubMed 31131499, whose abstract states the acetaldehyde condensation. All 92
  references remain cited.
- autocycle draws intermediate labels (upstream commit 195883d).

## Computed but NOT yet in the manuscript

Kept here so it cannot drift. Everything else found this round has been written in.

- **Reachability funnel figure** (`figures/output/Figure_reachability_funnel.png`). Exists and is
  current; not placed. It is a candidate to replace or accompany Figure 4, which is a placement
  decision rather than a writing one.
- **The twelve autocycle route figures** (`figures/routes/`). Current; not placed. Candidate
  replacement for the hand-drawn Figure 6, which would change that figure's panel structure.
- **Further targets are closed for now.** Glyoxylate is added. Apiose, ethylene glycol,
  threitol, glyceric acid and glycolic acid were checked against Crossref and PubMed and none has
  literature support as a prebiotic nucleic acid backbone: apiose returns only synthetic antiviral
  apiosyl nucleosides, and glycol nucleic acid is built on glycerol, already a target through ref
  80. Ribulose and fructose are supported by refs 31 and 74 but fructose sits at the 70th
  complexity percentile, a poor fit for a backbone target.
- **Figure S2** is now redundant with the numbers written into the discussion, and could be cut.

## Done: the Krishnamurthy papers are cited and the analysis is written up

Four references added (96 total, first-citation order preserved): Krishnamurthy & Liotta 82,
Sutton 85, Tabata 88 in the list, Cruz 89. The robustness result is a new paragraph after p98 and
SI Table 5. The Breslow attribution in p98 is qualified. Cruz supports the furanose/pyranose
limitation in the stereochemistry paragraph. Clay et al. 2022 is not cited: hydantoin is absent
from all five CRNRs, so it bears on nothing computed here.

## Offshoots this work could support

- **Stochastic kinetics on these networks.** Lauber et al. need two global parameters, not one per
  rule, and take dfG from eQuilibrator exactly as we do. The obstacle is that this repository holds
  MØD's output and not its grammar; the rules live in reac-space-exp, which is ours. With those,
  MØD's on-the-fly Gillespie simulation runs directly on this chemistry.
- **Kinetic ordering of the traced routes** once barriers exist, against the reachability ordering
  reported here.
- **Integer hyperflow on the same networks**, following Abel et al., to ask which of the traced
  routes are stoichiometrically realizable, closing the caveat rather than stating it.
- **The rule-dependence method itself** applied more widely: removing a disputed mechanism and
  re-searching is general, and the literature disagreement it settles here is not specific to sugars.

## Author-side

- ORCIDs (8 comments outstanding).
- TOC graphic.
- Conflict of interest and funding statements.
- Zenodo DOI for the archived repository.
- Delete the reviewer-suggestion paragraph (p4).
- Figure 5 content.
- Put `Jim_NA/` under version control.
