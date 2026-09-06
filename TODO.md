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

## Krishnamurthy glyoxylate papers, all verified against Crossref

Supplied 2026-09-07. Not yet cited; each would need renumbering, so they are held together.

- Krishnamurthy, R.; Liotta, C. L. The potential of glyoxylate as a prebiotic source molecule and
  a reactant in protometabolic pathways - the glyoxylose reaction. *Chem* **2023**, 9, 784-797.
  10.1016/j.chempr.2023.03.007. Proposes glyoxylate as an alternative source molecule to
  formaldehyde. Supports the new glyoxylate target far better than ref 81 alone.
- Sutton, Pulletikurti, Lin, Krishnamurthy, Liotta. Abiotic aldol reactions of formaldehyde with
  ketoses and aldoses. *Chem* **2025**, 11, 102553. 10.1016/j.chempr.2025.102553. By 13C NMR:
  formaldehyde aldol dominates, tetroses and pentoses are not observed, carbonyl migration is
  absent, and the Breslow autocatalytic pathway is doubted. **This bears directly on p98, which
  currently says our traced routes are the chain-growth steps "as Breslow described it".**
- Cruz, H.; Krishnamurthy, R. Selection of ribofuranose-isomer among pentoses by phosphorylation
  with diamidophosphate. *Angew. Chem. Int. Ed.* **2025**, 64. 10.1002/anie.202509810. Furanose
  against pyranose selection, a distinction this representation cannot make.
- Clay, Cooke, Kumar, Yadav, Krishnamurthy, Springsteen. A plausible prebiotic one-pot synthesis
  of orotate and pyruvate. *Angew. Chem. Int. Ed.* **2022**, 61. 10.1002/anie.202112572.
  Hydantoin plus glyoxylate; hydantoin is absent from all five CRNRs.

**Carbonyl-migration robustness is computed** (`scripts/rule_dependence.py`,
`ProcessedData/SI/rule_dependence.csv`) but not yet written into the manuscript. Removing the two
migration rules from the spontaneous network: the formose routes to ribose (2), threose (1) and
glycerol (2) are unchanged, so the headline result does not rest on carbonyl migration at all.
Formose deoxyribose and glyoxylate are lost, and every target in G, GA and PA is lost. FA keeps
all five, deoxyribose lengthening from 3 to 4.

## Author-side

- ORCIDs (8 comments outstanding).
- TOC graphic.
- Conflict of interest and funding statements.
- Zenodo DOI for the archived repository.
- Delete the reviewer-suggestion paragraph (p4).
- Figure 5 content.
- Put `Jim_NA/` under version control.
