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
- Sinks and hills as a table, not a figure. Only 9-77% of species classify, which a
  figure would hide.
- SI figure: chain depth against reaction count, worked on one convergent route.
- Confirm the Oro reference for the acetaldehyde + glyceraldehyde aldol to 2-deoxyribose
  from the primary source. The traced route recovers it unprompted, which is worth
  stating, but ref 64 has not been checked against the paper itself.

## Author-side

- ORCIDs (8 comments outstanding).
- TOC graphic.
- Conflict of interest and funding statements.
- Zenodo DOI for the archived repository.
- Delete the reviewer-suggestion paragraph (p4).
- Figure 5 content.
- Put `Jim_NA/` under version control.
