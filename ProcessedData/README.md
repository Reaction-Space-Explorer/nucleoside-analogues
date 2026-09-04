# ProcessedData

Everything here is derived from `OriginalData/`. Nothing is edited by hand.

Directory names are kept as published: the manuscript and its SI link to
`SpontaneousSMILES/` by URL, and `manifests/pathway_filename_map.csv` records
the one rename that has happened. Renaming a directory here breaks a citation,
so don't.

## Current pipeline

Regenerate with the scripts named; these are what the manuscript reports.

| Path | Contents | Produced by |
|---|---|---|
| `SI/SI_Table1_routes.csv` | Minimum-cost route count and critical-reaction count per target | `scripts/make_si_tables.py` |
| `SI/SI_Table2_steps.csv` | Chain depth and reaction count per target | `scripts/make_si_tables.py` |
| `SI/SI_Table3_pH_robustness.csv` | Spontaneity of each G3 network at pH 7.0, 7.4, 9.0, 11.0 | see `SI/README.md` |
| `SI/<network>_G3_energies_pH7.4.csv` | Per-reaction ΔrG′°, σ, estimability, both spontaneity calls | `scripts/compute_energies.py` |
| `SI/full/<network>_G<n>_energies_pH7.4.csv` | The same, to each network's deepest generation | `scripts/compute_energies.py` |
| `manifests/pathway_filename_map.csv` | 1,290 pathway files renamed off raw SMILES, with SHA1 | one-off, recorded |

## Inputs the pipeline reads

| Path | Contents |
|---|---|
| `RelsFiles/<network>/<network>G<n>ProcessedRels.tsv` | MØD output reshaped to one row per reaction: `Index, Reagents, Products, Rule`. Cumulative — generation *n* contains every earlier reaction |
| `MatchesFiles/<network>Matches.tsv` | Analogue library matched to network products on the InChIKey first block |
| `SpontaneousSMILES/<network>.tsv` | Species reachable by spontaneous reactions, per generation. **Linked from the SI** |

## Superseded, kept for provenance

These record how the published numbers were produced. They are not regenerated
and should not be used for new work: the energies predate the uncertainty and
unestimable handling, and the pathway files come from the threshold-limited
tree search that `hyperpath` replaces.

| Path | Contents | Superseded by |
|---|---|---|
| `RelsWithThermoFiles/` | Rels with a single `Energy Change` column, no uncertainty | `SI/<network>_G3_energies_pH7.4.csv` |
| `SpontaneousRelsWithThermoFiles/` | The above, filtered to ΔrG′° < 0 | as above, with a three-way call |
| `G3Pathways/`, `G3SinglePathways/` | Enumerated routes, one file per target molecule | `hyperpath.count_minimal_routes` |
| `G3TargetNucleosidePathways/` | Enumerations under a 10,000-route ceiling, hence lower bounds | `SI/SI_Table1_routes.csv` |
| `ShortestPathways/` | Shortest route per target, from the tree search | `SI/SI_Table2_steps.csv` |
| `PathwaySummaries/` | Whether any route was found per target | the `reachable` column of both SI tables |

## Analysis outputs

| Path | Contents |
|---|---|
| `ComplexityData/<network>ComplexityData.tsv` | Böttcher molecular complexity per matched species, by generation |
| `DescriptorsFiles/<network>Descriptors.tsv` | Physicochemical descriptors per species |
| `DatabaseMatches/<network>DatabaseMatches.tsv` | Overlap with HMDB, KEGG and ECMDB. The reference databases are not redistributable and are not vendored |
| `Nucleoside_Stereoisomers.tsv` | Enumerated stereoisomers of the analogue library. 96 MB — near GitHub's 100 MB file limit, so do not append to it |
| `G1SpontaneousMatches.tsv` | Generation-1 spontaneous matches |
