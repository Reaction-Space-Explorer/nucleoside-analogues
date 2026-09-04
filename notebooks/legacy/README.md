# Legacy notebooks

Salvaged from the archived
[ssiddhantsharma/nucleoside-analogues](https://github.com/ssiddhantsharma/nucleoside-analogues)
repository, which holds the 2021–2023 BMSIS Young Scientist Program work and is
now read-only. These notebooks existed nowhere else.

| Notebook | Produces |
|---|---|
| `DataPlots.ipynb` | Matches-vs-generation figures (manuscript Figure 4) |
| `DescriptorCalculator.ipynb` | Physicochemical descriptors → `../../figures/DescriptorPlots/` |
| `DatabaseMatchingScript.ipynb` | HMDB / ECMDB / KEGG overlap → `../../figures/DatabaseHitsPlots/` |
| `PlottingStructIsomers.ipynb` | Structural-isomer counts per network |
| `RemoveDegeneracy.ipynb` | Collapses stereoisomers to constitutional isomers |
| `TableMaker.ipynb`, `DatawarriorMatchingScript.ipynb`, `SolubilityCalculator.ipynb` | Exploratory |

**These are kept for provenance, not as the supported pipeline.** They address
data by relative path from their original working directories, carry no pinned
environment, and predate the corrections in `src/nucleoside_analogues/`. Use
them to see how a figure was made; use the package to make a new one.

`DatabaseMatchingScript.ipynb` additionally needs HMDB, ECMDB and KEGG
reference files that were never committed (and are not redistributable here);
only the computed matches in `../../ProcessedData/DatabaseMatches/` are present.
