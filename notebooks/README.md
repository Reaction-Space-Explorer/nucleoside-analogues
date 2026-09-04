# Notebooks

**These are a record of how the published results were produced. They are not
the supported pipeline** — use `src/nucleoside_analogues/` for new work, and
these to see how a number or a figure came about.

## Status: none of them run from a clone

This is pre-existing, not a consequence of the repository reorganisation. It
was verified against the state before that work: **0 of 21 relative data paths
resolved**, and eight of the nine pipeline notebooks address data through
Google Colab mount points such as

```
/content/drive/MyDrive/BMSIS /MinimalDirectory/ProcessedData/...
```

The `analysis/` notebooks likewise address data by relative path from their
original working directories in the archived repository, which had a different
layout from this one.

Making them runnable would mean rewriting their I/O layer. That is deliberately
not done here: editing the code would misrepresent what was actually executed,
and the package already supersedes their logic.

## `pipeline/` — the generation-3 analysis as originally run

| Notebook | Role | Superseded by |
|---|---|---|
| `PreppingAnalogues.ipynb` | `Cl` → `OH`/`NH2` substitution on the analogue library | — |
| `MatchingScript.ipynb` | InChIKey matching, library against network products | `matching` |
| `FindTargetSmiles.ipynb` | Locates target nucleosides in the networks | `matching` |
| `ProcessRels.ipynb` | Reshapes MØD long-format output | `rels.pivot_rels` |
| `GenThermodata.ipynb` | eQuilibrator ΔfG′° → ΔrG′° per reaction | — |
| `SpontaneousFilter.ipynb` | Keeps reactions with ΔrG′° < 0 | `rels.spontaneous_only` |
| `SpontVsNonSpont.ipynb` | Spontaneous vs total match counts | — |
| `PathwayScript.ipynb` | Exhaustive pathway enumeration | `hyperpath.count_minimal_routes` |
| `FindingShortestPathwayScript.ipynb` | Shortest-pathway tree search | `hyperpath.shortest_pathways` |

Three behaviours of these notebooks are worth knowing before you read their
output, and are why the package exists:

- **Pathway length is tree depth**, not a count of reactions. For a convergent
  synthesis the two differ, in 392 of the 1,181 deposited pathways (33%).
- **A timeout is reported identically to a genuine absence.** The `SIGALRM`
  handler swallows the exception, so a search that ran out of time and one that
  proved no pathway exists both emit "No spontaneous pathways".
- **`PathwayScript.ipynb` writes output files named after raw SMILES**
  (`{outdir}/*{smiles}{network}Pathways.txt`). Those filenames cannot be checked
  out on Windows. The 1,290 such files in this repository have been renamed —
  see `ProcessedData/manifests/pathway_filename_map.csv` — and CI now rejects
  any that reappear. **Re-running this notebook unmodified will recreate them.**

## `analysis/` — descriptor, plotting and database work

Salvaged from the archived
[ssiddhantsharma/nucleoside-analogues](https://github.com/ssiddhantsharma/nucleoside-analogues),
which holds the 2021–2023 BMSIS Young Scientist Program work and is now
read-only. These notebooks existed nowhere else: this repository's `Figures/`
directory was deleted in September 2024.

| Notebook | Produces |
|---|---|
| `DataPlots.ipynb` | Matches-vs-generation figures (manuscript Figure 4) |
| `DescriptorCalculator.ipynb` | Physicochemical descriptors → `figures/DescriptorPlots/` |
| `DatabaseMatchingScript.ipynb` | HMDB / ECMDB / KEGG overlap → `figures/DatabaseHitsPlots/` |
| `PlottingStructIsomers.ipynb` | Structural-isomer counts per network |
| `RemoveDegeneracy.ipynb` | Collapses stereoisomers to constitutional isomers |
| `TableMaker.ipynb`, `DatawarriorMatchingScript.ipynb`, `SolubilityCalculator.ipynb` | Exploratory |

`DatabaseMatchingScript.ipynb` additionally needs HMDB, ECMDB and KEGG reference
files that were never committed and are not redistributable here. Only the
computed matches are present, in `ProcessedData/DatabaseMatches/`.
