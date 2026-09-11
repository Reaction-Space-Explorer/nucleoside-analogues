# nucleoside-analogues

[![CI](https://github.com/Reaction-Space-Explorer/nucleoside-analogues/actions/workflows/ci.yml/badge.svg)](https://github.com/Reaction-Space-Explorer/nucleoside-analogues/actions/workflows/ci.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/code-BSD--3--Clause-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)

Which nucleic-acid-like molecules can prebiotic chemistry actually reach?

This repository intersects combinatorially enumerated **nucleoside analogue (NA)
libraries** with the products of in-silico **chemical reaction networks (CRNs)**,
then asks whether each matched structure is reachable by a sequence of
thermodynamically spontaneous reactions, and in how few steps.

Supports *"Computational Identification of Plausible Prebiotic Nucleic Acid
Analogue Synthetic Pathways"* (Ludwig, Sharma, Cruz-Simbron, Arya, Jain, Okoń,
Meringer & Cleaves).

## Quick start

```bash
git clone https://github.com/Reaction-Space-Explorer/nucleoside-analogues
cd nucleoside-analogues
uv sync --group dev && uv run pytest -q
```

No `uv`? `pip install -r requirements.txt`, or `docker build -t na . && docker run --rm na pytest -q`.

```python
from nucleoside_analogues.hyperpath import shortest_pathways
from nucleoside_analogues.rels import build_index, load_processed_rels, read_products, seeds_from_products

rels = load_processed_rels("ProcessedData/SpontaneousRelsWithThermoFiles/SpontaneousFormoseG3RelsWithThermo.tsv")
products = read_products("OriginalData/OriginalNetworkData/Products/formose_output.tsv")
result = shortest_pathways(build_index(rels), seeds_from_products(products))
result.cost["C(C(C(CO)O)O)=O"]   # threose: 1 step
```

One pass returns the exact optimum for *every* reachable species at once. The
largest network here (Formose at generation 6, 306,244 reactions over 117,874
molecules) solves in under half a second.

## Reproducing the reported results

Each script writes to `ProcessedData/SI/` and is safe to re-run. Every number in
the paper comes from one of them; none is entered by hand.

```bash
uv sync --extra thermo
uv run --extra thermo python scripts/compute_energies.py --workers 8   # ~40 core-hours; resumable
```

Energies first, then everything that reads them:

| Script | Produces |
|---|---|
| `make_si_tables.py` | SI Tables 1-3: routes, steps, pH robustness |
| `make_figure_data.py` | The funnel and sinks/hills tables behind Figures 3 and 4 |
| `make_pathway_energetics.py` | Route length against free energy, for Figure 6 |
| `route_commitment.py` | Reverse/forward flux per route, and each route's least committed step |
| `driving_force.py` | Max-min driving force per route, as a linear program |
| `hyperpath_benchmark.py` | Scale and running time of Algorithm 1 on the full network of each CRNR |
| `amino_acid_coverage.py` | Alpha-amino acids present in each CRNR, against the absent nucleobases |
| `rule_energetics.py` | ΔrG′° by reaction rule and pH, and how far each rule's estimates degenerate |
| `rule_dependence.py` | What survives removing carbonyl migration, the Cannizzaro reaction, or both |
| `matched_controls.py` | Each target against ~50 matched controls, with rank p-values and `control_statistics.csv` |
| `deposition_comparison.py` | Why reactions the earlier deposition called spontaneous are no longer |
| `ms_validation.py` | FT-ICR MS formulas recovered per network, over each network's own mass ceiling. Needs `MS_DATA` set to the unzipped `NHMFLMSData.11.16.20` directory; the archive is deposited beside the manuscript |
| `uronate_case.py` | Yi et al.'s uronate route to pentoses, and the estimator's coverage broken into its causes |
| `ketohexose_case.py` | The hexose branch point against Yi et al.: where the traced route matches their labelling and where minimum step count does not |
| `verify_matches.py` | Re-derive analogue matches and diff against the deposited set |
| `crosscheck_energies.py` | Recompute the G3 energies and diff against the deposited file, to show they are machine-independent |
| `match_formose_g6.py` | Extend Formose matching to generation six, verifying the deposited rows first |
| `bottcher_formose_g6.py` | Böttcher complexity for the extended Formose set |
| `database_matches.py` | Network products found in ChEBI and HMDB, by generation |
| `make_route_specs.py` | The traced routes as autocycle specs, in both bases |
| `kinetic_ordering.py` | Whether the route ordering survives as a kinetic one under Bell-Evans-Polanyi; a deposited negative |
| `descriptor_model.py` | A negative control, not used in the paper; see `ProcessedData/SI/README.md` |

`figures/algorithm1.tex` sets Algorithm 1 of the manuscript with `algorithm2e`;
build it with `xelatex algorithm1.tex` (needs the `algorithm2e` and `newtx`
TeX packages). The body is transcribed from the manuscript, so change it there
and rebuild rather than editing the PDF.

Then the figures:

```bash
uv sync --extra figures
for f in workflow overlap funnel pathway ms rule robustness depth bottcher ketohexose glyoxylate; do
  uv run --extra figures python figures/make_${f}_figure.py
done
uv run --extra figures python figures/make_figure6.py    # tiles the four Formose routes
```

`compute_energies.py` downloads eQuilibrator's ~1.3 GB compound cache on first
use. **No ChemAxon licence is needed**: `chemaxon_bypass` supplies the
protonation handling and dissociation constants come from `pka`. Most of the
runtime is compound resolution, not the energies; each worker resolves the
species in its own block, so eight workers buy rather less than eight times the
speed.

## What is here

| Path | Contents |
|---|---|
| `src/nucleoside_analogues/` | The supported pipeline |
| `scripts/` | Regenerate the deposited tables, energies and validations |
| `figures/` | `acs_style.py`, figure scripts, output in `figures/output/` |
| `tests/` | Chemistry invariants, experimental recall, pathway regression |
| `OriginalData/` | Raw MØD output and the CHO/CHNO analogue libraries. Five networks are used in the manuscript; HCN, Maillard and Urey-Miller are deposited but unused -- see [The three unused networks](#the-three-unused-networks) |
| `ProcessedData/` | Everything derived — see [`ProcessedData/README.md`](ProcessedData/README.md) |
| `notebooks/` | How the published results were produced — see [`notebooks/README.md`](notebooks/README.md) |

| Module | Role |
|---|---|
| `rels` | Read and reshape MØD output; build reaction indices |
| `hyperpath` | Exact shortest pathways, route counts, critical reactions |
| `matching` | Stereochemistry-flattened InChIKey matching |
| `invariants` | Atom/charge balance, generation monotonicity, motif screens |
| `descriptors` | Physicochemical descriptor panel |
| `pka` | pKa assignment and titratability |
| `thermo` | ΔrG′° with propagated uncertainty (needs the `thermo` extra) |
| `chemaxon_bypass` | Lets eQuilibrator build compounds without a ChemAxon licence |

`ProcessedData/` mixes current outputs with superseded ones kept for provenance;
its README says which is which. **Directory names there are load bearing** — the
SI cites `SpontaneousSMILES/` by URL.

## The three unused networks

`OriginalData/` carries three networks the manuscript does not use: Maillard
(glycine + glucose + water), Urey-Miller and HCN. They are deposited whole, and
Maillard in particular is now usable, so this records where it stands.

**Maillard could not be traced at all until recently, for a reason worth
knowing.** The original notebook held a hard-coded table of seed molecules
written in canonical SMILES, while the network itself writes its own. Lookup is
by exact string:

| | Network writes | Canonical is |
|---|---|---|
| glycine | `C(CN)(O)=O` | `NCC(=O)O` |
| glucose, open chain | `C(C(C(C(C(CO)O)O)O)O)=O` | `O=CC(O)C(O)C(O)C(O)CO` |

Same molecules, different strings, so the tracer never recognised its own
starting materials and reached three species -- the seeds alone -- out of 2,923.
`rels.seeds_from_products` fixes this by reading generation 0 from the deposited
product listing instead, and every script here uses it. With the seeds read that
way the spontaneous Maillard network reaches 2,128 species and all five targets.

**What is not done.** Maillard's deposited thermodynamics
(`ProcessedData/RelsWithThermoFiles/MaillardG3RelsWithThermo.tsv`) carry a point
dGr'o with no sigma, so its spontaneous set is the older hard dGr'o < 0 rule
rather than the 95% criterion used throughout this work. Anyone comparing
Maillard against the five must recompute its energies first:

    uv run --extra thermo python scripts/compute_energies.py --workers 8 Maillard

At 8,856 reactions against Formose's 306,244 this is short work. Until it is
done, any Maillard-versus-CRNR comparison is across two different spontaneity
filters and should not be reported.

**Depth.** Maillard and Urey-Miller were expanded only to generation three;
the five networks in the manuscript run to generations four through six. HCN
reaches generation six but builds no nucleobase, and none of the three is
matched against the analogue library in this work.


## Method

```
NA library (MOLGEN, CHO + CHNO)              CRN products (MØD)
   Cl → OH / NH2 substitution                  seeded expansion, G0…G6
              ↓                                          ↓
   collapse stereoisomers  ────►  match on InChIKey[:14]  ◄────┘
                                          ↓
                     eQuilibrator ΔfG′° → ΔrG′° ± σ per reaction
                                          ↓
                          three-way spontaneity call
                                          ↓
                    shortest hyperpath → steps, routes, critical reactions
```

Matching uses the first 14 InChIKey characters, which encode constitution and
charge but not stereochemistry; the networks do not track stereochemistry, so a
match means the network produced *some* stereoisomer of the target scaffold.

A reaction consumes and produces several species, so a network is a directed
hypergraph and a route is a hyperpath. Knuth's generalisation of Dijkstra's
algorithm solves the minimum-weight hyperpath exactly in one pass. Two
objectives measure different things and **reporting either as "steps" requires
saying which**: `chain` (`f = max`) is the longest chain of consecutive
reactions, and reproduces all 1,181 deposited generation-3 values; `reactions`
(`f = sum`) counts every use of a shared intermediate. They differ in 392 of
those 1,181 pathways. See `hyperpath.py` for the derivation and the reason
minimum-cost route counts are well defined where "total pathways" is not.

## Testing

```bash
uv run pytest -q && uv run ruff check src tests && uv run pyright
```

Three tiers: **structural invariants** (every SMILES parses, every reaction
conserves atoms and charge, no reagent from a later generation); **implausibility
screens** (peroxides, N–N and N–O bonds, orthoesters, strained rings — a
blocklist can only show a molecule is not obviously wrong); and **experimental
recall**, the only tier that validates rather than self-checks — the formose
network recovers 19 of 20 structures in the Omran/Decker set, all by generation
three. CI also rejects filenames Windows cannot check out.

## Known limitations

- **pH.** The two ammonia-seeded networks titrate broadly (86% and 57% of
  compounds between pH 7 and 11; 65% and 39% over the generation-3 subset). The
  carbon-only networks are nearly flat but not entirely: Glucose and Pyruvic acid
  each contain CO₂, whose second carbonate constant (10.33) is in range, and 11
  Pyruvic acid reactions change classification by pH 11. Only Formose has no
  titratable species. See `ProcessedData/SI/SI_Table3_pH_robustness.csv`.
- **Coverage, not precision, is the thermodynamic limit.** Every ΔrG′° carries a
  standard error and spontaneity requires the whole 95% interval below zero, but
  σ is small (median 1.6–2.4 kJ/mol) and that criterion agrees with a bare
  `ΔrG′° < 0` on every number reported. What bites is that 9–12% of reactions per
  network — **47% in Pyruvic acid** — fall outside both contribution bases and
  have no estimate at all. These are reported as unestimable, never as zero.
- **Protonation states are assigned, not computed.** `bypass_chemaxon=True` skips
  ChemAxon; cached compounds keep measured constants, the rest are matched
  against a SMARTS table in `pka`, cross-checked against dimorphite-dl in the
  tests.
- **Stereochemistry is discarded** throughout, so `n_chiral_centers`, `fcsp3_bm`
  and the 3D shape pair are computed on an arbitrary stereoisomer.
- **`exact_mass` is monoisotopic**, not average molecular weight — 180.042 against
  180.159 for aspirin. It was called `MW`, which reads as the latter.
- **Deposited pathways were optimised for chain length**, so their reaction counts
  are incidental rather than minimal.

## Related repositories

Two companion repositories in the same organisation:

| Repository | Role |
|---|---|
| [reac-space-exp](https://github.com/Reaction-Space-Explorer/reac-space-exp) | The MØD graph-grammar workflow that generated the networks used here. |
| [autocycle](https://github.com/Reaction-Space-Explorer/autocycle) | Publication figures for reaction cycles and synthetic routes, with 2D structures drawn inline at one constant bond length. |

[ssiddhantsharma/nucleoside-analogues](https://github.com/ssiddhantsharma/nucleoside-analogues)
is the archived 2021–2023 working repository, and the source of `figures/` and
`notebooks/analysis/`.

## Licence

Code under [BSD-3-Clause](LICENSE), data and figures under [CC BY 4.0](LICENSE-DATA).
`descriptors.py` adapts the panel from
[MycoPermeNet](https://github.com/Nevbarunegbe/Mycomembrane-permeability-project) (MIT).

## References

The manuscript carries the full bibliography. These five are the ones a reader
of this code will want:

- **The networks** — Arya *et al.*, *Chem. Sci.* **2022**, *13*, 4838.
  [10.1039/d2sc00256f](https://doi.org/10.1039/d2sc00256f)
- **The analogue space** — Cleaves *et al.*, *J. Chem. Inf. Model.* **2019**, *59*, 4266.
  [10.1021/acs.jcim.9b00632](https://doi.org/10.1021/acs.jcim.9b00632)
- **The shortest-path algorithm** — Knuth, *Inf. Process. Lett.* **1977**, *6*, 1.
  [10.1016/0020-0190(77)90002-3](https://doi.org/10.1016/0020-0190(77)90002-3)
- **Reaction networks as hypergraphs** — Gallo *et al.*, *Discrete Appl. Math.* **1993**, *42*, 177.
  [10.1016/0166-218X(93)90045-P](https://doi.org/10.1016/0166-218X(93)90045-P)
- **The free energies** — Noor *et al.*, *PLoS Comput. Biol.* **2013**, *9*, e1003098.
  [10.1371/journal.pcbi.1003098](https://doi.org/10.1371/journal.pcbi.1003098)
