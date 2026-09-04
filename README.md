# nucleoside-analogues

[![CI](https://github.com/Reaction-Space-Explorer/nucleoside-analogues/actions/workflows/ci.yml/badge.svg)](https://github.com/Reaction-Space-Explorer/nucleoside-analogues/actions/workflows/ci.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/code-BSD--3--Clause-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)

Which nucleic-acid-like molecules can prebiotic chemistry actually reach?

This repository intersects combinatorially enumerated **nucleoside analogue (NA)
libraries** with the molecular diversity produced by in-silico **chemical
reaction networks (CRNs)**, then asks whether each matched structure can be
reached by a sequence of thermodynamically spontaneous reactions — and in how
few steps.

It supports *"Computational Identification of Plausible Prebiotic Nucleic Acid
Analogue Synthetic Pathways"* (Ludwig, Sharma, Cruz-Simbron, Arya, Jain, Okoń,
Meringer & Cleaves).

---

## Quick start

```bash
git clone https://github.com/Reaction-Space-Explorer/nucleoside-analogues
cd nucleoside-analogues
uv sync --group dev
uv run pytest -q
```

No `uv`? `pip install -r requirements.txt`, or use the container:

```bash
docker build -t nucleoside-analogues .
docker run --rm nucleoside-analogues pytest -q
```

### Shortest pathways to every reachable molecule

```python
from nucleoside_analogues.hyperpath import shortest_pathways
from nucleoside_analogues.rels import build_index, load_processed_rels, read_products, seeds_from_products

rels = load_processed_rels(
    "ProcessedData/SpontaneousRelsWithThermoFiles/SpontaneousFormoseG3RelsWithThermo.tsv"
)
products = read_products("OriginalData/OriginalNetworkData/Products/formose_output.tsv")

result = shortest_pathways(build_index(rels), seeds_from_products(products))
result.cost["C(C(C(CO)O)O)=O"]   # threose: 1 step
```

One pass returns the exact optimum for *every* reachable species at once. The
largest network here (FormoseAmm at generation 4 — 145,820 reactions, 35,318
molecules) solves in about 0.25 s.

---

## What is here

| Path | Contents |
|---|---|
| `src/nucleoside_analogues/` | The supported pipeline (see below) |
| `tests/` | Chemistry invariants, experimental recall, pathway regression |
| `OriginalData/OriginalNetworkData/` | Raw MØD output: product listings and `Rels_N.tsv` reaction files |
| `OriginalData/OriginalNucleosideAnalogueData/` | CHO and CHNO analogue libraries |
| `ProcessedData/` | Matches, reshaped reactions, thermodynamics, traced pathways |
| `ProcessedData/manifests/` | Filename → SMILES map for the pathway outputs |
| `figures/` | 173 plots, salvaged from the archived working repository |
| `notebooks/pipeline/` | The generation-3 analysis as originally run |
| `notebooks/analysis/` | Descriptor, plotting and database notebooks, salvaged from the archive |

The notebooks are a record of how the published results were produced, not a
runnable pipeline: none of their data paths resolve from a clone, and most
address Google Colab mount points. See [`notebooks/README.md`](notebooks/README.md).

### The package

| Module | Role |
|---|---|
| `rels` | Read and reshape MØD output; build reaction indices |
| `hyperpath` | Exact shortest pathways as minimum-weight hyperpaths |
| `matching` | Stereochemistry-flattened InChIKey matching |
| `invariants` | Atom/charge balance, generation monotonicity, motif screens |
| `descriptors` | Physicochemical descriptor panel |
| `pka` | pKa assignment and titratability classification |
| `thermo` | ΔrG′° with propagated uncertainty (optional `thermo` extra) |

---

## Method

```
NA library (MOLGEN, CHO + CHNO)              CRN products (MØD)
   Cl → OH / NH2 substitution                  seeded expansion, G0…G6
              ↓                                          ↓
   collapse stereoisomers  ────►  match on InChIKey[:14]  ◄────┘
                                          ↓
                                  NA matches per network
                                          ↓
                     eQuilibrator ΔfG′° → ΔrG′° per reaction
                                          ↓
                          spontaneity filter (ΔrG′° < 0)
                                          ↓
                    shortest hyperpath → steps, ΔrG′°, route count
```

Matching uses the **first 14 characters of the InChIKey**, which encode
constitution and charge but not stereochemistry. The network expansions do not
track stereochemistry, so both sides are flattened to constitutional isomers. A
match therefore means the network produced *some* stereoisomer of the target
scaffold.

### Why hyperpaths

A reaction consumes several species and produces several species, so a reaction
network is a directed **hypergraph** and a synthetic route is a *hyperpath*. The
shortest one satisfies

```
d(m) = min over reactions r producing m of ( 1 + f{ d(x) : x in reagents(r) } )
```

with `d(seed) = 0`. For monotone non-decreasing `f` this is a superior function,
so Knuth's generalisation of Dijkstra's algorithm solves it exactly in one pass,
and cycles are excluded structurally.

Two objectives, measuring different things:

- **`chain`** (`f = max`) — longest chain of consecutive reactions. This is the
  derivation-tree depth, and is what the original notebooks reported as
  "Pathway Length". Reproduces all 1,181 deposited generation-3 values.
- **`reactions`** (`f = sum`) — total reactions, counting a shared intermediate
  once per use. An upper bound on the number of *distinct* reactions; minimising
  distinct reactions permits re-use of a sub-derivation, is not a superior
  function, and is not solved here.

**Reporting either number as "steps" requires saying which one.** They differ
for any convergent synthesis — in 392 of the 1,181 deposited pathways (33%).

### Counting routes

`count_minimal_routes` counts distinct *minimum-cost* derivations. This is well
defined because the cost ordering makes the minimum-cost sub-hypergraph acyclic.
"Total number of pathways" is **not** well defined without also fixing a length
bound and a convention for reordering independent branches; different
conventions disagree by many orders of magnitude on the same network.

---

## Testing

```bash
uv run pytest -q          # 87 tests against the deposited data
uv run ruff check src tests
uv run pyright
```

Three tiers, in increasing order of what they can establish:

1. **Structural invariants** — every SMILES parses, every reaction conserves
   atoms and charge, every species appears in the product listing, no reaction
   consumes a reagent from a later generation. All pass on the deposited data;
   they exist to fail loudly if a rule-set change alters the chemistry.
2. **Implausibility screens** — peroxides, N–N and N–O bonds, orthoesters,
   strained rings. A blocklist cannot show a molecule is real, only that it is
   not obviously wrong.
3. **Experimental recall** — the formose network recovers **19 of 21**
   structures in the Omran/Decker literature set, all within three generations.
   This is the only tier that validates rather than self-checks.

CI additionally rejects any filename Windows cannot check out.

---

## Known limitations

- **pH range.** The two ammonia-seeded networks titrate broadly: 86% of all
  compounds in Formose/NH₃ and 57% in Glucose/NH₃ change protonation between
  pH 7 and 11 (65% and 39% over the generation-3 subset the reported energies
  cover), through aliphatic amines with pKa ~9-11. The carbon-only networks are
  nearly flat, but not entirely: Glucose and Pyruvic acid each contain **one**
  titratable species, CO₂, whose second carbonate constant (10.33) falls in
  range. That is enough to move Pyruvic acid — 11 of its reactions change
  classification by pH 11 — so it is wrong to call those networks
  pH-independent. Formose alone has no titratable species at all. See
  `ProcessedData/SI/SI_Table3_pH_robustness.csv`; constants come from
  `nucleoside_analogues.pka` and are cross-checked against dimorphite-dl in the
  test suite.
- **Uncertainty is carried, but coverage is the real limit.** Every ΔrG′° has a
  standard error attached, and a reaction counts as spontaneous only where the
  whole 95% interval lies below zero. The errors are small — median σ is
  1.6-2.4 kJ/mol — so that criterion and a bare `ΔrG′° < 0` agree on every
  number reported here. What does bite is coverage: 9-12% of reactions per
  network, and **47% in Pyruvic acid**, fall outside the span of both the
  reactant- and group-contribution bases and have no estimate at all. Those are
  reported as unestimable, never as zero. Values are in **kJ/mol**.
- **Protonation states are assigned, not computed.** `bypass_chemaxon=True`
  skips ChemAxon's per-compound pKa determination, so the calculation needs no
  licence. Compounds already in the eQuilibrator cache keep their measured
  constants; the rest are matched against a small SMARTS table of literature
  values in `nucleoside_analogues.pka`.
- **Stereochemistry is discarded** throughout. Descriptors that depend on it —
  `n_chiral_centers`, `fcsp3_bm`, and the optional 3D shape pair — are computed
  on an arbitrary stereoisomer.
- **Deposited pathways were optimised for chain length, not reaction count**, so
  their reaction counts are incidental rather than minimal.
- The Urey-Miller and HCN networks in `OriginalData/` are **not used in the
  manuscript**. The network labelled Urey-Miller is seeded with HCN,
  cyanoacetylene, formaldehyde and ammonia rather than the CH₄/NH₃/H₂ spark
  mixture, and recovers none of the amino acids that define that experiment.

---

## Related repositories

| Repository | Role |
|---|---|
| [reac-space-exp](https://github.com/Reaction-Space-Explorer/reac-space-exp) | The MØD graph-grammar workflow that generated these networks (Arya et al., *Chem. Sci.* 2022) |
| [ssiddhantsharma/nucleoside-analogues](https://github.com/ssiddhantsharma/nucleoside-analogues) | Archived 2021–2023 working repository; source of `figures/` and `notebooks/legacy/` |
| [match-viz](https://github.com/Reaction-Space-Explorer/match-viz) | Interactive viewer for matched structures |

## Licence

Code and notebooks under [BSD-3-Clause](LICENSE); data and figures under
[CC BY 4.0](LICENSE-DATA).

`src/nucleoside_analogues/descriptors.py` adapts the descriptor panel from the
[MycoPermeNet project](https://github.com/Nevbarunegbe/Mycomembrane-permeability-project)
(MIT).

## References

- Arya, A. *et al.* An open source computational workflow for the discovery of
  autocatalytic networks in abiotic reactions. *Chem. Sci.* **2022**, *13*,
  4838–4853. <https://doi.org/10.1039/d2sc00256f>
- Cleaves, H. J. *et al.* One Among Millions: The Chemical Space of Nucleic
  Acid-Like Molecules. *J. Chem. Inf. Model.* **2019**, *59*, 4266–4277.
- Andersen, J. L.; Flamm, C.; Merkle, D.; Stadler, P. F. A Software Package for
  Chemically Inspired Graph Transformation. *Lect. Notes Comput. Sci.* **2016**,
  *9761*, 73–88.
- Knuth, D. E. A generalization of Dijkstra's algorithm. *Inf. Process. Lett.*
  **1977**, *6*, 1–5.
- Gallo, G.; Longo, G.; Pallottino, S.; Nguyen, S. Directed hypergraphs and
  applications. *Discrete Appl. Math.* **1993**, *42*, 177–201.
