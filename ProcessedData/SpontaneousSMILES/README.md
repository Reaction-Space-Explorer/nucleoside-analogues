# Spontaneously reachable matched analogues (SI File 1)

Regenerate with `uv run python scripts/spontaneous_smiles.py`. Nothing here is
edited by hand.

Each network is taken at the deepest generation it was generated to, and a
reaction counts as spontaneous only when the whole 95% interval of its
dGr'o lies below zero. Reactions returning a null estimate of exactly 0 +/- 0,
and those whose variance is the infinite-variance sentinel, carry no information
about their own sign and are excluded rather than assigned a value.

| File | Contents |
|---|---|
| `AllSpontaneousSmiles.tsv` | Every matched analogue reached by spontaneous reactions, all five CRNRs |
| `<Network>.tsv` | The same, per network |
| `AllSpontaneousSmiles_with_unestimable.tsv` | The same on the permissive basis, which additionally allows reactions with no usable estimate to pass |
| `SpontVsNonSpont.tsv` | Matched analogues per generation, spontaneous against not |
| `UniqueSmiles.tsv` | Those reached spontaneously in exactly one network |
| `UniqueSmiles.ipynb` | The 2023 notebook that built the first version of this listing, kept as provenance |

The counts here reproduce the `reachable_estimable_only` column of
`ProcessedData/SI/figure_funnel.csv`, which is what Figure 4 plots: 4,108 in F,
2,337 in FA, 2,233 in G, 2,396 in GA and 1 in PA, 11,075 in all.

An earlier version of this listing was built in 2023 from the generation-3
networks under a hard dGr'o < 0 test on the point estimate, and held 974
entries. It did not describe the method the manuscript reports.
