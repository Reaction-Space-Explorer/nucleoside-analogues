# Outstanding work

Kept in the repository so it survives between sessions. Ordered by what blocks
what, not by importance.

## Running now

- [ ] **Free energies at four pH values** (bizon, ~7 h, started). One resolution
      pass serves pH 7.0 / 7.4 / 9.0 / 11.0, since compound resolution dominates
      the cost and does not depend on pH. Verified before launch: pH 7.4
      reproduces the deposited file exactly and all four pH match SI Table 3.

## Blocked on that run

- [ ] Refresh **SI Table 3** at full depth rather than generation 3.
- [ ] **ΔrG distribution by pH and by reaction rule** (Chem. Sci. ESI SI12 style).
      The by-rule panel works now; the by-pH panel needs the run.
- [ ] **Max-min driving force**. Spike passed: the interface works, all sixteen
      routes are realizable, and MDF returns 42-106 kJ/mol on the seven routes
      tested. About a day. Needs pH-specific energies, so it follows the run.
- [ ] **Descriptors + SHAP**, with bootstrap chi-square and KS after Wolos.
      Needs reachability labels. Note: `formal_charge` and `abs_charge` are
      constant in all five networks and `NumAmideBonds` in four, so drop the
      constant features before fitting.

## Ready to start

- [ ] **autocycle for the reaction figures**. Write a converter from
      `hyperpath.trace()` output to autocycle's YAML so Figure 6 regenerates
      from the data rather than being drawn by hand.
- [ ] **Re-derive the database comparison**. ChEBI is open and scriptable
      (`ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi_lite_3_stars.sdf.gz`, 15 MB,
      3-star curated). **HMDB refuses scripted download (403 even with a browser
      user agent)** and must be fetched by hand once. **KEGG bulk data is
      licensed**, so it cannot be re-derived here; either drop it or keep the
      existing counts and say they are from the earlier release.
      The current `DatabaseMatches/` reconcile exactly against our own numbers
      but stop at Formose G5, one generation short of the network.
- [ ] **Sinks and hills**: computed, but only 9-77% of matched species classify
      because so many lack estimable reactions on one side. Report as a table
      with the classified fraction stated; it does not support a figure.
- [ ] Small SI figure: **chain depth versus reaction count** worked example, next
      to Algorithm S1. They differ in 392 of 1,181 deposited pathways and the
      distinction is currently prose only.

## Manuscript, needing the author

- [ ] **ORCIDs** — eight comments (ids 0-7) carry them, still unresolved.
- [ ] **TOC / abstract graphic** — required by ACS, does not exist.
- [ ] **Conflict of interest** and **funding** statements — both missing.
- [ ] **Zenodo deposit** and DOI, for the Data Availability section.
- [ ] Delete the **reviewer-suggestion paragraph** (p4): cover-letter material
      sitting in the manuscript body.
- [ ] Decide what **Figure 5** should show.
- [ ] Put `Jim_NA/` under version control. It holds Algorithm S1, three SI
      tables, four figures and a rewritten Methods, protected only by
      timestamped `.docx` backups.

## Done

Networks complete to their deepest generation (Formose reaches G6 after
`rels_6.txt` was recovered). Matches re-derived and identical. Chemistry audited
across 665,645 reactions, zero unbalanced. Böttcher values reproduce 400/400.
Descriptors verified against independent RDKit calls, with `MW` corrected to
`exact_mass` and an empty-SMILES bug fixed. SI Tables 1 and 2 at full depth.
Matched-control analysis answering the Abstract's question. Figures S1-S3 and
the reachability funnel. References at 92, all cited in ACS order.
