# TODO

## Running on bizon

- **Four-pH energies** (`compute_energies.py --workers 8`, pH 7.0/7.4/9.0/11.0). Four of five
  networks written: GlucoseAmm, FormoseAmm, Glucose, PyruvicAcid. **Formose G6 is last in the
  queue and still going**, about an hour and a quarter into its own pass; its four files exist but
  are empty, because each is written only when its network finishes. The local pH 7.4 copies are
  intact and committed, so nothing already computed is at risk.

## Blocked on that run

- **SI Table 3** still holds the generation-3 set (582 Formose reactions). Its caption says so and
  must change with it.
- dGr distribution by pH and by reaction rule.
- Max-min driving force on the traced routes. The spike passed at 42-106 kJ/mol on seven routes.
- Descriptors and SHAP, dropping the constant features.

## Open, needing a decision rather than work

- **JCIM Article rather than Letter.** Promote SI Extended Methods 1, 2 (with Algorithm S1), 4 and 6
  into the main Materials and Methods; keep 3 and 5 in the SI; promote Figure S4, the robustness
  grid, now that every reachable target survives removal of the disputed mechanism in every CRNR; split
  Results and Discussion; fix the heading styles, Materials and Methods being `normal` where
  Results and Discussion is `Heading 3`. Not started, awaiting a go.

- **Two orphans.** Figure S3 (the two route costs) and SI Table 4 (sinks and hills) carry captions
  but are referenced nowhere in the body, by any wording. Either cite them or cut them, as was done
  with the old Figure 5 and the Bottcher figure.

- **Which route figures, if any, go to the SI?** The glyoxylate pair is the strongest candidate:
  Formose reaches it through carbon chemistry, FormoseAmm through nitrogen, by transamination. The
  glucose routes are eleven steps and only legible at full width. The FormoseAmm duplicates of the
  four main-text panels would be the same picture twice.

## Open, needing something from outside

- **HMDB.** Download "Structures" (`structures.zip`) from <https://www.hmdb.ca/downloads> and unzip
  to `OriginalData/reference_databases/hmdb_structures.sdf`. The 403 is a Cloudflare challenge
  (`cf-mitigated: challenge`), not a user-agent block, so it needs a browser and cannot be scripted.
  `scripts/database_matches.py` picks the file up automatically and is otherwise finished. ChEBI is
  already re-derived; if HMDB never arrives, the deposited HMDB/KEGG/ECMDB counts stay and are flagged
  as from the original deposition, exactly as KEGG already is.
- **KEGG** cannot be re-derived; its bulk data is licensed. The deposited counts stay, flagged.

## Author-side

- ORCIDs. All eight authors' ORCIDs are already supplied in the docx comments and just need placing
  in the author block: ML 0009-0005-4426-4571, SS 0000-0003-1768-1802, RC 0000-0001-5942-9918,
  AA 0000-0001-6782-407X, DJ 0009-0005-5088-7222, RO 0009-0002-6595-3307, MM 0000-0001-8526-2429,
  HJC 0000-0003-4101-0654. Three further ORCIDs are supplied for people not on the author line --
  Alejandro Lozano Garcia, Jakob Anderson, Christopher Butch -- which needs a decision, not work.
- Comment 8 asks for Algorithm S1 to be reworked as pseudo-code; it bears on promoting it to the
  main Methods.
- TOC graphic.
- Conflict of interest and funding statements.
- Zenodo DOI for the archived repository.
- The cover letter is bundled in the same .docx (p0-p13) and asks for "rapid communication". If the
  paper goes to Article format the letter has to change with it, and most journals want the letter as
  a separate file anyway.
- Put `Jim_NA/` under version control. It holds the manuscript, the red-marked copy and INDEX.md,
  none of it versioned.

## Offshoots this work could support

- **Stochastic kinetics on these networks.** Lauber et al. need two global parameters, not one per
  rule, and take dfG from eQuilibrator exactly as we do. The obstacle is that this repository holds
  MØD's output and not its grammar; the rules live in reac-space-exp, which is ours.
- **Kinetic ordering of the traced routes** once real barriers exist. The parameter-free shortcut was
  tried and does not carry: see `scripts/kinetic_ordering.py`. Bell-Evans-Polanyi with global
  parameters does identify the rate-limiting step with the least committed one, which the paper now
  uses, but the ordering of whole routes is not separable from the drift of a maximum over more steps.
  Real barriers, per rule, are what would settle it.
- **Integer hyperflow**, following Abel et al. (ref 73), asking which traced routes are
  stoichiometrically realizable, closing the caveat rather than stating it.
- **The rule-dependence method itself** applied more widely: removing a disputed mechanism and
  re-searching is general, and the disagreement it settles here is not specific to sugars.

## Settled this round

Six Krishnamurthy/Liotta-school papers read in full. Three cited: Suarez-Marina 2019 as ref 37 (the
canonical monomers are not the favoured products of uncatalysed glycosylation -- the motivation for
the whole enumeration), Oro & Cox 1962 as ref 87 (the actual origin of the acetaldehyde route to
2-deoxyribose), and Yi et al. JACS Au 2023 as ref 96 (the uronate route to pentoses). The Trapp
citation, ref 81, is no longer stated as fact: the review at ref 44 flags exact-mass ion
chromatograms and unexcluded biological contamination, and the text now says so. Rule dependence
split and rewritten, both control bases disclosed, uronate route reported as present-but-unestimable.
103 references, order verified.

Yi et al. 2023 read and used, as ref 88. Carbonyl migration is established for tetroses without
Ca2+, with Ca2+ the switch between the enediol route our rules encode and the 1,2-hydride route, so
the rule-dependence removal is framed as a stress test rather than a rival account. Figure 7 added:
our traced route to the 3-ketohexose is their Scheme 5c exactly, found blind from glycolaldehyde
alone, while the 2-ketohexose exposes the metric -- a step shorter through the aldohexose they do not
observe, their own final migration excluded as a null estimate. The three estimator findings are now
one statement about isomerisation blindness in the Methods. Threose softened to the aldotetrose
constitution in the Figure 5 caption and p96, with their measured half-lives as the cost of
flattening.

## Settled by the science pass

Every quantitative claim in the manuscript was checked against the data. Five were wrong, all
in the same way: the number was computed ad hoc for a draft and deposited nowhere, so nothing
could catch it drifting.

- Figure S2's caption: formose-ammonia ceiling 315 Da (real 312.15), range 45 to 70% (real 50
  to 71 over each network's own ceiling), 16 of 511 formulas (real 16 of 492).
- The Methods claimed peak lists were exported from m/z 150 upward. The lowest exported peak is
  136.03 and six of the eight spectra start at 161.05; 155 is a floor adopted here.
- Figure 4's caption said 40 to 53%; the funnel gives 39.56 to 52.41.
- The permutation p was a Monte Carlo estimate quoted as 0.010. C(25,10) is 3,268,760, small
  enough to enumerate, and the exact value is 0.00742.

Correct and now reproducible rather than only asserted: every per-pair rank p (Phipson-Smyth,
0.020 to 0.042 with glyoxylate in F at 0.137), Fisher's 9.11e-11, the 44% to 82% unestimable
share against the earlier deposition, the 88% control percentile, Figure 2's 33% F-G overlap,
SI Table 3's pH drifts, the 112 route steps, the MDF range and its six-of-sixteen disagreement
with the least-committed step, and -30.3 +/- 2.2 kJ/mol per reaction.

Three scripts now emit what the manuscript quotes: ms_validation.py (ceiling columns, computed
from the products rather than hardcoded), matched_controls.py (p_value plus
control_statistics.csv) and the new deposition_comparison.py. Tests pin all of it.

Route length in Figure 6 is the longest chain, not the reaction count, which is why Glucose
deoxyribose is 25 reactions at depth <= 20. The text now says so.

CI never linted or format-checked figures/. It does now.

## Settled the round before, for the record

98 references, first-citation order verified. Glyoxylate added as a fifth target (Bean 81,
Krishnamurthy & Liotta 83). Sutton 86, Tabata 87, Lauber 88, Cruz 91, Abel 73. Figure 2 rebuilt
from the reaction listings after its Formose row proved stale; the reachability funnel became
Figure 4; Figure 7 uncrowded to two panels; Figure S2 (Böttcher) cut, its content now in the
discussion as percentiles and an AUROC, and the remaining SI figures renumbered S1-S4. Formose G6
matching closed, 3,404 matched species to 9,305. Apiose, ethylene glycol, threitol, glyceric and
glycolic acid were checked for backbone literature and have none, so they are not targets.

Figure 6 is now the four Formose spontaneous routes, tiled from the autocycle renderings so it
cannot drift from the traced data. `make_route_specs.py` takes a basis argument, and the 25
with-unestimable routes are generated too, so the three-reaction glucose derivation of deoxyribose
survives as a figure rather than only as a table row. Every claim in the new caption was checked
against the data: F and FA agree at all four depths, and two of the three steps in that shorter
glucose derivation have no free-energy estimate at all.
