# What this paper taught us

Written after the JCIM nucleoside-analogue manuscript, for the next one. Every
item below is here because it actually went wrong, and the incident is named so
the rule has weight rather than sounding like advice.

---

## 1. Numbers

**Every number in the text must be produced by a deposited script.** Not checked
against one — produced by one. A number that exists only in prose has nothing
holding it in place, and it will drift as the analysis moves under it.

This failed four times in one manuscript:

- The SI caption gave a network ceiling of 315 Da, a range of 45–70% and a count
  of 16 of 511. None of the three was produced by any script. The real values
  were 312.15 Da, 50–71%, and 16 of 492. All three were wrong.
- The permutation *p* was quoted as 0.010 from a Monte Carlo run nobody kept.
  The exact value over all 3,268,760 splits is 0.0074.
- A median was `sorted(x)[len(x)//2]`, the upper-middle order statistic, not the
  median. It moved eight of sixteen per-route values.

**Depositing a number catches errors in the act of depositing it.** Writing the
script that emits "eleven β-decarboxylations" revealed that eleven was the count
across five networks and the sentence said Formose. The draft was wrong; the act
of making it reproducible found it.

**A number copied from a docstring is not a deposited number.** The Methods
called the generation-4 FormoseAmm network, 145,820 reactions over 35,318
species, "the largest network considered here". The figure is real but was
lifted from a docstring describing the cost of the *rels reshape*, a different
operation; the largest network is Formose at generation six, 306,244 reactions
over 117,874 species. Check what the number in the source was measuring.

**Check superlatives against the paper's own tables.** That error needed no
measurement to catch: SI Table 2 already listed 306,244 reactions for Formose,
in the same document whose Methods called 145,820 the largest. Any claim of the
form largest, smallest, only or first should be re-derived from the tables
already in the file, mechanically, at the end of every draft.

**Sanity-check magnitudes against chemistry.** A decarboxylation came back at
−373 kJ/mol where a β-keto acid should be near −20. That mismatch was the signal
the estimate was unusable, and the reason to report it as unestimable rather
than favourable.

---

## 2. Test your own finding before you believe it

A correlation between route length and rate-limiting barrier looked strong
(ρ = +0.75, p = 0.0009) and would have supported "reachable in few steps means
fast". It does not survive the obvious confound: the least exergonic of *n*
steps drifts upward with *n*, so longer routes have weaker weakest steps by
construction. Against a length-preserving null it gives p = 0.24 on one
bookkeeping basis and p = 0.04 on another.

**A result that changes verdict with a defensible bookkeeping choice is not a
result.** It was deposited as a documented negative instead.

The same discipline killed a descriptor model that predicted reachability at
AUROC 0.86–0.93 — until the same descriptors predicted *what the estimator could
evaluate at all* at 0.93–0.96. It was reading the estimator's blind spot, not
chemistry.

**Before reporting a correlation, ask what would produce it if the hypothesis
were false, and simulate that.**

---

## 3. Verifying chemistry

**Check what a rule does, not what it is called.** The rule set had
`Keto-enol migration twice` and `Elimination + enol to keto` grouped as one
"carbonyl migration" family. Computing molecular formulas over every reaction
showed 10,671 of the first preserve formula and 19,131 of the second lose water.
They are different transformations. Split, the conclusion reversed: no target
depended on carbonyl migration at all, and every dependence attributed to it
belonged to the dehydration.

A test now asserts the classification by formula rather than by name.

**Know what your representation throws away.** These networks carry no
stereochemistry, so threose and erythrose are one species. A target named
"threose" is really the aldotetrose constitution, and the two diastereomers
differ about six-fold in half-life. Name nodes by constitution and say so.

**Match the feedstock before comparing to an experiment.** A published "no
pentoses observed" result looked like a contradiction until the feedstock was
checked: glycolaldehyde alone, where our network also has formaldehyde. Not a
discrepancy. The genuine comparison — same feedstock — was a different one.

**Estimators have structural blind spots, and they are describable.** Component
contribution sees a reaction only through the groups it changes, so
isomerisations are invisible to it: 5,429 of 5,466 null estimates in one network
were carbonyl migrations. Two published experimental routes are present in the
networks and excluded by the filter — one by a null estimate, one by unbounded
variance — for want of an estimate rather than for being unfavourable. Say this
plainly; it is a limitation with a mechanism, not a caveat.

**Check the claim in the direction it is written.** The SI said null estimates
arise where reagents and products share a decomposition, "as the keto-enol
migrations here do". The verified fact was the converse: 99.3% of null estimates
are keto-enol migrations. Only 50% of keto-enol migrations are null, so the
sentence overstated a real finding by reading a verified implication backwards.
Almost all A are B does not give almost all B are A, and the number that was
checked is rarely the one the sentence needs.

**Look up species by the spelling the data uses.** Lookup is by exact string.
A canonical SMILES is not the network's spelling, and a target silently reads as
unreachable. Build an explicit canonical → network-spelling map.

**RDKit:** `Chem.MolFromSmiles(x).GetRingInfo().NumRings()` returns 0 — the Mol
is freed before RingInfo is read. Bind the Mol first. This produced six false
"makes no ring" flags.

---

## 4. Document mechanics

These cost more time than the science.

**Renumbering must cover the whole file.** Reference-insertion passes stopped at
the reference list, so the SI's superscripts were never updated. Four insertions
later, eight citations were stale — Knuth cited as 57 when he is 59, eQuilibrator
as 69,70 when it is 70,71, a mass-spec acquisition as 47 when it is 49. They only
surfaced when SI sections were promoted into the main text.

**Resolve stale citations from what the sentence cites, not by arithmetic.** The
shifts were inconsistent (+1 in places, +2 in others) because they had
accumulated across several renumberings. Only content identified them.

**A token scheme must not contain re-matchable text.** Wrapping a search term in
sentinels — `\x00SI Table 2\x00` — still contains "SI Table 2", so the second
pass matched it again and collapsed two tables onto one number. Map to opaque
tokens, then resolve.

**Captions and the things they describe live in different places.** Rewriting a
caption for three panels left the old two-panel image embedded, contradicting
both the caption and the paper's conclusion. Renumbering a table's caption left
its body holding the superseded split. **Hash-compare embedded images against
their source files** after any figure regeneration.

**Orphan checks must not use a filter that hides orphans.** Treating any
paragraph starting "Figure 5" as a caption meant "Figure 5 shows…" was read as a
caption, not a reference, and a genuinely uncited figure passed as cited. Match
`^Figure N\.` with the period.

**Re-check first-citation order after any structural move.** Promoting a Methods
section moves its references thirty paragraphs earlier.

**Two copies diverge.** The moment a manuscript exists in both a `.docx` and
Google Docs, decide which is authoritative. Edits applied to one are lost by the
other.

---

## 5. Structure and placement

**Importance and placement should agree.** Audit this explicitly: the paper's
central quantitative result was an SI table while the main text gave it only in
prose; the strongest robustness result was an SI figure; and the more
informative of two workflow diagrams was in the SI while a concept cartoon held
the main-text slot.

**Cut what nothing depends on.** Two SI items had captions but no reference
anywhere in the body. An attempt to rescue one as evidence for another argument
failed on the arithmetic, which settled it. Keep the data deposited and record
why it was cut.

**After moving text, de-duplicate claim by claim.** Promoting an SI section into
the Methods left four statements of the same two cost definitions -- the Methods
prose, two consecutive paragraphs of the promoted text, and the caption of the
table that reports them -- and three separate citations of Knuth for the same
correctness result. The de-duplication pass after the move checked one claim,
the hypergraph definition, and stopped. Diff the moved text against its new
neighbours sentence by sentence, not once.

**Read the caption before deciding the prose is load-bearing.** The argument for
keeping one of those paragraphs was that the table reported both quantities and
the reader would not otherwise know what they were. The table's own caption
already defined both, in near-identical words, immediately above the numbers.
Whenever the case for keeping text is that a figure or table depends on it, read
that caption first; it is usually where the definition already lives.

**Check every figure and table is referenced, and every reference has a target.**
Both directions, every round.

---

## 6. Reading the literature

**Read the paper, not the abstract.** A citation was carried as established fact
when the authoritative review of that field says of it: claims rest on a highly
sensitive extracted ion chromatogram and biological contamination is not
excluded. The uncontested original for the same chemistry was not cited at all.

**Check whether the work you cite is contested by the people likely to review
you.** Five of six papers in one reading batch came from the same school; one
shared an author with us.

---

## 7. Process

**Do not let a tool's limitation look like a document defect.** A converter
flattened correctly formatted pseudocode and dropped table borders; the
manuscript was fine. Check the source before reporting a problem in it.

**Run the repo's own gate, not just the part of it you were thinking about.**
Six commits went out with CI red. The cause was `ruff format --check` on two new
scripts -- a wrapped `print`, nothing else -- while `pytest` passed locally every
time. The gate ran four steps and only one was being checked, so the failures
died in 28 seconds before the tests it was passing ever executed. Run the
workflow's commands, in its order, before pushing; and read the run afterwards
rather than assuming a local pass covers it.

**Scope discipline.** A request for a PDF became a 700 MB install and a mirror
hunt, after the requirement that motivated it had been dropped. Re-read what is
actually being asked before escalating the means.
